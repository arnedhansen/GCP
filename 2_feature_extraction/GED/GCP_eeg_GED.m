%% GCP Trial-Level GED
%
% Identical spatial filter approach as GCP_eeg_GED_subjects.m, but:
%   - Peak detection is done on INDIVIDUAL TRIALS rather than
%     condition-averaged spec
%
% Pipeline:
%   Phase 1 — Pool all conditions -> broadband GED -> common spatial filter
%   Phase 2 — For each condition & trial, narrowband scan (30-90 Hz),
%             compute per-trial power ratio, detrend, detect peaks.
%             Single-peak and dual-peak (hard 50 Hz) models per trial.
%   Method comparison A — Four benchmark branches:
%             1. Raw channel-space reference
%             2. Top GED component (highest eigenvalue)
%             3. Combined GED components before permutation/CV
%             4. Combined GED components after permutation+CV
%   Method comparison B — Three dual-peak assignment methods:
%             1. Hard 50 Hz boundary
%             2. Per-subject spectral-trough boundary
%             3. Gaussian Mixture Model on single-peak frequencies
%   Aggregation — mean & median peak frequency per condition per subject

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');

%%
nSubj = length(subjects);

% Time windows
baseline_window = [-1.5, -0.25];
stimulus_window = [0, 2.0];

% Gamma frequency range
gamma_range = [30, 90];

% Narrowband scanning parameters
scan_freqs = 30:1:90;
nFreqs     = length(scan_freqs);
scan_width = 3;

% GED parameters
lambda = 0.01;
ged_search_n = 10;          % search first N GED components
topo_display_n = 10;        % number of component topos to visualize
template_front_weight = 0.7; % anti-template weight for frontal channels
template_sigma_occ = 0.20;   % spatial smoothness for occipital template
template_sigma_front = 0.25; % spatial smoothness for frontal anti-template
score_w_corr = 2;           % composite-score weight for template correlation
score_w_gamma = 1.75;       % composite-score weight for pooled gamma evidence
score_w_eval = 1.5;         % composite-score weight for GED eigenvalue evidence
score_w_frontleak = 1;      % composite-score penalty for frontal leakage
viz_suppress_nonocc_outliers = false;  % visualization-only suppression/interpolation
viz_interp_k = 6;                   % nearest neighbors for channel interpolation
viz_nonocc_outlier_mult = 1.00;     % non-occipital outlier threshold multiplier (vs posterior pctl)
viz_topo_prctile = 99;              % robust percentile for color scaling
min_occfront_ratio = 0.50;          % hard eligibility threshold (occ/front ratio)
min_eigval_hard = 1.05;             % hard minimum GED eigenvalue
min_corr_hard = 0.1;               % hard minimum template correlation
min_gamma_hard = 0.05;              % hard minimum gamma evidence (log ratio)
cv_gain_min = 0.01;                 % minimum positive held-out gain to add component

% Multi-component selection (perm + CV) and weighted projection
comp_select_mode = 'perm_cv';
projection_mode = 'weighted_power';  % fixed in this implementation
n_perm = 200;                        % permutation count (balanced runtime/stability)
perm_alpha = 0.05;                   % family-wise alpha for max-eigen threshold
cv_folds = 10;                       % strict mode: K folds
cv_repeats = 2;                      % strict mode: repeated K-fold
random_seed = 13;                    % reproducible randomization
cv_min_train_trials = 10;            % minimum train trials per CV fold
cv_min_test_trials = 3;              % minimum test trials per CV fold

% Four-way benchmark config (ordered for plotting/metrics):
% raw -> top component -> pre-Perm+CV combined -> post-Perm+CV combined
benchmark_methods = {'raw', 'ged_top_eig', 'ged_preperm_hard_weighted', 'ged_permcv_weighted'};
nBenchmarkMethods = numel(benchmark_methods);
raw_reference_definition = 'posterior_roi';  % locked reference for raw benchmark branch
compute_detectability = true;
compute_separation = true;
compute_reliability = true;

% Detrending parameters for power-ratio spectrum
poly_order = 2;
centroid_freq_range = [40 80];
centroid_band_mask = scan_freqs >= centroid_freq_range(1) & scan_freqs <= centroid_freq_range(2);
centroid_posfrac_min = 0.20; % minimum positive-energy fraction for centroid validity
centroid_min_peak = 0.02;    % minimum positive peak in detrended spectrum

% Condition info
condNames  = {'c25', 'c50', 'c75', 'c100'};
trialCodes = [61, 62, 63, 64];
condLabels = {'25%', '50%', '75%', '100%'};

% Figure save directories
if ispc
    fig_save_dir = 'W:\Students\Arne\GCP\figures\eeg\ged';
else
    fig_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged';
end
if ~exist(fig_save_dir, 'dir'), mkdir(fig_save_dir); end
comp_sel_save_dir = fullfile(fig_save_dir, 'component_selection');
if ~exist(comp_sel_save_dir, 'dir'), mkdir(comp_sel_save_dir); end

method_base = fullfile(fig_save_dir, 'method_comparison');
method_dirs = {fullfile(method_base, 'hard_50Hz'), ...
    fullfile(method_base, 'spectral_trough'), ...
    fullfile(method_base, 'gmm')};
for mi = 1:3
    if ~exist(method_dirs{mi}, 'dir'), mkdir(method_dirs{mi}); end
end

%% Preallocate storage
all_trial_powratio     = cell(4, nSubj);
all_trial_peaks_single = cell(4, nSubj);
all_trial_peaks_low    = cell(4, nSubj);
all_trial_peaks_high   = cell(4, nSubj);
all_trial_centroid     = cell(4, nSubj);

all_trial_mean_single   = nan(4, nSubj);
all_trial_median_single = nan(4, nSubj);
all_trial_mean_low      = nan(4, nSubj);
all_trial_median_low    = nan(4, nSubj);
all_trial_mean_high     = nan(4, nSubj);
all_trial_median_high   = nan(4, nSubj);
all_trial_median_gap    = nan(4, nSubj);
all_trial_mean_centroid   = nan(4, nSubj);
all_trial_median_centroid = nan(4, nSubj);

all_trial_detrate_single = nan(4, nSubj);
all_trial_detrate_low    = nan(4, nSubj);
all_trial_detrate_high   = nan(4, nSubj);
all_trial_detrate_centroid = nan(4, nSubj);

all_topos       = cell(1, nSubj);
all_topo_labels = cell(1, nSubj);
all_eigenvalues = nan(1, nSubj);
all_selected_comp_idx  = nan(1, nSubj);
all_selected_comp_corr = nan(1, nSubj);
all_selected_comp_eval = nan(1, nSubj);
all_top5_corrs         = nan(5, nSubj);
all_top5_evals         = nan(5, nSubj);
all_top5_topos         = cell(1, nSubj);
all_simulated_templates = cell(1, nSubj);
all_selected_comp_indices_multi = cell(1, nSubj);
all_selected_comp_weights = cell(1, nSubj);
all_preperm_comp_indices_multi = cell(1, nSubj);
all_preperm_comp_weights = cell(1, nSubj);
all_component_selection_stats = cell(1, nSubj);
all_cv_curves = cell(1, nSubj);
all_perm_thresholds = nan(1, nSubj);

all_trial_powratio_bench = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_dt_bench = cell(nBenchmarkMethods, 4, nSubj);
benchmark_metric_detectability = nan(nBenchmarkMethods, 4, nSubj);
benchmark_metric_prominence = nan(nBenchmarkMethods, 4, nSubj);
benchmark_metric_separation_slope = nan(nBenchmarkMethods, nSubj);
benchmark_metric_separation_delta = nan(nBenchmarkMethods, nSubj);
benchmark_metric_reliability_trialcv = nan(nBenchmarkMethods, 4, nSubj);
benchmark_metric_reliability_subjspread = nan(nBenchmarkMethods, 4);

chanlocs_all = {};

%% Process each subject
for subj = 1:nSubj
    close all
    tic
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load dataEEG

    fsample = dataEEG_c25.fsample;

    trialIndices = { ...
        find(dataEEG_c25.trialinfo  == 61), ...
        find(dataEEG_c50.trialinfo  == 62), ...
        find(dataEEG_c75.trialinfo  == 63), ...
        find(dataEEG_c100.trialinfo == 64)};
    dataStructs = {dataEEG_c25, dataEEG_c50, dataEEG_c75, dataEEG_c100};

    if subj == 1
        chanlocs_all = dataEEG_c25.label;
    end

    nChans = length(dataEEG_c25.label);

    occ_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO)', 'once')), dataEEG_c25.label);
    occ_idx  = find(occ_mask);
    nOcc     = length(occ_idx);
    if nOcc == 0
        warning('No occipital channels matched for subject %s. Using all channels as template fallback.', subjects{subj});
        occ_idx = 1:nChans;
        nOcc = nChans;
    end
    front_mask = cellfun(@(l) ~isempty(regexp(l, '^(Fp|AF|F)', 'once')), dataEEG_c25.label);
    front_idx  = find(front_mask);
    if isempty(front_idx)
        warning('No frontal channels matched for subject %s. Frontal penalty disabled for this subject.', subjects{subj});
    end
    post_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|P)', 'once')), dataEEG_c25.label);
    post_idx  = find(post_mask);
    if isempty(post_idx)
        post_idx = occ_idx;
    end
    switch lower(raw_reference_definition)
        case 'posterior_roi'
            raw_w = zeros(nChans, 1);
            raw_w(post_idx) = 1;
            if sum(raw_w) > 0
                raw_w = raw_w / sum(raw_w);
            else
                raw_w = ones(nChans, 1) / nChans;
            end
        otherwise
            error('Unsupported raw_reference_definition: %s', raw_reference_definition);
    end
    nonocc_idx = setdiff(1:nChans, occ_idx);

    %% ================================================================
    %  PHASE 1: Build POOLED covariance across all conditions -> one GED
    %  ================================================================
    clc
    fprintf('Subject %s (%d/%d) — Phase 1: Full-head GED + occipital template (%d occ / %d ch)\n', ...
        subjects{subj}, subj, nSubj, nOcc, nChans);
    rng(random_seed + subj, 'twister');

    covStim_full = zeros(nChans);
    covBase_full = zeros(nChans);
    nTrials_total = 0;
    stim_cov_trials = {};
    base_cov_trials = {};

    dat_per_cond = cell(1, 4);

    for cond = 1:4
        dat    = dataStructs{cond};
        trlIdx = trialIndices{cond};
        if isempty(trlIdx), continue; end

        cfg = [];
        cfg.trials = trlIdx;
        dat = ft_selectdata(cfg, dat);
        dat_per_cond{cond} = dat;

        cfg_filt = [];
        cfg_filt.bpfilter   = 'yes';
        cfg_filt.bpfreq     = gamma_range;
        cfg_filt.bpfilttype = 'fir';
        cfg_filt.bpfiltord  = round(3 * fsample / gamma_range(1));
        dat_gamma = ft_preprocessing(cfg_filt, dat);

        cfg_t = [];
        cfg_t.latency = baseline_window;
        dat_base = ft_selectdata(cfg_t, dat_gamma);

        cfg_t.latency = stimulus_window;
        dat_stim = ft_selectdata(cfg_t, dat_gamma);

        nTrl = length(dat_stim.trial);
        for trl = 1:nTrl
            d = double(dat_stim.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            cov_stim_trl = (d * d') / size(d, 2);
            covStim_full = covStim_full + cov_stim_trl;
            stim_cov_trials{end + 1} = cov_stim_trl; %#ok<AGROW>

            d = double(dat_base.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            cov_base_trl = (d * d') / size(d, 2);
            covBase_full = covBase_full + cov_base_trl;
            base_cov_trials{end + 1} = cov_base_trl; %#ok<AGROW>
        end
        nTrials_total = nTrials_total + nTrl;
    end

    covStim_full = covStim_full / nTrials_total;
    covBase_full = covBase_full / nTrials_total;

    % Shrinkage regularization (full-head GED)
    covStim_reg = (1-lambda)*covStim_full + lambda*mean(diag(covStim_full))*eye(nChans);
    covBase_reg = (1-lambda)*covBase_full + lambda*mean(diag(covBase_full))*eye(nChans);

    % GED on full-head covariance
    [W_full, D_full] = eig(covStim_reg, covBase_reg);
    [evals_sorted, sortIdx] = sort(real(diag(D_full)), 'descend');
    W_full = W_full(:, sortIdx);

    nSearch = min(ged_search_n, size(W_full, 2));
    searchFilters = nan(nChans, nSearch);
    searchTopos = nan(nChans, nSearch);
    searchCorrs = nan(nSearch, 1);
    searchScores = nan(nSearch, 1);
    searchOccStrength = nan(nSearch, 1);
    searchFrontStrength = nan(nSearch, 1);
    searchOccFrontRatio = nan(nSearch, 1);
    searchGammaEvidence = nan(nSearch, 1);
    searchFrontLeak = nan(nSearch, 1);

    % Simulated signed occipital template with frontal anti-template.
    % Spatial smoothing uses layout coordinates when available.
    sim_template = zeros(nChans, 1);
    lay_labels = headmodel.layANThead.label;
    lay_pos = headmodel.layANThead.pos;
    chan_pos = nan(nChans, 2);
    for ch = 1:nChans
        li = find(strcmp(lay_labels, dataEEG_c25.label{ch}), 1, 'first');
        if ~isempty(li)
            chan_pos(ch, :) = lay_pos(li, :);
        end
    end
    has_pos = ~any(isnan(chan_pos), 2);
    occ_pos_idx = intersect(occ_idx, find(has_pos));
    front_pos_idx = intersect(front_idx, find(has_pos));
    if ~isempty(occ_pos_idx) && ~isempty(front_pos_idx)
        occ_ctr = mean(chan_pos(occ_pos_idx, :), 1);
        front_ctr = mean(chan_pos(front_pos_idx, :), 1);
        occ_d2 = sum((chan_pos - occ_ctr).^2, 2);
        front_d2 = sum((chan_pos - front_ctr).^2, 2);
        w_occ = exp(-occ_d2 / (2 * template_sigma_occ^2));
        w_front = exp(-front_d2 / (2 * template_sigma_front^2));
        sim_template = w_occ - template_front_weight * w_front;
        sim_template(~has_pos) = 0;
    else
        sim_template(occ_idx) = 1;
        sim_template(front_idx) = -template_front_weight;
    end
    if std(sim_template) > 0
        sim_template = (sim_template - mean(sim_template)) / std(sim_template);
    end

    % Full-head forward model for topoplot and component scoring
    covStim_full_reg = (1-lambda)*covStim_full + lambda*mean(diag(covStim_full))*eye(nChans);
    for ci = 1:nSearch
        w_ci = W_full(:, ci);
        topo_ci = covStim_full_reg * w_ci;
        r_ci = corr(topo_ci, sim_template, 'rows', 'complete');
        if ~isnan(r_ci) && r_ci < 0
            w_ci = -w_ci;
            topo_ci = -topo_ci;
            r_ci = -r_ci;
        end
        occ_strength = mean(abs(topo_ci(occ_idx)));
        if ~isempty(front_idx)
            front_strength = mean(abs(topo_ci(front_idx)));
            ratio_ci = occ_strength / max(front_strength, eps);
            front_leak_ci = front_strength / max(occ_strength, eps);
        else
            front_strength = 0;
            ratio_ci = Inf;
            front_leak_ci = 0;
        end
        % Pooled gamma evidence (all conditions combined).
        gamma_ev_ci = log( max((w_ci' * covStim_reg * w_ci), eps) / ...
                           max((w_ci' * covBase_reg * w_ci), eps) );

        searchFilters(:, ci) = w_ci;
        searchTopos(:, ci) = topo_ci;
        searchCorrs(ci) = r_ci;
        searchOccStrength(ci) = occ_strength;
        searchFrontStrength(ci) = front_strength;
        searchOccFrontRatio(ci) = ratio_ci;
        searchGammaEvidence(ci) = gamma_ev_ci;
        searchFrontLeak(ci) = front_leak_ci;
    end
    W_top = searchFilters(:, 1);

    % Candidate objectives (used for ranking after hard eligibility gates).
    corr_vec = searchCorrs;
    ratio_vec = searchOccFrontRatio;
    gamma_vec = searchGammaEvidence;
    eval_vec = log(max(evals_sorted(1:nSearch), eps));
    leak_vec = searchFrontLeak;

    valid_corr = isfinite(corr_vec);
    valid_gamma = isfinite(gamma_vec);
    valid_eval = isfinite(eval_vec);
    valid_leak = isfinite(leak_vec);

    if any(valid_corr), corr_mu = mean(corr_vec(valid_corr)); corr_sd = std(corr_vec(valid_corr));
    else, corr_mu = 0; corr_sd = 1; end
    if any(valid_gamma), gamma_mu = mean(gamma_vec(valid_gamma)); gamma_sd = std(gamma_vec(valid_gamma));
    else, gamma_mu = 0; gamma_sd = 1; end
    if any(valid_eval), eval_mu = mean(eval_vec(valid_eval)); eval_sd = std(eval_vec(valid_eval));
    else, eval_mu = 0; eval_sd = 1; end
    if any(valid_leak), leak_mu = mean(leak_vec(valid_leak)); leak_sd = std(leak_vec(valid_leak));
    else, leak_mu = 0; leak_sd = 1; end

    z_corr = (corr_vec - corr_mu) / max(corr_sd, eps);
    z_gamma = (gamma_vec - gamma_mu) / max(gamma_sd, eps);
    z_eval = (eval_vec - eval_mu) / max(eval_sd, eps);
    z_leak = (leak_vec - leak_mu) / max(leak_sd, eps);

    comp_score_raw = score_w_corr * z_corr + ...
                     score_w_gamma * z_gamma + ...
                     score_w_eval * z_eval - ...
                     score_w_frontleak * z_leak;
    comp_score = comp_score_raw;
    finite_metrics = isfinite(corr_vec) & isfinite(ratio_vec) & isfinite(gamma_vec) & ...
        isfinite(evals_sorted(1:nSearch));
    hard_eligible = finite_metrics & ...
        (evals_sorted(1:nSearch) >= min_eigval_hard) & ...
        (corr_vec >= min_corr_hard) & ...
        (ratio_vec >= min_occfront_ratio) & ...
        (gamma_vec >= min_gamma_hard);
    searchScores = comp_score;
    searchScores(~hard_eligible) = -Inf;
    if ~any(isfinite(searchScores))
        warning(['No components met hard thresholds for subject %s. ', ...
                 'Falling back to unconstrained composite ranking.'], subjects{subj});
        searchScores = comp_score;
    end
    [bestScore, bestIdx] = max(searchScores);
    if isempty(bestIdx) || isnan(bestScore)
        bestIdx = 1;
        bestScore = NaN;
    end

    candidate_table = struct();
    candidate_table.comp_idx = (1:nSearch)';
    candidate_table.eigenvalue = evals_sorted(1:nSearch);
    candidate_table.log_eigenvalue = eval_vec;
    candidate_table.corr = corr_vec;
    candidate_table.gamma = gamma_vec;
    candidate_table.ratio = ratio_vec;
    candidate_table.front_leak = leak_vec;
    candidate_table.score = searchScores;
    candidate_table.hard_eligible = hard_eligible;

    perm_threshold = NaN;
    perm_null_max = nan(n_perm, 1);
    perm_null_max_unconstrained = nan(n_perm, 1);
    perm_sig_mask = true(nSearch, 1);
    if strcmpi(comp_select_mode, 'perm_cv')
        nStimTrials = numel(stim_cov_trials);
        nBaseTrials = numel(base_cov_trials);
        nPool = nStimTrials + nBaseTrials;
        if nStimTrials > 0 && nBaseTrials > 0
            all_cov_trials = [stim_cov_trials, base_cov_trials];
            for pi = 1:n_perm
                pidx = randperm(nPool);
                stim_idx_perm = pidx(1:nStimTrials);
                base_idx_perm = pidx(nStimTrials+1:end);
                covStim_perm = zeros(nChans);
                covBase_perm = zeros(nChans);
                for ti = 1:numel(stim_idx_perm)
                    covStim_perm = covStim_perm + all_cov_trials{stim_idx_perm(ti)};
                end
                for ti = 1:numel(base_idx_perm)
                    covBase_perm = covBase_perm + all_cov_trials{base_idx_perm(ti)};
                end
                covStim_perm = covStim_perm / numel(stim_idx_perm);
                covBase_perm = covBase_perm / numel(base_idx_perm);
                covStim_perm_reg = (1-lambda)*covStim_perm + lambda*mean(diag(covStim_perm))*eye(nChans);
                covBase_perm_reg = (1-lambda)*covBase_perm + lambda*mean(diag(covBase_perm))*eye(nChans);
                [W_perm, D_perm] = eig(covStim_perm_reg, covBase_perm_reg);
                [eval_perm, perm_ord] = sort(real(diag(D_perm)), 'descend');
                W_perm = W_perm(:, perm_ord);
                nPermSearch = min(nSearch, numel(eval_perm));
                eval_perm = eval_perm(1:nPermSearch);
                perm_null_max_unconstrained(pi) = max(eval_perm);
                corr_perm = nan(nPermSearch, 1);
                ratio_perm = nan(nPermSearch, 1);
                gamma_perm = nan(nPermSearch, 1);
                for ci_perm = 1:nPermSearch
                    w_perm = W_perm(:, ci_perm);
                    topo_perm = covStim_perm_reg * w_perm;
                    r_perm = corr(topo_perm, sim_template, 'rows', 'complete');
                    if ~isnan(r_perm) && r_perm < 0
                        w_perm = -w_perm;
                        topo_perm = -topo_perm;
                        r_perm = -r_perm;
                    end
                    occ_perm = mean(abs(topo_perm(occ_idx)));
                    if ~isempty(front_idx)
                        front_perm = mean(abs(topo_perm(front_idx)));
                        ratio_perm(ci_perm) = occ_perm / max(front_perm, eps);
                    else
                        ratio_perm(ci_perm) = Inf;
                    end
                    gamma_perm(ci_perm) = log(max((w_perm' * covStim_perm_reg * w_perm), eps) / ...
                                             max((w_perm' * covBase_perm_reg * w_perm), eps));
                    corr_perm(ci_perm) = r_perm;
                end
                finite_perm = isfinite(eval_perm) & isfinite(corr_perm) & ...
                    isfinite(ratio_perm) & isfinite(gamma_perm);
                hard_perm = finite_perm & ...
                    (eval_perm >= min_eigval_hard) & ...
                    (corr_perm >= min_corr_hard) & ...
                    (ratio_perm >= min_occfront_ratio) & ...
                    (gamma_perm >= min_gamma_hard);
                if any(hard_perm)
                    perm_null_max(pi) = max(eval_perm(hard_perm));
                end
            end
            valid_null = perm_null_max(isfinite(perm_null_max));
            if numel(valid_null) >= max(20, round(0.10 * n_perm))
                perm_threshold = prctile(valid_null, 100*(1-perm_alpha));
            else
                warning('Too few valid permutation draws after hard gating for subject %s; relaxing null to unconstrained max eigen.', subjects{subj});
                valid_null_unconstrained = perm_null_max_unconstrained(isfinite(perm_null_max_unconstrained));
                if isempty(valid_null_unconstrained)
                    perm_threshold = Inf;
                else
                    perm_threshold = prctile(valid_null_unconstrained, 100*(1-perm_alpha));
                end
            end
            perm_sig_mask = hard_eligible & (evals_sorted(1:nSearch) > perm_threshold);
        else
            warning('Permutation testing skipped for subject %s due to insufficient trials.', subjects{subj});
            perm_sig_mask = false(nSearch, 1);
        end
    end
    all_perm_thresholds(subj) = perm_threshold;
    candidate_table.perm_sig = perm_sig_mask;

    candidate_pool_mask = hard_eligible & isfinite(searchScores);
    if strcmpi(comp_select_mode, 'perm_cv')
        candidate_pool_mask = candidate_pool_mask & perm_sig_mask;
    end
    candidate_pool_idx = find(candidate_pool_mask);
    if isempty(candidate_pool_idx)
        warning(['No components survived hard+permutation gating for subject %s. ', ...
                 'Falling back to single best ranked component.'], subjects{subj});
        candidate_pool_idx = bestIdx;
    else
        [~, pool_ord] = sort(evals_sorted(candidate_pool_idx), 'descend');
        candidate_pool_idx = candidate_pool_idx(pool_ord);
    end

    k_max = numel(candidate_pool_idx);
    if k_max < 1
        k_max = 1;
        candidate_pool_idx = bestIdx;
    end

    cv_curve_mean = nan(1, k_max);
    cv_comp_utility = nan(1, k_max);
    cv_fold_scores = nan(cv_repeats * cv_folds, k_max);
    cv_fold_sep = nan(cv_repeats * cv_folds, k_max);
    nCvRows = 0;
    if strcmpi(comp_select_mode, 'perm_cv') && nTrials_total >= cv_folds && ...
            nTrials_total >= (cv_min_train_trials + cv_min_test_trials) && k_max >= 1
        for rep = 1:cv_repeats
            trl_perm = randperm(nTrials_total);
            fold_id = zeros(1, nTrials_total);
            fold_edges = round(linspace(0, nTrials_total, cv_folds+1));
            for fi_cv = 1:cv_folds
                fold_range = (fold_edges(fi_cv)+1):fold_edges(fi_cv+1);
                if isempty(fold_range), continue; end
                fold_id(trl_perm(fold_range)) = fi_cv;
            end

            for fi_cv = 1:cv_folds
                test_mask = (fold_id == fi_cv);
                train_mask = ~test_mask;
                nTest = sum(test_mask);
                nTrain = sum(train_mask);
                if nTest < cv_min_test_trials || nTrain < cv_min_train_trials
                    continue;
                end

                covStim_train = zeros(nChans);
                covBase_train = zeros(nChans);
                covStim_test = zeros(nChans);
                covBase_test = zeros(nChans);
                train_idx = find(train_mask);
                test_idx = find(test_mask);
                for ti = train_idx
                    covStim_train = covStim_train + stim_cov_trials{ti};
                    covBase_train = covBase_train + base_cov_trials{ti};
                end
                for ti = test_idx
                    covStim_test = covStim_test + stim_cov_trials{ti};
                    covBase_test = covBase_test + base_cov_trials{ti};
                end
                covStim_train = covStim_train / numel(train_idx);
                covBase_train = covBase_train / numel(train_idx);
                covStim_test = covStim_test / numel(test_idx);
                covBase_test = covBase_test / numel(test_idx);

                covStim_train_reg = (1-lambda)*covStim_train + lambda*mean(diag(covStim_train))*eye(nChans);
                covBase_train_reg = (1-lambda)*covBase_train + lambda*mean(diag(covBase_train))*eye(nChans);
                covStim_test_reg = (1-lambda)*covStim_test + lambda*mean(diag(covStim_test))*eye(nChans);
                covBase_test_reg = (1-lambda)*covBase_test + lambda*mean(diag(covBase_test))*eye(nChans);

                [W_cv, D_cv] = eig(covStim_train_reg, covBase_train_reg);
                [~, cv_sort] = sort(real(diag(D_cv)), 'descend');
                W_cv = W_cv(:, cv_sort);
                nCvComp = min(nSearch, size(W_cv, 2));
                if nCvComp < 1
                    continue;
                end

                nCvRows = nCvRows + 1;
                for ki = 1:k_max
                    target_idx = candidate_pool_idx(ki);
                    target_topo = searchTopos(:, target_idx);
                    best_map_abs = -Inf;
                    best_map_idx = 1;
                    best_map_sign = 1;
                    for ci_cv = 1:nCvComp
                        topo_cv = covStim_train_reg * W_cv(:, ci_cv);
                        r_map = corr(topo_cv, target_topo, 'rows', 'complete');
                        if isfinite(r_map) && abs(r_map) > best_map_abs
                            best_map_abs = abs(r_map);
                            best_map_idx = ci_cv;
                            best_map_sign = sign(r_map);
                        end
                    end
                    w_cv = W_cv(:, best_map_idx);
                    if best_map_sign < 0
                        w_cv = -w_cv;
                    end
                    sep_cv = log(max(w_cv' * covStim_test_reg * w_cv, eps) / ...
                                 max(w_cv' * covBase_test_reg * w_cv, eps));
                    cv_fold_sep(nCvRows, ki) = sep_cv;
                    if ki == 1
                        cv_fold_scores(nCvRows, ki) = sep_cv;
                    else
                        cv_fold_scores(nCvRows, ki) = cv_fold_scores(nCvRows, ki-1) + sep_cv;
                    end
                end
            end
        end
    end

    if nCvRows > 0
        cv_curve_mean = nanmean(cv_fold_scores(1:nCvRows, :), 1);
        cv_gain_curve = [cv_curve_mean(1), diff(cv_curve_mean)];
        valid_gain_idx = find(cv_gain_curve > cv_gain_min);
        if isempty(valid_gain_idx)
            cv_k_opt = 1;
        else
            cv_k_opt = max(valid_gain_idx);
        end
        cv_k_opt = max(1, min(cv_k_opt, k_max));
        cv_comp_utility = nanmean(cv_fold_sep(1:nCvRows, 1:cv_k_opt), 1);
    else
        cv_k_opt = 1;
        cv_curve_mean = nan(1, k_max);
    end

    selected_idx = candidate_pool_idx(1:min(cv_k_opt, numel(candidate_pool_idx)));
    if isempty(selected_idx)
        selected_idx = bestIdx;
    end

    selected_weights = [];
    if strcmpi(projection_mode, 'weighted_power')
        if ~isempty(cv_comp_utility) && numel(cv_comp_utility) >= numel(selected_idx)
            selected_weights = cv_comp_utility(1:numel(selected_idx));
        else
            selected_weights = evals_sorted(selected_idx)';
        end
        selected_weights(~isfinite(selected_weights) | selected_weights < 0) = 0;
        if sum(selected_weights) <= 0
            selected_weights = evals_sorted(selected_idx)';
            selected_weights(~isfinite(selected_weights) | selected_weights < 0) = 0;
        end
        if sum(selected_weights) <= 0
            selected_weights = ones(1, numel(selected_idx));
        end
        selected_weights = selected_weights / sum(selected_weights);
    end

    preperm_idx = find(hard_eligible & isfinite(searchScores));
    if isempty(preperm_idx)
        warning(['No hard-eligible pre-Perm+CV components for subject %s. ', ...
                 'Falling back to best single component for pre-Perm+CV combined benchmark.'], subjects{subj});
        preperm_idx = bestIdx;
    else
        [~, pre_ord] = sort(evals_sorted(preperm_idx), 'descend');
        preperm_idx = preperm_idx(pre_ord);
    end
    preperm_weights = searchScores(preperm_idx)';
    preperm_weights(~isfinite(preperm_weights) | preperm_weights < 0) = 0;
    if sum(preperm_weights) <= 0
        preperm_weights = evals_sorted(preperm_idx)';
        preperm_weights(~isfinite(preperm_weights) | preperm_weights < 0) = 0;
    end
    if sum(preperm_weights) <= 0
        preperm_weights = ones(1, numel(preperm_idx));
    end
    preperm_weights = preperm_weights / sum(preperm_weights);

    bestIdx = selected_idx(1);
    bestScore = searchScores(bestIdx);
    bestCorr = searchCorrs(bestIdx);
    bestOcc = searchOccStrength(bestIdx);
    bestFront = searchFrontStrength(bestIdx);
    bestRatio = searchOccFrontRatio(bestIdx);
    bestGamma = searchGammaEvidence(bestIdx);
    bestLeak = searchFrontLeak(bestIdx);

    topComp = searchFilters(:, bestIdx);
    if numel(selected_idx) > 1
        topo_temp = searchTopos(:, selected_idx) * selected_weights(:);
    else
        topo_temp = covStim_full_reg * topComp;
    end

    finite_scores = find(isfinite(searchScores));
    if isempty(finite_scores)
        warning('No finite eligible component scores for subject %s; using unconstrained ordering for display.', subjects{subj});
        [~, topDispOrder] = sort(comp_score, 'descend');
    else
        [~, topDispOrder] = sort(searchScores, 'descend');
    end
    nDispTopo = min(topo_display_n, sum(isfinite(searchScores)));
    if nDispTopo < 1
        nDispTopo = 1;
    end
    dispTopoIdx = topDispOrder(1:nDispTopo);
    topoDispTopos = searchTopos(:, dispTopoIdx);
    topoDispCorrs = searchCorrs(dispTopoIdx);
    topoDispScores = searchScores(dispTopoIdx);

    nStore = min(5, nSearch);
    storeCompIdx = topDispOrder(1:nStore);
    storeCorrs = searchCorrs(storeCompIdx);
    storeEvals = evals_sorted(storeCompIdx);
    storeTopos = searchTopos(:, storeCompIdx);

    all_topos{subj}       = topo_temp;
    all_topo_labels{subj} = dataEEG_c25.label;
    all_eigenvalues(subj) = evals_sorted(bestIdx);
    all_selected_comp_idx(subj)  = bestIdx;
    all_selected_comp_corr(subj) = bestCorr;
    all_selected_comp_eval(subj) = evals_sorted(bestIdx);
    all_top5_corrs(1:nStore, subj) = storeCorrs;
    all_top5_evals(1:nStore, subj) = storeEvals;
    all_top5_topos{subj} = storeTopos;
    all_simulated_templates{subj} = sim_template;
    all_selected_comp_indices_multi{subj} = selected_idx;
    all_selected_comp_weights{subj} = selected_weights(:)';
    all_preperm_comp_indices_multi{subj} = preperm_idx;
    all_preperm_comp_weights{subj} = preperm_weights(:)';
    all_cv_curves{subj} = struct( ...
        'mean_curve', cv_curve_mean, ...
        'fold_scores', cv_fold_scores(1:max(nCvRows, 1), :), ...
        'fold_sep', cv_fold_sep(1:max(nCvRows, 1), :), ...
        'k_opt', cv_k_opt);
    all_component_selection_stats{subj} = struct( ...
        'mode', comp_select_mode, ...
        'projection_mode', projection_mode, ...
        'perm_threshold', perm_threshold, ...
        'perm_alpha', perm_alpha, ...
        'n_perm', n_perm, ...
        'n_cv_rows', nCvRows, ...
        'cv_folds', cv_folds, ...
        'cv_repeats', cv_repeats, ...
        'preperm_idx', preperm_idx, ...
        'preperm_weights', preperm_weights, ...
        'selected_idx', selected_idx, ...
        'selected_weights', selected_weights, ...
        'candidate_pool_idx', candidate_pool_idx, ...
        'candidate_table', candidate_table, ...
        'best_idx', bestIdx, ...
        'best_score', bestScore, ...
        'best_corr', bestCorr, ...
        'best_ratio', bestRatio, ...
        'best_gamma', bestGamma, ...
        'best_front', bestFront, ...
        'best_occ', bestOcc, ...
        'best_leak', bestLeak);

    %% Phase 1b: Component/template sanity figures
    cfg_topo = [];
    cfg_topo.layout    = headmodel.layANThead;
    cfg_topo.comment   = 'no';
    cfg_topo.marker    = 'off';
    cfg_topo.style     = 'straight';
    cfg_topo.gridscale = 300;
    cfg_topo.zlim      = 'maxabs';
    cfg_topo.colormap  = '*RdBu';
    cfg_topo.figure    = 'gcf';

    fig_sel = figure('Position', [0 0 1512 982], 'Color', 'w');
    title_str = sprintf(['Subject %s: RAW GED Top %d (best C%d, ', ...
        'score=%.3f, r=%.3f, occ/front=%.3f, gamma=%.3f, leak=%.3f)'], ...
        subjects{subj}, nDispTopo, bestIdx, bestScore, bestCorr, bestRatio, bestGamma, bestLeak);
    annotation(fig_sel, 'textbox', [0.01 0.965 0.98 0.03], ...
        'String', title_str, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 14, 'FontWeight', 'bold');
    nColsTopo = 5;
    nRowsTopo = ceil(nDispTopo / nColsTopo);
    for ci = 1:nDispTopo
        subplot(nRowsTopo, nColsTopo, ci);
        topo_data = [];
        topo_data.label  = all_topo_labels{subj};
        topo_plot_ci = topoDispTopos(:, ci);
        if viz_suppress_nonocc_outliers && ~isempty(nonocc_idx)
            post_abs = abs(topo_plot_ci(post_idx));
            post_abs = post_abs(isfinite(post_abs));
            if isempty(post_abs)
                post_abs = abs(topo_plot_ci(isfinite(topo_plot_ci)));
            end
            if isempty(post_abs)
                amp_thr_ci = 1;
            else
                amp_thr_ci = prctile(post_abs, viz_topo_prctile);
                if ~isfinite(amp_thr_ci) || amp_thr_ci <= 0
                    amp_thr_ci = max(post_abs);
                end
                if ~isfinite(amp_thr_ci) || amp_thr_ci <= 0
                    amp_thr_ci = 1;
                end
            end
            nonocc_out = false(nChans, 1);
            nonocc_out(nonocc_idx) = abs(topo_plot_ci(nonocc_idx)) > (viz_nonocc_outlier_mult * amp_thr_ci);
            if any(nonocc_out)
                valid_interp = find(~nonocc_out & isfinite(topo_plot_ci) & has_pos);
                out_idx = find(nonocc_out);
                for oi = out_idx(:)'
                    if has_pos(oi) && numel(valid_interp) >= 3
                        d = sqrt(sum((chan_pos(valid_interp, :) - chan_pos(oi, :)).^2, 2));
                        [d_sorted, d_ord] = sort(d, 'ascend');
                        k_use = min(viz_interp_k, numel(d_sorted));
                        nbr_idx = valid_interp(d_ord(1:k_use));
                        w = 1 ./ max(d_sorted(1:k_use), 1e-6);
                        topo_plot_ci(oi) = sum(w .* topo_plot_ci(nbr_idx)) / sum(w);
                    else
                        topo_plot_ci(oi) = sign(topo_plot_ci(oi)) * amp_thr_ci;
                    end
                end
            end
        end
        topo_data.avg    = topo_plot_ci;
        topo_data.dimord = 'chan';
        topo_abs_ci = abs(topo_data.avg(post_idx));
        topo_abs_ci = topo_abs_ci(isfinite(topo_abs_ci));
        if isempty(topo_abs_ci)
            topo_abs_ci = abs(topo_data.avg(:));
            topo_abs_ci = topo_abs_ci(isfinite(topo_abs_ci));
        end
        if isempty(topo_abs_ci)
            topo_clim_ci = 1;
        else
            topo_clim_ci = prctile(topo_abs_ci, viz_topo_prctile);
            if ~isfinite(topo_clim_ci) || topo_clim_ci <= 0
                topo_clim_ci = max(topo_abs_ci);
            end
            if ~isfinite(topo_clim_ci) || topo_clim_ci <= 0
                topo_clim_ci = 1;
            end
        end
        cfg_topo_ci = cfg_topo;
        cfg_topo_ci.zlim = [-topo_clim_ci topo_clim_ci];
        try
            ft_topoplotER(cfg_topo_ci, topo_data);
            colorbar;
        catch
            imagesc(topo_plot_ci); caxis([-topo_clim_ci topo_clim_ci]); colorbar;
        end
        comp_rank = dispTopoIdx(ci);
        title(sprintf('C%d: \\lambda=%.2f, score=%.2f, r=%.2f, ratio=%.2f, g=%.2f', ...
            comp_rank, evals_sorted(comp_rank), topoDispScores(ci), topoDispCorrs(ci), ...
            searchOccFrontRatio(comp_rank), searchGammaEvidence(comp_rank)), ...
            'FontSize', 7);
    end
    saveas(fig_sel, fullfile(comp_sel_save_dir, sprintf('GCP_eeg_GED_component_selection_subj%s_prePermCV.png', subjects{subj})));

    % Post-selection topoplots: components surviving permutation + CV cap.
    selTopoIdx = selected_idx(:)';
    nSelTopo = numel(selTopoIdx);
    nColsSel = 5;
    nRowsSel = ceil(max(nSelTopo, 1) / nColsSel);
    fig_post = figure('Position', [0 0 1512 982], 'Color', 'w');
    title_post = sprintf(['Subject %s: POST perm+CV (n=%d, k*=%d, perm thr=%.3f, ', ...
        'proj=%s)'], subjects{subj}, nSelTopo, cv_k_opt, perm_threshold, projection_mode);
    annotation(fig_post, 'textbox', [0.01 0.965 0.98 0.03], ...
        'String', title_post, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 14, 'FontWeight', 'bold');
    if nSelTopo == 0
        subplot(1, 1, 1);
        axis off;
        text(0.5, 0.5, 'No selected components', 'HorizontalAlignment', 'center');
    else
        for si = 1:nSelTopo
            subplot(nRowsSel, nColsSel, si);
            comp_rank = selTopoIdx(si);
            topo_data = [];
            topo_data.label = all_topo_labels{subj};
            topo_data.avg = searchTopos(:, comp_rank);
            topo_data.dimord = 'chan';
            topo_abs_ci = abs(topo_data.avg(post_idx));
            topo_abs_ci = topo_abs_ci(isfinite(topo_abs_ci));
            if isempty(topo_abs_ci)
                topo_abs_ci = abs(topo_data.avg(:));
                topo_abs_ci = topo_abs_ci(isfinite(topo_abs_ci));
            end
            if isempty(topo_abs_ci)
                topo_clim_ci = 1;
            else
                topo_clim_ci = prctile(topo_abs_ci, viz_topo_prctile);
                if ~isfinite(topo_clim_ci) || topo_clim_ci <= 0
                    topo_clim_ci = max(topo_abs_ci);
                end
                if ~isfinite(topo_clim_ci) || topo_clim_ci <= 0
                    topo_clim_ci = 1;
                end
            end
            cfg_topo_ci = cfg_topo;
            cfg_topo_ci.zlim = [-topo_clim_ci topo_clim_ci];
            try
                ft_topoplotER(cfg_topo_ci, topo_data);
                colorbar;
            catch
                imagesc(topo_data.avg); caxis([-topo_clim_ci topo_clim_ci]); colorbar;
            end
            w_show = NaN;
            if ~isempty(selected_weights) && numel(selected_weights) >= si
                w_show = selected_weights(si);
            end
            title(sprintf('C%d: \\lambda=%.2f, w=%.3f, r=%.2f, ratio=%.2f', ...
                comp_rank, evals_sorted(comp_rank), w_show, searchCorrs(comp_rank), ...
                searchOccFrontRatio(comp_rank)), 'FontSize', 8);
        end
    end
    saveas(fig_post, fullfile(comp_sel_save_dir, sprintf('GCP_eeg_GED_component_selection_subj%s_postPermCV.png', subjects{subj})));

    %% ================================================================
    %  PHASE 2: Per condition — trial-level narrowband scanning
    %  ================================================================
    W_sel = searchFilters(:, selected_idx);
    w_sel = all_selected_comp_weights{subj};
    W_pre = searchFilters(:, preperm_idx);
    w_pre = all_preperm_comp_weights{subj};
    if isempty(W_sel)
        warning('No selected filters for subject %s. Falling back to best single component.', subjects{subj});
        W_sel = zeros(0, 0);
    end
    if isempty(W_pre)
        warning('No pre-Perm+CV filters for subject %s. Falling back to best single component.', subjects{subj});
        W_pre = searchFilters(:, bestIdx);
    end
    if isempty(W_top)
        W_top = zeros(0, 1);
    end
    if isempty(w_sel) && ~isempty(W_sel)
        w_sel = ones(1, size(W_sel, 2)) / size(W_sel, 2);
    end
    if sum(w_sel) <= 0 && ~isempty(W_sel)
        w_sel = ones(1, size(W_sel, 2)) / size(W_sel, 2);
    end
    if ~isempty(W_sel)
        w_sel = w_sel(:)' / sum(w_sel);
    end
    if isempty(w_pre) && ~isempty(W_pre)
        w_pre = ones(1, size(W_pre, 2)) / size(W_pre, 2);
    end
    if sum(w_pre) <= 0 && ~isempty(W_pre)
        w_pre = ones(1, size(W_pre, 2)) / size(W_pre, 2);
    end
    if ~isempty(W_pre)
        w_pre = w_pre(:)' / sum(w_pre);
    end

    for cond = 1:4

        dat = dat_per_cond{cond};
        if isempty(dat), continue; end

        nTrl = length(dat.trial);
        powratio_methods = nan(nBenchmarkMethods, nTrl, nFreqs);

        for fi = 1:nFreqs
            clc
            fprintf('Subject    %s (%d/%d)\nCondition  %d/4\nFrequency  %d/%d\n', ...
                subjects{subj}, subj, nSubj, cond, fi, nFreqs);
            cf = scan_freqs(fi);
            bpfreq = [max(cf - scan_width/2, 1), cf + scan_width/2];

            cfg_filt = [];
            cfg_filt.bpfilter   = 'yes';
            cfg_filt.bpfreq     = bpfreq;
            cfg_filt.bpfilttype = 'fir';
            cfg_filt.bpfiltord  = round(3 * fsample / bpfreq(1));
            dat_nb = ft_preprocessing(cfg_filt, dat);

            cfg_t = [];
            cfg_t.latency = baseline_window;
            dat_base_nb = ft_selectdata(cfg_t, dat_nb);

            cfg_t.latency = stimulus_window;
            dat_stim_nb = ft_selectdata(cfg_t, dat_nb);

            for trl = 1:nTrl
                x_stim = double(dat_stim_nb.trial{trl});
                x_base = double(dat_base_nb.trial{trl});

                % Method 1: raw channel-space weighted power ratio.
                pow_stim_chan = mean(x_stim.^2, 2);
                pow_base_chan = mean(x_base.^2, 2);
                raw_pow_stim = sum(raw_w(:) .* pow_stim_chan(:));
                raw_pow_base = sum(raw_w(:) .* pow_base_chan(:));
                if isfinite(raw_pow_stim) && isfinite(raw_pow_base) && raw_pow_base > 0
                    powratio_methods(1, trl, fi) = raw_pow_stim / raw_pow_base;
                end

                % Method 2: top GED component by eigenvalue.
                if ~isempty(W_top)
                    comp_stim_top = W_top' * x_stim;
                    comp_base_top = W_top' * x_base;
                    pow_stim_top = mean(comp_stim_top.^2);
                    pow_base_top = mean(comp_base_top.^2);
                    if isfinite(pow_stim_top) && isfinite(pow_base_top) && pow_base_top > 0
                        powratio_methods(2, trl, fi) = pow_stim_top / pow_base_top;
                    end
                end

                % Method 3: pre-Perm+CV weighted GED combination.
                if ~isempty(W_pre)
                    comp_stim = W_pre' * x_stim;
                    comp_base = W_pre' * x_base;
                    pow_stim_vec = mean(comp_stim.^2, 2);
                    pow_base_vec = mean(comp_base.^2, 2);
                    valid_comp = isfinite(pow_stim_vec) & isfinite(pow_base_vec) & (pow_base_vec > 0);
                    if any(valid_comp)
                        ratio_vec = pow_stim_vec(valid_comp) ./ pow_base_vec(valid_comp);
                        w_use = w_pre(valid_comp);
                        if sum(w_use) <= 0
                            w_use = ones(1, numel(ratio_vec)) / numel(ratio_vec);
                        else
                            w_use = w_use / sum(w_use);
                        end
                        powratio_methods(3, trl, fi) = sum(ratio_vec(:)' .* w_use);
                    end
                end

                % Method 4: post-Perm+CV weighted GED combination.
                if ~isempty(W_sel)
                    comp_stim = W_sel' * x_stim;
                    comp_base = W_sel' * x_base;
                    pow_stim_vec = mean(comp_stim.^2, 2);
                    pow_base_vec = mean(comp_base.^2, 2);
                    valid_comp = isfinite(pow_stim_vec) & isfinite(pow_base_vec) & (pow_base_vec > 0);
                    if any(valid_comp)
                        ratio_vec = pow_stim_vec(valid_comp) ./ pow_base_vec(valid_comp);
                        w_use = w_sel(valid_comp);
                        if sum(w_use) <= 0
                            w_use = ones(1, numel(ratio_vec)) / numel(ratio_vec);
                        else
                            w_use = w_use / sum(w_use);
                        end
                        powratio_methods(4, trl, fi) = sum(ratio_vec(:)' .* w_use);
                    end
                end
            end
        end

        for mi = 1:nBenchmarkMethods
            pr_m = squeeze(powratio_methods(mi, :, :));
            all_trial_powratio_bench{mi, cond, subj} = pr_m;
            if ~isempty(pr_m)
                dt_m = nan(size(pr_m));
                for trl = 1:size(pr_m, 1)
                    p = polyfit(scan_freqs, pr_m(trl,:), poly_order);
                    dt_m(trl,:) = pr_m(trl,:) - polyval(p, scan_freqs);
                end
                all_trial_powratio_dt_bench{mi, cond, subj} = dt_m;
            end
        end

        % Keep legacy pipeline outputs based on post-Perm+CV weighted GED.
        powratio_trials = squeeze(powratio_methods(4, :, :));
        all_trial_powratio{cond, subj} = powratio_trials;

        %% Per-trial peak detection
        trl_peaks_single = nan(nTrl, 1);
        trl_peaks_low    = nan(nTrl, 1);
        trl_peaks_high   = nan(nTrl, 1);
        trl_centroid     = nan(nTrl, 1);

        for trl = 1:nTrl
            pr = powratio_trials(trl, :);
            if all(isnan(pr)), continue; end

            p = polyfit(scan_freqs, pr, poly_order);
            pr_dt = pr - polyval(p, scan_freqs);
            pr_dt_smooth = movmean(pr_dt, 5);

            % Stripe-center metric: spectral centroid of positive detrended mass in 40-80 Hz.
            pr_dt_band = pr_dt_smooth(centroid_band_mask);
            freq_band = scan_freqs(centroid_band_mask);
            w_pos = max(pr_dt_band, 0);
            pos_mass = sum(w_pos);
            total_abs_mass = sum(abs(pr_dt_band));
            if pos_mass > 0 && max(pr_dt_band) >= centroid_min_peak && ...
                    (pos_mass / max(total_abs_mass, eps)) >= centroid_posfrac_min
                trl_centroid(trl) = sum(freq_band .* w_pos) / pos_mass;
            end

            % Single-peak: tallest prominent peak
            [pks, locs] = findpeaks(pr_dt_smooth, scan_freqs, ...
                'MinPeakProminence', max(pr_dt_smooth) * 0.15, ...
                'MinPeakDistance', 5);

            if ~isempty(pks)
                [~, best_pk] = max(pks);
                trl_peaks_single(trl) = locs(best_pk);
            end

            % Dual-peak: hard 50 Hz boundary (primary method)
            [pks_all, locs_all] = findpeaks(pr_dt_smooth, scan_freqs, ...
                'MinPeakDistance', 5);
            pos_mask = pks_all > 0;
            pks_pos  = pks_all(pos_mask);
            locs_pos = locs_all(pos_mask);

            in_lo = locs_pos >= 30 & locs_pos <= 49;
            in_hi = locs_pos >= 50 & locs_pos <= 90;
            if any(in_lo)
                [~, bi] = max(pks_pos(in_lo));
                tmp = locs_pos(in_lo);
                trl_peaks_low(trl) = tmp(bi);
            end
            if any(in_hi)
                [~, bi] = max(pks_pos(in_hi));
                tmp = locs_pos(in_hi);
                trl_peaks_high(trl) = tmp(bi);
            end
        end

        all_trial_peaks_single{cond, subj} = trl_peaks_single;
        all_trial_peaks_low{cond, subj}    = trl_peaks_low;
        all_trial_peaks_high{cond, subj}   = trl_peaks_high;
        all_trial_centroid{cond, subj}     = trl_centroid;

        valid_s = ~isnan(trl_peaks_single);
        all_trial_mean_single(cond, subj)   = mean(trl_peaks_single(valid_s));
        all_trial_median_single(cond, subj) = median(trl_peaks_single(valid_s));
        all_trial_detrate_single(cond, subj) = sum(valid_s) / nTrl;

        valid_lo = ~isnan(trl_peaks_low);
        all_trial_mean_low(cond, subj)   = mean(trl_peaks_low(valid_lo));
        all_trial_median_low(cond, subj) = median(trl_peaks_low(valid_lo));
        all_trial_detrate_low(cond, subj) = sum(valid_lo) / nTrl;

        valid_hi = ~isnan(trl_peaks_high);
        all_trial_mean_high(cond, subj)   = mean(trl_peaks_high(valid_hi));
        all_trial_median_high(cond, subj) = median(trl_peaks_high(valid_hi));
        all_trial_detrate_high(cond, subj) = sum(valid_hi) / nTrl;
        valid_gap = valid_lo & valid_hi;
        if any(valid_gap)
            all_trial_median_gap(cond, subj) = median(trl_peaks_high(valid_gap) - trl_peaks_low(valid_gap));
        end

        valid_c = ~isnan(trl_centroid);
        all_trial_mean_centroid(cond, subj)   = mean(trl_centroid(valid_c));
        all_trial_median_centroid(cond, subj) = median(trl_centroid(valid_c));
        all_trial_detrate_centroid(cond, subj) = sum(valid_c) / nTrl;

    end % condition loop

    %% ================================================================
    %  PER-SUBJECT FIGURE (4 rows x 4 columns)
    %  ================================================================
    close all
    fig = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('Trial-Level GED: Subject %s', subjects{subj}), ...
        'FontSize', 18, 'FontWeight', 'bold');

    cmap_div = interp1([0 0.5 1], ...
        [0.17 0.27 0.53; 0.97 0.97 0.97; 0.70 0.09 0.17], linspace(0,1,256));

    % Pre-compute detrended matrices
    pr_dt_mats = cell(1, 4);
    for cond = 1:4
        pr_mat = all_trial_powratio{cond, subj};
        if ~isempty(pr_mat)
            nTrl = size(pr_mat, 1);
            dt = nan(size(pr_mat));
            for trl = 1:nTrl
                p = polyfit(scan_freqs, pr_mat(trl,:), poly_order);
                dt(trl,:) = pr_mat(trl,:) - polyval(p, scan_freqs);
            end
            pr_dt_mats{cond} = dt;
        end
    end
    row1_clim = ones(1, 4);
    for cond = 1:4
        if ~isempty(pr_dt_mats{cond})
            cond_vals = abs(pr_dt_mats{cond}(:));
            cond_vals = cond_vals(isfinite(cond_vals));
            if ~isempty(cond_vals)
                row1_clim(cond) = prctile(cond_vals, 95);
                if ~isfinite(row1_clim(cond)) || row1_clim(cond) <= 0
                    row1_clim(cond) = max(cond_vals);
                end
            end
        end
        if ~isfinite(row1_clim(cond)) || row1_clim(cond) <= 0
            row1_clim(cond) = 1;
        end
    end

    % --- Row 1: Heatmap of trial-level detrended power-ratio spectra ---
    for cond = 1:4
        subplot(4, 4, cond);
        if ~isempty(pr_dt_mats{cond})
            hold on;
            imagesc(scan_freqs, 1:size(pr_dt_mats{cond},1), pr_dt_mats{cond});
            colormap(gca, cmap_div);
            caxis([-row1_clim(cond) row1_clim(cond)]);
            cb = colorbar; cb.FontSize = 8;
            ctd = all_trial_centroid{cond, subj};
            if ~isempty(ctd)
                valid_ctd = ~isnan(ctd);
                tr_idx = find(valid_ctd);
                if ~isempty(tr_idx)
                    plot(ctd(valid_ctd), tr_idx, '.', 'Color', [0 0 0], 'MarkerSize', 7);
                end
                if numel(tr_idx) >= 2
                    plot(ctd(valid_ctd), tr_idx, '-', 'Color', [0 0 0], 'LineWidth', 0.6);
                end
            end
            xlabel('Freq [Hz]'); ylabel('Trial');
            set(gca, 'YDir', 'normal');
        end
        title(sprintf('%s Detrended', condLabels{cond}), 'FontSize', 11);
        set(gca, 'FontSize', 10); xlim([30 90]); box on;
    end

    % --- Row 2: Mean trial-level spectrum with single-peak markers ---
    for cond = 1:4
        subplot(4, 4, 4 + cond); hold on;
        if ~isempty(pr_dt_mats{cond})
            mu_dt = nanmean(pr_dt_mats{cond}, 1);
            nTrl = size(pr_dt_mats{cond}, 1);
            sem_dt = nanstd(pr_dt_mats{cond}, [], 1) / sqrt(nTrl);
            row2_abs = max(abs([mu_dt - sem_dt, mu_dt + sem_dt]), [], 'omitnan');
            if ~isfinite(row2_abs) || row2_abs <= 0
                row2_abs = 1;
            end
            faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
            patch([scan_freqs, fliplr(scan_freqs)], ...
                [mu_dt - sem_dt, fliplr(mu_dt + sem_dt)], ...
                colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            plot(scan_freqs, mu_dt, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);
            yline(0, 'k-', 'LineWidth', 0.5);
            mn_pf = all_trial_mean_single(cond, subj);
            md_pf = all_trial_median_single(cond, subj);
            if ~isnan(mn_pf)
                xline(mn_pf, '--', 'LineWidth', 2, 'Color', colors(cond,:));
                text(mn_pf + 1, max(mu_dt) * 0.9, ...
                    sprintf('mn:%.0f', mn_pf), 'FontSize', 9, ...
                    'Color', colors(cond,:), 'FontWeight', 'bold');
            end
            if ~isnan(md_pf)
                xline(md_pf, ':', 'LineWidth', 2, 'Color', colors(cond,:));
                text(md_pf + 1, max(mu_dt) * 0.7, ...
                    sprintf('md:%.0f', md_pf), 'FontSize', 9, ...
                    'Color', colors(cond,:));
            end
            ylim([-row2_abs row2_abs]);
        end
        xlabel('Freq [Hz]'); ylabel('\Delta PR');
        det = all_trial_detrate_single(cond, subj);
        title(sprintf('%s Single (det=%.0f%%)', condLabels{cond}, det*100), 'FontSize', 10);
        set(gca, 'FontSize', 10); xlim([30 90]);  box on;
    end

    % --- Row 3: Mean trial-level spectrum with dual-peak markers (50 Hz) ---
    for cond = 1:4
        subplot(4, 4, 8 + cond); hold on;
        if ~isempty(pr_dt_mats{cond})
            mu_dt = nanmean(pr_dt_mats{cond}, 1);
            nTrl = size(pr_dt_mats{cond}, 1);
            sem_dt = nanstd(pr_dt_mats{cond}, [], 1) / sqrt(nTrl);
            row3_abs = max(abs([mu_dt - sem_dt, mu_dt + sem_dt]), [], 'omitnan');
            if ~isfinite(row3_abs) || row3_abs <= 0
                row3_abs = 1;
            end
            faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
            patch([scan_freqs, fliplr(scan_freqs)], ...
                [mu_dt - sem_dt, fliplr(mu_dt + sem_dt)], ...
                colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            plot(scan_freqs, mu_dt, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);
            yline(0, 'k-', 'LineWidth', 0.5);
            xline(50, 'k:', 'LineWidth', 1, 'Alpha', 0.5);
            pf_lo = all_trial_median_low(cond, subj);
            pf_hi = all_trial_median_high(cond, subj);
            if ~isnan(pf_lo)
                xline(pf_lo, '--', 'LineWidth', 2, 'Color', [0 0 0.7]);
                text(pf_lo + 1, max(mu_dt) * 0.85, ...
                    sprintf('L:%.0f', pf_lo), 'FontSize', 9, ...
                    'Color', [0 0 0.7], 'FontWeight', 'bold');
            end
            if ~isnan(pf_hi)
                xline(pf_hi, '--', 'LineWidth', 2, 'Color', [0.7 0 0]);
                text(pf_hi + 1, max(mu_dt) * 0.65, ...
                    sprintf('H:%.0f', pf_hi), 'FontSize', 9, ...
                    'Color', [0.7 0 0], 'FontWeight', 'bold');
            end
            ylim([-row3_abs row3_abs]);
        end
        xlabel('Freq [Hz]'); ylabel('\Delta PR');
        det_lo = all_trial_detrate_low(cond, subj);
        det_hi = all_trial_detrate_high(cond, subj);
        title(sprintf('%s Dual (L:%.0f%% H:%.0f%%)', condLabels{cond}, det_lo*100, det_hi*100), ...
            'FontSize', 10);
        set(gca, 'FontSize', 10); xlim([30 90]);  box on;
    end

    % --- Row 4: Topoplot + histogram ---
    % Use the same occipital-channel definition as component selection.
    occ_highlight = dataEEG_c25.label(occ_idx);

    cfg_topo = [];
    cfg_topo.layout    = headmodel.layANThead;
    cfg_topo.comment   = 'no';
    cfg_topo.marker    = 'off';
    cfg_topo.style     = 'straight';
    cfg_topo.gridscale = 300;
    cfg_topo.zlim      = 'maxabs';
    cfg_topo.colormap  = '*RdBu';
    cfg_topo.figure    = 'gcf';
    cfg_topo.highlight          = {'on'};
    cfg_topo.highlightchannel   = {occ_highlight};
    cfg_topo.highlightsymbol    = {'.'};
    cfg_topo.highlightsize      = {12};
    cfg_topo.highlightcolor     = {[0 0 0]};

    subplot(4, 4, 13);
    if ~isempty(all_topos{subj})
        topo_data = [];
        topo_data.label  = all_topo_labels{subj};
        topo_plot_common = all_topos{subj};
        if viz_suppress_nonocc_outliers && ~isempty(nonocc_idx)
            post_abs = abs(topo_plot_common(post_idx));
            post_abs = post_abs(isfinite(post_abs));
            if isempty(post_abs)
                post_abs = abs(topo_plot_common(isfinite(topo_plot_common)));
            end
            if isempty(post_abs)
                amp_thr_common = 1;
            else
                amp_thr_common = prctile(post_abs, viz_topo_prctile);
                if ~isfinite(amp_thr_common) || amp_thr_common <= 0
                    amp_thr_common = max(post_abs);
                end
                if ~isfinite(amp_thr_common) || amp_thr_common <= 0
                    amp_thr_common = 1;
                end
            end
            nonocc_out = false(nChans, 1);
            nonocc_out(nonocc_idx) = abs(topo_plot_common(nonocc_idx)) > (viz_nonocc_outlier_mult * amp_thr_common);
            if any(nonocc_out)
                valid_interp = find(~nonocc_out & isfinite(topo_plot_common) & has_pos);
                out_idx = find(nonocc_out);
                for oi = out_idx(:)'
                    if has_pos(oi) && numel(valid_interp) >= 3
                        d = sqrt(sum((chan_pos(valid_interp, :) - chan_pos(oi, :)).^2, 2));
                        [d_sorted, d_ord] = sort(d, 'ascend');
                        k_use = min(viz_interp_k, numel(d_sorted));
                        nbr_idx = valid_interp(d_ord(1:k_use));
                        w = 1 ./ max(d_sorted(1:k_use), 1e-6);
                        topo_plot_common(oi) = sum(w .* topo_plot_common(nbr_idx)) / sum(w);
                    else
                        topo_plot_common(oi) = sign(topo_plot_common(oi)) * amp_thr_common;
                    end
                end
            end
        end
        topo_data.avg    = topo_plot_common;
        topo_data.dimord = 'chan';
        topo_abs_common = abs(topo_plot_common(post_idx));
        topo_abs_common = topo_abs_common(isfinite(topo_abs_common));
        if isempty(topo_abs_common)
            topo_abs_common = abs(topo_plot_common(isfinite(topo_plot_common)));
        end
        if isempty(topo_abs_common)
            topo_clim_common = 1;
        else
            topo_clim_common = prctile(topo_abs_common, viz_topo_prctile);
            if ~isfinite(topo_clim_common) || topo_clim_common <= 0
                topo_clim_common = max(topo_abs_common);
            end
            if ~isfinite(topo_clim_common) || topo_clim_common <= 0
                topo_clim_common = 1;
            end
        end
        cfg_topo_common = cfg_topo;
        cfg_topo_common.zlim = [-topo_clim_common topo_clim_common];
        try
            ft_topoplotER(cfg_topo_common, topo_data);
            cb = colorbar; cb.FontSize = 9;
        catch
            imagesc(topo_data.avg); caxis([-topo_clim_common topo_clim_common]); colorbar;
        end
        n_sel_show = numel(all_selected_comp_indices_multi{subj});
        title(sprintf('Weighted GED (%d comps, \\lambda=%.2f)', n_sel_show, all_eigenvalues(subj)), 'FontSize', 11);
    end

    subplot(4, 4, [14 15 16]); hold on;
    edges = 30:2:90;
    hist_mat = zeros(4, length(edges)-1);
    for cond = 1:4
        tpk = all_trial_peaks_single{cond, subj};
        if ~isempty(tpk)
            tpk = tpk(~isnan(tpk));
            hist_mat(cond,:) = histcounts(tpk, edges);
        end
    end
    centers = edges(1:end-1) + diff(edges)/2;
    bh = bar(centers, hist_mat', 'stacked', 'EdgeColor', 'none', 'BarWidth', 1);
    for cond = 1:4
        bh(cond).FaceColor = colors(cond,:);
    end
    xlabel('Peak Frequency [Hz]'); ylabel('Trial Count');
    title('Trial-Level Single Peak Distribution', 'FontSize', 12);
    legend(bh, condLabels, 'FontSize', 10, 'Location', 'best');
    set(gca, 'FontSize', 11); xlim([30 90]);  box on;

    saveas(fig, fullfile(fig_save_dir, sprintf('GCP_eeg_GED_subj%s.png', subjects{subj})));
    toc
end % subject loop

%% ====================================================================
%  FOUR-WAY BENCHMARK METRICS (raw vs top-eig GED vs pre/post combined GED)
%  ====================================================================
for mi = 1:nBenchmarkMethods
    for subj = 1:nSubj
        cond_medians = nan(4, 1);
        for cond = 1:4
            pr_dt = all_trial_powratio_dt_bench{mi, cond, subj};
            if isempty(pr_dt)
                continue;
            end
            nTrl = size(pr_dt, 1);
            peak_freq = nan(nTrl, 1);
            peak_prom = nan(nTrl, 1);
            for trl = 1:nTrl
                y = movmean(pr_dt(trl, :), 5);
                if all(isnan(y))
                    continue;
                end
                [pks, locs, ~, p] = findpeaks(y, scan_freqs, ...
                    'MinPeakDistance', 5, 'MinPeakProminence', max(y) * 0.15);
                if ~isempty(pks)
                    [~, bi] = max(pks);
                    peak_freq(trl) = locs(bi);
                    peak_prom(trl) = p(bi);
                end
            end

            if compute_detectability
                benchmark_metric_detectability(mi, cond, subj) = mean(~isnan(peak_freq));
                benchmark_metric_prominence(mi, cond, subj) = mean(peak_prom(~isnan(peak_prom)));
            end
            if compute_reliability
                vf = peak_freq(~isnan(peak_freq));
                if numel(vf) >= 2 && mean(vf) ~= 0
                    benchmark_metric_reliability_trialcv(mi, cond, subj) = std(vf) / abs(mean(vf));
                end
            end
            cond_medians(cond) = median(peak_freq(~isnan(peak_freq)));
        end

        if compute_separation
            vx = ~isnan(cond_medians);
            if sum(vx) >= 2
                p = polyfit(find(vx), cond_medians(vx)', 1);
                benchmark_metric_separation_slope(mi, subj) = p(1);
            end
            if ~isnan(cond_medians(1)) && ~isnan(cond_medians(4))
                benchmark_metric_separation_delta(mi, subj) = cond_medians(4) - cond_medians(1);
            end
        end
    end
end

for mi = 1:nBenchmarkMethods
    for cond = 1:4
        v = squeeze(benchmark_metric_reliability_trialcv(mi, cond, :));
        benchmark_metric_reliability_subjspread(mi, cond) = std(v(~isnan(v)));
    end
end

% Shared y-limits for benchmark lower-row panels (consistent across subjects).
single_vals = all_trial_median_single(isfinite(all_trial_median_single));
if isempty(single_vals)
    bench_single_ylim = [30 90];
else
    single_q = prctile(single_vals, [2 98]);
    bench_single_ylim = [max(30, single_q(1) - 1), min(90, single_q(2) + 1)];
end

gap_vals = all_trial_median_gap(isfinite(all_trial_median_gap));
if isempty(gap_vals)
    bench_gap_ylim = [-10 10];
else
    gap_q = prctile(gap_vals, [2 98]);
    gap_pad = max(1, 0.1 * (gap_q(2) - gap_q(1)));
    bench_gap_ylim = [gap_q(1) - gap_pad, gap_q(2) + gap_pad];
end

prom_vals = benchmark_metric_prominence(isfinite(benchmark_metric_prominence));
if isempty(prom_vals)
    bench_prom_ylim = [0 1];
else
    bench_prom_ylim = [0, max(prom_vals) * 1.15];
end

rel_vals = benchmark_metric_reliability_trialcv(isfinite(benchmark_metric_reliability_trialcv));
if isempty(rel_vals)
    bench_rel_ylim = [0 1];
else
    bench_rel_ylim = [0, max(rel_vals) * 1.15];
end

slope_vals = benchmark_metric_separation_slope(isfinite(benchmark_metric_separation_slope));
if isempty(slope_vals)
    bench_slope_absmax = 1;
else
    bench_slope_absmax = prctile(abs(slope_vals), 98);
    if ~isfinite(bench_slope_absmax) || bench_slope_absmax <= 0
        bench_slope_absmax = max(abs(slope_vals));
    end
    bench_slope_absmax = max(bench_slope_absmax * 1.1, 0.1);
end

delta_vals = benchmark_metric_separation_delta(isfinite(benchmark_metric_separation_delta));
if isempty(delta_vals)
    bench_delta_absmax = 1;
else
    bench_delta_absmax = prctile(abs(delta_vals), 98);
    if ~isfinite(bench_delta_absmax) || bench_delta_absmax <= 0
        bench_delta_absmax = max(abs(delta_vals));
    end
    bench_delta_absmax = max(bench_delta_absmax * 1.1, 0.1);
end

%% Subject-level benchmark figures
bench_method_labels = {'Raw', 'Top-Eig GED', 'Combined GED (pre-Perm+CV)', 'Combined GED (post-Perm+CV)'};
bench_method_colors = [0.2 0.2 0.2; 0.1 0.35 0.75; 0.85 0.55 0.10; 0.75 0.2 0.1];
for subj = 1:nSubj
    fig_bench_subj = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('Four-Way Spectrum Benchmark: Subject %s', subjects{subj}), ...
        'FontSize', 18, 'FontWeight', 'bold');

    for mi = 1:nBenchmarkMethods
        subplot(2, 4, mi); hold on;
        cond_line_handles = gobjects(1, 4);
        for cond = 1:4
            pr_dt = all_trial_powratio_dt_bench{mi, cond, subj};
            if isempty(pr_dt), continue; end
            mu = nanmean(pr_dt, 1);
            se = nanstd(pr_dt, [], 1) / sqrt(max(1, size(pr_dt, 1)));
            faceC = 0.8 * colors(cond,:) + 0.2 * [1 1 1];
            patch([scan_freqs, fliplr(scan_freqs)], ...
                [mu - se, fliplr(mu + se)], ...
                colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.25);
            cond_line_handles(cond) = plot(scan_freqs, movmean(mu, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2.0);
        end
        yline(0, 'k-', 'LineWidth', 0.5);
        title(bench_method_labels{mi}, 'FontSize', 12, 'Color', bench_method_colors(mi,:));
        xlabel('Frequency [Hz]'); ylabel('\Delta Power Ratio');
        xlim([30 90]);  box on; set(gca, 'FontSize', 10);
        if mi == 1
            valid_handles = isgraphics(cond_line_handles);
            if any(valid_handles)
                legend(cond_line_handles(valid_handles), condLabels(valid_handles), ...
                    'Location', 'best', 'FontSize', 9);
            end
        end
    end

    subplot(2, 4, 5); hold on;
    single_med = nan(1, 4);
    single_lo = nan(1, 4);
    single_hi = nan(1, 4);
    gap_med = nan(1, 4);
    gap_lo = nan(1, 4);
    gap_hi = nan(1, 4);
    for cond = 1:4
        sp = all_trial_peaks_single{cond, subj};
        sp = sp(~isnan(sp));
        if ~isempty(sp)
            single_med(cond) = median(sp);
            sq = prctile(sp, [25 75]);
            single_lo(cond) = single_med(cond) - sq(1);
            single_hi(cond) = sq(2) - single_med(cond);
        end

        lp = all_trial_peaks_low{cond, subj};
        hp = all_trial_peaks_high{cond, subj};
        vgap = ~isnan(lp) & ~isnan(hp);
        if any(vgap)
            g = hp(vgap) - lp(vgap);
            gap_med(cond) = median(g);
            gq = prctile(g, [25 75]);
            gap_lo(cond) = gap_med(cond) - gq(1);
            gap_hi(cond) = gq(2) - gap_med(cond);
        end
    end
    yyaxis left
    b4 = bar(1:4, single_med, 0.58, 'FaceColor', 'flat', 'EdgeColor', 'none');
    b4.CData = colors;
    errorbar(1:4, single_med, single_lo, single_hi, 'k', 'LineStyle', 'none', 'LineWidth', 1.1, 'CapSize', 5);
    ylabel('Single-peak median [Hz]');
    ylim(bench_single_ylim);
    yyaxis right
    errorbar(1:4, gap_med, gap_lo, gap_hi, 'ko-', 'LineWidth', 1.4, 'MarkerFaceColor', [0.2 0.2 0.2], 'CapSize', 5);
    ylabel('H-L separation [Hz]');
    ylim(bench_gap_ylim);
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 10);
    d_single = single_med(4) - single_med(1);
    d_gap = gap_med(4) - gap_med(1);
    text(0.02, 0.96, sprintf('\\DeltaSingle_{100-25}=%.2f Hz\n\\DeltaGap_{100-25}=%.2f Hz', d_single, d_gap), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
    title('Single median + H-L gap (post-Perm+CV)');
    box on;

    subplot(2, 4, 6); hold on;
    prom_vec = nan(nBenchmarkMethods, 1);
    rel_vec = nan(nBenchmarkMethods, 1);
    for mi = 1:nBenchmarkMethods
        prom_vec(mi) = nanmean(squeeze(benchmark_metric_prominence(mi, :, subj)));
        rel_vec(mi) = nanmean(squeeze(benchmark_metric_reliability_trialcv(mi, :, subj)));
    end
    yyaxis left
    b1 = bar(1:nBenchmarkMethods, prom_vec, 0.38, 'FaceColor', 'flat');
    for mi = 1:nBenchmarkMethods
        b1.CData(mi, :) = bench_method_colors(mi, :);
    end
    ylabel('Peak prominence');
    ylim(bench_prom_ylim);
    yyaxis right
    plot(1:nBenchmarkMethods, rel_vec, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', [0.2 0.2 0.2]);
    ylabel('Trial CV (lower better)');
    ylim(bench_rel_ylim);
    set(gca, 'XTick', 1:nBenchmarkMethods, 'XTickLabel', bench_method_labels, 'XTickLabelRotation', 20, 'FontSize', 10);
    text(0.02, 0.96, sprintf('\\DeltaProm_{Post-R}=%.2f\n\\DeltaCV_{Post-R}=%.3f', ...
        prom_vec(4)-prom_vec(1), rel_vec(4)-rel_vec(1)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
    title('Prominence + Reliability');
     box on;

    subplot(2, 4, 7); hold on;
    slope_vec = squeeze(benchmark_metric_separation_slope(:, subj));
    delta_vec = squeeze(benchmark_metric_separation_delta(:, subj));
    yyaxis left
    b2 = bar(1:nBenchmarkMethods, slope_vec, 0.38, 'FaceColor', 'flat');
    for mi = 1:nBenchmarkMethods
        b2.CData(mi, :) = bench_method_colors(mi, :);
    end
    ylabel('Condition slope [Hz/cond]');
    ylim([-bench_slope_absmax bench_slope_absmax]);
    yyaxis right
    plot(1:nBenchmarkMethods, delta_vec, 'ks--', 'LineWidth', 2, 'MarkerFaceColor', [0.2 0.2 0.2]);
    ylabel('\Delta median (100%-25%) [Hz]');
    ylim([-bench_delta_absmax bench_delta_absmax]);
    set(gca, 'XTick', 1:nBenchmarkMethods, 'XTickLabel', bench_method_labels, 'XTickLabelRotation', 20, 'FontSize', 10);
    text(0.02, 0.96, sprintf('Post-Perm+CV slope=%.2f\nPost-Perm+CV \\Delta=%.2f Hz', ...
        slope_vec(4), delta_vec(4)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
    title('Condition separation');
     box on;

    subplot(2, 4, 8); axis off;
    text(0.05, 0.9, sprintf(['Method order:\n1) %s\n2) %s\n3) %s\n4) %s'], ...
        bench_method_labels{1}, bench_method_labels{2}, bench_method_labels{3}, bench_method_labels{4}), ...
        'FontSize', 10, 'VerticalAlignment', 'top');

    saveas(fig_bench_subj, fullfile(fig_save_dir, sprintf('GCP_eeg_GED_benchmark_subj%s.png', subjects{subj})));
end

%% ====================================================================
%  CENTROID METRIC: Subject/group summaries and concordance
%  ====================================================================
fig_cent = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Trial-Level Gamma Centroid', ...
    'FontSize', 18, 'FontWeight', 'bold');

subplot(1, 2, 1); hold on;
for s = 1:nSubj
    yc = all_trial_median_centroid(:, s);
    if sum(~isnan(yc)) >= 2
        plot(1:4, yc, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    end
end
for c = 1:4
    vals = all_trial_median_centroid(c, :);
    vals = vals(~isnan(vals));
    if ~isempty(vals)
        xj = c + (rand(size(vals)) - 0.5) * 0.12;
        scatter(xj, vals, 110, colors(c,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
    end
end
mu_c = nanmean(all_trial_median_centroid, 2);
sem_c = nanstd(all_trial_median_centroid, [], 2) ./ sqrt(sum(~isnan(all_trial_median_centroid), 2));
errorbar(1:4, mu_c, sem_c, 'k', 'LineWidth', 2, 'CapSize', 10);
plot(1:4, mu_c, 'k-', 'LineWidth', 2.5);
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 13, 'Box', 'off');
xlim([0.3 4.7]); 
ylim([50 70]); 
ylabel('Centroid Frequency [Hz]');
title('Subject medians by condition', 'FontSize', 14, 'FontWeight', 'bold');

subplot(1, 2, 2); hold on;
y_all_c = [];
g_all_c = [];
for c = 1:4
    for s = 1:nSubj
        tc = all_trial_centroid{c, s};
        if ~isempty(tc)
            tc = tc(~isnan(tc));
            y_all_c = [y_all_c; tc(:)];
            g_all_c = [g_all_c; c * ones(length(tc), 1)];
        end
    end
end
if ~isempty(y_all_c)
    boxplot(y_all_c, g_all_c, 'Colors', 'k', 'Symbol', '', 'Widths', 0.15);
    for c = 1:4
        vals = y_all_c(g_all_c == c);
        xj = c + 0.15 + (rand(size(vals)) - 0.5) * 0.22;
        scatter(xj, vals, 10, colors(c,:), 'filled', 'MarkerFaceAlpha', 0.2);
    end
end
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 13, 'Box', 'off');
xlim([0.3 4.7]);
ylim(centroid_freq_range); 
ylabel('Centroid Frequency [Hz]');
title('All trials pooled', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig_cent, fullfile(fig_save_dir, 'GCP_eeg_GED_centroid_summary.png'));

%% ====================================================================
%  METHOD COMPARISON: Three dual-peak assignment methods
%  ====================================================================
fprintf('\nRunning method comparison...\n');

method_names = {'Hard 50 Hz', 'Spectral Trough', 'GMM'};

% Per-method storage: m_peaks_{low,high}{method}{cond, subj}
m_peaks_low  = cell(1, 3);
m_peaks_high = cell(1, 3);
for mi = 1:3
    m_peaks_low{mi}  = cell(4, nSubj);
    m_peaks_high{mi} = cell(4, nSubj);
end
m_median_low   = nan(4, nSubj, 3);
m_median_high  = nan(4, nSubj, 3);
m_detrate_low  = nan(4, nSubj, 3);
m_detrate_high = nan(4, nSubj, 3);

trough_freqs = nan(1, nSubj);
gmm_means    = nan(2, nSubj);
gmm_boundary = nan(1, nSubj);

for subj = 1:nSubj
    clc
    fprintf('Method comparison: Subject %s (%d/%d)\n', subjects{subj}, subj, nSubj);

    %% --- Method 1: Hard 50 Hz (already computed in Phase 2) ---
    for cond = 1:4
        m_peaks_low{1}{cond, subj}  = all_trial_peaks_low{cond, subj};
        m_peaks_high{1}{cond, subj} = all_trial_peaks_high{cond, subj};
    end

    %% --- Method 2: Per-subject spectral trough ---
    all_pr_dt = [];
    for cond = 1:4
        pr_mat = all_trial_powratio{cond, subj};
        if ~isempty(pr_mat)
            nTrl_m = size(pr_mat, 1);
            dt_mat = nan(size(pr_mat));
            for trl = 1:nTrl_m
                p = polyfit(scan_freqs, pr_mat(trl,:), poly_order);
                dt_mat(trl,:) = pr_mat(trl,:) - polyval(p, scan_freqs);
            end
            all_pr_dt = [all_pr_dt; dt_mat];
        end
    end

    mean_dt = nanmean(all_pr_dt, 1);
    mean_dt_smooth = movmean(mean_dt, 5);

    [avg_pks, avg_locs] = findpeaks(mean_dt_smooth, scan_freqs, 'MinPeakDistance', 5);
    avg_pos_mask = avg_pks > 0;
    avg_pks_pos  = avg_pks(avg_pos_mask);
    avg_locs_pos = avg_locs(avg_pos_mask);

    split_freq_m2 = 50;
    if length(avg_pks_pos) >= 2
        [~, sortI] = sort(avg_pks_pos, 'descend');
        pk1_loc = avg_locs_pos(sortI(1));
        pk2_loc = avg_locs_pos(sortI(2));
        lo_pk = min(pk1_loc, pk2_loc);
        hi_pk = max(pk1_loc, pk2_loc);
        lo_idx = find(scan_freqs == lo_pk);
        hi_idx = find(scan_freqs == hi_pk);
        if hi_idx > lo_idx
            [~, min_rel_idx] = min(mean_dt_smooth(lo_idx:hi_idx));
            split_freq_m2 = scan_freqs(lo_idx + min_rel_idx - 1);
        end
    end
    trough_freqs(subj) = split_freq_m2;

    for cond = 1:4
        pr_mat = all_trial_powratio{cond, subj};
        if isempty(pr_mat), continue; end
        nTrl_m = size(pr_mat, 1);
        trl_lo = nan(nTrl_m, 1);
        trl_hi = nan(nTrl_m, 1);

        for trl = 1:nTrl_m
            pr = pr_mat(trl, :);
            if all(isnan(pr)), continue; end
            p = polyfit(scan_freqs, pr, poly_order);
            pr_dt_smooth = movmean(pr - polyval(p, scan_freqs), 5);

            [pks_a, locs_a] = findpeaks(pr_dt_smooth, scan_freqs, 'MinPeakDistance', 5);
            pm = pks_a > 0;
            pp = pks_a(pm);
            lp = locs_a(pm);

            mask_lo = lp < split_freq_m2;
            mask_hi = lp >= split_freq_m2;
            if any(mask_lo)
                [~, bi] = max(pp(mask_lo));
                tmp = lp(mask_lo); trl_lo(trl) = tmp(bi);
            end
            if any(mask_hi)
                [~, bi] = max(pp(mask_hi));
                tmp = lp(mask_hi); trl_hi(trl) = tmp(bi);
            end
        end
        m_peaks_low{2}{cond, subj}  = trl_lo;
        m_peaks_high{2}{cond, subj} = trl_hi;
    end

    %% --- Method 3: GMM on single-peak frequencies ---
    all_single = [];
    for cond = 1:4
        sp = all_trial_peaks_single{cond, subj};
        if ~isempty(sp)
            all_single = [all_single; sp(~isnan(sp))];
        end
    end

    gmm_ok = false;
    if length(all_single) >= 20
        try
            gm = fitgmdist(all_single, 2, 'RegularizationValue', 0.1, ...
                'Options', statset('MaxIter', 500));
            [mu1, mu2] = deal(gm.mu(1), gm.mu(2));
            if abs(mu1 - mu2) >= 5
                gmm_ok = true;
                gmm_means(:, subj) = sort(gm.mu);
                gmm_boundary(subj) = mean(gm.mu);
            end
        catch
        end
    end

    for cond = 1:4
        sp = all_trial_peaks_single{cond, subj};
        if isempty(sp), continue; end
        nTrl_m = length(sp);
        trl_lo = nan(nTrl_m, 1);
        trl_hi = nan(nTrl_m, 1);

        if gmm_ok
            lo_comp_mu = min(gm.mu);
            for trl = 1:nTrl_m
                if isnan(sp(trl)), continue; end
                post = posterior(gm, sp(trl));
                [~, cluster] = max(post);
                if gm.mu(cluster) == lo_comp_mu
                    trl_lo(trl) = sp(trl);
                else
                    trl_hi(trl) = sp(trl);
                end
            end
        else
            for trl = 1:nTrl_m
                if isnan(sp(trl)), continue; end
                if sp(trl) < 50
                    trl_lo(trl) = sp(trl);
                else
                    trl_hi(trl) = sp(trl);
                end
            end
        end
        m_peaks_low{3}{cond, subj}  = trl_lo;
        m_peaks_high{3}{cond, subj} = trl_hi;
    end
end

% Compute aggregates for all methods
for mi = 1:3
    for subj = 1:nSubj
        for cond = 1:4
            lo = m_peaks_low{mi}{cond, subj};
            hi = m_peaks_high{mi}{cond, subj};
            if ~isempty(lo)
                valid_lo = ~isnan(lo);
                if any(valid_lo)
                    m_median_low(cond, subj, mi)  = median(lo(valid_lo));
                    m_detrate_low(cond, subj, mi) = sum(valid_lo) / length(lo);
                end
            end
            if ~isempty(hi)
                valid_hi = ~isnan(hi);
                if any(valid_hi)
                    m_median_high(cond, subj, mi)  = median(hi(valid_hi));
                    m_detrate_high(cond, subj, mi) = sum(valid_hi) / length(hi);
                end
            end
        end
    end
end

%% ====================================================================
%  GRAND-AVERAGE FOUR-WAY BENCHMARK
%  ====================================================================
fig_bench_group = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Four-Way Spectrum Benchmark: Grand Average', 'FontSize', 18, 'FontWeight', 'bold');

for mi = 1:nBenchmarkMethods
    subplot(2, 4, mi); hold on;
    cond_line_handles = gobjects(1, 4);
    panel_maxabs = 0;
    for cond = 1:4
        subj_curves = nan(nSubj, nFreqs);
        for s = 1:nSubj
            pr_dt = all_trial_powratio_dt_bench{mi, cond, s};
            if ~isempty(pr_dt)
                subj_mu = nanmean(pr_dt, 1);
                subj_curves(s, :) = normalize_maxabs_curve(subj_mu);
            end
        end
        mu = nanmean(subj_curves, 1);
        se = nanstd(subj_curves, [], 1) ./ sqrt(sum(~isnan(subj_curves(:,1))));
        panel_maxabs = max(panel_maxabs, max(abs([mu - se, mu + se]), [], 'omitnan'));
        faceC = 0.8 * colors(cond,:) + 0.2 * [1 1 1];
        patch([scan_freqs, fliplr(scan_freqs)], ...
            [mu - se, fliplr(mu + se)], ...
            colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.25);
        cond_line_handles(cond) = plot(scan_freqs, movmean(mu, 5), '-', ...
            'Color', colors(cond,:), 'LineWidth', 2.2);
    end
    yline(0, 'k-', 'LineWidth', 0.5);
    title(bench_method_labels{mi}, 'FontSize', 12, 'Color', bench_method_colors(mi,:));
    xlabel('Frequency [Hz]'); ylabel('\Delta Power Ratio');
    panel_maxabs = max(panel_maxabs, eps);
    ylim([-panel_maxabs panel_maxabs]);
    xlim([30 90]);  box on; set(gca, 'FontSize', 10);
    if mi == 1
        valid_handles = isgraphics(cond_line_handles);
        if any(valid_handles)
            legend(cond_line_handles(valid_handles), condLabels(valid_handles), ...
                'Location', 'best', 'FontSize', 9);
        end
    end
end

subplot(2, 4, 5); hold on;
single_mu = nanmean(all_trial_median_single, 2);
single_se = nanstd(all_trial_median_single, [], 2) ./ sqrt(sum(~isnan(all_trial_median_single), 2));
single_med = nanmedian(all_trial_median_single, 2);
gap_mu = nanmean(all_trial_median_gap, 2);
gap_se = nanstd(all_trial_median_gap, [], 2) ./ sqrt(sum(~isnan(all_trial_median_gap), 2));
gap_med = nanmedian(all_trial_median_gap, 2);
yyaxis left
b4g = bar(1:4, single_mu, 0.58, 'FaceColor', 'flat', 'EdgeColor', 'none');
b4g.CData = colors;
errorbar(1:4, single_mu, single_se, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);
scatter(1:4, single_med, 35, 'kd', 'filled');
ylabel('Single-peak mean [Hz]');
ylim(bench_single_ylim);
yyaxis right
errorbar(1:4, gap_mu, gap_se, 'ko-', 'LineWidth', 1.6, 'MarkerFaceColor', [0.2 0.2 0.2], 'CapSize', 6);
scatter(1:4, gap_med, 30, 'ks', 'filled');
ylabel('H-L separation mean [Hz]');
ylim(bench_gap_ylim);
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 10);
text(0.02, 0.96, sprintf('\\DeltaSingle_{100-25}=%.2f Hz\n\\DeltaGap_{100-25}=%.2f Hz', ...
    single_mu(4)-single_mu(1), gap_mu(4)-gap_mu(1)), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
title('Single mean + H-L gap (post-Perm+CV)');
box on;

subplot(2, 4, 6); hold on;
prom_group = squeeze(nanmean(benchmark_metric_prominence, 2));
prom_mu = nanmean(prom_group, 2);
prom_se = nanstd(prom_group, [], 2) ./ sqrt(sum(~isnan(prom_group), 2));
prom_med = nanmedian(prom_group, 2);
rel_group = squeeze(nanmean(benchmark_metric_reliability_trialcv, 2));
rel_mu = nanmean(rel_group, 2);
rel_se = nanstd(rel_group, [], 2) ./ sqrt(sum(~isnan(rel_group), 2));
rel_med = nanmedian(rel_group, 2);
yyaxis left
bar(1:nBenchmarkMethods, prom_mu, 0.38, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
errorbar(1:nBenchmarkMethods, prom_mu, prom_se, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);
scatter(1:nBenchmarkMethods, prom_med, 25, 'kd', 'filled');
ylabel('Peak prominence');
ylim(bench_prom_ylim);
yyaxis right
errorbar(1:nBenchmarkMethods, rel_mu, rel_se, 'ko-', 'LineWidth', 1.8, 'MarkerFaceColor', [0.2 0.2 0.2]);
scatter(1:nBenchmarkMethods, rel_med, 25, 'ks', 'filled');
ylabel('Trial CV (lower better)');
ylim(bench_rel_ylim);
set(gca, 'XTick', 1:nBenchmarkMethods, 'XTickLabel', bench_method_labels, 'XTickLabelRotation', 20, 'FontSize', 10);
text(0.02, 0.96, sprintf('\\DeltaProm_{Post-R}=%.2f\n\\DeltaCV_{Post-R}=%.3f', ...
    prom_mu(4)-prom_mu(1), rel_mu(4)-rel_mu(1)), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
title('Prominence + Reliability');
 box on;

subplot(2, 4, 7); hold on;
slope_mu = nanmean(benchmark_metric_separation_slope, 2);
slope_se = nanstd(benchmark_metric_separation_slope, [], 2) ./ sqrt(sum(~isnan(benchmark_metric_separation_slope), 2));
delta_mu = nanmean(benchmark_metric_separation_delta, 2);
delta_se = nanstd(benchmark_metric_separation_delta, [], 2) ./ sqrt(sum(~isnan(benchmark_metric_separation_delta), 2));
yyaxis left
bar(1:nBenchmarkMethods, slope_mu, 0.38, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
errorbar(1:nBenchmarkMethods, slope_mu, slope_se, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);
ylabel('Condition slope [Hz/cond]');
ylim([-bench_slope_absmax bench_slope_absmax]);
yyaxis right
errorbar(1:nBenchmarkMethods, delta_mu, delta_se, 'ks--', 'LineWidth', 1.8, 'MarkerFaceColor', [0.2 0.2 0.2]);
ylabel('\Delta median (100%-25%) [Hz]');
ylim([-bench_delta_absmax bench_delta_absmax]);
set(gca, 'XTick', 1:nBenchmarkMethods, 'XTickLabel', bench_method_labels, 'XTickLabelRotation', 20, 'FontSize', 10);
text(0.02, 0.96, sprintf('Post-Perm+CV slope=%.2f\nPost-Perm+CV \\Delta=%.2f Hz', ...
    slope_mu(4), delta_mu(4)), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
title('Condition separation');
 box on;

subplot(2, 4, 8); axis off;
text(0.05, 0.9, sprintf(['Method order:\n1) %s\n2) %s\n3) %s\n4) %s'], ...
    bench_method_labels{1}, bench_method_labels{2}, bench_method_labels{3}, bench_method_labels{4}), ...
    'FontSize', 10, 'VerticalAlignment', 'top');

saveas(fig_bench_group, fullfile(fig_save_dir, 'GCP_eeg_GED_benchmark_grandaverage.png'));

%% ====================================================================
%  STANDALONE CONDITION-SEPARATION METRICS (post-Perm+CV combined GED)
%  ====================================================================
close all
fig_cond_slope = figure('Position', [0 0 1512 982], 'Color', 'w');

post_idx = 4; % Combined GED (post-Perm+CV)
slope_post = benchmark_metric_separation_slope(post_idx, :);
delta_post = benchmark_metric_separation_delta(post_idx, :);

tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% ---------- Panel 1: condition slope ----------
nexttile; hold on;
valid_slope = isfinite(slope_post);
slope_vals = slope_post(valid_slope);
if ~isempty(slope_vals)
    boxplot(slope_vals(:), ones(numel(slope_vals), 1), 'Colors', 'k', ...
        'Symbol', '', 'Widths', 0.2);
    xj = 1 + (rand(numel(slope_vals), 1) - 0.5) * 0.10;
    scatter(xj, slope_vals(:), 250, [0.35 0.35 0.35], 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'k', 'LineWidth', 1);
end
yline(0, 'k--', 'LineWidth', 1.0);
xlim([0.45 1.45]);
ylim([-2 2]);
set(gca, 'XTick', 1, 'XTickLabel', {'Combined GED (post-Perm+CV)'}, ...
    'FontSize', 16, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
ylabel('Slope across contrast conditions [Hz/condition]', 'FontSize', 18, 'FontWeight', 'bold');
title('Contrast Condition Slope', 'FontSize', 20, 'FontWeight', 'bold');

% ---------- Panel 2: median shift (100%-25%) ----------
nexttile; hold on;
valid_delta = isfinite(delta_post);
delta_vals = delta_post(valid_delta);
if ~isempty(delta_vals)
    boxplot(delta_vals(:), ones(numel(delta_vals), 1), 'Colors', 'k', ...
        'Symbol', '', 'Widths', 0.2);
    xj = 1 + (rand(numel(delta_vals), 1) - 0.5) * 0.10;
    scatter(xj, delta_vals(:), 250, [0.35 0.35 0.35], 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'k', 'LineWidth', 1);
end
yline(0, 'k--', 'LineWidth', 1.0);
xlim([0.45 1.45]);
ylim([-6 6]);
set(gca, 'XTick', 1, 'XTickLabel', {'Combined GED (post-Perm+CV)'}, ...
    'FontSize', 16, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
ylabel('\Delta median (100% - 25%) [Hz]', 'FontSize', 18, 'FontWeight', 'bold');
title('Median Frequency Shift (100% - 25%)', 'FontSize', 20, 'FontWeight', 'bold');

saveas(fig_cond_slope, fullfile(fig_save_dir, 'GCP_eeg_GED_condition_slope.png'));

%% ====================================================================
%  SUMMARY DASHBOARD (backprojected combined-component data)
%  ====================================================================
fig_summary = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('GED Trials Summary Dashboard (Perm+CV backprojected)', ...
    'FontSize', 16, 'FontWeight', 'bold');

summary_metrics = { ...
    all_trial_median_single, ...
    all_trial_median_low, ...
    all_trial_median_high, ...
    all_trial_median_gap, ...
    squeeze(benchmark_metric_prominence(4, :, :)), ...
    squeeze(benchmark_metric_reliability_trialcv(4, :, :)), ...
    all_trial_median_centroid, ...
    all_trial_detrate_centroid * 100};
summary_names = {'Single median [Hz]', 'Low median [Hz]', 'High median [Hz]', ...
    'H-L separation [Hz]', 'Prominence', 'Trial CV', ...
    'Centroid median [Hz]', 'Centroid detectability [%]'};

for mi = 1:numel(summary_metrics)
    subplot(2, 4, mi); hold on;
    dat = summary_metrics{mi};
    mu = nanmean(dat, 2);
    sem = nanstd(dat, [], 2) ./ sqrt(sum(~isnan(dat), 2));
    med = nanmedian(dat, 2);
    for c = 1:4
        bar(c, mu(c), 0.6, 'FaceColor', colors(c,:), 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    end
    errorbar(1:4, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.1, 'CapSize', 5);
    scatter(1:4, med, 28, 'kd', 'filled');
    for s = 1:nSubj
        for c = 1:4
            if ~isnan(dat(c, s))
                scatter(c + (rand - 0.5) * 0.18, dat(c, s), 20, [0.5 0.5 0.5], ...
                    'filled', ...
                    'MarkerFaceAlpha', 0.5, ...
                    'MarkerEdgeColor', [1 1 1], ...    % white ring
                    'LineWidth', 0.5, ...              
                    'MarkerEdgeAlpha', 0.9);                       
            end
        end
    end
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 9);
    title(summary_names{mi}, 'FontSize', 10, 'FontWeight', 'bold');
    box on;
    if mi == 1
        ylim([45 70]);
    elseif mi == 2
        ylim([35 50]);
    elseif mi == 3
        ylim([50 75]);
    elseif mi == 4
        ylim([22.5 36]);
    elseif mi == 5
        ylim([0 10]);
    elseif mi == 7
        ylim([50 60]);
    end
end
saveas(fig_summary, fullfile(fig_save_dir, 'GCP_eeg_GED_metrics_summary.png'));

%% ====================================================================
%  GRAND AVERAGE: Single Peak (subject-level median, summary-style)
%  ====================================================================
close all
fprintf('\nCreating grand average figures...\n');

fig_box1 = figure('Position', [0 0 1512/2 982], 'Color', 'w');
hold on;

dat = all_trial_median_single;
mu = nanmean(dat, 2);
sem = nanstd(dat, [], 2) ./ sqrt(sum(~isnan(dat), 2));
med = nanmedian(dat, 2);

for c = 1:4
    bar(c, mu(c), 0.6, 'FaceColor', colors(c,:), 'EdgeColor', 'k', 'FaceAlpha', 0.7);
end
errorbar(1:4, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.1, 'CapSize', 5);
scatter(1:4, med, 60, 'kd', 'filled');

for s = 1:nSubj
    for c = 1:4
        if ~isnan(dat(c, s))
            scatter(c + (rand - 0.5) * 0.25, dat(c, s), 250, [0.5 0.5 0.5], ...
                'filled', ...
                'MarkerFaceAlpha', 0.5, ...
                'MarkerEdgeColor', [1 1 1], ...
                'LineWidth', 0.5, ...
                'MarkerEdgeAlpha', 0.9);
        end
    end
end

xlim([0.3 4.7]);
ylim([45 75]);
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 16, 'Box', 'off');
ylabel('Gamma Frequency [Hz]');
title('Gamma Frequency over Conditions', 'FontSize', 18, 'FontWeight', 'bold');

saveas(fig_box1, fullfile(fig_save_dir, 'GCP_eeg_GED_boxplot_SinglePeak.png'));

%% ====================================================================
%  GRAND AVERAGE: Single Peak all trials pooled — raincloud
%  ====================================================================
fig_trl1 = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

y_all = []; g_all = [];
for c = 1:4
    for s = 1:nSubj
        tpk = all_trial_peaks_single{c, s};
        if ~isempty(tpk)
            valid = tpk(~isnan(tpk));
            y_all = [y_all; valid(:)];
            g_all = [g_all; ones(length(valid), 1) * c];
        end
    end
end

for c = 1:4
    vals = y_all(g_all == c);
    if length(vals) >= 5
        [f_dens, xi] = ksdensity(vals);
        f_dens = f_dens / max(f_dens) * 0.35;
        patch(c - f_dens - 0.05, xi, colors(c,:), ...
            'FaceAlpha', 0.3, 'EdgeColor', colors(c,:), 'LineWidth', 1);
    end
end

if ~isempty(y_all)
    boxplot(y_all, g_all, 'Colors', 'k', 'Symbol', '', 'Widths', 0.15);
end
hold on;
for c = 1:4
    mask = g_all == c;
    vals = y_all(mask);
    xJit = c + 0.15 + (rand(size(vals)) - 0.5) * 0.25;
    scatter(xJit, vals, 12, colors(c,:), 'filled', 'MarkerFaceAlpha', 0.2);
end

xlim([0.3 4.7]); ylim([30 90]);
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 18, 'Box', 'off');
ylabel('Peak Gamma Frequency [Hz]');
title('Single Peak: All Trials (pooled across subjects)', ...
    'FontSize', 18, 'FontWeight', 'bold');

saveas(fig_trl1, fullfile(fig_save_dir, 'GCP_eeg_GED_boxplot_alltrials_SinglePeak.png'));

%% ====================================================================
%  GRAND AVERAGE: All-subjects subplot (mean trial spectra)
%  ====================================================================
nRows = ceil(nSubj / 5);
fig_all = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Trial-Level Mean Detrended Spectra: All Subjects (N=%d)', nSubj), ...
    'FontSize', 16, 'FontWeight', 'bold');
fig_all_legend_handles = gobjects(1, 4);

for s = 1:nSubj
    subplot(nRows, 5, s); hold on;
    subj_panel_maxabs = 0;
    for cond = 1:4
        pr_mat = all_trial_powratio{cond, s};
        if ~isempty(pr_mat)
            nTrl = size(pr_mat, 1);
            pr_dt_mat = nan(size(pr_mat));
            for trl = 1:nTrl
                p = polyfit(scan_freqs, pr_mat(trl,:), poly_order);
                pr_dt_mat(trl,:) = pr_mat(trl,:) - polyval(p, scan_freqs);
            end
            mu_dt = nanmean(pr_dt_mat, 1);
            mu_dt = normalize_maxabs_curve(mu_dt);
            subj_panel_maxabs = max(subj_panel_maxabs, max(abs(mu_dt), [], 'omitnan'));
            h_cond = plot(scan_freqs, movmean(mu_dt, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2);
            if s == 1
                fig_all_legend_handles(cond) = h_cond;
            end
            md_pf = all_trial_median_single(cond, s);
            if ~isnan(md_pf)
                xline(md_pf, '--', 'Color', colors(cond,:), 'LineWidth', 1.2);
            end
        end
    end
    yline(0, 'k-', 'LineWidth', 0.5);
    subj_panel_maxabs = max(subj_panel_maxabs, eps);
    xlabel('Freq [Hz]'); ylabel('\Delta PR');
    title(sprintf('Subj %s', subjects{s}), 'FontSize', 11);
    set(gca, 'FontSize', 10); xlim([30 90]); ylim([-subj_panel_maxabs subj_panel_maxabs]); box on;
    if s == 1
        valid_handles = isgraphics(fig_all_legend_handles);
        if any(valid_handles)
            legend(fig_all_legend_handles(valid_handles), condLabels(valid_handles), ...
                'FontSize', 8, 'Location', 'northwest');
        end
    end
end
saveas(fig_all, fullfile(fig_save_dir, 'GCP_eeg_GED_all_subjects.png'));

%% ====================================================================
%  DETECTION RATE FIGURE (primary: hard 50 Hz)
%  ====================================================================
fig_det = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Trial-Level Peak Detection Rate (mean across subjects)', ...
    'FontSize', 18, 'FontWeight', 'bold');

det_data   = {all_trial_detrate_single, all_trial_detrate_low, all_trial_detrate_high};
det_labels = {'Single Peak', 'Low Gamma', 'High Gamma'};

for di = 1:3
    subplot(1, 3, di); hold on;
    dr = det_data{di};

    mu_dr  = nanmean(dr, 2) * 100;
    sem_dr = nanstd(dr, [], 2) / sqrt(nSubj) * 100;

    b = bar(1:4, mu_dr, 0.6);
    b.FaceColor = 'flat';
    for c = 1:4
        b.CData(c,:) = colors(c,:);
    end
    errorbar(1:4, mu_dr, sem_dr, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 10);

    for s = 1:nSubj
        for c = 1:4
            if ~isnan(dr(c, s))
                scatter(c + (rand-0.5)*0.2, dr(c, s)*100, 30, [0.4 0.4 0.4], ...
                    'filled', 'MarkerFaceAlpha', 0.5);
            end
        end
    end

    ylim([0 105]);
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 14, 'Box', 'off');
    ylabel('Detection Rate [%]');
    title(det_labels{di}, 'FontSize', 16, 'FontWeight', 'bold');
    
end

saveas(fig_det, fullfile(fig_save_dir, 'GCP_eeg_GED_detection_rate.png'));

%% Save results
if ispc
    save_path = 'W:\Students\Arne\GCP\data\features\GCP_eeg_GED.mat';
else
    save_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_eeg_GED.mat';
end
save(save_path, ...
    'all_trial_powratio', ...
    'all_trial_peaks_single', 'all_trial_peaks_low', 'all_trial_peaks_high', ...
    'all_trial_centroid', ...
    'all_trial_mean_single', 'all_trial_median_single', ...
    'all_trial_mean_low', 'all_trial_median_low', ...
    'all_trial_mean_high', 'all_trial_median_high', 'all_trial_median_gap', ...
    'all_trial_mean_centroid', 'all_trial_median_centroid', ...
    'all_trial_detrate_single', 'all_trial_detrate_low', 'all_trial_detrate_high', ...
    'all_trial_detrate_centroid', ...
    'all_topos', 'all_topo_labels', 'all_eigenvalues', ...
    'all_selected_comp_idx', 'all_selected_comp_corr', 'all_selected_comp_eval', ...
    'all_selected_comp_indices_multi', 'all_selected_comp_weights', ...
    'all_preperm_comp_indices_multi', 'all_preperm_comp_weights', ...
    'all_component_selection_stats', 'all_cv_curves', 'all_perm_thresholds', ...
    'all_trial_powratio_bench', 'all_trial_powratio_dt_bench', ...
    'benchmark_methods', 'raw_reference_definition', ...
    'benchmark_metric_detectability', 'benchmark_metric_prominence', ...
    'benchmark_metric_separation_slope', 'benchmark_metric_separation_delta', ...
    'benchmark_metric_reliability_trialcv', 'benchmark_metric_reliability_subjspread', ...
    'all_top5_corrs', 'all_top5_evals', 'all_top5_topos', 'all_simulated_templates', ...
    'm_peaks_low', 'm_peaks_high', 'm_median_low', 'm_median_high', ...
    'm_detrate_low', 'm_detrate_high', ...
    'trough_freqs', 'gmm_means', 'gmm_boundary', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames');

clc
fprintf('Done.\n');

%% Quick trial-level GLMM check (single-peak gamma)
fprintf('\n============================================================\n');
fprintf('Trial-Level GLMM Check: GammaFrequency ~ Condition + (1|subjectID)\n');
fprintf('============================================================\n');

glmm_gamma = [];
glmm_cond  = {};
glmm_subj  = {};

for subj = 1:nSubj
    for cond = 1:4
        trl_freqs = all_trial_peaks_single{cond, subj};
        if isempty(trl_freqs)
            continue;
        end
        valid_mask = isfinite(trl_freqs);
        if any(valid_mask)
            n_valid = sum(valid_mask);
            glmm_gamma = [glmm_gamma; trl_freqs(valid_mask)];
            glmm_cond  = [glmm_cond; repmat(condLabels(cond), n_valid, 1)];
            glmm_subj  = [glmm_subj; repmat(subjects(subj), n_valid, 1)];
        end
    end
end

if isempty(glmm_gamma)
    fprintf('GLMM skipped: no valid trial-level gamma frequencies were found.\n');
else
    tbl_glmm = table( ...
        glmm_gamma, ...
        categorical(glmm_cond, condLabels, 'Ordinal', true), ...
        categorical(glmm_subj), ...
        'VariableNames', {'GammaFrequency', 'Condition', 'subjectID'});

    fprintf('Observations: %d\n', height(tbl_glmm));
    fprintf('Subjects: %d\n', numel(unique(tbl_glmm.subjectID)));
    fprintf('Condition counts:\n');
    cond_levels = categories(tbl_glmm.Condition);
    cond_counts = zeros(numel(cond_levels), 1);
    for ci = 1:numel(cond_levels)
        cond_counts(ci) = sum(tbl_glmm.Condition == cond_levels{ci});
    end
    disp(table(categorical(cond_levels, cond_levels, 'Ordinal', true), cond_counts, ...
        'VariableNames', {'Condition', 'nTrials'}));

    try
        glmm_gamma_model = fitglme( ...
            tbl_glmm, ...
            'GammaFrequency ~ Condition + (1|subjectID)', ...
            'Distribution', 'Normal', ...
            'Link', 'Identity', ...
            'FitMethod', 'Laplace');

        fprintf('\nModel fit summary:\n');
        disp(glmm_gamma_model);
        fprintf('\nFixed-effects coefficients:\n');
        disp(glmm_gamma_model.Coefficients);

        coef_tbl = glmm_gamma_model.Coefficients;
        coef_names = cellstr(string(coef_tbl.Name));
        is_cond_coef = startsWith(coef_names, 'Condition_');
        if any(is_cond_coef)
            cond_names = coef_names(is_cond_coef);
            p_raw = coef_tbl.pValue(is_cond_coef);
            p_fdr = bh_fdr_adjust(p_raw);
            sig_fdr = p_fdr < 0.05;

            fprintf('\nFDR-corrected p-values (Benjamini-Hochberg, q=0.05):\n');
            disp(table(cond_names, p_raw, p_fdr, sig_fdr, ...
                'VariableNames', {'Contrast', 'pRaw', 'pFDR', 'isSignificantFDR'}));
        else
            fprintf('\nNo condition contrasts found for FDR correction.\n');
        end
    catch ME
        fprintf('GLMM fit failed: %s\n', ME.message);
    end
end

function vec_out = normalize_maxabs_curve(vec_in)
vec_out = vec_in;
if isempty(vec_in)
    return;
end
scale = max(abs(vec_in(~isnan(vec_in))));
if ~isempty(scale) && isfinite(scale) && scale > 0
    vec_out = vec_in ./ scale;
end
end

function p_adj = bh_fdr_adjust(p_raw)
p = p_raw(:);
n = numel(p);
p_adj = nan(size(p));
if n == 0
    return;
end

[p_sorted, sort_idx] = sort(p);
rank_idx = (1:n)';
q_sorted = p_sorted .* n ./ rank_idx;
q_sorted = min(q_sorted, 1);

for k = n-1:-1:1
    q_sorted(k) = min(q_sorted(k), q_sorted(k+1));
end

p_adj(sort_idx) = q_sorted;
p_adj = reshape(p_adj, size(p_raw));
end
