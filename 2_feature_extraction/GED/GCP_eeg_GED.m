%% GCP Trial-Level GED
%
% Uses the same spatial filter approach as GCP_eeg_GED_subjects.m, with the
% key difference that peak detection is performed on individual trials
% rather than on condition-averaged spectra.
%
% Pipeline:
%   Phase 1 — Pool all conditions; broadband GED; derive common spatial filter.
%   Phase 2 — Per condition and trial: narrowband scan (30–90 Hz), power ratio,
%             polynomial detrending, peak detection. Single-peak model per trial.
%   Component comparison — Three branches: raw channel-space reference,
%             top GED component, and combined hard-eligible GED components.
%   Output — Mean and median peak frequency per condition per subject;
%             centroid metric; detectability, separation, reliability metrics.

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');

%%
nSubj = length(subjects);

% Time windows
baseline_window = [-1.5, -0.25];
full_window = [0, 2.0];
early_window = [0, 0.6];
late_window = [1.0, 2.0];
stimulus_window = full_window; % kept for backward compatibility with existing code paths

% Gamma frequency range
gamma_range = [30, 90];

% Narrowband scanning parameters
scan_freqs = 30:1:90;
nFreqs     = length(scan_freqs);
scan_width = 3;

% GED parameters
% Trial-limited per-subject setting:
% These values are tuned for stability when each subject has limited trial
% counts. If future data only add subjects (not trials/condition), these
% settings should generally remain unchanged.
lambda = 0.05;              % full window [0, 2 s]
lambda_full  = 0.05;        % full window [0, 2 s]
lambda_early = 0.10;        % early window [0, 0.6 s] — stronger regularization (fewer samples)
lambda_late  = 0.08;       % late window [1, 2 s] — moderate (1 s of data)
ged_search_n = 20;          % search first N GED components
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
viz_topo_prctile = 99.9;            % robust percentile for color scaling
min_occfront_ratio = 0.8;          % tightened hard threshold: occipital must exceed frontal
min_eigval_hard = 1.1;              % hard minimum GED eigenvalue
min_corr_hard = 0.1;                % tightened hard minimum template correlation
min_gamma_hard = 0.1;               % tightened hard minimum gamma evidence (log ratio)
exclude_dominant_outlier = true;    % remove dominant top outlier (lambda1/lambda2 + MAD)
outlier_ratio_thr = 8.0;            % lambda1/lambda2 threshold for dominant-outlier detection
outlier_mad_mult = 4.0;             % MAD multiplier on log-eigenvalue distance
outlier_min_rest = 0;               % minimum non-top eligible components for dominant-outlier exclusion (0 = no minimum)
threshold_inspection_targets = [1 3 5]; % requested valid-component targets for threshold inspection

random_seed = 13;                    % reproducible randomization

% Run mode: 'full' = full pipeline; 'component_check' = GED + component selection only (Phase 1)
run_mode = 'full';

% Sensitivity analysis settings
sensitivity_enable = true;
sensitivity_scale = [0.8 1.0 1.2]; % +/-20% for threshold sweeps
sensitivity_poly_orders = [];

% Three-way component selection config (ordered for plotting/metrics):
% raw -> top component -> combined hard-eligible weighted GED
benchmark_methods = {'raw', 'ged_top_eig', 'ged_combined_hard_weighted'};
nBenchmarkMethods = numel(benchmark_methods);
raw_reference_definition = 'posterior_roi';  % locked reference for raw reference branch
compute_detectability = true;
compute_separation = true;
compute_reliability = true;

% Detrending parameters for power-ratio spectrum
% Power ratios often have 1/f-like or exponential trends; log-domain detrending
% typically yields flatter residuals and more reliable peak detection.
poly_order = 3;                    % higher order captures more complex broadband trends
detrend_edge_exclude_n = 5;        % exclude edge bins from fit (reduces edge artifacts)
detrend_in_log = true;             % true: fit in log(pr), subtract; flatter residuals for peak finding
detrend_flat_edges = true;         % true: constrain residual to zero at first/last frequency
peak_min_prom_frac = 0.15;         % MinPeakProminence as fraction of current-spectrum max
peak_min_distance_hz = 5;          % MinPeakDistance in Hz
sensitivity_poly_orders = unique(max(1, round(poly_order * sensitivity_scale)));
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
    fig_save_dir_ged = 'W:\Students\Arne\GCP\figures\eeg\ged';
    fig_save_dir_component_comparison = 'W:\Students\Arne\GCP\figures\eeg\ged\component_comparison';
    fig_save_dir_component_selection = 'W:\Students\Arne\GCP\figures\eeg\ged\component_selection';
    fig_save_dir_subj = 'W:\Students\Arne\GCP\figures\eeg\ged\subj';
else
    fig_save_dir_ged = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged';
    fig_save_dir_component_comparison = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged/component_comparison';
    fig_save_dir_component_selection = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged/component_selection';
    fig_save_dir_subj = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged/subj';
end
if ~exist(fig_save_dir_ged, 'dir'), mkdir(fig_save_dir_ged); end
if ~exist(fig_save_dir_component_comparison, 'dir'), mkdir(fig_save_dir_component_comparison); end
if ~exist(fig_save_dir_component_selection, 'dir'), mkdir(fig_save_dir_component_selection); end
if ~exist(fig_save_dir_subj, 'dir'), mkdir(fig_save_dir_subj); end
comp_sel_save_dir = fig_save_dir_component_selection;

%% Preallocate storage
all_trial_powratio     = cell(4, nSubj);
all_trial_powratio_early = cell(4, nSubj);
all_trial_powratio_late  = cell(4, nSubj);
all_trial_peaks_single = cell(4, nSubj);
all_trial_peaks_single_early = cell(4, nSubj);
all_trial_peaks_single_late  = cell(4, nSubj);
all_trial_centroid     = cell(4, nSubj);

all_trial_mean_single   = nan(4, nSubj);
all_trial_median_single = nan(4, nSubj);
all_trial_mean_single_early   = nan(4, nSubj);
all_trial_median_single_early = nan(4, nSubj);
all_trial_mean_single_late    = nan(4, nSubj);
all_trial_median_single_late  = nan(4, nSubj);
all_trial_mean_centroid   = nan(4, nSubj);
all_trial_median_centroid = nan(4, nSubj);

all_trial_detrate_single = nan(4, nSubj);
all_trial_detrate_single_early = nan(4, nSubj);
all_trial_detrate_single_late  = nan(4, nSubj);
all_trial_detrate_centroid = nan(4, nSubj);
all_trial_gamma_power = nan(4, nSubj);
all_trial_gamma_power_early = nan(4, nSubj);
all_trial_gamma_power_late  = nan(4, nSubj);

all_topos       = cell(1, nSubj);
all_topos_early = cell(1, nSubj);
all_topos_late  = cell(1, nSubj);
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
all_component_selection_stats = cell(1, nSubj); % full window; kept for backward compat
all_component_selection_stats_full  = cell(1, nSubj);
all_component_selection_stats_early = cell(1, nSubj);
all_component_selection_stats_late  = cell(1, nSubj);
all_threshold_inspection = cell(1, nSubj);
warning_log_by_subj = cell(nSubj, 1);

all_trial_powratio_bench = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_bench_early = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_bench_late  = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_dt_bench = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_dt_bench_early = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_dt_bench_late  = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_components = cell(4, nSubj); % full window [component x trial x frequency]; kept for backward compat
all_trial_powratio_components_full  = cell(4, nSubj);
all_trial_powratio_components_early = cell(4, nSubj);
all_trial_powratio_components_late  = cell(4, nSubj);
benchmark_metric_detectability = nan(nBenchmarkMethods, 4, nSubj);
benchmark_metric_prominence = nan(nBenchmarkMethods, 4, nSubj);
benchmark_metric_separation_slope = nan(nBenchmarkMethods, nSubj);
benchmark_metric_separation_delta = nan(nBenchmarkMethods, nSubj);
benchmark_metric_reliability_trialcv = nan(nBenchmarkMethods, 4, nSubj);
benchmark_metric_reliability_subjspread = nan(nBenchmarkMethods, 4);
sensitivity_results = table();
sensitivity_results_full  = table();
sensitivity_results_early = table();
sensitivity_results_late  = table();
primary_slope_stats = struct();
primary_delta_stats = struct();

first_subj_datapath = fullfile(path, subjects{1}, 'eeg');
tmp_chanlocs = load(fullfile(first_subj_datapath, 'dataEEG.mat'), 'dataEEG_c25');
chanlocs_all = tmp_chanlocs.dataEEG_c25.label;

%% Process each subject
if strcmpi(run_mode, 'component_check')
    fprintf('Run mode: component_check — GED + component selection only (Phase 1)\n');
end
for subj = 1:nSubj
    close all
    tic
    datapath = fullfile(path, subjects{subj}, 'eeg');
    eeg_data = load(fullfile(datapath, 'dataEEG.mat'), ...
        'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    dataEEG_c25 = eeg_data.dataEEG_c25;
    dataEEG_c50 = eeg_data.dataEEG_c50;
    dataEEG_c75 = eeg_data.dataEEG_c75;
    dataEEG_c100 = eeg_data.dataEEG_c100;
    warning_log_subj = struct('subject', {}, 'code', {}, 'message', {}, 'metrics', {});

    fsample = dataEEG_c25.fsample;

    trialIndices = { ...
        find(dataEEG_c25.trialinfo  == 61), ...
        find(dataEEG_c50.trialinfo  == 62), ...
        find(dataEEG_c75.trialinfo  == 63), ...
        find(dataEEG_c100.trialinfo == 64)};
    dataStructs = {dataEEG_c25, dataEEG_c50, dataEEG_c75, dataEEG_c100};

    nChans = length(dataEEG_c25.label);

    occ_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO)', 'once')), dataEEG_c25.label);
    occ_idx  = find(occ_mask);
    nOcc     = length(occ_idx);
    if nOcc == 0
        msg = sprintf('No occipital channels matched for subject %s. Using all channels as template fallback.', subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'NO_OCC_CHANNELS', msg, struct('nChans', nChans));
        occ_idx = 1:nChans;
        nOcc = nChans;
    end
    front_mask = cellfun(@(l) ~isempty(regexp(l, '^(Fp|AF|F)', 'once')), dataEEG_c25.label);
    front_idx  = find(front_mask);
    if isempty(front_idx)
        msg = sprintf('No frontal channels matched for subject %s. Frontal penalty disabled for this subject.', subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'NO_FRONT_CHANNELS', msg, struct('nChans', nChans));
    end
    post_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|P)', 'once')), dataEEG_c25.label);
    post_idx  = find(post_mask);
    if isempty(post_idx)
        post_idx = occ_idx;
    end
    raw_w = [];
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
    if isempty(raw_w)
        error('raw_w was not initialized for subject %s.', subjects{subj});
    end
    nonocc_idx = setdiff(1:nChans, occ_idx);

    %% ================================================================
    %  PHASE 1: Build POOLED covariance per window -> three GEDs
    %  ================================================================
    clc
    fprintf('Subject %s (%d/%d) — Phase 1: Per-window GED (full, early, late) (%d occ / %d ch)\n', ...
        subjects{subj}, subj, nSubj, nOcc, nChans);
    rng(random_seed + subj, 'twister');

    stim_windows = {full_window, early_window, late_window};
    win_names   = {'full', 'early', 'late'};
    lambdas     = [lambda_full, lambda_early, lambda_late];

    covStim_full  = zeros(nChans);
    covStim_early = zeros(nChans);
    covStim_late  = zeros(nChans);
    covBase_full  = zeros(nChans);
    nTrials_total = 0;

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

        cfg_t.latency = full_window;
        dat_stim_full = ft_selectdata(cfg_t, dat_gamma);
        cfg_t.latency = early_window;
        dat_stim_early = ft_selectdata(cfg_t, dat_gamma);
        cfg_t.latency = late_window;
        dat_stim_late = ft_selectdata(cfg_t, dat_gamma);

        nTrl = length(dat_stim_full.trial);
        for trl = 1:nTrl
            d_base = double(dat_base.trial{trl});
            d_base = bsxfun(@minus, d_base, mean(d_base, 2));
            cov_base_trl = (d_base * d_base') / size(d_base, 2);
            covBase_full = covBase_full + cov_base_trl;

            d = double(dat_stim_full.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            cov_stim_trl = (d * d') / size(d, 2);
            covStim_full = covStim_full + cov_stim_trl;

            d = double(dat_stim_early.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            covStim_early = covStim_early + (d * d') / size(d, 2);

            d = double(dat_stim_late.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            covStim_late = covStim_late + (d * d') / size(d, 2);
        end
        nTrials_total = nTrials_total + nTrl;
    end

    if nTrials_total < 1
        error('No valid trials available for subject %s after trial selection.', subjects{subj});
    end
    covStim_full  = covStim_full / nTrials_total;
    covStim_early = covStim_early / nTrials_total;
    covStim_late  = covStim_late / nTrials_total;
    covBase_full  = covBase_full / nTrials_total;

    % Covariances per window (for loop below)
    covStim_per_win = {covStim_full, covStim_early, covStim_late};
    % covBase_reg: regularized once, shared by all three GEDs (full, early, late)
    covBase_reg = (1-lambda)*covBase_full + lambda*mean(diag(covBase_full))*eye(nChans);

    % Run GED + component selection per window
    searchFilters_full  = [];
    searchFilters_early = [];
    searchFilters_late  = [];
    W_top_full  = []; W_top_early  = []; W_top_late  = [];
    W_combined_full  = []; W_combined_early  = []; W_combined_late  = [];
    selected_idx_full  = []; selected_idx_early  = []; selected_idx_late  = [];
    w_combined_full  = []; w_combined_early  = []; w_combined_late  = [];
    all_searchTopos  = cell(1, 3);
    all_evals_sorted = cell(1, 3);
    all_bestIdx      = cell(1, 3);
    all_topo_temp    = cell(1, 3);

    % Simulated signed occipital template (same for all windows)
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

    for w = 1:3
        covStim_w = covStim_per_win{w};
        lam_w = lambdas(w);
        covStim_reg = (1-lam_w)*covStim_w + lam_w*mean(diag(covStim_w))*eye(nChans);

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

        % Forward model for topoplot and component scoring (window-specific)
        for ci = 1:nSearch
            w_ci = W_full(:, ci);
            topo_ci = covStim_reg * w_ci;
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

    % Candidate metrics (composite score favors gamma + template match after hard gates).
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
    comp_score = score_w_corr * z_corr + score_w_gamma * z_gamma + score_w_eval * z_eval - score_w_frontleak * z_leak;
    finite_metrics = isfinite(corr_vec) & isfinite(ratio_vec) & isfinite(gamma_vec) & ...
        isfinite(evals_sorted(1:nSearch));
    hard_eligible_raw = finite_metrics & ...
        (evals_sorted(1:nSearch) >= min_eigval_hard) & ...
        (corr_vec >= min_corr_hard) & ...
        (ratio_vec >= min_occfront_ratio) & ...
        (gamma_vec >= min_gamma_hard);
        if w == 1
            all_threshold_inspection{subj} = compute_threshold_inspection( ...
        evals_sorted(1:nSearch), corr_vec, ratio_vec, gamma_vec, ...
        min_eigval_hard, min_corr_hard, min_occfront_ratio, min_gamma_hard, ...
        threshold_inspection_targets, nSearch);
        end
        no_hard_threshold_match = ~any(hard_eligible_raw);
    hard_eligible = hard_eligible_raw;
    dominant_outlier_mask = false(nSearch, 1);
    if exclude_dominant_outlier
        [hard_eligible, dominant_outlier_idx] = exclude_dominant_top_outlier( ...
            evals_sorted(1:nSearch), hard_eligible_raw, outlier_ratio_thr, outlier_mad_mult, outlier_min_rest);
        if ~isempty(dominant_outlier_idx)
            dominant_outlier_mask(dominant_outlier_idx) = true;
            msg = sprintf(['Dominant top GED component excluded for subject %s (component C%d) ', ...
                           'before combined-component selection.'], subjects{subj}, dominant_outlier_idx);
            ratio12_val = NaN;
            if numel(evals_sorted) >= 2 && isfinite(evals_sorted(2)) && evals_sorted(2) > 0
                ratio12_val = evals_sorted(1) / evals_sorted(2);
            end
            warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'DOMINANT_OUTLIER_EXCLUDED', msg, ...
                struct('component_idx', dominant_outlier_idx, 'lambda1', evals_sorted(1), ...
                       'lambda2', evals_sorted(min(2, numel(evals_sorted))), 'lambda1_lambda2_ratio', ratio12_val));
        end
    end
    searchScores = comp_score;
    searchScores(~hard_eligible) = -Inf;
    if ~any(isfinite(searchScores))
        msg = sprintf(['No components met hard thresholds for subject %s. ', ...
                       'Falling back to unconstrained eigenvalue ranking.'], subjects{subj});
        top_fail_idx = NaN;
        top_fail_eig = NaN;
        top_fail_corr = NaN;
        top_fail_ratio = NaN;
        top_fail_gamma = NaN;
        top_fail_pass_eig = false;
        top_fail_pass_corr = false;
        top_fail_pass_ratio = false;
        top_fail_pass_gamma = false;
        finite_idx = find(finite_metrics);
        if ~isempty(finite_idx)
            [~, fail_ord] = sort(evals_sorted(finite_idx), 'descend');
            top_fail_idx = finite_idx(fail_ord(1));
            top_fail_eig = evals_sorted(top_fail_idx);
            top_fail_corr = corr_vec(top_fail_idx);
            top_fail_ratio = ratio_vec(top_fail_idx);
            top_fail_gamma = gamma_vec(top_fail_idx);
            top_fail_pass_eig = top_fail_eig >= min_eigval_hard;
            top_fail_pass_corr = top_fail_corr >= min_corr_hard;
            top_fail_pass_ratio = top_fail_ratio >= min_occfront_ratio;
            top_fail_pass_gamma = top_fail_gamma >= min_gamma_hard;
        end
        hard_metrics = struct( ...
            'n_search', nSearch, ...
            'n_finite_metrics', sum(finite_metrics), ...
            'n_pass_eig', sum(finite_metrics & (evals_sorted(1:nSearch) >= min_eigval_hard)), ...
            'n_pass_corr', sum(finite_metrics & (corr_vec >= min_corr_hard)), ...
            'n_pass_ratio', sum(finite_metrics & (ratio_vec >= min_occfront_ratio)), ...
            'n_pass_gamma', sum(finite_metrics & (gamma_vec >= min_gamma_hard)), ...
            'n_pass_all_raw', sum(hard_eligible_raw), ...
            'n_excluded_dominant_outlier', sum(dominant_outlier_mask), ...
            'top_fail_idx', top_fail_idx, ...
            'top_fail_eig', top_fail_eig, ...
            'top_fail_corr', top_fail_corr, ...
            'top_fail_ratio', top_fail_ratio, ...
            'top_fail_gamma', top_fail_gamma, ...
            'top_fail_pass_eig', top_fail_pass_eig, ...
            'top_fail_pass_corr', top_fail_pass_corr, ...
            'top_fail_pass_ratio', top_fail_pass_ratio, ...
            'top_fail_pass_gamma', top_fail_pass_gamma, ...
            'thr_eig', min_eigval_hard, ...
            'thr_corr', min_corr_hard, ...
            'thr_ratio', min_occfront_ratio, ...
            'thr_gamma', min_gamma_hard);
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'NO_HARD_ELIGIBLE_COMPONENTS', msg, hard_metrics);
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
    candidate_table.hard_eligible_raw = hard_eligible_raw;
    candidate_table.dominant_outlier = dominant_outlier_mask;
    candidate_table.hard_eligible = hard_eligible;

    combined_idx = find(hard_eligible & isfinite(searchScores));
    if isempty(combined_idx)
        msg = sprintf(['No hard-eligible finite-score components available for subject %s. ', ...
            'Falling back to single best component (unconstrained ranking).'], subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'NO_HARD_COMPONENTS', msg, ...
            struct('n_hard_eligible', sum(hard_eligible), 'n_finite_scores', sum(isfinite(searchScores)), ...
            'fallback_idx', bestIdx));
        combined_idx = bestIdx;
    end
    [~, combined_ord] = sort(searchScores(combined_idx), 'descend');
    combined_idx = combined_idx(combined_ord);

    combined_weights = searchScores(combined_idx)';
    combined_weights(~isfinite(combined_weights) | combined_weights <= 0) = 0;
    if sum(combined_weights) <= 0
        combined_weights = evals_sorted(combined_idx)';
        combined_weights(~isfinite(combined_weights) | combined_weights <= 0) = 0;
    end
    if sum(combined_weights) <= 0
        combined_weights = ones(1, numel(combined_idx));
    end
    combined_weights = combined_weights / sum(combined_weights);

    selected_idx = combined_idx;
    selected_weights = combined_weights;

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
            topo_temp = covStim_reg * topComp;
        end

    finite_scores = find(isfinite(searchScores));
    if isempty(finite_scores)
        msg = sprintf('No finite eligible component scores for subject %s; using unconstrained eigenvalue ordering for display.', subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'NO_FINITE_SCORES_FOR_DISPLAY', msg, ...
            struct('n_hard_eligible', sum(hard_eligible), 'n_search', nSearch));
        [~, topDispOrder] = sort(evals_sorted(1:nSearch), 'descend');
    else
        [~, topDispOrder] = sort(searchScores, 'descend');
    end

    nStore = min(5, nSearch);
    storeCompIdx = topDispOrder(1:nStore);
    storeCorrs = searchCorrs(storeCompIdx);
    storeEvals = evals_sorted(storeCompIdx);
    storeTopos = searchTopos(:, storeCompIdx);
    storeCorrsFixed = nan(5, 1);
    storeEvalsFixed = nan(5, 1);
    storeCorrsFixed(1:nStore) = storeCorrs(:);
    storeEvalsFixed(1:nStore) = storeEvals(:);

        % Store per-window filters for Phase 2
        if w == 1
            searchFilters_full = searchFilters;
            W_top_full = searchFilters(:, 1);
            W_combined_full = searchFilters(:, selected_idx);
            selected_idx_full = selected_idx;
            w_combined_full = selected_weights(:)';
        elseif w == 2
            searchFilters_early = searchFilters;
            W_top_early = searchFilters(:, 1);
            W_combined_early = searchFilters(:, selected_idx);
            selected_idx_early = selected_idx;
            w_combined_early = selected_weights(:)';
        else
            searchFilters_late = searchFilters;
            W_top_late = searchFilters(:, 1);
            W_combined_late = searchFilters(:, selected_idx);
            selected_idx_late = selected_idx;
            w_combined_late = selected_weights(:)';
        end
        all_topo_temp{w} = topo_temp;

        if w == 1
            all_topos{subj}       = topo_temp;
            all_topo_labels{subj} = dataEEG_c25.label;
            all_eigenvalues(subj) = evals_sorted(bestIdx);
            all_selected_comp_idx(subj)  = bestIdx;
            all_selected_comp_corr(subj) = bestCorr;
            all_selected_comp_eval(subj) = evals_sorted(bestIdx);
            all_top5_corrs(:, subj) = storeCorrsFixed;
            all_top5_evals(:, subj) = storeEvalsFixed;
            all_top5_topos{subj} = storeTopos;
            all_simulated_templates{subj} = sim_template;
            all_selected_comp_indices_multi{subj} = selected_idx;
            all_selected_comp_weights{subj} = selected_weights(:)';
        elseif w == 2
            all_topos_early{subj} = topo_temp;
        else
            all_topos_late{subj} = topo_temp;
        end
        comp_sel_struct = struct( ...
            'selection_mode', 'hard_eligible_weighted', ...
            'selected_idx', selected_idx, ...
            'selected_weights', selected_weights, ...
            'candidate_table', candidate_table, ...
            'best_idx', bestIdx, ...
            'best_score', bestScore, ...
            'best_corr', bestCorr, ...
            'best_ratio', bestRatio, ...
            'best_gamma', bestGamma, ...
            'best_front', bestFront, ...
            'best_occ', bestOcc, ...
            'best_leak', bestLeak, ...
            'no_hard_threshold_match', no_hard_threshold_match);
        if w == 1
            all_component_selection_stats{subj} = comp_sel_struct;
            all_component_selection_stats_full{subj} = comp_sel_struct;
        elseif w == 2
            all_component_selection_stats_early{subj} = comp_sel_struct;
        else
            all_component_selection_stats_late{subj} = comp_sel_struct;
        end

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

    % Combined-component topoplots: all hard-eligible selected components (ordered by composite score)
    [~, score_ord] = sort(searchScores(selected_idx), 'descend');
    selTopoIdx = selected_idx(score_ord)';
    selWeightsOrdered = selected_weights(score_ord);
    nSelTopo = numel(selTopoIdx);
    nColsSel = 5;
    nRowsSel = ceil(max(nSelTopo, 1) / nColsSel);
    % Scale figure height by number of rows: single-row figures use ~half height to reduce excess whitespace
    fig_height = min(982, round(491 * nRowsSel));  % 491 ≈ 982/2 per row, capped at 982
    fig_post = figure('Position', [0 0 1512 fig_height], 'Color', 'w');
    title_post = sprintf('Subject %s: Eligible Components (n=%d)', ...
        subjects{subj}, nSelTopo);
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
            if ~isempty(selWeightsOrdered) && numel(selWeightsOrdered) >= si
                w_show = selWeightsOrdered(si);
            end
            score_show = searchScores(comp_rank);
            g_log = searchGammaEvidence(comp_rank);
            g_pct = 100 * (exp(g_log) - 1);   % gamma amplification as % change from baseline
            title(sprintf('C%d:\\lambda=%.2f,score=%.2f,w=%.3f,r=%.2f,ratio=%.2f,g=%.0f%%', ...
                comp_rank, evals_sorted(comp_rank), score_show, w_show, searchCorrs(comp_rank), ...
                searchOccFrontRatio(comp_rank), g_pct), 'FontSize', 6);
        end
    end
        saveas(fig_post, fullfile(comp_sel_save_dir, sprintf('GCP_eeg_GED_component_selection_subj%s_%s.png', subjects{subj}, win_names{w})));
    end

    if strcmpi(run_mode, 'component_check')
        warning_log_by_subj{subj} = warning_log_subj;
        toc
        continue
    end

    %% ================================================================
    %  PHASE 2: Per condition — trial-level narrowband scanning
    %  ================================================================
    % Use window-specific filters for each power-ratio output
    for wi = 1:3
        if wi == 1
            W_comb = W_combined_full;
            w_comb = w_combined_full;
            W_t = W_top_full;
            sel_idx = selected_idx_full;
        elseif wi == 2
            W_comb = W_combined_early;
            w_comb = w_combined_early;
            W_t = W_top_early;
            sel_idx = selected_idx_early;
        else
            W_comb = W_combined_late;
            w_comb = w_combined_late;
            W_t = W_top_late;
            sel_idx = selected_idx_late;
        end
        if isempty(W_comb)
            msg = sprintf('No selected combined filters for subject %s (%s window).', subjects{subj}, win_names{wi});
            warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'EMPTY_W_COMBINED', msg, struct());
        end
        if isempty(W_t)
            W_t = zeros(nChans, 0);
        end
        if isempty(w_comb) && ~isempty(W_comb)
            w_comb = ones(1, size(W_comb, 2)) / size(W_comb, 2);
        end
        if sum(w_comb) <= 0 && ~isempty(W_comb)
            w_comb = ones(1, size(W_comb, 2)) / size(W_comb, 2);
        end
        if ~isempty(W_comb)
            w_comb = w_comb(:)' / sum(w_comb);
        end
        if wi == 1
            W_combined_full_norm = W_comb;
            w_combined_full_norm = w_comb;
            W_top_full_norm = W_t;
            selected_idx_full_norm = sel_idx;
        elseif wi == 2
            W_combined_early_norm = W_comb;
            w_combined_early_norm = w_comb;
            W_top_early_norm = W_t;
            selected_idx_early_norm = sel_idx;
        else
            W_combined_late_norm = W_comb;
            w_combined_late_norm = w_comb;
            W_top_late_norm = W_t;
            selected_idx_late_norm = sel_idx;
        end
    end
    % Require at least full-window filters
    if isempty(W_combined_full)
        msg = sprintf('No selected combined filters for subject %s (full window).', subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'EMPTY_W_COMBINED', msg, struct());
        error(msg);
    end

    % Per-window filters struct for Phase 2
    filters = struct('full', struct(), 'early', struct(), 'late', struct());
    filters.full.searchFilters = searchFilters_full;
    filters.full.W_top = W_top_full_norm;
    filters.full.W_combined = W_combined_full_norm;
    filters.full.selected_idx = selected_idx_full_norm;
    filters.full.w_combined = w_combined_full_norm;
    filters.early.searchFilters = searchFilters_early;
    filters.early.W_top = W_top_early_norm;
    filters.early.W_combined = W_combined_early_norm;
    filters.early.selected_idx = selected_idx_early_norm;
    filters.early.w_combined = w_combined_early_norm;
    filters.late.searchFilters = searchFilters_late;
    filters.late.W_top = W_top_late_norm;
    filters.late.W_combined = W_combined_late_norm;
    filters.late.selected_idx = selected_idx_late_norm;
    filters.late.w_combined = w_combined_late_norm;

    for cond = 1:4

        dat = dat_per_cond{cond};
        if isempty(dat), continue; end

        nTrl = length(dat.trial);
        powratio_methods_full = nan(nBenchmarkMethods, nTrl, nFreqs);
        powratio_methods_early = nan(nBenchmarkMethods, nTrl, nFreqs);
        powratio_methods_late = nan(nBenchmarkMethods, nTrl, nFreqs);
        nSearch_full = size(filters.full.searchFilters, 2);
        nSearch_early = size(filters.early.searchFilters, 2);
        nSearch_late = size(filters.late.searchFilters, 2);
        powratio_components       = nan(nSearch_full, nTrl, nFreqs);
        powratio_components_early = nan(nSearch_early, nTrl, nFreqs);
        powratio_components_late  = nan(nSearch_late, nTrl, nFreqs);

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

            for trl = 1:nTrl
                x_nb = double(dat_nb.trial{trl});
                t_nb = dat_nb.time{trl};
                idx_base = t_nb >= baseline_window(1) & t_nb <= baseline_window(2);
                idx_full = t_nb >= full_window(1) & t_nb <= full_window(2);
                idx_early = t_nb >= early_window(1) & t_nb <= early_window(2);
                idx_late = t_nb >= late_window(1) & t_nb <= late_window(2);

                x_base = x_nb(:, idx_base);
                x_full = x_nb(:, idx_full);
                x_early = x_nb(:, idx_early);
                x_late = x_nb(:, idx_late);

                nSearch_full = size(filters.full.searchFilters, 2);
                comp_base_all_full = filters.full.searchFilters(:, 1:nSearch_full)' * x_base;
                pow_base_all_full = mean(comp_base_all_full.^2, 2);
                ratio_all_full = nan(nSearch_full, 1);
                comp_stim_all_full = filters.full.searchFilters(:, 1:nSearch_full)' * x_full;
                pow_stim_all_full = mean(comp_stim_all_full.^2, 2);
                valid_all_full = isfinite(pow_stim_all_full) & isfinite(pow_base_all_full) & (pow_base_all_full > 0);
                ratio_all_full(valid_all_full) = pow_stim_all_full(valid_all_full) ./ pow_base_all_full(valid_all_full);
                powratio_components(:, trl, fi) = ratio_all_full;
                powratio_methods_full(:, trl, fi) = compute_method_ratios_from_components( ...
                    ratio_all_full, x_full, x_base, raw_w, filters.full.W_top, filters.full.W_combined, filters.full.selected_idx, filters.full.w_combined, nBenchmarkMethods);

                nSearch_early = size(filters.early.searchFilters, 2);
                comp_base_all_early = filters.early.searchFilters(:, 1:nSearch_early)' * x_base;
                pow_base_all_early = mean(comp_base_all_early.^2, 2);
                ratio_all_early = nan(nSearch_early, 1);
                comp_stim_all_early = filters.early.searchFilters(:, 1:nSearch_early)' * x_early;
                pow_stim_all_early = mean(comp_stim_all_early.^2, 2);
                valid_all_early = isfinite(pow_stim_all_early) & isfinite(pow_base_all_early) & (pow_base_all_early > 0);
                ratio_all_early(valid_all_early) = pow_stim_all_early(valid_all_early) ./ pow_base_all_early(valid_all_early);
                powratio_components_early(1:nSearch_early, trl, fi) = ratio_all_early;
                powratio_methods_early(:, trl, fi) = compute_method_ratios_from_components( ...
                    ratio_all_early, x_early, x_base, raw_w, filters.early.W_top, filters.early.W_combined, filters.early.selected_idx, filters.early.w_combined, nBenchmarkMethods);

                nSearch_late = size(filters.late.searchFilters, 2);
                comp_base_all_late = filters.late.searchFilters(:, 1:nSearch_late)' * x_base;
                pow_base_all_late = mean(comp_base_all_late.^2, 2);
                ratio_all_late = nan(nSearch_late, 1);
                comp_stim_all_late = filters.late.searchFilters(:, 1:nSearch_late)' * x_late;
                pow_stim_all_late = mean(comp_stim_all_late.^2, 2);
                valid_all_late = isfinite(pow_stim_all_late) & isfinite(pow_base_all_late) & (pow_base_all_late > 0);
                ratio_all_late(valid_all_late) = pow_stim_all_late(valid_all_late) ./ pow_base_all_late(valid_all_late);
                powratio_components_late(1:nSearch_late, trl, fi) = ratio_all_late;
                powratio_methods_late(:, trl, fi) = compute_method_ratios_from_components( ...
                    ratio_all_late, x_late, x_base, raw_w, filters.late.W_top, filters.late.W_combined, filters.late.selected_idx, filters.late.w_combined, nBenchmarkMethods);
            end
            clear dat_nb
        end
        all_trial_powratio_components{cond, subj} = powratio_components;          % full; backward compat
        all_trial_powratio_components_full{cond, subj}  = powratio_components;
        all_trial_powratio_components_early{cond, subj} = powratio_components_early;
        all_trial_powratio_components_late{cond, subj}  = powratio_components_late;

        for mi = 1:nBenchmarkMethods
            pr_m_full = squeeze(powratio_methods_full(mi, :, :));
            all_trial_powratio_bench{mi, cond, subj} = pr_m_full;
            if ~isempty(pr_m_full)
                dt_m_full = nan(size(pr_m_full));
                for trl = 1:size(pr_m_full, 1)
                    dt_m_full(trl,:) = detrend_power_ratio(pr_m_full(trl,:), scan_freqs, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
                end
                all_trial_powratio_dt_bench{mi, cond, subj} = dt_m_full;
            end

            pr_m_early = squeeze(powratio_methods_early(mi, :, :));
            all_trial_powratio_bench_early{mi, cond, subj} = pr_m_early;
            if ~isempty(pr_m_early)
                dt_m_early = nan(size(pr_m_early));
                for trl = 1:size(pr_m_early, 1)
                    dt_m_early(trl,:) = detrend_power_ratio(pr_m_early(trl,:), scan_freqs, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
                end
                all_trial_powratio_dt_bench_early{mi, cond, subj} = dt_m_early;
            end

            pr_m_late = squeeze(powratio_methods_late(mi, :, :));
            all_trial_powratio_bench_late{mi, cond, subj} = pr_m_late;
            if ~isempty(pr_m_late)
                dt_m_late = nan(size(pr_m_late));
                for trl = 1:size(pr_m_late, 1)
                    dt_m_late(trl,:) = detrend_power_ratio(pr_m_late(trl,:), scan_freqs, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
                end
                all_trial_powratio_dt_bench_late{mi, cond, subj} = dt_m_late;
            end
        end

        % Keep full-window outputs based on weighted combined GED branch.
        powratio_trials_full = squeeze(powratio_methods_full(3, :, :));
        powratio_trials_early = squeeze(powratio_methods_early(3, :, :));
        powratio_trials_late = squeeze(powratio_methods_late(3, :, :));
        all_trial_powratio{cond, subj} = powratio_trials_full;
        all_trial_powratio_early{cond, subj} = powratio_trials_early;
        all_trial_powratio_late{cond, subj} = powratio_trials_late;
        if ~isempty(powratio_trials_full)
            trl_gamma_power_full = mean(powratio_trials_full, 2, 'omitnan');
            all_trial_gamma_power(cond, subj) = mean(trl_gamma_power_full, 'omitnan');
        end
        if ~isempty(powratio_trials_early)
            trl_gamma_power_early = mean(powratio_trials_early, 2, 'omitnan');
            all_trial_gamma_power_early(cond, subj) = mean(trl_gamma_power_early, 'omitnan');
        end
        if ~isempty(powratio_trials_late)
            trl_gamma_power_late = mean(powratio_trials_late, 2, 'omitnan');
            all_trial_gamma_power_late(cond, subj) = mean(trl_gamma_power_late, 'omitnan');
        end

        %% Per-trial peak detection
        trl_peaks_single = nan(nTrl, 1);
        trl_centroid     = nan(nTrl, 1);

        for trl = 1:nTrl
            pr = powratio_trials_full(trl, :);
            if all(isnan(pr)), continue; end

            pr_dt = detrend_power_ratio(pr, scan_freqs, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
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
            mprom = max(0, max(pr_dt_smooth) * peak_min_prom_frac);
            [pks, locs] = findpeaks(pr_dt_smooth, scan_freqs, ...
                'MinPeakProminence', mprom, ...
                'MinPeakDistance', peak_min_distance_hz);

            if ~isempty(pks)
                [~, best_pk] = max(pks);
                trl_peaks_single(trl) = locs(best_pk);
            end
        end

        all_trial_peaks_single{cond, subj} = trl_peaks_single;
        all_trial_centroid{cond, subj}     = trl_centroid;

        valid_s = ~isnan(trl_peaks_single);
        all_trial_mean_single(cond, subj)   = mean(trl_peaks_single(valid_s));
        all_trial_median_single(cond, subj) = median(trl_peaks_single(valid_s));
        all_trial_detrate_single(cond, subj) = sum(valid_s) / nTrl;

        valid_c = ~isnan(trl_centroid);
        all_trial_mean_centroid(cond, subj)   = mean(trl_centroid(valid_c));
        all_trial_median_centroid(cond, subj) = median(trl_centroid(valid_c));
        all_trial_detrate_centroid(cond, subj) = sum(valid_c) / nTrl;

        % Time-split single-peak summaries.
        trl_peaks_single_early = detect_single_peaks_from_powratio( ...
            powratio_trials_early, scan_freqs, poly_order, detrend_edge_exclude_n, ...
            detrend_in_log, detrend_flat_edges, peak_min_prom_frac, peak_min_distance_hz);
        trl_peaks_single_late = detect_single_peaks_from_powratio( ...
            powratio_trials_late, scan_freqs, poly_order, detrend_edge_exclude_n, ...
            detrend_in_log, detrend_flat_edges, peak_min_prom_frac, peak_min_distance_hz);
        all_trial_peaks_single_early{cond, subj} = trl_peaks_single_early;
        all_trial_peaks_single_late{cond, subj} = trl_peaks_single_late;

        valid_s_early = ~isnan(trl_peaks_single_early);
        all_trial_mean_single_early(cond, subj) = mean(trl_peaks_single_early(valid_s_early));
        all_trial_median_single_early(cond, subj) = median(trl_peaks_single_early(valid_s_early));
        all_trial_detrate_single_early(cond, subj) = sum(valid_s_early) / nTrl;

        valid_s_late = ~isnan(trl_peaks_single_late);
        all_trial_mean_single_late(cond, subj) = mean(trl_peaks_single_late(valid_s_late));
        all_trial_median_single_late(cond, subj) = median(trl_peaks_single_late(valid_s_late));
        all_trial_detrate_single_late(cond, subj) = sum(valid_s_late) / nTrl;

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

    % Pre-compute detrended matrices (3 rows: heatmap, single-peak spectra, topoplot+histogram)
    pr_dt_mats = cell(1, 4);
    for cond = 1:4
        pr_mat = all_trial_powratio{cond, subj};
        if ~isempty(pr_mat)
            nTrl = size(pr_mat, 1);
            dt = nan(size(pr_mat));
            for trl = 1:nTrl
                dt(trl,:) = detrend_power_ratio(pr_mat(trl,:), scan_freqs, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
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
        subplot(3, 4, cond);
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
        subplot(3, 4, 4 + cond); hold on;
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

    % --- Row 3: Topoplot + histogram ---
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

    subplot(3, 4, 9);
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

    subplot(3, 4, [10 11 12]); hold on;
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

    saveas(fig, fullfile(fig_save_dir_subj, sprintf('GCP_eeg_GED_subj%s.png', subjects{subj})));
    warning_log_by_subj{subj} = warning_log_subj;
    toc
end % subject loop

warning_log = vertcat(warning_log_by_subj{:});
print_subject_warning_summary(warning_log);
print_threshold_inspection_summary(subjects, all_threshold_inspection, threshold_inspection_targets);

if strcmpi(run_mode, 'component_check')
    clc
    fprintf('Component check complete: GED components generated for all subjects.\n');
else
%% ====================================================================
%  THREE-WAY BENCHMARK METRICS (raw vs top-eig GED vs combined GED)
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
                mprom = max(0, max(y) * peak_min_prom_frac);
                [pks, locs, ~, p] = findpeaks(y, scan_freqs, ...
                    'MinPeakDistance', peak_min_distance_hz, ...
                    'MinPeakProminence', mprom);
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

% Shared y-limits for component selection lower-row panels (consistent across subjects).
single_vals = all_trial_median_single(isfinite(all_trial_median_single));
if isempty(single_vals)
    bench_single_ylim = [30 90];
else
    single_q = prctile(single_vals, [2 98]);
    bench_single_ylim = [max(30, single_q(1) - 1), min(90, single_q(2) + 1)];
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

%% Subject-level component comparison figures (raw, top GED, combined GED; edge flattening applied)
bench_method_labels = {'Raw', 'Top-Eig GED', 'Combined GED'};
bench_method_colors = [0.2 0.2 0.2; 0.1 0.35 0.75; 0.85 0.55 0.10];
for subj = 1:nSubj
    fig_bench_subj = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('Component Comparison: Subject %s', subjects{subj}), ...
        'FontSize', 18, 'FontWeight', 'bold');

    for mi = 1:nBenchmarkMethods
        subplot(2, 3, mi); hold on;
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

    subplot(2, 3, 4); hold on;
    single_med = nan(1, 4);
    single_lo = nan(1, 4);
    single_hi = nan(1, 4);
    for cond = 1:4
        sp = all_trial_peaks_single{cond, subj};
        sp = sp(~isnan(sp));
        if ~isempty(sp)
            single_med(cond) = median(sp);
            sq = prctile(sp, [25 75]);
            single_lo(cond) = single_med(cond) - sq(1);
            single_hi(cond) = sq(2) - single_med(cond);
        end
    end
    b4 = bar(1:4, single_med, 0.58, 'FaceColor', 'flat', 'EdgeColor', 'none');
    b4.CData = colors;
    errorbar(1:4, single_med, single_lo, single_hi, 'k', 'LineStyle', 'none', 'LineWidth', 1.1, 'CapSize', 5);
    ylabel('Single-peak median [Hz]');
    ylim(bench_single_ylim);
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 10);
    d_single = single_med(4) - single_med(1);
    text(0.02, 0.96, sprintf('\\DeltaSingle_{100-25}=%.2f Hz', d_single), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
    title('Single median (combined GED)');
    box on;

    subplot(2, 3, 5); hold on;
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
    text(0.02, 0.96, sprintf('\\DeltaProm_{Comb-R}=%.2f\n\\DeltaCV_{Comb-R}=%.3f', ...
        prom_vec(end)-prom_vec(1), rel_vec(end)-rel_vec(1)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
    title('Prominence + Reliability');
     box on;

    subplot(2, 3, 6); hold on;
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
    text(0.02, 0.96, sprintf('Combined slope=%.2f\nCombined \\Delta=%.2f Hz', ...
        slope_vec(end), delta_vec(end)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
    title('Condition separation');
     box on;

    saveas(fig_bench_subj, fullfile(fig_save_dir_component_comparison, sprintf('GCP_eeg_GED_component_comparison_subj%s.png', subjects{subj})));
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
saveas(fig_cent, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_centroid_summary.png'));

%% ====================================================================
%  GRAND-AVERAGE COMPONENT COMPARISON
%  ====================================================================
fig_bench_group = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Component Comparison: Grand Average', 'FontSize', 18, 'FontWeight', 'bold');

for mi = 1:nBenchmarkMethods
    subplot(2, 3, mi); hold on;
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

subplot(2, 3, 4); hold on;
single_mu = nanmean(all_trial_median_single, 2);
single_se = nanstd(all_trial_median_single, [], 2) ./ sqrt(sum(~isnan(all_trial_median_single), 2));
single_med = nanmedian(all_trial_median_single, 2);
b4g = bar(1:4, single_mu, 0.58, 'FaceColor', 'flat', 'EdgeColor', 'none');
b4g.CData = colors;
errorbar(1:4, single_mu, single_se, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);
scatter(1:4, single_med, 35, 'kd', 'filled');
ylabel('Single-peak mean [Hz]');
ylim(bench_single_ylim);
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 10);
text(0.02, 0.96, sprintf('\\DeltaSingle_{100-25}=%.2f Hz', single_mu(4)-single_mu(1)), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
title('Single mean (combined GED)');
box on;

subplot(2, 3, 5); hold on;
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
text(0.02, 0.96, sprintf('\\DeltaProm_{Comb-R}=%.2f\n\\DeltaCV_{Comb-R}=%.3f', ...
    prom_mu(end)-prom_mu(1), rel_mu(end)-rel_mu(1)), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
title('Prominence + Reliability');
 box on;

subplot(2, 3, 6); hold on;
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
text(0.02, 0.96, sprintf('Combined slope=%.2f\nCombined \\Delta=%.2f Hz', ...
    slope_mu(end), delta_mu(end)), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8);
title('Condition separation');
 box on;

saveas(fig_bench_group, fullfile(fig_save_dir_component_comparison, 'GCP_eeg_GED_component_comparison_grandaverage.png'));

%% ====================================================================
%  STANDALONE CONDITION-SEPARATION METRICS (combined GED)
%  ====================================================================
close all
fig_cond_slope = figure('Position', [0 0 1512 982], 'Color', 'w');

post_idx = 3; % Combined GED
slope_post = benchmark_metric_separation_slope(post_idx, :);
delta_post = benchmark_metric_separation_delta(post_idx, :);
primary_slope_stats = compute_one_sample_stats(slope_post);
primary_delta_stats = compute_one_sample_stats(delta_post);

fprintf('\n============================================================\n');
fprintf('Primary inference (combined GED, subject-level)\n');
fprintf('============================================================\n');
fprintf(['Slope: n=%d, mean=%.3f Hz/condition, t(%d)=%.3f, p=%.4f, ', ...
         '95%% CI [%.3f, %.3f], d=%.3f\n'], ...
    primary_slope_stats.n, primary_slope_stats.mean, primary_slope_stats.df, ...
    primary_slope_stats.tstat, primary_slope_stats.p, ...
    primary_slope_stats.ci_low, primary_slope_stats.ci_high, primary_slope_stats.cohens_d);
fprintf(['Delta (100%%-25%%): n=%d, mean=%.3f Hz, t(%d)=%.3f, p=%.4f, ', ...
         '95%% CI [%.3f, %.3f], d=%.3f\n'], ...
    primary_delta_stats.n, primary_delta_stats.mean, primary_delta_stats.df, ...
    primary_delta_stats.tstat, primary_delta_stats.p, ...
    primary_delta_stats.ci_low, primary_delta_stats.ci_high, primary_delta_stats.cohens_d);

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
set(gca, 'XTick', 1, 'XTickLabel', {'Combined GED'}, ...
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
set(gca, 'XTick', 1, 'XTickLabel', {'Combined GED'}, ...
    'FontSize', 16, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
ylabel('\Delta median (100% - 25%) [Hz]', 'FontSize', 18, 'FontWeight', 'bold');
title('Median Frequency Shift (100% - 25%)', 'FontSize', 20, 'FontWeight', 'bold');

saveas(fig_cond_slope, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_condition_slope.png'));

%% Mean gamma frequency shift bar plot
close all
fig_cond_shift_bar = figure('Position', [0 0 1512/2 982], 'Color', 'w');
valid_delta = isfinite(delta_post);
subj_idx = find(valid_delta);
delta_vals = delta_post(valid_delta);
bar(subj_idx, delta_vals, 0.75, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'k');
hold on;
yline(0, 'k--', 'LineWidth', 1.2);
xlim([0.5 max(subj_idx) + 0.5]);
ylim([-6 6]);
xticks(subj_idx);
if exist('subjects', 'var') == 1 && numel(subjects) >= max(subj_idx) && numel(subj_idx) <= 35
    xticklabels(subjects(subj_idx));
end
xlabel('Subject', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('\Delta mean frequency shift (100% - 25%) [Hz]', 'FontSize', 18, 'FontWeight', 'bold');
title('Gamma Frequency Shift', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
saveas(fig_cond_shift_bar, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_bar_GammaFreq.png'));

%% ====================================================================
%  SUMMARY DASHBOARD (backprojected combined-component data)
%  ====================================================================
fig_summary = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('GED Trials Summary Dashboard (component selection backprojected)', ...
    'FontSize', 16, 'FontWeight', 'bold');

summary_metrics = { ...
    all_trial_median_single, ...
    squeeze(benchmark_metric_prominence(3, :, :)), ...
    squeeze(benchmark_metric_reliability_trialcv(3, :, :)), ...
    all_trial_median_centroid, ...
    all_trial_gamma_power};
summary_names = {'Single median [Hz]', 'Prominence', 'Trial CV', ...
    'Centroid median [Hz]', 'Gamma Power over Conditions'};

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
        ylim([0 10]);
    elseif mi == 4
        ylim([50 60]);
    elseif mi == 5
        ylim([0 15]);
    end
end
saveas(fig_summary, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_metrics_summary.png'));

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

saveas(fig_box1, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_boxplot_GammaFreq.png'));

%% ====================================================================
%  GRAND AVERAGE: Single Peak with subject trajectories + ID labels
%  ====================================================================
fig_box1_traj = figure('Position', [0 0 1512/2 982], 'Color', 'w');
hold on;

dat = all_trial_median_single;
mu = nanmean(dat, 2);
sem = nanstd(dat, [], 2) ./ sqrt(sum(~isnan(dat), 2));
med = nanmedian(dat, 2);

for c = 1:4
    bar(c, mu(c), 0.6, 'FaceColor', colors(c,:), 'EdgeColor', 'k', 'FaceAlpha', 0.35);
end
errorbar(1:4, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.4, 'CapSize', 6);
scatter(1:4, med, 70, 'kd', 'filled');

for s = 1:nSubj
    subj_offset = (rand - 0.5) * 0.22;
    x_subj = (1:4) + subj_offset;
    y_subj = dat(:, s)';
    valid_subj = ~isnan(y_subj);

    if sum(valid_subj) >= 2
        plot(x_subj(valid_subj), y_subj(valid_subj), '-', ...
            'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
    end

    for c = 1:4
        if ~isnan(dat(c, s))
            scatter(x_subj(c), dat(c, s), 120, [0.45 0.45 0.45], ...
                'filled', ...
                'MarkerFaceAlpha', 0.55, ...
                'MarkerEdgeColor', [1 1 1], ...
                'LineWidth', 0.5, ...
                'MarkerEdgeAlpha', 0.9);
        end
    end

    if ~isnan(dat(4, s))
        subj_id_txt = string(subjects{s});
        text(x_subj(4) + 0.10, dat(4, s), subj_id_txt, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 9, ...
            'Color', [0.15 0.15 0.15], ...
            'BackgroundColor', [1 1 1], ...
            'Margin', 1);
    end
end

xlim([0.3 5.1]);
ylim([45 75]);
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 16, 'Box', 'off');
ylabel('Gamma Frequency [Hz]');
title('Gamma Frequency over Conditions (Subject Trajectories)', ...
    'FontSize', 18, 'FontWeight', 'bold');

saveas(fig_box1_traj, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_boxplot_GammaFreq_IDs.png'));

%% ====================================================================
%  GRAND AVERAGE: Single Peak (stats-style boxplot, non-baselined freq)
%  ====================================================================
fig_box1_statsstyle = figure('Position', [0 0 1512/2 982], 'Color', 'W');
hold on;

dat = all_trial_median_single;  % [condition x subject], absolute frequency (Hz)

% Subject-wise lines across contrast conditions
for s = 1:nSubj
    y_subj = dat(:, s);
    valid_subj = ~isnan(y_subj);
    if sum(valid_subj) >= 2
        plot(find(valid_subj), y_subj(valid_subj), '-', ...
            'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    end
end

% Boxplot for each condition
y_all = dat(:);
g_all = repmat((1:4)', nSubj, 1);
valid_all = ~isnan(y_all);
boxplot(y_all(valid_all), g_all(valid_all), 'Colors', [0.45 0.45 0.45], ...
    'Symbol', '', 'Widths', 0.5);

% Overlay jittered subject points
for c = 1:4
    y_c = dat(c, :);
    valid_c = ~isnan(y_c);
    x_jit = c + (rand(1, sum(valid_c)) - 0.5) * 0.10;
    scatter(x_jit, y_c(valid_c), 170, colors(c,:), 'filled', ...
        'MarkerEdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.7);
end

set(gca, 'XTick', 1:4, 'XTickLabel', strcat(condLabels, ' Contrast'), ...
    'FontSize', 15, 'Box', 'off');
xlim([0.5 4.5]);
%ylim([30 90]);
ylabel('Gamma Frequency [Hz]');
title('Gamma Frequency', 'FontSize', 30, 'FontWeight', 'bold');

saveas(fig_box1_statsstyle, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_boxplot_GammaFreq_statsStyle.png'));

%% ====================================================================
%  MAIN FIGURE: Gamma frequency over contrast (mean+-SEM + trajectories)
%  ====================================================================
fig_main_gamma = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

dat = all_trial_median_single;  % [condition x subject]
mu = nanmean(dat, 2);
sem = nanstd(dat, [], 2) ./ sqrt(sum(~isnan(dat), 2));

for s = 1:nSubj
    y_subj = dat(:, s);
    valid_subj = ~isnan(y_subj);
    if sum(valid_subj) >= 2
        plot(find(valid_subj), y_subj(valid_subj), '-', ...
            'Color', [0.78 0.78 0.78], 'LineWidth', 1.1);
    end
end

for c = 1:4
    y_c = dat(c, :);
    valid_c = ~isnan(y_c);
    x_jit = c + (rand(1, sum(valid_c)) - 0.5) * 0.12;
    scatter(x_jit, y_c(valid_c), 130, colors(c,:), 'filled', ...
        'MarkerEdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.7);
end

errorbar(1:4, mu, sem, 'k-o', 'LineWidth', 2.2, 'CapSize', 9, ...
    'MarkerFaceColor', [0.1 0.1 0.1], 'MarkerSize', 7);
set(gca, 'XTick', 1:4, 'XTickLabel', strcat(condLabels, ' Contrast'), ...
    'FontSize', 16, 'Box', 'off');
xlim([0.5 4.5]);
ylabel('Gamma Frequency [Hz]');
title('Gamma Frequency over Contrast', 'FontSize', 24, 'FontWeight', 'bold');
text(0.02, 0.98, sprintf('Linear trend (subject slope): p=%.4f, d=%.2f', ...
    primary_slope_stats.p, primary_slope_stats.cohens_d), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig_main_gamma, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_main_GammaFreq.png'));

%% ====================================================================
%  MAIN FIGURE: Gamma frequency over contrast by time window
%  ====================================================================
[gamma_slope_full, gamma_delta_full] = compute_condition_separation_from_matrix(all_trial_median_single);
[gamma_slope_early, gamma_delta_early] = compute_condition_separation_from_matrix(all_trial_median_single_early);
[gamma_slope_late, gamma_delta_late] = compute_condition_separation_from_matrix(all_trial_median_single_late);
fprintf('\nTime-split gamma separation (descriptive):\n');
fprintf('Full  slope=%.3f Hz/condition, delta(100-25)=%.3f Hz\n', ...
    mean(gamma_slope_full, 'omitnan'), mean(gamma_delta_full, 'omitnan'));
fprintf('Early slope=%.3f Hz/condition, delta(100-25)=%.3f Hz\n', ...
    mean(gamma_slope_early, 'omitnan'), mean(gamma_delta_early, 'omitnan'));
fprintf('Late  slope=%.3f Hz/condition, delta(100-25)=%.3f Hz\n', ...
    mean(gamma_slope_late, 'omitnan'), mean(gamma_delta_late, 'omitnan'));

fig_main_gamma_windows = figure('Position', [0 0 1512 982], 'Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot_gamma_window_panel(all_trial_median_single_early, condLabels, colors, nSubj);
ylabel('Gamma Frequency [Hz]');
title('Early (0-600 ms)', 'FontWeight', 'bold');

nexttile; hold on;
plot_gamma_window_panel(all_trial_median_single_late, condLabels, colors, nSubj);
title('Late (1000-2000 ms)', 'FontWeight', 'bold');

nexttile; hold on;
late_minus_early_subj = nanmean(all_trial_median_single_late - all_trial_median_single_early, 1);
valid_shift = isfinite(late_minus_early_subj);
if any(valid_shift)
    boxplot(late_minus_early_subj(valid_shift)', ones(sum(valid_shift), 1), ...
        'Colors', 'k', 'Symbol', '', 'Widths', 0.4);
    xj = 1 + (rand(1, sum(valid_shift)) - 0.5) * 0.20;
    scatter(xj, late_minus_early_subj(valid_shift), 70, [0.45 0.45 0.45], 'filled', ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.7);
    yline(0, 'k--', 'LineWidth', 1);
end
xlim([0.5 1.5]);
set(gca, 'XTick', 1, 'XTickLabel', {'Late-Early'}, 'Box', 'off', 'FontSize', 12);
ylabel('\Delta Gamma Frequency [Hz]');
title('Subject Shift', 'FontWeight', 'bold');
saveas(fig_main_gamma_windows, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_main_GammaFreq_timeSplit.png'));

%% ====================================================================
%  POWER INCREASE FIGURE: gamma-band power ratio over conditions
%  ====================================================================
fig_power_statsstyle = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

dat_power = all_trial_gamma_power; % mean trial-level stim/base gamma power ratio

for s = 1:nSubj
    y_subj = dat_power(:, s);
    valid_subj = ~isnan(y_subj);
    if sum(valid_subj) >= 2
        plot(find(valid_subj), y_subj(valid_subj), '-', ...
            'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    end
end

y_all = dat_power(:);
g_all = repmat((1:4)', nSubj, 1);
valid_all = ~isnan(y_all);
boxplot(y_all(valid_all), g_all(valid_all), 'Colors', [0.45 0.45 0.45], ...
    'Symbol', '', 'Widths', 0.5);

for c = 1:4
    y_c = dat_power(c, :);
    valid_c = ~isnan(y_c);
    x_jit = c + (rand(1, sum(valid_c)) - 0.5) * 0.10;
    scatter(x_jit, y_c(valid_c), 170, colors(c,:), 'filled', ...
        'MarkerEdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.7);
end
yline(1, 'k--', 'LineWidth', 1.2);
set(gca, 'XTick', 1:4, 'XTickLabel', strcat(condLabels, ' Contrast'), ...
    'FontSize', 15, 'Box', 'off');
xlim([0.5 4.5]);
ylabel('Gamma Power Ratio (Stimulus/Baseline)');
title('Gamma Power Increase', 'FontSize', 30, 'FontWeight', 'bold');
saveas(fig_power_statsstyle, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_boxplot_GammaPower_statsStyle.png'));

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

saveas(fig_trl1, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_boxplot_alltrials_GammaFreq.png'));

%% ====================================================================
%  GRAND AVERAGE: Mean detrended power spectrum
%  ====================================================================
close all
fig_grand_psd = figure('Position', [0 0 1512/2 982], 'Color', 'w');
hold on;
grand_line_handles = gobjects(1, 4);
grand_panel_maxabs = 0;

for cond = 1:4
    subj_curves = nan(nSubj, nFreqs);
    for s = 1:nSubj
        pr_mat = all_trial_powratio{cond, s};
        if isempty(pr_mat)
            continue;
        end
        nTrl = size(pr_mat, 1);
        pr_dt_mat = nan(size(pr_mat));
        for trl = 1:nTrl
            pr_dt_mat(trl,:) = detrend_power_ratio(pr_mat(trl,:), scan_freqs, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
        end
        subj_mu = nanmean(pr_dt_mat, 1);
        subj_curves(s, :) = normalize_maxabs_curve(subj_mu);
    end

    mu = nanmean(subj_curves, 1);
    n_valid_subj = sum(~isnan(subj_curves(:,1)));
    if n_valid_subj > 0
        se = nanstd(subj_curves, [], 1) ./ sqrt(n_valid_subj);
    else
        se = nan(size(mu));
    end

    if any(isfinite(mu))
        env_vals = [mu - se, mu + se];
        env_vals = env_vals(isfinite(env_vals));
        if ~isempty(env_vals)
            grand_panel_maxabs = max(grand_panel_maxabs, max(abs(env_vals)));
        end

        faceC = 0.8 * colors(cond,:) + 0.2 * [1 1 1];
        % SEB
        % patch([scan_freqs, fliplr(scan_freqs)], ...
        %    [mu - se, fliplr(mu + se)], ...
        %    colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.25);

        % Stronger edge-only smoothing for 30-39 Hz and 81-90 Hz.
        mu_plot = smooth_reflective_edges(mu, scan_freqs, [40 80], 5, 11);
        grand_line_handles(cond) = plot(scan_freqs, mu_plot, '-', ...
            'Color', colors(cond,:), 'LineWidth', 5);

        % Xlines for Median Gamma Freqs
        % md_pf = nanmedian(all_trial_median_single(cond, :));
        % if ~isnan(md_pf)
        %     xline(md_pf, '--', 'Color', colors(cond,:), 'LineWidth', 1.2, 'Alpha', 0.7);
        % end
    end
end

yline(0, 'k--', 'LineWidth', 0.5);
grand_panel_maxabs = max(grand_panel_maxabs, eps);
xlim([30 90]);
ylim([-grand_panel_maxabs grand_panel_maxabs]);
xlabel('Frequency [Hz]');
ylabel('\Delta Power Ratio (normalized)');
title(sprintf('Grand Average Detrended Power Spectrum'), 'FontSize', 25, 'FontWeight', 'bold');
set(gca, 'FontSize', 15);
valid_handles = isgraphics(grand_line_handles);
if any(valid_handles)
    legend(grand_line_handles(valid_handles), condLabels(valid_handles), ...
        'Location', 'best', 'FontSize', 15);
end
saveas(fig_grand_psd, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_grand_average_power_spectrum.png'));

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
                pr_dt_mat(trl,:) = detrend_power_ratio(pr_mat(trl,:), scan_freqs, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
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
saveas(fig_all, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_all_subjects.png'));

%% ====================================================================
%  DETECTION RATE FIGURE (single peak)
%  ====================================================================
fig_det = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Trial-Level Peak Detection Rate (mean across subjects)', ...
    'FontSize', 18, 'FontWeight', 'bold');

subplot(1, 1, 1); hold on;
dr = all_trial_detrate_single;

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
title('Single Peak', 'FontSize', 16, 'FontWeight', 'bold');

saveas(fig_det, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_detection_rate.png'));

%% ====================================================================
%  SENSITIVITY ANALYSIS (thresholds and detrending order)
%  ====================================================================
if sensitivity_enable
    sens_win_names = {'full', 'early', 'late'};
    sens_comp_stats = {all_component_selection_stats_full, all_component_selection_stats_early, all_component_selection_stats_late};
    sens_comp_cubes = {all_trial_powratio_components_full, all_trial_powratio_components_early, all_trial_powratio_components_late};

    for win_idx = 1:3
        wname = sens_win_names{win_idx};
        subj_stats_source = sens_comp_stats{win_idx};
        comp_cube_source = sens_comp_cubes{win_idx};

        fprintf('\n============================================================\n');
        fprintf('Sensitivity analysis (%s window): threshold robustness\n', wname);
        fprintf('============================================================\n');

        base_cfg = struct( ...
        'minEigThr', min_eigval_hard, ...
        'minCorrThr', min_corr_hard, ...
        'minOccFrontRatioThr', min_occfront_ratio, ...
        'minGammaThr', min_gamma_hard, ...
        'outlierRatioThr', outlier_ratio_thr, ...
        'outlierMadMult', outlier_mad_mult, ...
        'polyOrder', poly_order, ...
        'peakMinPromFrac', peak_min_prom_frac, ...
        'peakMinDistanceHz', peak_min_distance_hz);

    cfg_list = struct('sweepParam', {}, 'sweepScale', {}, ...
        'minEigThr', {}, 'minCorrThr', {}, 'minOccFrontRatioThr', {}, 'minGammaThr', {}, ...
        'outlierRatioThr', {}, 'outlierMadMult', {}, 'polyOrder', {}, ...
        'peakMinPromFrac', {}, 'peakMinDistanceHz', {});
    cfg_list(1).sweepParam = 'default';
    cfg_list(1).sweepScale = 1.0;
    cfg_list(1).minEigThr = base_cfg.minEigThr;
    cfg_list(1).minCorrThr = base_cfg.minCorrThr;
    cfg_list(1).minOccFrontRatioThr = base_cfg.minOccFrontRatioThr;
    cfg_list(1).minGammaThr = base_cfg.minGammaThr;
    cfg_list(1).outlierRatioThr = base_cfg.outlierRatioThr;
    cfg_list(1).outlierMadMult = base_cfg.outlierMadMult;
    cfg_list(1).polyOrder = base_cfg.polyOrder;
    cfg_list(1).peakMinPromFrac = base_cfg.peakMinPromFrac;
    cfg_list(1).peakMinDistanceHz = base_cfg.peakMinDistanceHz;

    sweep_params = {'minEigThr', 'minCorrThr', 'minOccFrontRatioThr', 'minGammaThr', ...
        'outlierRatioThr', 'outlierMadMult', 'peakMinPromFrac', 'peakMinDistanceHz'};
    for pi = 1:numel(sweep_params)
        pname = sweep_params{pi};
        for si = 1:numel(sensitivity_scale)
            sc = sensitivity_scale(si);
            if abs(sc - 1) < 1e-12
                continue;
            end
            cfg_i = base_cfg;
            cfg_i.(pname) = base_cfg.(pname) * sc;
            cfg_list(end + 1).sweepParam = pname; %#ok<AGROW>
            cfg_list(end).sweepScale = sc;
            cfg_list(end).minEigThr = cfg_i.minEigThr;
            cfg_list(end).minCorrThr = cfg_i.minCorrThr;
            cfg_list(end).minOccFrontRatioThr = cfg_i.minOccFrontRatioThr;
            cfg_list(end).minGammaThr = cfg_i.minGammaThr;
            cfg_list(end).outlierRatioThr = cfg_i.outlierRatioThr;
            cfg_list(end).outlierMadMult = cfg_i.outlierMadMult;
            cfg_list(end).polyOrder = cfg_i.polyOrder;
            cfg_list(end).peakMinPromFrac = cfg_i.peakMinPromFrac;
            cfg_list(end).peakMinDistanceHz = cfg_i.peakMinDistanceHz;
        end
    end
    for ip = 1:numel(sensitivity_poly_orders)
        po = sensitivity_poly_orders(ip);
        if po == poly_order
            continue;
        end
        cfg_i = base_cfg;
        cfg_i.polyOrder = po;
        cfg_list(end + 1).sweepParam = 'polyOrder'; %#ok<AGROW>
        cfg_list(end).sweepScale = po / max(poly_order, eps);
        cfg_list(end).minEigThr = cfg_i.minEigThr;
        cfg_list(end).minCorrThr = cfg_i.minCorrThr;
        cfg_list(end).minOccFrontRatioThr = cfg_i.minOccFrontRatioThr;
        cfg_list(end).minGammaThr = cfg_i.minGammaThr;
        cfg_list(end).outlierRatioThr = cfg_i.outlierRatioThr;
        cfg_list(end).outlierMadMult = cfg_i.outlierMadMult;
        cfg_list(end).polyOrder = cfg_i.polyOrder;
        cfg_list(end).peakMinPromFrac = cfg_i.peakMinPromFrac;
        cfg_list(end).peakMinDistanceHz = cfg_i.peakMinDistanceHz;
    end

        fprintf('Evaluating %d sensitivity configurations (%s window).\n', numel(cfg_list), wname);
        sens_rows = [];
        sens_sweep_param = cell(0, 1);
        sens_sweep_scale = [];
        for cfg_id = 1:numel(cfg_list)
            cfg_cur = cfg_list(cfg_id);

            slope_cfg = nan(1, nSubj);
            delta_cfg = nan(1, nSubj);

            for subj = 1:nSubj
                subj_stats = subj_stats_source{subj};
            if isempty(subj_stats) || ~isfield(subj_stats, 'candidate_table')
                continue;
            end
            cand = subj_stats.candidate_table;
            eig_vec = cand.eigenvalue(:);
            corr_vec = cand.corr(:);
            ratio_vec = cand.ratio(:);
            gamma_vec = cand.gamma(:);
            nComp = numel(eig_vec);
            finite_metrics = isfinite(eig_vec) & isfinite(corr_vec) & isfinite(ratio_vec) & isfinite(gamma_vec);

            hard_eligible_raw_cfg = finite_metrics & ...
                (eig_vec >= cfg_cur.minEigThr) & ...
                (corr_vec >= cfg_cur.minCorrThr) & ...
                (ratio_vec >= cfg_cur.minOccFrontRatioThr) & ...
                (gamma_vec >= cfg_cur.minGammaThr);

            hard_eligible_cfg = hard_eligible_raw_cfg;
            if exclude_dominant_outlier
                [hard_eligible_cfg, ~] = exclude_dominant_top_outlier( ...
                    eig_vec, hard_eligible_raw_cfg, cfg_cur.outlierRatioThr, cfg_cur.outlierMadMult, outlier_min_rest);
            end

            sel_idx = find(hard_eligible_cfg & isfinite(eig_vec));
            if isempty(sel_idx)
                finite_idx = find(finite_metrics & isfinite(eig_vec));
                if isempty(finite_idx)
                    continue;
                end
                [~, fin_ord] = sort(eig_vec(finite_idx), 'descend');
                sel_idx = finite_idx(fin_ord(1));
            end
            [~, sel_ord] = sort(eig_vec(sel_idx), 'descend');
            sel_idx = sel_idx(sel_ord);

            sel_w = eig_vec(sel_idx)';
            sel_w(~isfinite(sel_w) | sel_w <= 0) = 0;
            if sum(sel_w) <= 0
                sel_w = ones(1, numel(sel_idx));
            end
            sel_w = sel_w / sum(sel_w);

            cond_medians = nan(4, 1);
            for cond = 1:4
                comp_cube = comp_cube_source{cond, subj};
                if isempty(comp_cube)
                    continue;
                end
                nTrl = size(comp_cube, 2);
                trl_peaks = nan(nTrl, 1);
                for trl = 1:nTrl
                    pr_comp = squeeze(comp_cube(1:nComp, trl, :));
                    if isempty(pr_comp)
                        continue;
                    end
                    if isvector(pr_comp)
                        pr_comp = pr_comp(:)';
                    end
                    if size(pr_comp, 2) ~= nFreqs
                        pr_comp = pr_comp';
                    end
                    pr_sel = pr_comp(sel_idx, :);
                    if isempty(pr_sel)
                        continue;
                    end

                    num = nansum(pr_sel .* sel_w(:), 1);
                    den = nansum((~isnan(pr_sel)) .* sel_w(:), 1);
                    pr = num ./ den;
                    if all(~isfinite(pr))
                        continue;
                    end
                    pr_dt = detrend_power_ratio(pr, scan_freqs, cfg_cur.polyOrder, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
                    pr_dt_smooth = movmean(pr_dt, 5);

                    mprom = max(0, max(pr_dt_smooth) * cfg_cur.peakMinPromFrac);
                    [pks, locs] = findpeaks(pr_dt_smooth, scan_freqs, ...
                        'MinPeakProminence', mprom, ...
                        'MinPeakDistance', cfg_cur.peakMinDistanceHz);
                    if ~isempty(pks)
                        [~, best_pk] = max(pks);
                        trl_peaks(trl) = locs(best_pk);
                    end
                end
                cond_medians(cond) = median(trl_peaks(~isnan(trl_peaks)));
            end

            vx = ~isnan(cond_medians);
            if sum(vx) >= 2
                p_fit = polyfit(find(vx), cond_medians(vx)', 1);
                slope_cfg(subj) = p_fit(1);
            end
            if ~isnan(cond_medians(1)) && ~isnan(cond_medians(4))
                delta_cfg(subj) = cond_medians(4) - cond_medians(1);
            end
        end

        slope_stats_cfg = compute_one_sample_stats(slope_cfg);
        delta_stats_cfg = compute_one_sample_stats(delta_cfg);
        sens_sweep_param{end + 1, 1} = cfg_cur.sweepParam; %#ok<AGROW>
        sens_sweep_scale(end + 1, 1) = cfg_cur.sweepScale; %#ok<AGROW>
        sens_rows = [sens_rows; ...
            cfg_id, ...
            cfg_cur.minEigThr, cfg_cur.minCorrThr, cfg_cur.minOccFrontRatioThr, cfg_cur.minGammaThr, ...
            cfg_cur.outlierRatioThr, cfg_cur.outlierMadMult, cfg_cur.polyOrder, ...
            cfg_cur.peakMinPromFrac, cfg_cur.peakMinDistanceHz, ...
            slope_stats_cfg.n, slope_stats_cfg.mean, slope_stats_cfg.p, slope_stats_cfg.ci_low, slope_stats_cfg.ci_high, slope_stats_cfg.cohens_d, ...
            delta_stats_cfg.n, delta_stats_cfg.mean, delta_stats_cfg.p, delta_stats_cfg.ci_low, delta_stats_cfg.ci_high, delta_stats_cfg.cohens_d]; %#ok<AGROW>
        end

        if isempty(sens_rows)
            fprintf('Sensitivity analysis (%s window) skipped: no valid outputs.\n', wname);
            sens_results_win = table();
        else
        sens_results_win = array2table(sens_rows, 'VariableNames', ...
            {'cfgID', 'minEigThr', 'minCorrThr', 'minOccFrontRatioThr', 'minGammaThr', ...
            'outlierRatioThr', 'outlierMadMult', 'polyOrder', 'peakMinPromFrac', 'peakMinDistanceHz', ...
            'nSlope', 'meanSlope', 'pSlope', 'ciSlopeLow', 'ciSlopeHigh', 'dSlope', ...
            'nDelta', 'meanDelta', 'pDelta', 'ciDeltaLow', 'ciDeltaHigh', 'dDelta'});
            sens_results_win.sweepParam = string(sens_sweep_param);
            sens_results_win.sweepScale = sens_sweep_scale;
            sens_results_win.sigSlope = sens_results_win.pSlope < 0.05;
            sens_results_win.sigDelta = sens_results_win.pDelta < 0.05;

            fprintf('Sensitivity (%s): %d/%d significant slope, %d/%d delta (p<0.05)\n', ...
                wname, sum(sens_results_win.sigSlope), height(sens_results_win), ...
                sum(sens_results_win.sigDelta), height(sens_results_win));

            [~, sort_delta_idx] = sort(sens_results_win.pDelta, 'ascend');
            n_show = min(10, height(sens_results_win));
            disp(sens_results_win(sort_delta_idx(1:n_show), :));

            fig_sens = figure('Position', [0 0 1512 982], 'Color', 'w');
            tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
            nexttile; hold on;
            scatter(sens_results_win.cfgID, sens_results_win.meanSlope, 120, ...
                sens_results_win.pSlope, 'filled');
            yline(0, 'k--', 'LineWidth', 1.0);
            cb1 = colorbar; cb1.Label.String = 'p-value (slope)';
            xlabel('Configuration ID');
            ylabel('Mean slope [Hz/condition]');
            title(sprintf('Sensitivity (%s): slope', wname));
            box on;

            nexttile; hold on;
            scatter(sens_results_win.cfgID, sens_results_win.meanDelta, 120, ...
                sens_results_win.pDelta, 'filled');
            yline(0, 'k--', 'LineWidth', 1.0);
            cb2 = colorbar; cb2.Label.String = 'p-value (delta)';
            xlabel('Configuration ID');
            ylabel('Mean \Delta(100%-25%) [Hz]');
            title(sprintf('Sensitivity (%s): delta', wname));
            box on;
            saveas(fig_sens, fullfile(fig_save_dir_ged, sprintf('GCP_eeg_GED_sensitivity_summary_%s.png', wname)));
        end

        if win_idx == 1
            sensitivity_results = sens_results_win;
            sensitivity_results_full = sens_results_win;
        elseif win_idx == 2
            sensitivity_results_early = sens_results_win;
        else
            sensitivity_results_late = sens_results_win;
        end
    end
end

%% Save results
if ispc
    save_path = 'W:\Students\Arne\GCP\data\features\GCP_eeg_GED.mat';
else
    save_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_eeg_GED.mat';
end
save(save_path, ...
    'all_trial_powratio', ...
    'all_trial_powratio_early', 'all_trial_powratio_late', ...
    'all_trial_peaks_single', 'all_trial_centroid', ...
    'all_trial_peaks_single_early', 'all_trial_peaks_single_late', ...
    'all_trial_mean_single', 'all_trial_median_single', ...
    'all_trial_mean_single_early', 'all_trial_median_single_early', ...
    'all_trial_mean_single_late', 'all_trial_median_single_late', ...
    'all_trial_mean_centroid', 'all_trial_median_centroid', ...
    'all_trial_detrate_single', 'all_trial_detrate_single_early', 'all_trial_detrate_single_late', ...
    'all_trial_detrate_centroid', 'all_trial_gamma_power', 'all_trial_gamma_power_early', 'all_trial_gamma_power_late', ...
    'all_topos', 'all_topos_early', 'all_topos_late', 'all_topo_labels', 'all_eigenvalues', ...
    'all_selected_comp_idx', 'all_selected_comp_corr', 'all_selected_comp_eval', ...
    'all_selected_comp_indices_multi', 'all_selected_comp_weights', ...
    'all_component_selection_stats', 'all_component_selection_stats_full', 'all_component_selection_stats_early', 'all_component_selection_stats_late', 'warning_log', ...
    'all_trial_powratio_bench', 'all_trial_powratio_dt_bench', 'all_trial_powratio_components', ...
    'all_trial_powratio_components_full', 'all_trial_powratio_components_early', 'all_trial_powratio_components_late', ...
    'all_trial_powratio_bench_early', 'all_trial_powratio_bench_late', ...
    'all_trial_powratio_dt_bench_early', 'all_trial_powratio_dt_bench_late', ...
    'benchmark_methods', 'raw_reference_definition', ...
    'benchmark_metric_detectability', 'benchmark_metric_prominence', ...
    'benchmark_metric_separation_slope', 'benchmark_metric_separation_delta', ...
    'benchmark_metric_reliability_trialcv', 'benchmark_metric_reliability_subjspread', ...
    'primary_slope_stats', 'primary_delta_stats', 'sensitivity_results', ...
    'sensitivity_results_full', 'sensitivity_results_early', 'sensitivity_results_late', ...
    'all_top5_corrs', 'all_top5_evals', 'all_top5_topos', 'all_simulated_templates', ...
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
end % run_mode == 'full'

function method_ratios = compute_method_ratios_from_components(ratio_all, x_stim, x_base, raw_w, W_top, W_combined, selected_idx, w_combined, nBenchmarkMethods)
% ratio_all: power ratios indexed by component (1 = highest eigenvalue); must match W_top (col 1) and selected_idx
method_ratios = nan(nBenchmarkMethods, 1);

pow_stim_chan = mean(x_stim.^2, 2);
pow_base_chan = mean(x_base.^2, 2);
raw_pow_stim = sum(raw_w(:) .* pow_stim_chan(:));
raw_pow_base = sum(raw_w(:) .* pow_base_chan(:));
if isfinite(raw_pow_stim) && isfinite(raw_pow_base) && raw_pow_base > 0
    method_ratios(1) = raw_pow_stim / raw_pow_base;
end

if ~isempty(W_top)
    ratio_top = ratio_all(1);
    if isfinite(ratio_top)
        method_ratios(2) = ratio_top;
    end
end

if ~isempty(W_combined)
    ratio_vec = ratio_all(selected_idx);
    valid_comp = isfinite(ratio_vec);
    if any(valid_comp)
        w_use = w_combined(valid_comp);
        if sum(w_use) <= 0
            w_use = ones(1, sum(valid_comp)) / sum(valid_comp);
        else
            w_use = w_use / sum(w_use);
        end
        method_ratios(3) = sum(ratio_vec(valid_comp)' .* w_use);
    end
end
end

function trl_peaks_single = detect_single_peaks_from_powratio(powratio_trials, scan_freqs, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges, peak_min_prom_frac, peak_min_distance_hz)
nTrl = size(powratio_trials, 1);
trl_peaks_single = nan(nTrl, 1);
for trl = 1:nTrl
    pr = powratio_trials(trl, :);
    if all(isnan(pr))
        continue;
    end
    pr_dt = detrend_power_ratio(pr, scan_freqs, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
    pr_dt_smooth = movmean(pr_dt, 5);
    mprom = max(0, max(pr_dt_smooth) * peak_min_prom_frac);
    [pks, locs] = findpeaks(pr_dt_smooth, scan_freqs, ...
        'MinPeakProminence', mprom, ...
        'MinPeakDistance', peak_min_distance_hz);
    if ~isempty(pks)
        [~, best_pk] = max(pks);
        trl_peaks_single(trl) = locs(best_pk);
    end
end
end

function [slope_vec, delta_vec] = compute_condition_separation_from_matrix(dat)
nSubj_local = size(dat, 2);
slope_vec = nan(1, nSubj_local);
delta_vec = nan(1, nSubj_local);
for s = 1:nSubj_local
    y = dat(:, s);
    valid_y = isfinite(y);
    if sum(valid_y) >= 2
        p = polyfit(find(valid_y), y(valid_y)', 1);
        slope_vec(s) = p(1);
    end
    if isfinite(y(1)) && isfinite(y(4))
        delta_vec(s) = y(4) - y(1);
    end
end
end

function plot_gamma_window_panel(dat, condLabels, colors, nSubj)
mu = nanmean(dat, 2);
sem = nanstd(dat, [], 2) ./ sqrt(sum(~isnan(dat), 2));
for s = 1:nSubj
    y_subj = dat(:, s);
    valid_subj = ~isnan(y_subj);
    if sum(valid_subj) >= 2
        plot(find(valid_subj), y_subj(valid_subj), '-', 'Color', [0.80 0.80 0.80], 'LineWidth', 1);
    end
end
for c = 1:4
    y_c = dat(c, :);
    valid_c = ~isnan(y_c);
    x_jit = c + (rand(1, sum(valid_c)) - 0.5) * 0.12;
    scatter(x_jit, y_c(valid_c), 60, colors(c,:), 'filled', ...
        'MarkerEdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.6);
end
errorbar(1:4, mu, sem, 'k-o', 'LineWidth', 1.8, 'CapSize', 8, ...
    'MarkerFaceColor', [0.1 0.1 0.1], 'MarkerSize', 6);
set(gca, 'XTick', 1:4, 'XTickLabel', strcat(condLabels, ' Contrast'), 'FontSize', 11, 'Box', 'off');
xlim([0.5 4.5]);
end

function [eligible_mask_out, dominant_outlier_idx] = exclude_dominant_top_outlier(eigvals_sorted, eligible_mask_in, ratio_thr, mad_mult, min_rest_n)
eligible_mask_out = logical(eligible_mask_in(:));
dominant_outlier_idx = [];
if isempty(eigvals_sorted) || isempty(eligible_mask_out)
    return;
end

eigvals_sorted = eigvals_sorted(:);
if numel(eigvals_sorted) ~= numel(eligible_mask_out)
    return;
end
if nargin < 5 || isempty(min_rest_n)
    min_rest_n = 3;
end

eligible_idx = find(eligible_mask_out & isfinite(eigvals_sorted) & eigvals_sorted > 0);
if numel(eligible_idx) < 2
    return;
end

[~, elig_ord] = sort(eigvals_sorted(eligible_idx), 'descend');
eligible_ranked_idx = eligible_idx(elig_ord);
top_idx = eligible_ranked_idx(1);
lambda_top = eigvals_sorted(top_idx);
lambda_second = eigvals_sorted(eligible_ranked_idx(2));
rest_vals = eigvals_sorted(eligible_ranked_idx(2:end));
if numel(rest_vals) < min_rest_n || ~isfinite(lambda_top) || ~isfinite(lambda_second) || lambda_second <= 0
    return;
end

rest_log = log(rest_vals);
rest_log = rest_log(isfinite(rest_log));
if numel(rest_log) < min_rest_n
    return;
end

lambda_top_log = log(max(lambda_top, eps));
ratio12 = lambda_top / max(lambda_second, eps);
med_rest = median(rest_log);
mad_rest = mad(rest_log, 1);
if ~isfinite(mad_rest) || mad_rest <= eps
    mad_criterion = lambda_top_log > med_rest;
else
    mad_criterion = lambda_top_log > (med_rest + mad_mult * mad_rest);
end

if ratio12 >= ratio_thr && mad_criterion
    eligible_mask_out(top_idx) = false;
    dominant_outlier_idx = top_idx;
end
end

function warning_log = append_subject_warning(warning_log, subject_id, code, message, metrics)
if nargin < 5 || isempty(metrics)
    metrics = struct();
end
entry = struct('subject', subject_id, 'code', code, 'message', message, 'metrics', metrics);
if isempty(warning_log)
    warning_log = entry;
else
    warning_log(end+1) = entry;
end
end

function print_subject_warning_summary(warning_log)
fprintf('\n============================================================\n');
fprintf('Per-subject warning summary\n');
fprintf('============================================================\n');
if isempty(warning_log)
    fprintf('No warnings were logged across subjects.\n');
    return;
end

for wi = 1:numel(warning_log)
    w = warning_log(wi);
    fprintf('%3d) subject=%s | code=%s\n', wi, w.subject, w.code);
    fprintf('     %s\n', w.message);
    if strcmp(w.code, 'NO_HARD_ELIGIBLE_COMPONENTS')
        m = w.metrics;
        fprintf(['     hard-gate diagnostics: nSearch=%d, finite=%d, passEig=%d, passCorr=%d, ', ...
                 'passRatio=%d, passGamma=%d, passAllRaw=%d, excludedOutlier=%d\n'], ...
            m.n_search, m.n_finite_metrics, m.n_pass_eig, m.n_pass_corr, ...
            m.n_pass_ratio, m.n_pass_gamma, m.n_pass_all_raw, m.n_excluded_dominant_outlier);
        if isfield(m, 'top_fail_idx') && isfinite(m.top_fail_idx)
            fprintf(['     top failing component: C%d | eig=%.3f (thr=%.3f, pass=%d), ', ...
                     'corr=%.3f (thr=%.3f, pass=%d), ratio=%.3f (thr=%.3f, pass=%d), ', ...
                     'gamma=%.3f (thr=%.3f, pass=%d)\n'], ...
                m.top_fail_idx, m.top_fail_eig, m.thr_eig, m.top_fail_pass_eig, ...
                m.top_fail_corr, m.thr_corr, m.top_fail_pass_corr, ...
                m.top_fail_ratio, m.thr_ratio, m.top_fail_pass_ratio, ...
                m.top_fail_gamma, m.thr_gamma, m.top_fail_pass_gamma);
        end
    elseif strcmp(w.code, 'NO_HARD_COMPONENTS') || strcmp(w.code, 'TOO_FEW_HARD_COMPONENTS')
        m = w.metrics;
        if isfield(m, 'n_valid_components')
            fprintf('     combined-selection diagnostics: valid=%d, minRequired=%d\n', ...
                m.n_valid_components, m.min_required);
        else
            fprintf('     combined-selection diagnostics: nHardEligible=%d, nFiniteScores=%d\n', ...
                m.n_hard_eligible, m.n_finite_scores);
        end
    elseif strcmp(w.code, 'DOMINANT_OUTLIER_EXCLUDED')
        m = w.metrics;
        fprintf('     outlier diagnostics: C%d, lambda1=%.3f, lambda2=%.3f, ratio12=%.3f\n', ...
            m.component_idx, m.lambda1, m.lambda2, m.lambda1_lambda2_ratio);
    end
end

codes = {warning_log.code};
[u_codes, ~, ic] = unique(codes);
fprintf('\nWarning counts by code:\n');
for ci = 1:numel(u_codes)
    fprintf('  %s: %d\n', u_codes{ci}, sum(ic == ci));
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

function inspection = compute_threshold_inspection(eig_vec, corr_vec, ratio_vec, gamma_vec, ...
    thr_eig, thr_corr, thr_ratio, thr_gamma, target_counts, n_search)
inspection = struct();
inspection.n_search = n_search;
inspection.target_counts = target_counts(:)';
inspection.thr_default = [thr_eig, thr_corr, thr_ratio, thr_gamma];
inspection.base_pass_count = sum(isfinite(eig_vec) & isfinite(corr_vec) & isfinite(ratio_vec) & isfinite(gamma_vec) & ...
    eig_vec >= thr_eig & corr_vec >= thr_corr & ratio_vec >= thr_ratio & gamma_vec >= thr_gamma);

finite_mask = isfinite(eig_vec) & isfinite(corr_vec) & isfinite(ratio_vec) & isfinite(gamma_vec);
inspection.n_finite_metrics = sum(finite_mask);

scale_cap_eig = eig_vec ./ max(thr_eig, eps);
scale_cap_corr = corr_vec ./ max(thr_corr, eps);
scale_cap_ratio = ratio_vec ./ max(thr_ratio, eps);
scale_cap_gamma = gamma_vec ./ max(thr_gamma, eps);

component_scale_cap = min([scale_cap_eig(:), scale_cap_corr(:), scale_cap_ratio(:), scale_cap_gamma(:)], [], 2);
component_scale_cap(~finite_mask(:)) = -Inf;
component_scale_cap(component_scale_cap < 0) = -Inf;

valid_caps = component_scale_cap(isfinite(component_scale_cap));
valid_caps = valid_caps(valid_caps >= 0);
if isempty(valid_caps)
    max_scale_any = NaN;
else
    max_scale_any = max(valid_caps);
end
inspection.max_feasible_scale = max_scale_any;
inspection.max_pass_at_zero = sum(finite_mask & eig_vec >= 0 & corr_vec >= 0 & ratio_vec >= 0 & gamma_vec >= 0);

n_targets = numel(target_counts);
inspection.results = repmat(struct( ...
    'target_count', NaN, ...
    'joint_achievable', false, ...
    'joint_required_scale', NaN, ...
    'required_thr_eig', NaN, ...
    'required_thr_corr', NaN, ...
    'required_thr_ratio', NaN, ...
    'required_thr_gamma', NaN, ...
    'pass_count_at_required_scale', NaN, ...
    'by_criterion', struct()), 1, n_targets);

caps_sorted = sort(valid_caps, 'descend');
for ti = 1:n_targets
    k_target = target_counts(ti);
    inspection.results(ti).target_count = k_target;
    if k_target <= 0
        inspection.results(ti).joint_achievable = true;
        inspection.results(ti).joint_required_scale = 1;
        inspection.results(ti).required_thr_eig = thr_eig;
        inspection.results(ti).required_thr_corr = thr_corr;
        inspection.results(ti).required_thr_ratio = thr_ratio;
        inspection.results(ti).required_thr_gamma = thr_gamma;
        inspection.results(ti).pass_count_at_required_scale = inspection.base_pass_count;
    else
        if numel(caps_sorted) >= k_target
            req_scale = caps_sorted(k_target);
            req_scale = min(req_scale, 1);
            req_scale = max(req_scale, 0);
            pass_at_scale = sum(finite_mask & ...
                eig_vec >= (thr_eig * req_scale) & ...
                corr_vec >= (thr_corr * req_scale) & ...
                ratio_vec >= (thr_ratio * req_scale) & ...
                gamma_vec >= (thr_gamma * req_scale));

            inspection.results(ti).joint_achievable = pass_at_scale >= k_target;
            inspection.results(ti).joint_required_scale = req_scale;
            inspection.results(ti).required_thr_eig = thr_eig * req_scale;
            inspection.results(ti).required_thr_corr = thr_corr * req_scale;
            inspection.results(ti).required_thr_ratio = thr_ratio * req_scale;
            inspection.results(ti).required_thr_gamma = thr_gamma * req_scale;
            inspection.results(ti).pass_count_at_required_scale = pass_at_scale;
        end
    end

    by_criterion = struct();
    by_criterion.eig = compute_single_criterion_min_threshold( ...
        eig_vec, finite_mask & corr_vec >= thr_corr & ratio_vec >= thr_ratio & gamma_vec >= thr_gamma, ...
        thr_eig, k_target);
    by_criterion.corr = compute_single_criterion_min_threshold( ...
        corr_vec, finite_mask & eig_vec >= thr_eig & ratio_vec >= thr_ratio & gamma_vec >= thr_gamma, ...
        thr_corr, k_target);
    by_criterion.ratio = compute_single_criterion_min_threshold( ...
        ratio_vec, finite_mask & eig_vec >= thr_eig & corr_vec >= thr_corr & gamma_vec >= thr_gamma, ...
        thr_ratio, k_target);
    by_criterion.gamma = compute_single_criterion_min_threshold( ...
        gamma_vec, finite_mask & eig_vec >= thr_eig & corr_vec >= thr_corr & ratio_vec >= thr_ratio, ...
        thr_gamma, k_target);
    inspection.results(ti).by_criterion = by_criterion;
end
end

function print_threshold_inspection_summary(subjects, all_threshold_inspection, target_counts)
fprintf('\n============================================================\n');
fprintf('Threshold-inspection summary (hard-gate only)\n');
fprintf('============================================================\n');
fprintf(['For each subject, "scale" is the maximum common multiplier applied to all hard thresholds ', ...
         '(eig/corr/ratio/gamma) that still yields at least K valid components.\n']);
fprintf(['Criterion-specific minimums are also reported: each threshold is relaxed alone while ', ...
         'all other thresholds remain fixed.\n']);
fprintf('Current thresholds correspond to scale=1.000; lower scale means relaxed thresholds.\n');

if isempty(all_threshold_inspection)
    fprintf('No threshold-inspection results available.\n');
    return;
end

for subj = 1:numel(subjects)
    ins = all_threshold_inspection{subj};
    if isempty(ins)
        fprintf('  %s: no inspection data\n', subjects{subj});
        continue;
    end
    fprintf('  %s: basePass=%d/%d finite=%d\n', ...
        subjects{subj}, ins.base_pass_count, ins.n_search, ins.n_finite_metrics);
    for ti = 1:numel(target_counts)
        r = ins.results(ti);
        if r.joint_achievable
            fprintf(['    K=%d -> scale=%.3f | eig>=%.3f, corr>=%.3f, ratio>=%.3f, gamma>=%.3f ', ...
                     '(pass=%d)\n'], ...
                r.target_count, r.joint_required_scale, r.required_thr_eig, r.required_thr_corr, ...
                r.required_thr_ratio, r.required_thr_gamma, r.pass_count_at_required_scale);
        else
            fprintf('    K=%d -> not achievable within finite candidates (max nonnegative-pass=%d)\n', ...
                r.target_count, ins.max_pass_at_zero);
        end
        c = r.by_criterion;
        fprintf('      min-by-criterion: eig>=%.3f (%s), corr>=%.3f (%s), ratio>=%.3f (%s), gamma>=%.3f (%s)\n', ...
            c.eig.required_threshold, logical_str(c.eig.achievable), ...
            c.corr.required_threshold, logical_str(c.corr.achievable), ...
            c.ratio.required_threshold, logical_str(c.ratio.achievable), ...
            c.gamma.required_threshold, logical_str(c.gamma.achievable));
    end
end
end

function out = compute_single_criterion_min_threshold(metric_vec, base_mask, default_thr, k_target)
out = struct('achievable', false, 'required_threshold', NaN, 'pass_count', 0, 'pool_count', 0);
if k_target <= 0
    out.achievable = true;
    out.required_threshold = default_thr;
    out.pass_count = sum(base_mask & metric_vec >= default_thr);
    out.pool_count = sum(base_mask);
    return;
end

pool_vals = metric_vec(base_mask);
pool_vals = pool_vals(isfinite(pool_vals));
out.pool_count = numel(pool_vals);
if numel(pool_vals) < k_target
    return;
end

pool_sorted = sort(pool_vals, 'descend');
kth_val = pool_sorted(k_target);
out.required_threshold = min(default_thr, kth_val);
out.pass_count = sum(base_mask & metric_vec >= out.required_threshold);
out.achievable = out.pass_count >= k_target;
end

function txt = logical_str(tf)
if tf
    txt = 'ok';
else
    txt = 'NA';
end
end

function stats = compute_one_sample_stats(x)
x = x(:);
x = x(isfinite(x));
stats = struct('n', numel(x), 'mean', NaN, 'p', NaN, 'tstat', NaN, ...
    'df', NaN, 'ci_low', NaN, 'ci_high', NaN, 'cohens_d', NaN);
if isempty(x)
    return;
end
stats.mean = mean(x);
if numel(x) < 2
    return;
end
[~, p, ci, st] = ttest(x, 0, 'Alpha', 0.05);
stats.p = p;
stats.tstat = st.tstat;
stats.df = st.df;
stats.ci_low = ci(1);
stats.ci_high = ci(2);
sx = std(x);
if isfinite(sx) && sx > 0
    stats.cohens_d = stats.mean / sx;
end
end

function y = smooth_reflective(x, win)
y = x;
if nargin < 2 || isempty(win) || win <= 1 || isempty(x)
    return;
end
win = max(1, round(win));
if mod(win,2) == 0
    win = win + 1;
end
pad = floor(win/2);
if numel(x) <= pad
    y = movmean(x, min(win, numel(x)));
    return;
end
x_pad = [fliplr(x(1:pad)), x, fliplr(x(end-pad+1:end))];
y_pad = movmean(x_pad, win);
y = y_pad(pad+1:end-pad);
end

function y = smooth_reflective_edges(x, freqs, core_band, core_win, edge_win)
y = x;
if isempty(x) || isempty(freqs) || numel(x) ~= numel(freqs)
    return;
end
if nargin < 4 || isempty(core_win)
    core_win = 5;
end
if nargin < 5 || isempty(edge_win)
    edge_win = 11;
end
if nargin < 3 || isempty(core_band) || numel(core_band) ~= 2
    core_band = [40 80];
end

y_core = smooth_reflective(x, core_win);
y_edge = smooth_reflective(x, edge_win);
edge_mask = freqs < core_band(1) | freqs > core_band(2);
y(edge_mask) = y_edge(edge_mask);
y(~edge_mask) = y_core(~edge_mask);
end

function y_dt = detrend_power_ratio(y, x, ord, edge_exclude_n, in_log, flat_edges)
% Detrend power-ratio spectrum for peak detection.
% y: power ratio (positive values)
% x: frequency vector
% ord: polynomial order
% edge_exclude_n: bins to exclude from fit at each edge (default 5)
% in_log: if true, fit polynomial to log(y), return residual in log space (default true)
% flat_edges: if true, constrain residual to zero at first/last frequency (default true)
y_dt = y;
if isempty(y) || isempty(x) || numel(y) ~= numel(x)
    return;
end
if nargin < 4 || isempty(edge_exclude_n)
    edge_exclude_n = 5;
end
if nargin < 5
    in_log = true;
end
if nargin < 6
    flat_edges = true;
end
y_fit = y;
if in_log
    y_fit = log(max(y, eps));  % preserve shape; y(:) would force column and cause implicit expansion downstream
end
y_dt = detrend_poly_stable(y_fit, x, ord, edge_exclude_n, flat_edges);
end

function y_dt = detrend_poly_stable(y, x, ord, edge_exclude_n, flat_edges)
y_dt = y;
if isempty(y) || isempty(x) || numel(y) ~= numel(x)
    return;
end
if nargin < 4 || isempty(edge_exclude_n)
    edge_exclude_n = 0;
end
if nargin < 5
    flat_edges = true;
end
edge_exclude_n = max(0, round(edge_exclude_n));
n = numel(x);

fit_idx = (1 + edge_exclude_n):(n - edge_exclude_n);
if numel(fit_idx) < (ord + 1)
    fit_idx = 1:n;
end

% Restrict to valid indices for both y and x, and force same orientation
% to avoid implicit expansion when y is column and x is row (e.g. y 61x1, x 1x61)
fit_idx = fit_idx(fit_idx <= min(numel(y), numel(x)) & fit_idx >= 1);
yv = y(fit_idx);
xv = x(fit_idx);
valid_mask = isfinite(yv(:)) & isfinite(xv(:));
valid_fit = fit_idx(valid_mask);
if numel(valid_fit) < (ord + 1)
    valid_fit = find(isfinite(y(:)) & isfinite(x(:)));
end
if numel(valid_fit) < (ord + 1)
    return;
end

p = polyfit(x(valid_fit), y(valid_fit), ord);
baseline = polyval(p, x);
baseline = reshape(baseline, size(y));  % match y orientation to avoid implicit expansion
y_dt = y - baseline;

% Constrain residual to zero at edges
if flat_edges && n >= 2
    x1 = x(1);
    xn = x(n);
    if isfinite(x1) && isfinite(xn) && xn > x1
        d1 = y_dt(1);
        dn = y_dt(n);
        if isfinite(d1) && isfinite(dn)
            % Linear ramp from d1 at x1 to dn at xn; subtract so residual -> 0 at edges
            ramp = d1 + (dn - d1) * (x(:) - x1) / (xn - x1);
            ramp = reshape(ramp, size(y_dt));  % match orientation
            y_dt = y_dt - ramp;
        end
    end
end
end
