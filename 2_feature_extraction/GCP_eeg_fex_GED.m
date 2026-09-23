%% GCP Gamma Peak Frequency and Power with Generalized Eigendecomposition (GED)
%
% Broadband GED and component selection (per subject)
%   - Pool trials across contrast and compute gamma-band covariances
%      (30-90 Hz FIR) for the stimulus interval (0-2000 ms) and baseline.
%   - Solve GED (S_stim * w = lambda * S_base * w) with
%      regularisation and rank candidate components by eigenvalue.
%   - Retain candidates that pass SNR (lambda), spectral peak-form (PF),
%      posterior>frontal and posterior>temporal dominance (post_front,
%      post_temp), low whole-scalp single-channel peak fraction (rejects
%      ultra-focal EMG anywhere on the montage), and rising HF-slope EMG
%      exclusion. Combine all eligible components among the first 10
%      with eigenvalue-proportional weights.
%   - A signed occipital/frontal spatial template is used only to align
%      component polarity before scoring.
%
% Trial-level spectral scanning (per subject, condition, trial)
%   - Project each trial through the subject combined GED filter and
%      compute spectrum-based power scans on a 30-90 Hz grid in dB.
%   - Flag numerical-instability cases (near-floor baseline power across
%      many frequencies/components) and exclude unstable trials automatically.
%   - Detect per-trial peak frequency and define peak power as the mean power
%      within peak frequency +/- 5 Hz.
%
% Outputs
%   - Trial-level peak frequency/power (trials_peaks, trials_powratio,
%     trials_outlier_mask_power) for trial hypotheses / rainclouds.
%   - Condition-averaged spectral peaks (all_condition_peak_freq/power),
%     from mean trial spectra (no IQR trial exclusion), used by subject-level
%     master matrix, boxplots, and rainclouds.
%   - Topographies / spectra / Haufe / reconstructed freq for viz scripts.
%   - Optional GED-projected condition TFRs (toggle do_tfr below).
%   - Optional diagnostic/summary figures (toggle do_plots below).
%
% Helpers live in paths.code/2_feature_extraction/GED-helpers (added to path at startup).

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('GCP');
if ispc
    addpath('C:\Users\Administrator\Documents\GitHub\GCP\2_feature_extraction\GED-helpers')
else
    addpath('/Users/Arne/Documents/GitHub/GCP/2_feature_extraction/GED-helpers');
end
nSubj = length(subjects);
total_runtime_tic = tic;

%% Parameters

% Toggles
do_tfr = true;   % set false to skip GED-projected TFR feature extraction
do_plots = true; % set false to skip diagnostic/summary figure generation
topo_gridscale = 300; % topo interpolation density (was 300; raise for publication)

% Time intervals
baseline_window = [-1.5, -0.5];
stim_window = [0, 2.0];

% Gamma analysis (frequency grid and FieldTrip mtmfft multitaper bandwidth)
analysis_freq_range = [30 90];
scan_freq_step_hz = 1; % frequency grid step (Hz); powratio_*_freq_smooth_bins count bins on this grid
scan_freqs = analysis_freq_range(1):scan_freq_step_hz:analysis_freq_range(2);
nFreqs = length(scan_freqs);
mtmfft_tapsmofrq_hz = 3; % FieldTrip cfg.tapsmofrq for mtmfft (Hz)

% GED
lambda = 0.05;              % regularization
ged_search_n = 10;          % search first N GED components
min_eigval = 1.1;            % minimum GED eigenvalue (lambda >= 1.1)
min_powspctrm_form = 0.80;  % minimum PF (powspctrm-form) score for candidate eligibility
random_seed = 123;
% TFR (GED-projected multitaper)
tfr_foi = 30:1:90;
tfr_toi = -1.75:0.05:2.00;
tfr_win_sec = 0.50;
tfr_tapsmofrq = 5;
tfr_baseline_window = [baseline_window(1) + tfr_win_sec / 2, ...
    baseline_window(2) - tfr_win_sec / 2];
if tfr_baseline_window(1) > tfr_baseline_window(2)
    error('The TFR window is too long for the requested baseline interval.');
end
condCodes = [61, 62, 63, 64];

% Condition info
condNames  = {'c25', 'c50', 'c75', 'c100'};
condLabels = {'25%', '50%', '75%', '100%'};

% Figure save directories
gcp_root_path = paths.root;
gcp_feature_data_path = paths.features;
if ~exist(gcp_feature_data_path, 'dir')
    gcp_feature_data_path = gcp_root_path;
end
fig_save_dir_ged = fullfile(paths.figures, 'eeg', 'ged');
fig_save_dir_component_selection_base = fullfile(fig_save_dir_ged, 'component_selection');
if ~exist(fig_save_dir_ged, 'dir'), mkdir(fig_save_dir_ged); end
if ~exist(fig_save_dir_component_selection_base, 'dir'), mkdir(fig_save_dir_component_selection_base); end
fig_save_dir_component_selection_root = fig_save_dir_component_selection_base;
[component_parent_dir, component_leaf_dir] = fileparts(fig_save_dir_component_selection_root);
if any(strcmp(component_leaf_dir, subjects))
    % Defensive guard: avoid nested subject folders (e.g., .../601/602)
    fig_save_dir_component_selection_root = component_parent_dir;
end

%% Preallocate storage
trials_powratio = cell(4, nSubj);
trials_peaks = cell(4, nSubj);
trials_outlier_mask_power = cell(4, nSubj);
trials_centroid = cell(4, nSubj); % inprocess centroid figures only
trials_median_centroid = nan(4, nSubj); % inprocess centroid figures only
trials_gamma_power = nan(4, nSubj); % subject inclusion only

all_topos       = cell(1, nSubj);
all_topo_labels = cell(1, nSubj);
all_combined_spectrum = cell(1, nSubj);
all_combined_eigenvalue = nan(1, nSubj);
all_component_selection_stats  = cell(1, nSubj); % inprocess adequacy / TFR meta
all_combined_filter  = cell(1, nSubj); % TFR + Haufe; one filter across contrast
all_haufe_pattern = cell(4, nSubj);
freq_reconstructed_multicomp = cell(4, nSubj);
tfr_cond_trials = cell(4, nSubj); % returned by TFR helper; not saved
tfr_cond_avg = cell(4, nSubj);
ged_filter_meta = cell(1, nSubj); % returned by TFR helper; not saved
subject_runtime_seconds = nan(nSubj, 1);

all_condition_powspctrm = cell(4, nSubj);
freq_powspctrm_unsmoothed = cell(4, nSubj);
all_condition_peak_freq = nan(4, nSubj);
all_condition_peak_power = nan(4, nSubj);

%% Subject loop
for subj = 1:nSubj
    subj_runtime_tic = tic;
    fig_save_dir_component_selection = fullfile(fig_save_dir_component_selection_root, subjects{subj});
    if ~exist(fig_save_dir_component_selection, 'dir'), mkdir(fig_save_dir_component_selection); end
    datapath = fullfile(gcp_feature_data_path, subjects{subj}, 'eeg');
    eeg_data = load(fullfile(datapath, 'dataEEG.mat'), ...
        'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    dataEEG_c25 = eeg_data.dataEEG_c25;
    dataEEG_c50 = eeg_data.dataEEG_c50;
    dataEEG_c75 = eeg_data.dataEEG_c75;
    dataEEG_c100 = eeg_data.dataEEG_c100;
    dataStructs = {dataEEG_c25, dataEEG_c50, dataEEG_c75, dataEEG_c100};
    ref_dat = [];
    for c = 1:4
        if ~isempty(dataStructs{c}) && isstruct(dataStructs{c}) && isfield(dataStructs{c}, 'label')
            ref_dat = dataStructs{c};
            break
        end
    end
    if isempty(ref_dat)
        error('No valid EEG condition data for subject %s.', subjects{subj});
    end
    fsample = ref_dat.fsample;
    nChans = length(ref_dat.label);

    trialIndices = cell(1, 4);
    for c = 1:4
        if isempty(dataStructs{c}) || ~isstruct(dataStructs{c}) || ~isfield(dataStructs{c}, 'trialinfo')
            trialIndices{c} = [];
        else
            trialIndices{c} = find(dataStructs{c}.trialinfo == condCodes(c));
        end
    end

    % Find channels
    occ_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P10|P9)', 'once')), ref_dat.label);
    occ_idx  = find(occ_mask);
    front_mask = cellfun(@(l) ~isempty(regexp(l, '^(Fp|AF|F)', 'once')), ref_dat.label);
    front_idx  = find(front_mask);
    post_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P)', 'once')), ref_dat.label);
    post_idx  = find(post_mask);
    temp_mask = cellfun(@(l) ~isempty(regexp(l, '^(T|TP|FT)', 'once')), ref_dat.label);
    temp_idx = find(temp_mask);
    post_w = zeros(nChans, 1);
    post_w(post_idx) = 1;
    if sum(post_w) > 0
        post_w = post_w / sum(post_w);
    else
        post_w = ones(nChans, 1) / nChans;
    end

    %% Build pooled covariance (stimulus interval)
    clc; close all; fprintf('[GED] Subject GCP%s (%d/%d)\n', subjects{subj}, subj, nSubj);
    rng(random_seed + subj, 'twister');

    covStim  = zeros(nChans);
    covBase  = zeros(nChans);
    covStim_by_cond = cell(1, 4);
    nTrials_total = 0;
    nTrials_per_cond = zeros(1, 4);

    dat_per_cond = cell(1, 4);

    for cond = 1:4
        dat    = dataStructs{cond};
        trlIdx = trialIndices{cond};
        if isempty(trlIdx), continue; end

        dat = select_trials_by_index(dat, trlIdx);
        dat_per_cond{cond} = dat;

        cfg_filt = [];
        cfg_filt.bpfilter   = 'yes';
        cfg_filt.bpfreq     = analysis_freq_range;
        cfg_filt.bpfilttype = 'fir';
        cfg_filt.bpfiltord  = round(3 * fsample / analysis_freq_range(1));
        cfg_filt.feedback   = 'none';
        dat_gamma = ft_preprocessing(cfg_filt, dat);

        nTrl = numel(dat_gamma.trial);
        nTrials_per_cond(cond) = nTrl;
        if nTrl > 0
            cov_base = compute_pooled_covariance_window(dat_gamma, baseline_window);
            cov_stim = compute_pooled_covariance_window(dat_gamma, stim_window);
            covBase = covBase + cov_base * nTrl;
            covStim = covStim + cov_stim * nTrl;
            covStim_by_cond{cond} = cov_stim;
        end
        nTrials_total = nTrials_total + nTrl;
    end

    if nTrials_total < 1
        error('No valid trials available for subject %s after trial selection.', subjects{subj});
    end
    covStim  = covStim / nTrials_total;
    covBase  = covBase / nTrials_total;

    if do_plots
        plot_covariance_matrix_diagnostics( ...
            fig_save_dir_component_selection, subjects{subj}, ref_dat.label, ...
            covBase, covStim, lambda);
    end

    %% Simulated signed occipital template
    template_front_weight = 0.75; % anti-template weight for frontal channels
    template_sigma_occ = 0.12;   % spatial smoothness for occipital template
    template_sigma_front = 0.25; % spatial smoothness for frontal anti-template
    sim_template = zeros(nChans, 1);
    lay_labels = headmodel.layANThead.label;
    lay_pos = headmodel.layANThead.pos;
    chan_pos = nan(nChans, 2);
    for ch = 1:nChans
        li = find(strcmp(lay_labels, ref_dat.label{ch}), 1, 'first');
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

    %% GED + component selection (one filter, trials pooled across contrast)
    dat_search = dat_per_cond;
    lam_w = lambda;
    covStim_reg = (1-lam_w)*covStim + lam_w*mean(diag(covStim))*eye(nChans);
    covBase_reg = (1-lam_w)*covBase + lam_w*mean(diag(covBase))*eye(nChans);

    [W_ged, D_ged] = eig(covStim_reg, covBase_reg);
    [evals_sorted, sortIdx] = sort(real(diag(D_ged)), 'descend');
    W_ged = W_ged(:, sortIdx);

    nSearch = min(ged_search_n, size(W_ged, 2));
    searchFilters = nan(nChans, nSearch);
    searchTopos = nan(nChans, nSearch);
    searchCorrs = nan(nSearch, 1);
    searchPostStrength = nan(nSearch, 1);
    searchFrontStrength = nan(nSearch, 1);
    searchTempStrength = nan(nSearch, 1);
    searchPostFront = nan(nSearch, 1);
    searchPostTemp = nan(nSearch, 1);
    searchTopoPeakFrac = nan(nSearch, 1);
    searchEmgHfSlope = nan(nSearch, 1);
    searchEmgClass = repmat({'unassigned'}, nSearch, 1);
    searchMeanPrSpectrum = nan(nSearch, numel(scan_freqs));

    % Forward model for topoplot and component scoring.
    % Polarity/topo metrics first; proxy spectra are batched across
    % components in one mtmfft call (same per-segment spectra as before).
    for ci = 1:nSearch
        w_ci = W_ged(:, ci);
        topo_ci = covStim_reg * w_ci;
        r_ci = corr(topo_ci, sim_template, 'rows', 'complete');
        if ~isnan(r_ci) && r_ci < 0
            w_ci = -w_ci;
            topo_ci = -topo_ci;
            r_ci = -r_ci;
        end
        if ~isempty(post_idx)
            post_strength = mean(abs(topo_ci(post_idx)));
        else
            post_strength = 0;
        end
        if ~isempty(front_idx)
            front_strength = mean(abs(topo_ci(front_idx)));
            post_front_ci = post_strength / max(front_strength, eps);
        else
            front_strength = 0;
            post_front_ci = Inf;
        end
        if ~isempty(temp_idx)
            temp_strength = mean(abs(topo_ci(temp_idx)));
            post_temp_ci = post_strength / max(temp_strength, eps);
        else
            temp_strength = 0;
            post_temp_ci = Inf;
        end
        topo_abs = abs(topo_ci(:));
        topo_abs_sum = sum(topo_abs(isfinite(topo_abs)));
        if isfinite(topo_abs_sum) && topo_abs_sum > 0
            topo_peak_frac_ci = max(topo_abs(isfinite(topo_abs))) / topo_abs_sum;
        else
            topo_peak_frac_ci = NaN;
        end

        searchFilters(:, ci) = w_ci;
        searchTopos(:, ci) = topo_ci;
        searchCorrs(ci) = r_ci;
        searchPostStrength(ci) = post_strength;
        searchFrontStrength(ci) = front_strength;
        searchTempStrength(ci) = temp_strength;
        searchPostFront(ci) = post_front_ci;
        searchPostTemp(ci) = post_temp_ci;
        searchTopoPeakFrac(ci) = topo_peak_frac_ci;
    end

    cw_prefix = sprintf('[GED] Subject GCP%s (%d/%d)', ...
        subjects{subj}, subj, nSubj);
    proxies = estimate_components_artifact_proxies( ...
        searchFilters, dat_search, stim_window, baseline_window, ...
        fsample, scan_freqs, mtmfft_tapsmofrq_hz, cw_prefix);
    for ci = 1:nSearch
        searchEmgHfSlope(ci) = proxies(ci).hf_slope;
        searchMeanPrSpectrum(ci, :) = proxies(ci).mean_pr_spectrum(:)';
    end

    % Stage-1 gates: SNR, PF, post>front, post>temp, non-focal topo, rising HF EMG
    eval_raw_vec = evals_sorted(1:nSearch);
    post_front_vec = searchPostFront;
    post_temp_vec = searchPostTemp;
    topo_peak_frac_vec = searchTopoPeakFrac;
    emg_hf_slope_vec = searchEmgHfSlope;
    emg_hf_slope_vec(~isfinite(emg_hf_slope_vec)) = 0;
    [powspctrm_form_score_vec, ~] = compute_powspctrm_form_laplacian_score_from_spectra( ...
        searchMeanPrSpectrum, scan_freqs, analysis_freq_range);
    finite_metrics = isfinite(eval_raw_vec) & isfinite(post_front_vec) & ...
        isfinite(post_temp_vec) & isfinite(topo_peak_frac_vec) & ...
        isfinite(powspctrm_form_score_vec);
    pass_eig_gate = finite_metrics & (eval_raw_vec >= min_eigval);
    pass_peak_gate = finite_metrics & (powspctrm_form_score_vec >= min_powspctrm_form);
    pass_post_front_gate = finite_metrics & (post_front_vec > 1);
    pass_post_temp_gate = finite_metrics & (post_temp_vec > 1);
    fail_topo_peak_frac = finite_metrics & (topo_peak_frac_vec > 0.25);
    fail_emg_hf_slope = finite_metrics & (emg_hf_slope_vec > 0);
    for ci = 1:nSearch
        if fail_emg_hf_slope(ci) || fail_topo_peak_frac(ci)
            searchEmgClass{ci} = 'EMG';
        elseif ~(post_front_vec(ci) > 1)
            searchEmgClass{ci} = 'frontal';
        elseif ~(post_temp_vec(ci) > 1)
            searchEmgClass{ci} = 'temporal';
        else
            searchEmgClass{ci} = 'posterior';
        end
    end
    eligible = pass_eig_gate & pass_peak_gate & pass_post_front_gate & ...
        pass_post_temp_gate & ~fail_topo_peak_frac & ~fail_emg_hf_slope;
    no_threshold_match = ~any(eligible);
    selection_pool_mask = eligible;
    searchScores = eval_raw_vec;
    searchScores(~finite_metrics) = -Inf;
    searchScores(~selection_pool_mask) = -Inf;
    [bestScore, bestIdx] = max(searchScores);
    if isempty(bestIdx) || isnan(bestScore)
        bestIdx = 1;
        bestScore = NaN;
    end

    combined_idx = find(selection_pool_mask);
    if isempty(combined_idx)
        combined_weights = [];
    else
        [~, combined_ord] = sort(eval_raw_vec(combined_idx), 'descend');
        combined_idx = combined_idx(combined_ord);

        combined_weights = eval_raw_vec(combined_idx)';
        combined_weights(~isfinite(combined_weights) | combined_weights <= 0) = 0;
        if sum(combined_weights) <= 0
            combined_weights = ones(1, numel(combined_idx));
        end
        combined_weights = combined_weights / sum(combined_weights);
    end

    selected_idx = combined_idx;
    selected_weights = combined_weights;

    if isempty(selected_idx)
        bestIdx = NaN;
        bestScore = NaN;
        bestCorr = NaN;
        bestPost = NaN;
        bestFront = NaN;
        bestPostFront = NaN;
        bestPostTemp = NaN;
        bestTopoPeakFrac = NaN;
    else
        bestIdx = selected_idx(1);
        bestScore = searchScores(bestIdx);
        bestCorr = searchCorrs(bestIdx);
        bestPost = searchPostStrength(bestIdx);
        bestFront = searchFrontStrength(bestIdx);
        bestPostFront = post_front_vec(bestIdx);
        bestPostTemp = post_temp_vec(bestIdx);
        bestTopoPeakFrac = topo_peak_frac_vec(bestIdx);
    end
    all_topo_labels{subj} = ref_dat.label;
    all_component_selection_stats{subj} = struct( ...
        'selection_mode', 'fixed_preregistered_weighted', ...
        'selected_idx', selected_idx, ...
        'n_selected_ged_components', numel(selected_idx), ...
        'selected_weights', selected_weights, ...
        'best_idx', bestIdx, ...
        'best_score', bestScore, ...
        'best_corr', bestCorr, ...
        'best_post_front', bestPostFront, ...
        'best_post_temp', bestPostTemp, ...
        'best_topo_peak_frac', bestTopoPeakFrac, ...
        'best_front', bestFront, ...
        'best_post', bestPost, ...
        'post_front', post_front_vec, ...
        'post_temp', post_temp_vec, ...
        'topo_peak_frac', topo_peak_frac_vec, ...
        'emg_hf_slope', emg_hf_slope_vec, ...
        'emg_class', {searchEmgClass}, ...
        'eligible', eligible, ...
        'no_threshold_match', no_threshold_match);
    w_combined = selected_weights(:)';

    if isempty(selected_idx)
        W_combined = [];
        topo_temp = nan(nChans, 1);
    else
        W_combined = searchFilters(:, selected_idx);
        topo_temp = searchTopos(:, selected_idx) * w_combined(:);
    end
    all_topos{subj} = topo_temp;
    [sel_idx_spec, sel_w_spec] = sanitize_selected_components( ...
        selected_idx, w_combined, size(searchMeanPrSpectrum, 1));
    if isempty(sel_idx_spec) || isempty(searchMeanPrSpectrum)
        all_combined_spectrum{subj} = [];
        all_combined_eigenvalue(subj) = NaN;
    else
        all_combined_spectrum{subj} = sel_w_spec(:)' * searchMeanPrSpectrum(sel_idx_spec, :);
        evals_sel = evals_sorted(sel_idx_spec);
        evals_sel = evals_sel(:);
        evals_sel(~isfinite(evals_sel)) = NaN;
        all_combined_eigenvalue(subj) = sum(sel_w_spec(:) .* evals_sel);
    end

    cfg_topo = [];
    cfg_topo.layout    = headmodel.layANThead;
    cfg_topo.comment   = 'no';
    cfg_topo.marker    = 'off';
    cfg_topo.style     = 'straight';
    cfg_topo.gridscale = topo_gridscale;
    cfg_topo.zlim      = 'maxabs';
    cfg_topo.colormap  = '*RdBu';
    cfg_topo.figure    = 'gcf';
    if do_plots
        plot_selected_components( ...
            fig_save_dir_component_selection, subjects{subj}, scan_freqs, searchTopos, ...
            searchMeanPrSpectrum, evals_sorted(1:numel(eligible)), ...
            searchEmgClass, ...
            eligible, ...
            post_front_vec, post_temp_vec, topo_peak_frac_vec, emg_hf_slope_vec, ...
            cfg_topo, all_topo_labels{subj}, powspctrm_form_score_vec, ...
            selected_idx);
        plot_combined_topo_spectra( ...
            fig_save_dir_component_selection, subjects{subj}, scan_freqs, cfg_topo, all_topo_labels{subj}, ...
            searchTopos, searchMeanPrSpectrum, selected_idx, w_combined, ...
            analysis_freq_range);
    end

    adequate = ~isempty(selected_idx);
    W_comb = W_combined;
    w_comb = w_combined;
    if ~adequate
        W_comb = zeros(nChans, 0);
        w_comb = [];
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
    W_comb = normalize_filters_to_noise_metric(W_comb, covBase);
    all_combined_filter{subj} = build_combined_filter_vector(W_comb, w_comb);

    %% Condition-specific Haufe patterns and multicomponent reconstruction
    if adequate && ~isempty(W_comb)
        W_selected = W_comb;
        component_cov_pool = W_selected' * covStim * W_selected;
        reconstruction_patterns = covStim * W_selected * pinv(component_cov_pool);
        combined_filter = all_combined_filter{subj};
        for cond = 1:4
            cov_cond = covStim_by_cond{cond};
            dat_cond = dat_per_cond{cond};
            if isempty(cov_cond) || isempty(dat_cond)
                continue;
            end

            combined_variance_cond = combined_filter' * cov_cond * combined_filter;
            if isfinite(combined_variance_cond) && combined_variance_cond > eps
                all_haufe_pattern{cond, subj} = ...
                    (cov_cond * combined_filter) / combined_variance_cond;
            end

            nTrl_cond = numel(dat_cond.trial);
            dat_reconstructed_base = [];
            dat_reconstructed_base.label = dat_cond.label;
            dat_reconstructed_base.fsample = dat_cond.fsample;
            dat_reconstructed_base.trial = cell(1, nTrl_cond);
            dat_reconstructed_base.time = cell(1, nTrl_cond);
            dat_reconstructed_stim = dat_reconstructed_base;
            for trl = 1:nTrl_cond
                x = double(dat_cond.trial{trl});
                t = dat_cond.time{trl};
                idx_base = t >= baseline_window(1) & t <= baseline_window(2);
                idx_stim = t >= stim_window(1) & t <= stim_window(2);
                dat_reconstructed_base.trial{trl} = ...
                    reconstruction_patterns * (W_selected' * x(:, idx_base));
                dat_reconstructed_base.time{trl} = t(idx_base);
                dat_reconstructed_stim.trial{trl} = ...
                    reconstruction_patterns * (W_selected' * x(:, idx_stim));
                dat_reconstructed_stim.time{trl} = t(idx_stim);
            end

            cfg_freq = [];
            cfg_freq.method = 'mtmfft';
            cfg_freq.output = 'pow';
            cfg_freq.taper = 'dpss';
            cfg_freq.foi = scan_freqs;
            cfg_freq.tapsmofrq = mtmfft_tapsmofrq_hz;
            cfg_freq.pad = 'nextpow2';
            cfg_freq.keeptrials = 'yes';
            cfg_freq.feedback = 'none';
            freq_base = ft_freqanalysis(cfg_freq, dat_reconstructed_base);
            freq_stim = ft_freqanalysis(cfg_freq, dat_reconstructed_stim);

            freq_ratio_trials = freq_stim;
            freq_ratio_trials.powspctrm = 10 * log10( ...
                max(double(freq_stim.powspctrm), eps) ./ ...
                max(double(freq_base.powspctrm), eps));
            cfg_desc = [];
            cfg_desc.keeptrials = 'no';
            freq_ratio = ft_freqdescriptives(cfg_desc, freq_ratio_trials);
            freq_reconstructed_multicomp{cond, subj} = freq_ratio;
        end
    end

    %% Per-condition trial-level spectral scanning
    subj_powratio = cell(1, 4);
    subj_peaks = cell(1, 4);
    subj_centroid = cell(1, 4);
    subj_condition_avg = cell(1, 4);
    subj_condition_peak = nan(1, 4);
    for cond = 1:4
        dat = dat_per_cond{cond};
        if isempty(dat)
            continue;
        end

        nTrl = length(dat.trial);
        powratio_methods = nan(1, nTrl, nFreqs);
        nSearch_sel = size(W_comb, 2);
        powratio_components       = nan(nSearch_sel, nTrl, nFreqs);
        unstable_freq_counts = zeros(nTrl, 1);
        valid_freq_counts = zeros(nTrl, 1);

        % Baseline quality gate computed once per trial (not per frequency).
        baseline_power_raw = nan(nTrl, 1);
        baseline_power_comb = nan(nTrl, 1);
        trial_cache = cell(nTrl, 1);
        has_base = false(nTrl, 1);
        has_stim = false(nTrl, 1);
        for trl = 1:nTrl
            x = double(dat.trial{trl});
            t = dat.time{trl};
            idx_base = t >= baseline_window(1) & t <= baseline_window(2);
            idx_stim = t >= stim_window(1) & t <= stim_window(2);
            x_base = x(:, idx_base);
            x_stim = x(:, idx_stim);
            trial_cache{trl} = struct('x_base', x_base, 'x_stim', x_stim);
            has_base(trl) = ~isempty(x_base);
            has_stim(trl) = ~isempty(x_stim);
            if has_base(trl)
                pow_base_chan = mean(x_base.^2, 2);
                baseline_power_raw(trl) = sum(post_w(:) .* pow_base_chan(:));
                if adequate && ~isempty(W_comb)
                    x_base_proj = W_comb' * x_base;
                    baseline_power_comb(trl) = mean(x_base_proj(:).^2);
                end
            end
        end

        baseline_outlier_mad_mult = 3.5;    % robust log-power cutoff (MAD units) for baseline outliers
        ratio_floor_prctile = 20;           % robust baseline-power percentile used as floor anchor
        ratio_floor_frac = 0.25;            % floor is this fraction of the robust baseline-power anchor
        instability_near_floor_mult = 1.5;  % mark as near-floor when baseline <= this multiple of floor
        instability_trial_freq_frac_thr = 0.35; % exclude trial when unstable at >= this frequency fraction

        bad_base_raw = flag_unreliable_baseline_trials( ...
            baseline_power_raw, baseline_outlier_mad_mult);
        bad_base = bad_base_raw;
        if adequate
            bad_base = bad_base | flag_unreliable_baseline_trials( ...
                baseline_power_comb, baseline_outlier_mad_mult);
        end
        [base_floor, ~] = compute_baseline_floor_stats(baseline_power_comb, ratio_floor_prctile, ratio_floor_frac);

        if adequate
            trial_mask = has_base & has_stim & ~bad_base;
            [ratio_cube, ratio_trials_combined, near_floor_count, ...
                near_floor_count_combined, valid_freq_counts_combined] = ...
                compute_scan_ratio_for_window_and_combined_batch( ...
                trial_cache, W_comb, w_comb, 'x_stim', trial_mask, ...
                fsample, scan_freqs, mtmfft_tapsmofrq_hz, base_floor, instability_near_floor_mult);
            powratio_components = ratio_cube;
            for trl = 1:nTrl
                if ~isempty(ratio_trials_combined)
                    powratio_methods(1, trl, :) = ratio_trials_combined(trl, :);
                end
                if trl <= numel(valid_freq_counts_combined)
                    valid_freq_counts(trl) = valid_freq_counts_combined(trl);
                end
                if trl <= numel(near_floor_count_combined) && valid_freq_counts(trl) > 0
                    unstable_freq_counts(trl) = near_floor_count_combined(trl);
                elseif any(isfinite(powratio_components(:, trl, :)))
                    unstable_freq_counts(trl) = near_floor_count(trl);
                end
            end
        end
        unstable_trial_frac = unstable_freq_counts ./ max(valid_freq_counts, 1);
        trial_unstable = unstable_trial_frac >= instability_trial_freq_frac_thr;
        if any(trial_unstable)
            powratio_methods(:, trial_unstable, :) = NaN;
        end

        powratio_trials = squeeze(powratio_methods(1, :, :));
        trials_powratio{cond, subj} = powratio_trials;
        subj_powratio{cond} = powratio_trials;
        %% Per-trial peak detection
        trial_metric_outlier_iqr_mult = 1.5; % outlier threshold in IQR units around Q1/Q3
        [trl_peaks, trial_peak_power, trl_centroid] = ...
            compute_trial_peak_metrics_from_powratio( ...
            powratio_trials, scan_freqs, true(size(scan_freqs)), ...
            5, 5);

        trials_peaks{cond, subj} = trl_peaks;
        trials_centroid{cond, subj} = trl_centroid;
        subj_peaks{cond} = trl_peaks;
        subj_centroid{cond} = trl_centroid;

        % Trial-level metric outlier rejection (subject-condition specific).
        [outlier_mask_freq, ~] = detect_trial_metric_outliers_iqr( ...
            trl_peaks, trial_metric_outlier_iqr_mult);
        [outlier_mask_power, ~] = detect_trial_metric_outliers_iqr( ...
            trial_peak_power, trial_metric_outlier_iqr_mult);
        trl_peaks(outlier_mask_freq) = NaN;
        trial_peak_power(outlier_mask_power) = NaN;
        trials_peaks{cond, subj} = trl_peaks;
        trials_outlier_mask_power{cond, subj} = outlier_mask_power;
        subj_peaks{cond} = trl_peaks;

        % Condition-level spectra and peak metrics from trial-averaged spectra.
        % All valid (non-unstable) trial spectra enter the average; IQR on
        % trial peak metrics is retained only for exploratory trial summaries.
        cond_avg = compute_condition_average_powratio_ft(powratio_trials, scan_freqs);
        all_condition_powspctrm{cond, subj} = cond_avg;
        freq_powspctrm_unsmoothed{cond, subj} = ged_powcurve_to_freq_ft( ...
            cond_avg, scan_freqs, ref_dat);
        subj_condition_avg{cond} = cond_avg;

        [peak_hz, peak_power] = pick_tallest_peak(cond_avg, scan_freqs, 0, 5);
        all_condition_peak_freq(cond, subj) = peak_hz;
        all_condition_peak_power(cond, subj) = peak_power;
        subj_condition_peak(cond) = peak_hz;

        valid_c = isfinite(trl_centroid);
        trials_median_centroid(cond, subj) = median(trl_centroid(valid_c));

        % Peak power: highest dB value in the trial spectrum (inclusion).
        trials_gamma_power(cond, subj) = robust_trial_mean(trial_peak_power);

    end % condition loop

    if do_tfr && adequate && ~isempty(all_combined_filter{subj})
        fprintf('[GED TFR] Subject GCP%s (%d/%d)\n', subjects{subj}, subj, nSubj);
        stat_sel = all_component_selection_stats{subj};
        if isempty(stat_sel)
            stat_sel = struct();
        end
        [tfr_cond_trials(:, subj), tfr_cond_avg(:, subj), ged_filter_meta{subj}] = ...
            compute_ged_tfr_subject( ...
            dat_per_cond, all_combined_filter{subj}, all_topo_labels{subj}, ...
            subjects{subj}, subj, stat_sel, ...
            baseline_window, tfr_foi, tfr_toi, tfr_win_sec, tfr_tapsmofrq, condCodes);
    end

    % PER-SUBJECT FIGURES
    if do_plots
    close all
    cmap_div = interp1([0 0.5 1], ...
        [0.17 0.27 0.53; 0.97 0.97 0.97; 0.70 0.09 0.17], linspace(0,1,256));

    cfg_topo = [];
    cfg_topo.layout    = headmodel.layANThead;
    cfg_topo.comment   = 'no';
    cfg_topo.marker    = 'off';
    cfg_topo.style     = 'straight';
    cfg_topo.gridscale = topo_gridscale;
    cfg_topo.zlim      = 'maxabs';
    cfg_topo.colormap  = '*RdBu';
    cfg_topo.figure    = 'gcf';

    pr_source = subj_powratio;
    peaks_source = subj_peaks;
    centroid_source = subj_centroid;
    condavg_source = subj_condition_avg;
    condpeak_source = subj_condition_peak;
    topo_mat = searchTopos;
    selected_idx_plot = selected_idx;
    selected_w_plot = w_combined;
    eigvals_plot = evals_sorted;

        pr_raw_mats = cell(1, 4);
        row1_clim = ones(1, 4);
        for cond = 1:4
            pr_raw_mats{cond} = pr_source{cond};
            if ~isempty(pr_raw_mats{cond})
                cond_vals = abs(pr_raw_mats{cond}(:));
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

        fig = figure('Position', [0 0 1512 982], 'Color', 'w');
        sgtitle(sprintf('Trial-Level GED: Subject %s', subjects{subj}), ...
            'FontSize', 18, 'FontWeight', 'bold');

        % --- Row 1: Heatmap of trial-level spectra ---
        for cond = 1:4
            subplot(2, 4, cond);
            if ~isempty(pr_raw_mats{cond})
                hold on;
                imagesc(scan_freqs, 1:size(pr_raw_mats{cond},1), pr_raw_mats{cond});
                colormap(gca, cmap_div);
                caxis([-row1_clim(cond) row1_clim(cond)]);
                cb = colorbar; cb.FontSize = 8;
                ctd = centroid_source{cond};
                if ~isempty(ctd)
                    valid_ctd = ~isnan(ctd);
                    tr_idx = find(valid_ctd);
                    if ~isempty(tr_idx)
                        scatter(ctd(valid_ctd), tr_idx, 28, 'k', 'filled', ...
                            'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.75, ...
                            'MarkerFaceAlpha', 0.85);
                    end
                end
                xlabel('Freq [Hz]'); ylabel('Trial');
                set(gca, 'YDir', 'normal');
            end
            title(sprintf('%s Raw', condLabels{cond}), 'FontSize', 11);
            set(gca, 'FontSize', 10); xlim([30 90]); ylim([0 200]); box on;
        end

        % --- Row 2: Topoplot + histogram + combined spectra ---
        subplot(2, 4, 5);
        if ~isempty(topo_mat) && ~isempty(selected_idx_plot)
            topo_data = [];
            topo_data.label  = all_topo_labels{subj};
            w_plot_common = selected_w_plot(:);
            if numel(w_plot_common) ~= numel(selected_idx_plot) || ~any(isfinite(w_plot_common))
                w_plot_common = ones(numel(selected_idx_plot), 1);
            end
            w_plot_common(~isfinite(w_plot_common) | w_plot_common <= 0) = 0;
            if sum(w_plot_common) <= 0
                w_plot_common = ones(numel(selected_idx_plot), 1);
            end
            w_plot_common = w_plot_common / sum(w_plot_common);
            topo_plot_common = topo_mat(:, selected_idx_plot) * w_plot_common;
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
                topo_clim_common = prctile(topo_abs_common, 99.9);
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
            n_sel_show = numel(selected_idx_plot);
            lambda_show = NaN;
            if ~isempty(eigvals_plot) && ~isempty(selected_idx_plot)
                idx_show = selected_idx_plot(1);
                if idx_show >= 1 && idx_show <= numel(eigvals_plot)
                    lambda_show = eigvals_plot(idx_show);
                end
            end
            title(sprintf('Weighted GED (%d comps, \\lambda=%.2f)', n_sel_show, lambda_show), 'FontSize', 11);
        end

        subplot(2, 4, [6 7]); hold on;
        edges = 30:2:90;
        hist_mat = zeros(4, length(edges)-1);
        for cond = 1:4
            tpk = peaks_source{cond};
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
        title('Trial-Level Peak Distribution', 'FontSize', 12);
        legend(bh, condLabels, 'FontSize', 10, 'Location', 'best');
        set(gca, 'FontSize', 11); xlim([30 90]);  box on;

        subplot(2, 4, 8); hold on;
        for cond = 1:4
            avg_curve = condavg_source{cond};
            if isempty(avg_curve) || numel(avg_curve) ~= numel(scan_freqs)
                continue;
            end
            valid_curve = isfinite(avg_curve) & isfinite(scan_freqs);
            if sum(valid_curve) < 3
                continue;
            end
            x_curve = scan_freqs(valid_curve);
            y_curve = avg_curve(valid_curve);
            plot(x_curve, y_curve, '-', 'Color', colors(cond, :), 'LineWidth', 2.2, ...
                'DisplayName', condLabels{cond});
            peak_hz = condpeak_source(cond);
            if isfinite(peak_hz)
                xline(peak_hz, ':', 'Color', colors(cond, :), 'LineWidth', 1.2, ...
                    'HandleVisibility', 'off');
            end
        end
        yline(0, 'k--', 'LineWidth', 0.7, 'HandleVisibility', 'off');
        xlim([scan_freqs(1), scan_freqs(end)]);
        xlabel('Freq [Hz]');
        ylabel('Power [dB]');
        title('Condition-Averaged Spectra', 'FontSize', 11);
        set(gca, 'FontSize', 10, 'Box', 'on');
        peak_text_y = 0.95;
        peak_text_step = 0.08;
        for cond = 1:4
            peak_hz = condpeak_source(cond);
            if ~isfinite(peak_hz)
                continue;
            end
            text(0.98, peak_text_y - (cond - 1) * peak_text_step, sprintf('%.0f Hz', peak_hz), ...
                'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                'Color', colors(cond, :), 'FontSize', 10, 'FontWeight', 'bold');
        end
        legend('Location', 'southwest', 'FontSize', 9, 'Box', 'off');

        save_figure_png(fig, fullfile(fig_save_dir_component_selection, ...
            sprintf('GCP_eeg_GED_subj%s_trials_overview.png', subjects{subj})));
    end % do_plots (per-subject figures)
    subject_runtime_seconds(subj) = toc(subj_runtime_tic);
end % subject loop

% CENTROID METRIC: Subject/group summaries and concordance
if do_plots
fig_cent = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Trial-Level Gamma Centroid', ...
    'FontSize', 18, 'FontWeight', 'bold');

subplot(1, 2, 1); hold on;
for s = 1:nSubj
    yc = trials_median_centroid(:, s);
    if sum(~isnan(yc)) >= 2
        plot(1:4, yc, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    end
end
for c = 1:4
    vals = trials_median_centroid(c, :);
    vals = vals(~isnan(vals));
    if ~isempty(vals)
        xj = c + (rand(size(vals)) - 0.5) * 0.12;
        scatter(xj, vals, 110, colors(c,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
    end
end
med_c = nanmedian(trials_median_centroid, 2);
mad_c = nan(4, 1);
for c = 1:4
    mad_c(c) = robust_mad(trials_median_centroid(c, :));
end
errorbar(1:4, med_c, mad_c, 'k', 'LineWidth', 2, 'CapSize', 10);
plot(1:4, med_c, 'k-', 'LineWidth', 2.5);
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
        tc = trials_centroid{c, s};
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
ylim(analysis_freq_range);
ylabel('Centroid Frequency [Hz]');
title('All trials pooled', 'FontSize', 14, 'FontWeight', 'bold');
save_figure_png(fig_cent, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_centroid_summary.png'));

%% Standalone condition-separation metrics (combined GED)
close all
fig_cond_slope = figure('Position', [0 0 1512 982], 'Color', 'w');

slope_post = compute_condition_separation_from_matrix(all_condition_peak_freq);
delta_post = all_condition_peak_freq(4, :) - all_condition_peak_freq(1, :);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Panel 1: condition slope
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
ylim(compute_symmetric_ylim(slope_vals, 0.15, 0.4));
set(gca, 'XTick', 1, 'XTickLabel', {'Combined GED'}, ...
    'FontSize', 16, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
ylabel('Slope across contrast conditions [Hz/condition]', 'FontSize', 18, 'FontWeight', 'bold');
title('Contrast Condition Slope', 'FontSize', 20, 'FontWeight', 'bold');

% Panel 2: median shift (100%-25%)
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
ylim(compute_symmetric_ylim(delta_vals, 0.15, 1));
set(gca, 'XTick', 1, 'XTickLabel', {'Combined GED'}, ...
    'FontSize', 16, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
ylabel('\Delta median (100% - 25%) [Hz]', 'FontSize', 18, 'FontWeight', 'bold');
title('Median Frequency Shift (100% - 25%)', 'FontSize', 20, 'FontWeight', 'bold');

save_figure_png(fig_cond_slope, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_condition_slope.png'));

%% Mean gamma frequency shift bar plot
close all
fig_cond_shift_bar_freq = figure('Position', [0 0 1512 982], 'Color', 'w');
valid_delta = isfinite(delta_post);
subj_idx = find(valid_delta);
delta_vals = delta_post(valid_delta);
bar(subj_idx, delta_vals, 0.75, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'k');
hold on;
yline(0, 'k--', 'LineWidth', 1.2);
xlim([0.5 max(subj_idx) + 0.5]);
ylim(compute_symmetric_ylim(delta_vals, 0.15, 1));
xticks(subj_idx);
if exist('subjects', 'var') == 1 && numel(subjects) >= max(subj_idx) && numel(subj_idx) <= 35
    xticklabels(subjects(subj_idx));
end
xlabel('Subject', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('\Delta Frequency Shift (100% - 25%) [Hz]', 'FontSize', 18, 'FontWeight', 'bold');
title('Gamma Frequency Shift', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
save_figure_png(fig_cond_shift_bar_freq, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_bar_GammaFreq.png'));

%% Condition-shift figure: normalized gamma frequency trajectories
fig_condition_shift_freq = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

dat_freq = all_condition_peak_freq;  % [condition x subject], peak from condition-averaged spectra
dat_freq_shift = dat_freq - dat_freq(1, :);  % Anchor each subject at 25% condition

for s = 1:nSubj
    y_subj = dat_freq_shift(:, s);
    valid_subj = isfinite(y_subj);
    if sum(valid_subj) >= 2
        plot(find(valid_subj), y_subj(valid_subj), '-o', ...
            'Color', [0.75 0.75 0.75], 'LineWidth', 1.2, ...
            'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerSize', 6);
    end
end

med_freq_shift = nanmedian(dat_freq_shift, 2);
mad_freq_shift = nan(4, 1);
for c = 1:4
    mad_freq_shift(c) = robust_mad(dat_freq_shift(c, :));
end
errorbar(1:4, med_freq_shift, mad_freq_shift, '-o', ...
    'Color', colors(4, :), 'LineWidth', 3, 'CapSize', 10, ...
    'MarkerFaceColor', colors(4, :), 'MarkerSize', 8);

yline(0, 'k--', 'LineWidth', 1.5);
set(gca, 'XTick', 1:4, 'XTickLabel', strcat(condLabels, ' Contrast'), ...
    'FontSize', 18, 'Box', 'off');
xlim([0.5 4.5]);
ylim(compute_symmetric_ylim(dat_freq_shift(:), 0.15, 1));
xlabel('Contrast condition');
ylabel('\Delta Gamma Frequency [Hz]');
title('Gamma Peak Frequency Shift', ...
    'FontSize', 24, 'FontWeight', 'bold');

cond_shift_freq_path = fullfile(fig_save_dir_ged, 'GCP_eeg_GED_condition_shift_frequency.png');
save_figure_png(fig_condition_shift_freq, cond_shift_freq_path);

%% Frequency figure: gamma frequency over contrast
fig_main_gamma = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
plot_gamma_condition_panel(all_condition_peak_freq, condLabels, colors, nSubj);
ylabel('Gamma Frequency [Hz]');
title('Gamma Frequency', 'FontWeight', 'bold');
save_figure_png(fig_main_gamma, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_freq.png'));

%% Power figure: gamma power over contrast
fig_main_power = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
plot_gamma_condition_panel(all_condition_peak_power, condLabels, colors, nSubj);
ylabel('Gamma Peak Power [dB]');
title('Gamma Peak Power', 'FontWeight', 'bold');
save_figure_png(fig_main_power, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_power.png'));

%% Condition-shift figure: normalized peak power trajectories
fig_condition_shift_power = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

dat_power = all_condition_peak_power;
dat_power_shift = dat_power - dat_power(1, :);  % Anchor each subject at 25% condition

for s = 1:nSubj
    y_subj = dat_power_shift(:, s);
    valid_subj = isfinite(y_subj);
    if sum(valid_subj) >= 2
        plot(find(valid_subj), y_subj(valid_subj), '-o', ...
            'Color', [0.75 0.75 0.75], 'LineWidth', 1.2, ...
            'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerSize', 6);
    end
end

med_power_shift = nanmedian(dat_power_shift, 2);
mad_power_shift = nan(4, 1);
for c = 1:4
    mad_power_shift(c) = robust_mad(dat_power_shift(c, :));
end
errorbar(1:4, med_power_shift, mad_power_shift, '-o', ...
    'Color', colors(4, :), 'LineWidth', 3, 'CapSize', 10, ...
    'MarkerFaceColor', colors(4, :), 'MarkerSize', 8);

yline(0, 'k--', 'LineWidth', 1.5);
set(gca, 'XTick', 1:4, 'XTickLabel', strcat(condLabels, ' Contrast'), ...
    'FontSize', 18, 'Box', 'off');
xlim([0.5 4.5]);
ylim(compute_symmetric_ylim(dat_power_shift(:), 0.15, 0.2));
xlabel('Contrast condition');
ylabel('\Delta Peak Power [dB]');
title('Gamma Peak Power Shift', ...
    'FontSize', 24, 'FontWeight', 'bold');

cond_shift_power_path = fullfile(fig_save_dir_ged, 'GCP_eeg_GED_condition_shift_power.png');
save_figure_png(fig_condition_shift_power, cond_shift_power_path);

%% All-subjects combined topography + spectrum overview
cfg_topo_all = [];
cfg_topo_all.layout    = headmodel.layANThead;
cfg_topo_all.comment   = 'no';
cfg_topo_all.marker    = 'off';
cfg_topo_all.style     = 'straight';
cfg_topo_all.gridscale = topo_gridscale;
cfg_topo_all.zlim      = 'maxabs';
cfg_topo_all.colormap  = '*RdBu';
cfg_topo_all.figure    = 'gcf';
plot_all_subjects_combined_topo_spectra( ...
    fig_save_dir_component_selection_root, subjects, scan_freqs, analysis_freq_range, cfg_topo_all, ...
    all_topo_labels, all_topos, all_combined_spectrum, all_combined_eigenvalue);

plot_subjects_combined_topo_spectra_grid( ...
    fig_save_dir_component_selection_root, subjects, 1:min(5, nSubj), scan_freqs, analysis_freq_range, cfg_topo_all, ...
    all_topo_labels, all_topos, all_combined_spectrum, '01-05');
if nSubj >= 6
    plot_subjects_combined_topo_spectra_grid( ...
        fig_save_dir_component_selection_root, subjects, 6:min(10, nSubj), scan_freqs, analysis_freq_range, cfg_topo_all, ...
        all_topo_labels, all_topos, all_combined_spectrum, '06-10');
end
end % do_plots (group figures)

%% Save results
save_path = fullfile(gcp_root_path, 'data', 'features', 'GCP_eeg_GED.mat');
save(save_path, ...
    'trials_powratio', ...
    'trials_peaks', ...
    'trials_outlier_mask_power', ...
    'all_topos', 'all_topo_labels', ...
    'all_combined_filter', ...
    'all_combined_spectrum', ...
    'all_combined_eigenvalue', ...
    'all_haufe_pattern', 'freq_reconstructed_multicomp', ...
    'all_condition_powspctrm', ...
    'all_condition_peak_freq', ...
    'all_condition_peak_power', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames', ...
    'baseline_window', 'stim_window');

clc
powspctrm_save_path = fullfile(gcp_root_path, 'data', 'features', 'GCP_eeg_powspctrm_GED.mat');
% freq_* cells: one FieldTrip freq per (cond, subj), chan_freq, for ft_freqgrandaverage like AOC powl2{subj}.
% Do not add topolabel to those structs: FieldTrip ft_datatype sets iscomp if topolabel is present.
save(powspctrm_save_path, ...
    'freq_powspctrm_unsmoothed', ...
    'all_condition_peak_freq', ...
    'scan_freqs', 'condLabels', 'subjects');
fprintf('[GED] Subject-level freq (GED spectra) saved to: %s\n', powspctrm_save_path);

%% Save GED analysis cohort (subjects with valid gamma power)
SubjID = str2double(string(subjects(:)));
Include = any(isfinite(trials_gamma_power), 1)';
subject_inclusion = table(SubjID, Include, 'VariableNames', {'SubjID', 'Include'});
save(fullfile(paths.controls, 'GCP_subject_inclusion.mat'), 'subject_inclusion', '-v7.3');

%% GED-projected TFR (computed in the subject loop when do_tfr is true)
if do_tfr
    tfr_save_path = fullfile(gcp_feature_data_path, 'GCP_eeg_GED_TFR.mat');
    save(tfr_save_path, ...
        'tfr_cond_avg', ...
        'subjects', 'condNames', 'condLabels', ...
        'baseline_window', 'tfr_baseline_window', 'stim_window', ...
        'tfr_foi', 'tfr_toi', 'tfr_win_sec', 'tfr_tapsmofrq', ...
        '-v7.3');
    fprintf('[GED TFR] Saved GED-TFR data to: %s\n', tfr_save_path);
end

fprintf('[GED] DONE!\n');
fprintf('[GED] Feature extraction results saved to: %s\n', save_path);
for si = 1:nSubj
    if isfinite(subject_runtime_seconds(si))
        fprintf('[GED] Runtime Subject %s: %s\n', subjects{si}, format_runtime_hhmmss(subject_runtime_seconds(si)));
    else
        fprintf('[GED] Runtime Subject %s: n/a\n', subjects{si});
    end
end
fprintf('[GED] Runtime TOTAL: %s\n', format_runtime_hhmmss(toc(total_runtime_tic)));

