%% GCP Gamma Peak Frequency and Power with Generalized Eigendecomposition (GED)
%
% Broadband GED and component selection (per subject)
%   - Pool trials across conditions and compute gamma-band covariances
%      (30-90 Hz FIR) for three windows (full/early/late) and baseline.
%   - Solve window-specific GED (S_stim * w = lambda * S_base * w) with
%      regularisation and rank candidate components by eigenvalue.
%   - Retain candidates that pass SNR (lambda), spectral peak-form (PF),
%      occipital>frontal dominance (occdom), and transparent EMG gates
%      (temporal dominance or rising HF slope). Combine up to 5 eligible
%      components with eigenvalue-proportional weights.
%
% Trial-level spectral scanning (per subject, condition, trial)
%   - Project each trial to the combined GED component space and compute
%      spectrum-based power scans on a 30-90 Hz grid in dB.
%   - Flag numerical-instability cases (near-floor baseline power across
%      many frequencies/components) and exclude unstable trials automatically.
%   - Detect per-trial peak frequency and define peak power as the mean power
%      within peak frequency +/- 5 Hz.
%
% Outputs
%   - Trial-level peak frequency/power cell arrays (trials_peaks, ...).
%   - Condition-averaged spectral peaks (all_condition_peak_freq/power_*),
%     used by subject-level master matrix, boxplots, and rainclouds.
%   - Detectability, trial CV, peak power, and condition separation.
%   - Subject diagnostics (component selection, rejection reasons,
%     topographies, spectra), group summary figures.
%   - Optional GED-projected condition TFRs (toggle do_tfr below).
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
do_tfr = true; % set false to skip GED-projected TFR feature extraction

% Time windows
baseline_window = [-1.5, -0.5];
full_window = [0, 2.0];
early_window = [0, 1.0];
late_window = [1.0, 2.0];

% Gamma analysis (frequency grid and FieldTrip mtmfft multitaper bandwidth)
analysis_freq_range = [30 90];
scan_freq_step_hz = 1; % frequency grid step (Hz); powratio_*_freq_smooth_bins count bins on this grid
scan_freqs = analysis_freq_range(1):scan_freq_step_hz:analysis_freq_range(2);
nFreqs = length(scan_freqs);
mtmfft_tapsmofrq_hz = 3; % FieldTrip cfg.tapsmofrq for mtmfft (Hz)

% GED
lambda = 0.05;              % regularization
ged_search_n = 10;          % search first N GED components
min_eigval = 1.05;          % minimum GED eigenvalue (lambda >= 1.05)
min_powspctrm_form = 0.75;  % minimum PF (powspctrm-form) score for candidate eligibility
max_components_to_combine = 5; % top-K cap for lambda-weighted combination
random_seed = 123;
powratio_trial_freq_smooth_bins = 5;      % movmean length (frequency bins) on per-trial powratio for peak/centroid
powratio_condition_freq_smooth_bins = 1;  % movmean length on condition-mean powratio before condition-level peaks/plots
peak_power_halfwidth_hz = 5;  % peak power = mean power within peak_hz +/-

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
trials_powratio     = cell(4, nSubj);
trials_powratio_fullscan = cell(4, nSubj);
trials_powratio_early = cell(4, nSubj);
trials_powratio_late  = cell(4, nSubj);
trials_peaks = cell(4, nSubj);
trials_peaks_early = cell(4, nSubj);
trials_peaks_late  = cell(4, nSubj);
trials_outlier_mask_freq_full = cell(4, nSubj);
trials_outlier_mask_freq_early = cell(4, nSubj);
trials_outlier_mask_freq_late = cell(4, nSubj);
trials_outlier_mask_power_full = cell(4, nSubj);
trials_outlier_mask_power_early = cell(4, nSubj);
trials_outlier_mask_power_late = cell(4, nSubj);
trials_centroid     = cell(4, nSubj);

trials_mean   = nan(4, nSubj);
trials_median = nan(4, nSubj);
trials_mean_early   = nan(4, nSubj);
trials_median_early = nan(4, nSubj);
trials_mean_late    = nan(4, nSubj);
trials_median_late  = nan(4, nSubj);
trials_median_centroid_early = nan(4, nSubj);
trials_median_centroid_late  = nan(4, nSubj);
trials_trialcv_early       = nan(4, nSubj);
trials_trialcv_late        = nan(4, nSubj);
trials_mean_centroid   = nan(4, nSubj);
trials_median_centroid = nan(4, nSubj);

trials_gamma_power = nan(4, nSubj);
trials_gamma_power_early = nan(4, nSubj);
trials_gamma_power_late  = nan(4, nSubj);

all_topos       = cell(1, nSubj);
all_topos_early = cell(1, nSubj);
all_topos_late  = cell(1, nSubj);
all_topo_labels = cell(1, nSubj);
all_combined_spectrum_full = cell(1, nSubj);
all_combined_eigenvalue_full = nan(1, nSubj);
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
all_component_selection_stats_full  = cell(1, nSubj);
all_component_selection_stats_early = cell(1, nSubj);
all_component_selection_stats_late  = cell(1, nSubj);
all_combined_filter_full  = cell(1, nSubj);
all_combined_filter_early = cell(1, nSubj);
all_combined_filter_late  = cell(1, nSubj);
all_haufe_pattern_full = cell(4, nSubj);
all_haufe_patterns_multicomp_full = cell(4, nSubj);
all_reconstruction_patterns_full = cell(1, nSubj);
freq_reconstructed_multicomp_full = cell(4, nSubj);
subject_runtime_seconds = nan(nSubj, 1);

trials_powratio_components_full  = cell(4, nSubj);
trials_powratio_components_early = cell(4, nSubj);
trials_powratio_components_late  = cell(4, nSubj);
all_condition_powspctrm_full = cell(4, nSubj);
all_condition_powspctrm_early = cell(4, nSubj);
all_condition_powspctrm_late = cell(4, nSubj);
all_condition_powspctrm_full_unsmoothed = cell(4, nSubj);
all_condition_powspctrm_early_unsmoothed = cell(4, nSubj);
all_condition_powspctrm_late_unsmoothed = cell(4, nSubj);
freq_powspctrm_full = cell(4, nSubj);
freq_powspctrm_full_unsmoothed = cell(4, nSubj);
all_condition_peak_freq_full = nan(4, nSubj);
all_condition_peak_freq_early = nan(4, nSubj);
all_condition_peak_freq_late = nan(4, nSubj);
all_condition_peak_power_full = nan(4, nSubj);
all_condition_peak_power_early = nan(4, nSubj);
all_condition_peak_power_late = nan(4, nSubj);

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

    fsample = dataEEG_c25.fsample;

    trialIndices = { ...
        find(dataEEG_c25.trialinfo  == 61), ...
        find(dataEEG_c50.trialinfo  == 62), ...
        find(dataEEG_c75.trialinfo  == 63), ...
        find(dataEEG_c100.trialinfo == 64)};
    dataStructs = {dataEEG_c25, dataEEG_c50, dataEEG_c75, dataEEG_c100};

    nChans = length(dataEEG_c25.label);

    % Find channels
    occ_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P10|P9)', 'once')), dataEEG_c25.label);
    occ_idx  = find(occ_mask);
    nOcc     = length(occ_idx);
    front_mask = cellfun(@(l) ~isempty(regexp(l, '^(Fp|AF|F)', 'once')), dataEEG_c25.label);
    front_idx  = find(front_mask);
    post_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P)', 'once')), dataEEG_c25.label);
    post_idx  = find(post_mask);
    temp_mask = cellfun(@(l) ~isempty(regexp(l, '^(T|TP|FT)', 'once')), dataEEG_c25.label);
    temp_idx = find(temp_mask);
    post_w = zeros(nChans, 1);
    post_w(post_idx) = 1;
    if sum(post_w) > 0
        post_w = post_w / sum(post_w);
    else
        post_w = ones(nChans, 1) / nChans;
    end

    %% Build pooled covariance per window
    clc; close all; fprintf('[GED] Subject %s (%d/%d) GED (early, full, late) (%d occ channels)\n', subjects{subj}, subj, nSubj, nOcc);
    rng(random_seed + subj, 'twister');

    stim_windows = {full_window, early_window, late_window};
    win_names   = {'full', 'early', 'late'};
    win_names_cap = {'Full', 'Early', 'Late'};
    lambdas     = [lambda, lambda, lambda];

    covStim_full  = zeros(nChans);
    covStim_early = zeros(nChans);
    covStim_late  = zeros(nChans);
    covBase_full  = zeros(nChans);
    covStim_full_by_cond = cell(1, 4);
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
        cfg_filt.bpfreq     = analysis_freq_range;
        cfg_filt.bpfilttype = 'fir';
        cfg_filt.bpfiltord  = round(3 * fsample / analysis_freq_range(1));
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
        if nTrl > 0
            cfg_cov = [];
            cfg_cov.covariance = 'yes';
            cfg_cov.covariancewindow = 'all';
            cfg_cov.removemean = 'yes';

            tl_base = ft_timelockanalysis(cfg_cov, dat_base);
            tl_stim_full = ft_timelockanalysis(cfg_cov, dat_stim_full);
            tl_stim_early = ft_timelockanalysis(cfg_cov, dat_stim_early);
            tl_stim_late = ft_timelockanalysis(cfg_cov, dat_stim_late);

            covBase_full = covBase_full + double(tl_base.cov) * nTrl;
            covStim_full = covStim_full + double(tl_stim_full.cov) * nTrl;
            covStim_early = covStim_early + double(tl_stim_early.cov) * nTrl;
            covStim_late = covStim_late + double(tl_stim_late.cov) * nTrl;
            covStim_full_by_cond{cond} = double(tl_stim_full.cov);
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

    % Covariances per window
    covStim_per_win = {covStim_full, covStim_early, covStim_late};
    plot_covariance_matrix_diagnostics( ...
        fig_save_dir_component_selection, subjects{subj}, dataEEG_c25.label, ...
        covBase_full, covStim_per_win, win_names_cap, lambdas);

    %% Run GED + component selection per window
    searchFilters_full  = []; searchFilters_early = []; searchFilters_late  = [];
    selected_idx_full  = []; selected_idx_early  = []; selected_idx_late  = [];
    w_combined_full  = []; w_combined_early  = []; w_combined_late  = [];
    evals_sorted_full = []; evals_sorted_early = []; evals_sorted_late = [];
    searchCorrs_full = []; searchCorrs_early = []; searchCorrs_late = [];
    searchTopos_full = []; searchTopos_early = []; searchTopos_late = [];
    searchMeanPrSpectrum_full = []; searchMeanPrSpectrum_early = []; searchMeanPrSpectrum_late = [];
    searchEmgClass_full = {}; searchEmgClass_early = {}; searchEmgClass_late = {};
    eligible_full = []; eligible_early = []; eligible_late = [];
    occdom_full = []; occdom_early = []; occdom_late = [];
    emg_temp_full = []; emg_temp_early = []; emg_temp_late = [];
    emg_hf_slope_full = []; emg_hf_slope_early = []; emg_hf_slope_late = [];
    powspctrm_form_score_full = []; powspctrm_form_score_early = []; powspctrm_form_score_late = [];

    % Simulated signed occipital template (same for all windows)
    template_front_weight = 0.75; % anti-template weight for frontal channels
    template_sigma_occ = 0.12;   % spatial smoothness for occipital template
    template_sigma_front = 0.25; % spatial smoothness for frontal anti-template
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
        % Window-matched GED default: regularize baseline with the same lambda
        % used for the current stimulus window.
        covBase_reg = (1-lam_w)*covBase_full + lam_w*mean(diag(covBase_full))*eye(nChans);

        [W_full, D_full] = eig(covStim_reg, covBase_reg);
        [evals_sorted, sortIdx] = sort(real(diag(D_full)), 'descend');
        W_full = W_full(:, sortIdx);

        nSearch = min(ged_search_n, size(W_full, 2));
        searchFilters = nan(nChans, nSearch);
        searchTopos = nan(nChans, nSearch);
        searchCorrs = nan(nSearch, 1);
        searchOccStrength = nan(nSearch, 1);
        searchFrontStrength = nan(nSearch, 1);
        searchTempStrength = nan(nSearch, 1);
        searchOccdom = nan(nSearch, 1);
        searchEmgTemp = nan(nSearch, 1);
        searchEmgHfSlope = nan(nSearch, 1);
        searchEmgClass = repmat({'unassigned'}, nSearch, 1);
        searchMeanPrSpectrum = nan(nSearch, numel(scan_freqs));

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
                occdom_ci = occ_strength / max(front_strength, eps);
            else
                front_strength = 0;
                occdom_ci = Inf;
            end
            if ~isempty(temp_idx)
                temp_strength = mean(abs(topo_ci(temp_idx)));
            else
                temp_strength = 0;
            end
            emg_temp_ci = temp_strength / max(occ_strength, eps);
            proxy_ci = estimate_component_artifact_proxies( ...
                w_ci, dat_per_cond, stim_windows{w}, baseline_window, fsample, scan_freqs, mtmfft_tapsmofrq_hz);

            searchFilters(:, ci) = w_ci;
            searchTopos(:, ci) = topo_ci;
            searchCorrs(ci) = r_ci;
            searchOccStrength(ci) = occ_strength;
            searchFrontStrength(ci) = front_strength;
            searchTempStrength(ci) = temp_strength;
            searchOccdom(ci) = occdom_ci;
            searchEmgTemp(ci) = emg_temp_ci;
            searchEmgHfSlope(ci) = proxy_ci.hf_slope;
            searchMeanPrSpectrum(ci, :) = proxy_ci.mean_pr_spectrum(:)';
        end

        % Stage-1 gates: SNR, PF, occipital>frontal, not EMG
        eval_raw_vec = evals_sorted(1:nSearch);
        occdom_vec = searchOccdom;
        emg_temp_vec = searchEmgTemp;
        emg_hf_slope_vec = searchEmgHfSlope;
        emg_hf_slope_vec(~isfinite(emg_hf_slope_vec)) = 0;
        [powspctrm_form_score_vec, ~] = compute_powspctrm_form_laplacian_score_from_spectra( ...
            searchMeanPrSpectrum, scan_freqs, analysis_freq_range);
        finite_metrics = isfinite(eval_raw_vec) & isfinite(occdom_vec) & ...
            isfinite(emg_temp_vec) & isfinite(powspctrm_form_score_vec);
        pass_eig_gate = finite_metrics & (eval_raw_vec >= min_eigval);
        pass_peak_gate = finite_metrics & (powspctrm_form_score_vec >= min_powspctrm_form);
        pass_occdom_gate = finite_metrics & (occdom_vec > 1);
        fail_emg_temp = finite_metrics & (emg_temp_vec >= 1);
        fail_emg_hf_slope = finite_metrics & (emg_hf_slope_vec > 0);
        fail_emg = fail_emg_temp | fail_emg_hf_slope;
        for ci = 1:nSearch
            if fail_emg(ci)
                searchEmgClass{ci} = 'EMG';
            elseif ~(occdom_vec(ci) > 1)
                searchEmgClass{ci} = 'frontal';
            elseif searchTempStrength(ci) >= searchOccStrength(ci) && ...
                    searchTempStrength(ci) >= searchFrontStrength(ci)
                searchEmgClass{ci} = 'temporal';
            else
                searchEmgClass{ci} = 'occipital';
            end
        end
        eligible = pass_eig_gate & pass_peak_gate & pass_occdom_gate & ~fail_emg;
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
            combined_idx = combined_idx(1:min(max_components_to_combine, numel(combined_idx)));

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
            bestOcc = NaN;
            bestFront = NaN;
            bestRatio = NaN;
            bestLeak = NaN;
            topo_temp = nan(nChans, 1);
        else
            bestIdx = selected_idx(1);
            bestScore = searchScores(bestIdx);
            bestCorr = searchCorrs(bestIdx);
            bestOcc = searchOccStrength(bestIdx);
            bestFront = searchFrontStrength(bestIdx);
            bestRatio = occdom_vec(bestIdx);
            bestLeak = 1 / max(occdom_vec(bestIdx), eps);

            topComp = searchFilters(:, bestIdx);
            if numel(selected_idx) > 1
                topo_temp = searchTopos(:, selected_idx) * selected_weights(:);
            else
                topo_temp = covStim_reg * topComp;
            end
        end
        [~, topDispOrder] = sort(evals_sorted(1:nSearch), 'descend');
        nStore = min(5, nSearch);
        storeCompIdx = topDispOrder(1:nStore);
        storeCorrsTop = searchCorrs(storeCompIdx);
        storeEvalsTop = evals_sorted(storeCompIdx);
        storeTopos = searchTopos(:, storeCompIdx);
        storeCorrs = nan(5, 1);
        storeEvals = nan(5, 1);
        storeCorrs(1:nStore) = storeCorrsTop(:);
        storeEvals(1:nStore) = storeEvalsTop(:);

        % Store per-window filters for trial-level scanning
        if w == 1
            searchFilters_full = searchFilters;
            selected_idx_full = selected_idx;
            w_combined_full = selected_weights(:)';
            evals_sorted_full = evals_sorted;
            searchCorrs_full = searchCorrs;
            searchTopos_full = searchTopos;
            searchMeanPrSpectrum_full = searchMeanPrSpectrum;
            searchEmgClass_full = searchEmgClass;
            eligible_full = eligible;
            occdom_full = occdom_vec;
            emg_temp_full = emg_temp_vec;
            emg_hf_slope_full = emg_hf_slope_vec;
            powspctrm_form_score_full = powspctrm_form_score_vec;
        elseif w == 2
            searchFilters_early = searchFilters;
            selected_idx_early = selected_idx;
            w_combined_early = selected_weights(:)';
            evals_sorted_early = evals_sorted;
            searchCorrs_early = searchCorrs;
            searchTopos_early = searchTopos;
            searchMeanPrSpectrum_early = searchMeanPrSpectrum;
            searchEmgClass_early = searchEmgClass;
            eligible_early = eligible;
            occdom_early = occdom_vec;
            emg_temp_early = emg_temp_vec;
            emg_hf_slope_early = emg_hf_slope_vec;
            powspctrm_form_score_early = powspctrm_form_score_vec;
        else
            searchFilters_late = searchFilters;
            selected_idx_late = selected_idx;
            w_combined_late = selected_weights(:)';
            evals_sorted_late = evals_sorted;
            searchCorrs_late = searchCorrs;
            searchTopos_late = searchTopos;
            searchMeanPrSpectrum_late = searchMeanPrSpectrum;
            searchEmgClass_late = searchEmgClass;
            eligible_late = eligible;
            occdom_late = occdom_vec;
            emg_temp_late = emg_temp_vec;
            emg_hf_slope_late = emg_hf_slope_vec;
            powspctrm_form_score_late = powspctrm_form_score_vec;
        end

        if w == 1
            all_topos{subj}       = topo_temp;
            all_topo_labels{subj} = dataEEG_c25.label;
            if isempty(selected_idx) || ~isfinite(bestIdx)
                all_eigenvalues(subj) = NaN;
                all_selected_comp_idx(subj)  = NaN;
                all_selected_comp_eval(subj) = NaN;
            else
                all_eigenvalues(subj) = evals_sorted(bestIdx);
                all_selected_comp_idx(subj)  = bestIdx;
                all_selected_comp_eval(subj) = evals_sorted(bestIdx);
            end
            all_selected_comp_corr(subj) = bestCorr;
            all_top5_corrs(:, subj) = storeCorrs;
            all_top5_evals(:, subj) = storeEvals;
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
            'selection_mode', 'fixed_preregistered_weighted', ...
            'selected_idx', selected_idx, ...
            'n_selected_ged_components', numel(selected_idx), ...
            'selected_weights', selected_weights, ...
            'best_idx', bestIdx, ...
            'best_score', bestScore, ...
            'best_corr', bestCorr, ...
            'best_occdom', bestRatio, ...
            'best_front', bestFront, ...
            'best_occ', bestOcc, ...
            'best_front_leak', bestLeak, ...
            'occdom', occdom_vec, ...
            'emg_temp', emg_temp_vec, ...
            'emg_hf_slope', emg_hf_slope_vec, ...
            'emg_class', {searchEmgClass}, ...
            'eligible', eligible, ...
            'no_threshold_match', no_threshold_match);
        if w == 1
            all_component_selection_stats_full{subj} = comp_sel_struct;
        elseif w == 2
            all_component_selection_stats_early{subj} = comp_sel_struct;
        else
            all_component_selection_stats_late{subj} = comp_sel_struct;
        end
    end

    if isempty(selected_idx_full)
        W_combined_full = [];
        topo_temp_full = nan(nChans, 1);
    else
        W_combined_full = searchFilters_full(:, selected_idx_full);
        topo_temp_full = searchTopos_full(:, selected_idx_full) * w_combined_full(:);
    end
    if isempty(selected_idx_early)
        W_combined_early = [];
        topo_temp_early = nan(nChans, 1);
    else
        W_combined_early = searchFilters_early(:, selected_idx_early);
        topo_temp_early = searchTopos_early(:, selected_idx_early) * w_combined_early(:);
    end
    if isempty(selected_idx_late)
        W_combined_late = [];
        topo_temp_late = nan(nChans, 1);
    else
        W_combined_late = searchFilters_late(:, selected_idx_late);
        topo_temp_late = searchTopos_late(:, selected_idx_late) * w_combined_late(:);
    end
    all_topos{subj} = topo_temp_full;
    all_topos_early{subj} = topo_temp_early;
    all_topos_late{subj} = topo_temp_late;
    all_selected_comp_indices_multi{subj} = selected_idx_full;
    all_selected_comp_weights{subj} = w_combined_full(:)';
    [sel_idx_spec, sel_w_spec] = sanitize_selected_components( ...
        selected_idx_full, w_combined_full, size(searchMeanPrSpectrum_full, 1));
    if isempty(sel_idx_spec) || isempty(searchMeanPrSpectrum_full)
        all_combined_spectrum_full{subj} = [];
        all_combined_eigenvalue_full(subj) = NaN;
    else
        all_combined_spectrum_full{subj} = sel_w_spec(:)' * searchMeanPrSpectrum_full(sel_idx_spec, :);
        evals_sel = evals_sorted_full(sel_idx_spec);
        evals_sel = evals_sel(:);
        evals_sel(~isfinite(evals_sel)) = NaN;
        all_combined_eigenvalue_full(subj) = sum(sel_w_spec(:) .* evals_sel);
    end
    if isempty(selected_idx_full)
        all_selected_comp_idx(subj) = NaN;
        all_selected_comp_corr(subj) = NaN;
        all_selected_comp_eval(subj) = NaN;
        all_eigenvalues(subj) = NaN;
    else
        all_selected_comp_idx(subj) = selected_idx_full(1);
        all_selected_comp_corr(subj) = searchCorrs_full(selected_idx_full(1));
        all_selected_comp_eval(subj) = evals_sorted_full(selected_idx_full(1));
        all_eigenvalues(subj) = evals_sorted_full(selected_idx_full(1));
    end
    comp_stats_full = all_component_selection_stats_full{subj};
    comp_stats_early = all_component_selection_stats_early{subj};
    comp_stats_late = all_component_selection_stats_late{subj};

    comp_stats_full.selected_idx = selected_idx_full;
    comp_stats_full.selected_weights = w_combined_full;

    comp_stats_early.selected_idx = selected_idx_early;
    comp_stats_early.selected_weights = w_combined_early;

    comp_stats_late.selected_idx = selected_idx_late;
    comp_stats_late.selected_weights = w_combined_late;

    all_component_selection_stats_full{subj} = comp_stats_full;
    all_component_selection_stats_early{subj} = comp_stats_early;
    all_component_selection_stats_late{subj} = comp_stats_late;

    cfg_topo = [];
    cfg_topo.layout    = headmodel.layANThead;
    cfg_topo.comment   = 'no';
    cfg_topo.marker    = 'off';
    cfg_topo.style     = 'straight';
    cfg_topo.gridscale = 300;
    cfg_topo.zlim      = 'maxabs';
    cfg_topo.colormap  = '*RdBu';
    cfg_topo.figure    = 'gcf';
    plot_emg_exclusion_diagnostics( ...
        fig_save_dir_component_selection, subjects{subj}, 'full', scan_freqs, searchTopos_full, ...
        searchMeanPrSpectrum_full, evals_sorted_full(1:numel(eligible_full)), ...
        searchEmgClass_full, ...
        eligible_full, ...
        occdom_full, emg_temp_full, emg_hf_slope_full, ...
        cfg_topo, all_topo_labels{subj}, powspctrm_form_score_full, ...
        selected_idx_full);
    plot_emg_exclusion_diagnostics( ...
        fig_save_dir_component_selection, subjects{subj}, 'early', scan_freqs, searchTopos_early, ...
        searchMeanPrSpectrum_early, evals_sorted_early(1:numel(eligible_early)), ...
        searchEmgClass_early, ...
        eligible_early, ...
        occdom_early, emg_temp_early, emg_hf_slope_early, ...
        cfg_topo, all_topo_labels{subj}, powspctrm_form_score_early, ...
        selected_idx_early);
    plot_emg_exclusion_diagnostics( ...
        fig_save_dir_component_selection, subjects{subj}, 'late', scan_freqs, searchTopos_late, ...
        searchMeanPrSpectrum_late, evals_sorted_late(1:numel(eligible_late)), ...
        searchEmgClass_late, ...
        eligible_late, ...
        occdom_late, emg_temp_late, emg_hf_slope_late, ...
        cfg_topo, all_topo_labels{subj}, powspctrm_form_score_late, ...
        selected_idx_late);
    plot_combined_topo_spectra_windows( ...
        fig_save_dir_component_selection, subjects{subj}, scan_freqs, cfg_topo, all_topo_labels{subj}, ...
        searchTopos_full, searchMeanPrSpectrum_full, selected_idx_full, w_combined_full, ...
        searchTopos_early, searchMeanPrSpectrum_early, selected_idx_early, w_combined_early, ...
        searchTopos_late, searchMeanPrSpectrum_late, selected_idx_late, w_combined_late, ...
        analysis_freq_range);

    adequate_full = false;
    adequate_early = false;
    adequate_late = false;
    if ~isempty(all_component_selection_stats_full{subj}) && isfield(all_component_selection_stats_full{subj}, 'selected_idx')
        adequate_full = ~isempty(all_component_selection_stats_full{subj}.selected_idx);
    end
    if ~isempty(all_component_selection_stats_early{subj}) && isfield(all_component_selection_stats_early{subj}, 'selected_idx')
        adequate_early = ~isempty(all_component_selection_stats_early{subj}.selected_idx);
    end
    if ~isempty(all_component_selection_stats_late{subj}) && isfield(all_component_selection_stats_late{subj}, 'selected_idx')
        adequate_late = ~isempty(all_component_selection_stats_late{subj}.selected_idx);
    end
    if ~adequate_full
        all_selected_comp_idx(subj) = NaN;
        all_selected_comp_corr(subj) = NaN;
        all_selected_comp_eval(subj) = NaN;
        all_eigenvalues(subj) = NaN;
        all_selected_comp_indices_multi{subj} = NaN;
        all_selected_comp_weights{subj} = NaN;
    end

    %% Per-condition trial-level spectral scanning
    % Use window-specific filters for each dB-spectrum output
    for wi = 1:3
        if wi == 1
            W_comb = W_combined_full;
            w_comb = w_combined_full;
            sel_idx = selected_idx_full;
            window_adequate = adequate_full;
        elseif wi == 2
            W_comb = W_combined_early;
            w_comb = w_combined_early;
            sel_idx = selected_idx_early;
            window_adequate = adequate_early;
        else
            W_comb = W_combined_late;
            w_comb = w_combined_late;
            sel_idx = selected_idx_late;
            window_adequate = adequate_late;
        end
        if ~window_adequate
            W_comb = zeros(nChans, 0);
            w_comb = [];
            sel_idx = [];
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
        W_comb = normalize_filters_to_noise_metric(W_comb, covBase_full);
        if wi == 1
            W_combined_full_norm = W_comb;
            w_combined_full_norm = w_comb;
            selected_idx_full_norm = sel_idx;
        elseif wi == 2
            W_combined_early_norm = W_comb;
            w_combined_early_norm = w_comb;
            selected_idx_early_norm = sel_idx;
        else
            W_combined_late_norm = W_comb;
            w_combined_late_norm = w_comb;
            selected_idx_late_norm = sel_idx;
        end
    end
    % Per-window filters struct for trial-level scanning
    filters = struct('full', struct(), 'early', struct(), 'late', struct());
    filters.full.searchFilters = normalize_filters_to_noise_metric(searchFilters_full, covBase_full);
    filters.full.W_combined = W_combined_full_norm;
    filters.full.selected_idx = selected_idx_full_norm;
    filters.full.w_combined = w_combined_full_norm;
    filters.early.searchFilters = normalize_filters_to_noise_metric(searchFilters_early, covBase_full);
    filters.early.W_combined = W_combined_early_norm;
    filters.early.selected_idx = selected_idx_early_norm;
    filters.early.w_combined = w_combined_early_norm;
    filters.late.searchFilters = normalize_filters_to_noise_metric(searchFilters_late, covBase_full);
    filters.late.W_combined = W_combined_late_norm;
    filters.late.selected_idx = selected_idx_late_norm;
    filters.late.w_combined = w_combined_late_norm;

    % Persist the exact signed, noise-normalized channel-space filters used
    % for trial projection. Downstream analyses must reuse these vectors
    % rather than reconstructing the eigendecomposition.
    all_combined_filter_full{subj} = build_combined_filter_vector( ...
        filters.full.W_combined, filters.full.w_combined);
    all_combined_filter_early{subj} = build_combined_filter_vector( ...
        filters.early.W_combined, filters.early.w_combined);
    all_combined_filter_late{subj} = build_combined_filter_vector( ...
        filters.late.W_combined, filters.late.w_combined);

    %% Condition-specific Haufe patterns and multicomponent reconstruction
    if adequate_full && ~isempty(filters.full.W_combined)
        W_selected_full = filters.full.W_combined;
        component_cov_pool = W_selected_full' * covStim_full * W_selected_full;
        reconstruction_patterns = covStim_full * W_selected_full * pinv(component_cov_pool);
        all_reconstruction_patterns_full{subj} = reconstruction_patterns;

        combined_filter_full = all_combined_filter_full{subj};
        for cond = 1:4
            cov_cond = covStim_full_by_cond{cond};
            dat_cond = dat_per_cond{cond};
            if isempty(cov_cond) || isempty(dat_cond)
                continue;
            end

            component_cov_cond = W_selected_full' * cov_cond * W_selected_full;
            all_haufe_patterns_multicomp_full{cond, subj} = ...
                cov_cond * W_selected_full * pinv(component_cov_cond);

            combined_variance_cond = combined_filter_full' * cov_cond * combined_filter_full;
            if isfinite(combined_variance_cond) && combined_variance_cond > eps
                all_haufe_pattern_full{cond, subj} = ...
                    (cov_cond * combined_filter_full) / combined_variance_cond;
            end

            % Reconstruct the selected GED subspace at the sensors using one
            % pooled activation matrix so condition maps remain comparable.
            dat_reconstructed = dat_cond;
            for trl = 1:numel(dat_cond.trial)
                x = double(dat_cond.trial{trl});
                component_data = W_selected_full' * x;
                dat_reconstructed.trial{trl} = reconstruction_patterns * component_data;
            end

            cfg_select = [];
            cfg_select.latency = baseline_window;
            dat_reconstructed_base = ft_selectdata(cfg_select, dat_reconstructed);
            cfg_select.latency = full_window;
            dat_reconstructed_stim = ft_selectdata(cfg_select, dat_reconstructed);

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
            freq_reconstructed_multicomp_full{cond, subj} = freq_ratio;
        end
    end

    subj_powratio_fullscan = cell(1, 4);
    subj_powratio_early = cell(1, 4);
    subj_powratio_late = cell(1, 4);
    subj_peaks_full = cell(1, 4);
    subj_peaks_early = cell(1, 4);
    subj_peaks_late = cell(1, 4);
    subj_centroid_full = cell(1, 4);
    subj_centroid_early = cell(1, 4);
    subj_centroid_late = cell(1, 4);
    subj_condition_avg_full = cell(1, 4);
    subj_condition_avg_early = cell(1, 4);
    subj_condition_avg_late = cell(1, 4);
    subj_condition_peak_full = nan(1, 4);
    subj_condition_peak_early = nan(1, 4);
    subj_condition_peak_late = nan(1, 4);
    for cond = 1:4

        dat = dat_per_cond{cond};
        if isempty(dat), continue; end

        nTrl = length(dat.trial);
        powratio_methods_full = nan(1, nTrl, nFreqs);
        powratio_methods_early = nan(1, nTrl, nFreqs);
        powratio_methods_late = nan(1, nTrl, nFreqs);
        nSearch_full = size(filters.full.searchFilters, 2);
        nSearch_early = size(filters.early.searchFilters, 2);
        nSearch_late = size(filters.late.searchFilters, 2);
        powratio_components       = nan(nSearch_full, nTrl, nFreqs);
        powratio_components_early = nan(nSearch_early, nTrl, nFreqs);
        powratio_components_late  = nan(nSearch_late, nTrl, nFreqs);
        unstable_freq_counts_full = zeros(nTrl, 1);
        unstable_freq_counts_early = zeros(nTrl, 1);
        unstable_freq_counts_late = zeros(nTrl, 1);
        valid_freq_counts_full = zeros(nTrl, 1);
        valid_freq_counts_early = zeros(nTrl, 1);
        valid_freq_counts_late = zeros(nTrl, 1);

        % Baseline quality gate computed once per trial (not per frequency).
        baseline_power_raw = nan(nTrl, 1);
        baseline_power_comb_full = nan(nTrl, 1);
        baseline_power_comb_early = nan(nTrl, 1);
        baseline_power_comb_late = nan(nTrl, 1);
        trial_cache = cell(nTrl, 1);
        for trl = 1:nTrl
            x = double(dat.trial{trl});
            t = dat.time{trl};
            idx_base = t >= baseline_window(1) & t <= baseline_window(2);
            idx_full = t >= full_window(1) & t <= full_window(2);
            idx_early = t >= early_window(1) & t <= early_window(2);
            idx_late = t >= late_window(1) & t <= late_window(2);
            x_base = x(:, idx_base);
            x_full = x(:, idx_full);
            x_early = x(:, idx_early);
            x_late = x(:, idx_late);
            trial_cache{trl} = struct('x_base', x_base, 'x_full', x_full, 'x_early', x_early, 'x_late', x_late);
            if ~isempty(x_base)
                pow_base_chan = mean(x_base.^2, 2);
                baseline_power_raw(trl) = sum(post_w(:) .* pow_base_chan(:));
                if adequate_full && ~isempty(filters.full.W_combined)
                    x_base_full = filters.full.W_combined' * x_base;
                    baseline_power_comb_full(trl) = mean(x_base_full(:).^2);
                end
                if adequate_early && ~isempty(filters.early.W_combined)
                    x_base_early = filters.early.W_combined' * x_base;
                    baseline_power_comb_early(trl) = mean(x_base_early(:).^2);
                end
                if adequate_late && ~isempty(filters.late.W_combined)
                    x_base_late = filters.late.W_combined' * x_base;
                    baseline_power_comb_late(trl) = mean(x_base_late(:).^2);
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
        bad_base_full = bad_base_raw;
        bad_base_early = bad_base_raw;
        bad_base_late = bad_base_raw;
        if adequate_full
            bad_base_full = bad_base_full | flag_unreliable_baseline_trials( ...
                baseline_power_comb_full, baseline_outlier_mad_mult);
        end
        if adequate_early
            bad_base_early = bad_base_early | flag_unreliable_baseline_trials( ...
                baseline_power_comb_early, baseline_outlier_mad_mult);
        end
        if adequate_late
            bad_base_late = bad_base_late | flag_unreliable_baseline_trials( ...
                baseline_power_comb_late, baseline_outlier_mad_mult);
        end
        [base_floor_full, ~] = compute_baseline_floor_stats(baseline_power_comb_full, ratio_floor_prctile, ratio_floor_frac);
        [base_floor_early, ~] = compute_baseline_floor_stats(baseline_power_comb_early, ratio_floor_prctile, ratio_floor_frac);
        [base_floor_late, ~] = compute_baseline_floor_stats(baseline_power_comb_late, ratio_floor_prctile, ratio_floor_frac);

        has_base = false(nTrl, 1);
        has_full = false(nTrl, 1);
        has_early = false(nTrl, 1);
        has_late = false(nTrl, 1);
        for trl = 1:nTrl
            tc = trial_cache{trl};
            has_base(trl) = ~isempty(tc.x_base);
            has_full(trl) = ~isempty(tc.x_full);
            has_early(trl) = ~isempty(tc.x_early);
            has_late(trl) = ~isempty(tc.x_late);
        end

        if adequate_full
            trial_mask_full = has_base & has_full & ~bad_base_full;
            [ratio_cube_full, near_floor_count_full] = compute_scan_ratio_for_window_batch( ...
                trial_cache, filters.full.W_combined, 'x_full', trial_mask_full, ...
                fsample, scan_freqs, mtmfft_tapsmofrq_hz, base_floor_full, instability_near_floor_mult);
            powratio_components = ratio_cube_full;
            filter_vec_full = build_combined_filter_vector(filters.full.W_combined, filters.full.w_combined);
            [ratio_trials_full_combined, near_floor_count_full_combined, valid_freq_counts_full_combined] = ...
                compute_scan_ratio_for_combined_filter_batch( ...
                trial_cache, filter_vec_full, 'x_full', trial_mask_full, ...
                fsample, scan_freqs, mtmfft_tapsmofrq_hz, base_floor_full, instability_near_floor_mult);
            for trl = 1:nTrl
                ratio_mat_full = squeeze(powratio_components(:, trl, :));
                if isvector(ratio_mat_full)
                    ratio_mat_full = reshape(ratio_mat_full, size(powratio_components, 1), []);
                end
                if ~isempty(ratio_trials_full_combined)
                    powratio_methods_full(1, trl, :) = ratio_trials_full_combined(trl, :);
                end
                if trl <= numel(valid_freq_counts_full_combined)
                    valid_freq_counts_full(trl) = valid_freq_counts_full_combined(trl);
                end
                if trl <= numel(near_floor_count_full_combined) && valid_freq_counts_full(trl) > 0
                    unstable_freq_counts_full(trl) = near_floor_count_full_combined(trl);
                elseif any(any(isfinite(ratio_mat_full), 1))
                    unstable_freq_counts_full(trl) = near_floor_count_full(trl);
                end
            end
        end

        if adequate_early
            trial_mask_early = has_base & has_early & ~bad_base_early;
            [ratio_cube_early, near_floor_count_early] = compute_scan_ratio_for_window_batch( ...
                trial_cache, filters.early.W_combined, 'x_early', trial_mask_early, ...
                fsample, scan_freqs, mtmfft_tapsmofrq_hz, base_floor_early, instability_near_floor_mult);
            powratio_components_early = ratio_cube_early;
            filter_vec_early = build_combined_filter_vector(filters.early.W_combined, filters.early.w_combined);
            [ratio_trials_early_combined, near_floor_count_early_combined, valid_freq_counts_early_combined] = ...
                compute_scan_ratio_for_combined_filter_batch( ...
                trial_cache, filter_vec_early, 'x_early', trial_mask_early, ...
                fsample, scan_freqs, mtmfft_tapsmofrq_hz, base_floor_early, instability_near_floor_mult);
            for trl = 1:nTrl
                ratio_mat_early = squeeze(powratio_components_early(:, trl, :));
                if isvector(ratio_mat_early)
                    ratio_mat_early = reshape(ratio_mat_early, size(powratio_components_early, 1), []);
                end
                if ~isempty(ratio_trials_early_combined)
                    powratio_methods_early(1, trl, :) = ratio_trials_early_combined(trl, :);
                end
                if trl <= numel(valid_freq_counts_early_combined)
                    valid_freq_counts_early(trl) = valid_freq_counts_early_combined(trl);
                end
                if trl <= numel(near_floor_count_early_combined) && valid_freq_counts_early(trl) > 0
                    unstable_freq_counts_early(trl) = near_floor_count_early_combined(trl);
                elseif any(any(isfinite(ratio_mat_early), 1))
                    unstable_freq_counts_early(trl) = near_floor_count_early(trl);
                end
            end
        end

        if adequate_late
            trial_mask_late = has_base & has_late & ~bad_base_late;
            [ratio_cube_late, near_floor_count_late] = compute_scan_ratio_for_window_batch( ...
                trial_cache, filters.late.W_combined, 'x_late', trial_mask_late, ...
                fsample, scan_freqs, mtmfft_tapsmofrq_hz, base_floor_late, instability_near_floor_mult);
            powratio_components_late = ratio_cube_late;
            filter_vec_late = build_combined_filter_vector(filters.late.W_combined, filters.late.w_combined);
            [ratio_trials_late_combined, near_floor_count_late_combined, valid_freq_counts_late_combined] = ...
                compute_scan_ratio_for_combined_filter_batch( ...
                trial_cache, filter_vec_late, 'x_late', trial_mask_late, ...
                fsample, scan_freqs, mtmfft_tapsmofrq_hz, base_floor_late, instability_near_floor_mult);
            for trl = 1:nTrl
                ratio_mat_late = squeeze(powratio_components_late(:, trl, :));
                if isvector(ratio_mat_late)
                    ratio_mat_late = reshape(ratio_mat_late, size(powratio_components_late, 1), []);
                end
                if ~isempty(ratio_trials_late_combined)
                    powratio_methods_late(1, trl, :) = ratio_trials_late_combined(trl, :);
                end
                if trl <= numel(valid_freq_counts_late_combined)
                    valid_freq_counts_late(trl) = valid_freq_counts_late_combined(trl);
                end
                if trl <= numel(near_floor_count_late_combined) && valid_freq_counts_late(trl) > 0
                    unstable_freq_counts_late(trl) = near_floor_count_late_combined(trl);
                elseif any(any(isfinite(ratio_mat_late), 1))
                    unstable_freq_counts_late(trl) = near_floor_count_late(trl);
                end
            end
        end
        unstable_trial_frac_full = unstable_freq_counts_full ./ max(valid_freq_counts_full, 1);
        unstable_trial_frac_early = unstable_freq_counts_early ./ max(valid_freq_counts_early, 1);
        unstable_trial_frac_late = unstable_freq_counts_late ./ max(valid_freq_counts_late, 1);
        trial_unstable_full = unstable_trial_frac_full >= instability_trial_freq_frac_thr;
        trial_unstable_early = unstable_trial_frac_early >= instability_trial_freq_frac_thr;
        trial_unstable_late = unstable_trial_frac_late >= instability_trial_freq_frac_thr;
        if any(trial_unstable_full)
            powratio_methods_full(:, trial_unstable_full, :) = NaN;
            powratio_components(:, trial_unstable_full, :) = NaN;
        end
        if any(trial_unstable_early)
            powratio_methods_early(:, trial_unstable_early, :) = NaN;
            powratio_components_early(:, trial_unstable_early, :) = NaN;
        end
        if any(trial_unstable_late)
            powratio_methods_late(:, trial_unstable_late, :) = NaN;
            powratio_components_late(:, trial_unstable_late, :) = NaN;
        end
        powratio_methods_full_analysis = powratio_methods_full;
        powratio_methods_early_analysis = powratio_methods_early;
        powratio_methods_late_analysis = powratio_methods_late;
        powratio_components_analysis = powratio_components;
        powratio_components_early_analysis = powratio_components_early;
        powratio_components_late_analysis = powratio_components_late;

        trials_powratio_components_full{cond, subj}  = powratio_components_analysis;
        trials_powratio_components_early{cond, subj} = powratio_components_early_analysis;
        trials_powratio_components_late{cond, subj}  = powratio_components_late_analysis;

        % Keep full-window outputs based on weighted combined GED branch.
        powratio_trials_fullscan = squeeze(powratio_methods_full(1, :, :));
        powratio_trials_early_fullscan = squeeze(powratio_methods_early(1, :, :));
        powratio_trials_late_fullscan = squeeze(powratio_methods_late(1, :, :));
        powratio_trials_full = squeeze(powratio_methods_full_analysis(1, :, :));
        powratio_trials_early = squeeze(powratio_methods_early_analysis(1, :, :));
        powratio_trials_late = squeeze(powratio_methods_late_analysis(1, :, :));
        trials_powratio_fullscan{cond, subj} = powratio_trials_fullscan;
        trials_powratio{cond, subj} = powratio_trials_full;
        trials_powratio_early{cond, subj} = powratio_trials_early;
        trials_powratio_late{cond, subj} = powratio_trials_late;
        subj_powratio_fullscan{cond} = powratio_trials_fullscan;
        subj_powratio_early{cond} = powratio_trials_early;
        subj_powratio_late{cond} = powratio_trials_late;
        %% Per-trial peak detection
        trial_metric_outlier_iqr_mult = 1.5; % outlier threshold in IQR units around Q1/Q3
        [trl_peaks, trial_peak_power_full, trl_centroid] = ...
            compute_trial_peak_metrics_from_powratio_fullscan( ...
            powratio_trials_fullscan, scan_freqs, true(size(scan_freqs)), ...
            powratio_trial_freq_smooth_bins, peak_power_halfwidth_hz);

        trials_peaks{cond, subj} = trl_peaks;
        trials_centroid{cond, subj}     = trl_centroid;
        subj_peaks_full{cond} = trl_peaks;
        subj_centroid_full{cond} = trl_centroid;

        % Time-split peak summaries.
        [trl_peaks_early, trial_peak_power_early, trl_centroid_early] = ...
            compute_trial_peak_metrics_from_powratio_fullscan( ...
            powratio_trials_early_fullscan, scan_freqs, true(size(scan_freqs)), ...
            powratio_trial_freq_smooth_bins, peak_power_halfwidth_hz);
        [trl_peaks_late, trial_peak_power_late, trl_centroid_late] = ...
            compute_trial_peak_metrics_from_powratio_fullscan( ...
            powratio_trials_late_fullscan, scan_freqs, true(size(scan_freqs)), ...
            powratio_trial_freq_smooth_bins, peak_power_halfwidth_hz);
        % Trial-level metric outlier rejection (subject-condition specific).
        [outlier_mask_freq_full, ~] = detect_trial_metric_outliers_iqr( ...
            trl_peaks, trial_metric_outlier_iqr_mult);
        [outlier_mask_freq_early, ~] = detect_trial_metric_outliers_iqr( ...
            trl_peaks_early, trial_metric_outlier_iqr_mult);
        [outlier_mask_freq_late, ~] = detect_trial_metric_outliers_iqr( ...
            trl_peaks_late, trial_metric_outlier_iqr_mult);
        [outlier_mask_power_full, ~] = detect_trial_metric_outliers_iqr( ...
            trial_peak_power_full, trial_metric_outlier_iqr_mult);
        [outlier_mask_power_early, ~] = detect_trial_metric_outliers_iqr( ...
            trial_peak_power_early, trial_metric_outlier_iqr_mult);
        [outlier_mask_power_late, ~] = detect_trial_metric_outliers_iqr( ...
            trial_peak_power_late, trial_metric_outlier_iqr_mult);
        trl_peaks(outlier_mask_freq_full) = NaN;
        trl_peaks_early(outlier_mask_freq_early) = NaN;
        trl_peaks_late(outlier_mask_freq_late) = NaN;
        trial_peak_power_full(outlier_mask_power_full) = NaN;
        trial_peak_power_early(outlier_mask_power_early) = NaN;
        trial_peak_power_late(outlier_mask_power_late) = NaN;
        trials_peaks{cond, subj} = trl_peaks;
        trials_outlier_mask_freq_full{cond, subj} = outlier_mask_freq_full;
        trials_outlier_mask_freq_early{cond, subj} = outlier_mask_freq_early;
        trials_outlier_mask_freq_late{cond, subj} = outlier_mask_freq_late;
        trials_outlier_mask_power_full{cond, subj} = outlier_mask_power_full;
        trials_outlier_mask_power_early{cond, subj} = outlier_mask_power_early;
        trials_outlier_mask_power_late{cond, subj} = outlier_mask_power_late;
        trials_peaks_early{cond, subj} = trl_peaks_early;
        trials_peaks_late{cond, subj} = trl_peaks_late;
        subj_peaks_full{cond} = trl_peaks;
        subj_peaks_early{cond} = trl_peaks_early;
        subj_peaks_late{cond} = trl_peaks_late;
        subj_centroid_early{cond} = trl_centroid_early;
        subj_centroid_late{cond} = trl_centroid_late;
        outlier_rows_full = outlier_mask_freq_full | outlier_mask_power_full;
        outlier_rows_early = outlier_mask_freq_early | outlier_mask_power_early;
        outlier_rows_late = outlier_mask_freq_late | outlier_mask_power_late;
        powratio_trials_full_avg = powratio_trials_full;
        powratio_trials_early_avg = powratio_trials_early;
        powratio_trials_late_avg = powratio_trials_late;
        if ~isempty(powratio_trials_full_avg)
            powratio_trials_full_avg(outlier_rows_full, :) = NaN;
        end
        if ~isempty(powratio_trials_early_avg)
            powratio_trials_early_avg(outlier_rows_early, :) = NaN;
        end
        if ~isempty(powratio_trials_late_avg)
            powratio_trials_late_avg(outlier_rows_late, :) = NaN;
        end
        valid_s = ~isnan(trl_peaks);
        valid_s_early = ~isnan(trl_peaks_early);
        valid_s_late = ~isnan(trl_peaks_late);

        % Condition-level spectra and peak metrics from trial-averaged spectra
        % (FieldTrip averaging over the trial dimension).
        cond_avg_full = compute_condition_average_powratio_ft(powratio_trials_full_avg, scan_freqs);
        cond_avg_early = compute_condition_average_powratio_ft(powratio_trials_early_avg, scan_freqs);
        cond_avg_late = compute_condition_average_powratio_ft(powratio_trials_late_avg, scan_freqs);
        all_condition_powspctrm_full_unsmoothed{cond, subj} = cond_avg_full;
        all_condition_powspctrm_early_unsmoothed{cond, subj} = cond_avg_early;
        all_condition_powspctrm_late_unsmoothed{cond, subj} = cond_avg_late;
        cond_avg_full = movmean(cond_avg_full, max(1, round(powratio_condition_freq_smooth_bins)), 'omitnan');
        cond_avg_early = movmean(cond_avg_early, max(1, round(powratio_condition_freq_smooth_bins)), 'omitnan');
        cond_avg_late = movmean(cond_avg_late, max(1, round(powratio_condition_freq_smooth_bins)), 'omitnan');
        all_condition_powspctrm_full{cond, subj} = cond_avg_full;
        all_condition_powspctrm_early{cond, subj} = cond_avg_early;
        all_condition_powspctrm_late{cond, subj} = cond_avg_late;
        freq_powspctrm_full_unsmoothed{cond, subj} = ged_powcurve_to_freq_ft( ...
            all_condition_powspctrm_full_unsmoothed{cond, subj}, scan_freqs, dataEEG_c25);
        freq_powspctrm_full{cond, subj} = ged_powcurve_to_freq_ft( ...
            all_condition_powspctrm_full{cond, subj}, scan_freqs, dataEEG_c25);
        subj_condition_avg_full{cond} = cond_avg_full;
        subj_condition_avg_early{cond} = cond_avg_early;
        subj_condition_avg_late{cond} = cond_avg_late;

        [peak_full_hz, peak_full_power] = pick_tallest_peak(cond_avg_full, scan_freqs, 1, peak_power_halfwidth_hz);
        [peak_early_hz, peak_early_power] = pick_tallest_peak(cond_avg_early, scan_freqs, 1, peak_power_halfwidth_hz);
        [peak_late_hz, peak_late_power] = pick_tallest_peak(cond_avg_late, scan_freqs, 1, peak_power_halfwidth_hz);
        all_condition_peak_freq_full(cond, subj) = peak_full_hz;
        all_condition_peak_freq_early(cond, subj) = peak_early_hz;
        all_condition_peak_freq_late(cond, subj) = peak_late_hz;
        all_condition_peak_power_full(cond, subj) = peak_full_power;
        all_condition_peak_power_early(cond, subj) = peak_early_power;
        all_condition_peak_power_late(cond, subj) = peak_late_power;
        subj_condition_peak_full(cond) = peak_full_hz;
        subj_condition_peak_early(cond) = peak_early_hz;
        subj_condition_peak_late(cond) = peak_late_hz;

        trials_mean(cond, subj) = robust_trial_mean(trl_peaks(valid_s));
        trials_median(cond, subj) = median(trl_peaks(valid_s));
        valid_c = isfinite(trl_centroid);
        trials_mean_centroid(cond, subj) = robust_trial_mean(trl_centroid(valid_c));
        trials_median_centroid(cond, subj) = median(trl_centroid(valid_c));

        trials_mean_early(cond, subj) = robust_trial_mean(trl_peaks_early(valid_s_early));
        trials_median_early(cond, subj) = median(trl_peaks_early(valid_s_early));
        trials_mean_late(cond, subj) = robust_trial_mean(trl_peaks_late(valid_s_late));
        trials_median_late(cond, subj) = median(trl_peaks_late(valid_s_late));

        % Peak power: highest dB value in the smoothed trial spectrum.
        trials_gamma_power(cond, subj) = robust_trial_mean(trial_peak_power_full);
        trials_gamma_power_early(cond, subj) = robust_trial_mean(trial_peak_power_early);
        trials_gamma_power_late(cond, subj) = robust_trial_mean(trial_peak_power_late);

        % Time-split centroid summaries.
        valid_cent_early = isfinite(trl_centroid_early);
        if any(valid_cent_early)
            trials_median_centroid_early(cond, subj) = median(trl_centroid_early(valid_cent_early));
        end
        if sum(valid_s_early) >= 2
            vf_early = trl_peaks_early(valid_s_early);
            trials_trialcv_early(cond, subj) = std(vf_early) / abs(mean(vf_early));
        end

        valid_cent_late = isfinite(trl_centroid_late);
        if any(valid_cent_late)
            trials_median_centroid_late(cond, subj) = median(trl_centroid_late(valid_cent_late));
        end
        if sum(valid_s_late) >= 2
            vf_late = trl_peaks_late(valid_s_late);
            trials_trialcv_late(cond, subj) = std(vf_late) / abs(mean(vf_late));
        end

    end % condition loop

    %  PER-SUBJECT FIGURES (one each for full, early, late; raw spectra)
    close all
    cmap_div = interp1([0 0.5 1], ...
        [0.17 0.27 0.53; 0.97 0.97 0.97; 0.70 0.09 0.17], linspace(0,1,256));

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

    window_names = {'full', 'early', 'late'};
    for wi = 1:3
        if wi == 1
            pr_source = subj_powratio_fullscan;
            peaks_source = subj_peaks_full;
            centroid_source = subj_centroid_full;
            condavg_source = subj_condition_avg_full;
            condpeak_source = subj_condition_peak_full;
            topo_mat_window = searchTopos_full;
            selected_idx_window = selected_idx_full;
            selected_w_window = w_combined_full;
            eigvals_window = evals_sorted_full;
        elseif wi == 2
            pr_source = subj_powratio_early;
            peaks_source = subj_peaks_early;
            centroid_source = subj_centroid_early;
            condavg_source = subj_condition_avg_early;
            condpeak_source = subj_condition_peak_early;
            topo_mat_window = searchTopos_early;
            selected_idx_window = selected_idx_early;
            selected_w_window = w_combined_early;
            eigvals_window = evals_sorted_early;
        else
            pr_source = subj_powratio_late;
            peaks_source = subj_peaks_late;
            centroid_source = subj_centroid_late;
            condavg_source = subj_condition_avg_late;
            condpeak_source = subj_condition_peak_late;
            topo_mat_window = searchTopos_late;
            selected_idx_window = selected_idx_late;
            selected_w_window = w_combined_late;
            eigvals_window = evals_sorted_late;
        end

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
        sgtitle(sprintf('Trial-Level GED: Subject %s (%s)', subjects{subj}, window_names{wi}), ...
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
        if ~isempty(topo_mat_window) && ~isempty(selected_idx_window)
            topo_data = [];
            topo_data.label  = dataEEG_c25.label;
            w_plot_common = selected_w_window(:);
            if numel(w_plot_common) ~= numel(selected_idx_window) || ~any(isfinite(w_plot_common))
                w_plot_common = ones(numel(selected_idx_window), 1);
            end
            w_plot_common(~isfinite(w_plot_common) | w_plot_common <= 0) = 0;
            if sum(w_plot_common) <= 0
                w_plot_common = ones(numel(selected_idx_window), 1);
            end
            w_plot_common = w_plot_common / sum(w_plot_common);
            topo_plot_common = topo_mat_window(:, selected_idx_window) * w_plot_common;
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
            n_sel_show = numel(selected_idx_window);
            lambda_show = NaN;
            if ~isempty(eigvals_window) && ~isempty(selected_idx_window)
                idx_show = selected_idx_window(1);
                if idx_show >= 1 && idx_show <= numel(eigvals_window)
                    lambda_show = eigvals_window(idx_show);
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
            sprintf('GCP_eeg_GED_subj%s_trials_overview_%s.png', subjects{subj}, window_names{wi})));
    end
    subject_runtime_seconds(subj) = toc(subj_runtime_tic);
end % subject loop

% CENTROID METRIC: Subject/group summaries and concordance
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

slope_post = compute_condition_separation_from_matrix(all_condition_peak_freq_full);
delta_post = all_condition_peak_freq_full(4, :) - all_condition_peak_freq_full(1, :);
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

dat_freq = all_condition_peak_freq_full;  % [condition x subject], peak from condition-averaged spectra
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

%% Frequency figure: gamma frequency over contrast by time window
fig_main_gamma_windows = figure('Position', [0 0 1512 982], 'Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot_gamma_window_panel(all_condition_peak_freq_early, condLabels, colors, nSubj);
ylabel('Gamma Frequency [Hz]');
title('Early (0-500 ms)', 'FontWeight', 'bold');

nexttile; hold on;
plot_gamma_window_panel(all_condition_peak_freq_full, condLabels, colors, nSubj);
title('Full (0-2000 ms)', 'FontWeight', 'bold');

nexttile; hold on;
plot_gamma_window_panel(all_condition_peak_freq_late, condLabels, colors, nSubj);
title('Late (1000-2000 ms)', 'FontWeight', 'bold');
save_figure_png(fig_main_gamma_windows, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_freq_windows.png'));

%% Power figure: gamma power over contrast by time window
fig_main_power_windows = figure('Position', [0 0 1512 982], 'Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot_gamma_window_panel(all_condition_peak_power_early, condLabels, colors, nSubj);
ylabel('Gamma Peak Power [dB]');
title('Early (0-500 ms)', 'FontWeight', 'bold');

nexttile; hold on;
plot_gamma_window_panel(all_condition_peak_power_full, condLabels, colors, nSubj);
title('Full (0-2000 ms)', 'FontWeight', 'bold');

nexttile; hold on;
plot_gamma_window_panel(all_condition_peak_power_late, condLabels, colors, nSubj);
title('Late (1000-2000 ms)', 'FontWeight', 'bold');
save_figure_png(fig_main_power_windows, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_power_windows.png'));

%% Condition-shift figure: normalized peak power trajectories
fig_condition_shift_power = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

dat_power = all_condition_peak_power_full;
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

%% All-subjects FULL combined topography + spectrum overview
cfg_topo_all = [];
cfg_topo_all.layout    = headmodel.layANThead;
cfg_topo_all.comment   = 'no';
cfg_topo_all.marker    = 'off';
cfg_topo_all.style     = 'straight';
cfg_topo_all.gridscale = 300;
cfg_topo_all.zlim      = 'maxabs';
cfg_topo_all.colormap  = '*RdBu';
cfg_topo_all.figure    = 'gcf';
plot_all_subjects_full_combined_topo_spectra( ...
    fig_save_dir_component_selection_root, subjects, scan_freqs, analysis_freq_range, cfg_topo_all, ...
    all_topo_labels, all_topos, all_combined_spectrum_full, all_combined_eigenvalue_full);

%% Save results
save_path = fullfile(gcp_root_path, 'data', 'features', 'GCP_eeg_GED.mat');
save(save_path, ...
    'trials_powratio', ...
    'trials_powratio_fullscan', ...
    'trials_powratio_early', 'trials_powratio_late', ...
    'trials_peaks', 'trials_centroid', ...
    'trials_peaks_early', 'trials_peaks_late', ...
    'trials_mean', 'trials_median', ...
    'trials_mean_early', 'trials_median_early', ...
    'trials_mean_late', 'trials_median_late', ...
    'trials_median_centroid_early', 'trials_median_centroid_late', ...
    'trials_trialcv_early', 'trials_trialcv_late', ...
    'trials_mean_centroid', 'trials_median_centroid', ...
    'trials_gamma_power', 'trials_gamma_power_early', 'trials_gamma_power_late', ...
    'all_topos', 'all_topos_early', 'all_topos_late', 'all_topo_labels', ...
    'all_combined_spectrum_full', 'all_combined_eigenvalue_full', 'all_eigenvalues', ...
    'all_selected_comp_idx', 'all_selected_comp_corr', 'all_selected_comp_eval', ...
    'all_selected_comp_indices_multi', 'all_selected_comp_weights', ...
    'all_component_selection_stats_full', 'all_component_selection_stats_early', 'all_component_selection_stats_late', ...
    'all_combined_filter_full', 'all_combined_filter_early', 'all_combined_filter_late', ...
    'all_haufe_pattern_full', 'all_haufe_patterns_multicomp_full', ...
    'all_reconstruction_patterns_full', 'freq_reconstructed_multicomp_full', ...
    'trials_powratio_components_full', 'trials_powratio_components_early', 'trials_powratio_components_late', ...
    'all_condition_powspctrm_full', 'all_condition_powspctrm_early', 'all_condition_powspctrm_late', ...
    'all_condition_powspctrm_full_unsmoothed', 'all_condition_powspctrm_early_unsmoothed', 'all_condition_powspctrm_late_unsmoothed', ...
    'all_condition_peak_freq_full', 'all_condition_peak_freq_early', 'all_condition_peak_freq_late', ...
    'all_condition_peak_power_full', 'all_condition_peak_power_early', 'all_condition_peak_power_late', ...
    'trials_outlier_mask_freq_full', 'trials_outlier_mask_freq_early', 'trials_outlier_mask_freq_late', ...
    'trials_outlier_mask_power_full', 'trials_outlier_mask_power_early', 'trials_outlier_mask_power_late', ...
    'all_top5_corrs', 'all_top5_evals', 'all_top5_topos', 'all_simulated_templates', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames', ...
    'baseline_window', 'full_window', 'early_window', 'late_window');

clc
powspctrm_save_path = fullfile(gcp_root_path, 'data', 'features', 'GCP_eeg_powspctrm_GED.mat');
% freq_* cells: one FieldTrip freq per (cond, subj), chan_freq, for ft_freqgrandaverage like AOC powl2{subj}.
% Do not add topolabel to those structs: FieldTrip ft_datatype sets iscomp if topolabel is present.
save(powspctrm_save_path, ...
    'freq_powspctrm_full', 'freq_powspctrm_full_unsmoothed', ...
    'all_condition_peak_freq_full', 'scan_freqs', 'condLabels', 'subjects');
fprintf('[GED] Subject-level freq (full-window GED spectra) saved to: %s\n', powspctrm_save_path);

%% Save GED analysis cohort (subjects with valid gamma power)
SubjID = str2double(string(subjects(:)));
Include = any(isfinite(trials_gamma_power), 1)';
subject_inclusion = table(SubjID, Include, 'VariableNames', {'SubjID', 'Include'});
save(fullfile(paths.controls, 'GCP_subject_inclusion.mat'), 'subject_inclusion', '-v7.3');

%% GED-projected TFR
if do_tfr
    [tfr_cond_trials, tfr_cond_avg, ged_filter_meta] = compute_ged_tfr( ...
        subjects, paths, all_combined_filter_full, all_topo_labels, ...
        all_component_selection_stats_full, ...
        baseline_window, ...
        tfr_foi, tfr_toi, tfr_win_sec, tfr_tapsmofrq, ...
        condNames, condCodes);
    tfr_save_path = fullfile(gcp_feature_data_path, 'GCP_eeg_GED_TFR.mat');
    save(tfr_save_path, ...
        'tfr_cond_trials', 'tfr_cond_avg', ...
        'ged_filter_meta', ...
        'subjects', 'condNames', 'condLabels', ...
        'baseline_window', 'tfr_baseline_window', 'full_window', ...
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

