%% GCP Gamma Peak Frequency and Power with Generalized Eigendecomposition (GED)
%
% Broadband GED and component selection (per subject)
%   - Pool trials across conditions and compute gamma-band covariances
%      (30-90 Hz FIR) for three windows (full/early/late) and baseline.
%   - Solve window-specific GED (S_stim * w = lambda * S_base * w) with
%      regularisation and rank candidate components by eigenvalue.
%   - Score candidates using occipital evidence, spectral powspctrm form
%      quality, and artifact metrics. Exclude extreme outliers and build
%      eigenvalue-weighted combined component.
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
%   - Peak frequency and power summaries (full, early, late windows).
%   - Detectability, trial CV, peak power, and condition separation.
%   - Subject diagnostics (component selection, rejection reasons,
%     topographies, spectra), group summary figures.

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('GCP');
nSubj = length(subjects);
total_runtime_tic = tic;

%% Parameters

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
min_eigval = 1.1;           % minimum GED eigenvalue (lambda >= 1.1)
min_powspctrm_form = 0.8;        % minimum PF (powspctrm-form) score for candidate eligibility
random_seed = 123;
powratio_trial_freq_smooth_bins = 5;      % movmean length (frequency bins) on per-trial powratio for peak/centroid
powratio_condition_freq_smooth_bins = 1;  % movmean length on condition-mean powratio before condition-level peaks/plots
peak_power_halfwidth_hz = 5;  % peak power = mean power within peak_hz +/-

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
    nSearch_plan = min(ged_search_n, nChans);
    n_freq_calls_est = 3 * nSearch_plan + 3; % 3 windows of proxy spectra + up to 3 trial-scan batches
    ged_freq_progress_reset(subjects{subj}, subj, nSubj, n_freq_calls_est);

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
    rejection_flags_full = struct(); rejection_flags_early = struct(); rejection_flags_late = struct();
    warn_flags_full = struct(); warn_flags_early = struct(); warn_flags_late = struct();
    eligible_full = []; eligible_early = []; eligible_late = [];
    front_leak_full = []; front_leak_early = []; front_leak_late = [];
    temp_leak_full = []; temp_leak_early = []; temp_leak_late = [];
    combined_leak_full = []; combined_leak_early = []; combined_leak_late = [];
    lineharm_full = []; lineharm_early = []; lineharm_late = [];
    hf_slope_full = []; hf_slope_early = []; hf_slope_late = [];
    emg_artifact_score_full = []; emg_artifact_score_early = []; emg_artifact_score_late = [];
    extreme_component_outlier_full = []; extreme_component_outlier_early = []; extreme_component_outlier_late = [];
    powspctrm_form_score_full = []; powspctrm_form_score_early = []; powspctrm_form_score_late = [];
    powspctrm_form_deduct_full = []; powspctrm_form_deduct_early = []; powspctrm_form_deduct_late = [];
    thr_max_combined_leak_full = NaN; thr_max_combined_leak_early = NaN; thr_max_combined_leak_late = NaN;
    thr_max_hf_slope_full = NaN; thr_max_hf_slope_early = NaN; thr_max_hf_slope_late = NaN;
    thr_max_emg_score_full = NaN; thr_max_emg_score_early = NaN; thr_max_emg_score_late = NaN;
    thr_max_lineharm_full = NaN; thr_max_lineharm_early = NaN; thr_max_lineharm_late = NaN;

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
        ged_freq_progress_set_phase(sprintf('%s ProxySpectra Batch', win_names_cap{w}));
        searchFilters = nan(nChans, nSearch);
        searchTopos = nan(nChans, nSearch);
        searchCorrs = nan(nSearch, 1);
        searchOccStrength = nan(nSearch, 1);
        searchFrontStrength = nan(nSearch, 1);
        searchOccFrontRatio = nan(nSearch, 1);
        searchFrontLeak = nan(nSearch, 1);
        searchTempLeak = nan(nSearch, 1);
        searchLineHarmRatio = nan(nSearch, 1);
        searchHFSlope = nan(nSearch, 1);
        searchEmgClass = repmat({'unassigned'}, nSearch, 1);
        searchMeanPrSpectrum = nan(nSearch, numel(scan_freqs));
        searchTopoPosteriorConcentration = nan(nSearch, 1);

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
            temp_strength = mean(abs(topo_ci(temp_idx)));
            temp_leak_ci = temp_strength / max(occ_strength, eps);
            topo_abs = abs(topo_ci(:));
            topo_finite = topo_abs(isfinite(topo_abs));
            if isempty(topo_finite)
                topo_finite = 0;
            end
            posterior_mass = mean(topo_abs(post_idx));
            total_mass = mean(topo_finite);
            topo_posterior_concentration_ci = posterior_mass / max(total_mass, eps);
            % Lightweight artifact proxies from trial-level component spectra
            proxy_ci = estimate_component_artifact_proxies( ...
                w_ci, dat_per_cond, stim_windows{w}, baseline_window, fsample, scan_freqs, mtmfft_tapsmofrq_hz);

            searchFilters(:, ci) = w_ci;
            searchTopos(:, ci) = topo_ci;
            searchCorrs(ci) = r_ci;
            searchOccStrength(ci) = occ_strength;
            searchFrontStrength(ci) = front_strength;
            searchOccFrontRatio(ci) = ratio_ci;
            searchFrontLeak(ci) = front_leak_ci;
            searchTempLeak(ci) = temp_leak_ci;
            searchLineHarmRatio(ci) = proxy_ci.lineharm_ratio;
            searchHFSlope(ci) = proxy_ci.hf_slope;
            searchMeanPrSpectrum(ci, :) = proxy_ci.mean_pr_spectrum(:)';
            searchTopoPosteriorConcentration(ci) = topo_posterior_concentration_ci;
        end

        % Candidate components metrics
        max_combined_leak = 1.30;             % artifact guard: mean(front leak, temporal leak)
        max_lineharm_ratio = 0.60;            % artifact guard: line-harmonic dominance ratio
        max_hf_slope = -0.15;                 % informative only (diagnostics); not an eligibility gate
        max_emg_score = 0.85;                 % anchor: very high EMG score
        max_components_to_combine = 10;       % top-K cap for combined GED branch
        outlier_ratio_thr = 3.0;              % lambda1/lambda2 threshold for extreme-component outlier detection
        outlier_mad_mult = 4.0;               % MAD multiplier on log-eigenvalue distance
        topo_nonposterior_max = 0.28;         % posterior concentration below this is non-posterior
        occ_class_thr = 0.60;                 % occipital-evidence threshold for occipital class
        emg_class_thr = 0.50;                 % EMG-score threshold for EMG class
        min_occ_margin = 0.05;                % occipital-vs-EMG margin threshold
        corr_vec = searchCorrs;
        ratio_vec = searchOccFrontRatio;
        eval_raw_vec = evals_sorted(1:nSearch);
        leak_vec = searchFrontLeak;
        temp_leak_vec = searchTempLeak;
        combined_leak_vec = 0.5 * (leak_vec + temp_leak_vec);
        lineharm_vec = searchLineHarmRatio;
        hf_slope_vec = searchHFSlope;
        hf_slope_for_score = hf_slope_vec;
        hf_slope_for_score(~isfinite(hf_slope_for_score)) = 0;
        finite_metrics = isfinite(corr_vec) & isfinite(ratio_vec) & ...
            isfinite(eval_raw_vec) & isfinite(leak_vec) & isfinite(temp_leak_vec) & ...
            isfinite(combined_leak_vec) & isfinite(lineharm_vec);
        peak_bonus_vec = compute_peak_bonus_from_spectra( ...
            searchMeanPrSpectrum, scan_freqs, analysis_freq_range);
        [powspctrm_form_score_vec, powspctrm_form_deduct_vec] = compute_powspctrm_form_laplacian_score_from_spectra( ...
            searchMeanPrSpectrum, scan_freqs, analysis_freq_range);
        topo_posterior_vec = searchTopoPosteriorConcentration;
        topo_nonposterior_fail_vec = topo_posterior_vec < topo_nonposterior_max;
        occipital_evidence = 0.40 * normalize_robust(corr_vec) + ...
            0.25 * normalize_robust(ratio_vec) + ...
            0.35 * normalize_robust(topo_posterior_vec);
        emg_artifact_score = 0.30 * normalize_robust(leak_vec) + ...
            0.20 * normalize_robust(temp_leak_vec) + ...
            0.30 * normalize_robust(lineharm_vec) + ...
            0.20 * normalize_robust(max(hf_slope_for_score, 0));
        enforced_min_eigval = min_eigval;
        pass_eig_gate = finite_metrics & (eval_raw_vec >= enforced_min_eigval);
        pass_peak_gate = finite_metrics & (powspctrm_form_score_vec >= min_powspctrm_form);
        fail_leak = finite_metrics & (combined_leak_vec > max_combined_leak);
        fail_lineharm = finite_metrics & (lineharm_vec > max_lineharm_ratio);
        fail_hf_slope = false(nSearch, 1); % HF slope kept in EMG score only; no eligibility gate
        fail_emg_score = finite_metrics & (emg_artifact_score > max_emg_score);
        hard_artifact_flags = fail_leak | fail_emg_score;
        artifact_flags = fail_leak | fail_emg_score | ...
            topo_nonposterior_fail_vec;
        warn_any = fail_lineharm;
        warn_flags = struct( ...
            'combined_leak', false(nSearch, 1), ...
            'lineharm', fail_lineharm, ...
            'hf_slope', false(nSearch, 1), ...
            'emg_score', false(nSearch, 1), ...
            'any', warn_any);
        rejection_flags = struct( ...
            'combined_leak', fail_leak, ...
            'lineharm', false(nSearch, 1), ...
            'hf_slope', fail_hf_slope, ...
            'emg_score', fail_emg_score, ...
            'topo_nonposterior', topo_nonposterior_fail_vec);
        occ_minus_emg_vec = occipital_evidence - emg_artifact_score;
        for ci = 1:nSearch
            occ_minus_emg_ci = occ_minus_emg_vec(ci);
            if (occipital_evidence(ci) >= occ_class_thr) && ...
                    (emg_artifact_score(ci) < emg_class_thr) && ...
                    (occ_minus_emg_ci >= min_occ_margin)
                searchEmgClass{ci} = 'occipital';
            elseif (occipital_evidence(ci) < occ_class_thr) && ...
                    (emg_artifact_score(ci) >= emg_class_thr) && ...
                    (-occ_minus_emg_ci >= min_occ_margin)
                searchEmgClass{ci} = 'EMG';
            elseif (occipital_evidence(ci) >= occ_class_thr) && ...
                    (emg_artifact_score(ci) >= emg_class_thr)
                if abs(occ_minus_emg_ci) < min_occ_margin && ...
                        (topo_posterior_vec(ci) >= topo_nonposterior_max)
                    searchEmgClass{ci} = 'mixed';
                elseif occ_minus_emg_ci >= 0
                    searchEmgClass{ci} = 'occipital';
                else
                    searchEmgClass{ci} = 'EMG';
                end
            else
                searchEmgClass{ci} = 'mixed';
            end
        end
        extreme_component_outlier_mask = false(nSearch, 1);
        occipital_class_mask = cellfun(@(c) strcmpi(c, 'occipital'), searchEmgClass(:));
        mixed_class_mask = cellfun(@(c) strcmpi(c, 'mixed'), searchEmgClass(:));
        mixed_rescue_mask = mixed_class_mask & ...
            ~hard_artifact_flags;
        eligible_class_mask = (~artifact_flags & occipital_class_mask) | mixed_rescue_mask;
        raw_eligible_for_outlier = pass_eig_gate & pass_peak_gate & ...
            eligible_class_mask;
        [~, extreme_component_outlier_idx] = exclude_extreme_component_outlier( ...
            evals_sorted(1:nSearch), raw_eligible_for_outlier, outlier_ratio_thr, outlier_mad_mult);
        if ~isempty(extreme_component_outlier_idx)
            extreme_component_outlier_mask(extreme_component_outlier_idx) = true;
        end
        eligible = pass_eig_gate & pass_peak_gate & ...
            ~extreme_component_outlier_mask & eligible_class_mask;
        no_threshold_match = ~any(eligible);
        selection_pool_mask = eligible;
        searchScores = compute_calibrated_rank_aggregation_score( ...
            eval_raw_vec, powspctrm_form_score_vec, peak_bonus_vec, occipital_evidence, emg_artifact_score);
        searchScores(~finite_metrics) = -Inf;
        searchScores(~selection_pool_mask) = -Inf;
        [bestScore, bestIdx] = max(searchScores);
        if isempty(bestIdx) || isnan(bestScore)
            bestIdx = 1;
            bestScore = NaN;
        end

        combined_idx = find(selection_pool_mask & isfinite(searchScores));
        if isempty(combined_idx)
            combined_weights = [];
        else
            [~, combined_ord] = sort(searchScores(combined_idx), 'descend');
            if ~any(isfinite(searchScores(combined_idx)))
                [~, combined_ord] = sort(evals_sorted(combined_idx), 'descend');
            end
            combined_idx = combined_idx(combined_ord);
            combined_idx = combined_idx(1:min(max_components_to_combine, numel(combined_idx)));

            combined_weights = evals_sorted(combined_idx)';
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
            bestRatio = searchOccFrontRatio(bestIdx);
            bestLeak = searchFrontLeak(bestIdx);

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
            rejection_flags_full = rejection_flags;
            warn_flags_full = warn_flags;
            eligible_full = eligible;
            front_leak_full = leak_vec;
            temp_leak_full = temp_leak_vec;
            combined_leak_full = combined_leak_vec;
            lineharm_full = lineharm_vec;
            hf_slope_full = hf_slope_vec;
            emg_artifact_score_full = emg_artifact_score;
            extreme_component_outlier_full = extreme_component_outlier_mask;
            powspctrm_form_score_full = powspctrm_form_score_vec;
            powspctrm_form_deduct_full = powspctrm_form_deduct_vec;
            thr_max_combined_leak_full = max_combined_leak;
            thr_max_hf_slope_full = max_hf_slope;
            thr_max_emg_score_full = max_emg_score;
            thr_max_lineharm_full = max_lineharm_ratio;
        elseif w == 2
            searchFilters_early = searchFilters;
            selected_idx_early = selected_idx;
            w_combined_early = selected_weights(:)';
            evals_sorted_early = evals_sorted;
            searchCorrs_early = searchCorrs;
            searchTopos_early = searchTopos;
            searchMeanPrSpectrum_early = searchMeanPrSpectrum;
            searchEmgClass_early = searchEmgClass;
            rejection_flags_early = rejection_flags;
            warn_flags_early = warn_flags;
            eligible_early = eligible;
            front_leak_early = leak_vec;
            temp_leak_early = temp_leak_vec;
            combined_leak_early = combined_leak_vec;
            lineharm_early = lineharm_vec;
            hf_slope_early = hf_slope_vec;
            emg_artifact_score_early = emg_artifact_score;
            extreme_component_outlier_early = extreme_component_outlier_mask;
            powspctrm_form_score_early = powspctrm_form_score_vec;
            powspctrm_form_deduct_early = powspctrm_form_deduct_vec;
            thr_max_combined_leak_early = max_combined_leak;
            thr_max_hf_slope_early = max_hf_slope;
            thr_max_emg_score_early = max_emg_score;
            thr_max_lineharm_early = max_lineharm_ratio;
        else
            searchFilters_late = searchFilters;
            selected_idx_late = selected_idx;
            w_combined_late = selected_weights(:)';
            evals_sorted_late = evals_sorted;
            searchCorrs_late = searchCorrs;
            searchTopos_late = searchTopos;
            searchMeanPrSpectrum_late = searchMeanPrSpectrum;
            searchEmgClass_late = searchEmgClass;
            rejection_flags_late = rejection_flags;
            warn_flags_late = warn_flags;
            eligible_late = eligible;
            front_leak_late = leak_vec;
            temp_leak_late = temp_leak_vec;
            combined_leak_late = combined_leak_vec;
            lineharm_late = lineharm_vec;
            hf_slope_late = hf_slope_vec;
            emg_artifact_score_late = emg_artifact_score;
            extreme_component_outlier_late = extreme_component_outlier_mask;
            powspctrm_form_score_late = powspctrm_form_score_vec;
            powspctrm_form_deduct_late = powspctrm_form_deduct_vec;
            thr_max_combined_leak_late = max_combined_leak;
            thr_max_hf_slope_late = max_hf_slope;
            thr_max_emg_score_late = max_emg_score;
            thr_max_lineharm_late = max_lineharm_ratio;
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
            'best_ratio', bestRatio, ...
            'best_front', bestFront, ...
            'best_occ', bestOcc, ...
            'best_leak', bestLeak, ...
            'emg_artifact_score', emg_artifact_score, ...
            'occipital_evidence', occipital_evidence, ...
            'occ_minus_emg', occ_minus_emg_vec, ...
            'emg_class', {searchEmgClass}, ...
            'mixed_rescue_mask', mixed_rescue_mask, ...
            'reject_flags', artifact_flags, ...
            'warn_flags', warn_flags, ...
            'rejection_flags', rejection_flags, ...
            'no_threshold_match', no_threshold_match);
        if w == 1
            all_component_selection_stats_full{subj} = comp_sel_struct;
        elseif w == 2
            all_component_selection_stats_early{subj} = comp_sel_struct;
        else
            all_component_selection_stats_late{subj} = comp_sel_struct;
        end
    end

    W_combined_full = searchFilters_full(:, selected_idx_full);
    W_combined_early = searchFilters_early(:, selected_idx_early);
    W_combined_late = searchFilters_late(:, selected_idx_late);

    topo_temp_full = searchTopos_full(:, selected_idx_full) * w_combined_full(:);
    topo_temp_early = searchTopos_early(:, selected_idx_early) * w_combined_early(:);
    topo_temp_late = searchTopos_late(:, selected_idx_late) * w_combined_late(:);
    all_topos{subj} = topo_temp_full;
    all_topos_early{subj} = topo_temp_early;
    all_topos_late{subj} = topo_temp_late;
    all_selected_comp_indices_multi{subj} = selected_idx_full;
    all_selected_comp_weights{subj} = w_combined_full(:)';
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
        rejection_flags_full, ...
        front_leak_full, temp_leak_full, combined_leak_full, ...
        lineharm_full, hf_slope_full, emg_artifact_score_full, ...
        extreme_component_outlier_full, ...
        cfg_topo, all_topo_labels{subj}, powspctrm_form_score_full, powspctrm_form_deduct_full, ...
        thr_max_combined_leak_full, thr_max_hf_slope_full, ...
        thr_max_emg_score_full, thr_max_lineharm_full, ...
        selected_idx_full);
    plot_emg_exclusion_diagnostics( ...
        fig_save_dir_component_selection, subjects{subj}, 'early', scan_freqs, searchTopos_early, ...
        searchMeanPrSpectrum_early, evals_sorted_early(1:numel(eligible_early)), ...
        searchEmgClass_early, ...
        eligible_early, ...
        rejection_flags_early, ...
        front_leak_early, temp_leak_early, combined_leak_early, ...
        lineharm_early, hf_slope_early, emg_artifact_score_early, ...
        extreme_component_outlier_early, ...
        cfg_topo, all_topo_labels{subj}, powspctrm_form_score_early, powspctrm_form_deduct_early, ...
        thr_max_combined_leak_early, thr_max_hf_slope_early, ...
        thr_max_emg_score_early, thr_max_lineharm_early, ...
        selected_idx_early);
    plot_emg_exclusion_diagnostics( ...
        fig_save_dir_component_selection, subjects{subj}, 'late', scan_freqs, searchTopos_late, ...
        searchMeanPrSpectrum_late, evals_sorted_late(1:numel(eligible_late)), ...
        searchEmgClass_late, ...
        eligible_late, ...
        rejection_flags_late, ...
        front_leak_late, temp_leak_late, combined_leak_late, ...
        lineharm_late, hf_slope_late, emg_artifact_score_late, ...
        extreme_component_outlier_late, ...
        cfg_topo, all_topo_labels{subj}, powspctrm_form_score_late, powspctrm_form_deduct_late, ...
        thr_max_combined_leak_late, thr_max_hf_slope_late, ...
        thr_max_emg_score_late, thr_max_lineharm_late, ...
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
        trial_print_step = 25;
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
            ged_freq_progress_set_phase('Full TrialScan Batch');
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
                if mod(trl-1, trial_print_step) == 0 || trl == nTrl
                    clc
                    fprintf('[GED] Subject %s (%d/%d) Trial %d/%d Full Window\n', ...
                        subjects{subj}, subj, nSubj, trl, nTrl);
                end
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
            clc;
            fprintf('[GED] Subject %s (%d/%d) Trial 0/%d Early Window\n', ...
                subjects{subj}, subj, nSubj, nTrl);
            trial_mask_early = has_base & has_early & ~bad_base_early;
            ged_freq_progress_set_phase('Early TrialScan Batch');
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
                if mod(trl-1, trial_print_step) == 0 || trl == nTrl
                    fprintf('[GED] Subject %s (%d/%d) Trial %d/%d Early Window\n', ...
                        subjects{subj}, subj, nSubj, trl, nTrl);
                end
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
            clc;
            fprintf('[GED] Subject %s (%d/%d) Trial 0/%d Late Window\n', ...
                subjects{subj}, subj, nSubj, nTrl);
            trial_mask_late = has_base & has_late & ~bad_base_late;
            ged_freq_progress_set_phase('Late TrialScan Batch');
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
                if mod(trl-1, trial_print_step) == 0 || trl == nTrl
                    fprintf('[GED] Subject %s (%d/%d) Trial %d/%d Late Window\n', ...
                        subjects{subj}, subj, nSubj, trl, nTrl);
                end
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
    'all_topos', 'all_topos_early', 'all_topos_late', 'all_topo_labels', 'all_eigenvalues', ...
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

function [trl_peaks, trl_peak_power, trl_centroid] = ...
    compute_trial_peak_metrics_from_powratio_fullscan(powratio_trials_fullscan, scan_freqs_full, analysis_mask, ...
    smooth_n, peak_power_halfwidth_hz)
nTrl = size(powratio_trials_fullscan, 1);
trl_peaks = nan(nTrl, 1);
trl_peak_power = nan(nTrl, 1);
trl_centroid = nan(nTrl, 1);
scan_freqs_analysis = scan_freqs_full(analysis_mask);
centroid_band_mask = scan_freqs_full >= 30 & scan_freqs_full <= 90;
freq_band = scan_freqs_full(centroid_band_mask);
if ~isfinite(peak_power_halfwidth_hz) || peak_power_halfwidth_hz < 0
    peak_power_halfwidth_hz = 0;
end
for trl = 1:nTrl
    pr_full = powratio_trials_fullscan(trl, :);
    if all(~isfinite(pr_full))
        continue;
    end
    pr_proc = pr_full(analysis_mask);
    pr_proc = pr_proc(:)';
    valid = isfinite(pr_proc) & isfinite(scan_freqs_analysis);
    pr_proc = pr_proc(valid);
    x_use = scan_freqs_analysis(valid);
    if isempty(pr_proc)
        continue;
    end

    % Peak search over the full 30-90 Hz gamma band.
    [peak_hz, peak_power] = pick_tallest_peak(pr_proc, x_use, smooth_n, peak_power_halfwidth_hz);
    if isfinite(peak_hz)
        trl_peaks(trl) = peak_hz;
        trl_peak_power(trl) = peak_power;
    end

    pr_proc_full = movmean(pr_full, max(1, round(smooth_n)), 'omitnan');
    pr_proc_band = pr_proc_full(centroid_band_mask);
    w_pos = max(pr_proc_band, 0);
    pos_mass = sum(w_pos);
    if pos_mass > 0
        trl_centroid(trl) = sum(freq_band .* w_pos) / pos_mass;
    end
end
end

function [peak_hz, peak_power] = pick_tallest_peak(y, x, smooth_n, peak_power_halfwidth_hz)
peak_hz = NaN;
peak_power = NaN;
y = y(:)';
x = x(:)';
if ~isfinite(smooth_n) || smooth_n < 1
    smooth_n = 1;
end
if ~isfinite(peak_power_halfwidth_hz) || peak_power_halfwidth_hz < 0
    peak_power_halfwidth_hz = 0;
end
if isempty(y) || numel(y) ~= numel(x)
    return;
end
y = movmean(y, max(1, round(smooth_n)), 'omitnan');
core_mask = x >= 30 & x <= 90 & isfinite(y) & isfinite(x);
if any(core_mask)
    x_use = x(core_mask);
    y_use = y(core_mask);
else
    valid = isfinite(y) & isfinite(x);
    if ~any(valid)
        return;
    end
    x_use = x(valid);
    y_use = y(valid);
end
if numel(x_use) < 3
    return;
end

% Use local maxima only; if none pass criteria, keep NaN.
dx = diff(x_use);
dx = dx(isfinite(dx) & dx > 0);
if isempty(dx)
    return;
end
freq_step = median(dx);
if ~isfinite(freq_step) || freq_step <= 0
    return;
end
min_peak_width_hz = max(2 * freq_step, 2.0);
local_spread = iqr(y_use);
if ~isfinite(local_spread) || local_spread <= 0
    local_spread = std(y_use, 'omitnan');
end
if ~isfinite(local_spread) || local_spread <= 0
    local_spread = max(abs(y_use), [], 'omitnan');
end
if ~isfinite(local_spread) || local_spread <= 0
    local_spread = 1;
end
min_peak_prom_db = max(0.15, 0.25 * local_spread);

[pks, locs] = findpeaks(y_use, x_use, ...
    'MinPeakProminence', min_peak_prom_db, ...
    'MinPeakWidth', min_peak_width_hz, ...
    'SortStr', 'descend', ...
    'NPeaks', 1);
if isempty(pks) || isempty(locs) || ~isfinite(pks(1)) || ~isfinite(locs(1))
    return;
end

peak_hz = locs(1);
peak_power = pks(1);
if peak_power_halfwidth_hz > 0
    band_mask = abs(x_use - peak_hz) <= peak_power_halfwidth_hz;
    band_power = y_use(band_mask);
    band_power = band_power(isfinite(band_power));
    if ~isempty(band_power)
        peak_power = mean(band_power);
    end
end

end

function [ratio_db, near_floor_freq_mask, near_floor_row_mask] = compute_scan_ratio_from_timeseries(sig_stim, sig_base, fs, scan_freqs, tapsmofrq_hz, base_floor, near_floor_mult)
if ~iscell(sig_stim) && isvector(sig_stim)
    sig_stim = sig_stim(:)';
end
if ~iscell(sig_base) && isvector(sig_base)
    sig_base = sig_base(:)';
end
if iscell(sig_stim)
    nSig = numel(sig_stim);
else
    nSig = size(sig_stim, 1);
end
ratio_db = nan(nSig, numel(scan_freqs));
near_floor_freq_mask = false(1, numel(scan_freqs));
near_floor_row_mask = false(nSig, numel(scan_freqs));
if isempty(sig_stim) || isempty(sig_base) || fs <= 0
    return;
end
if iscell(sig_stim) ~= iscell(sig_base)
    return;
end
if iscell(sig_stim)
    if numel(sig_stim) ~= numel(sig_base)
        return;
    end
else
    if size(sig_stim, 1) ~= size(sig_base, 1)
        return;
    end
end
if nSig == 0
    return;
end
if ~isfinite(near_floor_mult) || near_floor_mult <= 0
    near_floor_mult = 1.5;
end
% Each signal is represented as one virtual channel (GED component time series).
[p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(sig_stim, sig_base, fs, scan_freqs, tapsmofrq_hz);
if isempty(p_stim_scan) || isempty(p_base_scan)
    return;
end
% Derive a baseline floor for instability detection only. Valid power ratios
% use linear stimulus and baseline power directly before conversion to dB.
floor_fallback = max(base_floor, eps);
for ri = 1:nSig
    p_stim_row = p_stim_scan(ri, :);
    p_base_row = p_base_scan(ri, :);
    valid_base = isfinite(p_base_row) & (p_base_row > 0);
    if ~any(valid_base)
        continue;
    end
    base_anchor = prctile(p_base_row(valid_base), 20);
    base_median = median(p_base_row(valid_base), 'omitnan');
    if ~isfinite(base_anchor) || base_anchor <= 0
        base_anchor = base_median;
    end
    if ~isfinite(base_anchor) || base_anchor <= 0
        base_anchor = floor_fallback;
    end
    floor_row = max(0.25 * base_anchor, eps);

    valid = isfinite(p_stim_row) & isfinite(p_base_row) & ...
        (p_stim_row > 0) & (p_base_row > 0);
    ratio_row = nan(1, numel(scan_freqs));
    ratio_row(valid) = 10 * log10(p_stim_row(valid) ./ p_base_row(valid));
    ratio_db(ri, :) = ratio_row;
    near_floor_row_mask(ri, valid) = p_base_row(valid) <= near_floor_mult * floor_row;
end

for fi = 1:numel(scan_freqs)
    valid_fi = isfinite(ratio_db(:, fi));
    if any(valid_fi)
        near_floor_freq_mask(fi) = mean(near_floor_row_mask(valid_fi, fi)) >= 0.5;
    end
end
end

function [ratio_cube, near_floor_freq_count_per_trial] = compute_scan_ratio_for_window_batch(trial_cache, search_filters, stim_field, trial_mask, fs, scan_freqs, tapsmofrq_hz, base_floor, near_floor_mult)
nTrl = numel(trial_cache);
nComp = size(search_filters, 2);
nFreq = numel(scan_freqs);
ratio_cube = nan(nComp, nTrl, nFreq);
near_floor_freq_count_per_trial = zeros(nTrl, 1);
if nComp == 0 || nTrl == 0 || ~any(trial_mask)
    return;
end

sig_stim_cells = cell(0, 1);
sig_base_cells = cell(0, 1);
row_comp_idx = zeros(0, 1);
row_trial_idx = zeros(0, 1);
for trl = 1:nTrl
    if ~trial_mask(trl)
        continue;
    end
    tc = trial_cache{trl};
    x_base = tc.x_base;
    x_stim = tc.(stim_field);
    if isempty(x_base) || isempty(x_stim)
        continue;
    end
    comp_base = search_filters' * x_base;
    comp_stim = search_filters' * x_stim;
    for ci = 1:nComp
        sig_base_cells{end+1, 1} = comp_base(ci, :);
        sig_stim_cells{end+1, 1} = comp_stim(ci, :);
        row_comp_idx(end+1, 1) = ci;
        row_trial_idx(end+1, 1) = trl;
    end
end
if isempty(sig_stim_cells)
    return;
end

[ratio_rows, ~, near_floor_row_mask] = compute_scan_ratio_from_timeseries( ...
    sig_stim_cells, sig_base_cells, fs, scan_freqs, tapsmofrq_hz, base_floor, near_floor_mult);
for ri = 1:size(ratio_rows, 1)
    ratio_cube(row_comp_idx(ri), row_trial_idx(ri), :) = ratio_rows(ri, :);
end
for trl = 1:nTrl
    rows = row_trial_idx == trl;
    if ~any(rows)
        continue;
    end
    near_floor_trial_mask = mean(near_floor_row_mask(rows, :), 1) >= 0.5;
    near_floor_freq_count_per_trial(trl) = sum(near_floor_trial_mask);
end
end

function [base_floor, base_median] = compute_baseline_floor_stats(base_power_vals, floor_prctile, floor_frac)
base_floor = eps;
base_median = eps;
if ~isfinite(floor_prctile)
    floor_prctile = 20;
end
if ~isfinite(floor_frac)
    floor_frac = 0.25;
end
if isempty(base_power_vals)
    return;
end
valid = isfinite(base_power_vals) & (base_power_vals > 0);
if ~any(valid)
    return;
end
vals = base_power_vals(valid);
base_median = median(vals, 'omitnan');
anchor = prctile(vals, floor_prctile);
if ~isfinite(anchor) || anchor <= 0
    anchor = base_median;
end
if ~isfinite(base_median) || base_median <= 0
    base_median = anchor;
end
base_floor = max(anchor * max(floor_frac, 0), eps);
base_median = max(base_median, base_floor);
end

function ylims = compute_symmetric_ylim(vals, pad_frac, min_half_range)
if ~isfinite(pad_frac) || pad_frac < 0
    pad_frac = 0.15;
end
if ~isfinite(min_half_range) || min_half_range <= 0
    min_half_range = 1;
end
vals = vals(:);
vals = vals(isfinite(vals));
if isempty(vals)
    half_range = min_half_range;
else
    max_abs = max(abs(vals));
    if ~isfinite(max_abs)
        max_abs = min_half_range;
    end
    half_range = max(max_abs * (1 + pad_frac), min_half_range);
end
ylims = [-half_range, half_range];
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
med = nanmedian(dat, 2);
madv = nan(4, 1);
for c = 1:4
    madv(c) = robust_mad(dat(c, :));
end
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
errorbar(1:4, med, madv, 'k-o', 'LineWidth', 1.8, 'CapSize', 8, ...
    'MarkerFaceColor', [0.1 0.1 0.1], 'MarkerSize', 6);
set(gca, 'XTick', 1:4, 'XTickLabel', strcat(condLabels, ' Contrast'), 'FontSize', 11, 'Box', 'off');
xlim([0.5 4.5]);
end

function z = normalize_robust(x)
x = x(:);
z = nan(size(x));
valid = isfinite(x);
if ~any(valid)
    return;
end
xv = x(valid);
xm = median(xv);
xs = mad(xv, 1);
if ~isfinite(xs) || xs <= eps
    xs = std(xv);
end
if ~isfinite(xs) || xs <= eps
    xs = 1;
end
z(valid) = (xv - xm) / xs;
z(valid) = max(min(z(valid), 3), -3);
end

function rank_vec = compute_descending_rank(x)
rank_vec = nan(size(x));
valid = isfinite(x);
if ~any(valid)
    return;
end
[~, ord] = sort(x(valid), 'descend');
idx_valid = find(valid);
rank_vals = nan(sum(valid), 1);
rank_vals(ord) = 1:numel(ord);
rank_vec(idx_valid) = rank_vals;
end

function score_vec = compute_calibrated_rank_aggregation_score( ...
    eval_raw_vec, powspctrm_form_vec, peak_bonus_vec, occipital_evidence_vec, emg_artifact_vec)
nComp = numel(eval_raw_vec);
score_vec = -Inf(nComp, 1);
if nComp == 0
    return;
end

z_eig = normalize_robust(log(max(eval_raw_vec(:), eps)));
z_pf = normalize_robust(powspctrm_form_vec(:));
z_pb = normalize_robust(peak_bonus_vec(:));
z_occ = normalize_robust(occipital_evidence_vec(:));
z_anti_emg = normalize_robust(-emg_artifact_vec(:));
metric_mat = [z_eig, z_pf, z_pb, z_occ, z_anti_emg];
finite_rows = all(isfinite(metric_mat), 2);
if ~any(finite_rows)
    return;
end

point_mat = nan(nComp, size(metric_mat, 2));
for mi = 1:size(metric_mat, 2)
    rank_m = compute_descending_rank(metric_mat(:, mi));
    valid_rank = isfinite(rank_m);
    if ~any(valid_rank)
        continue;
    end
    n_valid = sum(valid_rank);
    if n_valid <= 1
        point_mat(valid_rank, mi) = 1;
    else
        point_mat(valid_rank, mi) = (n_valid - rank_m(valid_rank)) / (n_valid - 1);
    end
end
score_vec = mean(point_mat, 2, 'omitnan');
score_vec(~finite_rows) = -Inf;
end

function plot_emg_exclusion_diagnostics(save_dir, subject_id, win_name, scan_freqs, searchTopos, ...
    searchMeanPrSpectrum, eigval_vec, emg_class, ...
    eligible, rejection_flags, ...
    front_leak_vec, temp_leak_vec, combined_leak_vec, lineharm_vec, hf_slope_vec, emg_score_vec, ...
    extreme_component_outlier, ...
    cfg_topo, topo_labels, ...
    powspctrm_form_score, powspctrm_form_multi_peak_deduction, ...
    max_combined_leak_thr, max_hf_slope_thr, max_emg_score_thr, max_lineharm_ratio_thr, ...
    selected_idx)
nComp = numel(eigval_vec);
if nComp < 1
    return;
end
selected_mask = false(nComp, 1);
if ~isempty(selected_idx)
    selected_idx = unique(selected_idx(:));
    selected_idx = selected_idx(isfinite(selected_idx));
    selected_idx = selected_idx(selected_idx >= 1 & selected_idx <= nComp);
    selected_idx = round(selected_idx);
    selected_mask(selected_idx) = true;
end
final_eligible_mask = logical(eligible(:));
final_rejected_mask = ~final_eligible_mask;
sel_idx = find(selected_mask);
rej_idx = find(final_rejected_mask & ~selected_mask);
[~, so] = sort(eigval_vec(sel_idx), 'descend');
sel_idx = sel_idx(so);
[~, ro] = sort(eigval_vec(rej_idx), 'descend');
rej_idx = rej_idx(ro);

% Top selected/rejected components in one figure (max 5 per row):
% row 1 = selected topographies, row 2 = selected spectra,
% row 3 = rejected topographies, row 4 = rejected spectra.
nCols = 5;
nShowSel = min(nCols, numel(sel_idx));
nShowRej = min(nCols, numel(rej_idx));
figSel = figure('Position', [0 0 1512 982], 'Color', 'w');
for k = 1:nCols
    % Row 1: selected topographies
    subplot(4, nCols, k);
    if k <= nShowSel
        ci = sel_idx(k);
        topo_data = [];
        topo_data.label = topo_labels;
        topo_data.avg = searchTopos(:, ci);
        topo_data.dimord = 'chan';
        cfg_ci = cfg_topo;
        topo_vals = topo_data.avg(isfinite(topo_data.avg));
        topo_clim_ci = max(abs(topo_vals));
        if ~isfinite(topo_clim_ci) || topo_clim_ci <= 0
            topo_clim_ci = 1;
        end
        cfg_ci.zlim = [-topo_clim_ci topo_clim_ci];
        try
            ft_topoplotER(cfg_ci, topo_data);
        catch
            imagesc(topo_data.avg(:)); axis tight;
            caxis([-topo_clim_ci topo_clim_ci]);
        end
        ttl = sprintf('C%d (%s)', ci, emg_class{ci});
        title(ttl, 'FontSize', 8, 'Interpreter', 'none');
    else
        axis off;
    end

    % Row 2: selected spectra
    subplot(4, nCols, nCols + k); hold on;
    if k <= nShowSel
        ci = sel_idx(k);
        spec_data = searchMeanPrSpectrum(ci, :);
        plot(scan_freqs, spec_data, '-', 'Color', [0 0 0], 'LineWidth', 1.4);
        yline(0, 'k--');
        spec_min = min(spec_data(isfinite(spec_data)));
        spec_max = max(spec_data(isfinite(spec_data)));
        if isfinite(spec_min) && isfinite(spec_max) && spec_min ~= spec_max
            spec_range = spec_max - spec_min;
            ylim([spec_min - 0.10 * spec_range, spec_max + 0.10 * spec_range]);
        end
        [info_lines, info_viol] = ...
            build_rejection_info_columns(ci, rejection_flags, front_leak_vec, temp_leak_vec, combined_leak_vec, ...
            lineharm_vec, hf_slope_vec, emg_score_vec, extreme_component_outlier, ...
            max_combined_leak_thr, max_hf_slope_thr, max_emg_score_thr, max_lineharm_ratio_thr);
        plot_rejection_info_text_columns(info_lines, info_viol);
        format_power_change_db_axis(gca);
        xlabel('Hz'); ylabel('Power [dB]');
        pf_deduct_ci = 0;
        if ~isempty(powspctrm_form_multi_peak_deduction) && numel(powspctrm_form_multi_peak_deduction) >= ci && ...
                isfinite(powspctrm_form_multi_peak_deduction(ci))
            pf_deduct_ci = powspctrm_form_multi_peak_deduction(ci);
        end
        title(sprintf('\\lambda=%.2f, PF=%.2f (%.2f)', ...
            eigval_vec(ci), powspctrm_form_score(ci), pf_deduct_ci), 'FontSize', 7);
        box on;
    else
        axis off;
    end

    % Row 3: rejected topographies
    subplot(4, nCols, 2 * nCols + k);
    if k <= nShowRej
        ci = rej_idx(k);
        topo_data = [];
        topo_data.label = topo_labels;
        topo_data.avg = searchTopos(:, ci);
        topo_data.dimord = 'chan';
        cfg_ci = cfg_topo;
        topo_vals = topo_data.avg(isfinite(topo_data.avg));
        topo_clim_ci = max(abs(topo_vals));
        if ~isfinite(topo_clim_ci) || topo_clim_ci <= 0
            topo_clim_ci = 1;
        end
        cfg_ci.zlim = [-topo_clim_ci topo_clim_ci];
        try
            ft_topoplotER(cfg_ci, topo_data);
        catch
            imagesc(topo_data.avg(:)); axis tight;
            caxis([-topo_clim_ci topo_clim_ci]);
        end
        ttl = sprintf('C%d (%s)', ci, emg_class{ci});
        title(ttl, 'FontSize', 8, 'Interpreter', 'none');
    else
        axis off;
    end

    % Row 4: rejected spectra
    subplot(4, nCols, 3 * nCols + k); hold on;
    if k <= nShowRej
        ci = rej_idx(k);
        spec_data = searchMeanPrSpectrum(ci, :);
        plot(scan_freqs, spec_data, '-', 'Color', [0 0 0], 'LineWidth', 1.4);
        yline(0, 'k--');
        spec_min = min(spec_data(isfinite(spec_data)));
        spec_max = max(spec_data(isfinite(spec_data)));
        if isfinite(spec_min) && isfinite(spec_max) && spec_min ~= spec_max
            spec_range = spec_max - spec_min;
            ylim([spec_min - 0.10 * spec_range, spec_max + 0.10 * spec_range]);
        end
        [info_lines, info_viol] = ...
            build_rejection_info_columns(ci, rejection_flags, front_leak_vec, temp_leak_vec, combined_leak_vec, ...
            lineharm_vec, hf_slope_vec, emg_score_vec, extreme_component_outlier, ...
            max_combined_leak_thr, max_hf_slope_thr, max_emg_score_thr, max_lineharm_ratio_thr);
        plot_rejection_info_text_columns(info_lines, info_viol);
        format_power_change_db_axis(gca);
        xlabel('Hz'); ylabel('Power [dB]');
        pf_deduct_ci = 0;
        if ~isempty(powspctrm_form_multi_peak_deduction) && numel(powspctrm_form_multi_peak_deduction) >= ci && ...
                isfinite(powspctrm_form_multi_peak_deduction(ci))
            pf_deduct_ci = powspctrm_form_multi_peak_deduction(ci);
        end
        title(sprintf('\\lambda=%.2f, PF=%.2f (%.2f)', ...
            eigval_vec(ci), powspctrm_form_score(ci), pf_deduct_ci), 'FontSize', 7);
        box on;
    else
        axis off;
    end
end
if nShowSel == 0
    annotation(figSel, 'textbox', [0.25 0.77 0.50 0.06], ...
        'String', 'NO ELIGIBLE COMPONENTS', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 16, 'Color', [0.85 0.05 0.05], ...
        'Interpreter', 'none');
end
annotation(figSel, 'textbox', [0.010 0.56 0.040 0.34], ...
    'String', 'Selected Components', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontWeight', 'bold', 'FontSize', 14, 'Rotation', 90, 'FitBoxToText', 'off');
annotation(figSel, 'textbox', [0.010 0.09 0.040 0.34], ...
    'String', 'Rejected Components', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontWeight', 'bold', 'FontSize', 14, 'Rotation', 90, 'FitBoxToText', 'off');
sgtitle(sprintf('Top 5 Selected and Rejected Components: %s (%s)', subject_id, win_name), ...
    'FontSize', 14, 'FontWeight', 'bold');
drawnow;
pause(0.05);
save_figure_png(figSel, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_topo_spectra_selected_%s.png', subject_id, win_name)));
close(figSel);

end

function [info_lines, info_viol] = ...
    build_rejection_info_columns(ci, rejection_flags, front_leak_vec, temp_leak_vec, combined_leak_vec, ...
    lineharm_vec, hf_slope_vec, emg_score_vec, extreme_component_outlier, ...
    max_combined_leak_thr, max_hf_slope_thr, max_emg_score_thr, max_lineharm_ratio_thr)
combined_leak_val = combined_leak_vec(ci);
extreme_outlier_fail = logical(extreme_component_outlier(ci));
info_lines = { ...
    sprintf('extreme outlier: %d', extreme_outlier_fail), ...
    sprintf('combined_leak gate: %d (%.2f <= %.2f)', ~get_flag_value(rejection_flags, 'combined_leak', ci), combined_leak_val, max_combined_leak_thr), ...
    sprintf('hf_slope (not gated): %.2f (ref %.2f)', hf_slope_vec(ci), max_hf_slope_thr), ...
    sprintf('emg_score gate: %d (%.2f <= %.2f)', ~get_flag_value(rejection_flags, 'emg_score', ci), emg_score_vec(ci), max_emg_score_thr), ...
    sprintf('topo_nonposterior gate: %d', ~get_flag_value(rejection_flags, 'topo_nonposterior', ci)), ...
    sprintf('lineharm warn: %d (%.2f > %.2f)', lineharm_vec(ci) > max_lineharm_ratio_thr, lineharm_vec(ci), max_lineharm_ratio_thr), ...
    sprintf('front/temp leak: %.2f/%.2f', front_leak_vec(ci), temp_leak_vec(ci))};
info_viol = [ ...
    extreme_outlier_fail, ...
    get_flag_value(rejection_flags, 'combined_leak', ci), ...
    get_flag_value(rejection_flags, 'hf_slope', ci), ...
    get_flag_value(rejection_flags, 'emg_score', ci), ...
    get_flag_value(rejection_flags, 'topo_nonposterior', ci), ...
    false, ...
    false];
end

function plot_covariance_matrix_diagnostics(save_dir, subject_id, chan_labels, covBase_full, covStim_per_win, win_names_cap, lambdas)
if isempty(covStim_per_win) || isempty(covBase_full)
    return;
end
nWins = min(3, numel(covStim_per_win));
if nWins < 1
    return;
end
if isempty(win_names_cap)
    win_names_cap = {'Full', 'Early', 'Late'};
end
if isempty(lambdas)
    lambdas = repmat(0.05, 1, nWins);
end

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
nMatCols = 4;
nPlotCols = 5;
tiledlayout(nWins, nPlotCols, 'Padding', 'compact', 'TileSpacing', 'compact');

nChans = size(covBase_full, 1);
if size(covBase_full, 2) ~= nChans
    close(fig);
    return;
end

% Build matrices first so column-wise shared color limits are used.
all_mats = cell(nWins, nMatCols);
qc_lines_by_win = cell(nWins, 1);
for wi = 1:nWins
    covStim_w = covStim_per_win{wi};
    if isempty(covStim_w) || any(size(covStim_w) ~= [nChans nChans])
        covStim_w = nan(nChans);
    end
    lam_w = lambdas(min(wi, numel(lambdas)));

    d_stim = diag(covStim_w);
    d_stim = d_stim(isfinite(d_stim));
    stim_diag_mean = mean(d_stim);
    if ~isfinite(stim_diag_mean)
        stim_diag_mean = 1;
    end
    d_base = diag(covBase_full);
    d_base = d_base(isfinite(d_base));
    base_diag_mean = mean(d_base);
    if ~isfinite(base_diag_mean)
        base_diag_mean = 1;
    end

    covStim_reg = (1 - lam_w) * covStim_w + lam_w * stim_diag_mean * eye(nChans);
    covBase_reg = (1 - lam_w) * covBase_full + lam_w * base_diag_mean * eye(nChans);
    covStim_reg = real(0.5 * (covStim_reg + covStim_reg'));
    covBase_reg = real(0.5 * (covBase_reg + covBase_reg'));
    covDiff = covStim_reg - covBase_reg;

    rc_base = rcond(covBase_reg);
    if ~isfinite(rc_base) || rc_base < 1e-12
        gedOp = pinv(covBase_reg) * covStim_reg;
    else
        gedOp = covBase_reg \ covStim_reg;
    end
    gedOp = real(gedOp);

    all_mats{wi, 1} = covStim_reg;
    all_mats{wi, 2} = covBase_reg;
    all_mats{wi, 3} = covDiff;
    all_mats{wi, 4} = gedOp;

    sym_stim = norm(covStim_reg - covStim_reg', 'fro') / max(norm(covStim_reg, 'fro'), eps);
    sym_base = norm(covBase_reg - covBase_reg', 'fro') / max(norm(covBase_reg, 'fro'), eps);
    eig_base = eig(covBase_reg);
    eig_base = real(eig_base(isfinite(eig_base)));
    if isempty(eig_base)
        min_eig_base = NaN;
    else
        min_eig_base = min(eig_base);
    end
    cond_base = cond(covBase_reg);
    try
        eigvals_ged = eig(covStim_reg, covBase_reg);
    catch
        eigvals_ged = eig(gedOp);
    end
    eigvals_ged = real(eigvals_ged(isfinite(eigvals_ged)));
    if isempty(eigvals_ged)
        top_ged = nan(1, 5);
    else
        eigvals_ged = sort(eigvals_ged, 'descend');
        top_ged = nan(1, 5);
        nTop = min(5, numel(eigvals_ged));
        top_ged(1:nTop) = eigvals_ged(1:nTop);
    end
    qc_lines_by_win{wi} = sprintf(['QC:\n', ...
        'sym(S): %.2e\n', ...
        'sym(R): %.2e\n', ...
        'minEig(R): %.3g\n', ...
        'cond(R): %.2e\n', ...
        'GED top5:\n', ...
        '1) %.3f\n', ...
        '2) %.3f\n', ...
        '3) %.3f\n', ...
        '4) %.3f\n', ...
        '5) %.3f'], ...
        sym_stim, sym_base, min_eig_base, cond_base, ...
        top_ged(1), top_ged(2), top_ged(3), top_ged(4), top_ged(5));
end

col_clims = ones(1, nMatCols);
for pi = 1:nMatCols
    col_vals = [];
    for wi = 1:nWins
        mat_vals = abs(all_mats{wi, pi}(:));
        mat_vals = mat_vals(isfinite(mat_vals));
        col_vals = [col_vals; mat_vals];
    end
    if ~isempty(col_vals)
        clim = prctile(col_vals, 99);
        if ~isfinite(clim) || clim <= 0
            clim = max(col_vals);
        end
        if isfinite(clim) && clim > 0
            col_clims(pi) = clim;
        end
    end
end

% Blue-white-red diverging map for signed matrix interpretation.
nC = 256;
nHalf = nC / 2;
t = linspace(0, 1, nHalf)';
low_col = [0.10 0.25 0.80];
mid_col = [1.00 1.00 1.00];
high_col = [0.80 0.15 0.15];
cmap_low = [ ...
    low_col(1) + (mid_col(1) - low_col(1)) * t, ...
    low_col(2) + (mid_col(2) - low_col(2)) * t, ...
    low_col(3) + (mid_col(3) - low_col(3)) * t];
cmap_high = [ ...
    mid_col(1) + (high_col(1) - mid_col(1)) * t, ...
    mid_col(2) + (high_col(2) - mid_col(2)) * t, ...
    mid_col(3) + (high_col(3) - mid_col(3)) * t];
colormap(fig, [cmap_low; cmap_high]);

for wi = 1:nWins
    mats = all_mats(wi, :);
    panel_titles = { ...
        sprintf('%s Stim (reg)', win_names_cap{wi}), ...
        sprintf('%s Base (reg)', win_names_cap{wi}), ...
        sprintf('%s Stim-Base', win_names_cap{wi}), ...
        sprintf('%s R^{-1}S', win_names_cap{wi})};

    for pi = 1:nMatCols
        nexttile;
        mat = mats{pi};
        if all(~isfinite(mat(:)))
            mat = zeros(nChans);
        end
        imagesc(mat);
        axis image;
        caxis([-col_clims(pi) col_clims(pi)]);
        colorbar;
        title(panel_titles{pi}, 'FontSize', 10, 'Interpreter', 'none');
        set(gca, 'FontSize', 8, 'YDir', 'normal');
        if nChans <= 80 && ~isempty(chan_labels) && numel(chan_labels) == nChans
            xticks(1:nChans); yticks(1:nChans);
            xticklabels(chan_labels); yticklabels(chan_labels);
            xtickangle(90);
        else
            xlabel('Channels');
            ylabel('Channels');
        end
    end

    nexttile;
    axis off;
    title(sprintf('%s Values', win_names_cap{wi}), 'FontSize', 10, 'Interpreter', 'none');
    text(0.01, 0.99, qc_lines_by_win{wi}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 8, 'Color', [0.05 0.05 0.05], 'Interpreter', 'none');
end

sgtitle(sprintf('GED Covariance Diagnostics: Subject %s | Shared color scale per column', subject_id), ...
    'FontSize', 14, 'FontWeight', 'bold');
save_figure_png(fig, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_covariance_matrix.png', subject_id)));
close(fig);
end

function plot_rejection_info_text_columns(info_lines, info_viol)
y0 = 1.52;
dy = 0.11;
criteria_font_size = 3.5;
if isempty(info_lines)
    return;
end
n_lines = numel(info_lines);
n_left = ceil(n_lines / 2);
x_cols = [0.02, 0.52];
for li = 1:n_lines
    if li <= n_left
        col_idx = 1;
        row_idx = li;
    else
        col_idx = 2;
        row_idx = li - n_left;
    end
    if info_viol(li)
        txt_col = [0.82 0.10 0.10];
    else
        txt_col = [0.10 0.10 0.10];
    end
    text(x_cols(col_idx), y0 - (row_idx - 1) * dy, info_lines{li}, ...
        'Units', 'normalized', 'Clipping', 'off', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'FontSize', criteria_font_size, 'Color', txt_col, 'Interpreter', 'none');
end
end

function val = get_flag_value(flags_struct, field_name, idx)
val = false;
if isfield(flags_struct, field_name)
    field_val = flags_struct.(field_name);
    if numel(field_val) >= idx && isfinite(field_val(idx))
        val = logical(field_val(idx));
    end
end
end

function plot_combined_topo_spectra_windows(save_dir, subject_id, scan_freqs, cfg_topo, topo_labels, ...
    searchTopos_full, searchMeanPrSpectrum_full, selected_idx_full, w_combined_full, ...
    searchTopos_early, searchMeanPrSpectrum_early, selected_idx_early, w_combined_early, ...
    searchTopos_late, searchMeanPrSpectrum_late, selected_idx_late, w_combined_late, ...
    analysis_freq_range)

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
win_names = {'early', 'full', 'late'};
all_topos = {searchTopos_early, searchTopos_full, searchTopos_late};
all_specs = {searchMeanPrSpectrum_early, searchMeanPrSpectrum_full, searchMeanPrSpectrum_late};
all_idx = {selected_idx_early, selected_idx_full, selected_idx_late};
all_w = {w_combined_early, w_combined_full, w_combined_late};

for wi = 1:3
    subplot(2, 3, wi);
    topo_mat = all_topos{wi};
    spec_mat = all_specs{wi};
    sel_idx = all_idx{wi};
    sel_w = all_w{wi};
    [sel_idx, sel_w] = sanitize_selected_components(sel_idx, sel_w, size(spec_mat, 1));
    if isempty(sel_idx) || isempty(topo_mat)
        axis off;
        text(0.5, 0.5, sprintf('No combined components (%s)', upper(win_names{wi})), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 11, 'Color', [0.7 0.1 0.1], 'Interpreter', 'none');
    else
        topo_vec = topo_mat(:, sel_idx) * sel_w(:);
        topo_data = [];
        topo_data.label = topo_labels;
        topo_data.avg = topo_vec;
        topo_data.dimord = 'chan';
        topo_vals = topo_vec(isfinite(topo_vec));
        topo_clim = max(abs(topo_vals));
        if ~isfinite(topo_clim) || topo_clim <= 0
            topo_clim = 1;
        end
        cfg_ci = cfg_topo;
        cfg_ci.zlim = [-topo_clim topo_clim];
        try
            ft_topoplotER(cfg_ci, topo_data);
        catch
            imagesc(topo_vec(:)); axis tight;
            caxis([-topo_clim topo_clim]); colorbar;
        end
        title(sprintf('Combined Topography (%s, n=%d)', upper(win_names{wi}), numel(sel_idx)), ...
            'FontSize', 11, 'Interpreter', 'none');
    end
    set(gca, 'FontSize', 10);

    subplot(2, 3, wi + 3); hold on;
    if isempty(sel_idx) || isempty(spec_mat)
        axis off;
        text(0.5, 0.5, sprintf('No combined spectrum (%s)', upper(win_names{wi})), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 11, 'Color', [0.7 0.1 0.1], 'Interpreter', 'none');
        continue;
    end
    spec_vec = sel_w(:)' * spec_mat(sel_idx, :);
    [pf_score, pf_peak_hz] = compute_combined_powspctrm_form_metrics( ...
        spec_vec, scan_freqs, analysis_freq_range);
    plot(scan_freqs, spec_vec, '-', 'Color', [0 0 0], 'LineWidth', 2.0);
    add_peak_point_overlay(scan_freqs, spec_vec, pf_peak_hz);
    yline(0, 'k--', 'LineWidth', 0.8);
    xlim([analysis_freq_range(1) analysis_freq_range(2)]);
    spec_finite = spec_vec(isfinite(spec_vec));
    if ~isempty(spec_finite)
        sp_min = min(spec_finite);
        sp_max = max(spec_finite);
        if isfinite(sp_min) && isfinite(sp_max) && sp_min < sp_max
            sp_range = sp_max - sp_min;
            ylim([sp_min - 0.12 * sp_range, sp_max + 0.20 * sp_range]);
        end
    end
    format_power_change_db_axis(gca);
    xlabel('Hz'); ylabel('Power [dB]');
    box on;
    text(0.02, 0.98, sprintf('PF = %.2f | peak = %.1f Hz', ...
        pf_score, pf_peak_hz), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 9, 'Interpreter', 'none', 'Color', [0.1 0.1 0.1]);
    title(sprintf('Combined Spectrum (%s)', upper(win_names{wi})), ...
        'FontSize', 11, 'Interpreter', 'none');
    set(gca, 'FontSize', 10);
end

sgtitle(sprintf('Combined GED Components: %s', subject_id), ...
    'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
save_figure_png(fig, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_topo_spectra_combined.png', subject_id)));
close(fig);
end

function add_peak_point_overlay(freq_axis, spec_vec, peak_hz)
if isempty(freq_axis) || isempty(spec_vec) || numel(freq_axis) ~= numel(spec_vec)
    return;
end
valid = isfinite(freq_axis) & isfinite(spec_vec);
x = freq_axis(valid);
y = spec_vec(valid);
if numel(x) < 2 || numel(y) < 2
    return;
end
[x, sort_idx] = sort(x(:)');
y = y(sort_idx);
baseline_y = min(y);
if ~isfinite(baseline_y)
    return;
end

if ~isfinite(peak_hz)
    [~, idx_peak] = max(y);
    peak_hz = x(idx_peak);
end

peak_power = interp1(x, y, peak_hz, 'linear', NaN);
if ~isfinite(peak_power)
    [~, idx_near] = min(abs(x - peak_hz));
    peak_hz = x(idx_near);
    peak_power = y(idx_near);
end
if isfinite(peak_hz) && isfinite(peak_power)
    plot([peak_hz peak_hz], [baseline_y peak_power], '--', ...
        'Color', [0.55 0.55 0.55], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    plot(peak_hz, peak_power, 'o', ...
        'Color', [0.55 0.55 0.55], 'MarkerFaceColor', [0.55 0.55 0.55], ...
        'MarkerSize', 4, 'HandleVisibility', 'off');
end
end

function [sel_idx, sel_w] = sanitize_selected_components(sel_idx, sel_w, nComp)
if isempty(nComp) || ~isfinite(nComp)
    nComp = 0;
end
sel_idx = sel_idx(:);
sel_idx = sel_idx(isfinite(sel_idx));
sel_idx = round(sel_idx);
sel_idx = sel_idx(sel_idx >= 1 & sel_idx <= nComp);
if isempty(sel_idx)
    sel_w = [];
    return;
end
if isempty(sel_w) || numel(sel_w) ~= numel(sel_idx)
    sel_w = ones(numel(sel_idx), 1);
else
    sel_w = sel_w(:);
end
sel_w(~isfinite(sel_w) | sel_w <= 0) = 0;
if sum(sel_w) <= 0
    sel_w = ones(numel(sel_idx), 1);
end
sel_w = sel_w / sum(sel_w);
end

function [pf_score, pf_peak_hz] = compute_combined_powspctrm_form_metrics( ...
    spec_vec, scan_freqs, analysis_freq_range)
pf_score = NaN;
pf_peak_hz = NaN;
if isempty(spec_vec) || isempty(scan_freqs) || numel(spec_vec) ~= numel(scan_freqs)
    return;
end
[pf_vec, ~] = compute_powspctrm_form_laplacian_score_from_spectra( ...
    spec_vec(:)', scan_freqs, analysis_freq_range);
if ~isempty(pf_vec) && isfinite(pf_vec(1))
    pf_score = pf_vec(1);
end
freq_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
if any(freq_mask)
    x = scan_freqs(freq_mask);
    y = spec_vec(freq_mask);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    if numel(y) >= 3
        peak_lo = max(30, analysis_freq_range(1));
        peak_hi = min(90, analysis_freq_range(2));
        inner_mask = x >= peak_lo & x <= peak_hi;
        if sum(inner_mask) >= 3
            x_peak = x(inner_mask);
            y_peak = y(inner_mask);
        else
            x_peak = x;
            y_peak = y;
        end
        [~, idx_max] = max(y_peak);
        pf_peak_hz = x_peak(idx_max);
    end
end
if ~isfinite(pf_peak_hz)
    if sum(freq_mask) >= 3
        if numel(y) >= 3
            peak_lo = max(30, analysis_freq_range(1));
            peak_hi = min(90, analysis_freq_range(2));
            inner_mask = x >= peak_lo & x <= peak_hi;
            if sum(inner_mask) >= 3
                x_peak = x(inner_mask);
                y_peak = y(inner_mask);
            else
                x_peak = x;
                y_peak = y;
            end
            [~, idx_max] = max(y_peak);
            pf_peak_hz = x_peak(idx_max);
        end
    end
end
end

function [pf_score_vec, pf_deduction_vec, diag] = compute_powspctrm_form_laplacian_score_from_spectra( ...
    mean_pr_spectrum, scan_freqs, analysis_freq_range)
rival_sep_min_hz = 10;
rival_height_ratio_min = 0.90;
deduct_per_rival = 0.05;
deduct_per_rival = max(0, min(1, deduct_per_rival));
nComp = size(mean_pr_spectrum, 1);
pf_score_vec = zeros(nComp, 1);
pf_deduction_vec = zeros(nComp, 1);
diag = struct( ...
    'laplace_r2', zeros(nComp, 1), ...
    'dominant_peak_hz', nan(nComp, 1), ...
    'n_rival_peaks', zeros(nComp, 1));
if isempty(mean_pr_spectrum) || isempty(scan_freqs)
    return;
end
freq_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
if ~any(freq_mask)
    return;
end
for ci = 1:nComp
    y = mean_pr_spectrum(ci, :);
    x = scan_freqs(freq_mask);
    y = y(freq_mask);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    if numel(y) < 7
        continue;
    end
    [x, sidx] = sort(x(:)');
    y = y(sidx);
    y_smooth = movmean(y, 5, 'omitnan');
    y_floor = prctile(y_smooth, 20);
    if ~isfinite(y_floor)
        y_floor = median(y_smooth, 'omitnan');
    end
    if ~isfinite(y_floor)
        y_floor = 0;
    end
    y_pos = max(y_smooth - y_floor, 0);
    y_scale = max(y_pos);
    if ~isfinite(y_scale) || y_scale <= eps
        continue;
    end
    y_norm = y_pos / y_scale;
    [dom_amp, dom_idx] = max(y_pos);
    if ~isfinite(dom_amp) || dom_amp <= eps || isempty(dom_idx)
        continue;
    end
    dom_hz = x(dom_idx);
    diag.dominant_peak_hz(ci) = dom_hz;

    p0 = [1.0, dom_hz, 6.0];
    lb = [0.20, min(x), 1.0];
    ub = [1.50, max(x), 20.0];
    obj = @(p) mean((y_norm - laplace_template_model(x, p)).^2, 'omitnan');
    opts = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 3000);
    p_est = fminsearch(@(pp) bounded_obj(pp, lb, ub, obj), p0, opts);
    p_est = min(max(p_est, lb), ub);
    y_hat = laplace_template_model(x, p_est);
    r2 = compute_r2_score(y_norm, y_hat);
    diag.laplace_r2(ci) = r2;

    robust_scale = robust_mad(y_pos);
    if ~isfinite(robust_scale) || robust_scale <= eps
        robust_scale = iqr(y_pos);
    end
    if ~isfinite(robust_scale) || robust_scale <= eps
        robust_scale = std(y_pos(isfinite(y_pos)));
    end
    if ~isfinite(robust_scale) || robust_scale <= eps
        robust_scale = 1;
    end
    rel_prom = 0.08 * max(y_pos);
    min_prom = max([0, rel_prom, 0.02, 0.15 * robust_scale]);
    [pks, locs] = findpeaks(y_pos, x, 'MinPeakProminence', min_prom, 'MinPeakDistance', 5);
    n_rivals = 0;
    if ~isempty(pks) && numel(pks) >= 2
        pks = pks(:);
        locs = locs(:);
        dom_mask = abs(locs - dom_hz) <= 1e-9;
        if ~any(dom_mask)
            [~, nearest_idx] = min(abs(locs - dom_hz));
            dom_mask = false(size(locs));
            dom_mask(nearest_idx) = true;
        end
        rival_mask = ~dom_mask & ...
            (abs(locs - dom_hz) >= rival_sep_min_hz) & ...
            (pks >= rival_height_ratio_min * dom_amp);
        n_rivals = sum(rival_mask);
    end
    diag.n_rival_peaks(ci) = n_rivals;
    deduct_val = min(r2, n_rivals * deduct_per_rival);
    pf_deduction_vec(ci) = deduct_val;
    pf_score_vec(ci) = max(0, min(1, r2 - deduct_val));
end
pf_score_vec(~isfinite(pf_score_vec)) = 0;
pf_deduction_vec(~isfinite(pf_deduction_vec)) = 0;
end

function yhat = laplace_template_model(x, p)
A = p(1);
mu = p(2);
b = p(3);
yhat = A .* exp(-abs(x - mu) ./ max(b, eps));
yhat = max(yhat, 0);
if max(yhat) > 0
    yhat = yhat ./ max(yhat);
end
end

function r2 = compute_r2_score(y, yhat)
if numel(y) ~= numel(yhat)
    r2 = 0;
    return;
end
valid = isfinite(y) & isfinite(yhat);
y = y(valid);
yhat = yhat(valid);
if numel(y) < 3
    r2 = 0;
    return;
end
sse = sum((y - yhat).^2);
sst = sum((y - mean(y)).^2);
if ~isfinite(sst) || sst <= eps
    r2 = 0;
else
    r2 = 1 - (sse / sst);
end
r2 = max(0, min(1, r2));
end

function v = bounded_obj(p, lb, ub, obj_fun)
p_clip = min(max(p, lb), ub);
v = obj_fun(p_clip);
if ~isfinite(v)
    v = 1e6;
end
end

function format_power_change_db_axis(axh)
if isempty(axh) || ~isgraphics(axh, 'axes')
    axh = gca;
end
yt = get(axh, 'YTick');
if isempty(yt)
    return;
end
yt_lbl = arrayfun(@(v) sprintf('%g', v), yt, 'UniformOutput', false);
set(axh, 'YTickLabel', yt_lbl);
end

function [peak_bonus_vec, peak_count_vec] = compute_peak_bonus_from_spectra( ...
    mean_pr_spectrum, scan_freqs, analysis_freq_range)
nComp = size(mean_pr_spectrum, 1);
peak_bonus_vec = zeros(nComp, 1);
peak_count_vec = zeros(nComp, 1);
peak_prom_abs_floor = 0.02;
if isempty(mean_pr_spectrum) || isempty(scan_freqs)
    return;
end
freq_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
for ci = 1:nComp
    y = mean_pr_spectrum(ci, :);
    if all(~isfinite(y))
        continue;
    end
    y_band = y(freq_mask);
    x_band = scan_freqs(freq_mask);
    valid = isfinite(y_band) & isfinite(x_band);
    y_band = y_band(valid);
    x_band = x_band(valid);
    if numel(y_band) < 5
        continue;
    end
    y_band = movmean(y_band, 3);
    y_shape = max(y_band - median(y_band), 0);
    peak_scale = max(y_shape);
    if ~isfinite(peak_scale) || peak_scale <= eps
        continue;
    end
    robust_scale = robust_mad(y_shape);
    if ~isfinite(robust_scale) || robust_scale <= eps
        robust_scale = iqr(y_shape);
    end
    if ~isfinite(robust_scale) || robust_scale <= eps
        robust_scale = std(y_shape(isfinite(y_shape)));
    end
    if ~isfinite(robust_scale) || robust_scale <= eps
        robust_scale = 1;
    end
    min_prom = max([0, 0.10 * peak_scale, peak_prom_abs_floor, 0.15 * robust_scale]);
    [pks, ~] = findpeaks(y_shape, x_band, 'MinPeakProminence', min_prom, 'MinPeakDistance', 5);
    if isempty(pks)
        continue;
    end
    pks = sort(pks(:), 'descend');
    n_keep = min(2, numel(pks));
    peak_count_vec(ci) = n_keep;
    amp_score = min(1, mean(pks(1:n_keep)) / max(peak_scale, eps));
    count_score = n_keep / 2;
    peak_bonus_vec(ci) = max(0, min(1, 0.65 * amp_score + 0.35 * count_score));
end
end

function [outlier_mask, stats] = detect_trial_metric_outliers_iqr(x, iqr_mult)
x = x(:);
outlier_mask = false(size(x));
stats = struct('n_finite', 0, 'q1', NaN, 'q3', NaN, 'iqr_val', NaN, 'lo', NaN, 'hi', NaN, 'n_outliers', 0);
valid = isfinite(x);
n_finite = sum(valid);
stats.n_finite = n_finite;
if n_finite == 0
    return;
end
q = quantile(x(valid), [0.25 0.75]);
q1 = q(1);
q3 = q(2);
iqr_val = q3 - q1;
stats.q1 = q1;
stats.q3 = q3;
stats.iqr_val = iqr_val;
if ~isfinite(q1) || ~isfinite(q3) || ~isfinite(iqr_val) || iqr_val <= eps
    return;
end
lo = q1 - iqr_mult * iqr_val;
hi = q3 + iqr_mult * iqr_val;
stats.lo = lo;
stats.hi = hi;
outlier_mask = valid & (x < lo | x > hi);
stats.n_outliers = sum(outlier_mask);
end

function m = robust_trial_mean(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    m = NaN;
    return;
end
if numel(x) < 4
    m = mean(x);
    return;
end
q = quantile(x, [0.25 0.75]);
iqr_val = q(2) - q(1);
if ~isfinite(iqr_val) || iqr_val <= eps
    m = median(x);
    return;
end
lo_thr = q(1) - 1.5 * iqr_val;
hi_thr = q(2) + 1.5 * iqr_val;
x_trim = x(x >= lo_thr & x <= hi_thr);
if isempty(x_trim)
    m = median(x);
else
    m = mean(x_trim);
end
end

function m = robust_mad(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    m = 0;
    return;
end
m = mad(x, 1);
if ~isfinite(m) || m <= eps
    m = std(x);
end
if ~isfinite(m) || m <= eps
    m = 0;
end
end

function proxy = estimate_component_artifact_proxies(filter_w, dat_per_cond, stim_window, base_window, fs, scan_freqs, tapsmofrq_hz)
proxy = struct( ...
    'lineharm_ratio', NaN, ...
    'stationarity_cv', NaN, ...
    'burst_ratio', NaN, ...
    'hf_slope', NaN, ...
    'mean_pr_spectrum', nan(1, numel(scan_freqs)), ...
    'n_trials_used', 0);
if isempty(filter_w) || isempty(dat_per_cond)
    return;
end

band_mask = scan_freqs >= 30 & scan_freqs <= 90;
hf_mask = scan_freqs >= 70 & scan_freqs <= min(110, max(scan_freqs));
harm_mask = ismember(round(scan_freqs), [50 60 100]);
if ~any(harm_mask)
    [~, idx50] = min(abs(scan_freqs - 50));
    harm_mask(idx50) = true;
end

trial_gamma = [];
trial_pr = [];
lineharm_acc = [];
trial_count = 0;
stim_sig_cells = {};
base_sig_cells = {};
for cond = 1:numel(dat_per_cond)
    dat = dat_per_cond{cond};
    if isempty(dat) || ~isfield(dat, 'trial')
        continue;
    end
    for trl = 1:numel(dat.trial)
        x = double(dat.trial{trl});
        t = dat.time{trl};
        if isempty(x) || isempty(t)
            continue;
        end
        x_proj = filter_w' * x;
        idx_base = t >= base_window(1) & t <= base_window(2);
        idx_stim = t >= stim_window(1) & t <= stim_window(2);
        if sum(idx_base) < 5 || sum(idx_stim) < 5
            continue;
        end
        base_sig_cells{end+1, 1} = x_proj(idx_base);
        stim_sig_cells{end+1, 1} = x_proj(idx_stim);
        trial_count = trial_count + 1;
    end
end
proxy.n_trials_used = trial_count;
if trial_count < 4
    return;
end

if ~isfinite(fs) || fs <= 0
    return;
end
trial_pr = compute_simple_power_ratio_scan_batch(stim_sig_cells, base_sig_cells, fs, scan_freqs, tapsmofrq_hz);
if isempty(trial_pr)
    return;
end
for ri = 1:size(trial_pr, 1)
    pr_full = trial_pr(ri, :);
    if all(~isfinite(pr_full)) || ~any(band_mask)
        continue;
    end
    pr_band = pr_full(band_mask);
    if all(~isfinite(pr_band))
        continue;
    end
    nonharm_val = nanmean(pr_full(band_mask & ~harm_mask));
    harm_val = nanmean(pr_full(band_mask & harm_mask));
    if ~isfinite(nonharm_val) || abs(nonharm_val) <= eps || ~isfinite(harm_val)
        continue;
    end
    trial_gamma(end+1, 1) = nanmean(pr_band);
    lineharm_acc(end+1, 1) = max(harm_val, 0) / max(abs(nonharm_val), eps);
end
if isempty(trial_gamma) || isempty(lineharm_acc)
    return;
end
proxy.mean_pr_spectrum = nanmean(trial_pr, 1);
proxy.lineharm_ratio = nanmean(lineharm_acc);
proxy.stationarity_cv = nanstd(trial_gamma) / max(abs(nanmean(trial_gamma)), eps);
burst_thr = prctile(abs(trial_gamma), 85);
if ~isfinite(burst_thr) || burst_thr <= 0
    burst_thr = median(abs(trial_gamma)) + mad(abs(trial_gamma), 1);
end
proxy.burst_ratio = mean(abs(trial_gamma) >= max(burst_thr, eps));

if sum(hf_mask) >= 3
    f_hf = scan_freqs(hf_mask);
    spec_hf = proxy.mean_pr_spectrum(hf_mask);
    xh_full = log(f_hf(:));
    yh_full = log(abs(spec_hf(:)) + eps);
    valid_hf = isfinite(xh_full) & isfinite(yh_full) & ...
        isfinite(spec_hf(:)) & (spec_hf(:) > -inf);
    n_hf = nnz(valid_hf);
    if n_hf >= 3
        xh = xh_full(valid_hf);
        yh = yh_full(valid_hf);
        p = polyfit(xh, yh, 1);
        if isfinite(p(1))
            proxy.hf_slope = p(1);
        end
    end
end

end

function bad_mask = flag_unreliable_baseline_trials(base_power_vals, mad_mult)
bad_mask = false(size(base_power_vals));
if isempty(base_power_vals)
    return;
end
valid = isfinite(base_power_vals) & (base_power_vals > 0);
if ~any(valid)
    return;
end
logp = log10(base_power_vals(valid));
center = median(logp, 'omitnan');
spread = robust_mad(logp);
if ~isfinite(spread) || spread <= eps
    spread = iqr(logp) / 1.349;
end
if ~isfinite(spread) || spread <= eps
    return;
end
zabs = abs((logp - center) / spread);
tmp = false(size(logp));
tmp(zabs > mad_mult) = true;
bad_mask(valid) = tmp;
end

function pr_mat = compute_simple_power_ratio_scan_batch(sig_stim_cells, sig_base_cells, fs, scan_freqs, tapsmofrq_hz)
% Component-selection proxy retained unchanged intentionally.
nRows = numel(sig_stim_cells);
pr_mat = nan(nRows, numel(scan_freqs));
if nRows == 0 || numel(sig_base_cells) ~= nRows || fs <= 0
    return;
end
[p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(sig_stim_cells, sig_base_cells, fs, scan_freqs, tapsmofrq_hz);
if isempty(p_stim_scan) || isempty(p_base_scan)
    return;
end
for ri = 1:nRows
    p_stim_row = p_stim_scan(ri, :);
    p_base_row = p_base_scan(ri, :);
    valid_base = isfinite(p_base_row) & (p_base_row > 0);
    if ~any(valid_base)
        continue;
    end
    base_anchor = prctile(p_base_row(valid_base), 20);
    base_median = median(p_base_row(valid_base), 'omitnan');
    if ~isfinite(base_anchor) || base_anchor <= 0
        base_anchor = base_median;
    end
    if ~isfinite(base_median) || base_median <= 0
        base_median = base_anchor;
    end
    base_floor = max(0.25 * base_anchor, eps);
    valid = isfinite(p_stim_row) & isfinite(p_base_row) & (p_base_row > 0);
    ratio_row = nan(1, numel(scan_freqs));
    ratio_row(valid) = 10 * log10((p_stim_row(valid) + base_floor) ./ (p_base_row(valid) + base_floor));
    pr_mat(ri, :) = ratio_row;
end
end

function p_scan = compute_scan_power_mtmfft_ft(sig, fs, scan_freqs, tapsmofrq_hz)
if ~iscell(sig) && isvector(sig)
    sig = sig(:)';
end
if iscell(sig)
    nSig = numel(sig);
else
    nSig = size(sig, 1);
end
p_scan = nan(nSig, numel(scan_freqs));
if isempty(sig) || fs <= 0 || isempty(scan_freqs) || ~isfinite(tapsmofrq_hz) || tapsmofrq_hz <= 0
    return;
end
if ~iscell(sig) && size(sig, 2) < 8
    return;
end

valid_rows = false(nSig, 1);
dat = [];
dat.label = {'GED'};
dat.fsample = fs;
dat.trial = {};
dat.time = {};
for si = 1:nSig
    if iscell(sig)
        x = double(sig{si});
    else
        x = double(sig(si, :));
    end
    if iscolumn(x)
        x = x';
    end
    if any(~isfinite(x))
        continue;
    end
    if numel(x) < 8
        continue;
    end
    x = x - mean(x);
    valid_rows(si) = true;
    dat.trial{end+1} = x;
    dat.time{end+1} = (0:(numel(x)-1)) / fs;
end
if ~any(valid_rows)
    return;
end

% Provide explicit sampleinfo to avoid FieldTrip reconstruction warnings.
nTrials_valid = numel(dat.trial);
dat.sampleinfo = zeros(nTrials_valid, 2);
sample_start = 1;
for ti = 1:nTrials_valid
    nSamp = numel(dat.trial{ti});
    dat.sampleinfo(ti, :) = [sample_start, sample_start + nSamp - 1];
    sample_start = sample_start + nSamp;
end

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'dpss';
cfg.foi = scan_freqs;
cfg.tapsmofrq = tapsmofrq_hz;
cfg.pad = 'nextpow2';
cfg.keeptrials = 'yes';
cfg.feedback = 'none';
progress_state = ged_freq_progress_next_call();
clc;
fprintf('[GED] Subject %s (%d/%d) %s %d/%d\n', ...
    progress_state.subject, progress_state.subj_idx, progress_state.subj_total, ...
    progress_state.phase, progress_state.call_idx, progress_state.call_total);
try
    freq = ft_freqanalysis(cfg, dat);
catch
    return;
end

pow = freq.powspctrm;
if ndims(pow) == 3
    pow = squeeze(pow(:, 1, :));
elseif isvector(pow)
    pow = pow(:)';
end
if isfield(freq, 'freq') && ~isempty(freq.freq)
    n_freq = numel(freq.freq);
else
    n_freq = numel(scan_freqs);
end
if isvector(pow) && numel(pow) == n_freq
    % Keep one-trial output as row [1 x nFreq].
    pow = reshape(pow, 1, []);
elseif size(pow, 1) == n_freq && size(pow, 2) ~= n_freq
    % Guard against squeeze() producing [nFreq x nTrials].
    pow = pow.';
end
if size(pow, 1) == 1 && sum(valid_rows) > 1
    pow = repmat(pow, sum(valid_rows), 1);
end

if size(pow, 2) ~= numel(scan_freqs) && isfield(freq, 'freq') && ~isempty(freq.freq)
    freq_axis = freq.freq(:)';
    pow_interp = nan(size(pow, 1), numel(scan_freqs));
    for ri = 1:size(pow, 1)
        row_pow = pow(ri, :);
        if numel(row_pow) ~= numel(freq_axis) && size(pow, 1) == numel(freq_axis) && size(pow, 2) == size(pow_interp, 1)
            row_pow = pow(:, ri).';
        end
        if numel(row_pow) ~= numel(freq_axis)
            continue;
        end
        pow_interp(ri, :) = interp1(freq_axis, row_pow, scan_freqs, 'linear', NaN);
    end
    pow = pow_interp;
end

p_scan(valid_rows, :) = pow;
end

function [p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(sig_stim, sig_base, fs, scan_freqs, tapsmofrq_hz)
if ~iscell(sig_stim) && isvector(sig_stim)
    sig_stim = sig_stim(:)';
end
if ~iscell(sig_base) && isvector(sig_base)
    sig_base = sig_base(:)';
end
nFreq = numel(scan_freqs);
if iscell(sig_stim)
    nStim = numel(sig_stim);
else
    nStim = size(sig_stim, 1);
end
p_stim_scan = nan(nStim, nFreq);
p_base_scan = nan(nStim, nFreq);
if isempty(sig_stim) || isempty(sig_base) || fs <= 0
    return;
end
if iscell(sig_stim) ~= iscell(sig_base)
    return;
end
if iscell(sig_stim)
    if numel(sig_stim) ~= numel(sig_base)
        return;
    end
    sig_all = [sig_stim(:); sig_base(:)];
else
    if size(sig_stim, 1) ~= size(sig_base, 1)
        return;
    end
    sig_all = [sig_stim; sig_base];
end
if isempty(sig_all)
    return;
end
p_all = compute_scan_power_mtmfft_ft(sig_all, fs, scan_freqs, tapsmofrq_hz);
if isempty(p_all) || size(p_all, 1) < 2 * nStim
    return;
end
p_stim_scan = p_all(1:nStim, :);
p_base_scan = p_all((nStim + 1):(2 * nStim), :);
end

function w_combined = build_combined_filter_vector(W_combined, w_components)
w_combined = [];
if isempty(W_combined)
    return;
end
nComp = size(W_combined, 2);
if nComp < 1
    return;
end
if isempty(w_components) || numel(w_components) ~= nComp
    w_components = ones(nComp, 1);
else
    w_components = w_components(:);
end
w_components(~isfinite(w_components) | w_components <= 0) = 0;
if sum(w_components) <= 0
    w_components = ones(nComp, 1);
end
w_components = w_components / sum(w_components);
w_combined = W_combined * w_components;
if ~all(isfinite(w_combined))
    w_combined = [];
end
end

function [ratio_trials, near_floor_freq_count_per_trial, valid_freq_count_per_trial] = ...
    compute_scan_ratio_for_combined_filter_batch(trial_cache, filter_vec, stim_field, trial_mask, fs, scan_freqs, tapsmofrq_hz, base_floor, near_floor_mult)
nTrl = numel(trial_cache);
ratio_trials = nan(nTrl, numel(scan_freqs));
near_floor_freq_count_per_trial = zeros(nTrl, 1);
valid_freq_count_per_trial = zeros(nTrl, 1);
if nTrl == 0 || isempty(filter_vec) || ~any(trial_mask)
    return;
end
filter_vec = filter_vec(:);
if ~all(isfinite(filter_vec))
    return;
end

sig_stim_cells = cell(0, 1);
sig_base_cells = cell(0, 1);
row_trial_idx = zeros(0, 1);
for trl = 1:nTrl
    if ~trial_mask(trl)
        continue;
    end
    tc = trial_cache{trl};
    x_base = tc.x_base;
    x_stim = tc.(stim_field);
    if isempty(x_base) || isempty(x_stim)
        continue;
    end
    sig_base_cells{end+1, 1} = (filter_vec' * x_base);
    sig_stim_cells{end+1, 1} = (filter_vec' * x_stim);
    row_trial_idx(end+1, 1) = trl;
end
if isempty(sig_stim_cells)
    return;
end

[ratio_rows, ~, near_floor_row_mask] = compute_scan_ratio_from_timeseries( ...
    sig_stim_cells, sig_base_cells, fs, scan_freqs, tapsmofrq_hz, base_floor, near_floor_mult);
for ri = 1:size(ratio_rows, 1)
    trl = row_trial_idx(ri);
    ratio_trials(trl, :) = ratio_rows(ri, :);
    valid_freq_count_per_trial(trl) = sum(isfinite(ratio_rows(ri, :)));
    near_floor_freq_count_per_trial(trl) = sum(near_floor_row_mask(ri, :));
end
end

function ged_freq_progress_reset(subject_label, subj_idx, subj_total, call_total_estimate)
global GED_FREQ_PROGRESS;
if ~isfinite(call_total_estimate) || call_total_estimate < 1
    call_total_estimate = 0;
end
GED_FREQ_PROGRESS = struct( ...
    'subject', subject_label, ...
    'subj_idx', subj_idx, ...
    'subj_total', subj_total, ...
    'call_idx', 0, ...
    'call_total', round(call_total_estimate), ...
    'phase', 'Spectral Batch');
end

function ged_freq_progress_set_phase(phase_label)
global GED_FREQ_PROGRESS;
if isempty(phase_label)
    return;
end
if isempty(GED_FREQ_PROGRESS) || ~isstruct(GED_FREQ_PROGRESS)
    GED_FREQ_PROGRESS = struct( ...
        'subject', 'NA', ...
        'subj_idx', 0, ...
        'subj_total', 0, ...
        'call_idx', 0, ...
        'call_total', 0, ...
        'phase', 'Spectral Batch');
end
GED_FREQ_PROGRESS.phase = phase_label;
end

function state = ged_freq_progress_next_call()
global GED_FREQ_PROGRESS;
if isempty(GED_FREQ_PROGRESS) || ~isstruct(GED_FREQ_PROGRESS)
    GED_FREQ_PROGRESS = struct( ...
        'subject', 'NA', ...
        'subj_idx', 0, ...
        'subj_total', 0, ...
        'call_idx', 0, ...
        'call_total', 0, ...
        'phase', 'Spectral Batch');
end
GED_FREQ_PROGRESS.call_idx = GED_FREQ_PROGRESS.call_idx + 1;
if GED_FREQ_PROGRESS.call_total < GED_FREQ_PROGRESS.call_idx
    GED_FREQ_PROGRESS.call_total = GED_FREQ_PROGRESS.call_idx;
end
state = GED_FREQ_PROGRESS;
end

function Wn = normalize_filters_to_noise_metric(W, covBase)
Wn = W;
if isempty(W) || isempty(covBase)
    return;
end
for ci = 1:size(W, 2)
    w = W(:, ci);
    if ~all(isfinite(w))
        continue;
    end
    denom = real(w' * covBase * w);
    if ~isfinite(denom) || denom <= eps
        continue;
    end
    Wn(:, ci) = w / sqrt(denom);
end
end

function [eligible_mask_out, extreme_component_outlier_idx] = exclude_extreme_component_outlier(eigvals_sorted, eligible_mask_in, ratio_thr, mad_mult)
eligible_mask_out = logical(eligible_mask_in(:));
extreme_component_outlier_idx = [];
if isempty(eigvals_sorted) || isempty(eligible_mask_out)
    return;
end

eigvals_sorted = eigvals_sorted(:);
if numel(eigvals_sorted) ~= numel(eligible_mask_out)
    return;
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
if isempty(rest_vals) || ~isfinite(lambda_top) || ~isfinite(lambda_second) || lambda_second <= 0
    return;
end

rest_log = log(rest_vals);
rest_log = rest_log(isfinite(rest_log));
if isempty(rest_log)
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
    extreme_component_outlier_idx = top_idx;
end
end

function save_figure_png(fig_handle, out_path)
if isempty(fig_handle) || ~ishandle(fig_handle)
    return;
end
if isempty(out_path)
    return;
end
out_dir = fileparts(out_path);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
drawnow;
pause(0.05);
set(fig_handle, 'PaperPositionMode', 'auto');
print(fig_handle, out_path, '-dpng', '-r300');
end

function freq = ged_powcurve_to_freq_ft(pow_vec, scan_freqs, data_recording)
% One-row FieldTrip freq struct (dimord chan_freq) for a subject-level GED spectrum in dB.
% Optional data_recording (e.g. preprocessed timelock/raw): attach elec/grad only.
% Do not set topolabel: ft_datatype treats any struct with topolabel as comp, so
% ft_checkdata then calls ft_datatype_comp and expects topo (see FieldTrip ft_datatype.m).
freq = struct();
freq.label = {'GED'};
freq.freq = scan_freqs(:)';
freq.dimord = 'chan_freq';
nf = numel(scan_freqs);
if isempty(scan_freqs)
    freq.powspctrm = nan(1, 0);
    return;
end
if isempty(pow_vec) || ~isnumeric(pow_vec) || numel(pow_vec) ~= nf
    freq.powspctrm = nan(1, nf);
else
    freq.powspctrm = reshape(double(pow_vec(:)).', 1, nf);
end
if nargin >= 3 && isstruct(data_recording)
    if isfield(data_recording, 'elec')
        el = data_recording.elec;
        if isstruct(el) && isfield(el, 'label') && ~isempty(el.label)
            freq.elec = el;
        end
    end
    if isfield(data_recording, 'grad')
        g = data_recording.grad;
        if isstruct(g) && isfield(g, 'label') && ~isempty(g.label)
            freq.grad = g;
        end
    end
end
end

function avg_curve = compute_condition_average_powratio_ft(pr_mat, scan_freqs)
avg_curve = nan(1, numel(scan_freqs));
if isempty(pr_mat) || isempty(scan_freqs)
    return;
end
if size(pr_mat, 2) ~= numel(scan_freqs)
    return;
end
valid_trials = any(isfinite(pr_mat), 2);
if ~any(valid_trials)
    return;
end
freq_dat = [];
freq_dat.label = {'GED'};
freq_dat.freq = scan_freqs(:)';
freq_dat.dimord = 'rpt_chan_freq';
freq_dat.powspctrm = nan(sum(valid_trials), 1, numel(scan_freqs));
trial_rows = find(valid_trials);
for ti = 1:numel(trial_rows)
    freq_dat.powspctrm(ti, 1, :) = pr_mat(trial_rows(ti), :);
end
cfg = [];
cfg.avgoverrpt = 'yes';
try
    freq_avg = ft_selectdata(cfg, freq_dat);
catch
    avg_curve = mean(pr_mat(valid_trials, :), 1, 'omitnan');
    return;
end
pow_avg = freq_avg.powspctrm;
if ndims(pow_avg) == 3
    pow_avg = squeeze(pow_avg(1, 1, :));
elseif ismatrix(pow_avg)
    pow_avg = squeeze(pow_avg);
end
avg_curve = pow_avg(:)';
if numel(avg_curve) ~= numel(scan_freqs)
    avg_curve = mean(pr_mat(valid_trials, :), 1, 'omitnan');
end
if all(~isfinite(avg_curve))
    avg_curve = mean(pr_mat(valid_trials, :), 1, 'omitnan');
end
end

function runtime_str = format_runtime_hhmmss(runtime_seconds)
if ~isfinite(runtime_seconds) || runtime_seconds < 0
    runtime_str = 'n/a';
    return;
end
runtime_seconds = round(runtime_seconds);
h = floor(runtime_seconds / 3600);
m = floor(mod(runtime_seconds, 3600) / 60);
s = mod(runtime_seconds, 60);
runtime_str = sprintf('%02d:%02d:%02d', h, m, s);
end
