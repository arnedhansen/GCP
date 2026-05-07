%% GCP Trial-Level Gamma Peak Frequency via Generalized Eigendecomposition (GED)
%
% Extracts trial-level gamma-band (30-90 Hz) peak frequencies from EEG data
% recorded during a visual grating task with four contrast conditions
% (25%, 50%, 75%, 100%). Generalized eigendecomposition is used to derive
% spatial filters that maximise stimulus-vs-baseline gamma power,
% replacing conventional channel-level analyses with a component-space
% approach that improves the signal-to-noise ratio of narrow-band
% gamma oscillations.
%
% Pipeline overview
% -----------------
% Broadband GED and component selection (per subject)
%   1. Pool trials across conditions and compute gamma-band covariances
%      (30-90 Hz FIR) for three stimulus windows (full/early/late) and
%      one shared pre-stimulus baseline.
%   2. Solve window-specific GED (S_stim * w = lambda * S_base * w) with
%      regularisation and rank candidate components by eigenvalue.
%   3. Score candidates using topography (occipital template agreement,
%      posterior concentration), spectral peak-form quality, and artifact
%      proxy metrics (front/temporal leak, line-harmonic ratio, stationarity,
%      burst ratio, HF slope, condition locking).
%   4. Apply adaptive hard rejection and dominant-outlier exclusion, then
%      build an eigenvalue-weighted combined occipital component set.
%
% Trial-level spectral scanning (per subject, condition, trial)
%   1. Project each trial to the combined GED component space and compute
%      spectrum-based power scans on a 30-90 Hz grid (default 0.5 Hz).
%      Scans are derived from spectral bins (FFT), not per-frequency FIR reruns.
%   2. Express power as robust dB ratio:
%      10*log10((p_stim + floor)/(p_base + floor)).
%   3. Flag numerical-instability cases (near-floor baseline power across
%      many frequencies/components) and exclude unstable trials automatically.
%   4. Detect per-trial spectral features (single-peak, low/high peaks,
%      centroid) and compute peak power.
%
% Outputs
% -------
%   - Per-condition/per-subject peak frequency and power summaries
%     (full, early, late windows).
%   - Combined-method benchmark summaries:
%     detectability, trial CV, peak power, and condition separation.
%   - Subject diagnostics (component selection, rejection reasons,
%     topographies, spectra), group summary figures, and GLMM check.
%   - Optional ground-truth simulation validation.
%   - Optional subject-level parallel execution (toggle in setup section).

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('GCP');
nSubj = length(subjects);

% Time windows
baseline_window = [-1.5, -0.25];
full_window = [0, 2.0];
early_window = [0, 0.6];
late_window = [1.0, 2.0];

% Gamma frequency range
gamma_range = [30, 90];

% Narrowband scanning parameters (direct analysis-band scan only)
analysis_freq_range = [30, 90];
scan_freq_step_hz = 0.5;            % analysis grid resolution (Hz)
scan_freqs = analysis_freq_range(1):scan_freq_step_hz:analysis_freq_range(2);
nFreqs = length(scan_freqs);
scan_width = 3;

%% GED parameters
lambda = 0.05;
ged_search_n = 20;          % search first N GED components
template_front_weight = 0.7; % anti-template weight for frontal channels
template_sigma_occ = 0.20;   % spatial smoothness for occipital template
template_sigma_front = 0.25; % spatial smoothness for frontal anti-template
min_eigval_hard = 1.1;              % hard minimum GED eigenvalue (lambda >= 1.1)
include_occipital_label_override = true; % enable lenient force-inclusion path for occipital-labeled components
min_peak_form_single_hard = 0.7;    % normal minimum PF score for candidate eligibility
max_minor_peaks_hard = 2;           % allow at most this many minor peaks for PF gate
max_minor_peak_rel_hard = 0.60;     % largest minor peak must stay below this fraction of dominant peak
occipital_pf_lenient_override = true;      % for occipital-labeled components, relax non-critical gates under lenient PF
occipital_pf_lenient_min = 0.60;           % lenient PF threshold for force-inclusion path
occipital_pf_fallback_min = 0.65;          % fallback PF threshold when regular pool is empty
max_frontleak_hard = 1.25;          % artifact guard: frontal leakage (front/occ)
max_templeak_hard = 1.35;           % artifact guard: temporal leakage (temp/occ)
min_corr_hard = -0.10;              % artifact guard: strongly anti-template components
max_lineharm_ratio_hard = 0.60;     % artifact guard: line-harmonic dominance ratio
max_stationarity_cv_hard = 1.20;    % artifact guard: trialwise gamma instability
max_burst_ratio_hard = 0.35;        % artifact guard: burst-dominated trials ratio
max_hf_slope_hard = -0.15;          % artifact guard: EMG-like high-frequency slope
min_condlock_rho_hard = 0.05;       % plausibility guard: weak condition locking
max_emg_score_hard = 0.85;          % adaptive guard anchor: very high EMG score
adaptive_class_quantile = 60;       % adaptive class split percentile (occ/eMG scores)
adaptive_artifact_mad_mult = 1.25;  % adaptive artifact bound scale (MAD units)
adaptive_min_eig_floor = 0.75;      % adaptive lower floor for GED eigenvalue gate
adaptive_occ_margin_floor = 0.05;   % adaptive floor for occ-vs-EMG separation margin
max_components_to_combine = 10;     % top-K cap for combined GED branch
artifact_proxy_max_trials = 24;     % speed cap for proxy estimation per component/window
exclude_dominant_outlier = true;    % remove dominant top outlier (lambda1/lambda2 + MAD)
outlier_ratio_thr = 8.0;            % lambda1/lambda2 threshold for dominant-outlier detection
outlier_mad_mult = 4.0;             % MAD multiplier on log-eigenvalue distance
outlier_min_rest = 0;               % minimum non-top eligible components for dominant-outlier exclusion (0 = no minimum)
baseline_outlier_exclusion_enable = true; % reject trials with unreliable baseline power estimates
baseline_outlier_mad_mult = 3.5;    % robust log-power cutoff (MAD units) for baseline outliers
baseline_outlier_min_trials = 12;   % minimum trial count before baseline outlier screening is applied
ratio_floor_prctile = 20;           % robust baseline-power percentile used as floor anchor
ratio_floor_frac = 0.25;            % floor is this fraction of the robust baseline-power anchor
instability_near_floor_mult = 1.5;  % mark as near-floor when baseline <= this multiple of floor
instability_component_frac_thr = 0.50; % unstable if >= this fraction of components are near-floor
instability_trial_freq_frac_thr = 0.35; % exclude trial when unstable at >= this frequency fraction
rank_stability_boot_reps = 200;     % bootstrap repetitions for rank-aggregation stability diagnostics
peak_bonus_rescue_min = 0.55;       % minimum peak-clarity to rescue flat-only exclusions
peak_bonus_min_peaks = 1;           % minimum number of detected peaks for rescue
hard_leak_severity_mult = 1.15;     % leak must exceed adaptive threshold by this factor for hard reject
hard_emg_severity_mult = 1.10;      % EMG score must exceed adaptive threshold by this factor for hard reject
occipital_pf_extreme_lineharm = 5.0;      % only reject occipital PF>=0.6 if lineharm exceeds this extreme threshold
occipital_pf_extreme_hf_slope = 2.0;      % only reject occipital PF>=0.6 if hf_slope exceeds this extreme threshold
occipital_pf_extreme_front_leak = 10.0;   % only reject occipital PF>=0.6 if front_leak exceeds this extreme threshold
occipital_pf_extreme_temp_leak = 10.0;    % only reject occipital PF>=0.6 if temp_leak exceeds this extreme threshold
topo_flat_std_min = 1e-6;           % near-zero std => degenerate/flat map
topo_flat_range_min = 1e-5;         % near-zero dynamic range => degenerate/flat map
topo_nonposterior_max = 0.28;       % posterior concentration below this is non-posterior
topo_fragmented_hotspot_rel = 8.0;  % high hotspot/median ratio => likely fragmented map
topo_occipital_override_min = 0.70; % strong posterior concentration override for class label
random_seed = 123;                    % reproducible randomization

% Raw dB spectrum peak detection parameters
centroid_freq_range = [40 80];
centroid_band_mask = scan_freqs >= centroid_freq_range(1) & scan_freqs <= centroid_freq_range(2);
centroid_posfrac_min = 0.20; % minimum positive-energy fraction for centroid validity
centroid_min_peak = 0.02;    % minimum positive peak in raw spectrum
trial_peak_smooth_n = 1;     % moving-average smoothing for trial-level peak detection
trial_metric_outlier_enable = true;      % apply trial-level outlier rejection on extracted frequency/power metrics
trial_metric_outlier_iqr_mult = 1.5;     % outlier threshold in IQR units around Q1/Q3
trial_metric_outlier_min_trials = 30;    % minimum finite trials to run IQR-based outlier rejection

% Condition info
condNames  = {'c25', 'c50', 'c75', 'c100'};
condLabels = {'25%', '50%', '75%', '100%'};

% Figure save directories
gcp_root_path = paths.root;
gcp_feature_data_path = fullfile(gcp_root_path, 'data', 'features');
if ~exist(gcp_feature_data_path, 'dir')
    gcp_feature_data_path = gcp_root_path;
end
fig_save_dir_ged = fullfile(gcp_root_path, 'figures', 'eeg', 'ged');
fig_save_dir_component_selection = fullfile(fig_save_dir_ged, 'component_selection');
fig_save_dir_emg_exclusion = fig_save_dir_component_selection;
fig_save_dir_powspctrm = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged/powspctrm';
if ~exist(fig_save_dir_ged, 'dir'), mkdir(fig_save_dir_ged); end
if ~exist(fig_save_dir_component_selection, 'dir'), mkdir(fig_save_dir_component_selection); end
if ~exist(fig_save_dir_powspctrm, 'dir'), mkdir(fig_save_dir_powspctrm); end
comp_sel_save_dir = fig_save_dir_component_selection;

%% Preallocate storage
all_trial_powratio     = cell(4, nSubj);
all_trial_powratio_fullscan = cell(4, nSubj);
all_trial_powratio_early = cell(4, nSubj);
all_trial_powratio_late  = cell(4, nSubj);
all_trial_powratio_plotstat = cell(4, nSubj);
all_trial_powratio_fullscan_plotstat = cell(4, nSubj);
all_trial_powratio_early_plotstat = cell(4, nSubj);
all_trial_powratio_late_plotstat  = cell(4, nSubj);
all_trial_peaks_single = cell(4, nSubj);
all_trial_peaks_low    = cell(4, nSubj);
all_trial_peaks_high   = cell(4, nSubj);
all_trial_peaks_single_early = cell(4, nSubj);
all_trial_peaks_single_late  = cell(4, nSubj);
all_trial_outlier_mask_freq_full = cell(4, nSubj);
all_trial_outlier_mask_freq_early = cell(4, nSubj);
all_trial_outlier_mask_freq_late = cell(4, nSubj);
all_trial_outlier_mask_power_full = cell(4, nSubj);
all_trial_outlier_mask_power_early = cell(4, nSubj);
all_trial_outlier_mask_power_late = cell(4, nSubj);
all_trial_centroid     = cell(4, nSubj);

all_trial_mean_single   = nan(4, nSubj);
all_trial_median_single = nan(4, nSubj);
all_trial_mean_low      = nan(4, nSubj);
all_trial_median_low    = nan(4, nSubj);
all_trial_mean_high     = nan(4, nSubj);
all_trial_median_high   = nan(4, nSubj);
all_trial_median_gap    = nan(4, nSubj);
all_trial_mean_single_early   = nan(4, nSubj);
all_trial_median_single_early = nan(4, nSubj);
all_trial_mean_single_late    = nan(4, nSubj);
all_trial_median_single_late  = nan(4, nSubj);
all_trial_median_low_early    = nan(4, nSubj);
all_trial_median_low_late     = nan(4, nSubj);
all_trial_median_high_early   = nan(4, nSubj);
all_trial_median_high_late    = nan(4, nSubj);
all_trial_median_gap_early    = nan(4, nSubj);
all_trial_median_gap_late     = nan(4, nSubj);
all_trial_median_centroid_early = nan(4, nSubj);
all_trial_median_centroid_late  = nan(4, nSubj);
all_trial_trialcv_early       = nan(4, nSubj);
all_trial_trialcv_late        = nan(4, nSubj);
all_trial_mean_centroid   = nan(4, nSubj);
all_trial_median_centroid = nan(4, nSubj);

all_trial_detrate_single = nan(4, nSubj);
all_trial_detrate_low    = nan(4, nSubj);
all_trial_detrate_high   = nan(4, nSubj);
all_trial_detrate_single_early = nan(4, nSubj);
all_trial_detrate_single_late  = nan(4, nSubj);
all_trial_detrate_centroid = nan(4, nSubj);
all_trial_gamma_power = nan(4, nSubj);
all_trial_gamma_power_early = nan(4, nSubj);
all_trial_gamma_power_late  = nan(4, nSubj);
all_trial_gamma_power_plotstat = nan(4, nSubj);
all_trial_gamma_power_early_plotstat = nan(4, nSubj);
all_trial_gamma_power_late_plotstat  = nan(4, nSubj);

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
warning_log_by_subj = cell(nSubj, 1);

all_trial_powratio_components_full  = cell(4, nSubj);
all_trial_powratio_components_early = cell(4, nSubj);
all_trial_powratio_components_late  = cell(4, nSubj);
benchmark_metric_detectability = nan(4, nSubj);
benchmark_metric_separation_slope = nan(nSubj, 1);
benchmark_metric_separation_delta = nan(nSubj, 1);
benchmark_metric_reliability_trialcv = nan(4, nSubj);
benchmark_metric_reliability_subjspread = nan(4, 1);
primary_slope_stats = struct();
primary_delta_stats = struct();

%% Subject loop
for subj = 1:nSubj
    close all
    tic
    comp_sel_save_dir = fullfile(fig_save_dir_component_selection, subjects{subj});
    if ~exist(comp_sel_save_dir, 'dir'), mkdir(comp_sel_save_dir); end
    fig_save_dir_emg_exclusion = comp_sel_save_dir;
    datapath = fullfile(gcp_feature_data_path, subjects{subj}, 'eeg');
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
    clc
    fprintf('Subject %s (%d/%d) GED (early, full, late) (%d occ channels)\n', ...
        subjects{subj}, subj, nSubj, nOcc);
    rng(random_seed + subj, 'twister');

    stim_windows = {full_window, early_window, late_window};
    win_names   = {'full', 'early', 'late'};
    lambdas     = [lambda, lambda, lambda];

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

    % Run GED + component selection per window
    searchFilters_full  = [];
    searchFilters_early = [];
    searchFilters_late  = [];
    W_top_full  = []; W_top_early  = []; W_top_late  = [];
    W_combined_full  = []; W_combined_early  = []; W_combined_late  = [];
    selected_idx_full  = []; selected_idx_early  = []; selected_idx_late  = [];
    w_combined_full  = []; w_combined_early  = []; w_combined_late  = [];
    evals_sorted_full = []; evals_sorted_early = []; evals_sorted_late = [];
    searchCorrs_full = []; searchCorrs_early = []; searchCorrs_late = [];
    searchTopos_full = []; searchTopos_early = []; searchTopos_late = [];
    searchMeanPrSpectrum_full = []; searchMeanPrSpectrum_early = []; searchMeanPrSpectrum_late = [];
    searchEmgClass_full = {}; searchEmgClass_early = {}; searchEmgClass_late = {};
    rejection_flags_full = struct(); rejection_flags_early = struct(); rejection_flags_late = struct();
    soft_warn_flags_full = struct(); soft_warn_flags_early = struct(); soft_warn_flags_late = struct();
    candidate_table_full = struct(); candidate_table_early = struct(); candidate_table_late = struct();

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
        searchOccFrontRatio = nan(nSearch, 1);
        searchFrontLeak = nan(nSearch, 1);
        searchTempStrength = nan(nSearch, 1);
        searchTempLeak = nan(nSearch, 1);
        searchLineHarmRatio = nan(nSearch, 1);
        searchStationarityCV = nan(nSearch, 1);
        searchBurstRatio = nan(nSearch, 1);
        searchHFSlope = nan(nSearch, 1);
        searchCondLockRho = nan(nSearch, 1);
        searchProxyUnknown = false(nSearch, 1);
        searchEmgClass = repmat({'unassigned'}, nSearch, 1);
        searchMeanPrSpectrum = nan(nSearch, numel(scan_freqs));
        searchTopoPosteriorConcentration = nan(nSearch, 1);
        searchTopoSpatialStd = nan(nSearch, 1);
        searchTopoDynamicRange = nan(nSearch, 1);
        searchTopoHotspotRel = nan(nSearch, 1);

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
        topo_spatial_std_ci = std(topo_ci(isfinite(topo_ci)));
        if ~isfinite(topo_spatial_std_ci)
            topo_spatial_std_ci = 0;
        end
        topo_dynamic_range_ci = max(topo_finite) - min(topo_finite);
        posterior_mass = mean(topo_abs(post_idx));
        total_mass = mean(topo_finite);
        topo_posterior_concentration_ci = posterior_mass / max(total_mass, eps);
        hotspot_thr_ci = prctile(topo_finite, 95);
        hotspot_mass_ci = mean(topo_finite(topo_finite >= hotspot_thr_ci));
        topo_hotspot_rel_ci = hotspot_mass_ci / max(median(topo_finite), eps);
        % Lightweight artifact proxies from trial-level component spectra.
        proxy_ci = estimate_component_artifact_proxies( ...
            w_ci, dat_per_cond, stim_windows{w}, baseline_window, scan_freqs, scan_width, artifact_proxy_max_trials);

        searchFilters(:, ci) = w_ci;
        searchTopos(:, ci) = topo_ci;
        searchCorrs(ci) = r_ci;
        searchOccStrength(ci) = occ_strength;
        searchFrontStrength(ci) = front_strength;
        searchOccFrontRatio(ci) = ratio_ci;
        searchFrontLeak(ci) = front_leak_ci;
        searchTempStrength(ci) = temp_strength;
        searchTempLeak(ci) = temp_leak_ci;
        searchLineHarmRatio(ci) = proxy_ci.lineharm_ratio;
        searchStationarityCV(ci) = proxy_ci.stationarity_cv;
        searchBurstRatio(ci) = proxy_ci.burst_ratio;
        searchHFSlope(ci) = proxy_ci.hf_slope;
        searchCondLockRho(ci) = proxy_ci.cond_lock_rho;
        searchProxyUnknown(ci) = proxy_ci.unknown_high_risk;
        searchMeanPrSpectrum(ci, :) = proxy_ci.mean_pr_spectrum(:)';
        searchTopoPosteriorConcentration(ci) = topo_posterior_concentration_ci;
        searchTopoSpatialStd(ci) = topo_spatial_std_ci;
        searchTopoDynamicRange(ci) = topo_dynamic_range_ci;
        searchTopoHotspotRel(ci) = topo_hotspot_rel_ci;
    end
    W_top = searchFilters(:, 1);

    % Candidate metrics (selection = eig + PF core gates + explicit artifact rejection).
    corr_vec = searchCorrs;
    ratio_vec = searchOccFrontRatio;
    eval_raw_vec = evals_sorted(1:nSearch);
    eval_vec = log(max(eval_raw_vec, eps));
    leak_vec = searchFrontLeak;
    temp_leak_vec = searchTempLeak;
    lineharm_vec = searchLineHarmRatio;
    stationarity_vec = searchStationarityCV;
    burst_vec = searchBurstRatio;
    hf_slope_vec = searchHFSlope;
    condlock_vec = searchCondLockRho;
    unknown_proxy_vec = logical(searchProxyUnknown(:));
    % Non-finite proxy values are explicitly marked as unknown/high-risk.
    unknown_proxy_vec = unknown_proxy_vec | ~isfinite(lineharm_vec) | ~isfinite(stationarity_vec) | ...
        ~isfinite(burst_vec) | ~isfinite(hf_slope_vec);
    finite_metrics = isfinite(corr_vec) & isfinite(ratio_vec) & ...
        isfinite(eval_raw_vec) & isfinite(leak_vec) & isfinite(temp_leak_vec) & ...
        isfinite(lineharm_vec) & isfinite(stationarity_vec) & isfinite(burst_vec) & ...
        isfinite(hf_slope_vec);
    [peak_bonus_vec, peak_count_vec] = compute_peak_bonus_from_spectra( ...
        searchMeanPrSpectrum, scan_freqs, analysis_freq_range);
    [peak_form_score_vec, peak_form_mode_vec, peak_form_diag] = compute_peak_form_template_score_from_spectra( ...
        searchMeanPrSpectrum, scan_freqs, analysis_freq_range);
    topo_posterior_vec = searchTopoPosteriorConcentration;
    topo_flat_fail_vec = (searchTopoSpatialStd < topo_flat_std_min) | (searchTopoDynamicRange < topo_flat_range_min);
    topo_nonposterior_fail_vec = topo_posterior_vec < topo_nonposterior_max;
    topo_fragmented_fail_vec = searchTopoHotspotRel > topo_fragmented_hotspot_rel;
    occipital_evidence = 0.40 * normalize_robust(corr_vec) + ...
        0.25 * normalize_robust(ratio_vec) + ...
        0.35 * normalize_robust(topo_posterior_vec);
    emg_artifact_score = 0.30 * normalize_robust(leak_vec) + ...
        0.20 * normalize_robust(temp_leak_vec) + ...
        0.20 * normalize_robust(lineharm_vec) + ...
        0.15 * normalize_robust(stationarity_vec) + ...
        0.15 * normalize_robust(max(hf_slope_vec, 0));
    default_thresholds = struct( ...
        'min_eigval', min_eigval_hard, ...
        'max_frontleak', max_frontleak_hard, ...
        'max_templeak', max_templeak_hard, ...
        'min_corr', min_corr_hard, ...
        'max_lineharm', max_lineharm_ratio_hard, ...
        'max_stationarity', max_stationarity_cv_hard, ...
        'max_burst_ratio', max_burst_ratio_hard, ...
        'max_hf_slope', max_hf_slope_hard, ...
        'min_condlock', min_condlock_rho_hard, ...
        'max_emg_score', max_emg_score_hard, ...
        'class_quantile', adaptive_class_quantile, ...
        'artifact_mad_mult', adaptive_artifact_mad_mult, ...
        'min_eig_floor', adaptive_min_eig_floor, ...
        'occ_margin_floor', adaptive_occ_margin_floor);
    adaptive_thr = build_adaptive_component_thresholds( ...
        eval_raw_vec, corr_vec, leak_vec, temp_leak_vec, lineharm_vec, stationarity_vec, ...
        burst_vec, hf_slope_vec, condlock_vec, occipital_evidence, emg_artifact_score, default_thresholds);
    enforced_min_eigval = adaptive_thr.min_eigval;
    pass_eig_gate = finite_metrics & (eval_raw_vec >= enforced_min_eigval);
    single_peak_mode_mask = strcmpi(peak_form_mode_vec, 'single') | strcmpi(peak_form_mode_vec, 'dominant');
    pass_single_peak_gate = finite_metrics & (peak_form_score_vec >= min_peak_form_single_hard) & ...
        (peak_form_diag.minor_peak_count <= max_minor_peaks_hard) & ...
        (peak_form_diag.minor_peak_relmax <= max_minor_peak_rel_hard);
    hard_eligible_raw = pass_eig_gate & pass_single_peak_gate;
    fail_front_leak = finite_metrics & (leak_vec > adaptive_thr.max_frontleak);
    fail_temp_leak = finite_metrics & (temp_leak_vec > adaptive_thr.max_templeak);
    fail_corr = finite_metrics & (corr_vec < adaptive_thr.min_corr);
    fail_lineharm = finite_metrics & (lineharm_vec > adaptive_thr.max_lineharm);
    fail_stationarity = finite_metrics & (stationarity_vec > adaptive_thr.max_stationarity);
    fail_burst = finite_metrics & (burst_vec > adaptive_thr.max_burst_ratio);
    fail_hf_slope = finite_metrics & (hf_slope_vec > adaptive_thr.max_hf_slope);
    fail_condlock = finite_metrics & (condlock_vec < adaptive_thr.min_condlock);
    fail_emg_score = finite_metrics & (emg_artifact_score > adaptive_thr.max_emg_score);
    severe_front_leak = finite_metrics & (leak_vec > adaptive_thr.max_frontleak * hard_leak_severity_mult);
    severe_temp_leak = finite_metrics & (temp_leak_vec > adaptive_thr.max_templeak * hard_leak_severity_mult);
    severe_emg_score = finite_metrics & (emg_artifact_score > adaptive_thr.max_emg_score * hard_emg_severity_mult);
    extreme_lineharm = finite_metrics & (lineharm_vec > occipital_pf_extreme_lineharm);
    extreme_hf_slope = finite_metrics & (hf_slope_vec > occipital_pf_extreme_hf_slope);
    extreme_front_leak = finite_metrics & (leak_vec > occipital_pf_extreme_front_leak);
    extreme_temp_leak = finite_metrics & (temp_leak_vec > occipital_pf_extreme_temp_leak);
    occ_minus_emg = occipital_evidence - emg_artifact_score;
    fail_occ_margin = finite_metrics & (emg_artifact_score >= adaptive_thr.emg_class_thr) & ...
        (occ_minus_emg < adaptive_thr.min_occ_margin);
    hard_reject_flags_raw = severe_front_leak | severe_temp_leak | fail_lineharm | fail_hf_slope | ...
        severe_emg_score | unknown_proxy_vec | topo_flat_fail_vec | ...
        topo_nonposterior_fail_vec | topo_fragmented_fail_vec;
    soft_warn_any = fail_front_leak | fail_temp_leak | fail_corr | fail_stationarity | ...
        fail_burst | fail_condlock | fail_emg_score | fail_occ_margin;
    soft_warn_flags = struct( ...
        'front_leak', fail_front_leak & ~severe_front_leak, ...
        'temp_leak', fail_temp_leak & ~severe_temp_leak, ...
        'corr', fail_corr, ...
        'stationarity', fail_stationarity, ...
        'burst', fail_burst, ...
        'condlock', fail_condlock, ...
        'emg_score', fail_emg_score & ~severe_emg_score, ...
        'occ_margin', fail_occ_margin, ...
        'any', soft_warn_any);
    artifact_flags = hard_reject_flags_raw;
    rejection_flags = struct( ...
        'unknown_proxy', unknown_proxy_vec, ...
        'front_leak', fail_front_leak, ...
        'temp_leak', fail_temp_leak, ...
        'corr', fail_corr, ...
        'lineharm', fail_lineharm, ...
        'stationarity', fail_stationarity, ...
        'burst', fail_burst, ...
        'hf_slope', fail_hf_slope, ...
        'condlock', fail_condlock, ...
        'emg_score', fail_emg_score, ...
        'occ_margin', fail_occ_margin, ...
        'severe_front_leak', severe_front_leak, ...
        'severe_temp_leak', severe_temp_leak, ...
        'severe_emg_score', severe_emg_score, ...
        'topo_flat', topo_flat_fail_vec, ...
        'topo_nonposterior', topo_nonposterior_fail_vec, ...
        'topo_fragmented', topo_fragmented_fail_vec);
    for ci = 1:nSearch
        if (topo_posterior_vec(ci) >= topo_occipital_override_min) && ...
                ~topo_flat_fail_vec(ci) && ~topo_fragmented_fail_vec(ci) && ...
                ~topo_nonposterior_fail_vec(ci)
            % Strongly posterior-constrained map: force occipital class label.
            searchEmgClass{ci} = 'occipital';
        elseif unknown_proxy_vec(ci)
            searchEmgClass{ci} = 'unclear';
        elseif (occipital_evidence(ci) >= adaptive_thr.occ_class_thr) && ...
                (emg_artifact_score(ci) < adaptive_thr.emg_class_thr) && ...
                (occ_minus_emg(ci) >= adaptive_thr.min_occ_margin)
            searchEmgClass{ci} = 'occipital';
        elseif (occipital_evidence(ci) < adaptive_thr.occ_class_thr) && ...
                (emg_artifact_score(ci) >= adaptive_thr.emg_class_thr) && ...
                (-occ_minus_emg(ci) >= adaptive_thr.min_occ_margin)
            searchEmgClass{ci} = 'EMG';
        elseif (occipital_evidence(ci) >= adaptive_thr.occ_class_thr) && ...
                (emg_artifact_score(ci) >= adaptive_thr.emg_class_thr)
            if abs(occ_minus_emg(ci)) < adaptive_thr.min_occ_margin && ...
                    (topo_posterior_vec(ci) >= topo_nonposterior_max)
                searchEmgClass{ci} = 'mixed';
            elseif occ_minus_emg(ci) >= 0
                searchEmgClass{ci} = 'occipital';
            else
                searchEmgClass{ci} = 'EMG';
            end
        else
            searchEmgClass{ci} = 'unclear';
        end
    end
    dominant_outlier_mask = false(nSearch, 1);
    if exclude_dominant_outlier
        [~, dominant_outlier_idx] = exclude_dominant_top_outlier( ...
            evals_sorted(1:nSearch), hard_eligible_raw, outlier_ratio_thr, outlier_mad_mult, outlier_min_rest);
        if ~isempty(dominant_outlier_idx)
            dominant_outlier_mask(dominant_outlier_idx) = true;
            artifact_flags(dominant_outlier_idx) = true;
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
    no_hard_threshold_match = ~any(hard_eligible_raw);
    if ~any(pass_single_peak_gate)
        pf_msg = sprintf('Subject %s has no components passing PF gate.', subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'NO_COMPONENTS_PASS_PF_GATE', pf_msg, ...
            struct('n_search', nSearch, 'n_pass_pf', sum(pass_single_peak_gate), ...
                   'thr_single_peak_pf', min_peak_form_single_hard, ...
                   'thr_max_minor_peaks', max_minor_peaks_hard, ...
                   'thr_max_minor_peak_rel', max_minor_peak_rel_hard));
    end
    hard_eligible = hard_eligible_raw & ~artifact_flags;
    hard_eligible(dominant_outlier_mask) = false;
    occipital_class_mask = cellfun(@(c) strcmpi(c, 'occipital'), searchEmgClass(:));
    very_good_pf_mask = finite_metrics & (peak_form_score_vec >= occipital_pf_lenient_min);
    pass_lenient_eig_gate = finite_metrics & (eval_raw_vec >= enforced_min_eigval);
    force_include_antiartifact_set1 = finite_metrics & ~unknown_proxy_vec & ...
        ~severe_front_leak & ~severe_temp_leak;
    force_include_antiartifact_set2 = finite_metrics & ~unknown_proxy_vec & ...
        ~fail_lineharm & ~fail_hf_slope & ~topo_fragmented_fail_vec;
    force_include_strong_ant_artifact_mask = force_include_antiartifact_set1 | force_include_antiartifact_set2;
    lenient_occipital_pf_mask = occipital_pf_lenient_override & occipital_class_mask & very_good_pf_mask & ...
        pass_lenient_eig_gate & ~dominant_outlier_mask & ~severe_emg_score & ...
        force_include_strong_ant_artifact_mask;
    force_include_occipital_mask = include_occipital_label_override & lenient_occipital_pf_mask;
    rescue_occipital_mask = false(nSearch, 1);
    % Lenient override intent: occipital + PF>=lenient threshold + non-severe EMG
    % remains selectable even if raw hard-reject criteria are triggered.
    hard_reject_flags = hard_reject_flags_raw & ~force_include_occipital_mask;
    artifact_flags = hard_reject_flags;
    selection_pool_mask = (hard_eligible | force_include_occipital_mask) & occipital_class_mask;
    [searchScores, rankagg_metrics, rankagg_stability] = compute_calibrated_rank_aggregation_score( ...
        eval_raw_vec, peak_form_score_vec, peak_bonus_vec, occipital_evidence, emg_artifact_score, rank_stability_boot_reps);
    searchScores(~finite_metrics) = -Inf;
    searchScores(~selection_pool_mask) = -Inf;
    if ~any(selection_pool_mask) && ~any(hard_eligible) && ~any(force_include_occipital_mask) && any(occipital_class_mask)
        rescue_base_mask = occipital_class_mask & finite_metrics & pass_lenient_eig_gate & ...
            ~dominant_outlier_mask & ~unknown_proxy_vec & ...
            ~severe_front_leak & ~severe_temp_leak & ~severe_emg_score;
        if any(rescue_base_mask)
            rescue_scores = rankagg_metrics.score_raw;
            rescue_scores(~rescue_base_mask) = -Inf;
            if ~any(isfinite(rescue_scores))
                rescue_scores = eval_raw_vec;
                rescue_scores(~rescue_base_mask) = -Inf;
            end
            if any(isfinite(rescue_scores))
                [~, rescue_ord] = sort(rescue_scores, 'descend');
                rescue_ord = rescue_ord(isfinite(rescue_scores(rescue_ord)));
                n_rescue = min(2, numel(rescue_ord));
                if n_rescue > 0
                    rescue_idx = rescue_ord(1:n_rescue);
                    rescue_occipital_mask(rescue_idx) = true;
                    selection_pool_mask = selection_pool_mask | rescue_occipital_mask;
                    searchScores = rankagg_metrics.score_raw;
                    searchScores(~finite_metrics) = -Inf;
                    searchScores(~selection_pool_mask) = -Inf;
                    msg = sprintf(['Rescue path activated for subject %s: selected %d finite non-severe ', ...
                        'occipital component(s) after hard/force pool collapse.'], subjects{subj}, n_rescue);
                    warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, ...
                        'RESCUE_OCCIPITAL_NON_SEVERE', msg, ...
                        struct('n_rescue', n_rescue, 'n_rescue_pool', sum(rescue_base_mask), ...
                               'n_occipital_class', sum(occipital_class_mask)));
                end
            end
        end
    end
    if ~any(isfinite(searchScores))
        msg = sprintf(['No components met eig/artifact criteria for subject %s. ', ...
                       'Falling back to unconstrained eigenvalue ranking.'], subjects{subj});
        top_fail_idx = NaN;
        top_fail_eig = NaN;
        top_fail_corr = NaN;
        top_fail_ratio = NaN;
        top_fail_pass_eig = false;
        top_fail_artifact = false;
        finite_idx = find(finite_metrics);
        if ~isempty(finite_idx)
            [~, fail_ord] = sort(evals_sorted(finite_idx), 'descend');
            top_fail_idx = finite_idx(fail_ord(1));
            top_fail_eig = evals_sorted(top_fail_idx);
            top_fail_corr = corr_vec(top_fail_idx);
            top_fail_ratio = ratio_vec(top_fail_idx);
            top_fail_pass_eig = top_fail_eig >= enforced_min_eigval;
            top_fail_artifact = artifact_flags(top_fail_idx);
        end
        hard_metrics = struct( ...
            'n_search', nSearch, ...
            'n_finite_metrics', sum(finite_metrics), ...
            'n_pass_eig', sum(pass_eig_gate), ...
            'n_pass_single_peak_pf', sum(pass_single_peak_gate), ...
            'n_artifact_flagged', sum(artifact_flags), ...
            'n_unknown_high_risk', sum(unknown_proxy_vec), ...
            'n_pass_all_raw', sum(hard_eligible_raw), ...
            'n_occipital_class', sum(occipital_class_mask), ...
            'n_single_peak_mode', sum(single_peak_mode_mask), ...
            'n_lenient_occipital_pf', sum(lenient_occipital_pf_mask), ...
            'n_force_include_occipital', sum(force_include_occipital_mask), ...
            'n_excluded_dominant_outlier', sum(dominant_outlier_mask), ...
            'top_fail_idx', top_fail_idx, ...
            'top_fail_eig', top_fail_eig, ...
            'top_fail_corr', top_fail_corr, ...
            'top_fail_ratio', top_fail_ratio, ...
            'top_fail_pass_eig', top_fail_pass_eig, ...
            'top_fail_artifact', top_fail_artifact, ...
            'thr_eig', enforced_min_eigval, ...
            'thr_single_peak_pf', min_peak_form_single_hard, ...
            'thr_corr', adaptive_thr.min_corr, ...
            'thr_front_leak', adaptive_thr.max_frontleak, ...
            'thr_temp_leak', adaptive_thr.max_templeak, ...
            'thr_lineharm', adaptive_thr.max_lineharm, ...
            'thr_stationarity', adaptive_thr.max_stationarity, ...
            'thr_burst_ratio', adaptive_thr.max_burst_ratio, ...
            'thr_hf_slope', adaptive_thr.max_hf_slope, ...
            'thr_condlock', adaptive_thr.min_condlock, ...
            'thr_emg_score', adaptive_thr.max_emg_score, ...
            'thr_occ_margin', adaptive_thr.min_occ_margin);
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'NO_HARD_ELIGIBLE_COMPONENTS', msg, hard_metrics);
        searchScores = rankagg_metrics.score_raw;
        searchScores(~finite_metrics) = -Inf;
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
    candidate_table.ratio = ratio_vec;
    candidate_table.front_leak = leak_vec;
    candidate_table.temp_leak = temp_leak_vec;
    candidate_table.lineharm_ratio = lineharm_vec;
    candidate_table.stationarity_cv = stationarity_vec;
    candidate_table.burst_ratio = burst_vec;
    candidate_table.hf_slope = hf_slope_vec;
    candidate_table.cond_lock_rho = condlock_vec;
    candidate_table.peak_clarity = peak_bonus_vec;
    candidate_table.peak_count = peak_count_vec;
    candidate_table.peak_form_score = peak_form_score_vec;
    candidate_table.peak_form_mode = peak_form_mode_vec;
    candidate_table.peak_form_basis = repmat({'raw_only_current'}, nSearch, 1);
    candidate_table.peak_form_best_single_similarity = peak_form_diag.best_single_similarity;
    candidate_table.peak_form_best_double_similarity = peak_form_diag.best_double_similarity;
    candidate_table.peak_form_similarity_raw = peak_form_diag.best_similarity_raw;
    candidate_table.peak_form_best_shift_hz = peak_form_diag.best_shift_hz;
    candidate_table.peak_form_best_center_hz = peak_form_diag.best_center_hz;
    candidate_table.peak_form_best_width_hz = peak_form_diag.best_width_hz;
    candidate_table.peak_form_best_separation_hz = peak_form_diag.best_separation_hz;
    candidate_table.peak_form_trough_depth = peak_form_diag.best_trough_depth;
    candidate_table.peak_form_minor_peak_count = peak_form_diag.minor_peak_count;
    candidate_table.peak_form_minor_peak_relmax = peak_form_diag.minor_peak_relmax;
    candidate_table.peak_form_edge_ratio = peak_form_diag.edge_ratio;
    candidate_table.peak_form_edge_run = peak_form_diag.edge_run_score;
    candidate_table.peak_form_edge_artifact_flag = peak_form_diag.edge_artifact_flag;
    candidate_table.peak_form_similarity_pre_penalty = peak_form_diag.pre_penalty_similarity;
    candidate_table.peak_form_total_penalty_raw = peak_form_diag.total_penalty_raw;
    candidate_table.peak_form_total_penalty_used = peak_form_diag.total_penalty_used;
    candidate_table.peak_form_dominant_penalty = peak_form_diag.dominant_penalty_tag;
    candidate_table.pf_quality_band = assign_pf_quality_band_labels(peak_form_score_vec);
    candidate_table.pf_rank = compute_descending_rank(peak_form_score_vec);
    candidate_table.occipital_evidence = occipital_evidence;
    candidate_table.topo_posterior_concentration = topo_posterior_vec;
    candidate_table.topo_spatial_std = searchTopoSpatialStd;
    candidate_table.topo_dynamic_range = searchTopoDynamicRange;
    candidate_table.topo_hotspot_rel = searchTopoHotspotRel;
    candidate_table.emg_artifact_score = emg_artifact_score;
    candidate_table.emg_class = searchEmgClass;
    candidate_table.unknown_high_risk = unknown_proxy_vec;
    candidate_table.score = searchScores;
    candidate_table.score_base = searchScores;
    candidate_table.score_final = searchScores;
    candidate_table.score_rankagg_raw = rankagg_metrics.score_raw;
    candidate_table.score_rankagg_eig = rankagg_metrics.eig_score;
    candidate_table.score_rankagg_peak_form = rankagg_metrics.peak_form_score;
    candidate_table.score_rankagg_peak_bonus = rankagg_metrics.peak_bonus_score;
    candidate_table.score_rankagg_occipital = rankagg_metrics.occipital_score;
    candidate_table.score_rankagg_anti_emg = rankagg_metrics.anti_emg_score;
    candidate_table.rank_stability_top1_freq = rankagg_stability.top1_freq;
    candidate_table.rank_stability_rank_mean = rankagg_stability.rank_mean;
    candidate_table.rank_stability_rank_std = rankagg_stability.rank_std;
    candidate_table.pass_eig_gate = pass_eig_gate;
    candidate_table.pass_single_peak_pf_gate = pass_single_peak_gate;
    candidate_table.single_peak_mode = single_peak_mode_mask;
    candidate_table.very_good_pf = very_good_pf_mask;
    candidate_table.pass_lenient_eig_gate = pass_lenient_eig_gate;
    candidate_table.force_include_antiartifact_set1 = force_include_antiartifact_set1;
    candidate_table.force_include_antiartifact_set2 = force_include_antiartifact_set2;
    candidate_table.force_include_strong_ant_artifact = force_include_strong_ant_artifact_mask;
    candidate_table.lenient_occipital_pf = lenient_occipital_pf_mask;
    candidate_table.force_include_occipital = force_include_occipital_mask;
    candidate_table.rescue_occipital_non_severe = rescue_occipital_mask;
    candidate_table.artifact_flag = artifact_flags;
    candidate_table.hard_reject_raw = hard_reject_flags_raw;
    candidate_table.hard_reject = hard_reject_flags;
    candidate_table.soft_warn = soft_warn_any;
    candidate_table.reject_reason = compute_primary_rejection_reason(rejection_flags);
    candidate_table.fail_unknown_proxy = rejection_flags.unknown_proxy;
    candidate_table.fail_front_leak = rejection_flags.front_leak;
    candidate_table.fail_temp_leak = rejection_flags.temp_leak;
    candidate_table.fail_corr = rejection_flags.corr;
    candidate_table.fail_lineharm = rejection_flags.lineharm;
    candidate_table.fail_stationarity = rejection_flags.stationarity;
    candidate_table.fail_burst = rejection_flags.burst;
    candidate_table.fail_hf_slope = rejection_flags.hf_slope;
    candidate_table.fail_condlock = rejection_flags.condlock;
    candidate_table.fail_emg_score = rejection_flags.emg_score;
    candidate_table.fail_occ_margin = rejection_flags.occ_margin;
    candidate_table.fail_severe_front_leak = rejection_flags.severe_front_leak;
    candidate_table.fail_severe_temp_leak = rejection_flags.severe_temp_leak;
    candidate_table.fail_severe_emg_score = rejection_flags.severe_emg_score;
    candidate_table.fail_topo_flat = rejection_flags.topo_flat;
    candidate_table.fail_topo_nonposterior = rejection_flags.topo_nonposterior;
    candidate_table.fail_topo_fragmented = rejection_flags.topo_fragmented;
    candidate_table.soft_warn_front_leak = soft_warn_flags.front_leak;
    candidate_table.soft_warn_temp_leak = soft_warn_flags.temp_leak;
    candidate_table.soft_warn_corr = soft_warn_flags.corr;
    candidate_table.soft_warn_stationarity = soft_warn_flags.stationarity;
    candidate_table.soft_warn_burst = soft_warn_flags.burst;
    candidate_table.soft_warn_condlock = soft_warn_flags.condlock;
    candidate_table.soft_warn_emg_score = soft_warn_flags.emg_score;
    candidate_table.soft_warn_occ_margin = soft_warn_flags.occ_margin;
    candidate_table.hard_eligible_raw = hard_eligible_raw;
    candidate_table.dominant_outlier = dominant_outlier_mask;
    candidate_table.hard_eligible = hard_eligible;
    candidate_table.thr_min_eigval = repmat(enforced_min_eigval, nSearch, 1);
    candidate_table.thr_min_single_peak_pf = repmat(min_peak_form_single_hard, nSearch, 1);
    candidate_table.thr_max_minor_peaks = repmat(max_minor_peaks_hard, nSearch, 1);
    candidate_table.thr_max_minor_peak_rel = repmat(max_minor_peak_rel_hard, nSearch, 1);
    candidate_table.thr_occ_class = repmat(adaptive_thr.occ_class_thr, nSearch, 1);
    candidate_table.thr_emg_class = repmat(adaptive_thr.emg_class_thr, nSearch, 1);
    candidate_table.thr_occ_margin = repmat(adaptive_thr.min_occ_margin, nSearch, 1);
    combined_idx = find(selection_pool_mask & isfinite(searchScores));
    fallback_occipital_idx = NaN;
    fallback_selected_mask = false(nSearch, 1);
    if isempty(combined_idx)
        msg = sprintf(['No occipital-labeled artifact-screened finite components available for subject %s. ', ...
            'This window will be marked as NaN for downstream metrics.'], subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'NO_OCCIPITAL_COMPONENTS', msg, ...
            struct('n_hard_eligible', sum(hard_eligible), 'n_occipital_labeled', sum(occipital_class_mask), ...
            'n_finite_scores', sum(isfinite(searchScores))));
    end
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
    candidate_table.fallback_selected = fallback_selected_mask;
    candidate_table.score_base = searchScores;
    candidate_table.score_final = searchScores;
    candidate_table.score = searchScores;

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

    finite_scores = find(isfinite(searchScores));
    if isempty(finite_scores)
        msg = sprintf('No finite eligible component scores for subject %s; using unconstrained eigenvalue ordering for display.', subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'NO_FINITE_SCORES_FOR_DISPLAY', msg, ...
            struct('n_hard_eligible', sum(hard_eligible), 'n_search', nSearch));
        [~, topDispOrder] = sort(evals_sorted(1:nSearch), 'descend');
    else
        [~, topDispOrder] = sort(evals_sorted(1:nSearch), 'descend');
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

        % Store per-window filters for trial-level scanning
        if w == 1
            searchFilters_full = searchFilters;
            W_top_full = searchFilters(:, 1);
            W_combined_full = searchFilters(:, selected_idx);
            selected_idx_full = selected_idx;
            w_combined_full = selected_weights(:)';
            evals_sorted_full = evals_sorted;
            searchCorrs_full = searchCorrs;
            searchTopos_full = searchTopos;
            searchMeanPrSpectrum_full = searchMeanPrSpectrum;
            searchEmgClass_full = searchEmgClass;
            rejection_flags_full = rejection_flags;
            soft_warn_flags_full = soft_warn_flags;
            candidate_table_full = candidate_table;
        elseif w == 2
            searchFilters_early = searchFilters;
            W_top_early = searchFilters(:, 1);
            W_combined_early = searchFilters(:, selected_idx);
            selected_idx_early = selected_idx;
            w_combined_early = selected_weights(:)';
            evals_sorted_early = evals_sorted;
            searchCorrs_early = searchCorrs;
            searchTopos_early = searchTopos;
            searchMeanPrSpectrum_early = searchMeanPrSpectrum;
            searchEmgClass_early = searchEmgClass;
            rejection_flags_early = rejection_flags;
            soft_warn_flags_early = soft_warn_flags;
            candidate_table_early = candidate_table;
        else
            searchFilters_late = searchFilters;
            W_top_late = searchFilters(:, 1);
            W_combined_late = searchFilters(:, selected_idx);
            selected_idx_late = selected_idx;
            w_combined_late = selected_weights(:)';
            evals_sorted_late = evals_sorted;
            searchCorrs_late = searchCorrs;
            searchTopos_late = searchTopos;
            searchMeanPrSpectrum_late = searchMeanPrSpectrum;
            searchEmgClass_late = searchEmgClass;
            rejection_flags_late = rejection_flags;
            soft_warn_flags_late = soft_warn_flags;
            candidate_table_late = candidate_table;
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
            'selection_mode', 'eig_artifact_screened_weighted', ...
            'selected_idx', selected_idx, ...
            'fallback_selected_idx', fallback_occipital_idx, ...
            'selected_weights', selected_weights, ...
            'candidate_table', candidate_table, ...
            'best_idx', bestIdx, ...
            'best_score', bestScore, ...
            'best_corr', bestCorr, ...
            'best_ratio', bestRatio, ...
            'best_front', bestFront, ...
            'best_occ', bestOcc, ...
            'best_leak', bestLeak, ...
            'emg_artifact_score', emg_artifact_score, ...
            'occipital_evidence', occipital_evidence, ...
            'emg_class', {searchEmgClass}, ...
            'unknown_high_risk', unknown_proxy_vec, ...
            'adaptive_thresholds', adaptive_thr, ...
            'hard_reject_flags', artifact_flags, ...
            'soft_warn_flags', soft_warn_flags, ...
            'rejection_flags', rejection_flags, ...
            'no_hard_threshold_match', no_hard_threshold_match);
        if w == 1
            all_component_selection_stats_full{subj} = comp_sel_struct;
        elseif w == 2
            all_component_selection_stats_early{subj} = comp_sel_struct;
        else
            all_component_selection_stats_late{subj} = comp_sel_struct;
        end

    % EMG exclusion diagnostics are generated after all windows are processed.
    end

    W_combined_full = searchFilters_full(:, selected_idx_full);
    W_combined_early = searchFilters_early(:, selected_idx_early);
    W_combined_late = searchFilters_late(:, selected_idx_late);
    W_top_full = searchFilters_full(:, 1);
    W_top_early = searchFilters_early(:, 1);
    W_top_late = searchFilters_late(:, 1);

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
    pf_benchmark = build_pf_benchmark_summary(candidate_table_full, candidate_table_early, candidate_table_late);
    all_component_selection_stats_full{subj}.selected_idx = selected_idx_full;
    all_component_selection_stats_full{subj}.selected_weights = w_combined_full;
    all_component_selection_stats_full{subj}.candidate_table = candidate_table_full;
    all_component_selection_stats_full{subj}.pf_benchmark = pf_benchmark.full;

    all_component_selection_stats_early{subj}.selected_idx = selected_idx_early;
    all_component_selection_stats_early{subj}.selected_weights = w_combined_early;
    all_component_selection_stats_early{subj}.candidate_table = candidate_table_early;
    all_component_selection_stats_early{subj}.pf_benchmark = pf_benchmark.early;

    all_component_selection_stats_late{subj}.selected_idx = selected_idx_late;
    all_component_selection_stats_late{subj}.selected_weights = w_combined_late;
    all_component_selection_stats_late{subj}.candidate_table = candidate_table_late;
    all_component_selection_stats_late{subj}.pf_benchmark = pf_benchmark.late;

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
        fig_save_dir_emg_exclusion, subjects{subj}, 'full', scan_freqs, searchTopos_full, ...
        searchMeanPrSpectrum_full, evals_sorted_full(1:numel(candidate_table_full.comp_idx)), ...
        searchEmgClass_full, ...
        candidate_table_full.hard_eligible, candidate_table_full.force_include_occipital, ...
        rejection_flags_full, ...
        candidate_table_full.front_leak, candidate_table_full.temp_leak, ...
        candidate_table_full.lineharm_ratio, candidate_table_full.hf_slope, ...
        cfg_topo, all_topo_labels{subj}, candidate_table_full.peak_form_score, ...
        candidate_table_full.peak_form_mode, candidate_table_full.peak_form_best_center_hz, ...
        candidate_table_full.peak_form_dominant_penalty, ...
        selected_idx_full, get_candidate_table_fallback_idx(candidate_table_full));
    plot_emg_exclusion_diagnostics( ...
        fig_save_dir_emg_exclusion, subjects{subj}, 'early', scan_freqs, searchTopos_early, ...
        searchMeanPrSpectrum_early, evals_sorted_early(1:numel(candidate_table_early.comp_idx)), ...
        searchEmgClass_early, ...
        candidate_table_early.hard_eligible, candidate_table_early.force_include_occipital, ...
        rejection_flags_early, ...
        candidate_table_early.front_leak, candidate_table_early.temp_leak, ...
        candidate_table_early.lineharm_ratio, candidate_table_early.hf_slope, ...
        cfg_topo, all_topo_labels{subj}, candidate_table_early.peak_form_score, ...
        candidate_table_early.peak_form_mode, candidate_table_early.peak_form_best_center_hz, ...
        candidate_table_early.peak_form_dominant_penalty, ...
        selected_idx_early, get_candidate_table_fallback_idx(candidate_table_early));
    plot_emg_exclusion_diagnostics( ...
        fig_save_dir_emg_exclusion, subjects{subj}, 'late', scan_freqs, searchTopos_late, ...
        searchMeanPrSpectrum_late, evals_sorted_late(1:numel(candidate_table_late.comp_idx)), ...
        searchEmgClass_late, ...
        candidate_table_late.hard_eligible, candidate_table_late.force_include_occipital, ...
        rejection_flags_late, ...
        candidate_table_late.front_leak, candidate_table_late.temp_leak, ...
        candidate_table_late.lineharm_ratio, candidate_table_late.hf_slope, ...
        cfg_topo, all_topo_labels{subj}, candidate_table_late.peak_form_score, ...
        candidate_table_late.peak_form_mode, candidate_table_late.peak_form_best_center_hz, ...
        candidate_table_late.peak_form_dominant_penalty, ...
        selected_idx_late, get_candidate_table_fallback_idx(candidate_table_late));
    plot_combined_topo_spectra_windows( ...
        fig_save_dir_emg_exclusion, subjects{subj}, scan_freqs, cfg_topo, all_topo_labels{subj}, ...
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
        msg = sprintf(['Subject %s excluded for FULL window downstream metrics: ', ...
            'no adequate GED components. FULL-window condition metrics set to NaN.'], subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, ...
            'WINDOW_EXCLUDED_NO_ADEQUATE_COMPONENTS', msg, struct('window', 'full'));
        all_selected_comp_idx(subj) = NaN;
        all_selected_comp_corr(subj) = NaN;
        all_selected_comp_eval(subj) = NaN;
        all_eigenvalues(subj) = NaN;
        all_selected_comp_indices_multi{subj} = NaN;
        all_selected_comp_weights{subj} = NaN;
    end
    if ~adequate_early
        msg = sprintf(['Subject %s excluded for EARLY window downstream metrics: ', ...
            'no adequate GED components. EARLY-window condition metrics set to NaN.'], subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, ...
            'WINDOW_EXCLUDED_NO_ADEQUATE_COMPONENTS', msg, struct('window', 'early'));
    end
    if ~adequate_late
        msg = sprintf(['Subject %s excluded for LATE window downstream metrics: ', ...
            'no adequate GED components. LATE-window condition metrics set to NaN.'], subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, ...
            'WINDOW_EXCLUDED_NO_ADEQUATE_COMPONENTS', msg, struct('window', 'late'));
    end

    %% Per-condition trial-level spectral scanning
    % Use window-specific filters for each dB-spectrum output
    for wi = 1:3
        if wi == 1
            W_comb = W_combined_full;
            w_comb = w_combined_full;
            W_t = W_top_full;
            sel_idx = selected_idx_full;
            window_adequate = adequate_full;
        elseif wi == 2
            W_comb = W_combined_early;
            w_comb = w_combined_early;
            W_t = W_top_early;
            sel_idx = selected_idx_early;
            window_adequate = adequate_early;
        else
            W_comb = W_combined_late;
            w_comb = w_combined_late;
            W_t = W_top_late;
            sel_idx = selected_idx_late;
            window_adequate = adequate_late;
        end
        if ~window_adequate
            W_comb = zeros(nChans, 0);
            w_comb = [];
            W_t = zeros(nChans, 0);
            sel_idx = [];
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
        W_t = normalize_filters_to_noise_metric(W_t, covBase_full);
        W_comb = normalize_filters_to_noise_metric(W_comb, covBase_full);
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
    % Per-window filters struct for trial-level scanning
    filters = struct('full', struct(), 'early', struct(), 'late', struct());
    filters.full.searchFilters = normalize_filters_to_noise_metric(searchFilters_full, covBase_full);
    filters.full.W_top = W_top_full_norm;
    filters.full.W_combined = W_combined_full_norm;
    filters.full.selected_idx = selected_idx_full_norm;
    filters.full.w_combined = w_combined_full_norm;
    filters.early.searchFilters = normalize_filters_to_noise_metric(searchFilters_early, covBase_full);
    filters.early.W_top = W_top_early_norm;
    filters.early.W_combined = W_combined_early_norm;
    filters.early.selected_idx = selected_idx_early_norm;
    filters.early.w_combined = w_combined_early_norm;
    filters.late.searchFilters = normalize_filters_to_noise_metric(searchFilters_late, covBase_full);
    filters.late.W_top = W_top_late_norm;
    filters.late.W_combined = W_combined_late_norm;
    filters.late.selected_idx = selected_idx_late_norm;
    filters.late.w_combined = w_combined_late_norm;

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

        bad_base_raw = flag_unreliable_baseline_trials( ...
            baseline_power_raw, baseline_outlier_mad_mult, baseline_outlier_min_trials, baseline_outlier_exclusion_enable);
        bad_base_full = bad_base_raw;
        bad_base_early = bad_base_raw;
        bad_base_late = bad_base_raw;
        if adequate_full
            bad_base_full = bad_base_full | flag_unreliable_baseline_trials( ...
                baseline_power_comb_full, baseline_outlier_mad_mult, baseline_outlier_min_trials, baseline_outlier_exclusion_enable);
        end
        if adequate_early
            bad_base_early = bad_base_early | flag_unreliable_baseline_trials( ...
                baseline_power_comb_early, baseline_outlier_mad_mult, baseline_outlier_min_trials, baseline_outlier_exclusion_enable);
        end
        if adequate_late
            bad_base_late = bad_base_late | flag_unreliable_baseline_trials( ...
                baseline_power_comb_late, baseline_outlier_mad_mult, baseline_outlier_min_trials, baseline_outlier_exclusion_enable);
        end
        [base_floor_raw, ~] = compute_baseline_floor_stats(baseline_power_raw, ratio_floor_prctile, ratio_floor_frac);
        [base_floor_full, ~] = compute_baseline_floor_stats(baseline_power_comb_full, ratio_floor_prctile, ratio_floor_frac);
        [base_floor_early, ~] = compute_baseline_floor_stats(baseline_power_comb_early, ratio_floor_prctile, ratio_floor_frac);
        [base_floor_late, ~] = compute_baseline_floor_stats(baseline_power_comb_late, ratio_floor_prctile, ratio_floor_frac);

        for trl = 1:nTrl
            tc = trial_cache{trl};
            x_base = tc.x_base;
            x_full = tc.x_full;
            x_early = tc.x_early;
            x_late = tc.x_late;
            if isempty(x_base)
                continue;
            end

            if adequate_full && ~bad_base_full(trl)
                comp_base = filters.full.searchFilters' * x_base;
                comp_stim = filters.full.searchFilters' * x_full;
                [ratio_mat_full, near_floor_mask_full] = compute_scan_ratio_from_timeseries( ...
                    comp_stim, comp_base, fsample, scan_freqs, scan_width, base_floor_full, instability_near_floor_mult);
                powratio_components(:, trl, :) = ratio_mat_full;
                if ~isempty(ratio_mat_full)
                    if ~isempty(filters.full.selected_idx)
                        w_use = filters.full.w_combined(:);
                        if numel(w_use) ~= numel(filters.full.selected_idx)
                            w_use = ones(numel(filters.full.selected_idx), 1) / max(numel(filters.full.selected_idx), 1);
                        end
                        w_use = w_use / max(sum(w_use), eps);
                        sel_ratio = ratio_mat_full(filters.full.selected_idx, :);
                        powratio_methods_full(1, trl, :) = (w_use' * sel_ratio);
                    end
                end
                valid_freq_counts_full(trl) = sum(any(isfinite(ratio_mat_full), 1));
                if any(any(isfinite(ratio_mat_full), 1))
                    unstable_freq_counts_full(trl) = sum(near_floor_mask_full);
                end
            end

            if adequate_early && ~bad_base_early(trl)
                comp_base = filters.early.searchFilters' * x_base;
                comp_stim = filters.early.searchFilters' * x_early;
                [ratio_mat_early, near_floor_mask_early] = compute_scan_ratio_from_timeseries( ...
                    comp_stim, comp_base, fsample, scan_freqs, scan_width, base_floor_early, instability_near_floor_mult);
                powratio_components_early(:, trl, :) = ratio_mat_early;
                if ~isempty(ratio_mat_early)
                    if ~isempty(filters.early.selected_idx)
                        w_use = filters.early.w_combined(:);
                        if numel(w_use) ~= numel(filters.early.selected_idx)
                            w_use = ones(numel(filters.early.selected_idx), 1) / max(numel(filters.early.selected_idx), 1);
                        end
                        w_use = w_use / max(sum(w_use), eps);
                        sel_ratio = ratio_mat_early(filters.early.selected_idx, :);
                        powratio_methods_early(1, trl, :) = (w_use' * sel_ratio);
                    end
                end
                valid_freq_counts_early(trl) = sum(any(isfinite(ratio_mat_early), 1));
                if any(any(isfinite(ratio_mat_early), 1))
                    unstable_freq_counts_early(trl) = sum(near_floor_mask_early);
                end
            end

            if adequate_late && ~bad_base_late(trl)
                comp_base = filters.late.searchFilters' * x_base;
                comp_stim = filters.late.searchFilters' * x_late;
                [ratio_mat_late, near_floor_mask_late] = compute_scan_ratio_from_timeseries( ...
                    comp_stim, comp_base, fsample, scan_freqs, scan_width, base_floor_late, instability_near_floor_mult);
                powratio_components_late(:, trl, :) = ratio_mat_late;
                if ~isempty(ratio_mat_late)
                    if ~isempty(filters.late.selected_idx)
                        w_use = filters.late.w_combined(:);
                        if numel(w_use) ~= numel(filters.late.selected_idx)
                            w_use = ones(numel(filters.late.selected_idx), 1) / max(numel(filters.late.selected_idx), 1);
                        end
                        w_use = w_use / max(sum(w_use), eps);
                        sel_ratio = ratio_mat_late(filters.late.selected_idx, :);
                        powratio_methods_late(1, trl, :) = (w_use' * sel_ratio);
                    end
                end
                valid_freq_counts_late(trl) = sum(any(isfinite(ratio_mat_late), 1));
                if any(any(isfinite(ratio_mat_late), 1))
                    unstable_freq_counts_late(trl) = sum(near_floor_mask_late);
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
        if any(trial_unstable_full) || any(trial_unstable_early) || any(trial_unstable_late)
            msg = sprintf(['Numerical instability exclusion for subject %s condition %s: ', ...
                'full=%d/%d, early=%d/%d, late=%d/%d trials removed (near-floor baseline).'], ...
                subjects{subj}, condLabels{cond}, ...
                sum(trial_unstable_full), nTrl, sum(trial_unstable_early), nTrl, sum(trial_unstable_late), nTrl);
            warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, ...
                'NUMERICAL_INSTABILITY_TRIAL_EXCLUSION', msg, struct( ...
                'condition', condLabels{cond}, ...
                'n_removed_full', sum(trial_unstable_full), ...
                'n_removed_early', sum(trial_unstable_early), ...
                'n_removed_late', sum(trial_unstable_late), ...
                'threshold_trial_freq_fraction', instability_trial_freq_frac_thr));
        end
        powratio_methods_full_analysis = powratio_methods_full;
        powratio_methods_early_analysis = powratio_methods_early;
        powratio_methods_late_analysis = powratio_methods_late;
        powratio_components_analysis = powratio_components;
        powratio_components_early_analysis = powratio_components_early;
        powratio_components_late_analysis = powratio_components_late;

        all_trial_powratio_components_full{cond, subj}  = powratio_components_analysis;
        all_trial_powratio_components_early{cond, subj} = powratio_components_early_analysis;
        all_trial_powratio_components_late{cond, subj}  = powratio_components_late_analysis;

        % Keep full-window outputs based on weighted combined GED branch.
        powratio_trials_fullscan = squeeze(powratio_methods_full(1, :, :));
        powratio_trials_early_fullscan = squeeze(powratio_methods_early(1, :, :));
        powratio_trials_late_fullscan = squeeze(powratio_methods_late(1, :, :));
        powratio_trials_full = squeeze(powratio_methods_full_analysis(1, :, :));
        powratio_trials_early = squeeze(powratio_methods_early_analysis(1, :, :));
        powratio_trials_late = squeeze(powratio_methods_late_analysis(1, :, :));
        powratio_trials_fullscan_plotstat = powratio_trials_fullscan;
        powratio_trials_full_plotstat = powratio_trials_full;
        powratio_trials_early_plotstat = powratio_trials_early;
        powratio_trials_late_plotstat = powratio_trials_late;
        all_trial_powratio_fullscan{cond, subj} = powratio_trials_fullscan;
        all_trial_powratio{cond, subj} = powratio_trials_full;
        all_trial_powratio_early{cond, subj} = powratio_trials_early;
        all_trial_powratio_late{cond, subj} = powratio_trials_late;
        all_trial_powratio_fullscan_plotstat{cond, subj} = powratio_trials_fullscan_plotstat;
        all_trial_powratio_plotstat{cond, subj} = powratio_trials_full_plotstat;
        all_trial_powratio_early_plotstat{cond, subj} = powratio_trials_early_plotstat;
        all_trial_powratio_late_plotstat{cond, subj} = powratio_trials_late_plotstat;
        %% Per-trial peak detection
        [trl_peaks_single, trial_peak_power_full, trl_peaks_low, trl_peaks_high, trl_centroid] = ...
            compute_trial_peak_metrics_from_powratio_fullscan( ...
            powratio_trials_fullscan, scan_freqs, true(size(scan_freqs)), ...
            centroid_band_mask, trial_peak_smooth_n, centroid_min_peak, centroid_posfrac_min);

        all_trial_peaks_single{cond, subj} = trl_peaks_single;
        all_trial_peaks_low{cond, subj}    = trl_peaks_low;
        all_trial_peaks_high{cond, subj}   = trl_peaks_high;
        all_trial_centroid{cond, subj}     = trl_centroid;

        valid_s = ~isnan(trl_peaks_single);
        all_trial_mean_single(cond, subj)   = robust_trial_mean(trl_peaks_single(valid_s));
        all_trial_median_single(cond, subj) = median(trl_peaks_single(valid_s));
        all_trial_detrate_single(cond, subj) = sum(valid_s) / nTrl;

        valid_lo = ~isnan(trl_peaks_low);
        all_trial_mean_low(cond, subj)   = robust_trial_mean(trl_peaks_low(valid_lo));
        all_trial_median_low(cond, subj) = median(trl_peaks_low(valid_lo));
        all_trial_detrate_low(cond, subj) = sum(valid_lo) / nTrl;

        valid_hi = ~isnan(trl_peaks_high);
        all_trial_mean_high(cond, subj)   = robust_trial_mean(trl_peaks_high(valid_hi));
        all_trial_median_high(cond, subj) = median(trl_peaks_high(valid_hi));
        all_trial_detrate_high(cond, subj) = sum(valid_hi) / nTrl;
        valid_gap = valid_lo & valid_hi;
        if any(valid_gap)
            all_trial_median_gap(cond, subj) = median(trl_peaks_high(valid_gap) - trl_peaks_low(valid_gap));
        end

        valid_c = ~isnan(trl_centroid);
        all_trial_mean_centroid(cond, subj)   = robust_trial_mean(trl_centroid(valid_c));
        all_trial_median_centroid(cond, subj) = median(trl_centroid(valid_c));
        all_trial_detrate_centroid(cond, subj) = sum(valid_c) / nTrl;

        % Time-split peak summaries.
        [trl_peaks_single_early, trial_peak_power_early, trl_peaks_low_early, trl_peaks_high_early, trl_centroid_early] = ...
            compute_trial_peak_metrics_from_powratio_fullscan( ...
            powratio_trials_early_fullscan, scan_freqs, true(size(scan_freqs)), ...
            centroid_band_mask, trial_peak_smooth_n, centroid_min_peak, centroid_posfrac_min);
        [trl_peaks_single_late, trial_peak_power_late, trl_peaks_low_late, trl_peaks_high_late, trl_centroid_late] = ...
            compute_trial_peak_metrics_from_powratio_fullscan( ...
            powratio_trials_late_fullscan, scan_freqs, true(size(scan_freqs)), ...
            centroid_band_mask, trial_peak_smooth_n, centroid_min_peak, centroid_posfrac_min);
        % Trial-level metric outlier rejection (subject-condition specific).
        if trial_metric_outlier_enable
            [outlier_mask_freq_full, ~] = detect_trial_metric_outliers_iqr( ...
                trl_peaks_single, trial_metric_outlier_iqr_mult, trial_metric_outlier_min_trials);
            [outlier_mask_freq_early, ~] = detect_trial_metric_outliers_iqr( ...
                trl_peaks_single_early, trial_metric_outlier_iqr_mult, trial_metric_outlier_min_trials);
            [outlier_mask_freq_late, ~] = detect_trial_metric_outliers_iqr( ...
                trl_peaks_single_late, trial_metric_outlier_iqr_mult, trial_metric_outlier_min_trials);
            [outlier_mask_power_full, ~] = detect_trial_metric_outliers_iqr( ...
                trial_peak_power_full, trial_metric_outlier_iqr_mult, trial_metric_outlier_min_trials);
            [outlier_mask_power_early, ~] = detect_trial_metric_outliers_iqr( ...
                trial_peak_power_early, trial_metric_outlier_iqr_mult, trial_metric_outlier_min_trials);
            [outlier_mask_power_late, ~] = detect_trial_metric_outliers_iqr( ...
                trial_peak_power_late, trial_metric_outlier_iqr_mult, trial_metric_outlier_min_trials);
        else
            outlier_mask_freq_full = false(size(trl_peaks_single));
            outlier_mask_freq_early = false(size(trl_peaks_single_early));
            outlier_mask_freq_late = false(size(trl_peaks_single_late));
            outlier_mask_power_full = false(size(trial_peak_power_full));
            outlier_mask_power_early = false(size(trial_peak_power_early));
            outlier_mask_power_late = false(size(trial_peak_power_late));
        end
        trl_peaks_single(outlier_mask_freq_full) = NaN;
        trl_peaks_single_early(outlier_mask_freq_early) = NaN;
        trl_peaks_single_late(outlier_mask_freq_late) = NaN;
        trial_peak_power_full(outlier_mask_power_full) = NaN;
        trial_peak_power_early(outlier_mask_power_early) = NaN;
        trial_peak_power_late(outlier_mask_power_late) = NaN;
        all_trial_peaks_single{cond, subj} = trl_peaks_single;
        all_trial_outlier_mask_freq_full{cond, subj} = outlier_mask_freq_full;
        all_trial_outlier_mask_freq_early{cond, subj} = outlier_mask_freq_early;
        all_trial_outlier_mask_freq_late{cond, subj} = outlier_mask_freq_late;
        all_trial_outlier_mask_power_full{cond, subj} = outlier_mask_power_full;
        all_trial_outlier_mask_power_early{cond, subj} = outlier_mask_power_early;
        all_trial_outlier_mask_power_late{cond, subj} = outlier_mask_power_late;
        all_trial_peaks_single_early{cond, subj} = trl_peaks_single_early;
        all_trial_peaks_single_late{cond, subj} = trl_peaks_single_late;

        valid_s_early = ~isnan(trl_peaks_single_early);
        all_trial_mean_single_early(cond, subj) = robust_trial_mean(trl_peaks_single_early(valid_s_early));
        all_trial_median_single_early(cond, subj) = median(trl_peaks_single_early(valid_s_early));
        all_trial_detrate_single_early(cond, subj) = sum(valid_s_early) / nTrl;

        valid_s_late = ~isnan(trl_peaks_single_late);
        all_trial_mean_single_late(cond, subj) = robust_trial_mean(trl_peaks_single_late(valid_s_late));
        all_trial_median_single_late(cond, subj) = median(trl_peaks_single_late(valid_s_late));
        all_trial_detrate_single_late(cond, subj) = sum(valid_s_late) / nTrl;

        % Peak power: highest dB value in the smoothed trial spectrum.
        all_trial_gamma_power(cond, subj) = robust_trial_mean(trial_peak_power_full);
        all_trial_gamma_power_early(cond, subj) = robust_trial_mean(trial_peak_power_early);
        all_trial_gamma_power_late(cond, subj) = robust_trial_mean(trial_peak_power_late);
        all_trial_gamma_power_plotstat(cond, subj) = robust_trial_mean(trial_peak_power_full);
        all_trial_gamma_power_early_plotstat(cond, subj) = robust_trial_mean(trial_peak_power_early);
        all_trial_gamma_power_late_plotstat(cond, subj) = robust_trial_mean(trial_peak_power_late);

        % Time-split dual-peak, centroid, and reliability summaries.
        valid_lo_early = ~isnan(trl_peaks_low_early);
        valid_hi_early = ~isnan(trl_peaks_high_early);
        if any(valid_lo_early)
            all_trial_median_low_early(cond, subj) = median(trl_peaks_low_early(valid_lo_early));
        end
        if any(valid_hi_early)
            all_trial_median_high_early(cond, subj) = median(trl_peaks_high_early(valid_hi_early));
        end
        valid_gap_early = valid_lo_early & valid_hi_early;
        if any(valid_gap_early)
            all_trial_median_gap_early(cond, subj) = median(trl_peaks_high_early(valid_gap_early) - trl_peaks_low_early(valid_gap_early));
        end
        valid_cent_early = isfinite(trl_centroid_early);
        if any(valid_cent_early)
            all_trial_median_centroid_early(cond, subj) = median(trl_centroid_early(valid_cent_early));
        end
        if sum(valid_s_early) >= 2
            vf_early = trl_peaks_single_early(valid_s_early);
            all_trial_trialcv_early(cond, subj) = std(vf_early) / abs(mean(vf_early));
        end

        valid_lo_late = ~isnan(trl_peaks_low_late);
        valid_hi_late = ~isnan(trl_peaks_high_late);
        if any(valid_lo_late)
            all_trial_median_low_late(cond, subj) = median(trl_peaks_low_late(valid_lo_late));
        end
        if any(valid_hi_late)
            all_trial_median_high_late(cond, subj) = median(trl_peaks_high_late(valid_hi_late));
        end
        valid_gap_late = valid_lo_late & valid_hi_late;
        if any(valid_gap_late)
            all_trial_median_gap_late(cond, subj) = median(trl_peaks_high_late(valid_gap_late) - trl_peaks_low_late(valid_gap_late));
        end
        valid_cent_late = isfinite(trl_centroid_late);
        if any(valid_cent_late)
            all_trial_median_centroid_late(cond, subj) = median(trl_centroid_late(valid_cent_late));
        end
        if sum(valid_s_late) >= 2
            vf_late = trl_peaks_single_late(valid_s_late);
            all_trial_trialcv_late(cond, subj) = std(vf_late) / abs(mean(vf_late));
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
            pr_source = all_trial_powratio_fullscan;
            peaks_single_source = all_trial_peaks_single;
            median_low_source = all_trial_median_low;
            median_high_source = all_trial_median_high;
            detrate_low_source = all_trial_detrate_low;
            detrate_high_source = all_trial_detrate_high;
            centroid_source = all_trial_centroid;
            topo_mat_window = searchTopos_full;
            selected_idx_window = selected_idx_full;
            selected_w_window = w_combined_full;
            eigvals_window = evals_sorted_full;
        elseif wi == 2
            pr_source = all_trial_powratio_early;
            peaks_single_source = all_trial_peaks_single_early;
            median_low_source = all_trial_median_low_early;
            median_high_source = all_trial_median_high_early;
            detrate_low_source = all_trial_detrate_single_early;
            detrate_high_source = all_trial_detrate_single_early;
            centroid_source = cell(size(all_trial_powratio_early));
            topo_mat_window = searchTopos_early;
            selected_idx_window = selected_idx_early;
            selected_w_window = w_combined_early;
            eigvals_window = evals_sorted_early;
        else
            pr_source = all_trial_powratio_late;
            peaks_single_source = all_trial_peaks_single_late;
            median_low_source = all_trial_median_low_late;
            median_high_source = all_trial_median_high_late;
            detrate_low_source = all_trial_detrate_single_late;
            detrate_high_source = all_trial_detrate_single_late;
            centroid_source = cell(size(all_trial_powratio_late));
            topo_mat_window = searchTopos_late;
            selected_idx_window = selected_idx_late;
            selected_w_window = w_combined_late;
            eigvals_window = evals_sorted_late;
        end

        pr_raw_mats = cell(1, 4);
        row1_clim = ones(1, 4);
        for cond = 1:4
            pr_raw_mats{cond} = pr_source{cond, subj};
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
            subplot(3, 4, cond);
            if ~isempty(pr_raw_mats{cond})
                hold on;
                imagesc(scan_freqs, 1:size(pr_raw_mats{cond},1), pr_raw_mats{cond});
                colormap(gca, cmap_div);
                caxis([-row1_clim(cond) row1_clim(cond)]);
                cb = colorbar; cb.FontSize = 8;
                ctd = centroid_source{cond, subj};
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
            title(sprintf('%s Raw', condLabels{cond}), 'FontSize', 11);
            set(gca, 'FontSize', 10); xlim([30 90]); box on;
        end

        % --- Row 2: Mean trial-level spectrum with dual-peak markers ---
        for cond = 1:4
            subplot(3, 4, 4 + cond); hold on;
            if ~isempty(pr_raw_mats{cond})
                mu_raw = nanmean(pr_raw_mats{cond}, 1);
                nTrl = size(pr_raw_mats{cond}, 1);
                sem_raw = nanstd(pr_raw_mats{cond}, [], 1) / sqrt(max(1, nTrl));
                row2_min = min(mu_raw - sem_raw, [], 'omitnan');
                row2_max = max(mu_raw + sem_raw, [], 'omitnan');
                row2_abs = max(abs([mu_raw - sem_raw, mu_raw + sem_raw]), [], 'omitnan');
                if ~isfinite(row2_abs) || row2_abs <= 0
                    row2_abs = 1;
                end
                if ~isfinite(row2_min)
                    row2_min = 0;
                end
                if ~isfinite(row2_max)
                    row2_max = 2;
                end
                faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
                patch([scan_freqs, fliplr(scan_freqs)], ...
                    [mu_raw - sem_raw, fliplr(mu_raw + sem_raw)], ...
                    colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
                plot(scan_freqs, mu_raw, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);
                yline(0, 'k-', 'LineWidth', 0.5);
                xline(50, 'k:', 'LineWidth', 1, 'Alpha', 0.5);
                pf_lo = median_low_source(cond, subj);
                pf_hi = median_high_source(cond, subj);
                if ~isnan(pf_lo)
                    xline(pf_lo, '--', 'LineWidth', 2, 'Color', [0 0 0.7]);
                    text(pf_lo + 1, max(mu_raw) * 0.85, ...
                        sprintf('L:%.0f', pf_lo), 'FontSize', 9, ...
                        'Color', [0 0 0.7], 'FontWeight', 'bold');
                end
                if ~isnan(pf_hi)
                    xline(pf_hi, '--', 'LineWidth', 2, 'Color', [0.7 0 0]);
                    text(pf_hi + 1, max(mu_raw) * 0.65, ...
                        sprintf('H:%.0f', pf_hi), 'FontSize', 9, ...
                        'Color', [0.7 0 0], 'FontWeight', 'bold');
                end
                y_lo = row2_min - 0.08 * abs(row2_max - row2_min);
                y_hi = row2_max + 0.08 * abs(row2_max - row2_min);
                if ~isfinite(y_lo) || ~isfinite(y_hi) || y_hi <= y_lo
                    y_lo = -10;
                    y_hi = 10;
                end
                ylim([y_lo y_hi]);
            end
            xlabel('Freq [Hz]');
            ylabel('Power [dB]');
            det_lo = detrate_low_source(cond, subj);
            det_hi = detrate_high_source(cond, subj);
            title(sprintf('%s Dual (L:%.0f%% H:%.0f%%)', condLabels{cond}, det_lo*100, det_hi*100), 'FontSize', 10);
            set(gca, 'FontSize', 10); xlim([30 90]);  box on;
        end

        % --- Row 3: Topoplot + histogram ---
        subplot(3, 4, 9);
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

        subplot(3, 4, [10 11 12]); hold on;
        edges = 30:2:90;
        hist_mat = zeros(4, length(edges)-1);
        for cond = 1:4
            tpk = peaks_single_source{cond, subj};
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

        save_figure_png(fig, fullfile(comp_sel_save_dir, ...
            sprintf('GCP_eeg_GED_subj%s_%s.png', subjects{subj}, window_names{wi})));
    end
    warning_log_by_subj{subj} = warning_log_subj;
    toc
end % subject loop

warning_log_cells = warning_log_by_subj;
for subj = 1:numel(warning_log_cells)
    warning_log_cells{subj} = warning_log_cells{subj}(:);
end
warning_log_cells = warning_log_cells(~cellfun(@isempty, warning_log_cells));
if isempty(warning_log_cells)
    warning_log = struct('subject', {}, 'code', {}, 'message', {}, 'metrics', {});
else
    warning_log = vertcat(warning_log_cells{:});
end
print_subject_warning_summary(warning_log);

%% Combined GED metrics (full window only)
benchmark_metric_detectability = nan(4, nSubj);
benchmark_metric_separation_slope = nan(nSubj, 1);
benchmark_metric_separation_delta = nan(nSubj, 1);
benchmark_metric_reliability_trialcv = nan(4, nSubj);
benchmark_metric_reliability_subjspread = nan(4, 1);

for subj = 1:nSubj
    cond_medians = nan(4, 1);
    for cond = 1:4
        peak_freq = all_trial_peaks_single{cond, subj};
        if isempty(peak_freq)
            continue;
        end
        benchmark_metric_detectability(cond, subj) = mean(isfinite(peak_freq));
        vf = peak_freq(isfinite(peak_freq));
        if numel(vf) >= 2 && mean(vf) ~= 0
            benchmark_metric_reliability_trialcv(cond, subj) = std(vf) / abs(mean(vf));
        end
        cond_medians(cond) = median(vf);
    end

    vx = ~isnan(cond_medians);
    if sum(vx) >= 2
        p = polyfit(find(vx), cond_medians(vx)', 1);
        benchmark_metric_separation_slope(subj) = p(1);
    end
    if ~isnan(cond_medians(1)) && ~isnan(cond_medians(4))
        benchmark_metric_separation_delta(subj) = cond_medians(4) - cond_medians(1);
    end
end

for cond = 1:4
    v = benchmark_metric_reliability_trialcv(cond, :);
    benchmark_metric_reliability_subjspread(cond, 1) = std(v(~isnan(v)));
end

% CENTROID METRIC: Subject/group summaries and concordance
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
save_figure_png(fig_cent, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_centroid_summary.png'));

%% Standalone condition-separation metrics (combined GED)
close all
fig_cond_slope = figure('Position', [0 0 1512 982], 'Color', 'w');

slope_post = benchmark_metric_separation_slope;
delta_post = benchmark_metric_separation_delta;
primary_slope_stats = compute_one_sample_stats(slope_post);
primary_delta_stats = compute_one_sample_stats(delta_post);

fprintf('Primary inference (combined GED, subject-level)\n');
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

save_figure_png(fig_cond_slope, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_condition_slope.png'));

%% Mean gamma frequency shift bar plot
close all
fig_cond_shift_bar = figure('Position', [0 0 1512 982], 'Color', 'w');
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
ylabel('\Delta Frequency Shift (100% - 25%) [Hz]', 'FontSize', 18, 'FontWeight', 'bold');
title('Gamma Frequency Shift', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
save_figure_png(fig_cond_shift_bar, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_bar_GammaFreq.png'));

%% Summary dashboard (backprojected combined-component data)
fig_summary = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('GED Trials Summary Dashboard (component selection backprojected)', ...
    'FontSize', 16, 'FontWeight', 'bold');
gamma_metric_full = all_trial_gamma_power;
gamma_metric_label = 'Peak Power [dB]';

summary_metrics = { ...
    all_trial_median_single, ...
    all_trial_median_low, ...
    all_trial_median_high, ...
    all_trial_median_gap, ...
    benchmark_metric_reliability_trialcv, ...
    all_trial_median_centroid, ...
    gamma_metric_full};
summary_names = {'Single median [Hz]', 'Low median [Hz]', 'High median [Hz]', ...
    'H-L separation [Hz]', 'Trial CV', 'Centroid median [Hz]', gamma_metric_label};

for mi = 1:numel(summary_metrics)
    subplot(2, 4, mi); hold on;
    dat = summary_metrics{mi};
    dat_plot = dat;
    % Figure-only outlier exclusion per metric and condition (subject-level).
    for c = 1:4
        cond_vals = dat_plot(c, :);
        valid_idx = find(isfinite(cond_vals));
        if numel(valid_idx) >= 4
            outlier_cond = isoutlier(cond_vals(valid_idx), 'median');
            dat_plot(c, valid_idx(outlier_cond)) = NaN;
        end
    end
    mu = nanmean(dat_plot, 2);
    sem = nanstd(dat_plot, [], 2) ./ sqrt(sum(~isnan(dat_plot), 2));
    med = nanmedian(dat_plot, 2);
    for c = 1:4
        bar(c, mu(c), 0.6, 'FaceColor', colors(c,:), 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    end
    errorbar(1:4, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.1, 'CapSize', 5);
    scatter(1:4, med, 28, 'kd', 'filled');
    for s = 1:nSubj
        for c = 1:4
            if ~isnan(dat_plot(c, s))
                scatter(c + (rand - 0.5) * 0.18, dat_plot(c, s), 20, [0.5 0.5 0.5], ...
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
    elseif mi == 6
        ylim([50 60]);
    elseif mi == 7
        ylim([-15 40]);
    end
end
apply_dynamic_summary_ylims();
save_figure_png(fig_summary, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_metrics_summary_full.png'));

fig_summary_early = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('GED Trials Summary Dashboard (component selection backprojected, early window)', ...
    'FontSize', 16, 'FontWeight', 'bold');
gamma_metric_early = all_trial_gamma_power_early;
summary_metrics_early = { ...
    all_trial_median_single_early, ...
    all_trial_median_low_early, ...
    all_trial_median_high_early, ...
    all_trial_median_gap_early, ...
    all_trial_trialcv_early, ...
    all_trial_median_centroid_early, ...
    gamma_metric_early};
for mi = 1:numel(summary_metrics_early)
    subplot(2, 4, mi); hold on;
    dat = summary_metrics_early{mi};
    dat_plot = dat;
    for c = 1:4
        cond_vals = dat_plot(c, :);
        valid_idx = find(isfinite(cond_vals));
        if numel(valid_idx) >= 4
            outlier_cond = isoutlier(cond_vals(valid_idx), 'median');
            dat_plot(c, valid_idx(outlier_cond)) = NaN;
        end
    end
    mu = nanmean(dat_plot, 2);
    sem = nanstd(dat_plot, [], 2) ./ sqrt(sum(~isnan(dat_plot), 2));
    med = nanmedian(dat_plot, 2);
    for c = 1:4
        bar(c, mu(c), 0.6, 'FaceColor', colors(c,:), 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    end
    errorbar(1:4, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.1, 'CapSize', 5);
    scatter(1:4, med, 28, 'kd', 'filled');
    for s = 1:nSubj
        for c = 1:4
            if ~isnan(dat_plot(c, s))
                scatter(c + (rand - 0.5) * 0.18, dat_plot(c, s), 20, [0.5 0.5 0.5], ...
                    'filled', ...
                    'MarkerFaceAlpha', 0.5, ...
                    'MarkerEdgeColor', [1 1 1], ...
                    'LineWidth', 0.5, ...
                    'MarkerEdgeAlpha', 0.9);
            end
        end
    end
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 9);
    title(summary_names{mi}, 'FontSize', 10, 'FontWeight', 'bold');
    box on;
end
apply_dynamic_summary_ylims();
save_figure_png(fig_summary_early, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_metrics_summary_early.png'));

fig_summary_late = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('GED Trials Summary Dashboard (component selection backprojected, late window)', ...
    'FontSize', 16, 'FontWeight', 'bold');
gamma_metric_late = all_trial_gamma_power_late;
summary_metrics_late = { ...
    all_trial_median_single_late, ...
    all_trial_median_low_late, ...
    all_trial_median_high_late, ...
    all_trial_median_gap_late, ...
    all_trial_trialcv_late, ...
    all_trial_median_centroid_late, ...
    gamma_metric_late};
for mi = 1:numel(summary_metrics_late)
    subplot(2, 4, mi); hold on;
    dat = summary_metrics_late{mi};
    dat_plot = dat;
    for c = 1:4
        cond_vals = dat_plot(c, :);
        valid_idx = find(isfinite(cond_vals));
        if numel(valid_idx) >= 4
            outlier_cond = isoutlier(cond_vals(valid_idx), 'median');
            dat_plot(c, valid_idx(outlier_cond)) = NaN;
        end
    end
    mu = nanmean(dat_plot, 2);
    sem = nanstd(dat_plot, [], 2) ./ sqrt(sum(~isnan(dat_plot), 2));
    med = nanmedian(dat_plot, 2);
    for c = 1:4
        bar(c, mu(c), 0.6, 'FaceColor', colors(c,:), 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    end
    errorbar(1:4, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.1, 'CapSize', 5);
    scatter(1:4, med, 28, 'kd', 'filled');
    for s = 1:nSubj
        for c = 1:4
            if ~isnan(dat_plot(c, s))
                scatter(c + (rand - 0.5) * 0.18, dat_plot(c, s), 20, [0.5 0.5 0.5], ...
                    'filled', ...
                    'MarkerFaceAlpha', 0.5, ...
                    'MarkerEdgeColor', [1 1 1], ...
                    'LineWidth', 0.5, ...
                    'MarkerEdgeAlpha', 0.9);
            end
        end
    end
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 9);
    title(summary_names{mi}, 'FontSize', 10, 'FontWeight', 'bold');
    box on;
end
apply_dynamic_summary_ylims();
save_figure_png(fig_summary_late, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_metrics_summary_late.png'));

%% Grand average gamma figures
close all
fprintf('\nCreating grand average figures...\n');

%% Grand average: single peak (stats-style boxplot, non-baselined frequency)
fig_box1_statsstyle = figure('Position', [0 0 1512 982], 'Color', 'w');
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

save_figure_png(fig_box1_statsstyle, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_boxplot_GammaFreq.png'));

%% Main figure: gamma frequency over contrast (mean+-SEM and trajectories)
close all
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
save_figure_png(fig_main_gamma, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_main_GammaFreq.png'));

%% Condition-shift figure: normalized gamma frequency trajectories
fig_condition_shift = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

dat_freq = all_trial_mean_single;  % [condition x subject], robust mean over trials
dat_shift = dat_freq - dat_freq(1, :);  % Anchor each subject at 25% condition

for s = 1:nSubj
    y_subj = dat_shift(:, s);
    valid_subj = isfinite(y_subj);
    if sum(valid_subj) >= 2
        plot(find(valid_subj), y_subj(valid_subj), '-o', ...
            'Color', [0.75 0.75 0.75], 'LineWidth', 1.2, ...
            'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerSize', 6);
    end
end

med_shift = nanmedian(dat_shift, 2);
mad_shift = nan(4, 1);
for c = 1:4
    mad_shift(c) = robust_mad(dat_shift(c, :));
end
errorbar(1:4, med_shift, mad_shift, '-o', ...
    'Color', colors(4, :), 'LineWidth', 2.8, 'CapSize', 10, ...
    'MarkerFaceColor', colors(4, :), 'MarkerSize', 8);

yline(0, 'k--', 'LineWidth', 1.5);
set(gca, 'XTick', 1:4, 'XTickLabel', strcat(condLabels, ' Contrast'), ...
    'FontSize', 18, 'Box', 'off');
xlim([0.5 4.5]);
ylim([-6 6]);
xlabel('Contrast condition');
ylabel('\Delta Gamma Frequency [Hz]');
title('Gamma Peak Frequency Shift', ...
    'FontSize', 24, 'FontWeight', 'bold');

cond_shift_path = fullfile(fig_save_dir_ged, 'GCP_eeg_GED_condition_shift.png');
save_figure_png(fig_condition_shift, cond_shift_path);

%% Main figure: gamma frequency over contrast by time window
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
save_figure_png(fig_main_gamma_windows, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_main_GammaFreq_timeSplit.png'));

%% Power figure: peak dB power over conditions
fig_power_statsstyle = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

dat_power = all_trial_gamma_power;
dat_power_plot = dat_power;

for s = 1:nSubj
    y_subj = dat_power_plot(:, s);
    valid_subj = ~isnan(y_subj);
    if sum(valid_subj) >= 2
        plot(find(valid_subj), y_subj(valid_subj), '-', ...
            'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    end
end

y_all = dat_power_plot(:);
g_all = repmat((1:4)', nSubj, 1);
valid_all = ~isnan(y_all);
boxplot(y_all(valid_all), g_all(valid_all), 'Colors', [0.45 0.45 0.45], ...
    'Symbol', '', 'Widths', 0.5);

for c = 1:4
    y_c = dat_power_plot(c, :);
    valid_c = ~isnan(y_c);
    x_jit = c + (rand(1, sum(valid_c)) - 0.5) * 0.10;
    scatter(x_jit, y_c(valid_c), 170, colors(c,:), 'filled', ...
        'MarkerEdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.7);
end
yline(0, 'k--', 'LineWidth', 1.2);
set(gca, 'XTick', 1:4, 'XTickLabel', strcat(condLabels, ' Contrast'), ...
    'FontSize', 15, 'Box', 'off');
xlim([0.5 4.5]);
ylabel('Peak Power [dB]');
title('Peak Power Increase', 'FontSize', 30, 'FontWeight', 'bold');
save_figure_png(fig_power_statsstyle, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_boxplot_power.png'));

%% Condition-shift figure: normalized peak power trajectories
fig_condition_shift_power = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

dat_power_shift = all_trial_gamma_power - all_trial_gamma_power(1, :);  % Anchor each subject at 25% condition

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
    'Color', colors(4, :), 'LineWidth', 2.8, 'CapSize', 10, ...
    'MarkerFaceColor', colors(4, :), 'MarkerSize', 8);

yline(0, 'k--', 'LineWidth', 1.5);
set(gca, 'XTick', 1:4, 'XTickLabel', strcat(condLabels, ' Contrast'), ...
    'FontSize', 18, 'Box', 'off');
xlim([0.5 4.5]);
ylim([-0.75 0.75]);
xlabel('Contrast condition');
ylabel('\Delta Peak Power [dB]');
title('Gamma Peak Power Shift', ...
    'FontSize', 24, 'FontWeight', 'bold');

cond_shift_power_path = fullfile(fig_save_dir_ged, 'GCP_eeg_GED_condition_shift_power.png');
save_figure_png(fig_condition_shift_power, cond_shift_power_path);

%% Save results
save_path = fullfile(gcp_root_path, 'data', 'features', 'GCP_eeg_GED.mat');
save(save_path, ...
    'all_trial_powratio', ...
    'all_trial_powratio_plotstat', ...
    'all_trial_powratio_fullscan', ...
    'all_trial_powratio_fullscan_plotstat', ...
    'all_trial_powratio_early', 'all_trial_powratio_late', ...
    'all_trial_powratio_early_plotstat', 'all_trial_powratio_late_plotstat', ...
    'all_trial_peaks_single', 'all_trial_peaks_low', 'all_trial_peaks_high', 'all_trial_centroid', ...
    'all_trial_peaks_single_early', 'all_trial_peaks_single_late', ...
    'all_trial_mean_single', 'all_trial_median_single', ...
    'all_trial_mean_low', 'all_trial_median_low', ...
    'all_trial_mean_high', 'all_trial_median_high', 'all_trial_median_gap', ...
    'all_trial_mean_single_early', 'all_trial_median_single_early', ...
    'all_trial_mean_single_late', 'all_trial_median_single_late', ...
    'all_trial_median_low_early', 'all_trial_median_low_late', ...
    'all_trial_median_high_early', 'all_trial_median_high_late', ...
    'all_trial_median_gap_early', 'all_trial_median_gap_late', ...
    'all_trial_median_centroid_early', 'all_trial_median_centroid_late', ...
    'all_trial_trialcv_early', 'all_trial_trialcv_late', ...
    'all_trial_mean_centroid', 'all_trial_median_centroid', ...
    'all_trial_detrate_single', 'all_trial_detrate_low', 'all_trial_detrate_high', 'all_trial_detrate_single_early', 'all_trial_detrate_single_late', ...
    'all_trial_detrate_centroid', 'all_trial_gamma_power', 'all_trial_gamma_power_early', 'all_trial_gamma_power_late', ...
    'all_trial_gamma_power_plotstat', 'all_trial_gamma_power_early_plotstat', 'all_trial_gamma_power_late_plotstat', ...
    'all_topos', 'all_topos_early', 'all_topos_late', 'all_topo_labels', 'all_eigenvalues', ...
    'all_selected_comp_idx', 'all_selected_comp_corr', 'all_selected_comp_eval', ...
    'all_selected_comp_indices_multi', 'all_selected_comp_weights', ...
    'all_component_selection_stats_full', 'all_component_selection_stats_early', 'all_component_selection_stats_late', ...
    'warning_log', ...
    'all_trial_powratio_components_full', 'all_trial_powratio_components_early', 'all_trial_powratio_components_late', ...
    'all_trial_outlier_mask_freq_full', 'all_trial_outlier_mask_freq_early', 'all_trial_outlier_mask_freq_late', ...
    'all_trial_outlier_mask_power_full', 'all_trial_outlier_mask_power_early', 'all_trial_outlier_mask_power_late', ...
    'benchmark_metric_detectability', ...
    'benchmark_metric_separation_slope', 'benchmark_metric_separation_delta', ...
    'benchmark_metric_reliability_trialcv', 'benchmark_metric_reliability_subjspread', ...
    'primary_slope_stats', 'primary_delta_stats', ...
    'all_top5_corrs', 'all_top5_evals', 'all_top5_topos', 'all_simulated_templates', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames');

% Save candidate tables as CSV files (one file per subject and window).
candidate_table_csv_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/controls';
if ~exist(candidate_table_csv_dir, 'dir')
    mkdir(candidate_table_csv_dir);
end
nSubj = numel(subjects);
for si = 1:nSubj
    sid = subjects{si};
    if ~isempty(all_component_selection_stats_full{si}) && isfield(all_component_selection_stats_full{si}, 'candidate_table')
        ct = all_component_selection_stats_full{si}.candidate_table;
        if ~isempty(fieldnames(ct))
            T = struct2table(ct);
        else
            T = table();
        end
    else
        T = table();
    end
    writetable(T, fullfile(candidate_table_csv_dir, sprintf('subj%s_full.csv', sid)));

    if ~isempty(all_component_selection_stats_early{si}) && isfield(all_component_selection_stats_early{si}, 'candidate_table')
        ct = all_component_selection_stats_early{si}.candidate_table;
        if ~isempty(fieldnames(ct))
            T = struct2table(ct);
        else
            T = table();
        end
    else
        T = table();
    end
    writetable(T, fullfile(candidate_table_csv_dir, sprintf('subj%s_early.csv', sid)));

    if ~isempty(all_component_selection_stats_late{si}) && isfield(all_component_selection_stats_late{si}, 'candidate_table')
        ct = all_component_selection_stats_late{si}.candidate_table;
        if ~isempty(fieldnames(ct))
            T = struct2table(ct);
        else
            T = table();
        end
    else
        T = table();
    end
    writetable(T, fullfile(candidate_table_csv_dir, sprintf('subj%s_late.csv', sid)));
end

clc
fprintf('Done.\n');

function [trl_peaks_single, trl_peak_power, trl_peaks_low, trl_peaks_high, trl_centroid] = ...
    compute_trial_peak_metrics_from_powratio_fullscan(powratio_trials_fullscan, scan_freqs_full, analysis_mask, ...
    centroid_band_mask, smooth_n, centroid_min_peak, centroid_posfrac_min)
nTrl = size(powratio_trials_fullscan, 1);
trl_peaks_single = nan(nTrl, 1);
trl_peak_power = nan(nTrl, 1);
trl_peaks_low = nan(nTrl, 1);
trl_peaks_high = nan(nTrl, 1);
trl_centroid = nan(nTrl, 1);
scan_freqs_analysis = scan_freqs_full(analysis_mask);
freq_band = scan_freqs_full(centroid_band_mask);
for trl = 1:nTrl
    pr_full = powratio_trials_fullscan(trl, :);
    if all(~isfinite(pr_full))
        continue;
    end
    pr_proc = pr_full(analysis_mask);
    pr_proc = movmean(pr_proc, max(1, round(smooth_n)), 'omitnan');
    pr_proc = pr_proc(:)';
    valid = isfinite(pr_proc) & isfinite(scan_freqs_analysis);
    pr_proc = pr_proc(valid);
    x_use = scan_freqs_analysis(valid);
    if isempty(pr_proc)
        continue;
    end

    % Tallest local maximum from findpeaks; fall back to global max if
    % findpeaks returns nothing (e.g. monotonic or boundary peak).
    [peak_hz, peak_power] = pick_tallest_peak(pr_proc, x_use);
    if isfinite(peak_hz)
        trl_peaks_single(trl) = peak_hz;
        trl_peak_power(trl) = peak_power;
    end

    low_mask = x_use >= 30 & x_use <= 49;
    if any(low_mask)
        [lo_hz, ~] = pick_tallest_peak(pr_proc(low_mask), x_use(low_mask));
        if isfinite(lo_hz)
            trl_peaks_low(trl) = lo_hz;
        end
    end
    high_mask = x_use >= 50 & x_use <= 90;
    if any(high_mask)
        [hi_hz, ~] = pick_tallest_peak(pr_proc(high_mask), x_use(high_mask));
        if isfinite(hi_hz)
            trl_peaks_high(trl) = hi_hz;
        end
    end

    pr_proc_full = movmean(pr_full, max(1, round(smooth_n)), 'omitnan');
    pr_proc_band = pr_proc_full(centroid_band_mask);
    w_pos = max(pr_proc_band, 0);
    pos_mass = sum(w_pos);
    total_abs_mass = sum(abs(pr_proc_band));
    if pos_mass > 0 && max(pr_proc_band) >= centroid_min_peak && ...
            (pos_mass / max(total_abs_mass, eps)) >= centroid_posfrac_min
        trl_centroid(trl) = sum(freq_band .* w_pos) / pos_mass;
    end
end
end

function [peak_hz, peak_power] = pick_tallest_peak(y, x)
peak_hz = NaN;
peak_power = NaN;
y = y(:)';
x = x(:)';
if isempty(y) || numel(y) ~= numel(x)
    return;
end
if numel(y) >= 3
    [pks, locs] = findpeaks(y, x);
    if ~isempty(pks)
        [peak_power, bi] = max(pks);
        peak_hz = locs(bi);
        return;
    end
end
[peak_power, idx] = max(y);
if isfinite(peak_power)
    peak_hz = x(idx);
else
    peak_power = NaN;
    peak_hz = NaN;
end
end

function [ratio_db, near_floor_freq_mask] = compute_scan_ratio_from_timeseries(sig_stim, sig_base, fs, scan_freqs, scan_width_hz, base_floor, near_floor_mult)
if isvector(sig_stim)
    sig_stim = sig_stim(:)';
end
if isvector(sig_base)
    sig_base = sig_base(:)';
end
nSig = size(sig_stim, 1);
ratio_db = nan(nSig, numel(scan_freqs));
near_floor_freq_mask = false(1, numel(scan_freqs));
if isempty(sig_stim) || isempty(sig_base) || fs <= 0 || size(sig_stim, 1) ~= size(sig_base, 1)
    return;
end
if nargin < 7 || ~isfinite(near_floor_mult) || near_floor_mult <= 0
    near_floor_mult = 1.5;
end
n_fft = 2^nextpow2(max([size(sig_stim, 2), size(sig_base, 2), 256]));
f_axis = fs * (0:(n_fft/2)) / n_fft;
px_stim = abs(fft(sig_stim, n_fft, 2)).^2;
px_base = abs(fft(sig_base, n_fft, 2)).^2;
px_stim = px_stim(:, 1:numel(f_axis));
px_base = px_base(:, 1:numel(f_axis));
floor_use = max(base_floor, eps);
for fi = 1:numel(scan_freqs)
    f0 = scan_freqs(fi);
    f_mask = f_axis >= max(0, f0 - scan_width_hz) & f_axis <= (f0 + scan_width_hz);
    if ~any(f_mask)
        continue;
    end
    p_stim = mean(px_stim(:, f_mask), 2, 'omitnan');
    p_base = mean(px_base(:, f_mask), 2, 'omitnan');
    valid = isfinite(p_stim) & isfinite(p_base) & (p_base > 0);
    ratio_col = nan(nSig, 1);
    ratio_col(valid) = 10 * log10((p_stim(valid) + floor_use) ./ (p_base(valid) + floor_use));
    ratio_db(:, fi) = ratio_col;
    if any(valid)
        near_floor_freq_mask(fi) = mean(p_base(valid) <= near_floor_mult * floor_use) >= 0.5;
    end
end
end

function [base_floor, base_median] = compute_baseline_floor_stats(base_power_vals, floor_prctile, floor_frac)
base_floor = eps;
base_median = eps;
if nargin < 2 || ~isfinite(floor_prctile)
    floor_prctile = 20;
end
if nargin < 3 || ~isfinite(floor_frac)
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

function apply_dynamic_summary_ylims()
ax = findall(gcf, 'Type', 'axes');
for ai = 1:numel(ax)
    this_ax = ax(ai);
    if numel(this_ax.XTick) ~= 4
        continue;
    end
    y_vals = [];
    ch = this_ax.Children;
    for ci = 1:numel(ch)
        if isprop(ch(ci), 'YData')
            yd = ch(ci).YData;
            y_vals = [y_vals; yd(:)]; %#ok<AGROW>
        end
    end
    y_vals = y_vals(isfinite(y_vals));
    if isempty(y_vals)
        continue;
    end
    q05 = prctile(y_vals, 5);
    q95 = prctile(y_vals, 95);
    if ~isfinite(q05) || ~isfinite(q95)
        continue;
    end
    if q95 <= q05
        y_min = min(y_vals);
        y_max = max(y_vals);
    else
        pad = 0.12 * (q95 - q05);
        y_min = q05 - pad;
        y_max = q95 + pad;
    end
    if ~isfinite(y_min) || ~isfinite(y_max) || y_max <= y_min
        continue;
    end
    title_str = lower(char(this_ax.Title.String));
    if contains(title_str, 'trial cv')
        y_min = max(0, min(y_min, 0));
        y_max = max(y_max, 0.15);
    elseif contains(title_str, 'peak power') && contains(title_str, 'norm')
        y_min = max(0, min(y_min, 0.7));
        y_max = max(y_max, 1.3);
    elseif contains(title_str, 'peak power')
        y_min = min(y_min, 0);
        y_max = max(y_max, 0);
    end
    ylim(this_ax, [y_min y_max]);
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

function labels = assign_pf_quality_band_labels(pf_scores)
labels = repmat({'fail'}, numel(pf_scores), 1);
if isempty(pf_scores)
    return;
end
valid = isfinite(pf_scores);
if ~any(valid)
    return;
end
vals = pf_scores(valid);
q_high = prctile(vals, 80);
q_pass = prctile(vals, 60);
q_border = prctile(vals, 40);
for i = 1:numel(pf_scores)
    s = pf_scores(i);
    if ~isfinite(s)
        labels{i} = 'fail';
    elseif s >= q_high
        labels{i} = 'high_pass';
    elseif s >= q_pass
        labels{i} = 'pass';
    elseif s >= q_border
        labels{i} = 'borderline';
    else
        labels{i} = 'fail';
    end
end
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

function [score_vec, metrics, stability] = compute_calibrated_rank_aggregation_score( ...
    eval_raw_vec, peak_form_vec, peak_bonus_vec, occipital_evidence_vec, emg_artifact_vec, n_boot)
nComp = numel(eval_raw_vec);
score_vec = -Inf(nComp, 1);
metrics = struct( ...
    'score_raw', nan(nComp, 1), ...
    'eig_score', nan(nComp, 1), ...
    'peak_form_score', nan(nComp, 1), ...
    'peak_bonus_score', nan(nComp, 1), ...
    'occipital_score', nan(nComp, 1), ...
    'anti_emg_score', nan(nComp, 1));
stability = struct( ...
    'top1_freq', nan(nComp, 1), ...
    'rank_mean', nan(nComp, 1), ...
    'rank_std', nan(nComp, 1));
if nargin < 6 || isempty(n_boot) || ~isfinite(n_boot) || n_boot < 1
    n_boot = 0;
end
if nComp == 0
    return;
end

z_eig = normalize_robust(log(max(eval_raw_vec(:), eps)));
z_pf = normalize_robust(peak_form_vec(:));
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
score_raw = mean(point_mat, 2, 'omitnan');
score_vec = score_raw;
score_vec(~finite_rows) = -Inf;
metrics.score_raw = score_raw;
metrics.eig_score = point_mat(:, 1);
metrics.peak_form_score = point_mat(:, 2);
metrics.peak_bonus_score = point_mat(:, 3);
metrics.occipital_score = point_mat(:, 4);
metrics.anti_emg_score = point_mat(:, 5);

if n_boot <= 0
    return;
end
n_metrics = size(point_mat, 2);
boot_rank = nan(nComp, n_boot);
top1_count = zeros(nComp, 1);
for bi = 1:n_boot
    metric_idx = randi(n_metrics, [1, n_metrics]);
    boot_score = mean(point_mat(:, metric_idx), 2, 'omitnan');
    if ~any(isfinite(boot_score))
        continue;
    end
    rank_b = compute_descending_rank(boot_score);
    boot_rank(:, bi) = rank_b;
    [~, top_idx] = max(boot_score);
    if ~isempty(top_idx) && isfinite(boot_score(top_idx))
        top1_count(top_idx) = top1_count(top_idx) + 1;
    end
end
stability.top1_freq = top1_count / n_boot;
stability.rank_mean = mean(boot_rank, 2, 'omitnan');
stability.rank_std = std(boot_rank, 0, 2, 'omitnan');
end

function pf_summary = build_pf_benchmark_summary(tbl_full, tbl_early, tbl_late)
pf_summary = struct();
pf_summary.full = summarize_pf_band_counts(tbl_full);
pf_summary.early = summarize_pf_band_counts(tbl_early);
pf_summary.late = summarize_pf_band_counts(tbl_late);
end

function out = summarize_pf_band_counts(tbl)
out = struct('n_total', 0, 'n_high_pass', 0, 'n_pass', 0, 'n_borderline', 0, 'n_fail', 0, 'pf_median', NaN);
if isempty(tbl) || ~isfield(tbl, 'peak_form_score')
    return;
end
scores = tbl.peak_form_score;
if isempty(scores)
    return;
end
labels = tbl.pf_quality_band;
if isempty(labels)
    labels = assign_pf_quality_band_labels(scores);
end
out.n_total = numel(scores);
out.n_high_pass = sum(strcmpi(labels, 'high_pass'));
out.n_pass = sum(strcmpi(labels, 'pass'));
out.n_borderline = sum(strcmpi(labels, 'borderline'));
out.n_fail = sum(strcmpi(labels, 'fail'));
out.pf_median = median(scores(isfinite(scores)));
end

function plot_emg_exclusion_diagnostics(save_dir, subject_id, win_name, scan_freqs, searchTopos, ...
    searchMeanPrSpectrum, eigval_vec, emg_class, ...
    hard_eligible, force_include_occipital, rejection_flags, ...
    front_leak_vec, temp_leak_vec, lineharm_vec, hf_slope_vec, ...
    cfg_topo, topo_labels, ...
    peak_form_score, peak_form_mode, peak_form_best_center_hz, peak_form_dominant_penalty, ...
    selected_idx, fallback_selected_idx)
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
occipital_class_mask = cellfun(@(c) strcmpi(c, 'occipital'), emg_class(:));
final_eligible_mask = (logical(hard_eligible(:)) | logical(force_include_occipital(:))) & occipital_class_mask;
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
        if isfinite(fallback_selected_idx) && (ci == fallback_selected_idx)
            text(1.02, 0.50, 'FALLBACK', 'Units', 'normalized', 'Clipping', 'off', ...
                'Color', [0.85 0.05 0.05], 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                'FontSize', 8, 'Interpreter', 'none');
        end
    else
        axis off;
    end

    % Row 2: selected spectra
    subplot(4, nCols, nCols + k); hold on;
    if k <= nShowSel
        ci = sel_idx(k);
        spec_data = searchMeanPrSpectrum(ci, :);
        h_raw = plot(scan_freqs, spec_data, '-', 'Color', [0 0 0], 'LineWidth', 1.4);
        yline(0, 'k--');
        spec_min = min(spec_data(isfinite(spec_data)));
        spec_max = max(spec_data(isfinite(spec_data)));
        if isfinite(spec_min) && isfinite(spec_max) && spec_min ~= spec_max
            spec_range = spec_max - spec_min;
            ylim([spec_min - 0.10 * spec_range, spec_max + 0.10 * spec_range]);
        end
        [info_lines, info_viol] = ...
            build_rejection_info_columns(ci, rejection_flags, front_leak_vec, temp_leak_vec, lineharm_vec, hf_slope_vec);
        plot_rejection_info_text_columns(info_lines, info_viol);
        format_power_change_db_axis(gca);
        xlabel('Hz'); ylabel('Power [dB]');
        title(sprintf('\\lambda=%.2f, PF=%.2f', ...
            eigval_vec(ci), peak_form_score(ci)), 'FontSize', 7);
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
        h_raw_r = plot(scan_freqs, spec_data, '-', 'Color', [0 0 0], 'LineWidth', 1.4);
        yline(0, 'k--');
        spec_min = min(spec_data(isfinite(spec_data)));
        spec_max = max(spec_data(isfinite(spec_data)));
        if isfinite(spec_min) && isfinite(spec_max) && spec_min ~= spec_max
            spec_range = spec_max - spec_min;
            ylim([spec_min - 0.10 * spec_range, spec_max + 0.10 * spec_range]);
        end
        [info_lines, info_viol] = ...
            build_rejection_info_columns(ci, rejection_flags, front_leak_vec, temp_leak_vec, lineharm_vec, hf_slope_vec);
        plot_rejection_info_text_columns(info_lines, info_viol);
        format_power_change_db_axis(gca);
        xlabel('Hz'); ylabel('Power [dB]');
        title(sprintf('\\lambda=%.2f, PF=%.2f', ...
            eigval_vec(ci), peak_form_score(ci)), 'FontSize', 7);
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
save_figure_png(figSel, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_topo_spectra_selected_%s.png', subject_id, win_name)));
close(figSel);

end

function [info_lines, info_viol] = ...
    build_rejection_info_columns(ci, rejection_flags, front_leak_vec, temp_leak_vec, lineharm_vec, hf_slope_vec)
info_lines = { ...
    sprintf('lineharm: %.2f', lineharm_vec(ci)), ...
    sprintf('hf_slope: %.2f', hf_slope_vec(ci)), ...
    sprintf('front_leak: %.2f', front_leak_vec(ci)), ...
    sprintf('temp_leak: %.2f', temp_leak_vec(ci))};
info_viol = [ ...
    get_flag_value(rejection_flags, 'lineharm', ci), ...
    get_flag_value(rejection_flags, 'hf_slope', ci), ...
    get_flag_value(rejection_flags, 'front_leak', ci), ...
    get_flag_value(rejection_flags, 'temp_leak', ci)];
end

function plot_rejection_info_text_columns(info_lines, info_viol)
y0 = 1.52;
dy = 0.11;
for li = 1:numel(info_lines)
    if info_viol(li)
        txt_col = [0.82 0.10 0.10];
    else
        txt_col = [0.10 0.10 0.10];
    end
    text(0.02, y0 - (li - 1) * dy, info_lines{li}, ...
        'Units', 'normalized', 'Clipping', 'off', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'FontSize', 6.5, 'Color', txt_col, 'Interpreter', 'none');
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
    [pf_score, pf_peak_hz] = compute_combined_peak_form_metrics( ...
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
if nargin < 3 || isempty(nComp) || ~isfinite(nComp)
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

function [pf_score, pf_peak_hz] = compute_combined_peak_form_metrics( ...
    spec_vec, scan_freqs, analysis_freq_range)
pf_score = NaN;
pf_peak_hz = NaN;
if isempty(spec_vec) || isempty(scan_freqs) || numel(spec_vec) ~= numel(scan_freqs)
    return;
end
[pf_vec, ~, pf_diag] = compute_peak_form_template_score_from_spectra( ...
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
        peak_lo = max(40, analysis_freq_range(1));
        peak_hi = min(80, analysis_freq_range(2));
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
            peak_lo = max(40, analysis_freq_range(1));
            peak_hi = min(80, analysis_freq_range(2));
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

function format_power_change_db_axis(axh)
if nargin < 1 || isempty(axh) || ~isgraphics(axh, 'axes')
    axh = gca;
end
yt = get(axh, 'YTick');
if isempty(yt)
    return;
end
yt_lbl = arrayfun(@(v) sprintf('%g', v), yt, 'UniformOutput', false);
set(axh, 'YTickLabel', yt_lbl);
end

function [edge_ratio, edge_run_score, edge_artifact_flag] = compute_edge_artifact_indicators( ...
    y_resid, x_band, analysis_freq_range, edge_ratio_hard, edge_run_hard)
edge_ratio = NaN;
edge_run_score = NaN;
edge_artifact_flag = false;
if isempty(y_resid) || isempty(x_band) || numel(y_resid) ~= numel(x_band)
    return;
end
x_lo = analysis_freq_range(1);
x_hi = analysis_freq_range(2);
edge_mask = (x_band <= (x_lo + 5)) | (x_band >= (x_hi - 5));
interior_mask = x_band >= (x_lo + 10) & x_band <= (x_hi - 10);
if sum(edge_mask) < 3 || sum(interior_mask) < 3
    return;
end
edge_amp = prctile(abs(y_resid(edge_mask)), 90);
interior_amp = prctile(abs(y_resid(interior_mask)), 90);
edge_ratio = edge_amp / max(interior_amp, eps);
left_mask = x_band <= (x_lo + 6);
right_mask = x_band >= (x_hi - 6);
left_run = 0;
right_run = 0;
if sum(left_mask) >= 4
    left_run = abs(mean(diff(y_resid(left_mask))));
end
if sum(right_mask) >= 4
    right_run = abs(mean(diff(y_resid(right_mask))));
end
edge_run_score = max(left_run, right_run);
edge_artifact_flag = isfinite(edge_ratio) && (edge_ratio > edge_ratio_hard) && ...
    isfinite(edge_run_score) && (edge_run_score > edge_run_hard);
end

function [peak_bonus_vec, peak_count_vec] = compute_peak_bonus_from_spectra( ...
    mean_pr_spectrum, scan_freqs, analysis_freq_range)
nComp = size(mean_pr_spectrum, 1);
peak_bonus_vec = zeros(nComp, 1);
peak_count_vec = zeros(nComp, 1);
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
    [pks, ~] = findpeaks(y_shape, x_band, 'MinPeakDistance', 5);
    if isempty(pks)
        pks = peak_scale;
    end
    pks = sort(pks(:), 'descend');
    n_keep = min(2, numel(pks));
    peak_count_vec(ci) = n_keep;
    amp_score = min(1, mean(pks(1:n_keep)) / max(peak_scale, eps));
    count_score = n_keep / 2;
    peak_bonus_vec(ci) = max(0, min(1, 0.65 * amp_score + 0.35 * count_score));
end
end

function [peak_form_score_vec, peak_form_mode_vec, diag] = compute_peak_form_template_score_from_spectra( ...
    mean_pr_spectrum, scan_freqs, analysis_freq_range)
shift_max_hz = 10;
single_widths_hz = [5 7 9 12];
double_separations_hz = [8 12 16 20];
double_widths_hz = [3 4 5 6];
min_trough_depth = 0.10;
min_similarity = 0.40;
smooth_n = 3;
peak_width_min_hz = 2.0;
peak_width_max_hz = 20.0;
edge_ratio_soft = 1.25;
edge_ratio_hard = 1.75;
edge_run_soft = 0.025;
edge_run_hard = 0.060;
nComp = size(mean_pr_spectrum, 1);
peak_form_score_vec = zeros(nComp, 1);
peak_form_mode_vec = repmat({'none'}, nComp, 1);
diag = struct( ...
    'best_single_similarity', nan(nComp, 1), ...
    'best_double_similarity', nan(nComp, 1), ...
    'best_similarity_raw', nan(nComp, 1), ...
    'pre_penalty_similarity', nan(nComp, 1), ...
    'total_penalty_raw', nan(nComp, 1), ...
    'total_penalty_used', nan(nComp, 1), ...
    'edge_penalty', nan(nComp, 1), ...
    'dominance_penalty', nan(nComp, 1), ...
    'hf_rise_penalty', nan(nComp, 1), ...
    'best_shift_hz', nan(nComp, 1), ...
    'best_center_hz', nan(nComp, 1), ...
    'best_width_hz', nan(nComp, 1), ...
    'best_separation_hz', nan(nComp, 1), ...
    'best_trough_depth', nan(nComp, 1), ...
    'minor_peak_count', nan(nComp, 1), ...
    'minor_peak_relmax', nan(nComp, 1), ...
    'edge_ratio', nan(nComp, 1), ...
    'edge_run_score', nan(nComp, 1), ...
    'edge_artifact_flag', false(nComp, 1), ...
    'roughness_ratio', nan(nComp, 1), ...
    'roughness_penalty', nan(nComp, 1), ...
    'dominant_penalty_tag', {repmat({'none'}, nComp, 1)});
if isempty(mean_pr_spectrum) || isempty(scan_freqs)
    return;
end
freq_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
if ~any(freq_mask)
    return;
end
smooth_n = max(1, round(smooth_n));
for ci = 1:nComp
    y = mean_pr_spectrum(ci, :);
    if all(~isfinite(y))
        continue;
    end
    x_band = scan_freqs(freq_mask);
    y_proc = y(freq_mask);
    y_proc = movmean(y_proc, smooth_n, 'omitnan');
    y_band = y_proc;
    valid = isfinite(x_band) & isfinite(y_band);
    x_band = x_band(valid);
    y_band = y_band(valid);
    y_resid = y_band;
    if numel(y_band) < 7
        continue;
    end
    y_shape = max(y_band - median(y_band(isfinite(y_band))), 0);
    if ~any(isfinite(y_shape))
        y_shape = zeros(size(y_band));
    end
    [single_best, single_meta] = evaluate_single_template_bank(y_shape, x_band, single_widths_hz, shift_max_hz);
    [double_best, double_meta] = evaluate_double_template_bank( ...
        y_shape, y_band, x_band, double_widths_hz, double_separations_hz, shift_max_hz, min_trough_depth);
    diag.best_single_similarity(ci) = single_best;
    diag.best_double_similarity(ci) = double_best;

    % Raw spectra can be globally offset, so peak shape is measured relative
    % to a local raw floor.
    y_floor = prctile(y_band(isfinite(y_band)), 20);
    if ~isfinite(y_floor)
        y_floor = median(y_band(isfinite(y_band)));
    end
    if ~isfinite(y_floor)
        y_floor = 0;
    end
    y_pos = max(y_band - y_floor, 0);
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
    [pks, locs, widths] = findpeaks(y_pos, x_band, 'MinPeakDistance', 5);
    if isempty(pks)
        % Fallback: keep best local maximum to avoid collapsing valid-but-noisy traces to PF=0.
        [dom_amp, dom_idx] = max(y_pos);
        if ~isfinite(dom_amp) || dom_amp <= eps
            continue;
        end
        dom_loc = x_band(dom_idx);
        widths = 6;
        pks = dom_amp;
        locs = dom_loc;
    end
    [~, dom_idx] = max(pks);
    dom_amp = pks(dom_idx);
    dom_width = widths(dom_idx);
    dom_loc = locs(dom_idx);
    minor_idx = setdiff(1:numel(pks), dom_idx);
    minor_count = numel(minor_idx);
    minor_relmax = 0;
    if ~isempty(minor_idx)
        minor_relmax = max(pks(minor_idx)) / max(dom_amp, eps);
    end
    diag.minor_peak_count(ci) = minor_count;
    diag.minor_peak_relmax(ci) = minor_relmax;

    width_low = peak_width_min_hz;
    width_high = peak_width_max_hz;
    width_score = 1;
    if isfinite(dom_width)
        if dom_width < width_low
            width_score = max(0.15, dom_width / max(width_low, eps));
        elseif dom_width > width_high
            width_score = max(0.20, width_high / max(dom_width, eps));
        end
    end
    snr_score = dom_amp / (dom_amp + robust_scale);
    amp_score = dom_amp / (dom_amp + 0.70 * robust_scale);
    shape_score = max(single_best, double_best);
    peak_core_hz = 6;
    peak_core_mask = abs(x_band - dom_loc) <= peak_core_hz;
    conc_ratio = sum(y_pos(peak_core_mask)) / max(sum(y_pos), eps);
    concentration_pen = 1;
    if isfinite(conc_ratio) && conc_ratio < 0.32
        concentration_pen = max(0.35, conc_ratio / 0.32);
    end
    dominant_quality = max(0, min(1, (0.40 * amp_score + 0.30 * snr_score + 0.30 * shape_score) * width_score * concentration_pen));

    best_pre_penalty = dominant_quality;
    dual_candidate = isfinite(double_best) && isfinite(single_best) && ...
        (double_best > (single_best + 0.05)) && isfinite(double_meta.trough_depth) && ...
        (double_meta.trough_depth >= max(0.05, 0.75 * min_trough_depth));
    if dual_candidate
        mode_raw = 'dual';
        best_shift = double_meta.shift_hz;
        best_width = double_meta.width_hz;
        best_sep = double_meta.sep_hz;
        best_trough = double_meta.trough_depth;
        if ~isempty(double_meta.centers_hz)
            best_centers = double_meta.centers_hz;
        else
            best_centers = dom_loc;
        end
    else
        mode_raw = 'single';
        best_shift = single_meta.shift_hz;
        best_width = dom_width;
        best_sep = double_meta.sep_hz;
        best_trough = double_meta.trough_depth;
        best_centers = dom_loc;
    end

    % Penalize EMG-like monotonic high-frequency rise in the upper gamma band.
    hf_pen = 1;
    hf_mask = x_band >= max(70, analysis_freq_range(2) - 15);
    if sum(hf_mask) >= 5
        hf_idx = find(hf_mask);
        hf_rho = corr((1:numel(hf_idx))', y_band(hf_idx)', 'rows', 'complete', 'type', 'Spearman');
        if isfinite(hf_rho) && hf_rho > 0.70
            hf_pen = max(0.65, 1 - 0.30 * (hf_rho - 0.70) / 0.30);
        end
    end

    minor_pen = 1;
    if minor_count > 0
        % Slightly stronger crowding penalty: messy spectra with several
        % minors should score lower, while one small shoulder remains tolerated.
        if minor_count >= 2
            minor_pen = minor_pen * max(0.75, 1 - 0.08 * (minor_count - 1));
        end

        % Strong penalty only when a secondary peak is large relative to the
        % dominant peak (C11/C14-like multi-peak pattern).
        strong_minor_rel_start = 0.55;
        strong_minor_rel_full = 0.85;
        if isfinite(minor_relmax) && minor_relmax > strong_minor_rel_start
            rel_span = max(strong_minor_rel_full - strong_minor_rel_start, eps);
            rel_frac = min(1, (minor_relmax - strong_minor_rel_start) / rel_span);
            minor_pen = minor_pen * max(0.30, 1 - 0.70 * rel_frac);
        end
    end

    [edge_ratio, edge_run_score, edge_flag] = compute_edge_artifact_indicators( ...
        y_resid, x_band, analysis_freq_range, edge_ratio_hard, edge_run_hard);
    diag.edge_ratio(ci) = edge_ratio;
    diag.edge_run_score(ci) = edge_run_score;
    diag.edge_artifact_flag(ci) = edge_flag;
    edge_pen = 1;
    if isfinite(edge_ratio) && edge_ratio > edge_ratio_soft
        edge_pen = edge_pen * max(0.80, 1 - 0.20 * min(1, (edge_ratio - edge_ratio_soft) / max(edge_ratio_hard - edge_ratio_soft, eps)));
    end
    if isfinite(edge_run_score) && edge_run_score > edge_run_soft
        edge_pen = edge_pen * max(0.80, 1 - 0.20 * min(1, (edge_run_score - edge_run_soft) / max(edge_run_hard - edge_run_soft, eps)));
    end
    if edge_flag
        edge_pen = edge_pen * 0.90;
    end

    dominance_pen = minor_pen;

    % Penalize rough/noisy spectra that can spuriously match broad templates.
    % Scale-invariant: normalize to unit range so absolute amplitude does not bias the ratio.
    roughness_ratio = NaN;
    roughness_pen = 1;
    y_finite = y_band(isfinite(y_band));
    if numel(y_finite) >= 3
        y_range = max(y_finite) - min(y_finite);
        if y_range > eps
            y_scaled = y_band / y_range;
        else
            y_scaled = y_band;
        end
        dy = diff(y_scaled);
        amp_scale = prctile(y_scaled, 75) - prctile(y_scaled, 25);
        if ~isfinite(amp_scale) || amp_scale <= eps
            amp_scale = std(y_scaled(isfinite(y_scaled)));
        end
        if ~isfinite(amp_scale) || amp_scale <= eps
            amp_scale = 1;
        end
        roughness_ratio = robust_mad(dy) / amp_scale;
        rough_ref = 0.32;
        rough_span = 0.58;
        rough_pen_floor = 0.38;
        rough_pen_max_loss = 0.62;
        if isfinite(roughness_ratio) && roughness_ratio > rough_ref
            loss_frac = min(1, (roughness_ratio - rough_ref) / rough_span);
            roughness_pen = max(rough_pen_floor, 1 - rough_pen_max_loss * loss_frac);
        end
    end

    penalty_raw = edge_pen * dominance_pen * hf_pen * roughness_pen;
    penalty_floor = 0.20;
    penalty_used = max(penalty_floor, penalty_raw);
    best_raw = best_pre_penalty * penalty_used;

    penalty_names = {'edge', 'dominance', 'hf_rise', 'roughness'};
    penalty_vals = [edge_pen, dominance_pen, hf_pen, roughness_pen];
    [worst_pen, worst_idx] = min(penalty_vals);
    if isfinite(worst_pen) && (worst_pen < 0.999)
        dominant_penalty_tag = penalty_names{worst_idx};
    else
        dominant_penalty_tag = 'none';
    end

    diag.best_similarity_raw(ci) = best_raw;
    diag.pre_penalty_similarity(ci) = best_pre_penalty;
    diag.total_penalty_raw(ci) = penalty_raw;
    diag.total_penalty_used(ci) = penalty_used;
    diag.edge_penalty(ci) = edge_pen;
    diag.dominance_penalty(ci) = dominance_pen;
    diag.hf_rise_penalty(ci) = hf_pen;
    diag.best_shift_hz(ci) = best_shift;
    if ~isempty(best_centers) && all(isfinite(best_centers))
        diag.best_center_hz(ci) = mean(best_centers);
    end
    diag.best_width_hz(ci) = best_width;
    diag.best_separation_hz(ci) = best_sep;
    diag.best_trough_depth(ci) = best_trough;
    diag.roughness_ratio(ci) = roughness_ratio;
    diag.roughness_penalty(ci) = roughness_pen;
    diag.dominant_penalty_tag{ci} = dominant_penalty_tag;

    if ~isfinite(best_raw)
        peak_form_score_vec(ci) = 0;
        peak_form_mode_vec{ci} = 'none';
    else
        score_mapped = best_raw;
        if isfinite(min_similarity) && (min_similarity > 0) && (best_raw < min_similarity)
            score_mapped = 0.75 * score_mapped;
        end
        peak_form_score_vec(ci) = max(0, min(1, score_mapped));
        peak_form_mode_vec{ci} = mode_raw;
    end
end
end

function [best_sim, meta] = evaluate_single_template_bank(y_shape, x_band, widths_hz, shift_max_hz)
best_sim = 0;
meta = struct('shift_hz', NaN, 'width_hz', NaN, 'centers_hz', []);
if isempty(widths_hz)
    return;
end
x_mid = mean(x_band);
x_min = min(x_band);
x_max = max(x_band);
for wi = 1:numel(widths_hz)
    w = widths_hz(wi);
    if ~isfinite(w) || w <= 0
        continue;
    end
    shift_vals = candidate_shift_values_hz(x_band, shift_max_hz);
    for si = 1:numel(shift_vals)
        shift_hz = shift_vals(si);
        center_hz = x_mid + shift_hz;
        if ~isfinite(center_hz) || center_hz < x_min || center_hz > x_max
            continue;
        end
        tpl = gaussian_template(x_band, center_hz, w);
        sim = safe_template_similarity(y_shape, tpl);
        if sim > best_sim
            best_sim = sim;
            meta.shift_hz = shift_hz;
            meta.width_hz = w;
            meta.centers_hz = center_hz;
        end
    end
end
end

function [best_sim, meta] = evaluate_double_template_bank(y_shape, y_band, x_band, widths_hz, separations_hz, shift_max_hz, min_trough_depth)
best_sim = 0;
meta = struct('shift_hz', NaN, 'width_hz', NaN, 'sep_hz', NaN, 'trough_depth', NaN, 'centers_hz', []);
if isempty(widths_hz) || isempty(separations_hz)
    return;
end
x_mid = mean(x_band);
x_min = min(x_band);
x_max = max(x_band);
shift_vals = candidate_shift_values_hz(x_band, shift_max_hz);
for wi = 1:numel(widths_hz)
    w = widths_hz(wi);
    if ~isfinite(w) || w <= 0
        continue;
    end
    for di = 1:numel(separations_hz)
        sep = separations_hz(di);
        if ~isfinite(sep) || sep <= 0
            continue;
        end
        for si = 1:numel(shift_vals)
            shift_hz = shift_vals(si);
            c1 = x_mid + shift_hz - sep / 2;
            c2 = x_mid + shift_hz + sep / 2;
            if ~isfinite(c1) || ~isfinite(c2) || c1 < x_min || c2 > x_max
                continue;
            end
            tpl = gaussian_template(x_band, c1, w) + gaussian_template(x_band, c2, w);
            sim = safe_template_similarity(y_shape, tpl);
            trough_depth = estimate_trough_depth(y_band, x_band, c1, c2);
            trough_scale = min(1, max(0, trough_depth) / max(min_trough_depth, eps));
            sim_adj = sim * trough_scale;
            if sim_adj > best_sim
                best_sim = sim_adj;
                meta.shift_hz = shift_hz;
                meta.width_hz = w;
                meta.sep_hz = sep;
                meta.trough_depth = trough_depth;
                meta.centers_hz = [c1 c2];
            end
        end
    end
end
end

function shifts = candidate_shift_values_hz(x_band, shift_max_hz)
if numel(x_band) >= 2
    df = median(diff(x_band));
else
    df = 1;
end
if ~isfinite(df) || df <= 0
    df = 1;
end
shift_max_hz = max(0, shift_max_hz);
max_possible = max(abs(x_band - mean(x_band)));
if ~isfinite(max_possible) || max_possible <= 0
    max_possible = shift_max_hz;
end
shift_limit = min(shift_max_hz, max_possible);
shifts = -shift_limit:df:shift_limit;
if isempty(shifts)
    shifts = 0;
end
if ~any(abs(shifts) < 1e-9)
    shifts = unique(sort([shifts, 0]));
end
end

function depth = estimate_trough_depth(y_band, x_band, c1, c2)
depth = 0;
if ~(isfinite(c1) && isfinite(c2))
    return;
end
if c1 > c2
    tmp = c1;
    c1 = c2;
    c2 = tmp;
end
[~, i1] = min(abs(x_band - c1));
[~, i2] = min(abs(x_band - c2));
if i1 == i2
    return;
end
idx_lo = min(i1, i2);
idx_hi = max(i1, i2);
y_pos = max(y_band(:), 0);
p1 = y_pos(i1);
p2 = y_pos(i2);
if idx_hi - idx_lo < 2 || ~isfinite(p1) || ~isfinite(p2)
    return;
end
valley = min(y_pos(idx_lo:idx_hi));
peak_ref = max(min(p1, p2), eps);
depth = max(0, min(1, 1 - valley / peak_ref));
end

function sim = safe_template_similarity(y_shape, tpl)
sim = 0;
if isempty(y_shape) || isempty(tpl)
    return;
end
ys = y_shape(:);
ys(~isfinite(ys)) = 0;
ys = max(ys, 0);
tpl = tpl(:);
tpl = max(tpl, 0);
if numel(tpl) ~= numel(ys)
    return;
end
mass = sum(ys);
if ~isfinite(mass) || mass <= eps
    return;
end
sim = sum(ys .* tpl) / mass;
if ~isfinite(sim)
    sim = 0;
end
sim = max(0, min(1, sim));
end

function g = gaussian_template(x, mu, sigma)
if ~isfinite(mu) || ~isfinite(sigma) || sigma <= 0
    g = zeros(size(x));
    return;
end
g = exp(-0.5 * ((x - mu) ./ sigma).^2);
end

function fallback_idx = get_candidate_table_fallback_idx(candidate_table)
fallback_idx = NaN;
if isfield(candidate_table, 'fallback_selected')
    idx = find(logical(candidate_table.fallback_selected), 1, 'first');
    if ~isempty(idx)
        fallback_idx = idx;
    end
end
end

function thr = build_adaptive_component_thresholds(eval_raw_vec, corr_vec, leak_vec, temp_leak_vec, ...
    lineharm_vec, stationarity_vec, burst_vec, hf_slope_vec, condlock_vec, occ_evidence, emg_score, defaults)
thr = struct();
thr.min_eigval = min(defaults.min_eigval, max(defaults.min_eig_floor, robust_percentile(eval_raw_vec, 35, defaults.min_eigval)));
thr.max_frontleak = max(defaults.max_frontleak, robust_upper_bound(leak_vec, defaults.max_frontleak, defaults.artifact_mad_mult));
thr.max_templeak = max(defaults.max_templeak, robust_upper_bound(temp_leak_vec, defaults.max_templeak, defaults.artifact_mad_mult));
thr.min_corr = min(defaults.min_corr, robust_lower_bound(corr_vec, defaults.min_corr, defaults.artifact_mad_mult));
thr.max_lineharm = max(defaults.max_lineharm, robust_upper_bound(lineharm_vec, defaults.max_lineharm, defaults.artifact_mad_mult));
thr.max_stationarity = max(defaults.max_stationarity, robust_upper_bound(stationarity_vec, defaults.max_stationarity, defaults.artifact_mad_mult));
thr.max_burst_ratio = max(defaults.max_burst_ratio, robust_upper_bound(burst_vec, defaults.max_burst_ratio, defaults.artifact_mad_mult));
thr.max_hf_slope = max(defaults.max_hf_slope, robust_upper_bound(hf_slope_vec, defaults.max_hf_slope, defaults.artifact_mad_mult));
thr.min_condlock = min(defaults.min_condlock, robust_lower_bound(condlock_vec, defaults.min_condlock, defaults.artifact_mad_mult));
thr.occ_class_thr = robust_percentile(occ_evidence, defaults.class_quantile, 0.6);
thr.emg_class_thr = robust_percentile(emg_score, defaults.class_quantile, 0.5);
thr.max_emg_score = max(defaults.max_emg_score, robust_upper_bound(emg_score, defaults.max_emg_score, defaults.artifact_mad_mult));
occ_minus_emg = occ_evidence(:) - emg_score(:);
margin_mad = robust_mad(occ_minus_emg);
thr.min_occ_margin = max(defaults.occ_margin_floor, 0.30 * margin_mad);
end

function val = robust_upper_bound(x, default_val, mad_mult)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    val = default_val;
    return;
end
val = median(x) + mad_mult * robust_mad(x);
if ~isfinite(val)
    val = default_val;
end
end

function val = robust_lower_bound(x, default_val, mad_mult)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    val = default_val;
    return;
end
val = median(x) - mad_mult * robust_mad(x);
if ~isfinite(val)
    val = default_val;
end
end

function val = robust_percentile(x, pctl, default_val)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    val = default_val;
    return;
end
val = prctile(x, pctl);
if ~isfinite(val)
    val = default_val;
end
end

function [outlier_mask, stats] = detect_trial_metric_outliers_iqr(x, iqr_mult, min_trials)
x = x(:);
outlier_mask = false(size(x));
stats = struct('n_finite', 0, 'q1', NaN, 'q3', NaN, 'iqr_val', NaN, 'lo', NaN, 'hi', NaN, 'n_outliers', 0);
valid = isfinite(x);
n_finite = sum(valid);
stats.n_finite = n_finite;
if n_finite < max(3, min_trials)
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

function reasons = compute_primary_rejection_reason(rejection_flags)
reason_order = {'unknown_proxy', 'front_leak', 'temp_leak', 'corr', ...
    'lineharm', 'stationarity', 'burst', 'hf_slope', 'condlock', 'emg_score', 'occ_margin', ...
    'topo_flat', 'topo_nonposterior', 'topo_fragmented'};
nComp = numel(rejection_flags.unknown_proxy);
reasons = repmat({'pass'}, nComp, 1);
for ci = 1:nComp
    for ri = 1:numel(reason_order)
        rname = reason_order{ri};
        if isfield(rejection_flags, rname) && rejection_flags.(rname)(ci)
            reasons{ci} = rname;
            break;
        end
    end
end
end

function proxy = estimate_component_artifact_proxies(filter_w, dat_per_cond, stim_window, base_window, scan_freqs, scan_width, max_trials)
proxy = struct( ...
    'lineharm_ratio', NaN, ...
    'stationarity_cv', NaN, ...
    'burst_ratio', NaN, ...
    'hf_slope', NaN, ...
    'cond_lock_rho', NaN, ...
    'unknown_high_risk', true, ...
    'mean_pr_spectrum', nan(1, numel(scan_freqs)), ...
    'n_trials_used', 0);
if isempty(filter_w) || isempty(dat_per_cond)
    return;
end
if nargin < 8 || isempty(max_trials)
    max_trials = inf;
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
trial_cond = [];
lineharm_acc = [];
trial_count = 0;
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
        pr_full = compute_simple_power_ratio_scan(x_proj(idx_stim), x_proj(idx_base), dat.fsample, scan_freqs, scan_width);
        if isempty(pr_full) || all(~isfinite(pr_full))
            continue;
        end
        pr_full = pr_full(:)';
        if ~any(band_mask)
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
        trial_gamma(end+1, 1) = nanmean(pr_band); %#ok<AGROW>
        trial_pr(end+1, :) = pr_full; %#ok<AGROW>
        trial_cond(end+1, 1) = cond; %#ok<AGROW>
        lineharm_acc(end+1, 1) = max(harm_val, 0) / max(abs(nonharm_val), eps); %#ok<AGROW>
        trial_count = trial_count + 1;
        if trial_count >= max_trials
            break;
        end
    end
    if trial_count >= max_trials
        break;
    end
end

proxy.n_trials_used = trial_count;
if trial_count < 4
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
    xh = log(scan_freqs(hf_mask));
    yh = log(abs(proxy.mean_pr_spectrum(hf_mask)) + eps);
    p = polyfit(xh(:), yh(:), 1);
    proxy.hf_slope = p(1);
end

cond_means = nan(4, 1);
for c = 1:4
    cond_means(c) = nanmean(trial_gamma(trial_cond == c));
end
valid_cond = isfinite(cond_means);
if sum(valid_cond) >= 3
    r = corr((1:4)', cond_means, 'rows', 'complete', 'type', 'Spearman');
    proxy.cond_lock_rho = r;
end

proxy.unknown_high_risk = ~(isfinite(proxy.lineharm_ratio) && isfinite(proxy.stationarity_cv) && ...
    isfinite(proxy.burst_ratio) && isfinite(proxy.hf_slope) && isfinite(proxy.cond_lock_rho));
end

function bad_mask = flag_unreliable_baseline_trials(base_power_vals, mad_mult, min_trials, enable_flag)
bad_mask = false(size(base_power_vals));
if nargin < 4 || ~enable_flag || isempty(base_power_vals)
    return;
end
valid = isfinite(base_power_vals) & (base_power_vals > 0);
if sum(valid) < max(3, min_trials)
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

function pr_scan = compute_simple_power_ratio_scan(x_stim, x_base, fs, scan_freqs, scan_width)
pr_scan = nan(size(scan_freqs));
if isempty(x_stim) || isempty(x_base) || fs <= 0
    return;
end
x_stim = x_stim(:)';
x_base = x_base(:)';
n_fft = 2^nextpow2(max([numel(x_stim), numel(x_base), 256]));
f_axis = fs * (0:(n_fft/2)) / n_fft;
px_stim = abs(fft(x_stim, n_fft)).^2;
px_base = abs(fft(x_base, n_fft)).^2;
px_stim = px_stim(1:numel(f_axis));
px_base = px_base(1:numel(f_axis));
base_band_scan = nan(size(scan_freqs));
for fi = 1:numel(scan_freqs)
    f0 = scan_freqs(fi);
    f_lo = max(0, f0 - scan_width);
    f_hi = f0 + scan_width;
    f_mask = f_axis >= f_lo & f_axis <= f_hi;
    if ~any(f_mask)
        continue;
    end
    p_stim = mean(px_stim(f_mask));
    p_base = mean(px_base(f_mask));
    base_band_scan(fi) = p_base;
    pr_scan(fi) = p_stim;
end
valid_base = isfinite(base_band_scan) & (base_band_scan > 0);
if ~any(valid_base)
    pr_scan(:) = NaN;
    return;
end
base_anchor = prctile(base_band_scan(valid_base), 20);
base_median = median(base_band_scan(valid_base), 'omitnan');
if ~isfinite(base_anchor) || base_anchor <= 0
    base_anchor = base_median;
end
if ~isfinite(base_median) || base_median <= 0
    base_median = base_anchor;
end
base_floor = max(0.25 * base_anchor, eps);
for fi = 1:numel(scan_freqs)
    p_stim = pr_scan(fi);
    p_base = base_band_scan(fi);
    if ~isfinite(p_stim) || ~isfinite(p_base) || p_base <= 0
        pr_scan(fi) = NaN;
        continue;
    end
    % dB ratio with robust additive floor for numerical stability.
    pr_scan(fi) = 10 * log10((p_stim + base_floor) / (p_base + base_floor));
end
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
    warning_log(end+1, 1) = entry;
end
end

function print_subject_warning_summary(warning_log)
fprintf('Per-subject warning summary\n');
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
        fprintf(['     selection diagnostics: nSearch=%d, finite=%d, passEig=%d, artifactFlagged=%d, ', ...
                 'unknownHighRisk=%d, passRaw=%d, excludedOutlier=%d\n'], ...
            m.n_search, m.n_finite_metrics, m.n_pass_eig, m.n_artifact_flagged, ...
            m.n_unknown_high_risk, m.n_pass_all_raw, m.n_excluded_dominant_outlier);
        if isfield(m, 'top_fail_idx') && isfinite(m.top_fail_idx)
            fprintf(['     top failing component: C%d | eig=%.3f (thr=%.3f, pass=%d), ', ...
                     'corr=%.3f (thr=%.3f), ratio=%.3f, artifact=%d\n'], ...
                m.top_fail_idx, m.top_fail_eig, m.thr_eig, m.top_fail_pass_eig, ...
                m.top_fail_corr, m.thr_corr, ...
                m.top_fail_ratio, m.top_fail_artifact);
        end
    elseif strcmp(w.code, 'NO_OCCIPITAL_COMPONENTS') || strcmp(w.code, 'TOO_FEW_HARD_COMPONENTS')
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

function save_figure_png(fig_handle, out_path)
if nargin < 1 || isempty(fig_handle) || ~ishandle(fig_handle)
    return;
end
if nargin < 2 || isempty(out_path)
    return;
end
drawnow;
pause(0.05);
try
    exportgraphics(fig_handle, out_path, 'Resolution', 300);
catch
    exportgraphics(fig_handle, out_path, 'Resolution', 600);
end
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
