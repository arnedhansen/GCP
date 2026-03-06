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
%             top GED component, and combined artifact-screened GED components.
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

% Gamma frequency range
gamma_range = [30, 90];

% Narrowband scanning parameters
% Detrending fit is computed on a wider band, then trimmed to analysis band.
scan_freqs_detrend = 10:1:110;
analysis_freq_range = [30, 90];
analysis_freq_mask_detrend = scan_freqs_detrend >= analysis_freq_range(1) & scan_freqs_detrend <= analysis_freq_range(2);
scan_freqs = scan_freqs_detrend(analysis_freq_mask_detrend);
nFreqs_detrend = length(scan_freqs_detrend);
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
viz_suppress_nonocc_outliers = false;  % visualization-only suppression/interpolation
viz_interp_k = 6;                   % nearest neighbors for channel interpolation
viz_nonocc_outlier_mult = 1.00;     % non-occipital outlier threshold multiplier (vs posterior pctl)
viz_topo_prctile = 99.9;            % robust percentile for color scaling
min_eigval_hard = 1.0;              % hard minimum GED eigenvalue (lambda > 1)
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
exclude_flat_topomap = true;        % remove spatially flat/one-color topographies
flat_topo_mad_ratio_thr = 0.18;     % relaxed flat-map guard: MAD relative to maxabs
flat_topo_iqr_ratio_thr = 0.24;     % relaxed flat-map guard: IQR relative to maxabs
flat_topo_span_ratio_thr = 0.25;    % relaxed flat-map guard: span relative to maxabs
flat_topo_std_ratio_thr = 0.08;     % relaxed flat-map guard: STD relative to maxabs
peak_form_weight = 1.25;            % dominant positive score term for shift-invariant gamma peak form
peak_bonus_weight = 0.00;           % legacy peak bonus weight (ranking disabled; retained for diagnostics/rescue)
peak_bonus_rescue_min = 0.55;       % minimum peak-clarity to rescue flat-only exclusions
peak_bonus_min_peaks = 1;           % minimum number of detected peaks for rescue
peak_form_shift_max_hz = 25;        % max template-shift alignment range for shift-invariant matching
peak_form_single_widths = [3 5 7 9];
peak_form_double_separations = [8 12 16 20];
peak_form_double_widths = [3 4 5 6];
peak_form_min_trough_depth = 0.10;  % minimum trough depth for double-peak plausibility
peak_form_min_similarity = 0.20;    % floor similarity before positive score is assigned
peak_form_smooth_n = 3;             % moving-average smoothing for template scoring
peak_form_enable_selfcheck = true;  % synthetic shift-invariance diagnostic
peak_form_selfcheck_centers = [40 50 60 70 80];
peak_form_selfcheck_width_hz = 5;
peak_form_selfcheck_tol = 0.03;
hard_leak_severity_mult = 1.15;     % leak must exceed adaptive threshold by this factor for hard reject
hard_emg_severity_mult = 1.10;      % EMG score must exceed adaptive threshold by this factor for hard reject
crosswin_w_topo = 0.50;             % cross-window matching weight for topography similarity
crosswin_w_spec = 0.35;             % cross-window matching weight for spectrum similarity
crosswin_w_feat = 0.15;             % cross-window matching weight for scalar-feature similarity
crosswin_match_min_conf = 0.65;     % minimum confidence required to accept a cross-window match
crosswin_consistency_min = 0.65;    % minimum mean confidence to treat a cross-window group as stable
crosswin_consistency_bonus = 0.12;  % additive score bonus for stable occipital-like cross-window groups

random_seed = 13;                    % reproducible randomization

if peak_form_enable_selfcheck
    [peak_form_shift_ok, peak_form_shift_stats] = validate_peak_form_shift_invariance( ...
        scan_freqs, analysis_freq_range, peak_form_shift_max_hz, peak_form_single_widths, ...
        peak_form_double_widths, peak_form_double_separations, peak_form_min_trough_depth, ...
        peak_form_min_similarity, peak_form_smooth_n, peak_form_selfcheck_centers, ...
        peak_form_selfcheck_width_hz, peak_form_selfcheck_tol);
    if ~peak_form_shift_ok
        warning('Peak-form shift-invariance check indicates center-frequency bias (max spread %.3f).', ...
            peak_form_shift_stats.single_spread);
    end
end

% Run mode: 'full' = full pipeline; 'component_check' = GED + component selection only (Phase 1)
run_mode = 'full';

% Three-way component selection config (ordered for plotting/metrics):
% raw -> top component -> combined artifact-screened weighted GED
benchmark_methods = {'raw', 'ged_top_eig', 'ged_combined_artifact_weighted'};
nBenchmarkMethods = numel(benchmark_methods);
raw_reference_definition = 'posterior_roi';  % locked reference for raw reference branch
compute_detectability = true;
compute_separation = true;
compute_reliability = true;

% Validation settings (reviewer-facing)
simulation_validation_enable = true;
simulation_reps = 40;
simulation_snr_levels = [0.4 0.8 1.2];
simulation_depth_levels = [0.6 1.0 1.4];
simulation_overlap_levels = [0.15 0.35 0.55];
simulation_artifact_levels = [0.0 0.25 0.5];

% Detrending parameters for power-ratio spectrum
poly_order = 2;                    % use linear-domain quadratic detrending (legacy-compatible)
detrend_edge_exclude_n = 0;        % fit on full 10-110 Hz support before trimming to 30-90 Hz
detrend_in_log = false;            % linear-domain detrending
detrend_flat_edges = false;        % no edge-flattening ramp
peak_min_prom_frac = 0.15;         % MinPeakProminence as fraction of current-spectrum max
peak_min_distance_hz = 5;          % MinPeakDistance in Hz
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
    gcp_root_path = 'W:\Students\Arne\GCP';
else
    gcp_root_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP';
end
fig_save_dir_ged = fullfile(gcp_root_path, 'figures', 'eeg', 'ged');
fig_save_dir_component_comparison = fullfile(fig_save_dir_ged, 'component_comparison');
fig_save_dir_component_selection = fullfile(fig_save_dir_ged, 'component_selection');
fig_save_dir_subj = fullfile(fig_save_dir_ged, 'subj');
fig_save_dir_emg_exclusion = fig_save_dir_component_selection;
if ~exist(fig_save_dir_ged, 'dir'), mkdir(fig_save_dir_ged); end
if ~exist(fig_save_dir_component_comparison, 'dir'), mkdir(fig_save_dir_component_comparison); end
if ~exist(fig_save_dir_component_selection, 'dir'), mkdir(fig_save_dir_component_selection); end
if ~exist(fig_save_dir_subj, 'dir'), mkdir(fig_save_dir_subj); end
comp_sel_save_dir = fig_save_dir_component_selection;

%% Preallocate storage
all_trial_powratio     = cell(4, nSubj);
all_trial_powratio_fullscan = cell(4, nSubj);
all_trial_powratio_early = cell(4, nSubj);
all_trial_powratio_late  = cell(4, nSubj);
all_trial_peaks_single = cell(4, nSubj);
all_trial_peaks_low    = cell(4, nSubj);
all_trial_peaks_high   = cell(4, nSubj);
all_trial_peaks_single_early = cell(4, nSubj);
all_trial_peaks_single_late  = cell(4, nSubj);
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
all_crosswin_match_stats = cell(1, nSubj);
warning_log_by_subj = cell(nSubj, 1);

all_trial_powratio_bench = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_bench_early = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_bench_late  = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_dt_bench = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_dt_bench_early = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_dt_bench_late  = cell(nBenchmarkMethods, 4, nSubj);
all_trial_powratio_components_full  = cell(4, nSubj);
all_trial_powratio_components_early = cell(4, nSubj);
all_trial_powratio_components_late  = cell(4, nSubj);
benchmark_metric_detectability = nan(nBenchmarkMethods, 4, nSubj);
benchmark_metric_prominence = nan(nBenchmarkMethods, 4, nSubj);
benchmark_metric_separation_slope = nan(nBenchmarkMethods, nSubj);
benchmark_metric_separation_delta = nan(nBenchmarkMethods, nSubj);
benchmark_metric_reliability_trialcv = nan(nBenchmarkMethods, 4, nSubj);
benchmark_metric_reliability_subjspread = nan(nBenchmarkMethods, 4);
simulation_validation_results = struct();
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
    temp_mask = cellfun(@(l) ~isempty(regexp(l, '^(T|TP|FT)', 'once')), dataEEG_c25.label);
    temp_idx = find(temp_mask);
    if isempty(temp_idx)
        temp_idx = setdiff(1:nChans, union(occ_idx, front_idx));
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
    searchOccFrontRatio_full = []; searchOccFrontRatio_early = []; searchOccFrontRatio_late = [];
    searchGammaEvidence_full = []; searchGammaEvidence_early = []; searchGammaEvidence_late = [];
    searchTopos_full = []; searchTopos_early = []; searchTopos_late = [];
    searchMeanPrSpectrum_full = []; searchMeanPrSpectrum_early = []; searchMeanPrSpectrum_late = [];
    gamma_vec_full = []; gamma_vec_early = []; gamma_vec_late = [];
    occipital_evidence_full = []; occipital_evidence_early = []; occipital_evidence_late = [];
    emg_artifact_score_full = []; emg_artifact_score_early = []; emg_artifact_score_late = [];
    searchEmgClass_full = {}; searchEmgClass_early = {}; searchEmgClass_late = {};
    unknown_proxy_full = []; unknown_proxy_early = []; unknown_proxy_late = [];
    hard_eligible_full = []; hard_eligible_early = []; hard_eligible_late = [];
    hard_reject_full = []; hard_reject_early = []; hard_reject_late = [];
    soft_warn_full = []; soft_warn_early = []; soft_warn_late = [];
    rejection_flags_full = struct(); rejection_flags_early = struct(); rejection_flags_late = struct();
    soft_warn_flags_full = struct(); soft_warn_flags_early = struct(); soft_warn_flags_late = struct();
    adaptive_thr_full = struct(); adaptive_thr_early = struct(); adaptive_thr_late = struct();
    component_signature_full = struct(); component_signature_early = struct(); component_signature_late = struct();
    candidate_table_full = struct(); candidate_table_early = struct(); candidate_table_late = struct();
    crosswin_id_full = []; crosswin_id_early = []; crosswin_id_late = [];
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
        searchGammaEvidence = nan(nSearch, 1);
        searchFrontLeak = nan(nSearch, 1);
        searchTempStrength = nan(nSearch, 1);
        searchTempLeak = nan(nSearch, 1);
        searchLineHarmRatio = nan(nSearch, 1);
        searchStationarityCV = nan(nSearch, 1);
        searchBurstRatio = nan(nSearch, 1);
        searchHFSlope = nan(nSearch, 1);
        searchCondLockRho = nan(nSearch, 1);
        searchProxyUnknown = false(nSearch, 1);
        searchOccipitalEvidence = nan(nSearch, 1);
        searchEmgArtifactScore = nan(nSearch, 1);
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
            ratio_ci = occ_strength / max(front_strength, eps);
            front_leak_ci = front_strength / max(occ_strength, eps);
        else
            front_strength = 0;
            ratio_ci = Inf;
            front_leak_ci = 0;
        end
        temp_strength = mean(abs(topo_ci(temp_idx)));
        temp_leak_ci = temp_strength / max(occ_strength, eps);
        % Pooled gamma evidence (all conditions combined).
        gamma_ev_ci = log( max((w_ci' * covStim_reg * w_ci), eps) / ...
                           max((w_ci' * covBase_reg * w_ci), eps) );
        % Lightweight artifact proxies from trial-level component spectra.
        proxy_ci = estimate_component_artifact_proxies( ...
            w_ci, dat_per_cond, stim_windows{w}, baseline_window, scan_freqs, scan_width, artifact_proxy_max_trials);

        searchFilters(:, ci) = w_ci;
        searchTopos(:, ci) = topo_ci;
        searchCorrs(ci) = r_ci;
        searchOccStrength(ci) = occ_strength;
        searchFrontStrength(ci) = front_strength;
        searchOccFrontRatio(ci) = ratio_ci;
        searchGammaEvidence(ci) = gamma_ev_ci;
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
    end
    W_top = searchFilters(:, 1);

    % Candidate metrics (selection = eig > 1 + explicit artifact rejection; no fixed gamma threshold).
    corr_vec = searchCorrs;
    ratio_vec = searchOccFrontRatio;
    gamma_vec = searchGammaEvidence;
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
        ~isfinite(burst_vec) | ~isfinite(hf_slope_vec) | ~isfinite(condlock_vec);
    finite_metrics = isfinite(corr_vec) & isfinite(ratio_vec) & isfinite(gamma_vec) & ...
        isfinite(eval_raw_vec) & isfinite(leak_vec) & isfinite(temp_leak_vec) & ...
        isfinite(lineharm_vec) & isfinite(stationarity_vec) & isfinite(burst_vec) & ...
        isfinite(hf_slope_vec) & isfinite(condlock_vec);
    [peak_bonus_vec, peak_count_vec] = compute_peak_bonus_from_spectra( ...
        searchMeanPrSpectrum, scan_freqs, analysis_freq_range, poly_order, detrend_edge_exclude_n, ...
        detrend_in_log, detrend_flat_edges, peak_min_prom_frac, peak_min_distance_hz);
    [peak_form_score_vec, peak_form_mode_vec, peak_form_diag] = compute_peak_form_template_score_from_spectra( ...
        searchMeanPrSpectrum, scan_freqs, analysis_freq_range, poly_order, detrend_edge_exclude_n, ...
        detrend_in_log, detrend_flat_edges, peak_form_shift_max_hz, peak_form_single_widths, ...
        peak_form_double_widths, peak_form_double_separations, peak_form_min_trough_depth, ...
        peak_form_min_similarity, peak_form_smooth_n);
    occipital_evidence = 0.40 * normalize_robust(corr_vec) + ...
        0.30 * normalize_robust(ratio_vec) + ...
        0.30 * normalize_robust(condlock_vec);
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
    searchOccipitalEvidence = occipital_evidence;
    searchEmgArtifactScore = emg_artifact_score;
    hard_eligible_raw = finite_metrics & (eval_raw_vec >= adaptive_thr.min_eigval);
    flat_topomap_mask = false(nSearch, 1);
    if exclude_flat_topomap
        flat_topomap_mask = detect_flat_topomap_outliers( ...
            searchTopos, flat_topo_mad_ratio_thr, flat_topo_iqr_ratio_thr, ...
            flat_topo_span_ratio_thr, flat_topo_std_ratio_thr);
    end
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
    occ_minus_emg = occipital_evidence - emg_artifact_score;
    fail_occ_margin = finite_metrics & (emg_artifact_score >= adaptive_thr.emg_class_thr) & ...
        (occ_minus_emg < adaptive_thr.min_occ_margin);
    peak_rescue_mask = false(nSearch, 1);
    if any(flat_topomap_mask)
        non_flat_fail = fail_front_leak | fail_temp_leak | fail_corr | fail_lineharm | ...
            fail_stationarity | fail_burst | fail_hf_slope | fail_condlock | fail_emg_score | ...
            fail_occ_margin | unknown_proxy_vec;
        peak_rescue_mask = flat_topomap_mask & finite_metrics & ~non_flat_fail & ...
            (peak_bonus_vec >= peak_bonus_rescue_min) & (peak_count_vec >= peak_bonus_min_peaks) & ...
            (occ_minus_emg >= adaptive_thr.min_occ_margin) & (corr_vec >= adaptive_thr.min_corr);
        flat_topomap_mask(peak_rescue_mask) = false;
    end
    hard_reject_flags = severe_front_leak | severe_temp_leak | fail_lineharm | fail_hf_slope | ...
        severe_emg_score | unknown_proxy_vec | flat_topomap_mask;
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
    artifact_flags = hard_reject_flags;
    rejection_flags = struct( ...
        'unknown_proxy', unknown_proxy_vec, ...
        'flat_topomap', flat_topomap_mask, ...
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
        'peak_rescue', peak_rescue_mask);
    for ci = 1:nSearch
        if unknown_proxy_vec(ci)
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
            searchEmgClass{ci} = 'mixed';
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
    hard_eligible = hard_eligible_raw & ~artifact_flags;
    hard_eligible(dominant_outlier_mask | flat_topomap_mask) = false;
    searchScores = eval_raw_vec + peak_form_weight * peak_form_score_vec + peak_bonus_weight * peak_bonus_vec;
    searchScores(~finite_metrics) = -Inf;
    searchScores(~hard_eligible) = -Inf;
    if ~any(isfinite(searchScores))
        msg = sprintf(['No components met eig/artifact criteria for subject %s. ', ...
                       'Falling back to unconstrained eigenvalue ranking.'], subjects{subj});
        top_fail_idx = NaN;
        top_fail_eig = NaN;
        top_fail_corr = NaN;
        top_fail_ratio = NaN;
        top_fail_gamma = NaN;
        top_fail_pass_eig = false;
        top_fail_artifact = false;
        finite_idx = find(finite_metrics);
        if ~isempty(finite_idx)
            [~, fail_ord] = sort(evals_sorted(finite_idx), 'descend');
            top_fail_idx = finite_idx(fail_ord(1));
            top_fail_eig = evals_sorted(top_fail_idx);
            top_fail_corr = corr_vec(top_fail_idx);
            top_fail_ratio = ratio_vec(top_fail_idx);
            top_fail_gamma = gamma_vec(top_fail_idx);
            top_fail_pass_eig = top_fail_eig >= adaptive_thr.min_eigval;
            top_fail_artifact = artifact_flags(top_fail_idx);
        end
        hard_metrics = struct( ...
            'n_search', nSearch, ...
            'n_finite_metrics', sum(finite_metrics), ...
            'n_pass_eig', sum(finite_metrics & (eval_raw_vec >= adaptive_thr.min_eigval)), ...
            'n_artifact_flagged', sum(artifact_flags), ...
            'n_unknown_high_risk', sum(unknown_proxy_vec), ...
            'n_pass_all_raw', sum(hard_eligible_raw), ...
            'n_excluded_dominant_outlier', sum(dominant_outlier_mask), ...
            'top_fail_idx', top_fail_idx, ...
            'top_fail_eig', top_fail_eig, ...
            'top_fail_corr', top_fail_corr, ...
            'top_fail_ratio', top_fail_ratio, ...
            'top_fail_gamma', top_fail_gamma, ...
            'top_fail_pass_eig', top_fail_pass_eig, ...
            'top_fail_artifact', top_fail_artifact, ...
            'thr_eig', adaptive_thr.min_eigval, ...
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
        searchScores = eval_raw_vec + peak_form_weight * peak_form_score_vec + peak_bonus_weight * peak_bonus_vec;
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
    candidate_table.gamma = gamma_vec;
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
    candidate_table.peak_form_best_single_similarity = peak_form_diag.best_single_similarity;
    candidate_table.peak_form_best_double_similarity = peak_form_diag.best_double_similarity;
    candidate_table.peak_form_similarity_raw = peak_form_diag.best_similarity_raw;
    candidate_table.peak_form_best_shift_hz = peak_form_diag.best_shift_hz;
    candidate_table.peak_form_best_width_hz = peak_form_diag.best_width_hz;
    candidate_table.peak_form_best_separation_hz = peak_form_diag.best_separation_hz;
    candidate_table.peak_form_trough_depth = peak_form_diag.best_trough_depth;
    candidate_table.occipital_evidence = occipital_evidence;
    candidate_table.emg_artifact_score = emg_artifact_score;
    candidate_table.emg_class = searchEmgClass;
    candidate_table.unknown_high_risk = unknown_proxy_vec;
    candidate_table.score = searchScores;
    candidate_table.score_base = searchScores;
    candidate_table.score_final = searchScores;
    candidate_table.consistency_bonus = zeros(nSearch, 1);
    candidate_table.artifact_flag = artifact_flags;
    candidate_table.hard_reject = hard_reject_flags;
    candidate_table.soft_warn = soft_warn_any;
    candidate_table.reject_reason = compute_primary_rejection_reason(rejection_flags);
    candidate_table.fail_unknown_proxy = rejection_flags.unknown_proxy;
    candidate_table.fail_flat_topomap = rejection_flags.flat_topomap;
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
    candidate_table.soft_warn_front_leak = soft_warn_flags.front_leak;
    candidate_table.soft_warn_temp_leak = soft_warn_flags.temp_leak;
    candidate_table.soft_warn_corr = soft_warn_flags.corr;
    candidate_table.soft_warn_stationarity = soft_warn_flags.stationarity;
    candidate_table.soft_warn_burst = soft_warn_flags.burst;
    candidate_table.soft_warn_condlock = soft_warn_flags.condlock;
    candidate_table.soft_warn_emg_score = soft_warn_flags.emg_score;
    candidate_table.soft_warn_occ_margin = soft_warn_flags.occ_margin;
    candidate_table.peak_rescued_from_flat = rejection_flags.peak_rescue;
    candidate_table.hard_eligible_raw = hard_eligible_raw;
    candidate_table.dominant_outlier = dominant_outlier_mask;
    candidate_table.hard_eligible = hard_eligible;
    candidate_table.crosswin_id = nan(nSearch, 1);
    candidate_table.crosswin_match_conf = nan(nSearch, 1);
    candidate_table.thr_min_eigval = repmat(adaptive_thr.min_eigval, nSearch, 1);
    candidate_table.thr_occ_class = repmat(adaptive_thr.occ_class_thr, nSearch, 1);
    candidate_table.thr_emg_class = repmat(adaptive_thr.emg_class_thr, nSearch, 1);
    candidate_table.thr_occ_margin = repmat(adaptive_thr.min_occ_margin, nSearch, 1);
    component_signature = build_component_signature_block( ...
        searchTopos, searchMeanPrSpectrum, corr_vec, ratio_vec, condlock_vec, ...
        lineharm_vec, hf_slope_vec, leak_vec, temp_leak_vec, eval_raw_vec);

    combined_idx = find(hard_eligible & isfinite(searchScores));
    if isempty(combined_idx)
        msg = sprintf(['No artifact-screened finite components available for subject %s. ', ...
            'Falling back to single best component (unconstrained ranking).'], subjects{subj});
        warning_log_subj = append_subject_warning(warning_log_subj, subjects{subj}, 'NO_HARD_COMPONENTS', msg, ...
            struct('n_hard_eligible', sum(hard_eligible), 'n_finite_scores', sum(isfinite(searchScores)), ...
            'fallback_idx', bestIdx));
        combined_idx = bestIdx;
    end
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

        % Store per-window filters for Phase 2
        if w == 1
            searchFilters_full = searchFilters;
            W_top_full = searchFilters(:, 1);
            W_combined_full = searchFilters(:, selected_idx);
            selected_idx_full = selected_idx;
            w_combined_full = selected_weights(:)';
            evals_sorted_full = evals_sorted;
            searchCorrs_full = searchCorrs;
            searchOccFrontRatio_full = searchOccFrontRatio;
            searchGammaEvidence_full = searchGammaEvidence;
            searchTopos_full = searchTopos;
            searchMeanPrSpectrum_full = searchMeanPrSpectrum;
            gamma_vec_full = gamma_vec;
            occipital_evidence_full = occipital_evidence;
            emg_artifact_score_full = emg_artifact_score;
            searchEmgClass_full = searchEmgClass;
            unknown_proxy_full = unknown_proxy_vec;
            hard_eligible_full = hard_eligible;
            hard_reject_full = artifact_flags;
            soft_warn_full = soft_warn_flags.any;
            rejection_flags_full = rejection_flags;
            soft_warn_flags_full = soft_warn_flags;
            adaptive_thr_full = adaptive_thr;
            component_signature_full = component_signature;
            candidate_table_full = candidate_table;
        elseif w == 2
            searchFilters_early = searchFilters;
            W_top_early = searchFilters(:, 1);
            W_combined_early = searchFilters(:, selected_idx);
            selected_idx_early = selected_idx;
            w_combined_early = selected_weights(:)';
            evals_sorted_early = evals_sorted;
            searchCorrs_early = searchCorrs;
            searchOccFrontRatio_early = searchOccFrontRatio;
            searchGammaEvidence_early = searchGammaEvidence;
            searchTopos_early = searchTopos;
            searchMeanPrSpectrum_early = searchMeanPrSpectrum;
            gamma_vec_early = gamma_vec;
            occipital_evidence_early = occipital_evidence;
            emg_artifact_score_early = emg_artifact_score;
            searchEmgClass_early = searchEmgClass;
            unknown_proxy_early = unknown_proxy_vec;
            hard_eligible_early = hard_eligible;
            hard_reject_early = artifact_flags;
            soft_warn_early = soft_warn_flags.any;
            rejection_flags_early = rejection_flags;
            soft_warn_flags_early = soft_warn_flags;
            adaptive_thr_early = adaptive_thr;
            component_signature_early = component_signature;
            candidate_table_early = candidate_table;
        else
            searchFilters_late = searchFilters;
            W_top_late = searchFilters(:, 1);
            W_combined_late = searchFilters(:, selected_idx);
            selected_idx_late = selected_idx;
            w_combined_late = selected_weights(:)';
            evals_sorted_late = evals_sorted;
            searchCorrs_late = searchCorrs;
            searchOccFrontRatio_late = searchOccFrontRatio;
            searchGammaEvidence_late = searchGammaEvidence;
            searchTopos_late = searchTopos;
            searchMeanPrSpectrum_late = searchMeanPrSpectrum;
            gamma_vec_late = gamma_vec;
            occipital_evidence_late = occipital_evidence;
            emg_artifact_score_late = emg_artifact_score;
            searchEmgClass_late = searchEmgClass;
            unknown_proxy_late = unknown_proxy_vec;
            hard_eligible_late = hard_eligible;
            hard_reject_late = artifact_flags;
            soft_warn_late = soft_warn_flags.any;
            rejection_flags_late = rejection_flags;
            soft_warn_flags_late = soft_warn_flags;
            adaptive_thr_late = adaptive_thr;
            component_signature_late = component_signature;
            candidate_table_late = candidate_table;
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
            'selection_mode', 'eig_artifact_screened_weighted', ...
            'selected_idx', selected_idx, ...
            'selected_weights', selected_weights, ...
            'candidate_table', candidate_table, ...
            'component_signature', component_signature, ...
            'best_idx', bestIdx, ...
            'best_score', bestScore, ...
            'best_corr', bestCorr, ...
            'best_ratio', bestRatio, ...
            'best_gamma', bestGamma, ...
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

    % Combined-component topoplots: selected components (ordered by eigenvalue).
    [~, score_ord] = sort(evals_sorted(selected_idx), 'descend');
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
        'VerticalAlignment', 'middle', 'FontSize', 18, 'FontWeight', 'bold');
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
            g_log = searchGammaEvidence(comp_rank);
            g_pct = 100 * (exp(g_log) - 1);   % gamma amplification as % change from baseline
            title(sprintf('C%d:\\lambda=%.2f,r=%.2f,ratio=%.2f,g=%.0f%%', ...
                comp_rank, evals_sorted(comp_rank), searchCorrs(comp_rank), ...
                searchOccFrontRatio(comp_rank), g_pct), 'FontSize', 6);
        end
    end
        saveas(fig_post, fullfile(comp_sel_save_dir, sprintf('GCP_eeg_GED_component_selection_subj%s_%s.png', subjects{subj}, win_names{w})));
        close(fig_post);
        % EMG exclusion diagnostics are generated after all windows are processed,
        % once cross-window IDs and consistency updates are available.
    end

    match_cfg = struct( ...
        'w_topo', crosswin_w_topo, ...
        'w_spec', crosswin_w_spec, ...
        'w_feat', crosswin_w_feat, ...
        'min_conf', crosswin_match_min_conf);
    crosswin_match = match_components_full_anchor( ...
        component_signature_full, component_signature_early, component_signature_late, match_cfg);
    crosswin_id_full = crosswin_match.id_full;
    crosswin_id_early = crosswin_match.id_early;
    crosswin_id_late = crosswin_match.id_late;

    [group_consistent_map, group_occipital_map, group_member_count, group_mean_conf] = ...
        compute_crosswindow_group_consistency( ...
            crosswin_match, candidate_table_full, candidate_table_early, candidate_table_late, ...
            crosswin_consistency_min);

    [candidate_table_full, selected_idx_full, w_combined_full] = ...
        apply_crosswindow_consistency_bonus(candidate_table_full, crosswin_id_full, ...
            group_consistent_map, group_occipital_map, crosswin_consistency_bonus, max_components_to_combine);
    [candidate_table_early, selected_idx_early, w_combined_early] = ...
        apply_crosswindow_consistency_bonus(candidate_table_early, crosswin_id_early, ...
            group_consistent_map, group_occipital_map, crosswin_consistency_bonus, max_components_to_combine);
    [candidate_table_late, selected_idx_late, w_combined_late] = ...
        apply_crosswindow_consistency_bonus(candidate_table_late, crosswin_id_late, ...
            group_consistent_map, group_occipital_map, crosswin_consistency_bonus, max_components_to_combine);

    if isempty(selected_idx_full)
        selected_idx_full = 1;
        w_combined_full = 1;
    end
    if isempty(selected_idx_early)
        selected_idx_early = 1;
        w_combined_early = 1;
    end
    if isempty(selected_idx_late)
        selected_idx_late = 1;
        w_combined_late = 1;
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
    all_topo_temp{1} = topo_temp_full;
    all_topo_temp{2} = topo_temp_early;
    all_topo_temp{3} = topo_temp_late;
    all_topos{subj} = topo_temp_full;
    all_topos_early{subj} = topo_temp_early;
    all_topos_late{subj} = topo_temp_late;
    all_selected_comp_indices_multi{subj} = selected_idx_full;
    all_selected_comp_weights{subj} = w_combined_full(:)';
    all_selected_comp_idx(subj) = selected_idx_full(1);
    all_selected_comp_corr(subj) = searchCorrs_full(selected_idx_full(1));
    all_selected_comp_eval(subj) = evals_sorted_full(selected_idx_full(1));
    all_eigenvalues(subj) = evals_sorted_full(selected_idx_full(1));

    sel_id_full = crosswin_id_full(selected_idx_full(1));
    sel_id_early = crosswin_id_early(selected_idx_early(1));
    sel_id_late = crosswin_id_late(selected_idx_late(1));
    sanity_report = struct( ...
        'subject', subjects{subj}, ...
        'selected_id_full', sel_id_full, ...
        'selected_id_early', sel_id_early, ...
        'selected_id_late', sel_id_late, ...
        'full_matches_early', isfinite(sel_id_full) && isfinite(sel_id_early) && (sel_id_full == sel_id_early), ...
        'full_matches_late', isfinite(sel_id_full) && isfinite(sel_id_late) && (sel_id_full == sel_id_late), ...
        'early_matches_late', isfinite(sel_id_early) && isfinite(sel_id_late) && (sel_id_early == sel_id_late), ...
        'selected_hard_reject_full', any(candidate_table_full.hard_reject(selected_idx_full)), ...
        'selected_hard_reject_early', any(candidate_table_early.hard_reject(selected_idx_early)), ...
        'selected_hard_reject_late', any(candidate_table_late.hard_reject(selected_idx_late)));
    crosswin_match.sanity_report = sanity_report;
    all_crosswin_match_stats{subj} = crosswin_match;

    candidate_table_full.crosswin_id = crosswin_id_full;
    candidate_table_early.crosswin_id = crosswin_id_early;
    candidate_table_late.crosswin_id = crosswin_id_late;
    candidate_table_full.crosswin_match_conf = crosswin_match.conf_full;
    candidate_table_early.crosswin_match_conf = crosswin_match.conf_early;
    candidate_table_late.crosswin_match_conf = crosswin_match.conf_late;
    candidate_table_full.crosswin_group_consistent = lookup_group_property(crosswin_id_full, group_consistent_map);
    candidate_table_early.crosswin_group_consistent = lookup_group_property(crosswin_id_early, group_consistent_map);
    candidate_table_late.crosswin_group_consistent = lookup_group_property(crosswin_id_late, group_consistent_map);
    candidate_table_full.crosswin_group_occipital = lookup_group_property(crosswin_id_full, group_occipital_map);
    candidate_table_early.crosswin_group_occipital = lookup_group_property(crosswin_id_early, group_occipital_map);
    candidate_table_late.crosswin_group_occipital = lookup_group_property(crosswin_id_late, group_occipital_map);
    candidate_table_full.crosswin_group_members = lookup_group_property(crosswin_id_full, group_member_count);
    candidate_table_early.crosswin_group_members = lookup_group_property(crosswin_id_early, group_member_count);
    candidate_table_late.crosswin_group_members = lookup_group_property(crosswin_id_late, group_member_count);
    candidate_table_full.crosswin_group_mean_conf = lookup_group_property(crosswin_id_full, group_mean_conf);
    candidate_table_early.crosswin_group_mean_conf = lookup_group_property(crosswin_id_early, group_mean_conf);
    candidate_table_late.crosswin_group_mean_conf = lookup_group_property(crosswin_id_late, group_mean_conf);

    all_component_selection_stats_full{subj}.selected_idx = selected_idx_full;
    all_component_selection_stats_full{subj}.selected_weights = w_combined_full;
    all_component_selection_stats_full{subj}.candidate_table = candidate_table_full;
    all_component_selection_stats_full{subj}.crosswin_id = crosswin_id_full;
    all_component_selection_stats_full{subj}.crosswin_match = crosswin_match;
    all_component_selection_stats_full{subj}.crosswin_group_consistent = group_consistent_map;

    all_component_selection_stats_early{subj}.selected_idx = selected_idx_early;
    all_component_selection_stats_early{subj}.selected_weights = w_combined_early;
    all_component_selection_stats_early{subj}.candidate_table = candidate_table_early;
    all_component_selection_stats_early{subj}.crosswin_id = crosswin_id_early;
    all_component_selection_stats_early{subj}.crosswin_match = crosswin_match;
    all_component_selection_stats_early{subj}.crosswin_group_consistent = group_consistent_map;

    all_component_selection_stats_late{subj}.selected_idx = selected_idx_late;
    all_component_selection_stats_late{subj}.selected_weights = w_combined_late;
    all_component_selection_stats_late{subj}.candidate_table = candidate_table_late;
    all_component_selection_stats_late{subj}.crosswin_id = crosswin_id_late;
    all_component_selection_stats_late{subj}.crosswin_match = crosswin_match;
    all_component_selection_stats_late{subj}.crosswin_group_consistent = group_consistent_map;

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
        searchMeanPrSpectrum_full, evals_sorted_full(1:numel(crosswin_id_full)), gamma_vec_full, ...
        occipital_evidence_full, emg_artifact_score_full, searchEmgClass_full, unknown_proxy_full, ...
        candidate_table_full.hard_eligible, hard_reject_full, soft_warn_full, rejection_flags_full, ...
        adaptive_thr_full, cfg_topo, all_topo_labels{subj}, candidate_table_full.peak_form_score, ...
        candidate_table_full.peak_form_mode, crosswin_id_full);
    plot_emg_exclusion_diagnostics( ...
        fig_save_dir_emg_exclusion, subjects{subj}, 'early', scan_freqs, searchTopos_early, ...
        searchMeanPrSpectrum_early, evals_sorted_early(1:numel(crosswin_id_early)), gamma_vec_early, ...
        occipital_evidence_early, emg_artifact_score_early, searchEmgClass_early, unknown_proxy_early, ...
        candidate_table_early.hard_eligible, hard_reject_early, soft_warn_early, rejection_flags_early, ...
        adaptive_thr_early, cfg_topo, all_topo_labels{subj}, candidate_table_early.peak_form_score, ...
        candidate_table_early.peak_form_mode, crosswin_id_early);
    plot_emg_exclusion_diagnostics( ...
        fig_save_dir_emg_exclusion, subjects{subj}, 'late', scan_freqs, searchTopos_late, ...
        searchMeanPrSpectrum_late, evals_sorted_late(1:numel(crosswin_id_late)), gamma_vec_late, ...
        occipital_evidence_late, emg_artifact_score_late, searchEmgClass_late, unknown_proxy_late, ...
        candidate_table_late.hard_eligible, hard_reject_late, soft_warn_late, rejection_flags_late, ...
        adaptive_thr_late, cfg_topo, all_topo_labels{subj}, candidate_table_late.peak_form_score, ...
        candidate_table_late.peak_form_mode, crosswin_id_late);

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
        powratio_methods_full = nan(nBenchmarkMethods, nTrl, nFreqs_detrend);
        powratio_methods_early = nan(nBenchmarkMethods, nTrl, nFreqs_detrend);
        powratio_methods_late = nan(nBenchmarkMethods, nTrl, nFreqs_detrend);
        nSearch_full = size(filters.full.searchFilters, 2);
        nSearch_early = size(filters.early.searchFilters, 2);
        nSearch_late = size(filters.late.searchFilters, 2);
        powratio_components       = nan(nSearch_full, nTrl, nFreqs_detrend);
        powratio_components_early = nan(nSearch_early, nTrl, nFreqs_detrend);
        powratio_components_late  = nan(nSearch_late, nTrl, nFreqs_detrend);

        for fi = 1:nFreqs_detrend
            clc
            fprintf('Subject    %s (%d/%d)\nCondition  %d/4\nFrequency  %d/%d\n', ...
                subjects{subj}, subj, nSubj, cond, fi, nFreqs_detrend);
            cf = scan_freqs_detrend(fi);
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
        powratio_methods_full_analysis = powratio_methods_full(:, :, analysis_freq_mask_detrend);
        powratio_methods_early_analysis = powratio_methods_early(:, :, analysis_freq_mask_detrend);
        powratio_methods_late_analysis = powratio_methods_late(:, :, analysis_freq_mask_detrend);
        powratio_components_analysis = powratio_components(:, :, analysis_freq_mask_detrend);
        powratio_components_early_analysis = powratio_components_early(:, :, analysis_freq_mask_detrend);
        powratio_components_late_analysis = powratio_components_late(:, :, analysis_freq_mask_detrend);

        all_trial_powratio_components_full{cond, subj}  = powratio_components_analysis;
        all_trial_powratio_components_early{cond, subj} = powratio_components_early_analysis;
        all_trial_powratio_components_late{cond, subj}  = powratio_components_late_analysis;

        for mi = 1:nBenchmarkMethods
            pr_m_full = squeeze(powratio_methods_full_analysis(mi, :, :));
            all_trial_powratio_bench{mi, cond, subj} = pr_m_full;
            if ~isempty(pr_m_full)
                dt_m_full = nan(size(pr_m_full));
                for trl = 1:size(pr_m_full, 1)
                    pr_full = squeeze(powratio_methods_full(mi, trl, :))';
                    dt_full = detrend_power_ratio(pr_full, scan_freqs_detrend, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
                    dt_m_full(trl,:) = dt_full(analysis_freq_mask_detrend);
                end
                all_trial_powratio_dt_bench{mi, cond, subj} = dt_m_full;
            end

            pr_m_early = squeeze(powratio_methods_early_analysis(mi, :, :));
            all_trial_powratio_bench_early{mi, cond, subj} = pr_m_early;
            if ~isempty(pr_m_early)
                dt_m_early = nan(size(pr_m_early));
                for trl = 1:size(pr_m_early, 1)
                    pr_full = squeeze(powratio_methods_early(mi, trl, :))';
                    dt_full = detrend_power_ratio(pr_full, scan_freqs_detrend, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
                    dt_m_early(trl,:) = dt_full(analysis_freq_mask_detrend);
                end
                all_trial_powratio_dt_bench_early{mi, cond, subj} = dt_m_early;
            end

            pr_m_late = squeeze(powratio_methods_late_analysis(mi, :, :));
            all_trial_powratio_bench_late{mi, cond, subj} = pr_m_late;
            if ~isempty(pr_m_late)
                dt_m_late = nan(size(pr_m_late));
                for trl = 1:size(pr_m_late, 1)
                    pr_full = squeeze(powratio_methods_late(mi, trl, :))';
                    dt_full = detrend_power_ratio(pr_full, scan_freqs_detrend, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
                    dt_m_late(trl,:) = dt_full(analysis_freq_mask_detrend);
                end
                all_trial_powratio_dt_bench_late{mi, cond, subj} = dt_m_late;
            end
        end

        % Keep full-window outputs based on weighted combined GED branch.
        powratio_trials_fullscan = squeeze(powratio_methods_full(3, :, :));
        powratio_trials_early_fullscan = squeeze(powratio_methods_early(3, :, :));
        powratio_trials_late_fullscan = squeeze(powratio_methods_late(3, :, :));
        powratio_trials_full = squeeze(powratio_methods_full_analysis(3, :, :));
        powratio_trials_early = squeeze(powratio_methods_early_analysis(3, :, :));
        powratio_trials_late = squeeze(powratio_methods_late_analysis(3, :, :));
        all_trial_powratio_fullscan{cond, subj} = powratio_trials_fullscan;
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
        trl_peaks_low    = nan(nTrl, 1);
        trl_peaks_high   = nan(nTrl, 1);
        trl_centroid     = nan(nTrl, 1);

        for trl = 1:nTrl
            pr_full = powratio_trials_fullscan(trl, :);
            if all(isnan(pr_full)), continue; end

            pr_dt_full = detrend_power_ratio(pr_full, scan_freqs_detrend, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
            pr_dt = pr_dt_full(analysis_freq_mask_detrend);
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

            % Dual-peak: hard 50 Hz boundary.
            [pks_all, locs_all] = findpeaks(pr_dt_smooth, scan_freqs, ...
                'MinPeakDistance', peak_min_distance_hz);
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

        % Time-split single-peak summaries.
        trl_peaks_single_early = detect_single_peaks_from_powratio_fullscan( ...
            powratio_trials_early_fullscan, scan_freqs_detrend, analysis_freq_mask_detrend, poly_order, detrend_edge_exclude_n, ...
            detrend_in_log, detrend_flat_edges, peak_min_prom_frac, peak_min_distance_hz);
        trl_peaks_single_late = detect_single_peaks_from_powratio_fullscan( ...
            powratio_trials_late_fullscan, scan_freqs_detrend, analysis_freq_mask_detrend, poly_order, detrend_edge_exclude_n, ...
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

    % Pre-compute detrended matrices (fit over 10-110 Hz, keep 30-90 Hz).
    pr_dt_mats = cell(1, 4);
    for cond = 1:4
        pr_mat_full = all_trial_powratio_fullscan{cond, subj};
        if ~isempty(pr_mat_full)
            nTrl = size(pr_mat_full, 1);
            dt = nan(nTrl, nFreqs);
            for trl = 1:nTrl
                dt_full = detrend_power_ratio(pr_mat_full(trl,:), scan_freqs_detrend, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
                dt(trl,:) = dt_full(analysis_freq_mask_detrend);
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

    % --- Row 2: Mean trial-level spectrum with dual-peak markers ---
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
            ylim([-row2_abs row2_abs]);
        end
        xlabel('Freq [Hz]'); ylabel('\Delta PR');
        det_lo = all_trial_detrate_low(cond, subj);
        det_hi = all_trial_detrate_high(cond, subj);
        title(sprintf('%s Dual (L:%.0f%% H:%.0f%%)', condLabels{cond}, det_lo*100, det_hi*100), 'FontSize', 10);
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
    fig_bench_subj = figure('Position', [0 0 1512/2 982], 'Color', 'w');
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
    all_trial_median_low, ...
    all_trial_median_high, ...
    all_trial_median_gap, ...
    squeeze(benchmark_metric_prominence(3, :, :)), ...
    squeeze(benchmark_metric_reliability_trialcv(3, :, :)), ...
    all_trial_median_centroid, ...
    all_trial_gamma_power};
summary_names = {'Single median [Hz]', 'Low median [Hz]', 'High median [Hz]', ...
    'H-L separation [Hz]', 'Prominence', 'Trial CV', ...
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
        ylim([35 50]);
    elseif mi == 3
        ylim([50 75]);
    elseif mi == 4
        ylim([22.5 36]);
    elseif mi == 5
        ylim([0 10]);
    elseif mi == 7
        ylim([50 60]);
    elseif mi == 8
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

dat_power = (all_trial_gamma_power - 1) * 100; % mean trial-level stim/base gamma power percentage change
dat_power_plot = dat_power;

% Figure-only outlier suppression (per condition): Tukey-style extreme outliers.
for c = 1:4
    vals = dat_power_plot(c, :);
    vals = vals(isfinite(vals));
    if numel(vals) < 4
        continue;
    end
    q = quantile(vals, [0.25 0.75]);
    iqr_val = q(2) - q(1);
    lo_thr = q(1) - 3 * iqr_val;
    hi_thr = q(2) + 3 * iqr_val;
    outlier_mask = dat_power_plot(c, :) < lo_thr | dat_power_plot(c, :) > hi_thr;
    dat_power_plot(c, outlier_mask) = NaN;
end

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
ylabel('Gamma Power Change from Baseline [%]');
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
        pr_mat_full = all_trial_powratio_fullscan{cond, s};
        if isempty(pr_mat_full)
            continue;
        end
        nTrl = size(pr_mat_full, 1);
        pr_dt_mat = nan(nTrl, nFreqs);
        for trl = 1:nTrl
            pr_dt_full = detrend_power_ratio(pr_mat_full(trl,:), scan_freqs_detrend, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
            pr_dt_mat(trl,:) = pr_dt_full(analysis_freq_mask_detrend);
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
        pr_mat_full = all_trial_powratio_fullscan{cond, s};
        if ~isempty(pr_mat_full)
            nTrl = size(pr_mat_full, 1);
            pr_dt_mat = nan(nTrl, nFreqs);
            for trl = 1:nTrl
                pr_dt_full = detrend_power_ratio(pr_mat_full(trl,:), scan_freqs_detrend, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
                pr_dt_mat(trl,:) = pr_dt_full(analysis_freq_mask_detrend);
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
            md_lo = all_trial_median_low(cond, s);
            if ~isnan(md_lo)
                xline(md_lo, ':', 'Color', [0 0 0.7], 'LineWidth', 1.1);
            end
            md_hi = all_trial_median_high(cond, s);
            if ~isnan(md_hi)
                xline(md_hi, ':', 'Color', [0.7 0 0], 'LineWidth', 1.1);
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
%  DETECTION RATE FIGURE (single + dual peaks)
%  ====================================================================
fig_det = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Trial-Level Peak Detection Rate (mean across subjects)', ...
    'FontSize', 18, 'FontWeight', 'bold');

det_data = {all_trial_detrate_single, all_trial_detrate_low, all_trial_detrate_high};
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

saveas(fig_det, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_detection_rate.png'));

if simulation_validation_enable
    fprintf('\n============================================================\n');
    fprintf('Ground-truth simulation validation (no CV)\n');
    fprintf('============================================================\n');
    simulation_validation_results = run_ground_truth_simulation_validation( ...
        headmodel, simulation_reps, simulation_snr_levels, simulation_depth_levels, ...
        simulation_overlap_levels, simulation_artifact_levels, min_eigval_hard, ...
        max_frontleak_hard, max_templeak_hard);
end

%% Save results
save_path = fullfile(gcp_root_path, 'data', 'features', 'GCP_eeg_GED.mat');
save(save_path, ...
    'all_trial_powratio', ...
    'all_trial_powratio_fullscan', ...
    'all_trial_powratio_early', 'all_trial_powratio_late', ...
    'all_trial_peaks_single', 'all_trial_peaks_low', 'all_trial_peaks_high', 'all_trial_centroid', ...
    'all_trial_peaks_single_early', 'all_trial_peaks_single_late', ...
    'all_trial_mean_single', 'all_trial_median_single', ...
    'all_trial_mean_low', 'all_trial_median_low', ...
    'all_trial_mean_high', 'all_trial_median_high', 'all_trial_median_gap', ...
    'all_trial_mean_single_early', 'all_trial_median_single_early', ...
    'all_trial_mean_single_late', 'all_trial_median_single_late', ...
    'all_trial_mean_centroid', 'all_trial_median_centroid', ...
    'all_trial_detrate_single', 'all_trial_detrate_low', 'all_trial_detrate_high', 'all_trial_detrate_single_early', 'all_trial_detrate_single_late', ...
    'all_trial_detrate_centroid', 'all_trial_gamma_power', 'all_trial_gamma_power_early', 'all_trial_gamma_power_late', ...
    'all_topos', 'all_topos_early', 'all_topos_late', 'all_topo_labels', 'all_eigenvalues', ...
    'all_selected_comp_idx', 'all_selected_comp_corr', 'all_selected_comp_eval', ...
    'all_selected_comp_indices_multi', 'all_selected_comp_weights', ...
    'all_component_selection_stats_full', 'all_component_selection_stats_early', 'all_component_selection_stats_late', ...
    'all_crosswin_match_stats', 'warning_log', ...
    'all_trial_powratio_bench', 'all_trial_powratio_dt_bench', ...
    'all_trial_powratio_components_full', 'all_trial_powratio_components_early', 'all_trial_powratio_components_late', ...
    'all_trial_powratio_bench_early', 'all_trial_powratio_bench_late', ...
    'all_trial_powratio_dt_bench_early', 'all_trial_powratio_dt_bench_late', ...
    'benchmark_methods', 'raw_reference_definition', ...
    'benchmark_metric_detectability', 'benchmark_metric_prominence', ...
    'benchmark_metric_separation_slope', 'benchmark_metric_separation_delta', ...
    'benchmark_metric_reliability_trialcv', 'benchmark_metric_reliability_subjspread', ...
    'primary_slope_stats', 'primary_delta_stats', ...
    'simulation_validation_results', ...
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

function trl_peaks_single = detect_single_peaks_from_powratio_fullscan(powratio_trials_fullscan, scan_freqs_full, analysis_mask, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges, peak_min_prom_frac, peak_min_distance_hz)
nTrl = size(powratio_trials_fullscan, 1);
trl_peaks_single = nan(nTrl, 1);
scan_freqs_analysis = scan_freqs_full(analysis_mask);
for trl = 1:nTrl
    pr_full = powratio_trials_fullscan(trl, :);
    if all(~isfinite(pr_full))
        continue;
    end
    pr_dt_full = detrend_power_ratio(pr_full, scan_freqs_full, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
    pr_dt = pr_dt_full(analysis_mask);
    pr_dt_smooth = movmean(pr_dt, 5);
    mprom = max(0, max(pr_dt_smooth) * peak_min_prom_frac);
    [pks, locs] = findpeaks(pr_dt_smooth, scan_freqs_analysis, ...
        'MinPeakProminence', mprom, ...
        'MinPeakDistance', peak_min_distance_hz);
    if ~isempty(pks)
        [~, best_pk] = max(pks);
        trl_peaks_single(trl) = locs(best_pk);
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

function cmap = emg_white_yellow_orange_red_colormap(n)
if nargin < 1 || isempty(n)
    n = 256;
end
anchors = [ ...
    1.00 1.00 1.00; ...
    1.00 1.00 0.65; ...
    1.00 0.65 0.20; ...
    0.78 0.00 0.00];
anchor_pos = [0.00, 0.10, 0.55, 1.00];
xi = linspace(0, 1, n);
cmap = interp1(anchor_pos, anchors, xi, 'linear');
cmap = max(min(cmap, 1), 0);
end

function plot_emg_exclusion_diagnostics(save_dir, subject_id, win_name, scan_freqs, searchTopos, ...
    searchMeanPrSpectrum, eval_vec, gamma_vec, occ_evidence, emg_score, emg_class, unknown_high_risk, ...
    hard_eligible, hard_reject_flags, soft_warn_flags, rejection_flags, adaptive_thr, cfg_topo, topo_labels, ...
    peak_form_score, peak_form_mode, crosswin_id)
nComp = numel(eval_vec);
if nComp < 1
    return;
end

figA = figure('Position', [0 0 1512 982]);
scatter_size = 60 + 80 * max(eval_vec - min(eval_vec), 0) / max(max(eval_vec) - min(eval_vec), eps);
g_pct = 100 * (exp(gamma_vec) - 1);
scatter(occ_evidence, emg_score, scatter_size, g_pct, 'filled'); hold on;
xline(adaptive_thr.occ_class_thr, 'k--', 'LineWidth', 1.0);
yline(adaptive_thr.emg_class_thr, 'k--', 'LineWidth', 1.0);
xline(0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
yline(0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
idx_unknown = find(unknown_high_risk);
if ~isempty(idx_unknown)
    scatter(occ_evidence(idx_unknown), emg_score(idx_unknown), 130, 'rx', 'LineWidth', 2.0);
end
colormap(gca, emg_white_yellow_orange_red_colormap(256));
g_finite = g_pct(isfinite(g_pct));
if isempty(g_finite)
    g_upper = 10;
else
    g_upper = prctile(g_finite, 95);
    if ~isfinite(g_upper) || g_upper < 10
        g_upper = max(10, max(g_finite));
    end
end
if ~isfinite(g_upper) || g_upper <= 0
    g_upper = 10;
end
caxis([0 g_upper]);
x_vals = occ_evidence(isfinite(occ_evidence));
y_vals = emg_score(isfinite(emg_score));
if isempty(x_vals)
    x_absmax = 1;
else
    x_absmax = max(abs(x_vals));
end
if isempty(y_vals)
    y_absmax = 1;
else
    y_absmax = max(abs(y_vals));
end
x_absmax = max([x_absmax, abs(adaptive_thr.occ_class_thr), 0.1]);
y_absmax = max([y_absmax, abs(adaptive_thr.emg_class_thr), 0.1]);
x_lim_absmax = 1.10 * x_absmax;
y_lim_absmax = 1.10 * y_absmax;
xlim([-x_lim_absmax x_lim_absmax]);
ylim([-y_lim_absmax y_lim_absmax]);
% Light-green bracket highlighting the selected-component region:
% occipital evidence >= threshold and EMG score <= threshold.
bracket_color = [0.70 0.90 0.70];
plot([adaptive_thr.occ_class_thr, x_lim_absmax], [adaptive_thr.emg_class_thr, adaptive_thr.emg_class_thr], ...
    '-', 'Color', bracket_color, 'LineWidth', 2.2);
plot([adaptive_thr.occ_class_thr, adaptive_thr.occ_class_thr], [-y_lim_absmax, adaptive_thr.emg_class_thr], ...
    '-', 'Color', bracket_color, 'LineWidth', 2.2);
for ci = 1:nComp
    if isfinite(occ_evidence(ci)) && isfinite(emg_score(ci))
        cid = NaN;
        if nargin >= 20 && ~isempty(crosswin_id) && numel(crosswin_id) >= ci
            cid = crosswin_id(ci);
        end
        if isfinite(cid)
            lbl = sprintf('C%d[ID%d]', ci, round(cid));
        else
            lbl = sprintf('C%d', ci);
        end
        text(occ_evidence(ci) + 0.015 * x_lim_absmax, emg_score(ci) + 0.015 * y_lim_absmax, ...
            lbl, 'FontSize', 7, 'Color', [0.15 0.15 0.15], 'Interpreter', 'none');
    end
end
cb = colorbar;
cb.Label.String = 'Gamma change [%]';
xlabel('OccipitalEvidenceScore');
ylabel('EMGArtifactScore');
title(sprintf('EMG-vs-Occipital Separation: %s (%s)', subject_id, win_name), 'FontSize', 14, 'FontWeight', 'bold');
box on;
saveas(figA, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_scatter_%s.png', subject_id, win_name)));
close(figA);

sel_idx = find(hard_eligible);
rej_idx = find(hard_reject_flags);
[~, so] = sort(eval_vec(sel_idx), 'descend');
sel_idx = sel_idx(so);
[~, ro] = sort(eval_vec(rej_idx), 'descend');
rej_idx = rej_idx(ro);

% Top 10 selected/rejected components in one figure:
% row 1 = selected topographies, row 2 = selected spectra,
% row 3 = rejected topographies, row 4 = rejected spectra.
nCols = 10;
nShowSel = min(nCols, numel(sel_idx));
nShowRej = min(nCols, numel(rej_idx));
figSel = figure('Position', [0 0 1512 982]);
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
        if nargin >= 20 && ~isempty(crosswin_id) && numel(crosswin_id) >= ci && isfinite(crosswin_id(ci))
            ttl = sprintf('Sel C%d [ID%d] (%s)', ci, round(crosswin_id(ci)), emg_class{ci});
        else
            ttl = sprintf('Sel C%d (%s)', ci, emg_class{ci});
        end
        title(ttl, 'FontSize', 8, 'Interpreter', 'none');
    else
        axis off;
    end

    % Row 2: selected spectra
    subplot(4, nCols, nCols + k); hold on;
    if k <= nShowSel
        ci = sel_idx(k);
        plot(scan_freqs, searchMeanPrSpectrum(ci, :), 'k-', 'LineWidth', 1.3);
        yline(0, 'k--');
        xlabel('Hz'); ylabel('PR');
        pf_mode = 'none';
        if ~isempty(peak_form_mode) && numel(peak_form_mode) >= ci && ~isempty(peak_form_mode{ci})
            pf_mode = peak_form_mode{ci};
        end
        title(sprintf('\\lambda=%.2f, g=%.0f%%, PF=%.2f (%s)', ...
            eval_vec(ci), 100 * (exp(gamma_vec(ci)) - 1), peak_form_score(ci), pf_mode), 'FontSize', 8);
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
        if nargin >= 20 && ~isempty(crosswin_id) && numel(crosswin_id) >= ci && isfinite(crosswin_id(ci))
            ttl = sprintf('Rej C%d [ID%d] (%s)', ci, round(crosswin_id(ci)), emg_class{ci});
        else
            ttl = sprintf('Rej C%d (%s)', ci, emg_class{ci});
        end
        title(ttl, 'FontSize', 8, 'Interpreter', 'none');
    else
        axis off;
    end

    % Row 4: rejected spectra
    subplot(4, nCols, 3 * nCols + k); hold on;
    if k <= nShowRej
        ci = rej_idx(k);
        plot(scan_freqs, searchMeanPrSpectrum(ci, :), 'r-', 'LineWidth', 1.3);
        yline(0, 'k--');
        xlabel('Hz'); ylabel('PR');
        pf_mode = 'none';
        if ~isempty(peak_form_mode) && numel(peak_form_mode) >= ci && ~isempty(peak_form_mode{ci})
            pf_mode = peak_form_mode{ci};
        end
        title(sprintf('\\lambda=%.2f, g=%.0f%%, PF=%.2f (%s)', ...
            eval_vec(ci), 100 * (exp(gamma_vec(ci)) - 1), peak_form_score(ci), pf_mode), 'FontSize', 8);
        box on;
    else
        axis off;
    end
end
sgtitle(sprintf('Top 10 Selected and Rejected Components: %s (%s)', subject_id, win_name), ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(figSel, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_topo_spectra_selected_%s.png', subject_id, win_name)));
close(figSel);

figC = figure('Position', [0 0 756 982]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
counts = [sum(hard_eligible), sum(logical(soft_warn_flags)), sum(logical(hard_reject_flags)), numel(eval_vec)];
cats = {'Selected', 'SoftWarn', 'HardRejected', 'Candidates'};
bar(counts, 0.6, 'FaceColor', [0.3 0.4 0.7]); hold on;
set(gca, 'XTick', 1:numel(cats), 'XTickLabel', cats, 'XTickLabelRotation', 20);
ylabel('Count');
title('Component counts');
box on;
nexttile;
reason_names = {'unknown_proxy', 'flat_topomap', 'front_leak', 'temp_leak', 'corr', ...
    'lineharm', 'stationarity', 'burst', 'hf_slope', 'condlock', 'emg_score', 'occ_margin'};
reason_counts = zeros(1, numel(reason_names));
for ri = 1:numel(reason_names)
    rname = reason_names{ri};
    if isfield(rejection_flags, rname)
        reason_counts(ri) = sum(hard_reject_flags & rejection_flags.(rname));
    end
end
bar(reason_counts, 0.65, 'FaceColor', [0.72 0.40 0.30]); hold on;
set(gca, 'XTick', 1:numel(reason_names), 'XTickLabel', reason_names, 'XTickLabelRotation', 35);
ylabel('Rejected components');
title('Rejection reasons (criterion-level)');
box on;
sgtitle(sprintf('EMG Exclusion Summary: %s (%s)', subject_id, win_name), 'FontSize', 14, 'FontWeight', 'bold');
saveas(figC, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_summary_%s.png', subject_id, win_name)));
close(figC);
end

function [peak_bonus_vec, peak_count_vec] = compute_peak_bonus_from_spectra( ...
    mean_pr_spectrum, scan_freqs, analysis_freq_range, poly_order, detrend_edge_exclude_n, ...
    detrend_in_log, detrend_flat_edges, peak_min_prom_frac, peak_min_distance_hz)
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
    y_dt = detrend_power_ratio(y, scan_freqs, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
    y_band = y_dt(freq_mask);
    x_band = scan_freqs(freq_mask);
    valid = isfinite(y_band) & isfinite(x_band);
    y_band = y_band(valid);
    x_band = x_band(valid);
    if numel(y_band) < 5
        continue;
    end
    y_band = movmean(y_band, 3);
    peak_scale = max(y_band) - median(y_band);
    if ~isfinite(peak_scale) || peak_scale <= eps
        continue;
    end
    min_prom = max(0, max(y_band) * peak_min_prom_frac);
    [~, ~, ~, prom] = findpeaks(y_band, x_band, ...
        'MinPeakProminence', min_prom, ...
        'MinPeakDistance', peak_min_distance_hz);
    if isempty(prom)
        continue;
    end
    prom = sort(prom(:), 'descend');
    n_keep = min(2, numel(prom));
    peak_count_vec(ci) = n_keep;
    prom_score = min(1, mean(prom(1:n_keep)) / max(peak_scale, eps));
    count_score = n_keep / 2;
    peak_bonus_vec(ci) = max(0, min(1, 0.65 * prom_score + 0.35 * count_score));
end
end

function [peak_form_score_vec, peak_form_mode_vec, diag] = compute_peak_form_template_score_from_spectra( ...
    mean_pr_spectrum, scan_freqs, analysis_freq_range, poly_order, detrend_edge_exclude_n, ...
    detrend_in_log, detrend_flat_edges, shift_max_hz, single_widths_hz, double_widths_hz, ...
    double_separations_hz, min_trough_depth, min_similarity, smooth_n)
nComp = size(mean_pr_spectrum, 1);
peak_form_score_vec = zeros(nComp, 1);
peak_form_mode_vec = repmat({'none'}, nComp, 1);
diag = struct( ...
    'best_single_similarity', nan(nComp, 1), ...
    'best_double_similarity', nan(nComp, 1), ...
    'best_similarity_raw', nan(nComp, 1), ...
    'best_shift_hz', nan(nComp, 1), ...
    'best_width_hz', nan(nComp, 1), ...
    'best_separation_hz', nan(nComp, 1), ...
    'best_trough_depth', nan(nComp, 1));
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
    y_dt = detrend_power_ratio(y, scan_freqs, poly_order, detrend_edge_exclude_n, detrend_in_log, detrend_flat_edges);
    x_band = scan_freqs(freq_mask);
    y_band = y_dt(freq_mask);
    valid = isfinite(x_band) & isfinite(y_band);
    x_band = x_band(valid);
    y_band = y_band(valid);
    if numel(y_band) < 7
        continue;
    end
    y_band = movmean(y_band, smooth_n);
    y_norm = normalize_positive_shape(y_band);
    if isempty(y_norm)
        continue;
    end

    [single_best, single_meta] = evaluate_single_template_bank(y_norm, x_band, single_widths_hz, shift_max_hz);
    [double_best, double_meta] = evaluate_double_template_bank( ...
        y_norm, y_band, x_band, double_widths_hz, double_separations_hz, shift_max_hz, min_trough_depth);

    diag.best_single_similarity(ci) = single_best;
    diag.best_double_similarity(ci) = double_best;

    mode_raw = 'single';
    best_raw = single_best;
    best_shift = single_meta.shift_hz;
    best_width = single_meta.width_hz;
    best_sep = NaN;
    best_trough = NaN;
    best_centers = single_meta.centers_hz;

    if double_best > best_raw
        mode_raw = 'double';
        best_raw = double_best;
        best_shift = double_meta.shift_hz;
        best_width = double_meta.width_hz;
        best_sep = double_meta.sep_hz;
        best_trough = double_meta.trough_depth;
        best_centers = double_meta.centers_hz;
    end

    % Penalize template matches that are forced to band edges.
    edge_margin_hz = 2;
    if ~isempty(best_centers) && any(best_centers <= (analysis_freq_range(1) + edge_margin_hz) | ...
            best_centers >= (analysis_freq_range(2) - edge_margin_hz))
        best_raw = 0.85 * best_raw;
    end

    % Penalize EMG-like monotonic high-frequency rise in the upper gamma band.
    hf_mask = x_band >= max(70, analysis_freq_range(2) - 15);
    if sum(hf_mask) >= 5
        hf_idx = find(hf_mask);
        hf_rho = corr((1:numel(hf_idx))', y_band(hf_idx)', 'rows', 'complete', 'type', 'Spearman');
        if isfinite(hf_rho) && hf_rho > 0.70
            best_raw = best_raw * max(0.65, 1 - 0.30 * (hf_rho - 0.70) / 0.30);
        end
    end

    diag.best_similarity_raw(ci) = best_raw;
    diag.best_shift_hz(ci) = best_shift;
    diag.best_width_hz(ci) = best_width;
    diag.best_separation_hz(ci) = best_sep;
    diag.best_trough_depth(ci) = best_trough;

    if ~isfinite(best_raw) || best_raw <= min_similarity
        peak_form_score_vec(ci) = 0;
        peak_form_mode_vec{ci} = 'none';
    else
        peak_form_score_vec(ci) = max(0, min(1, (best_raw - min_similarity) / max(1 - min_similarity, eps)));
        peak_form_mode_vec{ci} = mode_raw;
    end
end
end

function [best_sim, meta] = evaluate_single_template_bank(y_norm, x_band, widths_hz, shift_max_hz)
best_sim = 0;
meta = struct('shift_hz', NaN, 'width_hz', NaN, 'centers_hz', []);
if isempty(widths_hz)
    return;
end
x_mid = mean(x_band);
for wi = 1:numel(widths_hz)
    w = widths_hz(wi);
    if ~isfinite(w) || w <= 0
        continue;
    end
    shift_vals = candidate_shift_values_hz(x_band, shift_max_hz);
    for si = 1:numel(shift_vals)
        shift_hz = shift_vals(si);
        center_hz = x_mid + shift_hz;
        tpl = gaussian_template(x_band, center_hz, w);
        sim = safe_template_similarity(y_norm, tpl);
        if sim > best_sim
            best_sim = sim;
            meta.shift_hz = shift_hz;
            meta.width_hz = w;
            meta.centers_hz = center_hz;
        end
    end
end
end

function [best_sim, meta] = evaluate_double_template_bank(y_norm, y_band, x_band, widths_hz, separations_hz, shift_max_hz, min_trough_depth)
best_sim = 0;
meta = struct('shift_hz', NaN, 'width_hz', NaN, 'sep_hz', NaN, 'trough_depth', NaN, 'centers_hz', []);
if isempty(widths_hz) || isempty(separations_hz)
    return;
end
x_mid = mean(x_band);
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
            tpl = gaussian_template(x_band, c1, w) + gaussian_template(x_band, c2, w);
            sim = safe_template_similarity(y_norm, tpl);
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
shifts = -shift_max_hz:df:shift_max_hz;
if isempty(shifts)
    shifts = 0;
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

function y_norm = normalize_positive_shape(y)
y_norm = [];
if isempty(y)
    return;
end
y = y(:);
y = y - median(y(isfinite(y)));
y(~isfinite(y)) = 0;
y = max(y, 0);
y_mag = norm(y);
if ~isfinite(y_mag) || y_mag <= eps
    return;
end
y_norm = y / y_mag;
end

function sim = safe_template_similarity(y_norm, tpl)
sim = 0;
if isempty(y_norm) || isempty(tpl)
    return;
end
tpl = tpl(:);
tpl = max(tpl, 0);
tpl_mag = norm(tpl);
if ~isfinite(tpl_mag) || tpl_mag <= eps
    return;
end
tpl = tpl / tpl_mag;
if numel(tpl) ~= numel(y_norm)
    return;
end
sim = y_norm(:)' * tpl(:);
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

function [ok, stats] = validate_peak_form_shift_invariance( ...
    scan_freqs, analysis_freq_range, shift_max_hz, single_widths_hz, double_widths_hz, ...
    double_separations_hz, min_trough_depth, min_similarity, smooth_n, centers_hz, width_hz, spread_tol)
stats = struct( ...
    'single_scores', nan(numel(centers_hz), 1), ...
    'double_scores', nan(numel(centers_hz), 1), ...
    'single_spread', inf, ...
    'double_spread', inf);
ok = true;
if isempty(scan_freqs) || isempty(centers_hz)
    return;
end
analysis_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
x_band = scan_freqs(analysis_mask);
if isempty(x_band)
    return;
end
single_specs = nan(numel(centers_hz), numel(scan_freqs));
double_specs = nan(numel(centers_hz), numel(scan_freqs));
for ci = 1:numel(centers_hz)
    c = centers_hz(ci);
    single_specs(ci, :) = gaussian_template(scan_freqs, c, width_hz);
    sep_hz = 12;
    c1 = c - sep_hz / 2;
    c2 = c + sep_hz / 2;
    double_specs(ci, :) = gaussian_template(scan_freqs, c1, width_hz) + gaussian_template(scan_freqs, c2, width_hz);
end
[single_scores, ~] = compute_peak_form_template_score_from_spectra( ...
    single_specs, scan_freqs, analysis_freq_range, 0, 0, false, false, shift_max_hz, ...
    single_widths_hz, double_widths_hz, double_separations_hz, min_trough_depth, min_similarity, smooth_n);
[double_scores, ~] = compute_peak_form_template_score_from_spectra( ...
    double_specs, scan_freqs, analysis_freq_range, 0, 0, false, false, shift_max_hz, ...
    single_widths_hz, double_widths_hz, double_separations_hz, min_trough_depth, min_similarity, smooth_n);
stats.single_scores = single_scores;
stats.double_scores = double_scores;
stats.single_spread = max(single_scores) - min(single_scores);
stats.double_spread = max(double_scores) - min(double_scores);
ok = isfinite(stats.single_spread) && isfinite(stats.double_spread) && ...
    (stats.single_spread <= spread_tol) && (stats.double_spread <= spread_tol);
end

function signature = build_component_signature_block(searchTopos, searchMeanPrSpectrum, corr_vec, ratio_vec, ...
    condlock_vec, lineharm_vec, hf_slope_vec, leak_vec, temp_leak_vec, eval_raw_vec)
nComp = size(searchTopos, 2);
nFreq = size(searchMeanPrSpectrum, 2);
signature = struct();
signature.nComp = nComp;
signature.topo_norm = zeros(size(searchTopos));
signature.spec_norm = zeros(size(searchMeanPrSpectrum));
signature.scalar_raw = nan(nComp, 8);
signature.scalar_z = nan(nComp, 8);
if nComp < 1
    return;
end
for ci = 1:nComp
    topo_ci = searchTopos(:, ci);
    topo_finite = topo_ci(isfinite(topo_ci));
    if isempty(topo_finite)
        topo_finite = 0;
    end
    topo_ci = topo_ci - mean(topo_finite);
    topo_scale = norm(topo_ci(isfinite(topo_ci)));
    if ~isfinite(topo_scale) || topo_scale <= eps
        topo_scale = 1;
    end
    topo_ci(~isfinite(topo_ci)) = 0;
    signature.topo_norm(:, ci) = topo_ci / topo_scale;

    if nFreq > 0
        spec_ci = searchMeanPrSpectrum(ci, :);
        spec_ci = spec_ci(:)';
        if any(isfinite(spec_ci))
            spec_ci = spec_ci - nanmean(spec_ci);
            spec_scale = norm(spec_ci(isfinite(spec_ci)));
            if ~isfinite(spec_scale) || spec_scale <= eps
                spec_scale = 1;
            end
            spec_ci(~isfinite(spec_ci)) = 0;
            signature.spec_norm(ci, :) = spec_ci / spec_scale;
        end
    end
end
signature.scalar_raw(:, 1) = corr_vec(:);
signature.scalar_raw(:, 2) = ratio_vec(:);
signature.scalar_raw(:, 3) = condlock_vec(:);
signature.scalar_raw(:, 4) = lineharm_vec(:);
signature.scalar_raw(:, 5) = hf_slope_vec(:);
signature.scalar_raw(:, 6) = leak_vec(:);
signature.scalar_raw(:, 7) = temp_leak_vec(:);
signature.scalar_raw(:, 8) = eval_raw_vec(:);
for fi = 1:size(signature.scalar_raw, 2)
    signature.scalar_z(:, fi) = normalize_robust(signature.scalar_raw(:, fi));
end
end

function match = match_components_full_anchor(sig_full, sig_early, sig_late, cfg)
nFull = sig_full.nComp;
nEarly = sig_early.nComp;
nLate = sig_late.nComp;
id_full = (1:nFull)';
id_early = nan(nEarly, 1);
id_late = nan(nLate, 1);
conf_full = nan(nFull, 1);
conf_early = nan(nEarly, 1);
conf_late = nan(nLate, 1);

[map_full_to_early, conf_full_early, map_early_to_full, conf_early_full] = ...
    match_anchor_to_target(sig_full, sig_early, cfg);
[map_full_to_late, conf_full_late, map_late_to_full, conf_late_full] = ...
    match_anchor_to_target(sig_full, sig_late, cfg);

for ei = 1:nEarly
    fi = map_early_to_full(ei);
    if fi > 0
        id_early(ei) = id_full(fi);
        conf_early(ei) = conf_early_full(ei);
    end
end
for li = 1:nLate
    fi = map_late_to_full(li);
    if fi > 0
        id_late(li) = id_full(fi);
        conf_late(li) = conf_late_full(li);
    end
end
conf_full = max([conf_full_early(:), conf_full_late(:)], [], 2);

next_id = nFull + 1;
for ei = 1:nEarly
    if ~isfinite(id_early(ei))
        id_early(ei) = next_id;
        conf_early(ei) = 0;
        next_id = next_id + 1;
    end
end
for li = 1:nLate
    if ~isfinite(id_late(li))
        id_late(li) = next_id;
        conf_late(li) = 0;
        next_id = next_id + 1;
    end
end
for fi = 1:nFull
    if ~isfinite(conf_full(fi))
        conf_full(fi) = 0;
    end
end

match = struct();
match.id_full = id_full;
match.id_early = id_early;
match.id_late = id_late;
match.conf_full = conf_full;
match.conf_early = conf_early;
match.conf_late = conf_late;
match.full_to_early = map_full_to_early;
match.full_to_late = map_full_to_late;
match.full_to_early_conf = conf_full_early;
match.full_to_late_conf = conf_full_late;
match.cfg = cfg;
end

function [anchor_to_target, conf_anchor, target_to_anchor, conf_target] = match_anchor_to_target(sig_anchor, sig_target, cfg)
nAnchor = sig_anchor.nComp;
nTarget = sig_target.nComp;
anchor_to_target = zeros(nAnchor, 1);
conf_anchor = zeros(nAnchor, 1);
target_to_anchor = zeros(nTarget, 1);
conf_target = zeros(nTarget, 1);
if nAnchor < 1 || nTarget < 1
    return;
end
sim_mat = build_crosswindow_similarity_matrix(sig_anchor, sig_target, cfg);
sim_work = sim_mat;
assigned_anchor = false(nAnchor, 1);
assigned_target = false(nTarget, 1);
while true
    sim_work(assigned_anchor, :) = -Inf;
    sim_work(:, assigned_target) = -Inf;
    [best_val, lin_idx] = max(sim_work(:));
    if ~isfinite(best_val)
        break;
    end
    [ai, tj] = ind2sub(size(sim_work), lin_idx);
    assigned_anchor(ai) = true;
    assigned_target(tj) = true;
    if best_val >= cfg.min_conf
        anchor_to_target(ai) = tj;
        conf_anchor(ai) = best_val;
        target_to_anchor(tj) = ai;
        conf_target(tj) = best_val;
    end
end
end

function sim_mat = build_crosswindow_similarity_matrix(sig_a, sig_b, cfg)
nA = sig_a.nComp;
nB = sig_b.nComp;
sim_mat = nan(nA, nB);
for ai = 1:nA
    topo_a = sig_a.topo_norm(:, ai);
    spec_a = sig_a.spec_norm(ai, :);
    feat_a = sig_a.scalar_z(ai, :);
    for bj = 1:nB
        topo_b = sig_b.topo_norm(:, bj);
        spec_b = sig_b.spec_norm(bj, :);
        feat_b = sig_b.scalar_z(bj, :);
        sim_topo = corr(topo_a, topo_b, 'rows', 'complete');
        if ~isfinite(sim_topo), sim_topo = 0; end
        sim_topo = abs(sim_topo);
        sim_spec = corr(spec_a(:), spec_b(:), 'rows', 'complete');
        if ~isfinite(sim_spec), sim_spec = 0; end
        sim_spec = max(min((sim_spec + 1) / 2, 1), 0);
        feat_diff = feat_a - feat_b;
        feat_diff = feat_diff(isfinite(feat_diff));
        if isempty(feat_diff)
            sim_feat = 0;
        else
            feat_dist = sqrt(mean(feat_diff .^ 2));
            sim_feat = max(0, 1 - min(feat_dist / 3, 1));
        end
        sim_mat(ai, bj) = cfg.w_topo * sim_topo + cfg.w_spec * sim_spec + cfg.w_feat * sim_feat;
    end
end
sim_mat = max(min(sim_mat, 1), 0);
end

function [group_consistent, group_occipital, group_member_count, group_mean_conf] = ...
    compute_crosswindow_group_consistency(match, tbl_full, tbl_early, tbl_late, min_conf)
all_ids = unique([match.id_full(:); match.id_early(:); match.id_late(:)]);
all_ids = all_ids(isfinite(all_ids));
group_consistent = containers.Map('KeyType', 'double', 'ValueType', 'any');
group_occipital = containers.Map('KeyType', 'double', 'ValueType', 'any');
group_member_count = containers.Map('KeyType', 'double', 'ValueType', 'any');
group_mean_conf = containers.Map('KeyType', 'double', 'ValueType', 'any');
for ii = 1:numel(all_ids)
    gid = all_ids(ii);
    idx_f = find(match.id_full == gid);
    idx_e = find(match.id_early == gid);
    idx_l = find(match.id_late == gid);
    member_count = numel(idx_f) + numel(idx_e) + numel(idx_l);
    occ_votes = 0;
    hard_votes = 0;
    conf_vals = [];
    if ~isempty(idx_f)
        occ_votes = occ_votes + sum(strcmp(tbl_full.emg_class(idx_f), 'occipital'));
        hard_votes = hard_votes + sum(tbl_full.hard_reject(idx_f));
        conf_vals = [conf_vals; match.conf_full(idx_f)]; %#ok<AGROW>
    end
    if ~isempty(idx_e)
        occ_votes = occ_votes + sum(strcmp(tbl_early.emg_class(idx_e), 'occipital'));
        hard_votes = hard_votes + sum(tbl_early.hard_reject(idx_e));
        conf_vals = [conf_vals; match.conf_early(idx_e)]; %#ok<AGROW>
    end
    if ~isempty(idx_l)
        occ_votes = occ_votes + sum(strcmp(tbl_late.emg_class(idx_l), 'occipital'));
        hard_votes = hard_votes + sum(tbl_late.hard_reject(idx_l));
        conf_vals = [conf_vals; match.conf_late(idx_l)]; %#ok<AGROW>
    end
    conf_vals = conf_vals(isfinite(conf_vals) & conf_vals > 0);
    if isempty(conf_vals)
        mean_conf = 0;
    else
        mean_conf = mean(conf_vals);
    end
    occipital_like = (occ_votes >= max(1, ceil(0.5 * member_count))) && (hard_votes == 0);
    consistent = (member_count >= 2) && (mean_conf >= min_conf) && occipital_like;
    group_consistent(gid) = logical(consistent);
    group_occipital(gid) = logical(occipital_like);
    group_member_count(gid) = member_count;
    group_mean_conf(gid) = mean_conf;
end
end

function [candidate_table, selected_idx, selected_weights] = apply_crosswindow_consistency_bonus( ...
    candidate_table, crosswin_id, group_consistent_map, group_occipital_map, bonus_val, max_components_to_combine)
nComp = numel(candidate_table.comp_idx);
if isempty(crosswin_id) || numel(crosswin_id) ~= nComp
    crosswin_id = nan(nComp, 1);
end
score_final = candidate_table.score_base;
bonus_vec = zeros(nComp, 1);
for ci = 1:nComp
    gid = crosswin_id(ci);
    if ~isfinite(gid)
        continue;
    end
    if isKey(group_consistent_map, gid) && isKey(group_occipital_map, gid) && ...
            logical(group_consistent_map(gid)) && logical(group_occipital_map(gid)) && ...
            candidate_table.hard_eligible(ci) && candidate_table.soft_warn(ci)
        bonus_vec(ci) = bonus_val;
    end
end
score_final = score_final + bonus_vec;
eligible = find(candidate_table.hard_eligible & isfinite(score_final));
if isempty(eligible)
    finite_idx = find(isfinite(candidate_table.score_base));
    if isempty(finite_idx)
        selected_idx = 1;
    else
        [~, ord_fallback] = sort(candidate_table.score_base(finite_idx), 'descend');
        selected_idx = finite_idx(ord_fallback(1));
    end
else
    [~, ord] = sort(score_final(eligible), 'descend');
    selected_idx = eligible(ord);
    selected_idx = selected_idx(1:min(max_components_to_combine, numel(selected_idx)));
end
selected_weights = candidate_table.eigenvalue(selected_idx)';
selected_weights(~isfinite(selected_weights) | selected_weights <= 0) = 0;
if sum(selected_weights) <= 0
    selected_weights = ones(1, numel(selected_idx));
end
selected_weights = selected_weights / sum(selected_weights);
candidate_table.consistency_bonus = bonus_vec;
candidate_table.score_final = score_final;
candidate_table.score = score_final;
candidate_table.crosswin_id = crosswin_id;
end

function values = lookup_group_property(id_vec, prop_map)
values = nan(numel(id_vec), 1);
for ii = 1:numel(id_vec)
    gid = id_vec(ii);
    if ~isfinite(gid)
        continue;
    end
    if isKey(prop_map, gid)
        val = prop_map(gid);
        if islogical(val)
            values(ii) = double(val);
        else
            values(ii) = val;
        end
    end
end
end

function flat_mask = detect_flat_topomap_outliers(topo_mat, mad_ratio_thr, iqr_ratio_thr, span_ratio_thr, std_ratio_thr)
nComp = size(topo_mat, 2);
flat_mask = false(nComp, 1);
if isempty(topo_mat)
    return;
end
if nargin < 4 || isempty(span_ratio_thr)
    span_ratio_thr = 0.35;
end
if nargin < 5 || isempty(std_ratio_thr)
    std_ratio_thr = 0.12;
end
for ci = 1:nComp
    topo_vec = topo_mat(:, ci);
    topo_vec = topo_vec(isfinite(topo_vec));
    if isempty(topo_vec)
        flat_mask(ci) = true;
        continue;
    end
    amp_max = max(abs(topo_vec));
    if ~isfinite(amp_max) || amp_max <= eps
        flat_mask(ci) = true;
        continue;
    end
    mad_ratio = mad(topo_vec, 1) / amp_max;
    iqr_ratio = iqr(topo_vec) / amp_max;
    span_ratio = (max(topo_vec) - min(topo_vec)) / amp_max;
    if ~isfinite(mad_ratio), mad_ratio = 0; end
    if ~isfinite(iqr_ratio), iqr_ratio = 0; end
    if ~isfinite(span_ratio), span_ratio = 0; end
    std_ratio = std(topo_vec) / amp_max;
    if ~isfinite(std_ratio), std_ratio = 0; end
    flat_mask(ci) = (mad_ratio <= mad_ratio_thr) || (iqr_ratio <= iqr_ratio_thr) || ...
        (span_ratio <= span_ratio_thr) || (std_ratio <= std_ratio_thr);
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
reason_order = {'unknown_proxy', 'flat_topomap', 'front_leak', 'temp_leak', 'corr', ...
    'lineharm', 'stationarity', 'burst', 'hf_slope', 'condlock', 'emg_score', 'occ_margin'};
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
    pr_scan(fi) = (p_stim - p_base) / max(p_base, eps);
end
end

function sim_results = run_ground_truth_simulation_validation(headmodel, n_reps, snr_levels, depth_levels, overlap_levels, artifact_levels, min_eig_thr, max_frontleak_thr, max_templeak_thr)
rng(917, 'twister');
nChans = numel(headmodel.layANThead.label);
chan_labels = headmodel.layANThead.label(:);
occ_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO)', 'once')), chan_labels);
front_mask = cellfun(@(l) ~isempty(regexp(l, '^(Fp|AF|F)', 'once')), chan_labels);
temp_mask = cellfun(@(l) ~isempty(regexp(l, '^(T|TP|FT)', 'once')), chan_labels);
if ~any(occ_mask), occ_mask(:) = false; occ_mask(1:max(1, round(0.15*nChans))) = true; end
if ~any(front_mask), front_mask(:) = false; front_mask(1:max(1, round(0.15*nChans))) = true; end
if ~any(temp_mask), temp_mask = ~occ_mask & ~front_mask; end

sim_rows = [];
for snr = snr_levels
    for depth = depth_levels
        for overlap = overlap_levels
            for art = artifact_levels
                for rep = 1:n_reps
                    topo_gamma = zeros(nChans, 1);
                    topo_gamma(occ_mask) = 1;
                    topo_gamma(front_mask) = -0.2;
                    topo_gamma = topo_gamma + 0.08 * randn(nChans, 1);
                    topo_gamma = topo_gamma / max(norm(topo_gamma), eps);

                    topo_non = randn(nChans, 1);
                    topo_non = topo_non / max(norm(topo_non), eps);
                    topo_emg = zeros(nChans, 1);
                    topo_emg(front_mask | temp_mask) = 1;
                    topo_emg = topo_emg + 0.2 * randn(nChans, 1);
                    topo_emg = topo_emg / max(norm(topo_emg), eps);

                    topo_non = (1-overlap) * topo_non + overlap * topo_gamma;
                    topo_non = topo_non / max(norm(topo_non), eps);
                    topo_gamma = depth * topo_gamma;

                    cov_noise = eye(nChans);
                    cov_base = cov_noise + (0.8 * topo_non * topo_non') + (art * topo_emg * topo_emg');
                    cov_stim = cov_noise + ((0.8 + snr) * topo_gamma * topo_gamma') + ...
                        (0.8 * topo_non * topo_non') + (1.2 * art * topo_emg * topo_emg');

                    [W, D] = eig(cov_stim, cov_base);
                    [evals, ord] = sort(real(diag(D)), 'descend');
                    W = W(:, ord);
                    n_search = min(20, numel(evals));
                    ev = evals(1:n_search);
                    occ_strength = abs(W(occ_mask, 1:n_search));
                    front_strength = abs(W(front_mask, 1:n_search));
                    temp_strength = abs(W(temp_mask, 1:n_search));
                    occ_mean = mean(occ_strength, 1)';
                    front_leak = mean(front_strength, 1)' ./ max(occ_mean, eps);
                    temp_leak = mean(temp_strength, 1)' ./ max(occ_mean, eps);
                    eligible = isfinite(ev) & (ev >= min_eig_thr) & ...
                        (front_leak <= max_frontleak_thr) & (temp_leak <= max_templeak_thr);
                    sel_idx = find(eligible);
                    if isempty(sel_idx)
                        sel_idx = 1;
                    end
                    sel_idx = sel_idx(:)';
                    [~, ord_sel] = sort(ev(sel_idx), 'descend');
                    sel_idx = sel_idx(ord_sel);
                    sel_idx = sel_idx(1:min(10, numel(sel_idx)));
                    sel_w = ev(sel_idx)';
                    sel_w(sel_w <= 0 | ~isfinite(sel_w)) = 0;
                    if sum(sel_w) <= 0
                        sel_w = ones(size(sel_w));
                    end
                    sel_w = sel_w / sum(sel_w);
                    w_comb = W(:, sel_idx) * sel_w(:);

                    gamma_corr = abs(corr(w_comb, topo_gamma, 'rows', 'complete'));
                    emg_corr = abs(corr(w_comb, topo_emg, 'rows', 'complete'));
                    is_tp = gamma_corr >= 0.50;
                    is_fp = emg_corr >= 0.50 & gamma_corr < 0.50;
                    loc_err = 1 - gamma_corr;
                    sim_rows = [sim_rows; snr, depth, overlap, art, rep, is_tp, is_fp, loc_err]; %#ok<AGROW>
                end
            end
        end
    end
end

if isempty(sim_rows)
    sim_results = struct('table', table(), 'summary', struct());
    return;
end
sim_tbl = array2table(sim_rows, 'VariableNames', ...
    {'snr', 'depth', 'overlap', 'artifact', 'rep', 'tp', 'fp', 'locerr'});
sim_results = struct();
sim_results.table = sim_tbl;
sim_results.summary = struct( ...
    'sensitivity', mean(sim_tbl.tp), ...
    'false_positive_rate', mean(sim_tbl.fp), ...
    'specificity', 1 - mean(sim_tbl.fp), ...
    'mean_localization_error', mean(sim_tbl.locerr));
fprintf('Simulation summary: sens=%.3f, spec=%.3f, fpr=%.3f, locErr=%.3f\n', ...
    sim_results.summary.sensitivity, sim_results.summary.specificity, ...
    sim_results.summary.false_positive_rate, sim_results.summary.mean_localization_error);
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
        fprintf(['     selection diagnostics: nSearch=%d, finite=%d, passEig=%d, artifactFlagged=%d, ', ...
                 'unknownHighRisk=%d, passRaw=%d, excludedOutlier=%d\n'], ...
            m.n_search, m.n_finite_metrics, m.n_pass_eig, m.n_artifact_flagged, ...
            m.n_unknown_high_risk, m.n_pass_all_raw, m.n_excluded_dominant_outlier);
        if isfield(m, 'top_fail_idx') && isfinite(m.top_fail_idx)
            fprintf(['     top failing component: C%d | eig=%.3f (thr=%.3f, pass=%d), ', ...
                     'corr=%.3f (thr=%.3f), ratio=%.3f, gamma=%.3f, artifact=%d\n'], ...
                m.top_fail_idx, m.top_fail_eig, m.thr_eig, m.top_fail_pass_eig, ...
                m.top_fail_corr, m.thr_corr, ...
                m.top_fail_ratio, m.top_fail_gamma, m.top_fail_artifact);
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
