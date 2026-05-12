%% GCP Gamma Peak Frequency and Power with Generalized Eigendecomposition (GED)
%
% Registered report Stage 1 publication script.
%
% Scope
%   - Full stimulus window only analysis (0 to 2 s) with baseline contrast
%      against the prestimulus interval (-1.5 to -0.25 s).
%   - Subject specific GED with regularized covariance decomposition and
%      weighted multi component selection from occipital evidence, spectral
%      profile quality, and artifact indicators.
%   - Trial resolved spectral power ratio scans in dB across 30 to 90 Hz,
%      followed by trial peak frequency and peak power extraction.
%   - Condition wise subject summaries and bootstrap based reliability
%      assessment via peak frequency confidence interval width.
%
% Main outputs
%   - `data/features/GCP_eeg_GED.mat` with full window features and GED
%      component selection statistics (`all_component_selection_stats`).
%   - Subject reliability table in controls folder:
%      `GCP_eeg_GED_subject_reliability.csv`.
%   - Diagnostic and summary figures in `figures/eeg/ged`.

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('GCP');
nSubj = length(subjects);
total_runtime_tic = tic;

%% Parameters
baseline_window = [-1.5, -0.25];
full_window = [0, 2.0];
analysis_freq_range = [30 90];
scan_freq_step_hz = 1; % analysis grid resolution (Hz)
scan_freqs = analysis_freq_range(1):scan_freq_step_hz:analysis_freq_range(2);
nFreqs = numel(scan_freqs);
scan_width = 2.0; % spectral smoothing (Hz) for mtmfft

lambda = 0.05;
ged_search_n = 10;
min_eigval = 1.1;
min_powspctrm_form = 0.8;
random_seed = 123;
trial_peak_smooth_n = 10; % moving-average smoothing
peak_power_halfwidth_hz = 5; % peak power = mean power within peak_hz +/-

peak_bootstrap_reps = 1000; % bootstrap repetitions
peak_bootstrap_ci_prct = [2.5 97.5]; % bootstrap confidence interval percentile
reliability_ci_width_median_max_hz = 8; % median CI width maximum (Hz)
reliability_ci_width_cond_max_hz = 10; % condition-wise CI width maximum (Hz)
reliability_min_pass_conditions = 3; % minimum number of conditions to pass reliability criteria

condNames = {'c25', 'c50', 'c75', 'c100'};
condLabels = {'25%', '50%', '75%', '100%'}; 
condCodes = [61 62 63 64]; 

gcp_root_path = paths.root;
gcp_feature_data_path = paths.features;
if ~exist(gcp_feature_data_path, 'dir')
    gcp_feature_data_path = gcp_root_path;
end
fig_save_dir_ged = fullfile(paths.figures, 'eeg', 'ged');
fig_save_dir_component_selection = fullfile(fig_save_dir_ged, 'component_selection');
if ~exist(fig_save_dir_ged, 'dir'), mkdir(fig_save_dir_ged); end
if ~exist(fig_save_dir_component_selection, 'dir'), mkdir(fig_save_dir_component_selection); end

%% Storage
trials_powratio = cell(4, nSubj);
trials_powratio_plotstat = cell(4, nSubj);
trials_powratio_fullscan = cell(4, nSubj);
trials_powratio_fullscan_plotstat = cell(4, nSubj);
trials_peaks = cell(4, nSubj);
trials_centroid = cell(4, nSubj);
trials_outlier_mask_freq = cell(4, nSubj);
trials_outlier_mask_power = cell(4, nSubj);
trials_mean = nan(4, nSubj);
trials_median = nan(4, nSubj);
trials_trialcv = nan(4, nSubj);
trials_gamma_power = nan(4, nSubj);
trials_gamma_power_plotstat = nan(4, nSubj);
trials_mean_centroid = nan(4, nSubj);
trials_median_centroid = nan(4, nSubj);

all_topos = cell(1, nSubj);
all_topo_labels = cell(1, nSubj);
all_eigenvalues = nan(1, nSubj);
all_selected_comp_idx = nan(1, nSubj);
all_selected_comp_corr = nan(1, nSubj);
all_selected_comp_eval = nan(1, nSubj);
all_selected_comp_indices_multi = cell(1, nSubj);
all_selected_comp_weights = cell(1, nSubj);
all_component_selection_stats = cell(1, nSubj);
all_condition_powspctrm = cell(4, nSubj);
all_condition_peak_freq = nan(4, nSubj);
all_condition_peak_power = nan(4, nSubj);
trials_powratio_components = cell(4, nSubj);
trial_counts_initial_by_subj = zeros(nSubj, 1);
trial_counts_retained_by_subj = zeros(nSubj, 1);

subject_peak_boot_ci_width = nan(4, nSubj);
subject_peak_boot_ci_low = nan(4, nSubj);
subject_peak_boot_ci_high = nan(4, nSubj);
subject_peak_boot_n_valid = zeros(4, nSubj);
subject_reliability_median_ci_width = nan(1, nSubj);
subject_reliability_n_pass_conditions = zeros(1, nSubj);
subject_reliability_n_valid_conditions = zeros(1, nSubj);
subject_reliability_condition_pass = false(4, nSubj);
subject_reliability_pass = false(1, nSubj);
subject_reliability_reason = repmat({''}, 1, nSubj);
subject_runtime_seconds = nan(nSubj, 1);

%% Subject loop
for subj = 1:nSubj
    subj_runtime_tic = tic;
    rng(random_seed + subj, 'twister');

    subject_id = subjects{subj};
    fprintf('[GED] Subject %s (%d/%d) full window only
', subject_id, subj, nSubj);
    datapath = fullfile(gcp_feature_data_path, subject_id, 'eeg');
    eeg_data = load(fullfile(datapath, 'dataEEG.mat'), ...
        'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    dataStructs = {eeg_data.dataEEG_c25, eeg_data.dataEEG_c50, eeg_data.dataEEG_c75, eeg_data.dataEEG_c100};
    fsample = dataStructs{1}.fsample;
    labels = dataStructs{1}.label;
    nChans = numel(labels);
    nSearch = min(ged_search_n, nChans);

    occ_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P10|P9)', 'once')), labels);
    occ_idx = find(occ_mask);
    front_mask = cellfun(@(l) ~isempty(regexp(l, '^(Fp|AF|F)', 'once')), labels);
    front_idx = find(front_mask);
    post_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P)', 'once')), labels);
    post_idx = find(post_mask);
    temp_mask = cellfun(@(l) ~isempty(regexp(l, '^(T|TP|FT)', 'once')), labels);
    temp_idx = find(temp_mask);

    post_w = zeros(nChans, 1);
    post_w(post_idx) = 1;
    if sum(post_w) > 0
        post_w = post_w / sum(post_w);
    else
        post_w = ones(nChans, 1) / nChans;
    end

    covStim_full = zeros(nChans);
    covBase_full = zeros(nChans);
    nTrials_total = 0;
    dat_per_cond = cell(1, 4);

    for cond = 1:4
        dat = dataStructs{cond};
        trlIdx = find(dat.trialinfo == condCodes(cond));
        if isempty(trlIdx), continue; end

        cfg = [];
        cfg.trials = trlIdx;
        dat = ft_selectdata(cfg, dat);
        dat_per_cond{cond} = dat;

        cfg_filt = [];
        cfg_filt.bpfilter = 'yes';
        cfg_filt.bpfreq = analysis_freq_range;
        cfg_filt.bpfilttype = 'fir';
        cfg_filt.bpfiltord = round(3 * fsample / analysis_freq_range(1));
        dat_gamma = ft_preprocessing(cfg_filt, dat);

        cfg_t = [];
        cfg_t.latency = baseline_window;
        dat_base = ft_selectdata(cfg_t, dat_gamma);
        cfg_t.latency = full_window;
        dat_stim_full = ft_selectdata(cfg_t, dat_gamma);

        nTrl = numel(dat_stim_full.trial);
        if nTrl > 0
            cfg_cov = [];
            cfg_cov.covariance = 'yes';
            cfg_cov.covariancewindow = 'all';
            cfg_cov.removemean = 'yes';
            tl_base = ft_timelockanalysis(cfg_cov, dat_base);
            tl_stim_full = ft_timelockanalysis(cfg_cov, dat_stim_full);
            covBase_full = covBase_full + double(tl_base.cov) * nTrl;
            covStim_full = covStim_full + double(tl_stim_full.cov) * nTrl;
        end
        nTrials_total = nTrials_total + nTrl;
    end

    if nTrials_total < 1
        warning('No valid trials for %s. Skipping subject.', subject_id);
        continue;
    end
    covStim_full = covStim_full / nTrials_total;
    covBase_full = covBase_full / nTrials_total;

    covStim_reg = (1-lambda)*covStim_full + lambda*mean(diag(covStim_full))*eye(nChans);
    covBase_reg = (1-lambda)*covBase_full + lambda*mean(diag(covBase_full))*eye(nChans);
    [W_full, D_full] = eig(covStim_reg, covBase_reg);
    [evals_sorted_full, sortIdx] = sort(real(diag(D_full)), 'descend');
    W_full = W_full(:, sortIdx);

    template_front_weight = 0.75;
    template_sigma_occ = 0.12;
    template_sigma_front = 0.25;
    sim_template = zeros(nChans, 1);
    lay_labels = headmodel.layANThead.label;
    lay_pos = headmodel.layANThead.pos;
    chan_pos = nan(nChans, 2);
    for ch = 1:nChans
        li = find(strcmp(lay_labels, labels{ch}), 1, 'first');
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

    searchFilters_full = nan(nChans, nSearch);
    searchTopos_full = nan(nChans, nSearch);
    searchCorrs_full = nan(nSearch, 1);
    searchOccStrength = nan(nSearch, 1);
    searchFrontStrength = nan(nSearch, 1);
    searchOccFrontRatio = nan(nSearch, 1);
    searchFrontLeak = nan(nSearch, 1);
    searchTempLeak = nan(nSearch, 1);
    searchLineHarmRatio = nan(nSearch, 1);
    searchHFSlope = nan(nSearch, 1);
    searchMeanPrSpectrum_full = nan(nSearch, numel(scan_freqs));
    searchTopoPosteriorConcentration = nan(nSearch, 1);
    searchEmgClass_full = repmat({'unassigned'}, nSearch, 1);

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
        if isempty(topo_finite), topo_finite = 0; end
        posterior_mass = mean(topo_abs(post_idx));
        total_mass = mean(topo_finite);
        topo_posterior_concentration_ci = posterior_mass / max(total_mass, eps);

        proxy_ci = estimate_component_artifact_proxies( ...
            w_ci, dat_per_cond, full_window, baseline_window, fsample, scan_freqs, scan_width);

        searchFilters_full(:, ci) = w_ci;
        searchTopos_full(:, ci) = topo_ci;
        searchCorrs_full(ci) = r_ci;
        searchOccStrength(ci) = occ_strength;
        searchFrontStrength(ci) = front_strength;
        searchOccFrontRatio(ci) = ratio_ci;
        searchFrontLeak(ci) = front_leak_ci;
        searchTempLeak(ci) = temp_leak_ci;
        searchLineHarmRatio(ci) = proxy_ci.lineharm_ratio;
        searchHFSlope(ci) = proxy_ci.hf_slope;
        searchMeanPrSpectrum_full(ci, :) = proxy_ci.mean_pr_spectrum(:)';
        searchTopoPosteriorConcentration(ci) = topo_posterior_concentration_ci;
    end

    max_combined_leak = 1.30;
    max_emg_score = 0.85;
    max_components_to_combine = 10;
    corr_vec = searchCorrs_full;
    ratio_vec = searchOccFrontRatio;
    eval_raw_vec = evals_sorted_full(1:nSearch);
    leak_vec = searchFrontLeak;
    temp_leak_vec = searchTempLeak;
    combined_leak_vec = 0.5 * (leak_vec + temp_leak_vec);
    lineharm_vec = searchLineHarmRatio;
    hf_slope_vec = searchHFSlope;
    hf_slope_for_score = hf_slope_vec;
    hf_slope_for_score(~isfinite(hf_slope_for_score)) = 0;
    finite_metrics = isfinite(corr_vec) & isfinite(ratio_vec) & isfinite(eval_raw_vec) & ...
        isfinite(leak_vec) & isfinite(temp_leak_vec) & isfinite(combined_leak_vec) & isfinite(lineharm_vec);

    peak_bonus_vec = compute_peak_bonus_from_spectra(searchMeanPrSpectrum_full, scan_freqs, analysis_freq_range);
    [powspctrm_form_score_full, powspctrm_form_deduct_full] = ...
        compute_powspctrm_form_laplacian_score_from_spectra(searchMeanPrSpectrum_full, scan_freqs, analysis_freq_range);

    occipital_evidence = 0.40 * normalize_robust(corr_vec) + ...
        0.25 * normalize_robust(ratio_vec) + ...
        0.35 * normalize_robust(searchTopoPosteriorConcentration);
    emg_artifact_score_full = 0.30 * normalize_robust(leak_vec) + ...
        0.20 * normalize_robust(temp_leak_vec) + ...
        0.30 * normalize_robust(lineharm_vec) + ...
        0.20 * normalize_robust(max(hf_slope_for_score, 0));

    pass_eig_gate = finite_metrics & (eval_raw_vec >= min_eigval);
    pass_peak_gate = finite_metrics & (powspctrm_form_score_full >= min_powspctrm_form);
    artifact_flags = finite_metrics & ((combined_leak_vec > max_combined_leak) | (emg_artifact_score_full > max_emg_score));
    occipital_class_mask = occipital_evidence > emg_artifact_score_full;
    eligible_full = pass_eig_gate & pass_peak_gate & ~artifact_flags & occipital_class_mask;

    searchScores = compute_calibrated_rank_aggregation_score( ...
        eval_raw_vec, powspctrm_form_score_full, peak_bonus_vec, occipital_evidence, emg_artifact_score_full, 200);
    searchScores(~finite_metrics) = -Inf;
    searchScores(~eligible_full) = -Inf;

    selected_idx_full = find(eligible_full & isfinite(searchScores));
    if isempty(selected_idx_full)
        selected_idx_full = 1;
    else
        [~, ord] = sort(searchScores(selected_idx_full), 'descend');
        selected_idx_full = selected_idx_full(ord);
        selected_idx_full = selected_idx_full(1:min(max_components_to_combine, numel(selected_idx_full)));
    end

    w_combined_full = evals_sorted_full(selected_idx_full)';
    w_combined_full(~isfinite(w_combined_full) | w_combined_full <= 0) = 0;
    if sum(w_combined_full) <= 0
        w_combined_full = ones(1, numel(selected_idx_full));
    end
    w_combined_full = w_combined_full / sum(w_combined_full);

    W_combined_full = searchFilters_full(:, selected_idx_full);
    W_combined_full = normalize_filters_to_noise_metric(W_combined_full, covBase_full);
    filter_vec_full = build_combined_filter_vector(W_combined_full, w_combined_full);

    topo_temp_full = searchTopos_full(:, selected_idx_full) * w_combined_full(:);
    all_topos{subj} = topo_temp_full;
    all_topo_labels{subj} = labels;
    all_selected_comp_indices_multi{subj} = selected_idx_full;
    all_selected_comp_weights{subj} = w_combined_full(:)';
    all_selected_comp_idx(subj) = selected_idx_full(1);
    all_selected_comp_corr(subj) = searchCorrs_full(selected_idx_full(1));
    all_selected_comp_eval(subj) = evals_sorted_full(selected_idx_full(1));
    all_eigenvalues(subj) = evals_sorted_full(selected_idx_full(1));

    all_component_selection_stats{subj} = struct( ...
        'selection_mode', 'full_window_weighted', ...
        'selected_idx', selected_idx_full, ...
        'selected_weights', w_combined_full, ...
        'best_idx', selected_idx_full(1), ...
        'best_score', searchScores(selected_idx_full(1)), ...
        'best_corr', searchCorrs_full(selected_idx_full(1)), ...
        'eligible', eligible_full, ...
        'powspctrm_form_score', powspctrm_form_score_full, ...
        'powspctrm_form_deduct', powspctrm_form_deduct_full, ...
        'rejection_flags', artifact_flags, ...
        'emg_artifact_score', emg_artifact_score_full, ...
        'occipital_evidence', occipital_evidence);
    all_component_selection_stats{subj} = all_component_selection_stats{subj};

    trial_counts_initial_local = 0;
    trial_counts_retained_local = 0;

    for cond = 1:4
        dat = dat_per_cond{cond};
        if isempty(dat), continue; end

        nTrl = numel(dat.trial);
        trial_counts_initial_local = trial_counts_initial_local + nTrl;

        trial_cache = cell(nTrl, 1);
        baseline_power_raw = nan(nTrl, 1);
        baseline_power_comb_full = nan(nTrl, 1);
        for trl = 1:nTrl
            x = double(dat.trial{trl});
            t = dat.time{trl};
            idx_base = t >= baseline_window(1) & t <= baseline_window(2);
            idx_full = t >= full_window(1) & t <= full_window(2);
            x_base = x(:, idx_base);
            x_full = x(:, idx_full);
            trial_cache{trl} = struct('x_base', x_base, 'x_full', x_full);
            if ~isempty(x_base)
                pow_base_chan = mean(x_base.^2, 2);
                baseline_power_raw(trl) = sum(post_w(:) .* pow_base_chan(:));
                if ~isempty(W_combined_full)
                    x_base_full = W_combined_full' * x_base;
                    baseline_power_comb_full(trl) = mean(x_base_full(:).^2);
                end
            end
        end

        bad_base = flag_unreliable_baseline_trials(baseline_power_raw, 3.5);
        bad_base = bad_base | flag_unreliable_baseline_trials(baseline_power_comb_full, 3.5);
        [base_floor_full, ~] = compute_baseline_floor_stats(baseline_power_comb_full, 20, 0.25);

        has_base = false(nTrl, 1);
        has_full = false(nTrl, 1);
        for trl = 1:nTrl
            tc = trial_cache{trl};
            has_base(trl) = ~isempty(tc.x_base);
            has_full(trl) = ~isempty(tc.x_full);
        end
        trial_mask_full = has_base & has_full & ~bad_base;

        [powratio_components, near_floor_count_full] = compute_scan_ratio_for_window_batch( ...
            trial_cache, W_combined_full, 'x_full', trial_mask_full, ...
            fsample, scan_freqs, scan_width, base_floor_full, 1.5);
        [ratio_trials_full_combined, near_floor_count_full_combined, valid_freq_counts_full] = ...
            compute_scan_ratio_for_combined_filter_batch( ...
            trial_cache, filter_vec_full, 'x_full', trial_mask_full, fsample, scan_freqs, scan_width, base_floor_full, 1.5);

        unstable_frac = zeros(nTrl, 1);
        for trl = 1:nTrl
            if valid_freq_counts_full(trl) > 0
                unstable_frac(trl) = near_floor_count_full_combined(trl) / max(valid_freq_counts_full(trl), 1);
            elseif trl <= numel(near_floor_count_full)
                unstable_frac(trl) = near_floor_count_full(trl) / max(numel(scan_freqs), 1);
            end
        end
        reject_instability = unstable_frac >= 0.35;

        powratio_trials_fullscan = ratio_trials_full_combined;
        powratio_trials_fullscan(reject_instability, :) = NaN;
        trial_counts_retained_local = trial_counts_retained_local + sum(any(isfinite(powratio_trials_fullscan), 2));

        analysis_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
        [trl_peaks, trl_peak_power, trl_centroid] = compute_trial_peak_metrics_from_powratio_fullscan( ...
            powratio_trials_fullscan, scan_freqs, analysis_mask, trial_peak_smooth_n, peak_power_halfwidth_hz);

        [outlier_mask_freq, ~] = detect_trial_metric_outliers_iqr(trl_peaks, 3.0);
        [outlier_mask_power, ~] = detect_trial_metric_outliers_iqr(trl_peak_power, 3.0);
        valid_for_plot = ~(outlier_mask_freq | outlier_mask_power);

        trl_peaks_plot = trl_peaks;
        trl_peaks_plot(~valid_for_plot) = NaN;
        trl_peak_power_plot = trl_peak_power;
        trl_peak_power_plot(~valid_for_plot) = NaN;
        powratio_plot = powratio_trials_fullscan;
        powratio_plot(~valid_for_plot, :) = NaN;

        trials_powratio{cond, subj} = powratio_trials_fullscan;
        trials_powratio_plotstat{cond, subj} = powratio_plot;
        trials_powratio_fullscan{cond, subj} = powratio_trials_fullscan;
        trials_powratio_fullscan_plotstat{cond, subj} = powratio_plot;
        trials_powratio_components{cond, subj} = powratio_components;
        trials_peaks{cond, subj} = trl_peaks_plot;
        trials_centroid{cond, subj} = trl_centroid;
        trials_outlier_mask_freq{cond, subj} = outlier_mask_freq;
        trials_outlier_mask_power{cond, subj} = outlier_mask_power;

        trials_mean(cond, subj) = robust_trial_mean(trl_peaks_plot);
        trials_median(cond, subj) = median(trl_peaks_plot(isfinite(trl_peaks_plot)), 'omitnan');
        if ~any(isfinite(trl_peaks_plot)), trials_median(cond, subj) = NaN; end
        trials_trialcv(cond, subj) = robust_mad(trl_peaks_plot) / max(abs(robust_trial_mean(trl_peaks_plot)), eps);
        trials_gamma_power(cond, subj) = robust_trial_mean(trl_peak_power);
        trials_gamma_power_plotstat(cond, subj) = robust_trial_mean(trl_peak_power_plot);
        trials_mean_centroid(cond, subj) = robust_trial_mean(trl_centroid);
        trials_median_centroid(cond, subj) = median(trl_centroid(isfinite(trl_centroid)), 'omitnan');
        if ~any(isfinite(trl_centroid)), trials_median_centroid(cond, subj) = NaN; end

        avg_curve = compute_condition_average_powratio_ft(powratio_plot, scan_freqs);
        all_condition_powspctrm{cond, subj} = avg_curve;
        [peak_hz, peak_pow] = pick_tallest_peak(avg_curve, scan_freqs, trial_peak_smooth_n, peak_power_halfwidth_hz);
        all_condition_peak_freq(cond, subj) = peak_hz;
        all_condition_peak_power(cond, subj) = peak_pow;

        [ci_w, ci_l, ci_h, n_valid_boot] = compute_bootstrap_peak_precision_from_trials( ...
            powratio_plot, scan_freqs, trial_peak_smooth_n, peak_bootstrap_reps, peak_bootstrap_ci_prct);
        subject_peak_boot_ci_width(cond, subj) = ci_w;
        subject_peak_boot_ci_low(cond, subj) = ci_l;
        subject_peak_boot_ci_high(cond, subj) = ci_h;
        subject_peak_boot_n_valid(cond, subj) = n_valid_boot;
    end

    [median_ci_width, n_pass_cond, n_valid_cond, cond_pass_col, subj_pass, reason_str] = ...
        evaluate_criterion4_reliability_gate( ...
        subject_peak_boot_ci_width(:, subj)', 'FULL', condLabels, ...
        reliability_ci_width_median_max_hz, reliability_ci_width_cond_max_hz, reliability_min_pass_conditions);

    subject_reliability_median_ci_width(subj) = median_ci_width;
    subject_reliability_n_pass_conditions(subj) = n_pass_cond;
    subject_reliability_n_valid_conditions(subj) = n_valid_cond;
    subject_reliability_condition_pass(:, subj) = cond_pass_col;
    subject_reliability_pass(subj) = subj_pass;
    subject_reliability_reason{subj} = reason_str;

    trial_counts_initial_by_subj(subj) = trial_counts_initial_local;
    trial_counts_retained_by_subj(subj) = trial_counts_retained_local;
    subject_runtime_seconds(subj) = toc(subj_runtime_tic);
end

%% Group figures
fig_condition_avg_powspctrm = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
for cond = 1:4
    subj_curves = nan(nSubj, numel(scan_freqs));
    for s = 1:nSubj
        curv = all_condition_powspctrm{cond, s};
        if ~isempty(curv) && numel(curv) == numel(scan_freqs)
            subj_curves(s, :) = curv(:)';
        end
    end
    med_curve = nanmedian(subj_curves, 1);
    mad_curve = nan(1, numel(scan_freqs));
    for fi = 1:numel(scan_freqs)
        mad_curve(fi) = robust_mad(subj_curves(:, fi));
    end
    good = isfinite(scan_freqs) & isfinite(med_curve) & isfinite(mad_curve);
    if sum(good) < 3, continue; end
    x = scan_freqs(good);
    y = med_curve(good);
    e = mad_curve(good);
    fill([x, fliplr(x)], [y - e, fliplr(y + e)], colors(cond, :), ...
        'FaceAlpha', 0.14, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(x, y, '-', 'Color', colors(cond, :), 'LineWidth', 3.0, 'DisplayName', condLabels{cond});
    peak_freq_group = nanmedian(all_condition_peak_freq(cond, :));
    if isfinite(peak_freq_group)
        xline(peak_freq_group, ':', 'Color', colors(cond, :), 'LineWidth', 1.8, 'HandleVisibility', 'off');
    end
end
yline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlim([analysis_freq_range(1), analysis_freq_range(2)]);
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
title('GED Condition Averaged Power Spectra, Full Window', 'FontWeight', 'bold');
legend(condLabels, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
set(gca, 'Box', 'on');
save_figure_png(fig_condition_avg_powspctrm, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_condition_average_powspctrm_full.png'));

fig_trial_retention = figure('Position', [0 0 1512 982], 'Color', 'w');
bar(1:nSubj, [trial_counts_initial_by_subj trial_counts_retained_by_subj], 'grouped');
xticks(1:nSubj);
xticklabels(subjects);
xtickangle(45);
xlabel('Participant');
ylabel('Trial count');
title('Trial Retention, Full Window', 'FontWeight', 'bold');
legend({'Initial', 'Retained'}, 'Location', 'best');
grid on;
set(gca, 'Box', 'on');
save_figure_png(fig_trial_retention, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_trial_retention_by_subject_full.png'));

%% Save results
save_path = fullfile(gcp_root_path, 'data', 'features', 'GCP_eeg_GED.mat');
save(save_path, ...
    'trials_powratio', 'trials_powratio_plotstat', ...
    'trials_powratio_fullscan', 'trials_powratio_fullscan_plotstat', ...
    'trials_peaks', 'trials_centroid', ...
    'trials_mean', 'trials_median', 'trials_trialcv', ...
    'trials_gamma_power', 'trials_gamma_power_plotstat', ...
    'trials_mean_centroid', 'trials_median_centroid', ...
    'trials_outlier_mask_freq', 'trials_outlier_mask_power', ...
    'all_topos', 'all_topo_labels', 'all_eigenvalues', ...
    'all_selected_comp_idx', 'all_selected_comp_corr', 'all_selected_comp_eval', ...
    'all_selected_comp_indices_multi', 'all_selected_comp_weights', ...
    'all_component_selection_stats', 'all_component_selection_stats', ...
    'trials_powratio_components', 'all_condition_powspctrm', ...
    'all_condition_peak_freq', 'all_condition_peak_power', ...
    'trial_counts_initial_by_subj', 'trial_counts_retained_by_subj', ...
    'subject_peak_boot_ci_width', 'subject_peak_boot_ci_low', 'subject_peak_boot_ci_high', ...
    'subject_peak_boot_n_valid', 'subject_reliability_median_ci_width', ...
    'subject_reliability_n_pass_conditions', 'subject_reliability_n_valid_conditions', ...
    'subject_reliability_condition_pass', 'subject_reliability_pass', 'subject_reliability_reason', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames', ...
    'peak_bootstrap_reps', 'peak_bootstrap_ci_prct', ...
    'reliability_ci_width_median_max_hz', 'reliability_ci_width_cond_max_hz', ...
    'reliability_min_pass_conditions');

csv_output_dir = paths.controls;
if ~exist(csv_output_dir, 'dir')
    mkdir(csv_output_dir);
end
T_rel = table( ...
    subjects(:), ...
    subject_reliability_pass(:), ...
    subject_reliability_median_ci_width(:), ...
    subject_reliability_n_pass_conditions(:), ...
    subject_reliability_n_valid_conditions(:), ...
    subject_peak_boot_ci_width(1, :)', ...
    subject_peak_boot_ci_width(2, :)', ...
    subject_peak_boot_ci_width(3, :)', ...
    subject_peak_boot_ci_width(4, :)', ...
    subject_reliability_reason(:), ...
    'VariableNames', {'subject', 'criterion4_pass', ...
    'criterion4_median_ci_width_hz', 'criterion4_n_pass_conditions', ...
    'criterion4_n_valid_conditions', 'ci_width_25_hz', 'ci_width_50_hz', ...
    'ci_width_75_hz', 'ci_width_100_hz', 'criterion4_reason'});
writetable(T_rel, fullfile(csv_output_dir, 'GCP_eeg_GED_subject_reliability.csv'));

fprintf('[GED] DONE full window pipeline
');
fprintf('[GED] Feature extraction results saved to: %s
', save_path);
fprintf('[GED] Reliability summary saved to: %s
', csv_output_dir);
fprintf('[GED] Criterion 4 reliability pass FULL: %d/%d
', sum(subject_reliability_pass), nSubj);
for si = 1:nSubj
    if isfinite(subject_runtime_seconds(si))
        fprintf('[GED] Runtime Subject %s: %s
', subjects{si}, format_runtime_hhmmss(subject_runtime_seconds(si)));
    else
        fprintf('[GED] Runtime Subject %s: n/a
', subjects{si});
    end
end
fprintf('[GED] Runtime TOTAL: %s
', format_runtime_hhmmss(toc(total_runtime_tic)));

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

function [ratio_db, near_floor_freq_mask, near_floor_row_mask] = compute_scan_ratio_from_timeseries(sig_stim, sig_base, fs, scan_freqs, scan_width_hz, base_floor, near_floor_mult)
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
[p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(sig_stim, sig_base, fs, scan_freqs, scan_width_hz);
if isempty(p_stim_scan) || isempty(p_base_scan)
    return;
end
% IMPORTANT: derive stabilization floor from spectral baseline power per row.
% Using a global time-domain floor can over-dominate spectral bins and
% collapse ratios toward 0 dB.
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

    valid = isfinite(p_stim_row) & isfinite(p_base_row) & (p_base_row > 0);
    ratio_row = nan(1, numel(scan_freqs));
    ratio_row(valid) = 10 * log10((p_stim_row(valid) + floor_row) ./ (p_base_row(valid) + floor_row));
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

function [ratio_cube, near_floor_freq_count_per_trial] = compute_scan_ratio_for_window_batch(trial_cache, search_filters, stim_field, trial_mask, fs, scan_freqs, scan_width_hz, base_floor, near_floor_mult)
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
    sig_stim_cells, sig_base_cells, fs, scan_freqs, scan_width_hz, base_floor, near_floor_mult);
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

function [score_vec, metrics, stability] = compute_calibrated_rank_aggregation_score( ...
    eval_raw_vec, powspctrm_form_vec, peak_bonus_vec, occipital_evidence_vec, emg_artifact_vec, n_boot)
nComp = numel(eval_raw_vec);
score_vec = -Inf(nComp, 1);
metrics = struct( ...
    'score_raw', nan(nComp, 1), ...
    'eig_score', nan(nComp, 1), ...
    'powspctrm_form_score', nan(nComp, 1), ...
    'peak_bonus_score', nan(nComp, 1), ...
    'occipital_score', nan(nComp, 1), ...
    'anti_emg_score', nan(nComp, 1));
stability = struct( ...
    'top1_freq', nan(nComp, 1), ...
    'rank_mean', nan(nComp, 1), ...
    'rank_std', nan(nComp, 1));
if isempty(n_boot) || ~isfinite(n_boot) || n_boot < 1
    n_boot = 0;
end
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
score_raw = mean(point_mat, 2, 'omitnan');
score_vec = score_raw;
score_vec(~finite_rows) = -Inf;
metrics.score_raw = score_raw;
metrics.eig_score = point_mat(:, 1);
metrics.powspctrm_form_score = point_mat(:, 2);
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

function proxy = estimate_component_artifact_proxies(filter_w, dat_per_cond, stim_window, base_window, fs, scan_freqs, scan_width)
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
trial_pr = compute_simple_power_ratio_scan_batch(stim_sig_cells, base_sig_cells, fs, scan_freqs, scan_width);
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

function pr_mat = compute_simple_power_ratio_scan_batch(sig_stim_cells, sig_base_cells, fs, scan_freqs, scan_width)
nRows = numel(sig_stim_cells);
pr_mat = nan(nRows, numel(scan_freqs));
if nRows == 0 || numel(sig_base_cells) ~= nRows || fs <= 0
    return;
end
[p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(sig_stim_cells, sig_base_cells, fs, scan_freqs, scan_width);
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

function p_scan = compute_scan_power_mtmfft_ft(sig, fs, scan_freqs, scan_width_hz)
if ~iscell(sig) && isvector(sig)
    sig = sig(:)';
end
if iscell(sig)
    nSig = numel(sig);
else
    nSig = size(sig, 1);
end
p_scan = nan(nSig, numel(scan_freqs));
if isempty(sig) || fs <= 0 || isempty(scan_freqs) || ~isfinite(scan_width_hz) || scan_width_hz <= 0
    return;
end
if ~iscell(sig) && size(sig, 2) < 8
    return;
end

valid_rows = false(nSig, 1);
dat = [];
dat.label = {'comp'};
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
cfg.tapsmofrq = scan_width_hz;
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

function [p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(sig_stim, sig_base, fs, scan_freqs, scan_width_hz)
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
p_all = compute_scan_power_mtmfft_ft(sig_all, fs, scan_freqs, scan_width_hz);
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
    compute_scan_ratio_for_combined_filter_batch(trial_cache, filter_vec, stim_field, trial_mask, fs, scan_freqs, scan_width_hz, base_floor, near_floor_mult)
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
    sig_stim_cells, sig_base_cells, fs, scan_freqs, scan_width_hz, base_floor, near_floor_mult);
for ri = 1:size(ratio_rows, 1)
    trl = row_trial_idx(ri);
    ratio_trials(trl, :) = ratio_rows(ri, :);
    valid_freq_count_per_trial(trl) = sum(isfinite(ratio_rows(ri, :)));
    near_floor_freq_count_per_trial(trl) = sum(near_floor_row_mask(ri, :));
end
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
exportgraphics(fig_handle, out_path, 'Resolution', 300);
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
freq_dat.label = {'comp'};
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

function [median_ci_width, n_pass_cond, n_valid_cond, cond_pass_col, subj_pass, reason_str] = ...
    evaluate_criterion4_reliability_gate(ci_width_row, window_tag, cond_label_cell, ...
    median_max_hz, cond_max_hz, min_pass_cond)
median_ci_width = NaN;
n_pass_cond = 0;
n_valid_cond = 0;
cond_pass_col = false(numel(cond_label_cell), 1);
subj_pass = false;
reason_str = '';
if isempty(window_tag)
    window_tag = '?';
elseif ~ischar(window_tag)
    window_tag = char(window_tag);
end
if isempty(ci_width_row)
    reason_str = sprintf('%s: FAIL: no CI-width data.', window_tag);
    return;
end
ci_width_row = ci_width_row(:)';
n_cond = numel(cond_label_cell);
if numel(ci_width_row) < n_cond
    ci_width_row = [ci_width_row nan(1, n_cond - numel(ci_width_row))];
elseif numel(ci_width_row) > n_cond
    ci_width_row = ci_width_row(1:n_cond);
end
cond_pass = isfinite(ci_width_row) & (ci_width_row <= cond_max_hz);
cond_pass_col = cond_pass(:);
n_pass_cond = sum(cond_pass);
n_valid_cond = sum(isfinite(ci_width_row));
valid_ci = ci_width_row(isfinite(ci_width_row));
if ~isempty(valid_ci)
    median_ci_width = median(valid_ci);
end
subj_pass = isfinite(median_ci_width) && ...
    (median_ci_width <= median_max_hz) && ...
    (n_pass_cond >= min_pass_cond);
if subj_pass
    reason_str = sprintf('%s: PASS: median CI width %.2f Hz; reliable conditions %d/%d.', ...
        window_tag, median_ci_width, n_pass_cond, n_cond);
elseif n_valid_cond == 0
    reason_str = sprintf('%s: FAIL: no conditions with retained trials for bootstrap reliability.', ...
        window_tag);
elseif ~isfinite(median_ci_width) || median_ci_width > median_max_hz
    reason_str = sprintf(['%s: FAIL: median CI width %.2f Hz exceeds threshold %.2f Hz ', ...
        '(%d reliable conditions).'], ...
        window_tag, median_ci_width, median_max_hz, n_pass_cond);
else
    reason_str = sprintf(['%s: FAIL: only %d/%d conditions passed CI-width threshold %.2f Hz ', ...
        '(requires >=%d).'], ...
        window_tag, n_pass_cond, n_cond, cond_max_hz, min_pass_cond);
end
end

function [ci_width_hz, ci_low_hz, ci_high_hz, n_boot_valid] = ...
    compute_bootstrap_peak_precision_from_trials(pr_mat, scan_freqs, smooth_n, n_boot, ci_prct)
ci_width_hz = NaN;
ci_low_hz = NaN;
ci_high_hz = NaN;
n_boot_valid = 0;
if isempty(ci_prct) || numel(ci_prct) ~= 2
    ci_prct = [2.5 97.5];
end
if isempty(n_boot) || ~isfinite(n_boot) || n_boot < 1
    n_boot = 1000;
end
if ~isfinite(smooth_n) || smooth_n < 1
    smooth_n = 1;
end
if isempty(pr_mat) || isempty(scan_freqs) || size(pr_mat, 2) ~= numel(scan_freqs)
    return;
end
valid_trials = any(isfinite(pr_mat), 2);
trial_rows = find(valid_trials);
n_valid_trials = numel(trial_rows);
if n_valid_trials == 0
    return;
end
peak_hz_boot = nan(n_boot, 1);
for bi = 1:n_boot
    sample_idx = trial_rows(randi(n_valid_trials, n_valid_trials, 1));
    boot_mat = pr_mat(sample_idx, :);
    % Uniform trial mean only (matches compute_condition_average_powratio_ft fallback);
    % avoids ft_selectdata here — otherwise ~n_boot FieldTrip calls per condition/window.
    boot_use = any(isfinite(boot_mat), 2);
    if ~any(boot_use)
        continue;
    end
    boot_avg = mean(boot_mat(boot_use, :), 1, 'omitnan');
    boot_avg = movmean(boot_avg, max(1, round(smooth_n)), 'omitnan');
    [peak_hz_boot(bi), ~] = pick_tallest_peak(boot_avg, scan_freqs, 1, 0);
end
peak_hz_boot = peak_hz_boot(isfinite(peak_hz_boot));
n_boot_valid = numel(peak_hz_boot);
if n_boot_valid < max(20, ceil(0.20 * n_boot))
    return;
end
ci_bounds = prctile(peak_hz_boot, ci_prct);
if numel(ci_bounds) ~= 2 || any(~isfinite(ci_bounds))
    return;
end
ci_low_hz = ci_bounds(1);
ci_high_hz = ci_bounds(2);
ci_width_hz = ci_high_hz - ci_low_hz;
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
