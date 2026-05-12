%% GCP Gamma Peak Frequency and Power with Generalized Eigendecomposition (GED)
%
%   - Stimulus window only analysis (0 to 2 s) with baseline contrast
%      against the prestimulus interval (-1.5 to -0.25 s).
%   - Subject specific GED with regularized covariance decomposition and
%      weighted multi component selection from occipital evidence, spectral
%      profile quality, and artifact indicators.
%   - Trial resolved spectral power ratio scans in dB across 30 to 90 Hz,
%      followed by trial peak frequency and peak power extraction.
%   - Condition wise subject summaries and diagnostic figures.
%
% Main outputs
%   - `data/features/GCP_eeg_GED.mat` with stimulus window features and GED
%      component selection statistics (`all_component_selection_stats`).
%   - Diagnostic and summary figures in `figures/eeg/ged`.

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('GCP');
nSubj = length(subjects);
total_runtime_tic = tic;

%% Parameters

baseline_window = [-1.5, -0.25];
stim_window = [0, 2.0];

% Gamma analysis
analysis_freq_range = [30 90];
scan_freq_step_hz = 1; % analysis grid resolution (Hz)
scan_freqs = analysis_freq_range(1):scan_freq_step_hz:analysis_freq_range(2);
nFreqs = length(scan_freqs);
scan_width = 2.0; % spectral smoothing (Hz) for mtmfft

% GED
lambda = 0.05;              % regularization
ged_search_n = 10;          % search first N GED components
min_eigval = 1.1;           % minimum GED eigenvalue (lambda >= 1.1)
min_powspctrm_form = 0.8;   % minimum PF score for candidate eligibility
occ_class_thr = 0.60;
emg_class_thr = 0.50;
min_occ_margin = 0.05;
random_seed = 123;
trial_peak_smooth_n = 10;     % moving-average smoothing
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
fig_save_dir_component_selection_root = fullfile(fig_save_dir_ged, 'component_selection');
if ~exist(fig_save_dir_ged, 'dir'), mkdir(fig_save_dir_ged); end
if ~exist(fig_save_dir_component_selection_root, 'dir'), mkdir(fig_save_dir_component_selection_root); end

%% Preallocate storage
trials_powratio = cell(4, nSubj);
trials_powratio_plotstat = cell(4, nSubj);
trials_powratio_scan = cell(4, nSubj);
trials_powratio_scan_plotstat = cell(4, nSubj);
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
subject_runtime_seconds = nan(nSubj, 1);

trials_powratio_components = cell(4, nSubj);
all_condition_powspctrm = cell(4, nSubj);
all_condition_peak_freq = nan(4, nSubj);
all_condition_peak_power = nan(4, nSubj);
trial_counts_initial_by_subj = zeros(nSubj, 1);
trial_counts_retained_by_subj = zeros(nSubj, 1);

%% Subject loop
global GED_FREQ_PROGRESS;
GED_FREQ_PROGRESS = [];

for subj = 1:nSubj
    subj_runtime_tic = tic;
    rng(random_seed + subj, 'twister');

    subject_id = subjects{subj};
    fig_save_dir_component_selection = fullfile(fig_save_dir_component_selection_root, subject_id);
    if ~exist(fig_save_dir_component_selection, 'dir'), mkdir(fig_save_dir_component_selection); end

    clc;
    fprintf('[GED] Subject %s (%d/%d)\n', subject_id, subj, nSubj);

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

    covStim = zeros(nChans);
    covBase = zeros(nChans);
    nTrials_total = 0;
    dat_per_cond = cell(1, 4);
    valid_cond_mask = false(1, 4);

    for cond = 1:4
        dat = dataStructs{cond};
        trlIdx = find(dat.trialinfo == 60 + cond);
        if isempty(trlIdx), continue; end

        cfg = [];
        cfg.trials = trlIdx;
        dat = ft_selectdata(cfg, dat);
        dat_per_cond{cond} = dat;
        valid_cond_mask(cond) = ~isempty(dat.trial);

        cfg_filt = [];
        cfg_filt.bpfilter = 'yes';
        cfg_filt.bpfreq = analysis_freq_range;
        cfg_filt.bpfilttype = 'fir';
        cfg_filt.bpfiltord = round(3 * fsample / analysis_freq_range(1));
        dat_gamma = ft_preprocessing(cfg_filt, dat);

        cfg_t = [];
        cfg_t.latency = baseline_window;
        dat_base = ft_selectdata(cfg_t, dat_gamma);
        cfg_t.latency = stim_window;
        dat_stim = ft_selectdata(cfg_t, dat_gamma);

        nTrl = numel(dat_stim.trial);
        if nTrl > 0
            cfg_cov = [];
            cfg_cov.covariance = 'yes';
            cfg_cov.covariancewindow = 'all';
            cfg_cov.removemean = 'yes';
            tl_base = ft_timelockanalysis(cfg_cov, dat_base);
            tl_stim = ft_timelockanalysis(cfg_cov, dat_stim);
            covBase = covBase + double(tl_base.cov) * nTrl;
            covStim = covStim + double(tl_stim.cov) * nTrl;
            nTrials_total = nTrials_total + nTrl;
        end
    end

    if nTrials_total < 1
        warning('No valid trials for %s. Skipping subject.', subject_id);
        continue;
    end

    covStim = covStim / nTrials_total;
    covBase = covBase / nTrials_total;

    n_freq_calls_total = nSearch + 2 * sum(valid_cond_mask);
    ged_freq_progress_reset(subject_id, subj, nSubj, n_freq_calls_total);

    plot_covariance_matrix_diagnostics( ...
        fig_save_dir_component_selection, subject_id, labels, covBase, {covStim}, {'Stimulus'}, lambda);

    covStim_reg = (1 - lambda) * covStim + lambda * mean(diag(covStim)) * eye(nChans);
    covBase_reg = (1 - lambda) * covBase + lambda * mean(diag(covBase)) * eye(nChans);
    [W_ged, D_ged] = eig(covStim_reg, covBase_reg);
    [evals_sorted, sortIdx] = sort(real(diag(D_ged)), 'descend');
    W_ged = W_ged(:, sortIdx);

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

    ged_freq_progress_set_phase('ProxySpectra Batch');
    searchFilters = nan(nChans, nSearch);
    searchTopos = nan(nChans, nSearch);
    searchCorrs = nan(nSearch, 1);
    searchOccFrontRatio = nan(nSearch, 1);
    searchFrontLeak = nan(nSearch, 1);
    searchTempLeak = nan(nSearch, 1);
    searchLineHarmRatio = nan(nSearch, 1);
    searchHFSlope = nan(nSearch, 1);
    searchMeanPrSpectrum = nan(nSearch, numel(scan_freqs));
    searchTopoPosteriorConcentration = nan(nSearch, 1);
    searchEmgClass = repmat({'unassigned'}, nSearch, 1);

    for ci = 1:nSearch
        w_ci = W_ged(:, ci);
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
            w_ci, dat_per_cond, stim_window, baseline_window, fsample, scan_freqs, scan_width);

        searchFilters(:, ci) = w_ci;
        searchTopos(:, ci) = topo_ci;
        searchCorrs(ci) = r_ci;
        searchOccFrontRatio(ci) = ratio_ci;
        searchFrontLeak(ci) = front_leak_ci;
        searchTempLeak(ci) = temp_leak_ci;
        searchLineHarmRatio(ci) = proxy_ci.lineharm_ratio;
        searchHFSlope(ci) = proxy_ci.hf_slope;
        searchMeanPrSpectrum(ci, :) = proxy_ci.mean_pr_spectrum(:)';
        searchTopoPosteriorConcentration(ci) = topo_posterior_concentration_ci;
    end

    max_combined_leak = 1.30;
    max_lineharm_ratio = 0.60;
    max_hf_slope = -0.15;
    max_emg_score = 0.85;
    max_components_to_combine = 10;
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
    finite_metrics = isfinite(corr_vec) & isfinite(ratio_vec) & isfinite(eval_raw_vec) & ...
        isfinite(leak_vec) & isfinite(temp_leak_vec) & isfinite(combined_leak_vec) & isfinite(lineharm_vec);

    peak_bonus_vec = compute_peak_bonus_from_spectra(searchMeanPrSpectrum, scan_freqs, analysis_freq_range);
    [powspctrm_form_score, powspctrm_form_deduct] = ...
        compute_powspctrm_form_laplacian_score_from_spectra(searchMeanPrSpectrum, scan_freqs, analysis_freq_range);

    occipital_evidence = 0.40 * normalize_robust(corr_vec) + ...
        0.25 * normalize_robust(ratio_vec) + ...
        0.35 * normalize_robust(searchTopoPosteriorConcentration);
    emg_artifact_score = 0.30 * normalize_robust(leak_vec) + ...
        0.20 * normalize_robust(temp_leak_vec) + ...
        0.30 * normalize_robust(lineharm_vec) + ...
        0.20 * normalize_robust(max(hf_slope_for_score, 0));

    pass_eig_gate = finite_metrics & (eval_raw_vec >= min_eigval);
    pass_peak_gate = finite_metrics & (powspctrm_form_score >= min_powspctrm_form);
    artifact_flags = finite_metrics & ((combined_leak_vec > max_combined_leak) | (emg_artifact_score > max_emg_score));
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
            if abs(occ_minus_emg_ci) < min_occ_margin
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

    occipital_class_mask = cellfun(@(c) strcmpi(c, 'occipital'), searchEmgClass(:));
    mixed_class_mask = cellfun(@(c) strcmpi(c, 'mixed'), searchEmgClass(:));
    mixed_rescue_mask = mixed_class_mask & ...
        (occipital_evidence >= 0.50) & ...
        (occ_minus_emg_vec >= -0.05) & ...
        (searchTopoPosteriorConcentration >= 0.95);
    eligible = pass_eig_gate & pass_peak_gate & ~artifact_flags & ...
        (occipital_class_mask | mixed_rescue_mask);

    searchScores = compute_calibrated_rank_aggregation_score( ...
        eval_raw_vec, powspctrm_form_score, peak_bonus_vec, occipital_evidence, emg_artifact_score);
    searchScores(~finite_metrics) = -Inf;
    searchScores(~eligible) = -Inf;

    selected_idx = find(eligible & isfinite(searchScores));
    if isempty(selected_idx)
        selected_idx = 1;
    else
        [~, ord] = sort(searchScores(selected_idx), 'descend');
        selected_idx = selected_idx(ord);
        selected_idx = selected_idx(1:min(max_components_to_combine, numel(selected_idx)));
    end

    w_combined = evals_sorted(selected_idx)';
    w_combined(~isfinite(w_combined) | w_combined <= 0) = 0;
    if sum(w_combined) <= 0
        w_combined = ones(1, numel(selected_idx));
    end
    w_combined = w_combined / sum(w_combined);

    W_combined = searchFilters(:, selected_idx);
    W_combined = normalize_filters_to_noise_metric(W_combined, covBase);
    filter_vec = build_combined_filter_vector(W_combined, w_combined);
    topo_temp = searchTopos(:, selected_idx) * w_combined(:);

    all_topos{subj} = topo_temp;
    all_topo_labels{subj} = labels;
    all_selected_comp_indices_multi{subj} = selected_idx;
    all_selected_comp_weights{subj} = w_combined(:)';
    all_selected_comp_idx(subj) = selected_idx(1);
    all_selected_comp_corr(subj) = searchCorrs(selected_idx(1));
    all_selected_comp_eval(subj) = evals_sorted(selected_idx(1));
    all_eigenvalues(subj) = evals_sorted(selected_idx(1));

    all_component_selection_stats{subj} = struct( ...
        'selection_mode', 'stim_window_weighted', ...
        'selected_idx', selected_idx, ...
        'selected_weights', w_combined, ...
        'best_idx', selected_idx(1), ...
        'best_score', searchScores(selected_idx(1)), ...
        'best_corr', searchCorrs(selected_idx(1)), ...
        'eligible', eligible, ...
        'powspctrm_form_score', powspctrm_form_score, ...
        'powspctrm_form_deduct', powspctrm_form_deduct, ...
        'rejection_flags', artifact_flags, ...
        'emg_artifact_score', emg_artifact_score, ...
        'occipital_evidence', occipital_evidence, ...
        'occ_minus_emg', occ_minus_emg_vec, ...
        'emg_class', {searchEmgClass}, ...
        'mixed_rescue_mask', mixed_rescue_mask);

    occ_highlight = labels(occ_idx);
    cfg_topo = [];
    cfg_topo.layout = headmodel.layANThead;
    cfg_topo.comment = 'no';
    cfg_topo.marker = 'off';
    cfg_topo.style = 'straight';
    cfg_topo.gridscale = 300;
    cfg_topo.zlim = 'maxabs';
    cfg_topo.colormap = '*RdBu';
    cfg_topo.figure = 'gcf';
    cfg_topo.highlight = {'on'};
    cfg_topo.highlightchannel = {occ_highlight};
    cfg_topo.highlightsymbol = {'.'};
    cfg_topo.highlightsize = {12};
    cfg_topo.highlightcolor = {[0 0 0]};

    plot_emg_exclusion_diagnostics( ...
        fig_save_dir_component_selection, subject_id, scan_freqs, searchTopos, ...
        searchMeanPrSpectrum, evals_sorted(1:nSearch), searchEmgClass, ...
        eligible, artifact_flags, ...
        searchFrontLeak, searchTempLeak, combined_leak_vec, ...
        lineharm_vec, hf_slope_vec, emg_artifact_score, ...
        false(nSearch, 1), ...
        cfg_topo, labels, powspctrm_form_score, powspctrm_form_deduct, ...
        max_combined_leak, max_hf_slope, max_emg_score, max_lineharm_ratio, ...
        selected_idx);

    trial_counts_initial_local = 0;
    trial_counts_retained_local = 0;
    subj_powratio_scan = cell(1, 4);
    subj_peaks = cell(1, 4);
    subj_centroid = cell(1, 4);
    subj_condition_avg = cell(1, 4);
    subj_condition_peak = nan(1, 4);

    cmap_div = interp1([0 0.5 1], ...
        [0.17 0.27 0.53; 0.97 0.97 0.97; 0.70 0.09 0.17], linspace(0, 1, 256));

    for cond = 1:4
        dat = dat_per_cond{cond};
        if isempty(dat), continue; end

        nTrl = numel(dat.trial);
        trial_counts_initial_local = trial_counts_initial_local + nTrl;

        trial_cache = cell(nTrl, 1);
        baseline_power_raw = nan(nTrl, 1);
        baseline_power_comb = nan(nTrl, 1);
        for trl = 1:nTrl
            x = double(dat.trial{trl});
            t = dat.time{trl};
            idx_base = t >= baseline_window(1) & t <= baseline_window(2);
            idx_stim = t >= stim_window(1) & t <= stim_window(2);
            x_base = x(:, idx_base);
            x_stim = x(:, idx_stim);
            trial_cache{trl} = struct('x_base', x_base, 'x_stim', x_stim);
            if ~isempty(x_base)
                pow_base_chan = mean(x_base.^2, 2);
                baseline_power_raw(trl) = sum(post_w(:) .* pow_base_chan(:));
                if ~isempty(W_combined)
                    x_base_comb = W_combined' * x_base;
                    baseline_power_comb(trl) = mean(x_base_comb(:).^2);
                end
            end
        end

        bad_base = flag_unreliable_baseline_trials(baseline_power_raw, 3.5);
        bad_base = bad_base | flag_unreliable_baseline_trials(baseline_power_comb, 3.5);
        [base_floor, ~] = compute_baseline_floor_stats(baseline_power_comb, 20, 0.25);

        has_base = false(nTrl, 1);
        has_stim = false(nTrl, 1);
        for trl = 1:nTrl
            tc = trial_cache{trl};
            has_base(trl) = ~isempty(tc.x_base);
            has_stim(trl) = ~isempty(tc.x_stim);
        end
        trial_mask = has_base & has_stim & ~bad_base;

        ged_freq_progress_set_phase(sprintf('TrialScan Batch %s', condLabels{cond}));
        [powratio_components, near_floor_count] = compute_scan_ratio_for_window_batch( ...
            trial_cache, W_combined, 'x_stim', trial_mask, ...
            fsample, scan_freqs, scan_width, base_floor, 1.5);
        [ratio_trials_combined, near_floor_count_combined, valid_freq_counts] = ...
            compute_scan_ratio_for_combined_filter_batch( ...
            trial_cache, filter_vec, 'x_stim', trial_mask, fsample, scan_freqs, scan_width, base_floor, 1.5);

        unstable_frac = zeros(nTrl, 1);
        for trl = 1:nTrl
            if valid_freq_counts(trl) > 0
                unstable_frac(trl) = near_floor_count_combined(trl) / max(valid_freq_counts(trl), 1);
            elseif trl <= numel(near_floor_count)
                unstable_frac(trl) = near_floor_count(trl) / max(numel(scan_freqs), 1);
            end
        end
        reject_instability = unstable_frac >= 0.35;

        powratio_trials_scan = ratio_trials_combined;
        powratio_trials_scan(reject_instability, :) = NaN;
        trial_counts_retained_local = trial_counts_retained_local + sum(any(isfinite(powratio_trials_scan), 2));

        analysis_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
        [trl_peaks, trl_peak_power, trl_centroid] = compute_trial_peak_metrics_from_powratio_scan( ...
            powratio_trials_scan, scan_freqs, analysis_mask, trial_peak_smooth_n, peak_power_halfwidth_hz);

        [outlier_mask_freq, ~] = detect_trial_metric_outliers_iqr(trl_peaks, 3.0);
        [outlier_mask_power, ~] = detect_trial_metric_outliers_iqr(trl_peak_power, 3.0);
        valid_for_plot = ~(outlier_mask_freq | outlier_mask_power);

        trl_peaks_plot = trl_peaks;
        trl_peaks_plot(~valid_for_plot) = NaN;
        trl_peak_power_plot = trl_peak_power;
        trl_peak_power_plot(~valid_for_plot) = NaN;
        powratio_plot = powratio_trials_scan;
        powratio_plot(~valid_for_plot, :) = NaN;

        trials_powratio{cond, subj} = powratio_trials_scan;
        trials_powratio_plotstat{cond, subj} = powratio_plot;
        trials_powratio_scan{cond, subj} = powratio_trials_scan;
        trials_powratio_scan_plotstat{cond, subj} = powratio_plot;
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

        subj_powratio_scan{cond} = powratio_trials_scan;
        subj_peaks{cond} = trl_peaks_plot;
        subj_centroid{cond} = trl_centroid;
        subj_condition_avg{cond} = avg_curve;
        subj_condition_peak(cond) = peak_hz;
    end

    fig_trials = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('Trial-Level GED: Subject %s', subject_id), ...
        'FontSize', 18, 'FontWeight', 'bold');

    pr_raw_mats = cell(1, 4);
    row1_clim = ones(1, 4);
    for cond = 1:4
        pr_raw_mats{cond} = subj_powratio_scan{cond};
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

    for cond = 1:4
        subplot(2, 4, cond);
        if ~isempty(pr_raw_mats{cond})
            hold on;
            imagesc(scan_freqs, 1:size(pr_raw_mats{cond}, 1), pr_raw_mats{cond});
            colormap(gca, cmap_div);
            caxis([-row1_clim(cond) row1_clim(cond)]);
            cb = colorbar;
            cb.FontSize = 8;
            ctd = subj_centroid{cond};
            if ~isempty(ctd)
                valid_ctd = ~isnan(ctd);
                tr_idx = find(valid_ctd);
                if ~isempty(tr_idx)
                    scatter(ctd(valid_ctd), tr_idx, 28, 'k', 'filled', ...
                        'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.75, ...
                        'MarkerFaceAlpha', 0.85);
                end
            end
            xlabel('Freq [Hz]');
            ylabel('Trial');
            set(gca, 'YDir', 'normal');
        end
        title(sprintf('%s Raw', condLabels{cond}), 'FontSize', 11);
        set(gca, 'FontSize', 10);
        xlim([30 90]);
        ylim([0 200]);
        box on;
    end

    subplot(2, 4, 5);
    if ~isempty(searchTopos) && ~isempty(selected_idx)
        topo_data = [];
        topo_data.label = labels;
        topo_plot = searchTopos(:, selected_idx) * w_combined(:);
        topo_data.avg = topo_plot;
        topo_data.dimord = 'chan';
        topo_abs = abs(topo_plot(post_idx));
        topo_abs = topo_abs(isfinite(topo_abs));
        if isempty(topo_abs)
            topo_abs = abs(topo_plot(isfinite(topo_plot)));
        end
        if isempty(topo_abs)
            topo_clim = 1;
        else
            topo_clim = prctile(topo_abs, 99.9);
            if ~isfinite(topo_clim) || topo_clim <= 0
                topo_clim = max(topo_abs);
            end
            if ~isfinite(topo_clim) || topo_clim <= 0
                topo_clim = 1;
            end
        end
        cfg_topo_plot = cfg_topo;
        cfg_topo_plot.zlim = [-topo_clim topo_clim];
        try
            ft_topoplotER(cfg_topo_plot, topo_data);
            cb = colorbar;
            cb.FontSize = 9;
        catch
            imagesc(topo_data.avg);
            caxis([-topo_clim topo_clim]);
            colorbar;
        end
        lambda_show = evals_sorted(selected_idx(1));
        title(sprintf('Weighted GED (%d comps, \\lambda=%.2f)', numel(selected_idx), lambda_show), 'FontSize', 11);
    end

    subplot(2, 4, [6 7]); hold on;
    edges = 30:2:90;
    hist_mat = zeros(4, numel(edges) - 1);
    for cond = 1:4
        tpk = subj_peaks{cond};
        if ~isempty(tpk)
            tpk = tpk(~isnan(tpk));
            hist_mat(cond, :) = histcounts(tpk, edges);
        end
    end
    centers = edges(1:end-1) + diff(edges) / 2;
    bh = bar(centers, hist_mat', 'stacked', 'EdgeColor', 'none', 'BarWidth', 1);
    for cond = 1:4
        bh(cond).FaceColor = colors(cond, :);
    end
    xlabel('Peak Frequency [Hz]');
    ylabel('Trial Count');
    title('Trial-Level Peak Distribution', 'FontSize', 12);
    legend(bh, condLabels, 'FontSize', 10, 'Location', 'best');
    set(gca, 'FontSize', 11);
    xlim([30 90]);
    box on;

    subplot(2, 4, 8); hold on;
    for cond = 1:4
        avg_curve = subj_condition_avg{cond};
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
        peak_hz = subj_condition_peak(cond);
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
        peak_hz = subj_condition_peak(cond);
        if ~isfinite(peak_hz)
            continue;
        end
        text(0.98, peak_text_y - (cond - 1) * peak_text_step, sprintf('%.0f Hz', peak_hz), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'Color', colors(cond, :), 'FontSize', 10, 'FontWeight', 'bold');
    end
    legend('Location', 'southwest', 'FontSize', 9, 'Box', 'off');
    save_figure_png(fig_trials, fullfile(fig_save_dir_component_selection, ...
        sprintf('GCP_eeg_GED_subj%s_trials_overview.png', subject_id)));
    close(fig_trials);

    trial_counts_initial_by_subj(subj) = trial_counts_initial_local;
    trial_counts_retained_by_subj(subj) = trial_counts_retained_local;
    subject_runtime_seconds(subj) = toc(subj_runtime_tic);
end

%% Condition-averaged spectra
fig_condition_avg_powspctrm = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
panel_min = inf;
panel_max = -inf;
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
    if sum(good) < 3
        continue;
    end
    x = scan_freqs(good);
    y = med_curve(good);
    e = mad_curve(good);
    fill([x, fliplr(x)], [y - e, fliplr(y + e)], colors(cond, :), ...
        'FaceAlpha', 0.14, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(x, y, '-', 'Color', colors(cond, :), 'LineWidth', 3.2, 'DisplayName', condLabels{cond});
    peak_freq_group = nanmedian(all_condition_peak_freq(cond, :));
    if isfinite(peak_freq_group)
        xline(peak_freq_group, ':', 'Color', colors(cond, :), 'LineWidth', 2.0, 'HandleVisibility', 'off');
    end
    panel_min = min(panel_min, min(y - e));
    panel_max = max(panel_max, max(y + e));
end
yline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlim([analysis_freq_range(1), analysis_freq_range(2)]);
if isfinite(panel_min) && isfinite(panel_max) && panel_max > panel_min
    yr = panel_max - panel_min;
    ylim([panel_min - 0.10 * yr, panel_max + 0.12 * yr]);
end
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
title('GED Condition-Averaged Power Spectra', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'FontSize', 11, 'Box', 'on');
legend(condLabels, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
save_figure_png(fig_condition_avg_powspctrm, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_condition_average_powspctrm.png'));

%% Trial retention figure
fig_trial_retention = figure('Position', [0 0 1512 982], 'Color', 'w');
bar(1:nSubj, [trial_counts_initial_by_subj trial_counts_retained_by_subj], 'grouped');
ax = gca;
if numel(ax.Children) >= 2
    ax.Children(2).FaceColor = [0.85 0.85 0.85];
    ax.Children(2).EdgeColor = [0.45 0.45 0.45];
    ax.Children(1).FaceColor = [0.20 0.55 0.85];
    ax.Children(1).EdgeColor = [0.10 0.30 0.55];
end
xlabel('Participant');
ylabel('Trial count');
title('Trial Retention by Participant', 'FontSize', 18, 'FontWeight', 'bold');
xticks(1:nSubj);
xticklabels(subjects);
xtickangle(45);
box on;
grid on;
legend({'Initial', 'Retained'}, 'Location', 'best');
save_figure_png(fig_trial_retention, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_trial_retention_by_subject.png'));

%% Centroid summary
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
        scatter(xj, vals, 110, colors(c, :), 'filled', ...
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
            y_all_c = [y_all_c; tc(:)]; %#ok<AGROW>
            g_all_c = [g_all_c; c * ones(length(tc), 1)]; %#ok<AGROW>
        end
    end
end
if ~isempty(y_all_c)
    boxplot(y_all_c, g_all_c, 'Colors', 'k', 'Symbol', '', 'Widths', 0.15);
    for c = 1:4
        vals = y_all_c(g_all_c == c);
        xj = c + 0.15 + (rand(size(vals)) - 0.5) * 0.22;
        scatter(xj, vals, 10, colors(c, :), 'filled', 'MarkerFaceAlpha', 0.2);
    end
end
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 13, 'Box', 'off');
xlim([0.3 4.7]);
ylim(analysis_freq_range);
ylabel('Centroid Frequency [Hz]');
title('All trials pooled', 'FontSize', 14, 'FontWeight', 'bold');
save_figure_png(fig_cent, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_centroid_summary.png'));

%% Condition separation summary
fig_cond_slope = figure('Position', [0 0 1512 982], 'Color', 'w');
slope_post = compute_condition_separation_from_matrix(all_condition_peak_freq);
delta_post = all_condition_peak_freq(4, :) - all_condition_peak_freq(1, :);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
valid_slope = isfinite(slope_post);
slope_vals = slope_post(valid_slope);
if ~isempty(slope_vals)
    boxplot(slope_vals(:), ones(numel(slope_vals), 1), 'Colors', 'k', 'Symbol', '', 'Widths', 0.2);
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

nexttile; hold on;
valid_delta = isfinite(delta_post);
delta_vals = delta_post(valid_delta);
if ~isempty(delta_vals)
    boxplot(delta_vals(:), ones(numel(delta_vals), 1), 'Colors', 'k', 'Symbol', '', 'Widths', 0.2);
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
if exist('subjects', 'var') == 1 && ~isempty(subj_idx) && numel(subjects) >= max(subj_idx) && numel(subj_idx) <= 35
    xticklabels(subjects(subj_idx));
end
xlabel('Subject', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('\Delta Frequency Shift (100% - 25%) [Hz]', 'FontSize', 18, 'FontWeight', 'bold');
title('Gamma Frequency Shift', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
save_figure_png(fig_cond_shift_bar_freq, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_bar_GammaFreq.png'));

%% Condition shift figure: frequency
fig_condition_shift_freq = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
dat_freq_shift = all_condition_peak_freq - all_condition_peak_freq(1, :);
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
title('Gamma Peak Frequency Shift', 'FontSize', 24, 'FontWeight', 'bold');
save_figure_png(fig_condition_shift_freq, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_condition_shift_frequency.png'));

%% Condition shift figure: power
fig_condition_shift_power = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
dat_power_shift = all_condition_peak_power - all_condition_peak_power(1, :);
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
title('Gamma Peak Power Shift', 'FontSize', 24, 'FontWeight', 'bold');
save_figure_png(fig_condition_shift_power, fullfile(fig_save_dir_ged, 'GCP_eeg_GED_condition_shift_power.png'));

%% Save results
save_path = fullfile(gcp_root_path, 'data', 'features', 'GCP_eeg_GED.mat');
save(save_path, ...
    'trials_powratio', 'trials_powratio_plotstat', ...
    'trials_powratio_scan', 'trials_powratio_scan_plotstat', ...
    'trials_peaks', 'trials_centroid', ...
    'trials_mean', 'trials_median', 'trials_trialcv', ...
    'trials_gamma_power', 'trials_gamma_power_plotstat', ...
    'trials_mean_centroid', 'trials_median_centroid', ...
    'trials_outlier_mask_freq', 'trials_outlier_mask_power', ...
    'all_topos', 'all_topo_labels', 'all_eigenvalues', ...
    'all_selected_comp_idx', 'all_selected_comp_corr', 'all_selected_comp_eval', ...
    'all_selected_comp_indices_multi', 'all_selected_comp_weights', ...
    'all_component_selection_stats', ...
    'trials_powratio_components', 'all_condition_powspctrm', ...
    'all_condition_peak_freq', 'all_condition_peak_power', ...
    'trial_counts_initial_by_subj', 'trial_counts_retained_by_subj', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames');

clc;
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
    compute_trial_peak_metrics_from_powratio_scan(powratio_trials_scan, scan_freqs, analysis_mask, ...
    smooth_n, peak_power_halfwidth_hz)
nTrl = size(powratio_trials_scan, 1);
trl_peaks = nan(nTrl, 1);
trl_peak_power = nan(nTrl, 1);
trl_centroid = nan(nTrl, 1);
scan_freqs_analysis = scan_freqs(analysis_mask);
centroid_band_mask = scan_freqs >= 30 & scan_freqs <= 90;
freq_band = scan_freqs(centroid_band_mask);
if ~isfinite(peak_power_halfwidth_hz) || peak_power_halfwidth_hz < 0
    peak_power_halfwidth_hz = 0;
end
for trl = 1:nTrl
    pr_scan = powratio_trials_scan(trl, :);
    if all(~isfinite(pr_scan))
        continue;
    end
    pr_proc = pr_scan(analysis_mask);
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

    pr_proc_scan = movmean(pr_scan, max(1, round(smooth_n)), 'omitnan');
    pr_proc_band = pr_proc_scan(centroid_band_mask);
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
score_raw = mean(point_mat, 2, 'omitnan');
score_vec = score_raw;
score_vec(~finite_rows) = -Inf;
end

function plot_emg_exclusion_diagnostics(save_dir, subject_id, scan_freqs, searchTopos, ...
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
sgtitle(sprintf('Top 5 Selected and Rejected Components: %s', subject_id), ...
    'FontSize', 14, 'FontWeight', 'bold');
drawnow;
pause(0.05);
save_figure_png(figSel, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_topo_spectra_selected.png', subject_id)));
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

function plot_covariance_matrix_diagnostics(save_dir, subject_id, chan_labels, covBase, covStim_per_win, win_names_cap, lambdas)
if isempty(covStim_per_win) || isempty(covBase)
    return;
end
nWins = min(3, numel(covStim_per_win));
if nWins < 1
    return;
end
if isempty(win_names_cap)
    win_names_cap = {'Stimulus'};
end
if isempty(lambdas)
    lambdas = repmat(0.05, 1, nWins);
end

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
nMatCols = 4;
nPlotCols = 5;
tiledlayout(nWins, nPlotCols, 'Padding', 'compact', 'TileSpacing', 'compact');

nChans = size(covBase, 1);
if size(covBase, 2) ~= nChans
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
    d_base = diag(covBase);
    d_base = d_base(isfinite(d_base));
    base_diag_mean = mean(d_base);
    if ~isfinite(base_diag_mean)
        base_diag_mean = 1;
    end

    covStim_reg = (1 - lam_w) * covStim_w + lam_w * stim_diag_mean * eye(nChans);
    covBase_reg = (1 - lam_w) * covBase + lam_w * base_diag_mean * eye(nChans);
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
    if nWins == 1
        panel_titles = {'Stim (reg)', 'Base (reg)', 'Stim-Base', 'R^{-1}S'};
    end

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
    if nWins == 1
        title('Values', 'FontSize', 10, 'Interpreter', 'none');
    else
        title(sprintf('%s Values', win_names_cap{wi}), 'FontSize', 10, 'Interpreter', 'none');
    end
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

function reason = build_rejection_reason_short(ci, emg_class, rejection_flags, extreme_component_outlier)
reasons = {};

if numel(emg_class) >= ci && ~isempty(emg_class{ci})
    cls = lower(strtrim(emg_class{ci}));
    if strcmp(cls, 'mixed') || strcmp(cls, 'emg')
        reasons{end + 1} = cls; %#ok<AGROW>
    end
end

if numel(extreme_component_outlier) >= ci && logical(extreme_component_outlier(ci))
    reasons{end + 1} = 'extreme_outlier'; %#ok<AGROW>
end
if get_flag_value(rejection_flags, 'combined_leak', ci)
    reasons{end + 1} = 'combined_leak'; %#ok<AGROW>
end
if get_flag_value(rejection_flags, 'emg_score', ci)
    reasons{end + 1} = 'emg_score'; %#ok<AGROW>
end
if get_flag_value(rejection_flags, 'topo_nonposterior', ci)
    reasons{end + 1} = 'topo_nonposterior'; %#ok<AGROW>
end

if isempty(reasons)
    reason = 'criteria_not_met';
    return;
end

n_show = min(2, numel(reasons));
reason = strjoin(reasons(1:n_show), ',');
if numel(reasons) > n_show
    reason = sprintf('%s,+%d', reason, numel(reasons) - n_show);
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

function [powspctrm_form_score_vec, powspctrm_form_mode_vec, diag] = compute_powspctrm_form_template_score_from_spectra( ...
    mean_pr_spectrum, scan_freqs, analysis_freq_range)
multi_peak_sep_min_hz = 15;
multi_peak_height_ratio_min = 0.80;
multi_peak_penalty_mult = 0.62;
multi_peak_penalty_mult = max(0, min(1, multi_peak_penalty_mult));
shift_max_hz = 10;
single_widths_hz = [5 7 9 12];
double_separations_hz = [8 12 16 20];
double_widths_hz = [3 4 5 6];
min_trough_depth = 0.10;
min_similarity = 0.40;
smooth_n = 5;
peak_width_min_hz = 2.0;
peak_width_max_hz = 20.0;
edge_ratio_soft = 1.25;
edge_ratio_limit = 1.75;
edge_run_soft = 0.025;
edge_run_limit = 0.060;
powspctrm_form_prom_abs_floor = 0.02;
nComp = size(mean_pr_spectrum, 1);
powspctrm_form_score_vec = zeros(nComp, 1);
powspctrm_form_mode_vec = repmat({'none'}, nComp, 1);
diag = struct( ...
    'best_single_similarity', nan(nComp, 1), ...
    'best_double_similarity', nan(nComp, 1), ...
    'best_similarity_raw', nan(nComp, 1), ...
    'pre_penalty_similarity', nan(nComp, 1), ...
    'total_penalty_raw', nan(nComp, 1), ...
    'total_penalty_used', nan(nComp, 1), ...
    'edge_penalty', nan(nComp, 1), ...
    'dominance_penalty', nan(nComp, 1), ...
    'dominance_rival_count_major', nan(nComp, 1), ...
    'dominance_rival_count_minor', nan(nComp, 1), ...
    'dominance_rival_amp_ratio_sum', nan(nComp, 1), ...
    'hf_rise_penalty', nan(nComp, 1), ...
    'best_shift_hz', nan(nComp, 1), ...
    'best_center_hz', nan(nComp, 1), ...
    'best_width_hz', nan(nComp, 1), ...
    'best_separation_hz', nan(nComp, 1), ...
    'best_trough_depth', nan(nComp, 1), ...
    'edge_ratio', nan(nComp, 1), ...
    'edge_run_score', nan(nComp, 1), ...
    'edge_artifact_flag', false(nComp, 1), ...
    'roughness_ratio', nan(nComp, 1), ...
    'roughness_penalty', nan(nComp, 1), ...
    'multi_peak_penalty', nan(nComp, 1), ...
    'multi_peak_flag', false(nComp, 1), ...
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
    rel_prom = 0.08 * max(y_pos);
    min_prom = max([0, rel_prom, powspctrm_form_prom_abs_floor, 0.15 * robust_scale]);
    [pks, locs, widths] = findpeaks(y_pos, x_band, 'MinPeakProminence', min_prom, 'MinPeakDistance', 5);
    if isempty(pks)
        continue;
    end
    pks = pks(:);
    locs = locs(:);
    widths = widths(:);
    [~, dom_idx] = max(pks);
    dom_amp = pks(dom_idx);
    dom_width = widths(dom_idx);
    dom_loc = locs(dom_idx);

    multi_peak_pen = 1;
    multi_peak_hit = false;
    if numel(pks) >= 2 && isfinite(dom_amp) && dom_amp > eps
        sep_ok = abs(locs - dom_loc) >= multi_peak_sep_min_hz;
        height_ok = pks >= multi_peak_height_ratio_min * dom_amp;
        rival_mask = ((1:numel(pks))' ~= dom_idx) & sep_ok & height_ok;
        if any(rival_mask)
            multi_peak_hit = true;
            n_rivals = sum(rival_mask);
            rival_amp_ratio = pks(rival_mask) ./ max(dom_amp, eps);
            amp_excess = max(0, mean(rival_amp_ratio) - multi_peak_height_ratio_min);
            amp_scale = amp_excess / max(1 - multi_peak_height_ratio_min, eps);
            extra_loss = 0.18 * max(0, n_rivals - 1) + 0.15 * min(1, amp_scale);
            multi_peak_pen = max(0.15, multi_peak_penalty_mult - extra_loss);
        end
    end

    dominance_pen = 1;
    n_rival_major = 0;
    n_rival_minor = 0;
    rival_amp_ratio_sum = 0;
    if numel(pks) >= 2 && isfinite(dom_amp) && dom_amp > eps
        sep_from_dom = abs(locs - dom_loc);
        is_rival = ((1:numel(pks))' ~= dom_idx);
        rival_ratios = pks ./ max(dom_amp, eps);
        rival_sep = sep_from_dom;

        major_ratio_min = 0.55;
        minor_ratio_min = 0.30;
        major_sep_min_hz = 4.0;
        minor_sep_min_hz = 2.5;

        major_mask = is_rival & (rival_sep >= major_sep_min_hz) & (rival_ratios >= major_ratio_min);
        minor_mask = is_rival & (rival_sep >= minor_sep_min_hz) & (rival_ratios >= minor_ratio_min);

        n_rival_major = sum(major_mask);
        n_rival_minor = sum(minor_mask);
        if any(major_mask)
            rival_amp_ratio_sum = sum(rival_ratios(major_mask));
        end

        clutter_load = 0;
        clutter_load = clutter_load + 0.35 * max(0, n_rival_major - 1);
        clutter_load = clutter_load + 0.12 * max(0, n_rival_minor - 2);
        clutter_load = clutter_load + 0.45 * max(0, rival_amp_ratio_sum - 0.85);

        max_loss = 0.55;
        dominance_pen = max(0.45, 1 - min(max_loss, clutter_load));
    end

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

    hf_pen = 1;
    hf_mask = x_band >= max(70, analysis_freq_range(2) - 15);
    if sum(hf_mask) >= 5
        hf_idx = find(hf_mask);
        hf_rho = corr((1:numel(hf_idx))', y_band(hf_idx)', 'rows', 'complete', 'type', 'Spearman');
        if isfinite(hf_rho) && hf_rho > 0.70
            hf_pen = max(0.65, 1 - 0.30 * (hf_rho - 0.70) / 0.30);
        end
    end

    [edge_ratio, edge_run_score, edge_flag] = compute_edge_artifact_indicators( ...
        y_resid, x_band, analysis_freq_range, edge_ratio_limit, edge_run_limit);
    diag.edge_ratio(ci) = edge_ratio;
    diag.edge_run_score(ci) = edge_run_score;
    diag.edge_artifact_flag(ci) = edge_flag;
    edge_pen = 1;
    if isfinite(edge_ratio) && edge_ratio > edge_ratio_soft
        edge_pen = edge_pen * max(0.80, 1 - 0.20 * min(1, (edge_ratio - edge_ratio_soft) / max(edge_ratio_limit - edge_ratio_soft, eps)));
    end
    if isfinite(edge_run_score) && edge_run_score > edge_run_soft
        edge_pen = edge_pen * max(0.80, 1 - 0.20 * min(1, (edge_run_score - edge_run_soft) / max(edge_run_limit - edge_run_soft, eps)));
    end
    if edge_flag
        edge_pen = edge_pen * 0.90;
    end

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
        rough_ref = 0.40;
        rough_span = 0.70;
        rough_pen_floor = 0.55;
        rough_pen_max_loss = 0.45;
        if isfinite(roughness_ratio) && roughness_ratio > rough_ref
            loss_frac = min(1, (roughness_ratio - rough_ref) / rough_span);
            roughness_pen = max(rough_pen_floor, 1 - rough_pen_max_loss * loss_frac);
        end
    end

    penalty_raw = edge_pen * dominance_pen * hf_pen * roughness_pen * multi_peak_pen;
    penalty_floor = 0.20;
    penalty_used = max(penalty_floor, penalty_raw);
    best_raw = best_pre_penalty * penalty_used;

    penalty_names = {'edge', 'dominance', 'hf_rise', 'roughness', 'multi_peak'};
    penalty_vals = [edge_pen, dominance_pen, hf_pen, roughness_pen, multi_peak_pen];
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
    diag.dominance_rival_count_major(ci) = n_rival_major;
    diag.dominance_rival_count_minor(ci) = n_rival_minor;
    diag.dominance_rival_amp_ratio_sum(ci) = rival_amp_ratio_sum;
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
    diag.multi_peak_penalty(ci) = multi_peak_pen;
    diag.multi_peak_flag(ci) = multi_peak_hit;
    diag.dominant_penalty_tag{ci} = dominant_penalty_tag;

    if ~isfinite(best_raw)
        powspctrm_form_score_vec(ci) = 0;
        powspctrm_form_mode_vec{ci} = 'none';
    else
        score_mapped = best_raw;
        if isfinite(min_similarity) && (min_similarity > 0) && (best_raw < min_similarity)
            score_mapped = 0.75 * score_mapped;
        end
        powspctrm_form_score_vec(ci) = max(0, min(1, score_mapped));
        powspctrm_form_mode_vec{ci} = mode_raw;
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

function [edge_ratio, edge_run_score, edge_artifact_flag] = compute_edge_artifact_indicators( ...
    y_resid, x_band, analysis_freq_range, edge_ratio_limit, edge_run_limit)
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
edge_artifact_flag = isfinite(edge_ratio) && (edge_ratio > edge_ratio_limit) && ...
    isfinite(edge_run_score) && (edge_run_score > edge_run_limit);
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
    pr_scan = trial_pr(ri, :);
    if all(~isfinite(pr_scan)) || ~any(band_mask)
        continue;
    end
    pr_band = pr_scan(band_mask);
    if all(~isfinite(pr_band))
        continue;
    end
    nonharm_val = nanmean(pr_scan(band_mask & ~harm_mask));
    harm_val = nanmean(pr_scan(band_mask & harm_mask));
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
    xh_all = log(f_hf(:));
    yh_all = log(abs(spec_hf(:)) + eps);
    valid_hf = isfinite(xh_all) & isfinite(yh_all) & ...
        isfinite(spec_hf(:)) & (spec_hf(:) > -inf);
    n_hf = nnz(valid_hf);
    if n_hf >= 3
        xh = xh_all(valid_hf);
        yh = yh_all(valid_hf);
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

function pr_scan = compute_simple_power_ratio_scan(x_stim, x_base, fs, scan_freqs, scan_width)
pr_scan = nan(size(scan_freqs));
if isempty(x_stim) || isempty(x_base) || fs <= 0
    return;
end
x_stim = x_stim(:)';
x_base = x_base(:)';
[p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(x_stim, x_base, fs, scan_freqs, scan_width);
if isempty(p_stim_scan) || isempty(p_base_scan)
    return;
end
pr_scan = p_stim_scan(1, :);
base_band_scan = p_base_scan(1, :);
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

function ged_freq_progress_reset(subject_label, subj_idx, subj_total, call_total_estimate)
global GED_FREQ_PROGRESS;
if isempty(subject_label)
    subject_label = 'NA';
end
if ~isfinite(subj_idx) || subj_idx < 0
    subj_idx = 0;
end
if ~isfinite(subj_total) || subj_total < 0
    subj_total = 0;
end
if ~isfinite(call_total_estimate) || call_total_estimate < 1
    call_total_estimate = 0;
end
GED_FREQ_PROGRESS = struct( ...
    'subject', subject_label, ...
    'subj_idx', round(subj_idx), ...
    'subj_total', round(subj_total), ...
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

function y = smooth_reflective(x, win)
y = x;
if isempty(win) || win <= 1 || isempty(x)
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
