%% GCP Trial-Level Narrowband Scanning GED (Common Spatial Filter)
%
% Identical spatial filter approach as GCP_eeg_GED_subjects.m, but:
%   - Peak detection is done on INDIVIDUAL TRIALS rather than
%     condition-averaged spectra
%
% Pipeline:
%   Phase 1 — Pool all conditions -> broadband GED -> common spatial filter
%   Phase 2 — For each condition & trial, narrowband scan (30-90 Hz),
%             compute per-trial power ratio, detrend, detect peaks.
%             Single-peak and dual-peak (hard 50 Hz) models per trial.
%   Method comparison — Three dual-peak assignment methods:
%             1. Hard 50 Hz boundary
%             2. Per-subject spectral-trough boundary
%             3. Gaussian Mixture Model on single-peak frequencies
%   Aggregation — mean & median peak frequency per condition per subject

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');

nSubj = length(subjects);

% Time windows
baseline_window = [-1.5, -0.25];
stimulus_window = [0.3, 2.0];

% Gamma frequency range
gamma_range = [30, 90];

% Narrowband scanning parameters
scan_freqs = 30:1:90;
nFreqs     = length(scan_freqs);
scan_width = 4;

% GED parameters
lambda = 0.01;
ged_search_n = 10;          % search first N GED components
topo_display_n = 10;        % number of component topos to visualize
template_front_weight = 0.7; % anti-template weight for frontal channels
template_sigma_occ = 0.20;   % spatial smoothness for occipital template
template_sigma_front = 0.25; % spatial smoothness for frontal anti-template
score_w_corr = 1.00;        % composite-score weight for template correlation
score_w_gamma = 1.00;       % composite-score weight for pooled gamma evidence
score_w_frontleak = 1.00;   % composite-score penalty for frontal leakage

% Detrending parameters for power-ratio spectrum
poly_order = 2;

% Condition info
condNames  = {'c25', 'c50', 'c75', 'c100'};
trialCodes = [61, 62, 63, 64];
condLabels = {'25%', '50%', '75%', '100%'};

% Figure save directories
if ispc
    fig_save_dir = 'W:\Students\Arne\GCP\figures\eeg\ged\trials';
else
    fig_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged/trials';
end
if ~exist(fig_save_dir, 'dir'), mkdir(fig_save_dir); end

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

all_trial_mean_single   = nan(4, nSubj);
all_trial_median_single = nan(4, nSubj);
all_trial_mean_low      = nan(4, nSubj);
all_trial_median_low    = nan(4, nSubj);
all_trial_mean_high     = nan(4, nSubj);
all_trial_median_high   = nan(4, nSubj);

all_trial_detrate_single = nan(4, nSubj);
all_trial_detrate_low    = nan(4, nSubj);
all_trial_detrate_high   = nan(4, nSubj);

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

chanlocs_all = {};

%% Process each subject
for subj = 1:nSubj
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

    occ_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I)', 'once')), dataEEG_c25.label);
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

    %% ================================================================
    %  PHASE 1: Build POOLED covariance across all conditions -> one GED
    %  ================================================================
    clc
    fprintf('Subject %s (%d/%d) — Phase 1: Full-head GED + occipital template (%d occ / %d ch)\n', ...
        subjects{subj}, subj, nSubj, nOcc, nChans);

    covStim_full = zeros(nChans);
    covBase_full = zeros(nChans);
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

        cfg_t.latency = stimulus_window;
        dat_stim = ft_selectdata(cfg_t, dat_gamma);

        nTrl = length(dat_stim.trial);
        for trl = 1:nTrl
            d = double(dat_stim.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            covStim_full = covStim_full + (d * d') / size(d, 2);

            d = double(dat_base.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            covBase_full = covBase_full + (d * d') / size(d, 2);
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

    % Pareto-filter candidates on three objectives (maximize corr, gamma evidence, ratio),
    % then select by composite score z(corr) + z(gamma) - z(frontLeak).
    corr_vec = searchCorrs;
    ratio_vec = searchOccFrontRatio;
    gamma_vec = searchGammaEvidence;
    leak_vec = searchFrontLeak;

    pareto_mask = true(nSearch, 1);
    for i = 1:nSearch
        for j = 1:nSearch
            if i == j, continue; end
            better_or_equal = (corr_vec(j) >= corr_vec(i)) && ...
                              (gamma_vec(j) >= gamma_vec(i)) && ...
                              (ratio_vec(j) >= ratio_vec(i));
            strictly_better = (corr_vec(j) > corr_vec(i)) || ...
                              (gamma_vec(j) > gamma_vec(i)) || ...
                              (ratio_vec(j) > ratio_vec(i));
            if better_or_equal && strictly_better
                pareto_mask(i) = false;
                break;
            end
        end
    end

    valid_corr = isfinite(corr_vec);
    valid_gamma = isfinite(gamma_vec);
    valid_leak = isfinite(leak_vec);

    if any(valid_corr), corr_mu = mean(corr_vec(valid_corr)); corr_sd = std(corr_vec(valid_corr));
    else, corr_mu = 0; corr_sd = 1; end
    if any(valid_gamma), gamma_mu = mean(gamma_vec(valid_gamma)); gamma_sd = std(gamma_vec(valid_gamma));
    else, gamma_mu = 0; gamma_sd = 1; end
    if any(valid_leak), leak_mu = mean(leak_vec(valid_leak)); leak_sd = std(leak_vec(valid_leak));
    else, leak_mu = 0; leak_sd = 1; end

    z_corr = (corr_vec - corr_mu) / max(corr_sd, eps);
    z_gamma = (gamma_vec - gamma_mu) / max(gamma_sd, eps);
    z_leak = (leak_vec - leak_mu) / max(leak_sd, eps);

    comp_score = score_w_corr * z_corr + score_w_gamma * z_gamma - score_w_frontleak * z_leak;
    comp_score(~pareto_mask) = -Inf;
    if ~any(isfinite(comp_score))
        warning('Pareto filter removed all candidates for subject %s; falling back to unconstrained score.', subjects{subj});
        comp_score = score_w_corr * z_corr + score_w_gamma * z_gamma - score_w_frontleak * z_leak;
    end

    searchScores = comp_score;
    [bestScore, bestIdx] = max(searchScores);
    if isempty(bestIdx) || isnan(bestScore)
        bestIdx = 1;
        bestScore = NaN;
    end
    bestCorr = searchCorrs(bestIdx);
    bestOcc = searchOccStrength(bestIdx);
    bestFront = searchFrontStrength(bestIdx);
    bestRatio = searchOccFrontRatio(bestIdx);
    bestGamma = searchGammaEvidence(bestIdx);
    bestLeak = searchFrontLeak(bestIdx);

    topComp = searchFilters(:, bestIdx);
    topo_temp = covStim_full_reg * topComp;

    [~, topDispOrder] = sort(searchScores, 'descend');
    nDispTopo = min(topo_display_n, nSearch);
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
    all_eigenvalues(subj) = evals_sorted(1);
    all_selected_comp_idx(subj)  = bestIdx;
    all_selected_comp_corr(subj) = bestCorr;
    all_selected_comp_eval(subj) = evals_sorted(bestIdx);
    all_top5_corrs(1:nStore, subj) = storeCorrs;
    all_top5_evals(1:nStore, subj) = storeEvals;
    all_top5_topos{subj} = storeTopos;
    all_simulated_templates{subj} = sim_template;

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
    sgtitle(sprintf(['Subject %s: GED Top %d (selected C%d, ', ...
        'score=%.3f, r=%.3f, occ=%.3f, front=%.3f, occ/front=%.3f, gamma=%.3f, leak=%.3f)'], ...
        subjects{subj}, nDispTopo, bestIdx, bestScore, bestCorr, bestOcc, bestFront, bestRatio, bestGamma, bestLeak), ...
        'FontSize', 14, 'FontWeight', 'bold');
    nColsTopo = 5;
    nRowsTopo = ceil(nDispTopo / nColsTopo);
    for ci = 1:nDispTopo
        subplot(nRowsTopo, nColsTopo, ci);
        topo_data = [];
        topo_data.label  = all_topo_labels{subj};
        topo_data.avg    = topoDispTopos(:, ci);
        topo_data.dimord = 'chan';
        try
            ft_topoplotER(cfg_topo, topo_data);
            colorbar;
        catch
            imagesc(topoDispTopos(:, ci)); colorbar;
        end
        comp_rank = dispTopoIdx(ci);
        title(sprintf('C%d: score=%.2f, r=%.2f, ratio=%.2f, g=%.2f', ...
            comp_rank, topoDispScores(ci), topoDispCorrs(ci), ...
            searchOccFrontRatio(comp_rank), searchGammaEvidence(comp_rank)), ...
            'FontSize', 7);
    end
    saveas(fig_sel, fullfile(fig_save_dir, sprintf('GCP_eeg_GED_trials_top5_selection_subj%s.png', subjects{subj})));

    %% ================================================================
    %  PHASE 2: Per condition — trial-level narrowband scanning
    %  ================================================================
    for cond = 1:4

        dat = dat_per_cond{cond};
        if isempty(dat), continue; end

        nTrl = length(dat.trial);
        powratio_trials = nan(nTrl, nFreqs);

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
                comp_stim = topComp' * double(dat_stim_nb.trial{trl});
                comp_base = topComp' * double(dat_base_nb.trial{trl});
                pow_stim = mean(comp_stim.^2);
                pow_base = mean(comp_base.^2);
                if pow_base > 0
                    powratio_trials(trl, fi) = pow_stim / pow_base;
                end
            end
        end

        all_trial_powratio{cond, subj} = powratio_trials;

        %% Per-trial peak detection
        trl_peaks_single = nan(nTrl, 1);
        trl_peaks_low    = nan(nTrl, 1);
        trl_peaks_high   = nan(nTrl, 1);

        for trl = 1:nTrl
            pr = powratio_trials(trl, :);
            if all(isnan(pr)), continue; end

            p = polyfit(scan_freqs, pr, poly_order);
            pr_dt = pr - polyval(p, scan_freqs);
            pr_dt_smooth = movmean(pr_dt, 5);

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
    all_vals = cell2mat(cellfun(@(x) x(:), pr_dt_mats(~cellfun(@isempty, pr_dt_mats)), ...
        'UniformOutput', false)');
    if ~isempty(all_vals)
        cl_shared = prctile(abs(all_vals(~isnan(all_vals))), 98);
    else
        cl_shared = 1;
    end

    % --- Row 1: Heatmap of trial-level detrended power-ratio spectra ---
    for cond = 1:4
        subplot(4, 4, cond);
        if ~isempty(pr_dt_mats{cond})
            imagesc(scan_freqs, 1:size(pr_dt_mats{cond},1), pr_dt_mats{cond});
            colormap(gca, cmap_div);
            caxis([-cl_shared cl_shared]);
            cb = colorbar; cb.FontSize = 8;
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
        end
        xlabel('Freq [Hz]'); ylabel('\Delta PR');
        det = all_trial_detrate_single(cond, subj);
        title(sprintf('%s Single (det=%.0f%%)', condLabels{cond}, det*100), 'FontSize', 10);
        set(gca, 'FontSize', 10); xlim([30 90]); grid on; box on;
    end

    % --- Row 3: Mean trial-level spectrum with dual-peak markers (50 Hz) ---
    for cond = 1:4
        subplot(4, 4, 8 + cond); hold on;
        if ~isempty(pr_dt_mats{cond})
            mu_dt = nanmean(pr_dt_mats{cond}, 1);
            nTrl = size(pr_dt_mats{cond}, 1);
            sem_dt = nanstd(pr_dt_mats{cond}, [], 1) / sqrt(nTrl);
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
        end
        xlabel('Freq [Hz]'); ylabel('\Delta PR');
        det_lo = all_trial_detrate_low(cond, subj);
        det_hi = all_trial_detrate_high(cond, subj);
        title(sprintf('%s Dual (L:%.0f%% H:%.0f%%)', condLabels{cond}, det_lo*100, det_hi*100), ...
            'FontSize', 10);
        set(gca, 'FontSize', 10); xlim([30 90]); grid on; box on;
    end

    % --- Row 4: Topoplot + histogram ---
    occ_highlight = dataEEG_c25.label(cellfun(@(l) ...
        ~isempty(regexp(l, '[OI]', 'once')), dataEEG_c25.label));

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
        topo_data.avg    = all_topos{subj};
        topo_data.dimord = 'chan';
        try
            ft_topoplotER(cfg_topo, topo_data);
            cb = colorbar; cb.FontSize = 9;
        catch
            imagesc(topo_data.avg); colorbar;
        end
        title(sprintf('Common Filter (\\lambda=%.1f)', all_eigenvalues(subj)), 'FontSize', 11);
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
    legend(condLabels, 'FontSize', 10, 'Location', 'best');
    set(gca, 'FontSize', 11); xlim([30 90]); grid on; box on;

    saveas(fig, fullfile(fig_save_dir, sprintf('GCP_eeg_GED_trials_subj%s.png', subjects{subj})));
    toc
end % subject loop

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
%  METHOD COMPARISON FIGURES
%  ====================================================================
close all

for mi = 1:3
    save_dir = method_dirs{mi};

    %% --- Subject-level median dual-peak: raincloud ---
    fig_box = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('%s: Dual Peak Frequency (Median)', method_names{mi}), ...
        'FontSize', 18, 'FontWeight', 'bold');

    dual_data   = {m_median_low(:,:,mi), m_median_high(:,:,mi)};
    dual_titles = {'Low Gamma', 'High Gamma'};
    dual_ylims  = {[30 55], [45 90]};

    for di = 1:2
        subplot(1, 2, di); hold on;
        peak_data = dual_data{di};

        for s = 1:nSubj
            pf = peak_data(:, s);
            if sum(~isnan(pf)) >= 2
                plot(1:4, pf, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
            end
        end

        for c = 1:4
            pf = peak_data(c, :);
            pf = pf(~isnan(pf));
            if length(pf) >= 2
                [f_dens, xi] = ksdensity(pf);
                f_dens = f_dens / max(f_dens) * 0.3;
                patch(c - f_dens - 0.05, xi, colors(c,:), ...
                    'FaceAlpha', 0.3, 'EdgeColor', colors(c,:), 'LineWidth', 1);
            end
        end

        y_box = peak_data(:);
        g_box = repelem((1:4)', nSubj, 1);
        valid_box = ~isnan(y_box);
        if any(valid_box)
            boxplot(y_box(valid_box), g_box(valid_box), 'Colors', 'k', ...
                'Symbol', '', 'Widths', 0.15);
        end

        hold on;
        for c = 1:4
            pf = peak_data(c, :);
            pf = pf(~isnan(pf));
            xJit = c + 0.15 + (rand(size(pf)) - 0.5) * 0.1;
            scatter(xJit, pf, 200, colors(c,:), 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end

        xlim([0.3 4.7]); ylim(dual_ylims{di});
        set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 16, 'Box', 'off');
        ylabel('Peak Gamma Frequency [Hz]');
        title(dual_titles{di}, 'FontSize', 16, 'FontWeight', 'bold');
    end

    saveas(fig_box, fullfile(save_dir, 'boxplot_DualGamma.png'));

    %% --- All-trials pooled dual-peak: raincloud ---
    fig_trl = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('%s: All Trials (pooled)', method_names{mi}), ...
        'FontSize', 18, 'FontWeight', 'bold');

    dual_trl_data   = {m_peaks_low{mi}, m_peaks_high{mi}};
    dual_trl_titles = {'Low Gamma', 'High Gamma'};
    dual_trl_ylims  = {[30 55], [45 90]};

    for di = 1:2
        subplot(1, 2, di); hold on;

        y_all = []; g_all = [];
        for c = 1:4
            for s = 1:nSubj
                tpk = dual_trl_data{di}{c, s};
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
            vals = y_all(g_all == c);
            xJit = c + 0.15 + (rand(size(vals)) - 0.5) * 0.25;
            scatter(xJit, vals, 12, colors(c,:), 'filled', 'MarkerFaceAlpha', 0.2);
        end

        xlim([0.3 4.7]); ylim(dual_trl_ylims{di});
        set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 16, 'Box', 'off');
        ylabel('Peak Gamma Frequency [Hz]');
        title(dual_trl_titles{di}, 'FontSize', 16, 'FontWeight', 'bold');
    end

    saveas(fig_trl, fullfile(save_dir, 'boxplot_alltrials_DualGamma.png'));

    %% --- Detection rate ---
    fig_det_m = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('%s: Detection Rate', method_names{mi}), ...
        'FontSize', 18, 'FontWeight', 'bold');

    det_data_m   = {all_trial_detrate_single, m_detrate_low(:,:,mi), m_detrate_high(:,:,mi)};
    det_labels_m = {'Single Peak', 'Low Gamma', 'High Gamma'};

    for di = 1:3
        subplot(1, 3, di); hold on;
        dr = det_data_m{di};

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
        title(det_labels_m{di}, 'FontSize', 16, 'FontWeight', 'bold');
        grid on;
    end

    saveas(fig_det_m, fullfile(save_dir, 'detection_rate.png'));
    close all
end

%% ====================================================================
%  GRAND AVERAGE: Single Peak boxplots (mean & median) — raincloud
%  ====================================================================
fprintf('\nCreating grand average figures...\n');

fig_box1 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Trial-Level Peak Frequency: Single Peak', 'FontSize', 18, 'FontWeight', 'bold');

agg_data   = {all_trial_mean_single, all_trial_median_single};
agg_titles = {'Mean over Trials', 'Median over Trials'};

for ai = 1:2
    subplot(1, 2, ai); hold on;
    peak_data = agg_data{ai};

    for s = 1:nSubj
        pf = peak_data(:, s);
        if sum(~isnan(pf)) >= 2
            plot(1:4, pf, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
        end
    end

    for c = 1:4
        pf = peak_data(c, :);
        pf = pf(~isnan(pf));
        if length(pf) >= 2
            [f_dens, xi] = ksdensity(pf);
            f_dens = f_dens / max(f_dens) * 0.3;
            patch(c - f_dens - 0.05, xi, colors(c,:), ...
                'FaceAlpha', 0.3, 'EdgeColor', colors(c,:), 'LineWidth', 1);
        end
    end

    y_box = peak_data(:);
    g_box = repelem((1:4)', nSubj, 1);
    valid_box = ~isnan(y_box);
    if any(valid_box)
        boxplot(y_box(valid_box), g_box(valid_box), 'Colors', 'k', ...
            'Symbol', '', 'Widths', 0.15);
    end

    hold on;
    for c = 1:4
        pf = peak_data(c, :);
        pf = pf(~isnan(pf));
        xJit = c + 0.15 + (rand(size(pf)) - 0.5) * 0.1;
        scatter(xJit, pf, 200, colors(c,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end

    xlim([0.3 4.7]); ylim([30 90]);
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 16, 'Box', 'off');
    ylabel('Peak Gamma Frequency [Hz]');
    title(agg_titles{ai}, 'FontSize', 16, 'FontWeight', 'bold');
end

saveas(fig_box1, fullfile(fig_save_dir, 'GCP_eeg_GED_trials_boxplot_SinglePeak.png'));

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

saveas(fig_trl1, fullfile(fig_save_dir, 'GCP_eeg_GED_trials_boxplot_alltrials_SinglePeak.png'));

%% ====================================================================
%  GRAND AVERAGE: All-subjects subplot (mean trial spectra)
%  ====================================================================
nRows = ceil(nSubj / 5);
fig_all = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Trial-Level Mean Detrended Spectra: All Subjects (N=%d)', nSubj), ...
    'FontSize', 16, 'FontWeight', 'bold');

for s = 1:nSubj
    subplot(nRows, 5, s); hold on;
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
            plot(scan_freqs, movmean(mu_dt, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2);
            md_pf = all_trial_median_single(cond, s);
            if ~isnan(md_pf)
                xline(md_pf, '--', 'Color', colors(cond,:), 'LineWidth', 1.2);
            end
        end
    end
    yline(0, 'k-', 'LineWidth', 0.5);
    xlabel('Freq [Hz]'); ylabel('\Delta PR');
    title(sprintf('Subj %s', subjects{s}), 'FontSize', 11);
    set(gca, 'FontSize', 10); xlim([30 90]); grid on; box on;
    if s == 1
        legend(condLabels, 'FontSize', 8, 'Location', 'best');
    end
end
saveas(fig_all, fullfile(fig_save_dir, 'GCP_eeg_GED_trials_all_subjects.png'));

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
    grid on;
end

saveas(fig_det, fullfile(fig_save_dir, 'GCP_eeg_GED_trials_detection_rate.png'));

%% Save results
if ispc
    save_path = 'W:\Students\Arne\GCP\data\features\GCP_eeg_GED_trials.mat';
else
    save_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_eeg_GED_trials.mat';
end
save(save_path, ...
    'all_trial_powratio', ...
    'all_trial_peaks_single', 'all_trial_peaks_low', 'all_trial_peaks_high', ...
    'all_trial_mean_single', 'all_trial_median_single', ...
    'all_trial_mean_low', 'all_trial_median_low', ...
    'all_trial_mean_high', 'all_trial_median_high', ...
    'all_trial_detrate_single', 'all_trial_detrate_low', 'all_trial_detrate_high', ...
    'all_topos', 'all_topo_labels', 'all_eigenvalues', ...
    'all_selected_comp_idx', 'all_selected_comp_corr', 'all_selected_comp_eval', ...
    'all_top5_corrs', 'all_top5_evals', 'all_top5_topos', 'all_simulated_templates', ...
    'm_peaks_low', 'm_peaks_high', 'm_median_low', 'm_median_high', ...
    'm_detrate_low', 'm_detrate_high', ...
    'trough_freqs', 'gmm_means', 'gmm_boundary', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames');

clc
fprintf('Done.\n');
