%% GCP Trial-Level Narrowband Scanning GED (Common Spatial Filter)
%
% Identical spatial filter approach as GCP_eeg_ged.m, but:
%   - Uses only POSTERIOR channels (labels containing O, I, or P, excluding
%     those with F or T) for the GED
%   - Peak detection is done on INDIVIDUAL TRIALS rather than
%     condition-averaged spectra
%
% Pipeline:
%   Phase 1 — Pool all conditions -> broadband GED -> common spatial filter
%             (posterior channels only)
%   Phase 2 — For each condition & trial, narrowband scan (30-90 Hz),
%             compute per-trial power ratio, detrend, detect peaks.
%             Both single-peak and dual-peak models are applied per trial.
%   Aggregation — mean & median peak frequency per condition per subject
%
% Rationale: averaging spectra before peak detection can smear out peaks
% when gamma frequency fluctuates trial-to-trial. Detecting peaks per
% trial and then averaging preserves this variability.

clear; close all; clc

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

% Detrending parameters for power-ratio spectrum
poly_order = 2;

% Condition info
condNames  = {'c25', 'c50', 'c75', 'c100'};
trialCodes = [61, 62, 63, 64];
condLabels = {'25%', '50%', '75%', '100%'};

% Figure save directory
if ispc
    fig_save_dir = 'W:\Students\Arne\GCP\figures\eeg\ged\trials';
else
    fig_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged/trials';
end
if ~exist(fig_save_dir, 'dir'), mkdir(fig_save_dir); end

%% Preallocate storage
% Trial-level power-ratio matrices (trials x freqs)
all_trial_powratio = cell(4, nSubj);

% Trial-level peak frequencies (vectors, one entry per trial)
all_trial_peaks_single = cell(4, nSubj);
all_trial_peaks_low    = cell(4, nSubj);
all_trial_peaks_high   = cell(4, nSubj);

% Aggregated: mean & median per condition per subject
all_trial_mean_single   = nan(4, nSubj);
all_trial_median_single = nan(4, nSubj);
all_trial_mean_low      = nan(4, nSubj);
all_trial_median_low    = nan(4, nSubj);
all_trial_mean_high     = nan(4, nSubj);
all_trial_median_high   = nan(4, nSubj);

% Detection rates (fraction of trials with a valid peak)
all_trial_detrate_single = nan(4, nSubj);
all_trial_detrate_low    = nan(4, nSubj);
all_trial_detrate_high   = nan(4, nSubj);

% Common filter info (one per subject)
all_topos       = cell(1, nSubj);
all_topo_labels = cell(1, nSubj);
all_eigenvalues = nan(1, nSubj);

% Channel labels
chanlocs_all = {};

%% Process each subject
for subj = 1:nSubj

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

    % Select posterior channels (containing O, I, or P, but not F or T)
    all_labels = dataEEG_c25.label;
    occ_idx = find(cellfun(@(l) ...
        ~isempty(regexp(l, '[OIP]', 'once')) && isempty(regexp(l, '[FT]', 'once')), ...
        all_labels));
    occ_labels = all_labels(occ_idx);
    nChans = length(occ_idx);

    %% ================================================================
    %  PHASE 1: Build POOLED covariance across all conditions -> one GED
    %  ================================================================
    clc
    fprintf('Subject %s (%d/%d) — Phase 1: Pooled GED (posterior, %d ch)\n', ...
        subjects{subj}, subj, nSubj, nChans);

    covStim_pooled = zeros(nChans);
    covBase_pooled = zeros(nChans);
    nTrials_total  = 0;

    dat_per_cond = cell(1, 4);

    for cond = 1:4
        dat    = dataStructs{cond};
        trlIdx = trialIndices{cond};
        if isempty(trlIdx), continue; end

        cfg = [];
        cfg.trials  = trlIdx;
        cfg.channel = occ_labels;
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
            covStim_pooled = covStim_pooled + (d * d') / size(d, 2);

            d = double(dat_base.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            covBase_pooled = covBase_pooled + (d * d') / size(d, 2);
        end
        nTrials_total = nTrials_total + nTrl;
    end

    covStim_pooled = covStim_pooled / nTrials_total;
    covBase_pooled = covBase_pooled / nTrials_total;

    % Shrinkage regularization
    covStim_pooled = (1-lambda)*covStim_pooled + lambda*mean(diag(covStim_pooled))*eye(nChans);
    covBase_pooled = (1-lambda)*covBase_pooled + lambda*mean(diag(covBase_pooled))*eye(nChans);

    % GED on pooled covariance
    [W, D] = eig(covStim_pooled, covBase_pooled);
    [evals_sorted, sortIdx] = sort(real(diag(D)), 'descend');
    W = W(:, sortIdx);

    topComp = W(:, 1);
    all_eigenvalues(subj) = evals_sorted(1);

    % Sign correction
    topo_temp = covStim_pooled * topComp;
    [~, mxI] = max(abs(topo_temp));
    if topo_temp(mxI) < 0, topComp = -topComp; end

    all_topos{subj}     = covStim_pooled * topComp;
    all_topo_labels{subj} = occ_labels;

    %% ================================================================
    %  PHASE 2: Per condition — trial-level narrowband scanning
    %  ================================================================
    for cond = 1:4

        dat = dat_per_cond{cond};
        if isempty(dat), continue; end

        nTrl = length(dat.trial);

        %% Build trial x freq power-ratio matrix
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

            % Detrend
            p = polyfit(scan_freqs, pr, poly_order);
            pr_dt = pr - polyval(p, scan_freqs);
            pr_dt_smooth = movmean(pr_dt, 5);

            %% Single-peak detection
            [pks, locs] = findpeaks(pr_dt_smooth, scan_freqs, ...
                'MinPeakProminence', max(pr_dt_smooth) * 0.15, ...
                'MinPeakDistance', 5);

            if ~isempty(pks)
                [~, best_pk] = max(pks);
                trl_peaks_single(trl) = locs(best_pk);
            end

            %% Dual-peak detection
            [pks_all, locs_all] = findpeaks(pr_dt_smooth, scan_freqs, ...
                'MinPeakDistance', 5);
            pos_mask  = pks_all > 0;
            pks_pos   = pks_all(pos_mask);
            locs_pos  = locs_all(pos_mask);
            nPosPeaks = length(pks_pos);

            if nPosPeaks >= 2
                in_lo = locs_pos >= 30 & locs_pos <= 55;
                in_hi = locs_pos >= 50 & locs_pos <= 90;
                lo_freq = NaN; hi_freq = NaN;
                if any(in_lo)
                    [~, bi] = max(pks_pos(in_lo));
                    tmp = locs_pos(in_lo);
                    lo_freq = tmp(bi);
                end
                if any(in_hi)
                    [~, bi] = max(pks_pos(in_hi));
                    tmp = locs_pos(in_hi);
                    hi_freq = tmp(bi);
                end
                % Resolve overlap-zone conflict
                if ~isnan(lo_freq) && ~isnan(hi_freq) && lo_freq == hi_freq
                    lo_pks = pks_pos(in_lo);
                    lo_lcs = locs_pos(in_lo);
                    hi_pks = pks_pos(in_hi);
                    hi_lcs = locs_pos(in_hi);
                    alt_lo_mask = lo_lcs ~= lo_freq;
                    alt_hi_mask = hi_lcs ~= hi_freq;
                    has_alt_lo = any(alt_lo_mask);
                    has_alt_hi = any(alt_hi_mask);
                    if has_alt_lo && has_alt_hi
                        [~, ai] = max(lo_pks(alt_lo_mask));
                        tmp2 = lo_lcs(alt_lo_mask);
                        lo_freq = tmp2(ai);
                    elseif has_alt_lo
                        [~, ai] = max(lo_pks(alt_lo_mask));
                        tmp2 = lo_lcs(alt_lo_mask);
                        lo_freq = tmp2(ai);
                    elseif has_alt_hi
                        [~, ai] = max(hi_pks(alt_hi_mask));
                        tmp2 = hi_lcs(alt_hi_mask);
                        hi_freq = tmp2(ai);
                    end
                end
                trl_peaks_low(trl)  = lo_freq;
                trl_peaks_high(trl) = hi_freq;
            elseif nPosPeaks == 1
                the_peak = locs_pos(1);
                if the_peak >= 45 && the_peak <= 65
                    trl_peaks_low(trl)  = the_peak;
                    trl_peaks_high(trl) = the_peak;
                elseif the_peak < 45
                    trl_peaks_low(trl) = the_peak;
                else
                    trl_peaks_high(trl) = the_peak;
                end
            end
        end

        %% Store trial-level peaks
        all_trial_peaks_single{cond, subj} = trl_peaks_single;
        all_trial_peaks_low{cond, subj}    = trl_peaks_low;
        all_trial_peaks_high{cond, subj}   = trl_peaks_high;

        %% Aggregate: mean, median, detection rate
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
    sgtitle(sprintf('Trial-Level GED (posterior): Subject %s', subjects{subj}), ...
        'FontSize', 20, 'FontWeight', 'bold');

    cmap_div = interp1([0 0.5 1], ...
        [0.17 0.27 0.53; 0.97 0.97 0.97; 0.70 0.09 0.17], linspace(0,1,256));

    % --- Row 1: Heatmap of trial-level detrended power-ratio spectra ---
    % Pre-compute detrended matrices and find a common color limit
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
        title(sprintf('%s Detrended', condLabels{cond}), 'FontSize', 12);
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
        title(sprintf('%s Single (det=%.0f%%)', condLabels{cond}, det*100), 'FontSize', 11);
        set(gca, 'FontSize', 10); xlim([30 90]); grid on; box on;
    end

    % --- Row 3: Mean trial-level spectrum with dual-peak markers ---
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
            xline(50, 'k:', 'LineWidth', 0.8, 'Alpha', 0.4);
            xline(55, 'k:', 'LineWidth', 0.8, 'Alpha', 0.4);
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
            'FontSize', 11);
        set(gca, 'FontSize', 10); xlim([30 90]); grid on; box on;
    end

    % --- Row 4: Topoplot + histogram of trial-level peak frequencies ---
    cfg_topo = [];
    cfg_topo.layout    = headmodel.layANThead;
    cfg_topo.comment   = 'no';
    cfg_topo.marker    = 'off';
    cfg_topo.style     = 'straight';
    cfg_topo.gridscale = 300;
    cfg_topo.zlim      = 'maxabs';
    cfg_topo.colormap  = '*RdBu';
    cfg_topo.figure    = 'gcf';

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
        title(sprintf('Common Filter (\\lambda=%.1f)', all_eigenvalues(subj)), 'FontSize', 12);
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
    title('Trial-Level Single Peak Distribution', 'FontSize', 14);
    legend(condLabels, 'FontSize', 11, 'Location', 'best');
    set(gca, 'FontSize', 12); xlim([30 90]); grid on; box on;

    saveas(fig, fullfile(fig_save_dir, sprintf('GED_trials_subj%s.png', subjects{subj})));

end % subject loop

%% ====================================================================
%  GRAND AVERAGE FIGURE 1: Boxplots — Single Peak (mean & median)
%  ====================================================================
fprintf('\nCreating grand average figures...\n');
close all

fig_box1 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Trial-Level Peak Frequency: Single Peak', 'FontSize', 26, 'FontWeight', 'bold');

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

    y_box = peak_data(:);
    g_box = repelem((1:4)', nSubj, 1);
    valid_box = ~isnan(y_box);
    if any(valid_box)
        boxplot(y_box(valid_box), g_box(valid_box), 'Colors', 'k', 'Symbol', '');
    end

    hold on;
    for c = 1:4
        pf = peak_data(c, :);
        pf = pf(~isnan(pf));
        xJit = c + (rand(size(pf)) - 0.5) * 0.1;
        scatter(xJit, pf, 250, colors(c,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end

    xlim([0.5 4.5]); ylim([30 90]);
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 18, 'Box', 'off');
    ylabel('Peak Gamma Frequency [Hz]');
    title(agg_titles{ai}, 'FontSize', 22, 'FontWeight', 'bold');
end

saveas(fig_box1, fullfile(fig_save_dir, 'GED_trials_boxplot_SinglePeak.png'));

%% ====================================================================
%  GRAND AVERAGE FIGURE 2: Boxplots — Dual Peak (median)
%  ====================================================================
fig_box2 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Trial-Level Peak Frequency: Dual Peak (Median)', 'FontSize', 26, 'FontWeight', 'bold');

dual_data   = {all_trial_median_low, all_trial_median_high};
dual_titles = {'Low Gamma (Median)', 'High Gamma (Median)'};
dual_ylims  = {[30 65], [40 90]};

for di = 1:2
    subplot(1, 2, di); hold on;
    peak_data = dual_data{di};

    for s = 1:nSubj
        pf = peak_data(:, s);
        if sum(~isnan(pf)) >= 2
            plot(1:4, pf, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
        end
    end

    y_box = peak_data(:);
    g_box = repelem((1:4)', nSubj, 1);
    valid_box = ~isnan(y_box);
    if any(valid_box)
        boxplot(y_box(valid_box), g_box(valid_box), 'Colors', 'k', 'Symbol', '');
    end

    hold on;
    for c = 1:4
        pf = peak_data(c, :);
        pf = pf(~isnan(pf));
        xJit = c + (rand(size(pf)) - 0.5) * 0.1;
        scatter(xJit, pf, 250, colors(c,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end

    xlim([0.5 4.5]); ylim(dual_ylims{di});
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 18, 'Box', 'off');
    ylabel('Peak Gamma Frequency [Hz]');
    title(dual_titles{di}, 'FontSize', 22, 'FontWeight', 'bold');
end

saveas(fig_box2, fullfile(fig_save_dir, 'GED_trials_boxplot_DualGamma.png'));

%% ====================================================================
%  TRIAL-LEVEL BOXPLOTS: All individual trials pooled across subjects
%  ====================================================================

% --- Single Peak: all trials ---
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

if ~isempty(y_all)
    boxplot(y_all, g_all, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
end
hold on;
for c = 1:4
    mask = g_all == c;
    vals = y_all(mask);
    xJit = c + (rand(size(vals)) - 0.5) * 0.35;
    scatter(xJit, vals, 15, colors(c,:), 'filled', 'MarkerFaceAlpha', 0.25);
end

xlim([0.5 4.5]); ylim([30 90]);
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 20, 'Box', 'off');
ylabel('Peak Gamma Frequency [Hz]');
title('Single Peak: All Trials (pooled across subjects)', ...
    'FontSize', 26, 'FontWeight', 'bold');

saveas(fig_trl1, fullfile(fig_save_dir, 'GED_trials_boxplot_alltrials_SinglePeak.png'));

% --- Dual Peak: Low + High gamma, all trials ---
fig_trl2 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Dual Peak: All Trials (pooled across subjects)', ...
    'FontSize', 26, 'FontWeight', 'bold');

dual_trl_data   = {all_trial_peaks_low, all_trial_peaks_high};
dual_trl_titles = {'Low Gamma', 'High Gamma'};
dual_trl_ylims  = {[30 65], [40 90]};

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

    if ~isempty(y_all)
        boxplot(y_all, g_all, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
    end
    hold on;
    for c = 1:4
        mask = g_all == c;
        vals = y_all(mask);
        xJit = c + (rand(size(vals)) - 0.5) * 0.35;
        scatter(xJit, vals, 15, colors(c,:), 'filled', 'MarkerFaceAlpha', 0.25);
    end

    xlim([0.5 4.5]); ylim(dual_trl_ylims{di});
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 18, 'Box', 'off');
    ylabel('Peak Gamma Frequency [Hz]');
    title(dual_trl_titles{di}, 'FontSize', 22, 'FontWeight', 'bold');
end

saveas(fig_trl2, fullfile(fig_save_dir, 'GED_trials_boxplot_alltrials_DualGamma.png'));

%% ====================================================================
%  GRAND AVERAGE FIGURE 3: All-subjects subplot (mean trial spectra)
%  ====================================================================
nRows = ceil(nSubj / 5);
fig_all = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Trial-Level Mean Detrended Spectra: All Subjects (N=%d)', nSubj), ...
    'FontSize', 18, 'FontWeight', 'bold');

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
saveas(fig_all, fullfile(fig_save_dir, 'GED_trials_all_subjects.png'));

%% ====================================================================
%  DETECTION RATE FIGURE
%  ====================================================================
fig_det = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Trial-Level Peak Detection Rate (mean across subjects)', ...
    'FontSize', 22, 'FontWeight', 'bold');

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
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 16, 'Box', 'off');
    ylabel('Detection Rate [%]');
    title(det_labels{di}, 'FontSize', 20, 'FontWeight', 'bold');
    grid on;
end

saveas(fig_det, fullfile(fig_save_dir, 'GED_trials_detection_rate.png'));

%% Save results
if ispc
    save_path = 'W:\Students\Arne\GCP\data\features\ged_gamma_peaks_trials.mat';
else
    save_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/ged_gamma_peaks_trials.mat';
end
save(save_path, ...
    'all_trial_powratio', ...
    'all_trial_peaks_single', 'all_trial_peaks_low', 'all_trial_peaks_high', ...
    'all_trial_mean_single', 'all_trial_median_single', ...
    'all_trial_mean_low', 'all_trial_median_low', ...
    'all_trial_mean_high', 'all_trial_median_high', ...
    'all_trial_detrate_single', 'all_trial_detrate_low', 'all_trial_detrate_high', ...
    'all_topos', 'all_topo_labels', 'all_eigenvalues', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames');

clc
fprintf('Done.\n');
