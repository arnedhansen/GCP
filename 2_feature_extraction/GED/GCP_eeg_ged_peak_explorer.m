%% GCP GED Peak Explorer — Assumption-Free Peak Counting
%
% Runs the same pooled GED + narrowband scanning as the main trial-level
% script, but strips away all dual-peak models and boundary logic.
% Instead, it simply detects ALL peaks in each trial's detrended
% power-ratio spectrum (using findpeaks) and records:
%   - How many peaks per trial
%   - At which frequencies they occur
%
% Visualisation:
%   Per subject  — (1) heatmap of detrended spectra with peak markers,
%                  (2) histogram of all detected peak frequencies,
%                  (3) peak-count distribution across trials
%   Grand avg    — (1) pooled peak-frequency histogram (all subjects),
%                  (2) peak-count distribution (all subjects),
%                  (3) 2-D density heatmap (frequency x subject),
%                  (4) per-subject overview spectra

clear; close all; clc

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');

nSubj = length(subjects);

baseline_window = [-1.5, -0.25];
stimulus_window = [0.3, 2.0];

gamma_range = [30, 90];

scan_freqs = 30:1:90;
nFreqs     = length(scan_freqs);
scan_width = 4;

lambda     = 0.01;
poly_order = 2;

condNames  = {'c25', 'c50', 'c75', 'c100'};
trialCodes = [61, 62, 63, 64];
condLabels = {'25%', '50%', '75%', '100%'};

if ispc
    fig_save_dir = 'W:\Students\Arne\GCP\figures\eeg\ged\peak_explorer';
else
    fig_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged/peak_explorer';
end
if ~exist(fig_save_dir, 'dir'), mkdir(fig_save_dir); end

%% Preallocate
all_trial_powratio = cell(4, nSubj);
all_peak_freqs     = cell(4, nSubj);   % each cell = cell array of peak-freq vectors per trial
all_peak_counts    = cell(4, nSubj);   % each cell = vector of #peaks per trial

all_topos       = cell(1, nSubj);
all_topo_labels = cell(1, nSubj);
all_eigenvalues = nan(1, nSubj);

%% ====================================================================
%  PROCESS EACH SUBJECT
%  ====================================================================
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

    nChans = length(dataEEG_c25.label);

    occ_mask = cellfun(@(l) ~isempty(regexp(l, '[OI]', 'once')), dataEEG_c25.label);
    occ_idx  = find(occ_mask);
    nOcc     = length(occ_idx);

    %% Phase 1: Pooled GED — common spatial filter (occipital)
    clc
    fprintf('Subject %s (%d/%d) — Phase 1: Occipital GED (%d occ / %d ch)\n', ...
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

    % Extract occipital submatrices for GED
    covStim_occ = covStim_full(occ_idx, occ_idx);
    covBase_occ = covBase_full(occ_idx, occ_idx);

    % Shrinkage regularization (occipital)
    covStim_occ = (1-lambda)*covStim_occ + lambda*mean(diag(covStim_occ))*eye(nOcc);
    covBase_occ = (1-lambda)*covBase_occ + lambda*mean(diag(covBase_occ))*eye(nOcc);

    % GED on occipital covariance
    [W_occ, D_occ] = eig(covStim_occ, covBase_occ);
    [evals_sorted, sortIdx] = sort(real(diag(D_occ)), 'descend');
    W_occ = W_occ(:, sortIdx);
    w_occ = W_occ(:, 1);
    all_eigenvalues(subj) = evals_sorted(1);

    % Zero-pad to full channel space
    topComp = zeros(nChans, 1);
    topComp(occ_idx) = w_occ;

    % Full-head forward model for topoplot
    covStim_full_reg = (1-lambda)*covStim_full + lambda*mean(diag(covStim_full))*eye(nChans);
    topo_temp = covStim_full_reg * topComp;
    [~, mxI] = max(abs(topo_temp));
    if topo_temp(mxI) < 0, topComp = -topComp; end

    all_topos{subj}       = covStim_full_reg * topComp;
    all_topo_labels{subj} = dataEEG_c25.label;

    %% Phase 2: Narrowband scanning + assumption-free peak detection
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

        %% Assumption-free peak detection: find ALL peaks per trial
        trial_peak_freqs  = cell(nTrl, 1);
        trial_peak_counts = nan(nTrl, 1);

        for trl = 1:nTrl
            pr = powratio_trials(trl, :);
            if all(isnan(pr))
                trial_peak_counts(trl) = 0;
                trial_peak_freqs{trl}  = [];
                continue;
            end

            p = polyfit(scan_freqs, pr, poly_order);
            pr_dt = pr - polyval(p, scan_freqs);
            pr_dt_smooth = movmean(pr_dt, 5);

            [pks, locs] = findpeaks(pr_dt_smooth, scan_freqs, ...
                'MinPeakDistance', 5);

            pos_mask = pks > 0;
            locs = locs(pos_mask);
            pks  = pks(pos_mask);

            trial_peak_freqs{trl}  = locs(:)';
            trial_peak_counts(trl) = length(locs);
        end

        all_peak_freqs{cond, subj}  = trial_peak_freqs;
        all_peak_counts{cond, subj} = trial_peak_counts;

    end % condition loop

    %% ================================================================
    %  PER-SUBJECT FIGURE (3 rows x 4 columns + summary row)
    %  ================================================================
    close all
    fig = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('Peak Explorer: Subject %s  (\\lambda_1=%.2f)', ...
        subjects{subj}, all_eigenvalues(subj)), ...
        'FontSize', 16, 'FontWeight', 'bold');

    cmap_div = interp1([0 0.5 1], ...
        [0.17 0.27 0.53; 0.97 0.97 0.97; 0.70 0.09 0.17], linspace(0,1,256));

    % Compute detrended matrices + shared colour limits
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
    all_vals = cell2mat(cellfun(@(x) x(:), ...
        pr_dt_mats(~cellfun(@isempty, pr_dt_mats)), 'UniformOutput', false)');
    if ~isempty(all_vals)
        cl_shared = prctile(abs(all_vals(~isnan(all_vals))), 98);
    else
        cl_shared = 1;
    end

    % --- Row 1: Heatmaps with peak markers ---
    for cond = 1:4
        subplot(4, 4, cond); hold on;
        if ~isempty(pr_dt_mats{cond})
            imagesc(scan_freqs, 1:size(pr_dt_mats{cond},1), pr_dt_mats{cond});
            colormap(gca, cmap_div);
            caxis([-cl_shared cl_shared]);
            set(gca, 'YDir', 'normal');

            tpf = all_peak_freqs{cond, subj};
            for trl = 1:length(tpf)
                if ~isempty(tpf{trl})
                    plot(tpf{trl}, repmat(trl, size(tpf{trl})), ...
                        'k.', 'MarkerSize', 4);
                end
            end
        end
        xlabel('Freq [Hz]'); ylabel('Trial');
        nTrl_c = size(pr_dt_mats{cond}, 1);
        title(sprintf('%s  (n=%d)', condLabels{cond}, nTrl_c), 'FontSize', 11);
        set(gca, 'FontSize', 10); xlim([30 90]); box on;
    end

    % --- Row 2: Mean spectrum + individual trial spectra (faint) ---
    % Compute robust shared y-limits across conditions from smoothed trial data
    all_smooth_vals = [];
    for cond = 1:4
        if ~isempty(pr_dt_mats{cond})
            nTrl_tmp = size(pr_dt_mats{cond}, 1);
            for trl = 1:nTrl_tmp
                all_smooth_vals = [all_smooth_vals, movmean(pr_dt_mats{cond}(trl,:), 5)];
            end
        end
    end
    if ~isempty(all_smooth_vals)
        yl_lo = prctile(all_smooth_vals(~isnan(all_smooth_vals)), 1);
        yl_hi = prctile(all_smooth_vals(~isnan(all_smooth_vals)), 99);
        yl_pad = (yl_hi - yl_lo) * 0.15;
        spec_ylim = [yl_lo - yl_pad, yl_hi + yl_pad];
    else
        spec_ylim = [-1 1];
    end

    for cond = 1:4
        subplot(4, 4, 4 + cond); hold on;
        if ~isempty(pr_dt_mats{cond})
            nTrl = size(pr_dt_mats{cond}, 1);
            for trl = 1:nTrl
                plot(scan_freqs, movmean(pr_dt_mats{cond}(trl,:), 5), ...
                    '-', 'Color', [colors(cond,:) 0.08], 'LineWidth', 0.5);
            end
            mu_dt = nanmean(pr_dt_mats{cond}, 1);
            plot(scan_freqs, movmean(mu_dt, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2.5);
            yline(0, 'k-', 'LineWidth', 0.5);
        end
        xlabel('Freq [Hz]'); ylabel('\Delta PR');
        title(sprintf('%s Spectra', condLabels{cond}), 'FontSize', 11);
        set(gca, 'FontSize', 10); xlim([30 90]); ylim(spec_ylim); grid on; box on;
    end

    % --- Row 3: Peak-frequency histograms per condition ---
    edges = 30:2:90;
    centers = edges(1:end-1) + diff(edges)/2;
    for cond = 1:4
        subplot(4, 4, 8 + cond); hold on;
        tpf = all_peak_freqs{cond, subj};
        all_pf = [];
        if ~isempty(tpf)
            all_pf = [tpf{:}];
        end
        if ~isempty(all_pf)
            histogram(all_pf, edges, 'FaceColor', colors(cond,:), ...
                'EdgeColor', 'w', 'FaceAlpha', 0.8);
        end
        xlabel('Peak Freq [Hz]'); ylabel('Count');
        nPeaks = length(all_pf);
        nTrl_c = length(tpf);
        title(sprintf('%s: %d peaks / %d trials', condLabels{cond}, nPeaks, nTrl_c), ...
            'FontSize', 10);
        set(gca, 'FontSize', 10); xlim([30 90]); box on;
    end

    % --- Row 4: Peak-count distribution (stacked) + combined histogram ---
    subplot(4, 4, 13); hold on;
    max_pk = 0;
    for cond = 1:4
        pc = all_peak_counts{cond, subj};
        if ~isempty(pc), max_pk = max(max_pk, max(pc)); end
    end
    pk_edges = -0.5:1:(max_pk + 1.5);
    pk_centers = 0:max_pk+1;
    pk_count_mat = zeros(4, length(pk_centers));
    for cond = 1:4
        pc = all_peak_counts{cond, subj};
        if ~isempty(pc)
            pk_count_mat(cond,:) = histcounts(pc, pk_edges);
        end
    end
    bh_pk = bar(pk_centers, pk_count_mat', 'stacked', 'EdgeColor', 'w', 'BarWidth', 0.85);
    for cond = 1:4
        bh_pk(cond).FaceColor = colors(cond,:);
    end
    xlabel('# Peaks per Trial'); ylabel('Trial Count');
    title('Peak Count Distribution', 'FontSize', 11);
    legend(condLabels, 'FontSize', 8, 'Location', 'best');
    set(gca, 'FontSize', 10); box on;

    subplot(4, 4, [14 15 16]); hold on;
    hist_mat = zeros(4, length(edges)-1);
    for cond = 1:4
        tpf = all_peak_freqs{cond, subj};
        if ~isempty(tpf)
            all_pf = [tpf{:}];
            hist_mat(cond,:) = histcounts(all_pf, edges);
        end
    end
    bh = bar(centers, hist_mat', 'stacked', 'EdgeColor', 'none', 'BarWidth', 1);
    for cond = 1:4
        bh(cond).FaceColor = colors(cond,:);
    end
    xlabel('Peak Frequency [Hz]'); ylabel('Count (all peaks)');
    title('All Detected Peaks — Stacked by Condition', 'FontSize', 12);
    legend(condLabels, 'FontSize', 10, 'Location', 'best');
    set(gca, 'FontSize', 11); xlim([30 90]); grid on; box on;

    saveas(fig, fullfile(fig_save_dir, ...
        sprintf('GCP_peak_explorer_subj%s.png', subjects{subj})));

end % subject loop

%% ====================================================================
%  GRAND AVERAGE FIGURES
%  ====================================================================
fprintf('\nCreating grand-average figures...\n');
close all

%% --- Figure 1: Pooled peak-frequency histogram (all subjects, all trials) ---
fig1 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('All Detected Peaks — Pooled Across %d Subjects', nSubj), ...
    'FontSize', 18, 'FontWeight', 'bold');

edges = 30:2:90;
centers = edges(1:end-1) + diff(edges)/2;

% (a) Stacked by condition
subplot(2, 2, 1); hold on;
hist_mat = zeros(4, length(edges)-1);
total_peaks = zeros(1, 4);
total_trials = zeros(1, 4);
for cond = 1:4
    for s = 1:nSubj
        tpf = all_peak_freqs{cond, s};
        if ~isempty(tpf)
            all_pf = [tpf{:}];
            hist_mat(cond,:) = hist_mat(cond,:) + histcounts(all_pf, edges);
            total_peaks(cond) = total_peaks(cond) + length(all_pf);
            total_trials(cond) = total_trials(cond) + length(tpf);
        end
    end
end
bh = bar(centers, hist_mat', 'stacked', 'EdgeColor', 'none', 'BarWidth', 1);
for cond = 1:4
    bh(cond).FaceColor = colors(cond,:);
end
xlabel('Peak Frequency [Hz]'); ylabel('Count');
title('Stacked by Condition', 'FontSize', 14);
legend(condLabels, 'FontSize', 12, 'Location', 'best');
set(gca, 'FontSize', 12); xlim([30 90]); grid on; box on;

% (b) Collapsed across conditions
subplot(2, 2, 2); hold on;
all_peaks_pooled = [];
for cond = 1:4
    for s = 1:nSubj
        tpf = all_peak_freqs{cond, s};
        if ~isempty(tpf)
            all_pf = [tpf{:}];
            all_peaks_pooled = [all_peaks_pooled, all_pf];
        end
    end
end
histogram(all_peaks_pooled, edges, 'FaceColor', [0.3 0.3 0.7], ...
    'EdgeColor', 'w', 'FaceAlpha', 0.8);
xlabel('Peak Frequency [Hz]'); ylabel('Count');
title(sprintf('All Conditions Collapsed  (N=%d peaks)', length(all_peaks_pooled)), ...
    'FontSize', 14);
set(gca, 'FontSize', 12); xlim([30 90]); grid on; box on;

% (c) Peak count distribution
subplot(2, 2, 3); hold on;
all_counts_pooled = [];
for cond = 1:4
    for s = 1:nSubj
        pc = all_peak_counts{cond, s};
        if ~isempty(pc)
            all_counts_pooled = [all_counts_pooled; pc(:)];
        end
    end
end
max_pk = max(all_counts_pooled);
pk_edges = -0.5:1:(max_pk + 1.5);
histogram(all_counts_pooled, pk_edges, 'FaceColor', [0.4 0.4 0.4], ...
    'EdgeColor', 'w', 'FaceAlpha', 0.8);
xlabel('# Peaks per Trial'); ylabel('Trial Count');
mu_cnt = mean(all_counts_pooled);
md_cnt = median(all_counts_pooled);
title(sprintf('Peak Count per Trial  (mean=%.1f, median=%.0f)', mu_cnt, md_cnt), ...
    'FontSize', 14);
set(gca, 'FontSize', 12); box on;

% (d) Summary table as text
subplot(2, 2, 4); axis off;
txt = sprintf('Summary\n-------\n');
for cond = 1:4
    txt = [txt, sprintf('%s:  %d peaks across %d trials  (%.1f peaks/trial)\n', ...
        condLabels{cond}, total_peaks(cond), total_trials(cond), ...
        total_peaks(cond)/max(total_trials(cond),1))];
end
txt = [txt, sprintf('\nTotal:  %d peaks across %d trials\n', ...
    sum(total_peaks), sum(total_trials))];
txt = [txt, sprintf('Overall: %.1f peaks/trial\n', sum(total_peaks)/max(sum(total_trials),1))];
text(0.1, 0.5, txt, 'FontSize', 14, 'FontName', 'FixedWidth', ...
    'VerticalAlignment', 'middle');

saveas(fig1, fullfile(fig_save_dir, 'GCP_peak_explorer_grand_histogram.png'));

%% --- Figure 2: 2-D density — frequency x subject ---
fig2 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Peak Frequency Density per Subject (all conditions pooled)', ...
    'FontSize', 16, 'FontWeight', 'bold');

freq_bins = 30:1:90;
density_mat = nan(nSubj, length(freq_bins)-1);
for s = 1:nSubj
    pf_subj = [];
    for cond = 1:4
        tpf = all_peak_freqs{cond, s};
        if ~isempty(tpf)
            pf_subj = [pf_subj, [tpf{:}]];
        end
    end
    if ~isempty(pf_subj)
        density_mat(s,:) = histcounts(pf_subj, freq_bins, 'Normalization', 'probability');
    end
end

subplot(1, 1, 1);
imagesc(freq_bins(1:end-1) + 0.5, 1:nSubj, density_mat);
colormap(parula);
cb = colorbar; cb.Label.String = 'P(peak at freq)'; cb.FontSize = 12;
xlabel('Frequency [Hz]'); ylabel('Subject');
set(gca, 'YTick', 1:nSubj, 'YTickLabel', subjects, 'FontSize', 11);
title('Normalized Peak Density', 'FontSize', 14);
set(gca, 'YDir', 'normal'); box on;

saveas(fig2, fullfile(fig_save_dir, 'GCP_peak_explorer_density_heatmap.png'));

%% --- Figure 3: Per-subject overview (mean spectra + detected peaks) ---
nRows = ceil(nSubj / 5);
fig3 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Mean Detrended Spectra + Detected Peaks  (N=%d)', nSubj), ...
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
        end
    end
    yline(0, 'k-', 'LineWidth', 0.5);

    % Mark most common peak frequencies (mode per 2 Hz bin)
    pf_subj = [];
    for cond = 1:4
        tpf = all_peak_freqs{cond, s};
        if ~isempty(tpf)
            pf_subj = [pf_subj, [tpf{:}]];
        end
    end
    if ~isempty(pf_subj)
        yl = ylim;
        [~, mode_freq] = max(histcounts(pf_subj, 30:2:90));
        mode_center = 30 + (mode_freq - 1) * 2 + 1;
        xline(mode_center, 'k--', 'LineWidth', 1.5);
        text(mode_center + 1, yl(2) * 0.85, sprintf('%d', mode_center), ...
            'FontSize', 8, 'FontWeight', 'bold');
    end

    xlabel('Hz'); ylabel('\Delta PR');
    title(sprintf('S%s', subjects{s}), 'FontSize', 10);
    set(gca, 'FontSize', 9); xlim([30 90]); grid on; box on;
    if s == 1
        legend(condLabels, 'FontSize', 7, 'Location', 'best');
    end
end

saveas(fig3, fullfile(fig_save_dir, 'GCP_peak_explorer_all_subjects.png'));

%% --- Figure 4: Per-condition peak frequency distributions (violin-style) ---
fig4 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Peak Frequency Distributions by Condition (per subject)', ...
    'FontSize', 18, 'FontWeight', 'bold');

for cond = 1:4
    subplot(2, 2, cond); hold on;

    for s = 1:nSubj
        tpf = all_peak_freqs{cond, s};
        if isempty(tpf), continue; end
        pf = [tpf{:}];
        if length(pf) >= 3
            [f_dens, xi] = ksdensity(pf, 'Support', [29.5, 90.5]);
            f_dens = f_dens / max(f_dens) * 0.4;
            y_pos = s;
            patch(xi, y_pos + f_dens, colors(cond,:), ...
                'FaceAlpha', 0.4, 'EdgeColor', colors(cond,:), 'LineWidth', 0.8);
            plot([median(pf) median(pf)], [y_pos y_pos + 0.35], ...
                'k-', 'LineWidth', 2);
        end
    end

    xlabel('Peak Frequency [Hz]'); ylabel('Subject');
    set(gca, 'YTick', 1:nSubj, 'YTickLabel', subjects, 'FontSize', 10);
    xlim([30 90]); ylim([0 nSubj + 1]);
    title(sprintf('%s  (%d total peaks)', condLabels{cond}, total_peaks(cond)), ...
        'FontSize', 14);
    grid on; box on;
end

saveas(fig4, fullfile(fig_save_dir, 'GCP_peak_explorer_violin_by_cond.png'));

%% Save data
if ispc
    save_path = 'W:\Students\Arne\GCP\data\features\GCP_eeg_GED_peak_explorer.mat';
else
    save_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_eeg_GED_peak_explorer.mat';
end
save(save_path, ...
    'all_trial_powratio', 'all_peak_freqs', 'all_peak_counts', ...
    'all_topos', 'all_topo_labels', 'all_eigenvalues', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames');

clc
fprintf('Peak Explorer — Done.\n');
fprintf('Figures saved to: %s\n', fig_save_dir);
