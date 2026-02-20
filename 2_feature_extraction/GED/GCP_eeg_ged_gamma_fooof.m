%% GCP Narrowband Scanning GED + Frequency Sliding for Gamma Peak Detection
%
% Two complementary methods for robust individual gamma peak identification:
%
% METHOD 1: Narrowband Scanning GED (Cohen, 2021, J Neurosci Methods)
%   - Loop over candidate frequencies (30-90 Hz in 1 Hz steps)
%   - At each frequency: narrowband filter, build stimulus & baseline
%     covariance matrices, run GED
%   - The largest eigenvalue at each frequency = multivariate SNR
%   - Detrend eigenvalue spectrum (2nd-order polynomial, analogous to
%     FOOOF aperiodic removal) to remove broadband trends
%   - Peak of detrended spectrum = individual gamma peak frequency
%
% METHOD 2: Frequency Sliding (Cohen, 2014, J Neurosci)
%   - Apply GED spatial filter (from peak frequency) to broadband gamma data
%   - Hilbert transform -> phase angle time series -> temporal derivative
%   - Instantaneous frequency distribution during stimulus period
%   - Mode of distribution = individual gamma peak frequency
%
% Both methods are visualized side-by-side for comparison.

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
scan_width = 4; % +-2 Hz bandwidth at each step

% GED parameters
lambda = 0.01; % shrinkage regularization

% Detrending parameters for eigenvalue spectrum
poly_order = 2; % 2nd-order polynomial to capture broadband trend

% Condition info
condNames  = {'c25', 'c50', 'c75', 'c100'};
trialCodes = [61, 62, 63, 64];
condLabels = {'25%', '50%', '75%', '100%'};

% Figure save directory
if ispc
    fig_save_dir = 'W:\Students\Arne\GCP\figures\eeg\ged';
else
    fig_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged';
end
if ~exist(fig_save_dir, 'dir'), mkdir(fig_save_dir); end

%% Preallocate storage
% Method 1: Narrowband scanning GED
all_eigspec       = cell(4, nSubj);  % raw eigenvalue spectrum
all_eigspec_dt    = cell(4, nSubj);  % detrended eigenvalue spectrum
all_ged_topos     = cell(4, nSubj);  % topography at peak frequency
all_ged_peakfreq  = nan(4, nSubj);   % peak freq from detrended eigspec
all_ged_peakpow   = nan(4, nSubj);   % peak eigenvalue (detrended)

% Method 2: Frequency sliding
all_fs_hist_freqs  = cell(4, nSubj);
all_fs_hist_counts = cell(4, nSubj);
all_fs_peakfreq    = nan(4, nSubj);
all_fs_median      = nan(4, nSubj);

% Channel labels
chanlocs_all = {};

%% Process each subject
for subj = 1:nSubj

    fprintf('\n========================================\n');
    fprintf('Processing Subject %d/%d: %s\n', subj, nSubj, subjects{subj});
    fprintf('========================================\n');

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

    %% Process each condition
    for cond = 1:4

        fprintf('  Condition %d/4: %s\n', cond, condLabels{cond});

        dat    = dataStructs{cond};
        trlIdx = trialIndices{cond};
        if isempty(trlIdx)
            warning('No trials for condition %d, subject %s', cond, subjects{subj});
            continue;
        end

        cfg = [];
        cfg.trials = trlIdx;
        dat = ft_selectdata(cfg, dat);
        nChans = length(dat.label);

        %% ============================================================
        %  METHOD 1: Narrowband Scanning GED
        %  ============================================================

        % --- Pass 1: Compute eigenvalue at each candidate frequency ---
        eigspec = nan(1, length(scan_freqs));

        for fi = 1:length(scan_freqs)
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
            dat_base = ft_selectdata(cfg_t, dat_nb);

            cfg_t.latency = stimulus_window;
            dat_stim = ft_selectdata(cfg_t, dat_nb);

            % Covariance matrices (trial-averaged)
            covBase = zeros(nChans);
            for trl = 1:length(dat_base.trial)
                d = double(dat_base.trial{trl});
                d = bsxfun(@minus, d, mean(d, 2));
                covBase = covBase + (d * d') / size(d, 2);
            end
            covBase = covBase / length(dat_base.trial);

            covStim = zeros(nChans);
            for trl = 1:length(dat_stim.trial)
                d = double(dat_stim.trial{trl});
                d = bsxfun(@minus, d, mean(d, 2));
                covStim = covStim + (d * d') / size(d, 2);
            end
            covStim = covStim / length(dat_stim.trial);

            % Shrinkage regularization
            covBase = (1-lambda)*covBase + lambda*mean(diag(covBase))*eye(nChans);
            covStim = (1-lambda)*covStim + lambda*mean(diag(covStim))*eye(nChans);

            % GED: largest eigenvalue = max SNR at this frequency
            [~, D] = eig(covStim, covBase);
            eigspec(fi) = max(real(diag(D)));
        end

        all_eigspec{cond, subj} = eigspec;

        % --- Detrend eigenvalue spectrum ---
        % Fit polynomial to capture broadband trend (analogous to FOOOF aperiodic fit)
        p = polyfit(scan_freqs, eigspec, poly_order);
        eigspec_trend = polyval(p, scan_freqs);
        eigspec_dt = eigspec - eigspec_trend;

        all_eigspec_dt{cond, subj} = eigspec_dt;

        % --- Find peak in detrended spectrum ---
        eigspec_dt_smooth = movmean(eigspec_dt, 5);
        [pks, locs] = findpeaks(eigspec_dt_smooth, scan_freqs, ...
            'MinPeakProminence', max(eigspec_dt_smooth) * 0.15, ...
            'MinPeakDistance', 5);

        if ~isempty(pks)
            [~, best_pk] = max(pks);
            ged_peak_freq = locs(best_pk);
        else
            [~, mi] = max(eigspec_dt_smooth);
            ged_peak_freq = scan_freqs(mi);
        end
        ged_peak_pow = eigspec_dt(scan_freqs == ged_peak_freq);

        all_ged_peakfreq(cond, subj) = ged_peak_freq;
        all_ged_peakpow(cond, subj)  = ged_peak_pow;

        % --- Pass 2: Re-run GED at peak frequency for spatial filter ---
        bpfreq_peak = [max(ged_peak_freq - scan_width/2, 1), ged_peak_freq + scan_width/2];
        cfg_filt = [];
        cfg_filt.bpfilter   = 'yes';
        cfg_filt.bpfreq     = bpfreq_peak;
        cfg_filt.bpfilttype = 'fir';
        cfg_filt.bpfiltord  = round(3 * fsample / bpfreq_peak(1));
        dat_peak = ft_preprocessing(cfg_filt, dat);

        cfg_t = [];
        cfg_t.latency = baseline_window;
        dat_base_pk = ft_selectdata(cfg_t, dat_peak);
        cfg_t.latency = stimulus_window;
        dat_stim_pk = ft_selectdata(cfg_t, dat_peak);

        covBase_pk = zeros(nChans);
        for trl = 1:length(dat_base_pk.trial)
            d = double(dat_base_pk.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            covBase_pk = covBase_pk + (d * d') / size(d, 2);
        end
        covBase_pk = covBase_pk / length(dat_base_pk.trial);

        covStim_pk = zeros(nChans);
        for trl = 1:length(dat_stim_pk.trial)
            d = double(dat_stim_pk.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            covStim_pk = covStim_pk + (d * d') / size(d, 2);
        end
        covStim_pk = covStim_pk / length(dat_stim_pk.trial);

        covBase_pk = (1-lambda)*covBase_pk + lambda*mean(diag(covBase_pk))*eye(nChans);
        covStim_pk = (1-lambda)*covStim_pk + lambda*mean(diag(covStim_pk))*eye(nChans);

        [W_pk, D_pk] = eig(covStim_pk, covBase_pk);
        [~, topIdx] = max(real(diag(D_pk)));
        topComp = W_pk(:, topIdx);

        % Sign correction
        topo_temp = covStim_pk * topComp;
        [~, mxI] = max(abs(topo_temp));
        if topo_temp(mxI) < 0, topComp = -topComp; end

        all_ged_topos{cond, subj} = covStim_pk * topComp;

        %% ============================================================
        %  METHOD 2: Frequency Sliding on GED Component
        %  ============================================================
        cfg_broad = [];
        cfg_broad.bpfilter   = 'yes';
        cfg_broad.bpfreq     = gamma_range;
        cfg_broad.bpfilttype = 'fir';
        cfg_broad.bpfiltord  = round(3 * fsample / gamma_range(1));
        dat_broad = ft_preprocessing(cfg_broad, dat);

        cfg_t = [];
        cfg_t.latency = stimulus_window;
        dat_stim_broad = ft_selectdata(cfg_t, dat_broad);

        all_inst_freq = [];
        for trl = 1:length(dat_stim_broad.trial)
            comp_ts = topComp' * double(dat_stim_broad.trial{trl});
            analytic = hilbert(comp_ts);
            phase_angles = angle(analytic);
            dphi = diff(unwrap(phase_angles));
            inst_freq = dphi * fsample / (2 * pi);
            inst_freq = medfilt1(inst_freq, 10);
            valid = inst_freq >= gamma_range(1) & inst_freq <= gamma_range(2);
            all_inst_freq = [all_inst_freq, inst_freq(valid)];
        end

        % Histogram of instantaneous frequencies
        bin_edges   = gamma_range(1):1:gamma_range(2);
        bin_centers = bin_edges(1:end-1) + 0.5;
        counts = histcounts(all_inst_freq, bin_edges);
        counts = counts / sum(counts);

        all_fs_hist_freqs{cond, subj}  = bin_centers;
        all_fs_hist_counts{cond, subj} = counts;

        [~, mode_idx] = max(counts);
        all_fs_peakfreq(cond, subj) = bin_centers(mode_idx);
        all_fs_median(cond, subj)   = median(all_inst_freq);

    end % condition loop

    %% ================================================================
    %  PER-SUBJECT FIGURE
    %  ================================================================
    close all
    fig = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('GED Gamma Peak Detection — Subject %s', subjects{subj}), ...
        'FontSize', 20, 'FontWeight', 'bold');

    % --- Row 1: Raw + detrended eigenvalue spectra (Method 1) ---
    for cond = 1:4
        subplot(4, 4, cond); hold on;
        if ~isempty(all_eigspec{cond, subj})
            plot(scan_freqs, all_eigspec{cond, subj}, '-', ...
                'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            plot(scan_freqs, movmean(all_eigspec{cond, subj}, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2);
            % Polynomial trend
            p = polyfit(scan_freqs, all_eigspec{cond, subj}, poly_order);
            plot(scan_freqs, polyval(p, scan_freqs), 'k--', 'LineWidth', 1.5);
        end
        xlabel('Freq [Hz]'); ylabel('\lambda');
        title(sprintf('%s — Raw Eigspec', condLabels{cond}), 'FontSize', 12);
        set(gca, 'FontSize', 11); xlim([30 90]); grid on; box on;
    end

    % --- Row 2: Detrended eigenvalue spectra with peak ---
    for cond = 1:4
        subplot(4, 4, 4 + cond); hold on;
        if ~isempty(all_eigspec_dt{cond, subj})
            plot(scan_freqs, all_eigspec_dt{cond, subj}, '-', ...
                'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            plot(scan_freqs, movmean(all_eigspec_dt{cond, subj}, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2.5);
            yline(0, 'k-', 'LineWidth', 0.5);
            xline(all_ged_peakfreq(cond, subj), '--', 'LineWidth', 2, ...
                'Color', colors(cond,:));
            text(all_ged_peakfreq(cond, subj) + 1, ...
                max(all_eigspec_dt{cond, subj}) * 0.9, ...
                sprintf('%.0f Hz', all_ged_peakfreq(cond, subj)), ...
                'FontSize', 12, 'Color', colors(cond,:), 'FontWeight', 'bold');
        end
        xlabel('Freq [Hz]'); ylabel('\Delta\lambda');
        title(sprintf('%s — Detrended', condLabels{cond}), 'FontSize', 12);
        set(gca, 'FontSize', 11); xlim([30 90]); grid on; box on;
    end

    % --- Row 3: Frequency sliding distributions (Method 2) ---
    for cond = 1:4
        subplot(4, 4, 8 + cond); hold on;
        if ~isempty(all_fs_hist_counts{cond, subj})
            bar(all_fs_hist_freqs{cond, subj}, all_fs_hist_counts{cond, subj}, 1, ...
                'FaceColor', colors(cond,:), 'FaceAlpha', 0.6, 'EdgeColor', 'none');
            xline(all_fs_peakfreq(cond, subj), '--', 'LineWidth', 2, ...
                'Color', colors(cond,:));
            text(all_fs_peakfreq(cond, subj) + 1, ...
                max(all_fs_hist_counts{cond, subj}) * 0.9, ...
                sprintf('%.0f Hz', all_fs_peakfreq(cond, subj)), ...
                'FontSize', 12, 'Color', colors(cond,:), 'FontWeight', 'bold');
        end
        xlabel('Freq [Hz]'); ylabel('P(f)');
        title(sprintf('%s — Freq Sliding', condLabels{cond}), 'FontSize', 12);
        set(gca, 'FontSize', 11); xlim([30 90]); grid on; box on;
    end

    % --- Row 4: Topographies at GED peak frequency ---
    cfg_topo = [];
    cfg_topo.layout    = headmodel.layANThead;
    cfg_topo.comment   = 'no';
    cfg_topo.marker    = 'off';
    cfg_topo.style     = 'straight';
    cfg_topo.gridscale = 300;
    cfg_topo.zlim      = 'maxabs';
    cfg_topo.colormap  = '*RdBu';
    cfg_topo.figure    = 'gcf';

    for cond = 1:4
        subplot(4, 4, 12 + cond);
        if ~isempty(all_ged_topos{cond, subj})
            topo_data = [];
            topo_data.label  = chanlocs_all;
            topo_data.avg    = all_ged_topos{cond, subj};
            topo_data.dimord = 'chan';
            try
                ft_topoplotER(cfg_topo, topo_data);
                cb = colorbar;
                cb.FontSize = 9;
            catch
                imagesc(topo_data.avg); colorbar;
            end
            title(sprintf('%s — Topo @ %.0f Hz', condLabels{cond}, ...
                all_ged_peakfreq(cond, subj)), 'FontSize', 12);
        end
    end

    saveas(fig, fullfile(fig_save_dir, sprintf('GED_scan_subj%s.png', subjects{subj})));

end % subject loop

%% ====================================================================
%  GRAND AVERAGE FIGURE 1: Detrended eigenvalue spectra + peak scatter
%  ====================================================================
fprintf('\nCreating grand average figures...\n');
close all

fig_ga1 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Narrowband Scanning GED — Grand Average (N=%d)', nSubj), ...
    'FontSize', 20, 'FontWeight', 'bold');

% --- Left: Raw eigenvalue spectra ---
subplot(2, 2, 1); hold on;
for cond = 1:4
    eigspec_mat = nan(nSubj, length(scan_freqs));
    for s = 1:nSubj
        if ~isempty(all_eigspec{cond, s})
            eigspec_mat(s,:) = all_eigspec{cond, s};
        end
    end
    mu  = nanmean(eigspec_mat, 1);
    sem = nanstd(eigspec_mat, [], 1) / sqrt(sum(~isnan(eigspec_mat(:,1))));
    mu_s  = movmean(mu, 5);
    sem_s = movmean(sem, 5);
    faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
    patch([scan_freqs, fliplr(scan_freqs)], ...
          [mu_s - sem_s, fliplr(mu_s + sem_s)], ...
          colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(scan_freqs, mu_s, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);
end
xlabel('Frequency [Hz]'); ylabel('Eigenvalue (\lambda)');
title('Raw Eigenvalue Spectra (mean +/- SEM)', 'FontSize', 14);
legend(condLabels, 'Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 13); xlim([30 90]); grid on; box on;

% --- Right: Detrended eigenvalue spectra ---
subplot(2, 2, 2); hold on;
for cond = 1:4
    eigdt_mat = nan(nSubj, length(scan_freqs));
    for s = 1:nSubj
        if ~isempty(all_eigspec_dt{cond, s})
            eigdt_mat(s,:) = all_eigspec_dt{cond, s};
        end
    end
    mu  = nanmean(eigdt_mat, 1);
    sem = nanstd(eigdt_mat, [], 1) / sqrt(sum(~isnan(eigdt_mat(:,1))));
    mu_s  = movmean(mu, 5);
    sem_s = movmean(sem, 5);
    faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
    patch([scan_freqs, fliplr(scan_freqs)], ...
          [mu_s - sem_s, fliplr(mu_s + sem_s)], ...
          colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(scan_freqs, mu_s, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);
end
yline(0, 'k-', 'LineWidth', 0.5);
xlabel('Frequency [Hz]'); ylabel('\Delta\lambda (detrended)');
title('Detrended Eigenvalue Spectra (mean +/- SEM)', 'FontSize', 14);
legend(condLabels, 'Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 13); xlim([30 90]); grid on; box on;

% --- Bottom left: Peak frequency scatter (GED) ---
subplot(2, 2, 3); hold on;
for cond = 1:4
    pf = all_ged_peakfreq(cond, :);
    pf = pf(~isnan(pf));
    scatter(ones(size(pf)) * cond + 0.15*(rand(size(pf))-0.5), pf, ...
        60, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.6);
    mu_pf  = nanmean(all_ged_peakfreq(cond,:));
    sem_pf = nanstd(all_ged_peakfreq(cond,:)) / sqrt(sum(~isnan(all_ged_peakfreq(cond,:))));
    errorbar(cond, mu_pf, sem_pf, 'k', 'LineWidth', 2, 'CapSize', 10);
    plot(cond, mu_pf, 'kd', 'MarkerSize', 12, ...
        'MarkerFaceColor', colors(cond,:), 'LineWidth', 1.5);
end
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 13);
ylabel('Peak Gamma Frequency [Hz]');
title('Individual Peak Frequencies (Scanning GED)', 'FontSize', 14);
ylim([30 90]); grid on; box on;

% --- Bottom right: Peak frequency scatter (Freq Sliding) ---
subplot(2, 2, 4); hold on;
for cond = 1:4
    pf = all_fs_peakfreq(cond, :);
    pf = pf(~isnan(pf));
    scatter(ones(size(pf)) * cond + 0.15*(rand(size(pf))-0.5), pf, ...
        60, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.6);
    mu_pf  = nanmean(all_fs_peakfreq(cond,:));
    sem_pf = nanstd(all_fs_peakfreq(cond,:)) / sqrt(sum(~isnan(all_fs_peakfreq(cond,:))));
    errorbar(cond, mu_pf, sem_pf, 'k', 'LineWidth', 2, 'CapSize', 10);
    plot(cond, mu_pf, 'kd', 'MarkerSize', 12, ...
        'MarkerFaceColor', colors(cond,:), 'LineWidth', 1.5);
end
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 13);
ylabel('Peak Gamma Frequency [Hz]');
title('Individual Peak Frequencies (Freq Sliding)', 'FontSize', 14);
ylim([30 90]); grid on; box on;

saveas(fig_ga1, fullfile(fig_save_dir, 'GED_scan_grand_average.png'));

%% ====================================================================
%  GRAND AVERAGE FIGURE 2: Method comparison
%  ====================================================================
fig_ga2 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Method Comparison: Scanning GED vs Frequency Sliding (N=%d)', nSubj), ...
    'FontSize', 20, 'FontWeight', 'bold');

% --- Left: Correlation between methods ---
subplot(1, 2, 1); hold on;
all_ged_flat = all_ged_peakfreq(:);
all_fs_flat  = all_fs_peakfreq(:);
valid = ~isnan(all_ged_flat) & ~isnan(all_fs_flat);
for cond = 1:4
    scatter(all_ged_peakfreq(cond,:), all_fs_peakfreq(cond,:), ...
        80, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.7);
end
plot([30 90], [30 90], 'k--', 'LineWidth', 1.5);
if sum(valid) > 2
    [r, p_val] = corr(all_ged_flat(valid), all_fs_flat(valid), 'type', 'Spearman');
    text(32, 85, sprintf('r_s = %.2f, p = %.3f', r, p_val), 'FontSize', 14);
end
xlabel('GED Peak Freq [Hz]'); ylabel('Freq Sliding Peak Freq [Hz]');
title('Agreement Between Methods', 'FontSize', 15);
set(gca, 'FontSize', 14); xlim([30 90]); ylim([30 90]);
axis square; grid on; box on;
legend(condLabels, 'Location', 'southeast', 'FontSize', 12);

% --- Right: Frequency sliding grand average distributions ---
subplot(1, 2, 2); hold on;
bin_centers_common = gamma_range(1)+0.5 : 1 : gamma_range(2)-0.5;
for cond = 1:4
    counts_mat = nan(nSubj, length(bin_centers_common));
    for s = 1:nSubj
        if ~isempty(all_fs_hist_counts{cond, s})
            counts_mat(s,:) = interp1(all_fs_hist_freqs{cond, s}, ...
                all_fs_hist_counts{cond, s}, bin_centers_common, 'linear', 0);
        end
    end
    mu = nanmean(counts_mat, 1);
    sem = nanstd(counts_mat, [], 1) / sqrt(sum(~isnan(counts_mat(:,1))));
    mu_s = movmean(mu, 3);
    sem_s = movmean(sem, 3);
    faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
    patch([bin_centers_common, fliplr(bin_centers_common)], ...
          [mu_s - sem_s, fliplr(mu_s + sem_s)], ...
          colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bin_centers_common, mu_s, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);
end
xlabel('Frequency [Hz]'); ylabel('P(f)');
title('Grand Average Freq Sliding Distributions', 'FontSize', 15);
legend(condLabels, 'Location', 'northeast', 'FontSize', 12);
set(gca, 'FontSize', 14); xlim([30 90]); grid on; box on;

saveas(fig_ga2, fullfile(fig_save_dir, 'GED_method_comparison.png'));

%% ====================================================================
%  GRAND AVERAGE FIGURE 3: All subjects subplot
%  ====================================================================
nRows = ceil(nSubj / 5);
fig_all = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Detrended Eigenvalue Spectra — All Subjects (N=%d)', nSubj), ...
    'FontSize', 18, 'FontWeight', 'bold');

for s = 1:nSubj
    subplot(nRows, 5, s); hold on;
    for cond = 1:4
        if ~isempty(all_eigspec_dt{cond, s})
            plot(scan_freqs, movmean(all_eigspec_dt{cond, s}, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2);
            xline(all_ged_peakfreq(cond, s), '--', 'Color', colors(cond,:), ...
                'LineWidth', 1.2);
        end
    end
    yline(0, 'k-', 'LineWidth', 0.5);
    xlabel('Freq [Hz]'); ylabel('\Delta\lambda');
    title(sprintf('Subj %s', subjects{s}), 'FontSize', 11);
    set(gca, 'FontSize', 10); xlim([30 90]); grid on; box on;
    if s == 1
        legend(condLabels, 'FontSize', 8, 'Location', 'best');
    end
end
saveas(fig_all, fullfile(fig_save_dir, 'GED_scan_all_subjects.png'));

%% Save results
if ispc
    save_path = 'W:\Students\Arne\GCP\data\features\ged_gamma_peaks.mat';
else
    save_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/ged_gamma_peaks.mat';
end
save(save_path, ...
    'all_eigspec', 'all_eigspec_dt', 'all_ged_topos', ...
    'all_ged_peakfreq', 'all_ged_peakpow', ...
    'all_fs_hist_freqs', 'all_fs_hist_counts', 'all_fs_peakfreq', 'all_fs_median', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames');

fprintf('\n========================================\n');
fprintf('GED Gamma Peak Detection Complete!\n');
fprintf('Results saved to:\n  %s\n', save_path);
fprintf('Figures saved to:\n  %s\n', fig_save_dir);
fprintf('========================================\n');
