%% GCP Narrowband Scanning GED + Frequency Sliding for Gamma Peak Detection
%
% Two complementary methods for robust individual gamma peak identification:
%
% METHOD 1: Narrowband Scanning GED (Cohen, 2021, J Neurosci Methods)
%   - Loop over candidate frequencies (30-90 Hz in 1 Hz steps)
%   - At each frequency: narrowband filter, build stimulus & baseline
%     covariance matrices, run GED
%   - The largest eigenvalue at each frequency = multivariate SNR
%   - Peak of eigenvalue-vs-frequency curve = individual gamma peak freq
%   - Leverages spatial structure across all channels for denoising
%
% METHOD 2: Frequency Sliding (Cohen, 2014, J Neurosci)
%   - Apply GED spatial filter (from best frequency) to broadband gamma data
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

% Condition info
condNames  = {'c25', 'c50', 'c75', 'c100'};
trialCodes = [61, 62, 63, 64];
condLabels = {'25%', '50%', '75%', '100%'};

% Figure save directory
if ispc
    fig_save_dir = 'W:\Students\Arne\GCP\figures\eeg\ged';
    headmodel_path = 'W:\Students\Arne\MA\headmodel\ant128lay.mat';
else
    fig_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged';
    headmodel_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/ant128lay.mat';
end
if ~exist(fig_save_dir, 'dir'), mkdir(fig_save_dir); end
load(headmodel_path);

%% Preallocate storage
% Method 1: Narrowband scanning GED
all_eigspec      = cell(4, nSubj);   % eigenvalue spectrum per freq [cond x subj]
all_ged_topos    = cell(4, nSubj);   % topography at peak frequency
all_ged_peakfreq = nan(4, nSubj);    % peak freq from eigenvalue spectrum
all_ged_peakpow  = nan(4, nSubj);    % peak eigenvalue

% Method 2: Frequency sliding
all_fs_hist_freqs  = cell(4, nSubj); % histogram bin centers
all_fs_hist_counts = cell(4, nSubj); % histogram counts
all_fs_peakfreq    = nan(4, nSubj);  % peak freq from freq sliding
all_fs_median      = nan(4, nSubj);  % median instantaneous frequency
all_fs_timeseries  = cell(4, nSubj); % raw freq sliding time series (one trial)

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
        eigspec = nan(1, length(scan_freqs));
        best_eigval = -Inf;
        best_topo   = [];

        for fi = 1:length(scan_freqs)
            cf = scan_freqs(fi);
            bpfreq = [max(cf - scan_width/2, 1), cf + scan_width/2];

            % Narrowband filter
            cfg_filt = [];
            cfg_filt.bpfilter   = 'yes';
            cfg_filt.bpfreq     = bpfreq;
            cfg_filt.bpfilttype = 'fir';
            cfg_filt.bpfiltord  = round(3 * fsample / bpfreq(1));
            dat_nb = ft_preprocessing(cfg_filt, dat);

            % Extract baseline and stimulus windows
            cfg_t = [];
            cfg_t.latency = baseline_window;
            dat_base = ft_selectdata(cfg_t, dat_nb);

            cfg_t.latency = stimulus_window;
            dat_stim = ft_selectdata(cfg_t, dat_nb);

            % Compute covariance matrices (average across trials)
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

            % GED
            [W, D] = eig(covStim, covBase);
            evals = real(diag(D));
            [max_eval, maxIdx] = max(evals);
            eigspec(fi) = max_eval;

            % Track the best frequency for topography and freq sliding
            if max_eval > best_eigval
                best_eigval  = max_eval;
                best_W       = W;
                best_covStim = covStim;
                best_maxIdx  = maxIdx;
            end
        end

        all_eigspec{cond, subj} = eigspec;

        % Find peak of eigenvalue spectrum (smoothed)
        eigspec_smooth = movmean(eigspec, 5);
        [~, peak_idx] = max(eigspec_smooth);
        ged_peak_freq = scan_freqs(peak_idx);
        ged_peak_pow  = eigspec(peak_idx);

        all_ged_peakfreq(cond, subj) = ged_peak_freq;
        all_ged_peakpow(cond, subj)  = ged_peak_pow;

        % Topography at best frequency
        topComp = best_W(:, best_maxIdx);
        topo_temp = best_covStim * topComp;
        [~, mxI] = max(abs(topo_temp));
        if topo_temp(mxI) < 0, topComp = -topComp; end
        all_ged_topos{cond, subj} = best_covStim * topComp;

        %% ============================================================
        %  METHOD 2: Frequency Sliding on GED Component
        %  ============================================================
        %  Use the spatial filter from the best narrowband GED,
        %  but apply it to broadband gamma-filtered data.

        % Broadband gamma filter (plateau-shaped via wider bandwidth)
        cfg_broad = [];
        cfg_broad.bpfilter   = 'yes';
        cfg_broad.bpfreq     = gamma_range;
        cfg_broad.bpfilttype = 'fir';
        cfg_broad.bpfiltord  = round(3 * fsample / gamma_range(1));
        dat_broad = ft_preprocessing(cfg_broad, dat);

        % Project broadband data onto GED spatial filter
        cfg_t = [];
        cfg_t.latency = stimulus_window;
        dat_stim_broad = ft_selectdata(cfg_t, dat_broad);

        % Compute frequency sliding per trial, then pool
        all_inst_freq = [];
        example_fs = [];
        for trl = 1:length(dat_stim_broad.trial)
            comp_ts = topComp' * double(dat_stim_broad.trial{trl});

            % Hilbert transform
            analytic = hilbert(comp_ts);
            phase_angles = angle(analytic);

            % Temporal derivative of phase -> instantaneous frequency
            dphi = diff(unwrap(phase_angles));

            % Scale to Hz
            inst_freq = dphi * fsample / (2 * pi);

            % Median filter to remove phase-slip artifacts (kernel = 10 samples)
            inst_freq = medfilt1(inst_freq, 10);

            % Keep only values within gamma range (physiological constraint)
            valid = inst_freq >= gamma_range(1) & inst_freq <= gamma_range(2);
            all_inst_freq = [all_inst_freq, inst_freq(valid)];

            if trl == 1
                example_fs = inst_freq;
            end
        end

        all_fs_timeseries{cond, subj} = example_fs;

        % Compute histogram of instantaneous frequencies
        bin_edges   = gamma_range(1):1:gamma_range(2);
        bin_centers = bin_edges(1:end-1) + 0.5;
        counts = histcounts(all_inst_freq, bin_edges);
        counts = counts / sum(counts); % normalize to probability

        all_fs_hist_freqs{cond, subj}  = bin_centers;
        all_fs_hist_counts{cond, subj} = counts;

        % Peak = mode of distribution
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

    % --- Row 1: Eigenvalue spectra (Method 1) ---
    for cond = 1:4
        subplot(3, 4, cond); hold on;
        if ~isempty(all_eigspec{cond, subj})
            plot(scan_freqs, all_eigspec{cond, subj}, '-', ...
                'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            plot(scan_freqs, movmean(all_eigspec{cond, subj}, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2.5);
            xline(all_ged_peakfreq(cond, subj), '--', 'LineWidth', 2, ...
                'Color', colors(cond,:));
            text(all_ged_peakfreq(cond, subj) + 1, ...
                max(all_eigspec{cond, subj}) * 0.95, ...
                sprintf('%.0f Hz', all_ged_peakfreq(cond, subj)), ...
                'FontSize', 12, 'Color', colors(cond,:), 'FontWeight', 'bold');
        end
        xlabel('Frequency [Hz]'); ylabel('Eigenvalue (\lambda)');
        title(sprintf('%s — Scanning GED', condLabels{cond}), 'FontSize', 13);
        set(gca, 'FontSize', 12); xlim([30 90]); grid on; box on;
    end

    % --- Row 2: Frequency sliding distribution (Method 2) ---
    for cond = 1:4
        subplot(3, 4, 4 + cond); hold on;
        if ~isempty(all_fs_hist_counts{cond, subj})
            bar(all_fs_hist_freqs{cond, subj}, all_fs_hist_counts{cond, subj}, 1, ...
                'FaceColor', colors(cond,:), 'FaceAlpha', 0.6, 'EdgeColor', 'none');
            xline(all_fs_peakfreq(cond, subj), '--', 'LineWidth', 2, ...
                'Color', colors(cond,:));
            text(all_fs_peakfreq(cond, subj) + 1, ...
                max(all_fs_hist_counts{cond, subj}) * 0.95, ...
                sprintf('%.0f Hz', all_fs_peakfreq(cond, subj)), ...
                'FontSize', 12, 'Color', colors(cond,:), 'FontWeight', 'bold');
        end
        xlabel('Frequency [Hz]'); ylabel('P(f)');
        title(sprintf('%s — Freq Sliding', condLabels{cond}), 'FontSize', 13);
        set(gca, 'FontSize', 12); xlim([30 90]); grid on; box on;
    end

    % --- Row 3: Topographies at GED peak freq ---
    for cond = 1:4
        subplot(3, 4, 8 + cond);
        if ~isempty(all_ged_topos{cond, subj})
            topo_data = [];
            topo_data.label  = chanlocs_all;
            topo_data.avg    = all_ged_topos{cond, subj};
            topo_data.dimord = 'chan';
            cfg = [];
            cfg.layout   = layANThead;
            cfg.comment  = 'no';
            cfg.marker   = 'off';
            cfg.style    = 'straight';
            cfg.gridscale = 300;
            cfg.zlim     = 'maxabs';
            cfg.colormap = 'parula';
            try
                ft_topoplotER(cfg, topo_data);
            catch
                imagesc(topo_data.avg); colorbar;
            end
            title(sprintf('%s — Topo at %.0f Hz', condLabels{cond}, ...
                all_ged_peakfreq(cond, subj)), 'FontSize', 13);
        end
    end

    saveas(fig, fullfile(fig_save_dir, sprintf('GED_scan_subj%s.png', subjects{subj})));

end % subject loop

%% ====================================================================
%  GRAND AVERAGE FIGURE 1: Eigenvalue spectra overlay (all conditions)
%  ====================================================================
fprintf('\nCreating grand average figures...\n');
close all

fig_ga1 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Narrowband Scanning GED — Grand Average (N=%d)', nSubj), ...
    'FontSize', 20, 'FontWeight', 'bold');

% --- Left panel: Eigenvalue spectra per condition ---
subplot(1, 2, 1); hold on;
for cond = 1:4
    % Collect all eigenvalue spectra for this condition
    eigspec_mat = nan(nSubj, length(scan_freqs));
    for s = 1:nSubj
        if ~isempty(all_eigspec{cond, s})
            eigspec_mat(s,:) = all_eigspec{cond, s};
        end
    end
    mu  = nanmean(eigspec_mat, 1);
    sem = nanstd(eigspec_mat, [], 1) / sqrt(sum(~isnan(eigspec_mat(:,1))));

    mu_smooth  = movmean(mu, 5);
    sem_smooth = movmean(sem, 5);

    faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
    patch([scan_freqs, fliplr(scan_freqs)], ...
          [mu_smooth - sem_smooth, fliplr(mu_smooth + sem_smooth)], ...
          colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(scan_freqs, mu_smooth, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);

    [~, pi] = max(mu_smooth);
    xline(scan_freqs(pi), '--', 'Color', colors(cond,:), 'LineWidth', 1.5);
end
xlabel('Frequency [Hz]'); ylabel('Eigenvalue (\lambda)');
title('Eigenvalue Spectra (mean +/- SEM)', 'FontSize', 15);
legend(condLabels, 'Location', 'northeast', 'FontSize', 13);
set(gca, 'FontSize', 14); xlim([30 90]); grid on; box on;

% --- Right panel: Individual peak frequencies (scatter + boxplot style) ---
subplot(1, 2, 2); hold on;
positions = [1 2 3 4];
for cond = 1:4
    pf = all_ged_peakfreq(cond, :);
    pf = pf(~isnan(pf));
    scatter(ones(size(pf)) * positions(cond) + 0.15*(rand(size(pf))-0.5), pf, ...
        60, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.6);
    plot(positions(cond), nanmean(all_ged_peakfreq(cond,:)), 'kd', ...
        'MarkerSize', 12, 'MarkerFaceColor', colors(cond,:), 'LineWidth', 1.5);
    % Error bar for SEM
    mu_pf  = nanmean(all_ged_peakfreq(cond,:));
    sem_pf = nanstd(all_ged_peakfreq(cond,:)) / sqrt(sum(~isnan(all_ged_peakfreq(cond,:))));
    errorbar(positions(cond), mu_pf, sem_pf, 'k', 'LineWidth', 2, 'CapSize', 10);
end
set(gca, 'XTick', positions, 'XTickLabel', condLabels, 'FontSize', 14);
ylabel('Peak Gamma Frequency [Hz]');
title('Individual Peak Frequencies (GED)', 'FontSize', 15);
ylim([30 90]); grid on; box on;

saveas(fig_ga1, fullfile(fig_save_dir, 'GED_scan_grand_average.png'));

%% ====================================================================
%  GRAND AVERAGE FIGURE 2: Method comparison (GED vs Freq Sliding)
%  ====================================================================
fig_ga2 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Method Comparison: Scanning GED vs Frequency Sliding (N=%d)', nSubj), ...
    'FontSize', 20, 'FontWeight', 'bold');

% --- Top left: GED peak freq per condition ---
subplot(2, 2, 1); hold on;
for cond = 1:4
    pf = all_ged_peakfreq(cond, :);
    pf = pf(~isnan(pf));
    scatter(ones(size(pf))*cond + 0.15*(rand(size(pf))-0.5), pf, ...
        50, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.5);
    errorbar(cond, nanmean(pf), nanstd(pf)/sqrt(numel(pf)), ...
        'k', 'LineWidth', 2, 'CapSize', 10);
    plot(cond, nanmean(pf), 'kd', 'MarkerSize', 12, ...
        'MarkerFaceColor', colors(cond,:), 'LineWidth', 1.5);
end
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 14);
ylabel('Peak Frequency [Hz]'); title('Scanning GED', 'FontSize', 15);
ylim([30 90]); grid on; box on;

% --- Top right: Freq sliding peak freq per condition ---
subplot(2, 2, 2); hold on;
for cond = 1:4
    pf = all_fs_peakfreq(cond, :);
    pf = pf(~isnan(pf));
    scatter(ones(size(pf))*cond + 0.15*(rand(size(pf))-0.5), pf, ...
        50, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.5);
    errorbar(cond, nanmean(pf), nanstd(pf)/sqrt(numel(pf)), ...
        'k', 'LineWidth', 2, 'CapSize', 10);
    plot(cond, nanmean(pf), 'kd', 'MarkerSize', 12, ...
        'MarkerFaceColor', colors(cond,:), 'LineWidth', 1.5);
end
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 14);
ylabel('Peak Frequency [Hz]'); title('Frequency Sliding', 'FontSize', 15);
ylim([30 90]); grid on; box on;

% --- Bottom left: Correlation between methods ---
subplot(2, 2, 3); hold on;
all_ged_flat = all_ged_peakfreq(:);
all_fs_flat  = all_fs_peakfreq(:);
valid = ~isnan(all_ged_flat) & ~isnan(all_fs_flat);
for cond = 1:4
    scatter(all_ged_peakfreq(cond,:), all_fs_peakfreq(cond,:), ...
        80, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.7);
end
plot([30 90], [30 90], 'k--', 'LineWidth', 1.5);
if sum(valid) > 2
    [r, p] = corr(all_ged_flat(valid), all_fs_flat(valid), 'type', 'Spearman');
    text(32, 85, sprintf('r_s = %.2f, p = %.3f', r, p), 'FontSize', 14);
end
xlabel('GED Peak Freq [Hz]'); ylabel('Freq Sliding Peak Freq [Hz]');
title('Agreement Between Methods', 'FontSize', 15);
set(gca, 'FontSize', 14); xlim([30 90]); ylim([30 90]);
axis square; grid on; box on;
legend(condLabels, 'Location', 'southeast', 'FontSize', 12);

% --- Bottom right: Frequency sliding grand average distributions ---
subplot(2, 2, 4); hold on;
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
    mu_smooth = movmean(mu, 3);
    plot(bin_centers_common, mu_smooth, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);
end
xlabel('Frequency [Hz]'); ylabel('P(f)');
title('Grand Average Freq Sliding Distributions', 'FontSize', 15);
legend(condLabels, 'Location', 'northeast', 'FontSize', 13);
set(gca, 'FontSize', 14); xlim([30 90]); grid on; box on;

saveas(fig_ga2, fullfile(fig_save_dir, 'GED_method_comparison.png'));

%% ====================================================================
%  GRAND AVERAGE FIGURE 3: All subjects subplot (like your current fig)
%  ====================================================================
nRows = ceil(nSubj / 5);
fig_all = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Scanning GED Eigenvalue Spectra — All Subjects (N=%d)', nSubj), ...
    'FontSize', 18, 'FontWeight', 'bold');

for s = 1:nSubj
    subplot(nRows, 5, s); hold on;
    for cond = 1:4
        if ~isempty(all_eigspec{cond, s})
            plot(scan_freqs, movmean(all_eigspec{cond, s}, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2);
            xline(all_ged_peakfreq(cond, s), '--', 'Color', colors(cond,:), ...
                'LineWidth', 1.2);
        end
    end
    xlabel('Freq [Hz]'); ylabel('\lambda');
    title(sprintf('Subj %s', subjects{s}), 'FontSize', 11);
    set(gca, 'FontSize', 10); xlim([30 90]); grid on; box on;
    if s == 1
        legend(condLabels, 'FontSize', 8, 'Location', 'northeast');
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
    'all_eigspec', 'all_ged_topos', 'all_ged_peakfreq', 'all_ged_peakpow', ...
    'all_fs_hist_freqs', 'all_fs_hist_counts', 'all_fs_peakfreq', 'all_fs_median', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames');

fprintf('\n========================================\n');
fprintf('GED Gamma Peak Detection Complete!\n');
fprintf('Results saved to:\n  %s\n', save_path);
fprintf('Figures saved to:\n  %s\n', fig_save_dir);
fprintf('========================================\n');
