%% GCP GED + Spectral Centroid: Hybrid Approach (Common Spatial Filter)
%
% Combines the strengths of both methods:
%   - GED spatial filtering for optimal gamma signal isolation
%     (data-driven weighted combination of all channels)
%   - Spectral centroid for robust frequency estimation
%     (center of mass, not peak-picking)
%
% Uses a COMMON spatial filter across all contrast conditions:
%   Phase 1 (per subject): Pool covariance matrices across all 4 conditions
%           and run GED once -> one spatial filter per subject
%   Phase 2 (per condition): Apply the common filter to each condition
%           separately for power spectrum and centroid computation
%
% This guarantees the same gamma source is compared across contrast levels,
% avoiding the topography instability seen with per-condition GED.

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

% GED parameters
lambda = 0.01; % shrinkage regularization

% Multitaper parameters for power spectrum
taper_smoothing = 3; % ±3 Hz smoothing (dpss)

% Condition info
condNames  = {'c25', 'c50', 'c75', 'c100'};
trialCodes = [61, 62, 63, 64];
condLabels = {'25%', '50%', '75%', '100%'};

% Figure save directory
if ispc
    fig_save_dir = 'W:\Students\Arne\GCP\figures\eeg\ged_centroid';
    feat_root    = 'W:\Students\Arne\GCP\data\features';
else
    fig_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged_centroid';
    feat_root    = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features';
end
if ~exist(fig_save_dir, 'dir'), mkdir(fig_save_dir); end

%% Preallocate
centroid_freq  = nan(4, nSubj);
mean_power     = nan(4, nSubj);
ged_eigenvalue = nan(1, nSubj);  % single eigenvalue per subject (pooled)

all_powspec_stim = cell(4, nSubj);
all_powspec_base = cell(4, nSubj);
all_powspec_diff = cell(4, nSubj);
all_topos        = cell(1, nSubj);  % one topo per subject (common filter)

freq_common = [];

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

    nChans = length(dataEEG_c25.label);

    %% ================================================================
    %  PHASE 1: Build POOLED covariance across all conditions -> one GED
    %  ================================================================
    clc
    fprintf('Subject %s (%d/%d) — Phase 1: Pooled GED\n', ...
        subjects{subj}, subj, nSubj);

    covStim_pooled = zeros(nChans);
    covBase_pooled = zeros(nChans);
    nTrials_total  = 0;

    dat_per_cond = cell(1, 4);  % store selected data for Phase 2

    for cond = 1:4
        dat    = dataStructs{cond};
        trlIdx = trialIndices{cond};
        if isempty(trlIdx), continue; end

        cfg = [];
        cfg.trials = trlIdx;
        dat = ft_selectdata(cfg, dat);
        dat_per_cond{cond} = dat;

        % Broadband gamma filter
        cfg_filt = [];
        cfg_filt.bpfilter   = 'yes';
        cfg_filt.bpfreq     = gamma_range;
        cfg_filt.bpfilttype = 'fir';
        cfg_filt.bpfiltord  = round(3 * fsample / gamma_range(1));
        dat_gamma = ft_preprocessing(cfg_filt, dat);

        % Stimulus and baseline segments
        cfg_t = [];
        cfg_t.latency = baseline_window;
        dat_base = ft_selectdata(cfg_t, dat_gamma);

        cfg_t.latency = stimulus_window;
        dat_stim = ft_selectdata(cfg_t, dat_gamma);

        % Accumulate covariance (weighted by trial count)
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

    % GED on pooled covariance -> common spatial filter
    [W, D] = eig(covStim_pooled, covBase_pooled);
    [evals_sorted, sortIdx] = sort(real(diag(D)), 'descend');
    W = W(:, sortIdx);

    topComp = W(:, 1);
    ged_eigenvalue(subj) = evals_sorted(1);

    % Sign correction
    topo_temp = covStim_pooled * topComp;
    [~, mxI] = max(abs(topo_temp));
    if topo_temp(mxI) < 0, topComp = -topComp; end

    all_topos{subj} = covStim_pooled * topComp;

    %% ================================================================
    %  PHASE 2: Apply common filter to each condition separately
    %  ================================================================
    for cond = 1:4
        fprintf('Subject %s (%d/%d) — Phase 2: Condition %d/4\n', ...
            subjects{subj}, subj, nSubj, cond);

        dat = dat_per_cond{cond};
        if isempty(dat), continue; end

        % Select stimulus and baseline from UNFILTERED data
        cfg_t_raw = [];
        cfg_t_raw.latency = stimulus_window;
        dat_stim_raw = ft_selectdata(cfg_t_raw, dat);

        cfg_t_raw.latency = baseline_window;
        dat_base_raw = ft_selectdata(cfg_t_raw, dat);

        % Per-trial multitaper power spectra
        nTrl_stim = length(dat_stim_raw.trial);
        nTrl_base = length(dat_base_raw.trial);

        stim_dur = diff(stimulus_window);
        base_dur = diff(baseline_window);
        freq_res = 1 / min(stim_dur, base_dur);

        freqs_fft = 0 : freq_res : fsample/2;
        gamma_idx = freqs_fft >= gamma_range(1) & freqs_fft <= gamma_range(2);
        freqs_gamma = freqs_fft(gamma_idx);

        if isempty(freq_common)
            freq_common = freqs_gamma;
        end

        nTapers = max(1, floor(2 * taper_smoothing * min(stim_dur, base_dur) - 1));

        % Stimulus spectra through common spatial filter
        pow_stim_trials = zeros(nTrl_stim, sum(gamma_idx));
        for trl = 1:nTrl_stim
            comp_ts = topComp' * double(dat_stim_raw.trial{trl});
            n = length(comp_ts);
            tapers = dpss(n, taper_smoothing * n / fsample, nTapers);
            P_trl = zeros(1, length(freqs_fft));
            for ti = 1:nTapers
                Y = fft(comp_ts .* tapers(:, ti)', length(freqs_fft)*2 - 1);
                P_trl = P_trl + abs(Y(1:length(freqs_fft))).^2;
            end
            P_trl = P_trl / nTapers;
            pow_stim_trials(trl, :) = P_trl(gamma_idx);
        end
        pow_stim_mean = mean(pow_stim_trials, 1);

        % Baseline spectra through common spatial filter
        pow_base_trials = zeros(nTrl_base, sum(gamma_idx));
        for trl = 1:nTrl_base
            comp_ts = topComp' * double(dat_base_raw.trial{trl});
            n = length(comp_ts);
            tapers = dpss(n, taper_smoothing * n / fsample, nTapers);
            P_trl = zeros(1, length(freqs_fft));
            for ti = 1:nTapers
                Y = fft(comp_ts .* tapers(:, ti)', length(freqs_fft)*2 - 1);
                P_trl = P_trl + abs(Y(1:length(freqs_fft))).^2;
            end
            P_trl = P_trl / nTapers;
            pow_base_trials(trl, :) = P_trl(gamma_idx);
        end
        pow_base_mean = mean(pow_base_trials, 1);

        % Baseline correction (dB)
        pow_diff = 10 * log10(pow_stim_mean ./ pow_base_mean);

        all_powspec_stim{cond, subj} = pow_stim_mean;
        all_powspec_base{cond, subj} = pow_base_mean;
        all_powspec_diff{cond, subj} = pow_diff;

        % Spectral centroid on rectified baseline-corrected spectrum
        pow_rect = max(pow_diff, 0);
        if sum(pow_rect) > 0
            centroid_freq(cond, subj) = sum(freqs_gamma .* pow_rect) / sum(pow_rect);
        else
            centroid_freq(cond, subj) = NaN;
        end

        mean_power(cond, subj) = mean(pow_diff);

    end % condition loop

    %% ================================================================
    %  Per-subject figure
    %  ================================================================
    close all
    fig = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('GED + Centroid (Common Filter): Subject %s', subjects{subj}), ...
        'FontSize', 20, 'FontWeight', 'bold');

    % Row 1: Power spectra with centroid markers (4 conditions)
    for cond = 1:4
        subplot(2, 4, cond); hold on;
        if ~isempty(all_powspec_diff{cond, subj})
            plot(freq_common, all_powspec_diff{cond, subj}, '-', ...
                'Color', colors(cond,:), 'LineWidth', 2.5);
            yline(0, 'k-', 'LineWidth', 0.5);
            if ~isnan(centroid_freq(cond, subj))
                xline(centroid_freq(cond, subj), '--', 'LineWidth', 2, ...
                    'Color', colors(cond,:));
                text(centroid_freq(cond, subj) + 1, ...
                    max(all_powspec_diff{cond, subj}) * 0.85, ...
                    sprintf('%.1f Hz', centroid_freq(cond, subj)), ...
                    'FontSize', 12, 'Color', colors(cond,:), 'FontWeight', 'bold');
            end
        end
        xlabel('Freq [Hz]'); ylabel('Power [dB]');
        title(sprintf('%s Power Spectrum', condLabels{cond}), 'FontSize', 12);
        set(gca, 'FontSize', 11); xlim(gamma_range); grid on; box on;
    end

    % Row 2, left: Common topography
    cfg_topo = [];
    cfg_topo.layout    = headmodel.layANThead;
    cfg_topo.comment   = 'no';
    cfg_topo.marker    = 'off';
    cfg_topo.style     = 'straight';
    cfg_topo.gridscale = 300;
    cfg_topo.zlim      = 'maxabs';
    cfg_topo.colormap  = '*RdBu';
    cfg_topo.figure    = 'gcf';

    subplot(2, 4, 5);
    if ~isempty(all_topos{subj})
        topo_data = [];
        topo_data.label  = chanlocs_all;
        topo_data.avg    = all_topos{subj};
        topo_data.dimord = 'chan';
        try
            ft_topoplotER(cfg_topo, topo_data);
            cb = colorbar; cb.FontSize = 9;
        catch
            imagesc(topo_data.avg); colorbar;
        end
        title(sprintf('Common Filter (\\lambda=%.1f)', ged_eigenvalue(subj)), 'FontSize', 12);
    end

    % Row 2, center-right: All 4 conditions overlaid
    subplot(2, 4, [6 7 8]); hold on;
    hl_subj = gobjects(1, 4);
    for cond = 1:4
        if ~isempty(all_powspec_diff{cond, subj})
            hl_subj(cond) = plot(freq_common, all_powspec_diff{cond, subj}, '-', ...
                'Color', colors(cond,:), 'LineWidth', 2.5);
            xline(centroid_freq(cond, subj), '--', 'Color', colors(cond,:), 'LineWidth', 1.5);
        end
    end
    yline(0, 'k-', 'LineWidth', 0.5);
    xlabel('Freq [Hz]'); ylabel('Power [dB]');
    title('All Conditions (Same Spatial Filter)', 'FontSize', 14);
    legend(hl_subj, condLabels, 'FontSize', 11, 'Location', 'best');
    set(gca, 'FontSize', 12); xlim(gamma_range); grid on; box on;

    saveas(fig, fullfile(fig_save_dir, sprintf('GED_centroid_subj%s.png', subjects{subj})));

end % subject loop

%% ====================================================================
%  Grand average power spectrum
%  ====================================================================
fprintf('\nCreating grand average figures...\n');
close all

fig_ga = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

hl = gobjects(1, 4);
for cond = 1:4
    pow_mat = nan(nSubj, length(freq_common));
    for s = 1:nSubj
        if ~isempty(all_powspec_diff{cond, s})
            pow_mat(s, :) = all_powspec_diff{cond, s};
        end
    end
    mu  = nanmean(pow_mat, 1);
    sem = nanstd(pow_mat, [], 1) / sqrt(sum(~isnan(pow_mat(:,1))));

    faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
    patch([freq_common, fliplr(freq_common)], [mu - sem, fliplr(mu + sem)], ...
        colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    hl(cond) = plot(freq_common, mu, '-', 'Color', colors(cond,:), 'LineWidth', 3);

    ga_centroid = nanmean(centroid_freq(cond, :));
    xline(ga_centroid, '--', 'Color', colors(cond,:), 'LineWidth', 2);
end

yline(0, 'k--', 'LineWidth', 1);
set(gca, 'FontSize', 16);
xlim(gamma_range);
xlabel('Frequency [Hz]', 'FontSize', 20);
ylabel('Power [dB]', 'FontSize', 20);
legend(hl, condLabels, 'FontSize', 14, 'Location', 'best');
title(sprintf('GED + Centroid: Grand Average Power Spectrum (N=%d)', nSubj), ...
    'FontSize', 22, 'FontWeight', 'bold');

saveas(fig_ga, fullfile(fig_save_dir, 'GED_centroid_grand_average.png'));

%% ====================================================================
%  All subjects subplot
%  ====================================================================
close all
nCols = 5;
nRows = ceil(nSubj / nCols);
fig_all = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('GED + Centroid: All Subjects (N=%d)', nSubj), ...
    'FontSize', 18, 'FontWeight', 'bold');

for s = 1:nSubj
    subplot(nRows, nCols, s); hold on;
    for cond = 1:4
        if ~isempty(all_powspec_diff{cond, s})
            plot(freq_common, all_powspec_diff{cond, s}, '-', ...
                'Color', colors(cond,:), 'LineWidth', 2);
            xline(centroid_freq(cond, s), '--', 'Color', colors(cond,:), ...
                'LineWidth', 1.2);
        end
    end
    yline(0, 'k-', 'LineWidth', 0.5);
    xlabel('Freq [Hz]'); ylabel('Power [dB]');
    title(sprintf('Subj %s', subjects{s}), 'FontSize', 11);
    set(gca, 'FontSize', 10); xlim(gamma_range); grid on; box on;
    if s == 1
        legend(condLabels, 'FontSize', 8, 'Location', 'best');
    end
end
saveas(fig_all, fullfile(fig_save_dir, 'GED_centroid_all_subjects.png'));

%% ====================================================================
%  Boxplot: Centroid frequency
%  ====================================================================
close all
fig_box1 = figure('Position', [0 0 1200 982], 'Color', 'w');
hold on;

for s = 1:nSubj
    cf = centroid_freq(:, s);
    if sum(~isnan(cf)) >= 2
        plot(1:4, cf, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    end
end

y_box = centroid_freq(:);
g_box = repelem((1:4)', nSubj, 1);
valid = ~isnan(y_box);
boxplot(y_box(valid), g_box(valid), 'Colors', 'k', 'Symbol', '');

hold on;
for c = 1:4
    cf = centroid_freq(c, :);
    cf = cf(~isnan(cf));
    xJit = c + (rand(size(cf)) - 0.5) * 0.1;
    scatter(xJit, cf, 250, colors(c,:), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

xlim([0.5 4.5]); ylim([30 90]);
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 20, 'Box', 'off');
ylabel('Gamma Centroid Frequency [Hz]');
title('GED + Centroid Frequency', 'FontSize', 30, 'FontWeight', 'bold');

saveas(fig_box1, fullfile(fig_save_dir, 'GED_centroid_boxplot_freq.png'));

%% ====================================================================
%  Boxplot: Mean gamma power
%  ====================================================================
close all
fig_box2 = figure('Position', [0 0 1200 982], 'Color', 'w');
hold on;

for s = 1:nSubj
    mp = mean_power(:, s);
    if sum(~isnan(mp)) >= 2
        plot(1:4, mp, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    end
end

y_box2 = mean_power(:);
g_box2 = repelem((1:4)', nSubj, 1);
valid2 = ~isnan(y_box2);
boxplot(y_box2(valid2), g_box2(valid2), 'Colors', 'k', 'Symbol', '');

hold on;
for c = 1:4
    mp = mean_power(c, :);
    mp = mp(~isnan(mp));
    xJit = c + (rand(size(mp)) - 0.5) * 0.1;
    scatter(xJit, mp, 250, colors(c,:), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
xlim([0.5 4.5]);
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 20, 'Box', 'off');
ylabel('Mean Gamma Power [dB]');
title('GED + Mean Gamma Power (30-90 Hz)', 'FontSize', 30, 'FontWeight', 'bold');

saveas(fig_box2, fullfile(fig_save_dir, 'GED_centroid_boxplot_power.png'));

%% ====================================================================
%  Strip chart: GED eigenvalue (one value per subject, pooled filter)
%  ====================================================================
close all
fig_box3 = figure('Position', [0 0 800 982], 'Color', 'w');
hold on;

ev_valid = ged_eigenvalue(~isnan(ged_eigenvalue));
xJit = 1 + (rand(size(ev_valid)) - 0.5) * 0.3;
scatter(xJit, ev_valid, 300, [0.3 0.3 0.3], 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
boxplot(ev_valid, 'Colors', 'k', 'Symbol', '', 'Positions', 1);

mu_ev  = nanmean(ged_eigenvalue);
sem_ev = nanstd(ged_eigenvalue) / sqrt(sum(~isnan(ged_eigenvalue)));
errorbar(1.4, mu_ev, sem_ev, 'k', 'LineWidth', 2, 'CapSize', 12);
plot(1.4, mu_ev, 'kd', 'MarkerSize', 14, 'MarkerFaceColor', [0.5 0.5 0.5], 'LineWidth', 1.5);

xlim([0.4 1.8]);
set(gca, 'XTick', 1, 'XTickLabel', {'Pooled'}, 'FontSize', 20, 'Box', 'off');
ylabel('Top Eigenvalue (\lambda)');
title('Pooled Gamma SNR (GED Eigenvalue)', 'FontSize', 26, 'FontWeight', 'bold');

saveas(fig_box3, fullfile(fig_save_dir, 'GED_centroid_boxplot_eigenvalue.png'));

%% ====================================================================
%  Save results
%  ====================================================================
save(fullfile(feat_root, 'ged_centroid.mat'), ...
    'centroid_freq', 'mean_power', 'ged_eigenvalue', ...
    'all_powspec_stim', 'all_powspec_base', 'all_powspec_diff', ...
    'all_topos', 'freq_common', ...
    'subjects', 'condLabels', 'condNames', 'chanlocs_all');

eeg_ged_centroid_data = [];
for subj_i = 1:nSubj
    subject_id = repmat(str2double(subjects{subj_i}), 4, 1);
    conditions = (1:4)';
    cfreqs     = centroid_freq(:, subj_i);
    mpows      = mean_power(:, subj_i);
    evals      = repmat(ged_eigenvalue(subj_i), 4, 1);

    subj_data = struct('ID', num2cell(subject_id), ...
        'Condition',        num2cell(conditions), ...
        'GEDCentroidFreq',  num2cell(cfreqs), ...
        'GEDMeanPower',     num2cell(mpows), ...
        'GEDEigenvalue',    num2cell(evals));

    savepath = fullfile(feat_root, subjects{subj_i}, 'eeg');
    if ~exist(savepath, 'dir'), mkdir(savepath); end
    save(fullfile(savepath, 'eeg_ged_centroid_subj.mat'), 'subj_data');

    eeg_ged_centroid_data = [eeg_ged_centroid_data; subj_data];
end
save(fullfile(feat_root, 'eeg_ged_centroid_matrix.mat'), 'eeg_ged_centroid_data');

%% Print summary
clc
fprintf('=== GED + Centroid (Common Filter) Summary (N=%d) ===\n\n', nSubj);
fprintf('Pooled eigenvalue: %.2f ± %.2f (mean ± SEM)\n\n', ...
    nanmean(ged_eigenvalue), nanstd(ged_eigenvalue)/sqrt(sum(~isnan(ged_eigenvalue))));
fprintf('%-10s  %12s  %12s\n', 'Condition', 'Centroid [Hz]', 'Power [dB]');
fprintf('%s\n', repmat('-', 1, 38));
for cond = 1:4
    fprintf('%-10s  %8.1f ± %-4.1f  %8.2f ± %-4.2f\n', ...
        condLabels{cond}, ...
        nanmean(centroid_freq(cond,:)), nanstd(centroid_freq(cond,:))/sqrt(sum(~isnan(centroid_freq(cond,:)))), ...
        nanmean(mean_power(cond,:)),    nanstd(mean_power(cond,:))/sqrt(sum(~isnan(mean_power(cond,:)))));
end
fprintf('\nDone.\n');
