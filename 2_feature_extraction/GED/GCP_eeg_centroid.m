%% GCP Gamma Spectral Centroid Analysis
%
% Computes the spectral centroid (center of gravity) of the gamma band
% power spectrum as a more robust alternative to peak-picking via findpeaks.
%
% Spectral centroid = sum(f .* P(f)) / sum(P(f))
%
% Uses the already-computed FOOOF'd, baselined, smoothed power spectra
% from GCP_eeg_fex.m. Only the positive (above-baseline) portion of the
% spectrum contributes — this gives the center of mass of the gamma power
% increase relative to baseline.
%
% Extracted features per subject per condition:
%   - Gamma spectral centroid frequency [Hz]
%   - Mean gamma power [dB] (averaged over 30-90 Hz, occipital channels)

clear; close all; clc

%% Setup
startup
[subjects, path, colors, ~] = setup('GCP');

condLabels = {'25%', '50%', '75%', '100%'};
condNames  = {'c25', 'c50', 'c75', 'c100'};
nSubj = length(subjects);

gamma_range = [30 90];

if ispc
    fig_save_dir = 'W:\Students\Arne\GCP\figures\eeg\centroid';
    feat_root    = 'W:\Students\Arne\GCP\data\features';
else
    fig_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/centroid';
    feat_root    = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features';
end
if ~exist(fig_save_dir, 'dir'), mkdir(fig_save_dir); end

%% Define occipital channels (from first subject)
datapath = strcat(path, subjects{1}, '/eeg');
cd(datapath);
load('power_spectra.mat', 'pow_c25_fooof_bl_smooth');
occ_channels = {};
for i = 1:length(pow_c25_fooof_bl_smooth.label)
    label = pow_c25_fooof_bl_smooth.label{i};
    if contains(label, 'O') || contains(label, 'I')
        occ_channels{end+1} = label;
    end
end
fprintf('Using %d occipital channels.\n', numel(occ_channels));

%% Preallocate
centroid_freq = nan(4, nSubj);
mean_power    = nan(4, nSubj);

all_spectra   = cell(4, nSubj);  % occipital-averaged power spectra
freq_vec      = [];               % common frequency vector

%% Extract features
for subj = 1:nSubj
    clc
    fprintf('Spectral centroid: Subject %s (%d/%d)\n', subjects{subj}, subj, nSubj);

    load(fullfile(feat_root, subjects{subj}, 'eeg', 'power_spectra.mat'), ...
        'pow_c25_fooof_bl_smooth', 'pow_c50_fooof_bl_smooth', ...
        'pow_c75_fooof_bl_smooth', 'pow_c100_fooof_bl_smooth');

    pows = {pow_c25_fooof_bl_smooth, pow_c50_fooof_bl_smooth, ...
            pow_c75_fooof_bl_smooth, pow_c100_fooof_bl_smooth};

    for cond = 1:4
        pow = pows{cond};

        chan_idx = ismember(pow.label, occ_channels);
        freq_idx = pow.freq >= gamma_range(1) & pow.freq <= gamma_range(2);
        freqs = pow.freq(freq_idx);

        if isempty(freq_vec)
            freq_vec = freqs;
        end

        spectrum = mean(pow.powspctrm(chan_idx, freq_idx), 1);
        all_spectra{cond, subj} = spectrum;

        % Mean gamma power (full band, including negative values)
        mean_power(cond, subj) = mean(spectrum);

        % Spectral centroid on rectified spectrum (only above-baseline power)
        spectrum_rect = max(spectrum, 0);
        if sum(spectrum_rect) > 0
            centroid_freq(cond, subj) = sum(freqs .* spectrum_rect) / sum(spectrum_rect);
        else
            centroid_freq(cond, subj) = NaN;
        end
    end
end

%% ====================================================================
%  Per-subject control figures
%  ====================================================================
for subj = 1:nSubj
    close all
    fig = figure('Position', [0 0 1512 982], 'Color', 'w');
    hold on;

    for cond = 1:4
        spectrum = all_spectra{cond, subj};
        plot(freq_vec, spectrum, '-', 'Color', colors(cond,:), 'LineWidth', 3);
        xline(centroid_freq(cond, subj), '--', 'Color', colors(cond,:), 'LineWidth', 2);
    end

    yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5);
    max_abs = max(abs([all_spectra{:, subj}]));
    ylim([-max_abs*1.3, max_abs*1.3]);
    xlim(gamma_range);

    for cond = 1:4
        text(centroid_freq(cond, subj) + 0.5, -max_abs*(0.4 + 0.12*(cond-1)), ...
            sprintf('%.1f Hz', centroid_freq(cond, subj)), ...
            'FontSize', 16, 'Color', colors(cond,:), 'FontWeight', 'bold');
    end

    set(gca, 'FontSize', 18);
    xlabel('Frequency [Hz]', 'FontSize', 22);
    ylabel('Power [dB]', 'FontSize', 22);
    legend(condLabels, 'FontSize', 16, 'Location', 'best');
    title(sprintf('Subject %s: FOOOFed Power Spectrum + Centroid', subjects{subj}), ...
        'FontSize', 24, 'FontWeight', 'bold');

    saveas(fig, fullfile(fig_save_dir, sprintf('GCP_centroid_subj%s.png', subjects{subj})));
end

%% ====================================================================
%  Grand average power spectrum with centroid markers
%  ====================================================================
close all
fig_ga = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

hl = gobjects(1, 4);
for cond = 1:4
    spec_mat = cell2mat(all_spectra(cond, :)');
    mu  = mean(spec_mat, 1);
    sem = std(spec_mat, [], 1) / sqrt(nSubj);

    faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
    patch([freq_vec, fliplr(freq_vec)], [mu - sem, fliplr(mu + sem)], ...
        colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    hl(cond) = plot(freq_vec, mu, '-', 'Color', colors(cond,:), 'LineWidth', 3);

    ga_centroid = nanmean(centroid_freq(cond, :));
    xline(ga_centroid, '--', 'Color', colors(cond,:), 'LineWidth', 2);
end

yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5);
set(gca, 'FontSize', 18);
xlim(gamma_range);
xlabel('Frequency [Hz]', 'FontSize', 22);
ylabel('Power [dB]', 'FontSize', 22);
legend(hl, condLabels, 'FontSize', 16, 'Location', 'best');
title(sprintf('Grand Average Power Spectrum + Centroid (N=%d)', nSubj), ...
    'FontSize', 26, 'FontWeight', 'bold');

saveas(fig_ga, fullfile(fig_save_dir, 'GCP_centroid_grand_average.png'));

%% ====================================================================
%  All subjects subplot
%  ====================================================================
close all
nCols = 5;
nRows = ceil(nSubj / nCols);
fig_all = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Spectral Centroid: All Subjects (N=%d)', nSubj), ...
    'FontSize', 18, 'FontWeight', 'bold');

for subj = 1:nSubj
    subplot(nRows, nCols, subj); hold on;
    for cond = 1:4
        plot(freq_vec, all_spectra{cond, subj}, '-', ...
            'Color', colors(cond,:), 'LineWidth', 2);
        xline(centroid_freq(cond, subj), '--', 'Color', colors(cond,:), 'LineWidth', 1.2);
    end
    yline(0, 'k-', 'LineWidth', 0.5);
    xlabel('Freq [Hz]'); ylabel('Power');
    title(sprintf('Subj %s', subjects{subj}), 'FontSize', 11);
    set(gca, 'FontSize', 10); xlim(gamma_range); grid on; box on;
    if subj == 1
        legend(condLabels, 'FontSize', 8, 'Location', 'best');
    end
end
saveas(fig_all, fullfile(fig_save_dir, 'GCP_centroid_all_subjects.png'));

%% ====================================================================
%  Boxplot: Spectral centroid frequency
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
title('Spectral Centroid Frequency', 'FontSize', 30, 'FontWeight', 'bold');

saveas(fig_box1, fullfile(fig_save_dir, 'GCP_boxplot_centroid_freq.png'));

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
title('Mean Gamma Power (30-90 Hz)', 'FontSize', 30, 'FontWeight', 'bold');

saveas(fig_box2, fullfile(fig_save_dir, 'GCP_boxplot_mean_gamma_power.png'));

%% ====================================================================
%  Save results
%  ====================================================================
save(fullfile(feat_root, 'gamma_centroid.mat'), ...
    'centroid_freq', 'mean_power', 'all_spectra', 'freq_vec', ...
    'subjects', 'condLabels', 'condNames', 'occ_channels');

eeg_centroid_data = [];
for subj = 1:nSubj
    subject_id = repmat(str2double(subjects{subj}), 4, 1);
    conditions = (1:4)';
    cfreqs = centroid_freq(:, subj);
    mpows  = mean_power(:, subj);

    subj_data = struct('ID', num2cell(subject_id), ...
        'Condition',        num2cell(conditions), ...
        'CentroidFreq',     num2cell(cfreqs), ...
        'MeanGammaPower',   num2cell(mpows));

    savepath = fullfile(feat_root, subjects{subj}, 'eeg');
    if ~exist(savepath, 'dir'), mkdir(savepath); end
    save(fullfile(savepath, 'eeg_centroid_subj.mat'), 'subj_data');

    eeg_centroid_data = [eeg_centroid_data; subj_data];
end
save(fullfile(feat_root, 'eeg_centroid_matrix.mat'), 'eeg_centroid_data');

%% Print summary
clc
fprintf('=== Spectral Centroid Summary (N=%d) ===\n\n', nSubj);
fprintf('%-10s  %12s  %12s\n', 'Condition', 'Centroid [Hz]', 'Power [dB]');
fprintf('%s\n', repmat('-', 1, 38));
for cond = 1:4
    fprintf('%-10s  %8.1f ± %-4.1f  %8.4f ± %-6.4f\n', ...
        condLabels{cond}, ...
        nanmean(centroid_freq(cond,:)), nanstd(centroid_freq(cond,:))/sqrt(sum(~isnan(centroid_freq(cond,:)))), ...
        nanmean(mean_power(cond,:)), nanstd(mean_power(cond,:))/sqrt(sum(~isnan(mean_power(cond,:)))));
end
fprintf('\nDone.\n');
