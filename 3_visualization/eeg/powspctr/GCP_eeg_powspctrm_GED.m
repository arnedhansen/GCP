%% GCP GED Power Spectrum Visualizations
%
%  Grand average and single subjects (smoothed and unsmoothed).
%  Data: GCP_eeg_powspctrm_GED.mat (FieldTrip freq per subject).

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('GCP');
data_path = fullfile(paths.features, 'GCP_eeg_powspctrm_GED.mat');
fig_dir = fullfile(paths.figures, 'eeg', 'powspctrm');

%% Load data
dat = load(data_path);
scan_freqs = dat.scan_freqs;
analysis_freq_range = [scan_freqs(1), scan_freqs(end)];
condLabels = dat.condLabels;
nCond = numel(condLabels);
nSubj = numel(subjects);
peak_hz = dat.all_condition_peak_freq_full;

%% Grand average powerspectrum (smoothed)
close all
fig_grand = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
for cond = 1:nCond
    cfg = [];
    cfg.keepindividual = 'yes';
    cfg.parameter = 'powspctrm';
    ga = ft_freqgrandaverage(cfg, dat.freq_powspctrm_full{cond, :}{:});
    nf = numel(ga.freq);
    pow = reshape(double(ga.powspctrm), size(ga.powspctrm, 1), nf);
    n_fin = sum(isfinite(pow), 1);
    ym = mean(pow, 1, 'omitnan');
    ye = std(pow, 0, 1, 'omitnan') ./ sqrt(max(1, n_fin));
    fx = ga.freq(:).';
    lineProps = {'-', 'Color', colors(cond, :), 'LineWidth', 4, ...
        'DisplayName', condLabels{cond}};
    seb = shadedErrorBar(fx, ym, ye, 'lineProps', lineProps, ...
        'transparent', true, 'patchSaturation', 0.15);
    set(seb.edge, 'Visible', 'off');
end
xlim(analysis_freq_range);
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
title('');
legend('Location', 'best', 'Box', 'off');
set(gca, 'Box', 'off', 'FontSize', 20);
drawnow; pause(0.05);
exportgraphics(fig_grand, fullfile(fig_dir, 'GCP_eeg_GED_powspctrm_grand_average_smoothed.png'), 'Resolution', 600);

%% Single subjects powerspectra overview (smoothed)
close all
fig_subj = figure('Position', [0 0 1512 982], 'Color', 'w');
nRows = ceil(nSubj / 5);
for subj = 1:nSubj
    subplot(nRows, 5, subj);
    hold on;
    panel_min = inf;
    panel_max = -inf;
    for cond = 1:nCond
        fq = dat.freq_powspctrm_full{cond, subj};
        if isempty(fq)
            text(31, 0.1, 'NO ELIGIBLE GED COMPONENTS', 'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold');
            continue;
        end
        curv = squeeze(fq.powspctrm(1, :));
        if numel(curv) ~= numel(scan_freqs)
            text(31, 0.1, 'NO ELIGIBLE GED COMPONENTS', 'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold');
            continue;
        end
        plot(scan_freqs, curv(:)', '-', 'Color', colors(cond, :), 'LineWidth', 3);
        panel_min = min(panel_min, min(curv, [], 'omitnan'));
        panel_max = max(panel_max, max(curv, [], 'omitnan'));
        pk = peak_hz(cond, subj);
        if isfinite(pk)
            xline(pk, ':', 'Color', colors(cond, :), 'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end
    yline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    xlim(analysis_freq_range);
    xticks(30:10:90);
    if isfinite(panel_min) && isfinite(panel_max) && panel_max > panel_min
        ylim([0, panel_max * 1.1]);
    end
    peak_text_y = 0.95;
    peak_text_step = min(0.05, 0.95 / max(nCond, 1));
    for cond = 1:nCond
        pk = peak_hz(cond, subj);
        if ~isfinite(pk)
            text(31, 0.1, 'NO ELIGIBLE GED COMPONENTS', 'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold');
            continue;
        end
        text(0.98, peak_text_y - (cond - 1) * peak_text_step, sprintf('%.0f Hz', pk), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'Color', colors(cond, :), 'FontSize', 9, 'FontWeight', 'bold');
    end
    title(sprintf('Participant %d', subj));
    set(gca, 'FontSize', 15, 'Box', 'on');
    if subj == 1 || subj == 6
        legend(condLabels, 'Location', 'best', 'Box', 'off');
        ylabel('Power [dB]');
    end
    xlabel('Frequency [Hz]');
end
drawnow; pause(0.05);
exportgraphics(fig_subj, fullfile(fig_dir, 'GCP_eeg_GED_powspctrm_overview_subjects_smoothed.png'), 'Resolution', 600);

fprintf('Saved full window GED power spectrum figures (smoothed) to: %s\n', fig_dir);

%% Grand average powerspectrum (unsmoothed)
close all
fig_grand = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
for cond = 1:nCond
    cfg = [];
    cfg.keepindividual = 'yes';
    cfg.parameter = 'powspctrm';
    ga = ft_freqgrandaverage(cfg, dat.freq_powspctrm_full_unsmoothed{cond, :}{:});
    nf = numel(ga.freq);
    pow = reshape(double(ga.powspctrm), size(ga.powspctrm, 1), nf);
    n_fin = sum(isfinite(pow), 1);
    ym = mean(pow, 1, 'omitnan');
    ye = std(pow, 0, 1, 'omitnan') ./ sqrt(max(1, n_fin));
    fx = ga.freq(:).';
    lineProps = {'-', 'Color', colors(cond, :), 'LineWidth', 4, ...
        'DisplayName', condLabels{cond}};
    seb = shadedErrorBar(fx, ym, ye, 'lineProps', lineProps, ...
        'transparent', true, 'patchSaturation', 0.15);
    set(seb.edge, 'Visible', 'off');
end
xlim(analysis_freq_range);
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
title('');
legend('Location', 'best', 'Box', 'off');
set(gca, 'Box', 'off', 'FontSize', 20);
drawnow; pause(0.05);
exportgraphics(fig_grand, fullfile(fig_dir, 'GCP_eeg_GED_powspctrm_grand_average_unsmoothed.png'), 'Resolution', 600);

%% Single subjects powerspectra overview (unsmoothed)
close all
fig_subj = figure('Position', [0 0 1512 982], 'Color', 'w');
nRows = ceil(nSubj / 5);
for subj = 1:nSubj
    subplot(nRows, 5, subj);
    hold on;
    panel_min = inf;
    panel_max = -inf;
    for cond = 1:nCond
        fq = dat.freq_powspctrm_full_unsmoothed{cond, subj};
        if isempty(fq)
            text(31, 0.1, 'NO ELIGIBLE GED COMPONENTS', 'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold');
            continue;
        end
        curv = squeeze(fq.powspctrm(1, :));
        if numel(curv) ~= numel(scan_freqs)
            text(31, 0.1, 'NO ELIGIBLE GED COMPONENTS', 'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold');
            continue;
        end
        plot(scan_freqs, curv(:)', '-', 'Color', colors(cond, :), 'LineWidth', 3);
        panel_min = min(panel_min, min(curv, [], 'omitnan'));
        panel_max = max(panel_max, max(curv, [], 'omitnan'));
        pk = peak_hz(cond, subj);
        if isfinite(pk)
            xline(pk, ':', 'Color', colors(cond, :), 'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end
    yline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    xlim(analysis_freq_range);
    xticks(30:10:90);
    if isfinite(panel_min) && isfinite(panel_max) && panel_max > panel_min
        ylim([0, panel_max * 1.1]);
    end
    peak_text_y = 0.95;
    peak_text_step = min(0.05, 0.95 / max(nCond, 1));
    for cond = 1:nCond
        pk = peak_hz(cond, subj);
        if ~isfinite(pk)
            text(31, 0.1, 'NO ELIGIBLE GED COMPONENTS', 'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold');
            continue;
        end
        text(0.98, peak_text_y - (cond - 1) * peak_text_step, sprintf('%.0f Hz', pk), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'Color', colors(cond, :), 'FontSize', 9, 'FontWeight', 'bold');
    end
    title(sprintf('Participant %d', subj));
    set(gca, 'FontSize', 15, 'Box', 'on');
    if subj == 1 || subj == 6
        legend(condLabels, 'Location', 'best', 'Box', 'off');
        ylabel('Power [dB]');
    end
    xlabel('Frequency [Hz]');
end
drawnow; pause(0.05);
exportgraphics(fig_subj, fullfile(fig_dir, 'GCP_eeg_GED_powspctrm_overview_subjects_unsmoothed.png'), 'Resolution', 600);

fprintf('Saved full window GED power spectrum figures (unsmoothed) to: %s\n', fig_dir);
