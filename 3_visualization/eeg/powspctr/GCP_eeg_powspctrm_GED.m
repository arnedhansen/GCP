%% GCP GED Power Spectrum Visualizations
%
%  Grand average and single subjects (smoothed and unsmoothed).

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
data_path = fullfile(paths.features, 'GCP_eeg_powspctrm_GED.mat');
fig_dir = fullfile(paths.figures, 'eeg', 'powspctrm');

%% Load data
dat = load(data_path);
nSubj = numel(subjects);
subj_idx = arrayfun(@(s) find(strcmp(dat.subjects, subjects{s}), 1), 1:nSubj);
scan_freqs = dat.scan_freqs;
analysis_freq_range = [scan_freqs(1), scan_freqs(end)];
condLabels = dat.condLabels;
nCond = numel(condLabels);
peak_hz = dat.all_condition_peak_freq_full;

%% Grand average powerspectrum (smoothed)
close all
fig_grand = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
for cond = 1:nCond
    cfg = [];
    cfg.keepindividual = 'yes';
    cfg.parameter = 'powspctrm';
    pow_per_subj = cell(1, numel(subj_idx));
    for i = 1:numel(subj_idx)
        si = subj_idx(i);
        if si <= size(dat.freq_powspctrm_full, 2)
            pow_per_subj{i} = dat.freq_powspctrm_full{cond, si};
        end
    end
    pow_per_subj = pow_per_subj(~cellfun(@isempty, pow_per_subj));
    if isempty(pow_per_subj)
        warning('GCP_eeg_powspctrm_GED:NoData', ...
            'No non-empty freq structs for condition %d (%s); skipping grand average trace.', ...
            cond, condLabels{cond});
        continue;
    end
    ga = ft_freqgrandaverage(cfg, pow_per_subj{:});
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
for i = 1:nSubj
    subj = subj_idx(i);
    subplot(nRows, 5, i);
    hold on;
    panel_min = inf;
    panel_max = -inf;
    plotted_any = false;
    for cond = 1:nCond
        fq = dat.freq_powspctrm_full{cond, subj};
        if isempty(fq)
            continue;
        end
        curv = squeeze(fq.powspctrm(1, :));
        if numel(curv) ~= numel(scan_freqs)
            continue;
        end
        plot(scan_freqs, curv(:)', '-', 'Color', colors(cond, :), 'LineWidth', 3);
        plotted_any = true;
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
            continue;
        end
        text(0.98, peak_text_y - (cond - 1) * peak_text_step, sprintf('%.0f Hz', pk), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'Color', colors(cond, :), 'FontSize', 9, 'FontWeight', 'bold');
    end
    title(sprintf('Participant %s', subjects{i}));
    set(gca, 'FontSize', 15, 'Box', 'on');
    if ~plotted_any
        text(31, 0.1, 'NO ELIGIBLE GED COMPONENTS', 'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold');
        leg_h = gobjects(nCond, 1);
        for cond = 1:nCond
            leg_h(cond) = patch(NaN, NaN, colors(cond, :), 'EdgeColor', 'none');
        end
        legend(leg_h, condLabels, 'Location', 'best', 'Box', 'off');
    end
    if i == 1 || i == 6
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
    pow_per_subj = cell(1, numel(subj_idx));
    for i = 1:numel(subj_idx)
        si = subj_idx(i);
        if si <= size(dat.freq_powspctrm_full_unsmoothed, 2)
            pow_per_subj{i} = dat.freq_powspctrm_full_unsmoothed{cond, si};
        end
    end
    pow_per_subj = pow_per_subj(~cellfun(@isempty, pow_per_subj));
    if isempty(pow_per_subj)
        warning('GCP_eeg_powspctrm_GED:NoData', ...
            'No non-empty freq structs for condition %d (%s); skipping grand average trace.', ...
            cond, condLabels{cond});
        continue;
    end
    ga = ft_freqgrandaverage(cfg, pow_per_subj{:});
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
for i = 1:nSubj
    subj = subj_idx(i);
    subplot(nRows, 5, i);
    hold on;
    panel_min = inf;
    panel_max = -inf;
    plotted_any = false;
    for cond = 1:nCond
        fq = dat.freq_powspctrm_full_unsmoothed{cond, subj};
        if isempty(fq)
            continue;
        end
        curv = squeeze(fq.powspctrm(1, :));
        if numel(curv) ~= numel(scan_freqs)
            continue;
        end
        plot(scan_freqs, curv(:)', '-', 'Color', colors(cond, :), 'LineWidth', 3);
        plotted_any = true;
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
            continue;
        end
        text(0.98, peak_text_y - (cond - 1) * peak_text_step, sprintf('%.0f Hz', pk), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'Color', colors(cond, :), 'FontSize', 9, 'FontWeight', 'bold');
    end
    title(sprintf('Participant %s', subjects{i}));
    set(gca, 'FontSize', 15, 'Box', 'on');
    if ~plotted_any
        text(31, 0.1, 'NO ELIGIBLE GED COMPONENTS', 'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold');
        leg_h = gobjects(nCond, 1);
        for cond = 1:nCond
            leg_h(cond) = patch(NaN, NaN, colors(cond, :), 'EdgeColor', 'none');
        end
        legend(leg_h, condLabels, 'Location', 'best', 'Box', 'off');
    end
    if i == 1 || i == 6
        ylabel('Power [dB]');
    end
    xlabel('Frequency [Hz]');
end
drawnow; pause(0.05);
exportgraphics(fig_subj, fullfile(fig_dir, 'GCP_eeg_GED_powspctrm_overview_subjects_unsmoothed.png'), 'Resolution', 600);

fprintf('Saved full window GED power spectrum figures (unsmoothed) to: %s\n', fig_dir);
