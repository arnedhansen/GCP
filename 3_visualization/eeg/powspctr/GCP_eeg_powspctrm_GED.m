%% GCP GED Power Spectrum Visualizations, full window only

startup
[subjects, paths, colors] = setup('GCP');

data_path = fullfile(paths.features, 'GCP_eeg_GED.mat');
fig_dir = fullfile(paths.figures, 'eeg', 'powspctrm', 'ged');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

if ~isfile(data_path)
    error('GED feature file not found: %s', data_path);
end
dat = load(data_path);

required_vars = {'all_condition_powspctrm', 'all_condition_peak_freq', 'scan_freqs', 'condLabels'};
for vi = 1:numel(required_vars)
    if ~isfield(dat, required_vars{vi})
        error('Variable "%s" is missing in %s', required_vars{vi}, data_path);
    end
end

scan_freqs = dat.scan_freqs;
analysis_freq_range = [scan_freqs(1), scan_freqs(end)];
condLabels = dat.condLabels;
nCond = numel(condLabels);

if isfield(dat, 'subjects') && ~isempty(dat.subjects)
    subj_labels = dat.subjects;
else
    subj_labels = subjects;
end
nSubj = numel(subj_labels);

fig_grand = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
for cond = 1:nCond
    subj_curves = nan(nSubj, numel(scan_freqs));
    for s = 1:nSubj
        curv = dat.all_condition_powspctrm{cond, s};
        if ~isempty(curv) && numel(curv) == numel(scan_freqs)
            subj_curves(s, :) = curv(:)';
        end
    end
    med_curve = nanmedian(subj_curves, 1);
    mad_curve = nan(1, numel(scan_freqs));
    for fi = 1:numel(scan_freqs)
        mad_curve(fi) = mad(subj_curves(:, fi), 1);
    end
    good = isfinite(scan_freqs) & isfinite(med_curve) & isfinite(mad_curve);
    if sum(good) < 3
        continue;
    end
    x = scan_freqs(good);
    y = med_curve(good);
    e = mad_curve(good);
    fill([x, fliplr(x)], [y - e, fliplr(y + e)], colors(cond, :), ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(x, y, '-', 'Color', colors(cond, :), 'LineWidth', 3.0, 'DisplayName', condLabels{cond});
    peak_hz = nanmedian(dat.all_condition_peak_freq(cond, :));
    if isfinite(peak_hz)
        xline(peak_hz, ':', 'Color', colors(cond, :), 'LineWidth', 1.6, 'HandleVisibility', 'off');
    end
end
yline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlim(analysis_freq_range);
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
title('GED Condition Averaged Spectra, Full Window', 'FontWeight', 'bold');
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
set(gca, 'Box', 'on', 'FontSize', 12);
exportgraphics(fig_grand, fullfile(fig_dir, 'GCP_eeg_GED_powspctrm_full_group.png'), 'Resolution', 600);

fig_subj = figure('Position', [0 0 1512 982], 'Color', 'w');
nRows = ceil(nSubj / 5);
for s = 1:nSubj
    subplot(nRows, 5, s);
    hold on;
    panel_min = inf;
    panel_max = -inf;
    for cond = 1:nCond
        curv = dat.all_condition_powspctrm{cond, s};
        if isempty(curv) || numel(curv) ~= numel(scan_freqs)
            continue;
        end
        plot(scan_freqs, curv(:)', '-', 'Color', colors(cond, :), 'LineWidth', 1.8);
        panel_min = min(panel_min, min(curv, [], 'omitnan'));
        panel_max = max(panel_max, max(curv, [], 'omitnan'));
        pk = dat.all_condition_peak_freq(cond, s);
        if isfinite(pk)
            xline(pk, ':', 'Color', colors(cond, :), 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
    end
    yline(0, 'k--', 'LineWidth', 0.6, 'HandleVisibility', 'off');
    xlim(analysis_freq_range);
    if isfinite(panel_min) && isfinite(panel_max) && panel_max > panel_min
        yr = panel_max - panel_min;
        ylim([panel_min - 0.08 * yr, panel_max + 0.12 * yr]);
    end
    title(sprintf('Subj %s', subj_labels{s}), 'Interpreter', 'none', 'FontSize', 10);
    set(gca, 'FontSize', 9, 'Box', 'on');
end
legend(condLabels, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
sgtitle('GED Full Window Spectra by Subject', 'FontWeight', 'bold');
exportgraphics(fig_subj, fullfile(fig_dir, 'GCP_eeg_GED_powspctrm_full_subjects.png'), 'Resolution', 600);

fprintf('Saved full window GED power spectrum figures to: %s\n', fig_dir);
