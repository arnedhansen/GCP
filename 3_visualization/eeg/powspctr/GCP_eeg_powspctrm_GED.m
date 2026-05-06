%% GCP GED Power Spectrum Visualizations
%
% Recreates GED power-spectrum figures from the saved feature file.
% Includes:
%   1) Grand-average raw spectra across subjects (condition-wise)
%   2) All-subject raw spectra subplots
%   3) Final combined GED branch spectra (early/full/late windows)
%
% Data source:
%   data/features/GCP_eeg_GED.mat
%
% Outputs:
%   1) Subject-wise figures with 3 window panels
%      - GCP_eeg_GED_powspctrm_subj<id>.png
%   2) Window-wise all-subject overview figures
%      - GCP_eeg_GED_powspctrm_early_ALL.png
%      - GCP_eeg_GED_powspctrm_full_ALL.png
%      - GCP_eeg_GED_powspctrm_late_ALL.png

%% Setup
startup
[subjects, paths, colors] = setup('GCP');

%% Paths
data_path = fullfile(paths.features, 'GCP_eeg_GED.mat');
fig_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/powspctrm';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Load GED feature file
if ~isfile(data_path)
    error('GED feature file not found: %s', data_path);
end
S = load(data_path);

required_vars = { ...
    'all_trial_powratio', ...
    'all_trial_powratio_bench', ...
    'all_trial_powratio_bench_early', ...
    'all_trial_powratio_bench_late', ...
    'all_trial_median_single', ...
    'benchmark_methods'};
for vi = 1:numel(required_vars)
    if ~isfield(S, required_vars{vi})
        error('Variable "%s" is missing in %s', required_vars{vi}, data_path);
    end
end

if isfield(S, 'scan_freqs') && ~isempty(S.scan_freqs)
    scan_freqs = S.scan_freqs;
else
    scan_freqs = 30:90;
end
analysis_freq_range = [scan_freqs(1), scan_freqs(end)];

if isfield(S, 'condLabels') && ~isempty(S.condLabels)
    condLabels = S.condLabels;
else
    condLabels = {'25%', '50%', '75%', '100%'};
end
nCond = numel(condLabels);

if isfield(S, 'subjects') && ~isempty(S.subjects)
    subj_labels = S.subjects;
else
    subj_labels = subjects;
end
nSubj = min(numel(subj_labels), numel(subjects));
subj_labels = subj_labels(1:nSubj);

%% Grand-average raw spectrum
close all
fig_grand_psd = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

grand_line_handles = gobjects(1, nCond);
smooth_n = 3; % mild, uniform smoothing across full spectrum
plot_order = nCond:-1:1; % 100%, 75%, 50%, 25%
for oi = 1:numel(plot_order)
    cond = plot_order(oi);
    subj_curves = nan(nSubj, numel(scan_freqs));
    for s = 1:nSubj
        pr_mat_full = S.all_trial_powratio{cond, s};
        if ~isempty(pr_mat_full)
            subj_curves(s, :) = nanmean(pr_mat_full, 1);
        end
    end

    mu = nanmean(subj_curves, 1);
    n_valid_subj = sum(~isnan(subj_curves(:, 1)));
    if n_valid_subj > 0
        se = nanstd(subj_curves, [], 1) ./ sqrt(n_valid_subj);
    else
        se = nan(size(mu));
    end

    if any(isfinite(mu))
        mu_plot = movmean(mu, smooth_n, 'omitnan');
        se_plot = movmean(se, smooth_n, 'omitnan');
        good = isfinite(scan_freqs) & isfinite(mu_plot) & isfinite(se_plot);
        if sum(good) < 3
            continue;
        end

        hse = shadedErrorBar(scan_freqs(good), mu_plot(good), se_plot(good));
        set(hse.mainLine, 'Color', colors(cond, :), 'LineWidth', 5);
        set(hse.patch, 'FaceColor', colors(cond, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        set(hse.edge, 'Visible', 'off');
        grand_line_handles(cond) = hse.mainLine;

        if any(isfinite(mu_plot))
            [~, peak_idx] = max(mu_plot);
            if ~isempty(peak_idx) && isfinite(scan_freqs(peak_idx))
                xline(scan_freqs(peak_idx), '--', 'Color', colors(cond, :), 'LineWidth', 2.4);
            end
        end
    end
end
xlim([30 90]);
ylim([0 2.5]);
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
title('Grand Average Power Spectrum (dB)', 'FontSize', 25, 'FontWeight', 'bold');
set(gca, 'FontSize', 20);
valid_handles = isgraphics(grand_line_handles);
if any(valid_handles)
    legend_order = plot_order(isgraphics(grand_line_handles(plot_order)));
    legend(grand_line_handles(legend_order), condLabels(legend_order), ...
        'Location', 'best', 'FontSize', 20, 'Box','off');
end
drawnow
exportgraphics(fig_grand_psd, fullfile(fig_dir, 'GCP_eeg_GED_grand_average_power_spectrum.png'), 'Resolution', 600);

%% All-subject raw spectra subplots
nRows = ceil(nSubj / 5);
fig_all = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Trial-Level Mean Raw Spectra: All Subjects (N=%d)', nSubj), ...
    'FontSize', 16, 'FontWeight', 'bold');
fig_all_legend_handles = gobjects(1, nCond);

for s = 1:nSubj
    subplot(nRows, 5, s); hold on;
    subj_panel_maxabs = 0;
    subj_panel_min = inf;
    subj_panel_max = -inf;
    peak_freq_txt = strings(nCond, 1);

    for cond = 1:nCond
        pr_mat_full = S.all_trial_powratio{cond, s};
        if ~isempty(pr_mat_full)
            mu_raw = nanmean(pr_mat_full, 1);
            subj_panel_maxabs = max(subj_panel_maxabs, max(abs(mu_raw), [], 'omitnan'));
            subj_panel_min = min(subj_panel_min, min(mu_raw, [], 'omitnan'));
            subj_panel_max = max(subj_panel_max, max(mu_raw, [], 'omitnan'));
            h_cond = plot(scan_freqs, movmean(mu_raw, 5), '-', ...
                'Color', colors(cond, :), 'LineWidth', 2);
            if s == 1
                fig_all_legend_handles(cond) = h_cond;
            end
            md_pf = S.all_trial_median_single(cond, s);
            if ~isnan(md_pf)
                xline(md_pf, '--', 'Color', colors(cond, :), 'LineWidth', 1.2);
                peak_freq_txt(cond) = sprintf('%s: %.0fHz', condLabels{cond}, md_pf);
            else
                peak_freq_txt(cond) = sprintf('%s: n/a', condLabels{cond});
            end
        end
    end

    yline(0, 'k-', 'LineWidth', 0.5);
    xlabel('Freq [Hz]');
    ylabel('Power [dB]');
    title(sprintf('Subj %s', subj_labels{s}), 'FontSize', 11, 'Interpreter', 'none');

    if ~isfinite(subj_panel_min), subj_panel_min = -10; end
    if ~isfinite(subj_panel_max), subj_panel_max = 10; end
    y_lo = subj_panel_min - 0.05 * abs(subj_panel_max - subj_panel_min);
    y_hi = subj_panel_max + 0.05 * abs(subj_panel_max - subj_panel_min);
    if ~isfinite(y_lo) || ~isfinite(y_hi) || y_hi <= y_lo
        y_lo = -10;
        y_hi = 10;
    end
    set(gca, 'FontSize', 10); xlim([30 90]); ylim([y_lo y_hi]); box on;

    text(0.97, 0.97, strjoin(cellstr(peak_freq_txt), newline), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'Color', 'k', 'FontSize', 8);

    if s == 1
        valid_handles = isgraphics(fig_all_legend_handles);
        if any(valid_handles)
            legend(fig_all_legend_handles(valid_handles), condLabels(valid_handles), ...
                'FontSize', 8, 'Location', 'northwest');
        end
    end
end
exportgraphics(fig_all, fullfile(fig_dir, 'GCP_eeg_GED_all_subjects_power_spectrum.png'), 'Resolution', 600);

%% Final combined GED method index
combined_method_idx = find(strcmpi(S.benchmark_methods, 'ged_combined_artifact_weighted'), 1, 'first');
if isempty(combined_method_idx)
    nBenchmarkMethods = numel(S.benchmark_methods);
    combined_method_idx = min(3, nBenchmarkMethods);
end

time_windows = {'early', 'full', 'late'};
time_cells = {S.all_trial_powratio_bench_early, S.all_trial_powratio_bench, S.all_trial_powratio_bench_late};

%% Subject-wise figures (three windows per subject)
for s = 1:nSubj
    fig_ps = figure('Position', [0 0 1512 982], 'Color', 'w');
    for wi = 1:3
        subplot(1, 3, wi); hold on;
        panel_min = inf;
        panel_max = -inf;

        for cond = 1:nCond
            pr_mat = time_cells{wi}{combined_method_idx, cond, s};
            if isempty(pr_mat)
                continue;
            end

            mu = mean(pr_mat, 1, 'omitnan');
            sem = std(pr_mat, 0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(pr_mat), 1)));
            good = isfinite(mu) & isfinite(sem);
            if sum(good) < 3
                continue;
            end

            x = scan_freqs(good);
            y = mu(good);
            e = sem(good);
            patch([x, fliplr(x)], [y - e, fliplr(y + e)], colors(cond, :), ...
                'FaceAlpha', 0.18, 'EdgeColor', 'none');
            plot(x, y, '-', 'Color', colors(cond, :), 'LineWidth', 2.0);

            panel_min = min(panel_min, min(y - e));
            panel_max = max(panel_max, max(y + e));
        end

        yline(0, 'k--', 'LineWidth', 0.8);
        xlim(analysis_freq_range);
        if isfinite(panel_min) && isfinite(panel_max) && panel_max > panel_min
            yr = panel_max - panel_min;
            ylim([panel_min - 0.08 * yr, panel_max + 0.12 * yr]);
        end
        xlabel('Frequency [Hz]');
        ylabel('Power [dB]');
        title(sprintf('%s window', upper(time_windows{wi})), 'Interpreter', 'none');
        set(gca, 'FontSize', 11, 'Box', 'on');
    end

    legend(condLabels, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
    sgtitle(sprintf('Final Combined GED Power Spectra: Subject %s', subj_labels{s}), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');

    out_subj = fullfile(fig_dir, sprintf('GCP_eeg_GED_powspctrm_subj%s.png', subj_labels{s}));
    exportgraphics(fig_ps, out_subj, 'Resolution', 600);
end

%% Window-wise all-subject figures
for wi = 1:3
    fig_all_win = figure('Position', [0 0 1512 982], 'Color', 'w');
    for s = 1:nSubj
        subplot(ceil(nSubj / 5), 5, s); hold on;
        panel_min = inf;
        panel_max = -inf;

        for cond = 1:nCond
            pr_mat = time_cells{wi}{combined_method_idx, cond, s};
            if isempty(pr_mat)
                continue;
            end

            mu = mean(pr_mat, 1, 'omitnan');
            if ~any(isfinite(mu))
                continue;
            end

            plot(scan_freqs, mu, '-', 'Color', colors(cond, :), 'LineWidth', 1.8);
            panel_min = min(panel_min, min(mu, [], 'omitnan'));
            panel_max = max(panel_max, max(mu, [], 'omitnan'));
        end

        yline(0, 'k--', 'LineWidth', 0.6);
        xlim(analysis_freq_range);
        if isfinite(panel_min) && isfinite(panel_max) && panel_max > panel_min
            yr = panel_max - panel_min;
            ylim([panel_min - 0.08 * yr, panel_max + 0.12 * yr]);
        end
        title(sprintf('Subj %s', subj_labels{s}), 'FontSize', 10, 'Interpreter', 'none');
        xlabel('Hz');
        ylabel('dB');
        set(gca, 'FontSize', 9, 'Box', 'on');
    end

    sgtitle(sprintf('Final Combined GED Power Spectra (%s, all subjects)', upper(time_windows{wi})), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
    out_all = fullfile(fig_dir, sprintf('GCP_eeg_GED_powspctrm_%s_ALL.png', time_windows{wi}));
    exportgraphics(fig_all_win, out_all, 'Resolution', 600);
end

clc
fprintf('Saved final GED power-spectrum figures in:\n  %s\n', fig_dir);

function y_out = smooth_reflective_edges(y_in, freqs, central_range, central_n, edge_n)
if nargin < 5
    edge_n = central_n;
end

y_out = y_in;
if isempty(y_in) || numel(y_in) ~= numel(freqs)
    return;
end

idx_center = freqs >= central_range(1) & freqs <= central_range(2);
idx_edges = ~idx_center;

if any(idx_center)
    y_out(idx_center) = movmean(y_in(idx_center), central_n, 'omitnan');
end
if any(idx_edges)
    y_out(idx_edges) = movmean(y_in(idx_edges), edge_n, 'omitnan');
end
end
