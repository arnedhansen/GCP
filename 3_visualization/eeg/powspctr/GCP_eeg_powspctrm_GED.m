%% GCP GED Final Power Spectrum Visualizations
%
% Recreates the final GED power-spectrum figures from the saved feature file.
% Only the final combined GED branch is visualized (early/full/late windows).
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
fig_dir = fullfile(paths.figures, 'eeg', 'powspctrm', 'ged');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Load GED feature file
if ~isfile(data_path)
    error('GED feature file not found: %s', data_path);
end
S = load(data_path);

required_vars = { ...
    'all_trial_powratio_bench', ...
    'all_trial_powratio_bench_early', ...
    'all_trial_powratio_bench_late', ...
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
