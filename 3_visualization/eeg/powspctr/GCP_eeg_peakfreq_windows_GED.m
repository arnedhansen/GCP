%% GCP GED Peak Gamma Frequency Across Windows
%
% Visualizes participant-level shifts in median peak gamma frequency across
% GED windows (early, late, full), separated by condition.
%
% Data source:
%   data/features/GCP_eeg_GED.mat
% Required variables:
%   all_trial_median_single_early, all_trial_median_single_late,
%   all_trial_median_single, condLabels

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
dat = load(data_path);

required_vars = { ...
    'all_trial_median_single_early', ...
    'all_trial_median_single_late', ...
    'all_trial_median_single'};
for vi = 1:numel(required_vars)
    if ~isfield(dat, required_vars{vi})
        error('Variable "%s" is missing in %s', required_vars{vi}, data_path);
    end
end

if isfield(dat, 'condLabels')
    condLabels = dat.condLabels;
else
    condLabels = {'25%', '50%', '75%', '100%'};
end

nCond = size(dat.all_trial_median_single, 1);
nSubj_data = size(dat.all_trial_median_single, 2);
nSubj = min(numel(subjects), nSubj_data);

window_labels = {'Early (0-0.6 s)', 'Late (1-2 s)', 'Full (0-2 s)'};
xw = 1:3;

fig_peak_windows = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Peak Gamma Frequency Across GED Windows (Subject Trajectories)', ...
    'FontSize', 18, 'FontWeight', 'bold');

for cond = 1:nCond
    subplot(2, 2, cond); hold on;

    peak_mat = [dat.all_trial_median_single_early(cond, 1:nSubj); ...
                dat.all_trial_median_single_late(cond, 1:nSubj); ...
                dat.all_trial_median_single(cond, 1:nSubj)]';

    % Participant-level trajectories
    for s = 1:nSubj
        y = peak_mat(s, :);
        if any(isfinite(y))
            plot(xw, y, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
            scatter(xw(isfinite(y)), y(isfinite(y)), 18, [0.65 0.65 0.65], ...
                'filled', 'MarkerFaceAlpha', 0.65, 'MarkerEdgeAlpha', 0.65);
        end
    end

    % Group mean +/- SEM
    mu_w = mean(peak_mat, 1, 'omitnan');
    sem_w = std(peak_mat, 0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(peak_mat), 1)));
    errorbar(xw, mu_w, sem_w, '-o', ...
        'Color', colors(cond, :), ...
        'MarkerFaceColor', colors(cond, :), ...
        'LineWidth', 2.5, 'CapSize', 8);

    n_valid = sum(any(isfinite(peak_mat), 2));
    title(sprintf('%s (N=%d)', condLabels{cond}, n_valid), 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Median peak gamma [Hz]');
    set(gca, 'XTick', xw, 'XTickLabel', window_labels, 'XTickLabelRotation', 15, ...
        'FontSize', 11, 'Box', 'on');
    xlim([0.7 3.3]);
    ylim([30 90]);
    grid on;
end

out_path = fullfile(fig_dir, 'GCP_eeg_GED_peak_gamma_windows_by_subject.png');
exportgraphics(fig_peak_windows, out_path, 'Resolution', 600);

clc
fprintf('Saved figure:\n  %s\n', out_path);
