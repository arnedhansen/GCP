%% GCP GED Backprojected Power Spectra Visualization
%
% Visualizes trial-level GED backprojected power-ratio spectra from
% 2_feature_extraction/GED/GCP_eeg_GED.m outputs.
%
% Included visualizations:
%   1) Grand-average spectra with SEM:
%      - pre-Perm+CV weighted GED (method 3)
%      - post-Perm+CV weighted GED (method 4)
%      - raw and detrended variants
%   2) Subject-level detrended spectra (post-Perm+CV).

%% Setup
startup
[subjects, ~, colors] = setup('GCP');

%% Paths
if ispc
    data_path = 'W:\Students\Arne\GCP\data\features\GCP_eeg_GED.mat';
    fig_dir = 'W:\Students\Arne\GCP\figures\eeg\powspctrm\ged';
else
    data_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_eeg_GED.mat';
    fig_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/powspctrm/ged';
end
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Load GED features
if ~isfile(data_path)
    error('GED feature file not found: %s', data_path);
end
S = load(data_path);

required_vars = {'all_trial_powratio_bench', 'all_trial_powratio_dt_bench', 'benchmark_methods'};
for vi = 1:numel(required_vars)
    if ~isfield(S, required_vars{vi})
        error('Variable "%s" is missing in %s', required_vars{vi}, data_path);
    end
end

if isfield(S, 'scan_freqs')
    scan_freqs = S.scan_freqs;
else
    scan_freqs = 30:90;
end

if isfield(S, 'condLabels')
    condLabels = S.condLabels;
else
    condLabels = {'25%', '50%', '75%', '100%'};
end

nSubj = numel(subjects);
nCond = 4;
nFreq = numel(scan_freqs);

idx_pre = find(strcmpi(S.benchmark_methods, 'ged_preperm_hard_weighted'), 1, 'first');
idx_post = find(strcmpi(S.benchmark_methods, 'ged_permcv_weighted'), 1, 'first');

if isempty(idx_pre) || isempty(idx_post)
    error(['Could not identify benchmark methods "ged_preperm_hard_weighted" and ', ...
        '"ged_permcv_weighted" in benchmark_methods.']);
end

%% Aggregate subject means (condition x subject x frequency)
mu_pre_raw = nan(nCond, nSubj, nFreq);
mu_post_raw = nan(nCond, nSubj, nFreq);
mu_pre_dt = nan(nCond, nSubj, nFreq);
mu_post_dt = nan(nCond, nSubj, nFreq);

for subj = 1:nSubj
    for cond = 1:nCond
        pr_pre = S.all_trial_powratio_bench{idx_pre, cond, subj};
        pr_post = S.all_trial_powratio_bench{idx_post, cond, subj};
        dt_pre = S.all_trial_powratio_dt_bench{idx_pre, cond, subj};
        dt_post = S.all_trial_powratio_dt_bench{idx_post, cond, subj};

        if ~isempty(pr_pre)
            mu_pre_raw(cond, subj, :) = mean(pr_pre, 1, 'omitnan');
        end
        if ~isempty(pr_post)
            mu_post_raw(cond, subj, :) = mean(pr_post, 1, 'omitnan');
        end
        if ~isempty(dt_pre)
            mu_pre_dt(cond, subj, :) = mean(dt_pre, 1, 'omitnan');
        end
        if ~isempty(dt_post)
            mu_post_dt(cond, subj, :) = mean(dt_post, 1, 'omitnan');
        end
    end
end

%% Grand-average plot (raw + detrended, pre + post)
fig_grand = figure('Position', [0 0 1512 982]);
set(fig_grand, 'Color', 'w');
sgtitle('GED Backprojected Spectra (Grand Average)', 'FontSize', 18, 'FontWeight', 'bold');

subplot(2, 2, 1);
plot_group_curves(scan_freqs, mu_pre_raw, colors, condLabels, 'Pre-Perm+CV weighted (raw)', 'Power ratio');

subplot(2, 2, 2);
plot_group_curves(scan_freqs, mu_post_raw, colors, condLabels, 'Post-Perm+CV weighted (raw)', 'Power ratio');

subplot(2, 2, 3);
plot_group_curves(scan_freqs, mu_pre_dt, colors, condLabels, 'Pre-Perm+CV weighted (detrended)', '\Delta power ratio');

subplot(2, 2, 4);
plot_group_curves(scan_freqs, mu_post_dt, colors, condLabels, 'Post-Perm+CV weighted (detrended)', '\Delta power ratio');

saveas(fig_grand, fullfile(fig_dir, 'GCP_eeg_powspctrm_GED_grand_average.png'));

%% Subject-level detrended post-Perm+CV plot
nCols = 5;
nRows = ceil(nSubj / nCols);
fig_subj = figure('Position', [0 0 1512 982]);
set(fig_subj, 'Color', 'w');
sgtitle('GED Backprojected Spectra: Subject Means (Post-Perm+CV, detrended)', ...
    'FontSize', 16, 'FontWeight', 'bold');

legend_handles = gobjects(1, nCond);
for subj = 1:nSubj
    subplot(nRows, nCols, subj); hold on;
    y_abs = 0;
    for cond = 1:nCond
        y = squeeze(mu_post_dt(cond, subj, :))';
        if all(~isfinite(y))
            continue;
        end
        h = plot(scan_freqs, movmean(y, 5), '-', 'Color', colors(cond, :), 'LineWidth', 1.8);
        if subj == 1
            legend_handles(cond) = h;
        end
        y_abs = max(y_abs, max(abs(y), [], 'omitnan'));
    end
    yline(0, 'k-', 'LineWidth', 0.5);
    xlim([scan_freqs(1), scan_freqs(end)]);
    if y_abs > 0 && isfinite(y_abs)
        ylim([-y_abs, y_abs]);
    end
    set(gca, 'FontSize', 9, 'Box', 'on');
    title(sprintf('Subj %s', subjects{subj}), 'FontSize', 10);
    xlabel('Hz'); ylabel('\DeltaPR');
end

valid_h = isgraphics(legend_handles);
if any(valid_h)
    legend(legend_handles(valid_h), condLabels(valid_h), 'Location', 'bestoutside');
end

saveas(fig_subj, fullfile(fig_dir, 'GCP_eeg_powspctrm_GED_subjects_postPermCV_detrended.png'));

clc
fprintf('Saved GED backprojected power-spectrum figures:\n');
fprintf('  %s\n', fullfile(fig_dir, 'GCP_eeg_powspctrm_GED_grand_average.png'));
fprintf('  %s\n', fullfile(fig_dir, 'GCP_eeg_powspctrm_GED_subjects_postPermCV_detrended.png'));

%% Local function
function plot_group_curves(freqs, data_csf, colors, condLabels, ttl, ylab)
    hold on;
    nCond = size(data_csf, 1);
    h = gobjects(1, nCond);
    for cond = 1:nCond
        y_subj = squeeze(data_csf(cond, :, :)); % subjects x freqs
        mu = mean(y_subj, 1, 'omitnan');
        sem = std(y_subj, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(y_subj), 1));
        valid = isfinite(mu) & isfinite(sem);
        if ~any(valid)
            continue;
        end
        x = freqs(valid);
        y1 = mu(valid) - sem(valid);
        y2 = mu(valid) + sem(valid);
        fill([x, fliplr(x)], [y1, fliplr(y2)], colors(cond, :), ...
            'FaceAlpha', 0.16, 'EdgeColor', 'none');
        h(cond) = plot(freqs, movmean(mu, 5), '-', ...
            'Color', colors(cond, :), 'LineWidth', 2.4);
    end
    yline(0, 'k-', 'LineWidth', 0.5);
    xlim([freqs(1), freqs(end)]);
    xlabel('Frequency [Hz]');
    ylabel(ylab);
    title(ttl, 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'FontSize', 11, 'Box', 'on');
    valid_h = isgraphics(h);
    if any(valid_h)
        legend(h(valid_h), condLabels(valid_h), 'Location', 'best', 'FontSize', 9);
    end
end
