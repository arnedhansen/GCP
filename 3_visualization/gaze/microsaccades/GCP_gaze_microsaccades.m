%% GCP Microsaccades

% Visualizations:
%   Boxplots of microsaccade rate

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load gaze data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    cd(datapath);
    load('ms_rate.mat');
    ms_lc(subj) = lc_msrate;
    ms_hc(subj) = hc_msrate;
end

%% Plot Euclidean deviations BOXPLOTS
dataMSrate = [ms_lc', ms_hc'];
conditions = {'Low Contrast', 'High Contrast'};
close all
figure;
set(gcf, 'Position', [250, 200, 1000, 600], 'Color', 'w');
colors = {'b', 'r'};
hold on;

% Boxplots
boxplot(dataMSrate, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
for subj = 1:length(subjects)
    plot(1:length(conditions), dataMSrate(subj, :), '-o', 'Color', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'w');
end

% Scatter plot for individual points
scatterHandles = gobjects(1, length(conditions));
for condIdx = 1:length(conditions)
    scatterHandles(condIdx) = scatter(repelem(condIdx, length(subjects)), dataMSrate(:, condIdx), 100, colors{condIdx}, 'filled', 'MarkerEdgeColor', 'k');
end

box off
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xlabel('Conditions', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Mean Microsaccade Rate [Microsaccades / Second]', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions, 'FontSize', 15);
set(gca, 'LineWidth', 1.5);
set(gca, 'XLim', [0.5 length(conditions) + 0.5]);
set(gca, 'YLim', [min(dataMSrate(:)) * 0.85 max(dataMSrate(:)) * 1.15]);

legend(scatterHandles, conditions, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
title('Sternberg Mean Microsaccade Rate', 'FontName', 'Arial', 'FontSize', 25);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/microsaccades/GCP_gaze_microsaccades_boxplot.png');

% Stats
means = mean(dataMSrate, 'omitnan');
Stds = std(dataMSrate, 'omitnan');

% Display the results
disp('Means:');
disp(means);
disp('Stds:');
disp(Stds);

%% Plot Euclidean deviations BARPLOT for PERCENTAGE CHANGE
close all
% Calculate Percentage Changes
percent_change = ((ms_hc - ms_lc) ./ ms_lc) * 100;

% Plot Percentage Changes BARPLOT
figure;
set(gcf, 'Position', [300, 250, 900, 500], 'Color', 'w');
hold on;

% Bar plot
bar(1:length(subjects), percent_change, 'FaceColor', 'black', 'EdgeColor', 'k', 'LineWidth', 1.2);
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');

% Customisation
xlabel('Subjects', 'FontName', 'Arial', 'FontSize', 16);
ylabel('Percentage Change [%]', 'FontName', 'Arial', 'FontSize', 16);
title('Percentage Change in Microsaccade Rate (HC vs LC)', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XTick', 1:length(subjects), 'XTickLabel', subjects, 'FontSize', 12, 'XTickLabelRotation', 45);
set(gca, 'LineWidth', 1.5);
set(gca, 'YLim', [min(percent_change) * 1.15, max(percent_change) * 1.15]);

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/microsaccades/GCP_gaze_percentage_change_barplot.png');
hold off;

% Display the results
disp('Percentage Changes:');
disp(percent_change);

