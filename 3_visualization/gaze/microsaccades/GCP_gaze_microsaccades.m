%% GCP Microsaccades

% Visualizations:
%   Boxplots of microsaccade rate
%   Barplots of microsaccade rate percentage difference

%% Setup
clear
clc
close all

% path to features
path_feat = '/Volumes/methlab/Students/Arne/GCP/data/features/';
load(fullfile(path_feat,'gaze_matrix_bl.mat'),'gaze_data_bl');

% extract unique subjects and number of loads
allIDs     = [gaze_data_bl.ID];
subjIDs    = unique(allIDs);
nSub       = numel(subjIDs);
conditions = {'c25','c50','c75','c100'};
nConds     = numel(conditions);

% preallocate matrix: rows=subjects, cols=conditions
pctMS = nan(nSub, nConds);

% fill matrix
for i = 1:numel(gaze_data_bl)
    entry   = gaze_data_bl(i);
    subjIdx = find(subjIDs == entry.ID);
    condIdx = entry.Condition;    % assuming Condition 1→c25, …,4→c100
    pctMS(subjIdx, condIdx) = entry.PctMSRate;
end

%% 1) Boxplots of % Change by Load Condition
close all
figure;
set(gcf, 'Position', [200, 150, 1000, 600], 'Color', 'w');
hold on;

% draw boxplot
boxplot(pctMS, 'Labels', conditions, 'Widths', 0.6);

% connect individual subjects
for s = 1:nSub
    plot(1:nConds, pctMS(s,:), '-o', ...
         'Color', [0.5 0.5 0.5], 'MarkerFaceColor', 'w');
end

% zero‐line
yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

% labels & styling
xlabel('Load Condition', 'FontSize', 16);
ylabel('Microsaccade Rate % Change from Baseline', 'FontSize', 16);
title('Microsaccade‐Rate Percentage Change by Load', 'FontSize', 18);
set(gca, 'FontSize', 14, 'LineWidth', 1.5);

% save
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/microsaccades/GCP_gaze_pctchange_boxplot.png');
hold off;

%% 2) Barplot of Mean % Change with SEM
means = mean(pctMS, 'omitnan');
sems  = std(pctMS,  'omitnan') ./ sqrt(sum(~isnan(pctMS),1));

figure;
set(gcf, 'Position', [250, 200, 900, 500], 'Color', 'w');
hold on;

% bar + errorbars
b = bar(1:nConds, means, 'FaceColor', 'flat');
errorbar(1:nConds, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1.2);

% zero‐line
yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

% labels & styling
set(gca, 'XTick', 1:nConds, 'XTickLabel', conditions, 'FontSize', 14);
xlabel('Load Condition', 'FontSize', 16);
ylabel('Mean Microsaccade Rate % Change', 'FontSize', 16);
title('Mean % Change in Microsaccade Rate per Load', 'FontSize', 18);
set(gca, 'LineWidth', 1.5);

% save
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/microsaccades/GCP_gaze_pctchange_barplot.png');
hold off;

%% Display Summary Statistics
fprintf('\nMean %% Change (PctMSRate) by Condition:\n');
for k = 1:nConds
    fprintf('  %s: Mean = %.2f%%, SEM = %.2f%%\n', ...
            conditions{k}, means(k), sems(k));
end
