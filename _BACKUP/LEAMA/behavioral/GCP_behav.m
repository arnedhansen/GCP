%% Visualization of behavioral data for GCP

% Visualizations:
%   Boxplots of accuracy

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load RT and Acc data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/behavioral');
    cd(datapath);
    load('acc.mat');
    load('rt.mat');
    % acc1(subj) = l1acc;
    % acc2(subj) = l2acc;
    % acc3(subj) = l3acc;
    % acc4(subj) = l4acc;
    % acc5(subj) = l5acc;
    % acc6(subj) = l6acc;
    % acc7(subj) = l7acc;
    % acc8(subj) = l8acc;
    % rt1(subj) = l1rt*1000;
    % rt2(subj) = l2rt*1000;
    % rt3(subj) = l3rt*1000;
    % rt4(subj) = l4rt*1000;
    % rt5(subj) = l5rt*1000;
    % rt6(subj) = l6rt*1000;
    % rt7(subj) = l7rt*1000;
    % rt8(subj) = l8rt*1000;
    hcacc(subj) = hc_acc;
    lcacc(subj) = lc_acc;
end

%% Plot Accuracy BOXPLOTS
dataAcc = [hcacc', lcacc'];
conditions = {'High Contrast', 'Low Contrast'};
close all

figure;
set(gcf, 'Position', [0, 0, 1000, 600], 'Color', 'w');
colors = {'r', 'b'};
hold on;

% Boxplots
boxplot(dataAcc, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
for subj = 1:length(subjects)
    plot(1:length(conditions), dataAcc(subj, :), '-o', 'Color', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'w');
end

% Scatter plot for individual points
scatterHandles = gobjects(1, length(conditions));
for condIdx = 1:length(conditions)
    jitter = 0.05 * randn(1, length(subjects));
    scatterHandles(condIdx) = scatter(repelem(condIdx, length(subjects)) + jitter, dataAcc(:, condIdx), ...
                                      100, colors{condIdx}, 'filled', 'MarkerEdgeColor', 'k');
end

yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xlabel('Conditions', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Accuracy [%]', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions, 'FontSize', 15);
set(gca, 'LineWidth', 1.5);
set(gca, 'XLim', [0.5 length(conditions) + 0.5]);
set(gca, 'YLim', [min(dataAcc(:)) * 0.9 101]);

legend(scatterHandles, conditions, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
title('Grating Task Accuracy', 'FontName', 'Arial', 'FontSize', 25);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/behavioral/GCP_acc_boxplot.png');

% Stats
means = mean(dataAcc, 'omitnan');
stds = std(dataAcc, 'omitnan');
save /Volumes/methlab/Students/Arne/GCP/data/features/accuracy means stds

% Display the results
disp('Means:');
disp(means);
disp('Stds:');
disp(stds);
