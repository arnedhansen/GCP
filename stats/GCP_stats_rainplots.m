%% RAINPLOTs for GCP data
clear
clc
close all

% Load data
data = readtable('/Volumes/methlab/Students/Arne/GCP/data/features/merged_data.csv');

% Split and combine data for all variables by condition
acc_dat = {data.Accuracy(data.Condition == 1), data.Accuracy(data.Condition == 2)};
rt_dat = {data.ReactionTime(data.Condition == 1), data.ReactionTime(data.Condition == 2)};
gazedev_dat = {data.GazeDeviation(data.Condition == 1), data.GazeDeviation(data.Condition == 2)};
pups_dat = {data.PupilSize(data.Condition == 1), data.PupilSize(data.Condition == 2)};
ms_dat = {data.MSRate(data.Condition == 1), data.MSRate(data.Condition == 2)};
pow_dat = {data.GammaPower(data.Condition == 1), data.GammaPower(data.Condition == 2)};
freq_dat = {data.GammaFreq(data.Condition == 1), data.GammaFreq(data.Condition == 2)};

%% Set up the figure
close all

% Choose data to plot
variables = {acc_dat, rt_dat, gazedev_dat, pups_dat, ms_dat, pow_dat, freq_dat};
labels = {'Accuracy', 'Reaction Time', 'Gaze Deviation', 'Pupil Size', 'Microsaccade Rate', 'Gamma Power', 'Gamma Frequency'};
for var = 1:length(variables)
    close all
    % Extract the current data variable
    dat = variables{var};
    % dat = pow_dat
    max_vals = max(cellfun(@(x) max(abs(x)), dat));

    % Define colours for conditions
    [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
    colours = cb([5, 6], :); % cb pastel colours (5 = bright blue; 6 = orange)

    % Create raincloud plot
    figure;
    set(gcf, 'Position', [300, 250, 1600, 900], 'Color', 'w');
    % = rm_raincloud(data, colours, add_boxplot, plot_top_to_bottom, raindrop_size, plot_mean_dots, connecting_lines)
    h = rm_raincloud(dat, colours, 1, 0, 500, 0, 1);

    % Adjust plot aesthetics
    xline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'Alpha', 0.25);
    if var == 1
        xlim([92.5 100])
    end
    %xlim([-max_vals*1.25 max_vals*2]) % for y-axis
    xlabel(labels(var), 'FontSize', 20);
    title([labels(var)], 'FontSize', 25);
    saveas(gcf, strcat('/Volumes/methlab/Students/Arne/GCP/figures/stats/GCP_stats_rainclouds_', labels{var}, '.png'))
end
