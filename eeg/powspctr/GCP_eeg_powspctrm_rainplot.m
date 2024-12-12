%% RAINPLOT for GCP Power and Frequency data
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
for i = 1 %:length(variables)
    % Extract the current data variable
    dat = variables{i};
    dat = pow_dat
    max_vals = max(cellfun(@(x) max(abs(x)), dat));

    % Define colours for conditions
    [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
    colours = cb([5, 6], :); % Two colours: one for Condition 1 (low contrast), one for Condition 2 (high contrast)
    %colours = [0 0 1; 1 0 0];  % Blue for low contrast, Red for high contrast

    % Create raincloud plot
    figure;
    set(gcf, 'Position', [300, 250, 1600, 900], 'Color', 'w');
    % = rm_raincloud(data, colours, plot_top_to_bottom, raindrop_size)
    h = rm_raincloud(dat, colours, 0, 500);

    % Adjust plot aesthetics
    set(gca, 'YTick', [0, 1.8], 'YTickLabel', {'High Contrast', 'Low Contrast'}, 'FontSize', 20);
    %xlim([-max_vals*1.25 max_vals*2.5]) % for y-axis
    xlim([-1 2.25])
    xlabel('Gamma Power [db]', 'FontSize', 20);
    title(['Gamma Power' newline 'relative to baseline'], 'FontSize', 25);
end
