%% GCP Stats Overview
% Boxplots for all variables
% Percentage change barplots
% Gamma Differences vs. Gaze Differences

%% Load data
clc
clear
close all
data = readtable('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/merged_data.csv');
data.ReactionTime = data.ReactionTime .* 1000;
data.PupilSize = [];
variables = {'Accuracy', 'ReactionTime', 'GazeDeviation', 'GazeSTD', 'MSRate', 'GammaPower', 'GammaFreq'};
save_names = {'acc', 'rt', 'gazedev', 'gazestd', 'ms', 'pow', 'freq'};
colors = color_def('GCP');

%% BOXPLOTS for each variable
y_axis_labels = {'Accuracy [%]', 'Reaction Time [ms]', 'Gaze Deviation [px]', 'Gaze STD [px]', 'Microsaccade Rate [%]', 'Gamma Power [%]', 'Gamma Frequency [Hz]'};

% Unique subject identifiers
subjects = unique(data.ID);
font_size = 20;

for i = 1:length(variables)
    close all
    figure;
    set(gcf, 'Position', [100, 200, 1000, 800], 'Color', 'w');
    hold on;

    % Set axis limits
    ylim([min(data.(variables{i})) max(data.(variables{i}))])
    xlim([0.5 4.5])

    % Create boxplot
    boxplot(data.(variables{i}), data.Condition, 'Labels', {'25% Contrast', '50% Contrast', '75% Contrast', '100% Contrast'}, 'Colors', 'k');
    set(gca, 'FontSize', 30);

    % Define jitter
    jitterAmount = 0.05;
    for cond = 1:4
        cond_data = data(data.Condition == cond, :);
        x_jittered{cond} = cond_data.Condition + (rand(size(cond_data.(variables{i}))) - 0.5) * jitterAmount;
    end

    % Connect points with a line
    for subj = 1:length(subjects)
        % Extract data for this subject
        subj_data = data(data.ID == subjects(subj), :);

        if height(subj_data) == 2
            x = [x_jittered{1}'; x_jittered{2}'];
            x = x(:, subj);
            y = subj_data.(variables{i});
            plot(x, y, '-o', 'Color', [0.5, 0.5, 0.5, 0.5], 'LineWidth', 2);
        end
    end

    % Connect points for each subject
    for subj = 1:length(subjects)
        subj_data = data(data.ID == subjects(subj), :);
        if height(subj_data) == 4  % ensure this subject has all four conditions
            % Compute jittered x-positions for this subject
            x = subj_data.Condition + (rand(size(subj_data.Condition)) - 0.5) * jitterAmount;
            y = subj_data.(variables{i});
            % Sort by condition so lines go in order
            [x_sorted, sort_idx] = sort(x);
            y_sorted = y(sort_idx);
            % Plot thick, semi-transparent grey line with circle markers
            h = plot(x_sorted, y_sorted, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3);
            h.Color(4) = 0.2;   % set the line1s transparency after plotting
        end
    end

    % Scatter individual data points
    for cond = 1:4
        cond_data = data(data.Condition == cond, :);
        scatter(x_jittered{cond}, cond_data.(variables{i}), 300, 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', colors(cond, :), 'jitter', 'off', 'SizeData', 300);
    end

    % Add title and labels
    ylabel(y_axis_labels{i}, "FontSize", 30);
    set(gca, 'Box', 'off');
    ax = gca;
    ax.FontSize = font_size;
    title(variables{i}, "FontSize", 40);
    if i == 6
        title('Gamma Power', "FontSize", 40)
    elseif i == 7
        title('Gamma Frequency', "FontSize", 40)
    end

    % Add customization to microsaccade plot
    if i == 5
        yline(0, '--')
        box on
        title('Microsaccade Rate', "FontSize", 40);
        ylim([-60 60])
    end

    % Save individual subplot
    saveas(gcf, strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/stats/overview/GCP_stats_boxplots_', save_names{i}, '.png'));
end
