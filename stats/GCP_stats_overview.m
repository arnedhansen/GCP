% Load your data into MATLAB, assuming you have converted it to a table format
clc
clear
close all
data = readtable('/Volumes/methlab/Students/Arne/GCP/data/features/merged_data.csv');
data.ReactionTime = data.ReactionTime .* 1000;
variables = {'Accuracy', 'ReactionTime', 'GazeDeviation', 'MSRate', 'GammaPower', 'GammaFreq'};

% Split the data by contrast condition
low_contrast = data(data.Condition == 1, :);
high_contrast = data(data.Condition == 2, :);

%% BOXPLOTS
close all
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
y_axis_labels = {'Accuracy [%]', 'Reaction Time [ms]', 'Gaze Deviation [px]', 'Microsaccade Rate [ms/s]', 'Gamma Power [dB]', 'Gamma Frequency [Hz]'};

% Unique subject identifiers
subjects = unique(data.ID);

for i = 1:length(variables)
    subplot(2, 3, i);
    hold on;

    % Set axis limits
    ylim([min(data.(variables{i})) max(data.(variables{i}))])
    xlim([0.5 2.5])

    % Create boxplot
    boxplot(data.(variables{i}), data.Condition, 'Labels', {'Low Contrast', 'High Contrast'});

    % Overlay individual data points and connect them
    for subj = 1:length(subjects)
        % Extract data for this subject
        subj_data = data(data.ID == subjects(subj), :);

        if height(subj_data) == 2 % Ensure the subject has data for both conditions
            % X coordinates: condition indices
            x = subj_data.Condition;
            % Y coordinates: variable values
            y = subj_data.(variables{i});

            % Connect points with a line
            plot(x, y, '-o', 'Color', [0.5, 0.5, 0.5, 0.5], 'LineWidth', 0.5);
        end
    end

    % Scatter individual data points with jitter
    jitterAmount = 0.001;  % Amount of jitter in the x-direction
    for cond = 1:2  % Two conditions (Low Contrast = 1, High Contrast = 2)
        % Filter data based on condition
        cond_data = data(data.Condition == cond, :);

        % Determine color based on condition
        if cond == 1
            markerColor = [0, 0, 1];  % Blue for Low Contrast
        else
            markerColor = [1, 0, 0];  % Red for High Contrast
        end

        % Add jitter and scatter the points
        x_jittered = cond_data.Condition + (rand(size(cond_data.(variables{i}))) - 0.5) * jitterAmount;  % Add random jitter to x-axis
        
        % Clear any previous scatter points (in case of overlaid dots)
        scatter(x_jittered, cond_data.(variables{i}), 36, 'MarkerEdgeColor', markerColor, ...
            'MarkerFaceColor', markerColor, 'jitter', 'off', 'SizeData', 38);
    end

    % Add title and labels
    title(variables{i}, "FontSize", 20);
    ylabel(y_axis_labels{i}, "FontSize", 15);
    hold off;
end
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/stats/GCP_stats_overview_boxplots.png');

%% PERCENTAGE CHANGE BARPLOTS
close all
% Preallocate percentage change matrix
percent_change = zeros(length(subjects), length(variables));

% Loop through each variable to calculate percentage change
for i = 1:length(variables)
    for subj = 1:length(subjects)
        % Extract data for this subject
        subj_data = data(data.ID == subjects(subj), :);

        if height(subj_data) == 2 % Ensure the subject has data for both conditions
            % Low and high contrast values
            low_value = subj_data{subj_data.Condition == 1, variables{i}};
            high_value = subj_data{subj_data.Condition == 2, variables{i}};

            % Calculate percentage change ((high - low) / low) * 100
            percent_change(subj, i) = ((high_value - low_value) / low_value) * 100;
        else
            percent_change(subj, i) = NaN; % Handle missing data
        end
    end
end

% Plot percentage change bar plots
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

for i = 1:length(variables)
    subplot(2, 3, i);
    hold on;

    % Bar plot for each participant
    bar(1:length(subjects), percent_change(:, i), 'FaceColor', 'k', 'EdgeColor', 'none');

    % Formatting
    xlim([0.5, length(subjects) + 0.5]);
    abw = max(abs([min(percent_change(:, i), [], 'omitnan'), max(percent_change(:, i), [], 'omitnan')]));
    ylim([-abw*1.25 abw*1.25]);
    if i == 5
        ylim([-100 100])
    end
    xticks(1:length(subjects));
    xticklabels(subjects);
    xlabel('Subjects');
    ylabel('% Change');
    title(variables{i}, 'FontSize', 20);
    hold off;
end
sgtitle('Percentage Change (HC - LC)', 'FontSize', 24);
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/stats/GCP_stats_overview_barplots_percentage_change.png');

%% Plot DIFFERENCES in GAMMA POWER and FREQUENCY against MICROSACCADES and GAZE DEVIATION
close all
% Calculate the differences in Gamma Power and Gamma Frequency
gamma_power_diff = percent_change(:, 5);
gamma_freq_diff = percent_change(:, 6);

% Extract the corresponding values for Gaze Deviation and Microsaccade Rate
gaze_deviation_diff = percent_change(:, 3);
ms_rate_diff = percent_change(:, 4);

% Set up the figure with the same aesthetics as before
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

% Plot Gamma Power difference vs Gaze Deviation difference
subplot(2, 2, 1);
hold on;
scatter(gaze_deviation_diff, gamma_power_diff, 36, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'SizeData', 100);
xlabel('Gaze Deviation Difference [%]', 'FontSize', 15);
ylabel('Gamma Power Difference [%]', 'FontSize', 15);

% Calculate max_abs_range for this subplot
max_abs_range = max(abs([gamma_power_diff; gaze_deviation_diff]), [], 'all');
xlim([-50 50]);
ylim([-100 100]);

% Dashed lines at x = 0 and y = 0
xline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'Alpha', 0.25);
yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'Alpha', 0.25);

title('Gamma Power vs Gaze Deviation', 'FontSize', 20);
hold off;

% Plot Gamma Power difference vs Microsaccade Rate difference
subplot(2, 2, 2);
hold on;
scatter(ms_rate_diff, gamma_power_diff, 36, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'SizeData', 100);
xlabel('Microsaccade Rate Difference [%]', 'FontSize', 15);
ylabel('Gamma Power Difference [%]', 'FontSize', 15);

% Calculate max_abs_range for this subplot
max_abs_range = max(abs([gamma_power_diff; ms_rate_diff]), [], 'all');
xlim([-25 25]);
ylim([-100 100]);

% Dashed lines at x = 0 and y = 0
xline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'Alpha', 0.25);
yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'Alpha', 0.25);

title('Gamma Power vs Microsaccade Rate', 'FontSize', 20);
hold off;

% Plot Gamma Frequency difference vs Gaze Deviation difference
subplot(2, 2, 3);
hold on;
scatter(gaze_deviation_diff, gamma_freq_diff, 36, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'SizeData', 100);
xlabel('Gaze Deviation Difference [%]', 'FontSize', 15);
ylabel('Gamma Frequency Difference [%]', 'FontSize', 15);

% Calculate max_abs_range for this subplot
max_abs_range = max(abs([gamma_freq_diff; gaze_deviation_diff]), [], 'all');
max_abs_range = 45;
xlim([-max_abs_range*1.25 max_abs_range*1.25]);
ylim([-max_abs_range*1.25 max_abs_range*1.25]);

% Dashed lines at x = 0 and y = 0
xline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'Alpha', 0.25);
yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'Alpha', 0.25);

title('Gamma Frequency vs Gaze Deviation', 'FontSize', 20);
hold off;

% Plot Gamma Frequency difference vs Microsaccade Rate difference
subplot(2, 2, 4);
hold on;
scatter(ms_rate_diff, gamma_freq_diff, 36, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'SizeData', 100);
xlabel('Microsaccade Rate Difference [%]', 'FontSize', 15);
ylabel('Gamma Frequency Difference [%]', 'FontSize', 15);

% Calculate max_abs_range for this subplot
max_abs_range = max(abs([gamma_freq_diff; ms_rate_diff]), [], 'all');
max_abs_range = 45;
xlim([-max_abs_range*1.25 max_abs_range*1.25]);
ylim([-max_abs_range*1.25 max_abs_range*1.25]);

% Dashed lines at x = 0 and y = 0
xline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'Alpha', 0.25);
yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'Alpha', 0.25);

title('Gamma Frequency vs Microsaccade Rate', 'FontSize', 20);
hold off;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/stats/GCP_stats_overview_associations.png');

%% CORRELATION matrix
% close all
% % Calculate the correlation matrix between variables
% corr_matrix = corr(table2array(data(:, variables)), 'Rows', 'pairwise');
% 
% % Correlation matrix
% % Generate heatmap
% figure('Color', 'white');
% set(gcf, 'Position', [200, 400, 1200, 1500]);
% h = heatmap(corr_matrix);
% 
% % Define and add color map
% num_colors = 256;
% cpoints = [0, 0, 0, 1;      % blue at the lowest value
%            0.46, 1, 1, 1;  % transition to white starts
%            0.54, 1, 1, 1;  % transition from white ends
%            1, 1, 0, 0];    % red at the highest value
% 
% % Preallocate the colormap array.
% cmap = zeros(num_colors, 3);
% 
% % Linearly interpolate to fill in the colormap.
% for i=1:size(cmap,1)
%     % Normalize the current index to the 0-1 range based on the colormap size.
%     val = (i-1)/(num_colors-1);
%     % Find the first point where the value is greater than the current normalized value.
%     idx = find(cpoints(:,1) >= val, 1);
%     if idx == 1
%         % Use the first colour for values below the first point.
%         cmap(i,:) = cpoints(1,2:4);
%     else
%         % Linearly interpolate between the two bounding colours.
%         range = cpoints(idx,1) - cpoints(idx-1,1);
%         frac = (val - cpoints(idx-1,1)) / range;
%         cmap(i,:) = (1-frac)*cpoints(idx-1,2:4) + frac*cpoints(idx,2:4);
%     end
% end
% colormap(h, cmap);
% 
% % Customization
% h.Title = 'Correlation Heatmap';
% h.XDisplayLabels = variables; % Assuming varNames contains variable names excluding 'Date', 'Weekday', and 'Overall Score'
% h.YDisplayLabels = variables;
% h.ColorLimits = [-1, 1]; % Color limits to ensure proper color mapping
% h.FontSize = 12; % Increase the size of the x and y axis tick labels
% hTitle.FontSize = 25;
% saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/stats/GCP_stats_overview_correlation_matrix.png');
