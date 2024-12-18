%% GCP Microsaccades

% Visualizations:
%   Boxplots of microsaccade rate
%   Barplots of microsaccade rate percentage difference
%   Histograms of microsaccade timing

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load microsaccade data
ms_details_lc = {};
ms_details_hc = {};
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    cd(datapath);

    % Load ms_rate.mat
    load('ms_rate.mat');
    ms_lc(subj) = lc_msrate;
    ms_hc(subj) = hc_msrate;

    % Load ms_data.mat
    load('ms_data.mat');
    ms_details_lc{subj} = ms_data_lc;
    ms_details_hc{subj} = ms_data_hc;
end

%% Plot Microsaccades BOXPLOTS
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

%% Plot Microsaccades BARPLOT for PERCENTAGE CHANGE
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

%% Plot Microsaccade TIMING Data
close all
figure('Name', 'Microsaccade Timing Distribution', 'NumberTitle', 'off');
set(gcf, 'Position', [100, 200, 1000, 800], 'Color', 'w'); % Adjust figure size
set(gca, 'FontSize', 20);
hold on;

% Define colours for the GCP project
colors = color_def('GCP'); % Get colours: Beige and Purple

% Collect all timings for low and high conditions
ms_timings_lc = [];
ms_timings_hc = [];
for subj = 1:length(ms_details_lc)
    if ~isempty(ms_details_lc{subj})
        for trl = 1:length(ms_details_lc{subj})
            ms_timings_lc = [ms_timings_lc; ms_details_lc{subj}{trl}.Onset]; %#ok<AGROW>
        end
    end
    if ~isempty(ms_details_hc{subj})
        for trl = 1:length(ms_details_hc{subj})
            ms_timings_hc = [ms_timings_hc; ms_details_hc{subj}{trl}.Onset]; %#ok<AGROW>
        end
    end
end

% Convert to ms
ms_timings_lc = ms_timings_lc * 2;
ms_timings_hc = ms_timings_hc * 2;

% Histogram for LOW contrast
histogram(ms_timings_lc, 'BinWidth', 10, ... 
    'FaceColor', colors(1, :), 'EdgeColor', 'none', ...
    'Normalization', 'probability', 'FaceAlpha', 0.85, ...
    'DisplayName', 'Low Condition');

% Histogram for HIGH contrast
histogram(ms_timings_hc, 'BinWidth', 10, ...
    'FaceColor', colors(2, :), 'EdgeColor', 'none', ...
    'Normalization', 'probability', 'FaceAlpha', 0.5, ...
    'DisplayName', 'High Condition');


% Finalise plot
xlabel('Time [ms]');
ylabel('Density');
title('Microsaccade Timing Distribution');
legend('Location', 'best');
xlim([0 2000])
hold off;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/microsaccades/GCP_gaze_microsaccades_timing_histogram.png');

%% 3. Microsaccade Peak Velocities Distribution
close all
figure('Name', 'Microsaccade Timing Distribution', 'NumberTitle', 'off');
set(gcf, 'Position', [100, 200, 1000, 800], 'Color', 'w'); % Adjust figure size
set(gca, 'FontSize', 20);
hold on;

% Define colours for the GCP project
colors = color_def('GCP'); % Get colours: Beige and Purple

% Load data
ms_peak_lc = [];
ms_peak_hc = [];
for subj = 1:length(ms_details_lc)
    if ~isempty(ms_details_lc{subj})
        for trl = 1:length(ms_details_lc{subj})
            ms_peak_lc = [ms_peak_lc; ms_details_lc{subj}{trl}.Peak]; %#ok<AGROW>
        end
    end
    if ~isempty(ms_details_hc{subj})
        for trl = 1:length(ms_details_hc{subj})
            ms_peak_hc = [ms_peak_hc; ms_details_hc{subj}{trl}.Peak]; %#ok<AGROW>
        end
    end
end

% Convert to ms
ms_peak_lc = ms_peak_lc * 2;
ms_peak_hc = ms_peak_hc * 2;

% Histogram for LOW contrast
histogram(ms_peak_lc, 'BinWidth', 10, ... 
    'FaceColor', colors(1, :), 'EdgeColor', 'none', ...
    'Normalization', 'probability', 'FaceAlpha', 0.85, ...
    'DisplayName', 'Low Condition');

% Histogram for HIGH contrast
histogram(ms_peak_hc, 'BinWidth', 10, ...
    'FaceColor', colors(2, :), 'EdgeColor', 'none', ...
    'Normalization', 'probability', 'FaceAlpha', 0.5, ...
    'DisplayName', 'High Condition');

% Finalise plot
xlabel('Time [ms]');
ylabel('Density');
title('Microsaccade Peak Velocity Distribution');
legend('Location', 'best');
xlim([0 2000])
hold off;

% Save figure
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/microsaccades/GCP_gaze_microsaccades_peak_velocity_histogram.png');
