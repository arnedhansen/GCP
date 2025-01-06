%% GCP 3d scatter plot for gaze metrics
%   gaze deviation
%   saccades
%   microsaccades


%% als GA
figure; scatter(gazeSDy(gazeSDy<10)+gazeSDx(gazeSDy<10), gazeDev(gazeSDy<10));lsline



%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Initialise variables for metrics
devs_lc = zeros(1, length(subjects));
devs_hc = zeros(1, length(subjects));
ms_lc = zeros(1, length(subjects));
ms_hc = zeros(1, length(subjects));
saccades_lc = zeros(1, length(subjects));
saccades_hc = zeros(1, length(subjects));

%% Load gaze deviation and microsaccade data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    cd(datapath);
    
    % Load gaze deviation data
    load('gaze_dev.mat');
    devs_lc(subj) = lc_gdev;
    devs_hc(subj) = hc_gdev;
    
    % Load microsaccade rate
    load('ms_rate.mat');
    ms_lc(subj) = lc_msrate;
    ms_hc(subj) = hc_msrate;
    
    % Load saccade data
    load('gaze_metrics.mat'); % Assuming saccade data is in this file
    saccades_lc(subj) = lc_saccades;
    saccades_hc(subj) = hc_saccades;
end

%% Combine metrics for plotting
gaze_dev = [devs_lc, devs_hc];
microsaccades = [ms_lc, ms_hc];
saccades = [saccades_lc, saccades_hc];

% Create labels for conditions (Low Contrast = 1, High Contrast = 2)
conditions = [ones(1, length(subjects)), 2 * ones(1, length(subjects))];

%% Create 3D Scatter Plot with Microsaccade Plane
figure;
set(gcf, 'Position', [0, 0, 2000, 1200], 'Color', 'w');
scatter3(gaze_dev, saccades, microsaccades, 250, conditions, 'filled');
xlabel('Gaze Deviation [px]');
ylabel('Saccades');
zlabel('Microsaccades [ms/s]');
title('3D Scatter Plot of Gaze Metrics with Microsaccade Plane');
colors = color_def('GCP');
colormap([colors(1, :); colors(2, :)]); % Colour map for conditions

% Define the plane for microsaccades
[X, Y] = meshgrid(linspace(min(gaze_dev), max(gaze_dev), 50), ...
                  linspace(min(saccades), max(saccades), 50));
Z = mean(microsaccades) * ones(size(X)); % Example: mean microsaccade value

% Add the plane
hold on;
surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'cyan');
grid on;

legend({'Data Points', 'Microsaccade Plane'}, 'Location', 'northeastoutside');

%%
close all
figure;
subplot(2, 3, 1)
scatter(microsaccades, gaze_dev);lsline
subplot(2, 3, 2)
scatter(microsaccades, saccades);lsline
subplot(2, 3, 3)
scatter(saccades, gaze_dev);lsline
subplot(2, 3, 4)
scatter(microsaccades, gaze_devSDx+SDy);lsline % als GA
subplot(2, 3, 5)
scatter(microsaccades, saccades);lsline
subplot(2, 3, 6)
scatter(gazeSD+++, gaze_dev);lsline


%% mit gaze std machen
