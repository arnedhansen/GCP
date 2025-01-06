%% GCP 3d scatter plot for gaze metrics
%   gaze deviation
%   saccades
%   microsaccades

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

%% Load gaze data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    cd(datapath);
    
    % Load gaze deviation (euclidean distances) data
    load('gaze_dev.mat');
    devs_lc(subj) = lc_gdev;
    devs_hc(subj) = hc_gdev;
    
    % Load gaze standard deviation data
    load('gaze_std.mat');
    stdx_lc(subj) = lc_gSDx;
    stdy_lc(subj) = lc_gSDy;
    stdx_hc(subj) = hc_gSDx;
    stdy_hc(subj) = hc_gSDy;

    % Load microsaccade rate
    load('ms_rate.mat');
    ms_lc(subj) = lc_msrate;
    ms_hc(subj) = hc_msrate;
    
    % Load saccade data
    load('gaze_metrics.mat'); % Assuming saccade data is in this file
    saccades_lc(subj) = lc_saccades;
    saccades_hc(subj) = hc_saccades;
end

%% Aggregate measures
gaze_dev = [devs_lc, devs_hc];
gazeSDx = [stdx_lc, stdx_hc];
gazeSDy = [stdy_lc, stdy_hc];
microsaccades = [ms_lc, ms_hc];
saccades = [saccades_lc, saccades_hc];

%% Relation of gaze standard deviation with gaze deviation from euclidean distances
close all
figure; 
set(gcf, 'Position', [0, 0, 800, 600], 'Color', 'w');
xlabel('Gaze Standard Deviation [px]');
ylabel('Gaze Deviation [px]');
title('Scatter of gaze standard deviation with gaze deviation from euclidean distances');
% scatter(stdy_lc(stdy_lc<10)+stdx_lc(stdy_lc<10), devs_lc(stdy_lc<10));
% scatter(stdy_lc(stdy_lc<10), devs_lc(stdy_lc<10));
scatter(stdy_lc, devs_lc);
lsline;

%% Create 3D Scatter Plot with Microsaccade Plane
conditions = [ones(1, length(subjects)), 2 * ones(1, length(subjects))];
close all
figure;
set(gcf, 'Position', [0, 0, 2000, 1200], 'Color', 'w');
scatter3(gaze_dev, saccades, microsaccades, 250, conditions, 'filled');
xlabel('Gaze Deviation [px]');
ylabel('Saccades');
zlabel('Microsaccades [ms/s]');
title('3D Scatter Plot of Gaze Metrics');
colors = color_def('GCP');
colormap([colors(1, :); colors(2, :)]); 
grid on;
ax = gca;
ax.FontSize = 16;

%% Gaze STD vs. Gaze Deviation
close all
figure;
colors = color_def('GCP'); 
colormap([colors(1, :); colors(2, :)]);
marker_size = 150; 
font_size = 18;
set(gcf, 'Position', [0, 0, 1200, 800], 'Color', 'w');
scatter(gazeSDx + gazeSDy, gaze_dev, marker_size, conditions, 'filled');
lsline;
title('Gaze STD vs Gaze Deviation', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('Gaze SD [px]', 'FontSize', font_size);
ylabel('Gaze Deviation [px]', 'FontSize', font_size);
ax = gca;
ax.FontSize = font_size;
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/relations/GCP_gaze_metrics_relations_STDvsDEV.png');

%% Subplots for gaze deviation from euclidean distances
close all
figure;
set(gcf, 'Position', [0, 0, 1600, 900], 'Color', 'w');
colormap([colors(1, :); colors(2, :)]);

% Subplot 1: Microsaccades vs Gaze Deviation
subplot(1, 3, 1);
scatter(microsaccades, gaze_dev, marker_size, conditions, 'filled');
lsline;
title('Microsaccades vs Gaze Deviation', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('Microsaccades [ms/s]', 'FontSize', font_size);
ylabel('Gaze Deviation [px]', 'FontSize', font_size);
ax = gca;
ax.FontSize = font_size;

% Subplot 2: Microsaccades vs Saccades
subplot(1, 3, 2);
scatter(microsaccades, saccades, marker_size, conditions, 'filled');
lsline;
title('Microsaccades vs Saccades', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('Microsaccades [ms/s]', 'FontSize', font_size);
ylabel('Saccades [count]', 'FontSize', font_size);
ax = gca;
ax.FontSize = font_size;

% Subplot 3: Saccades vs Gaze Deviation
subplot(1, 3, 3);
scatter(saccades, gaze_dev, marker_size, conditions, 'filled');
lsline;
title('Saccades vs Gaze Deviation', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('Saccades [count]', 'FontSize', font_size);
ylabel('Gaze Deviation [px]', 'FontSize', font_size);
ax = gca;
ax.FontSize = font_size;

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/relations/GCP_gaze_metrics_relations_DEV.png');

%% Subplots for gaze standard deviation
close all
figure;
set(gcf, 'Position', [0, 0, 1600, 900], 'Color', 'w');
colormap([colors(1, :); colors(2, :)]);

% Subplot 1: Microsaccades vs Gaze STD
subplot(1, 3, 1);
scatter(microsaccades, gazeSDx + gazeSDy, marker_size, conditions, 'filled');
lsline;
title('Microsaccades vs Gaze STD', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('Microsaccades [ms/s]', 'FontSize', font_size);
ylabel('Gaze Deviation SD [px]', 'FontSize', font_size);
ax = gca;
ax.FontSize = font_size;

% Subplot 2: Microsaccades vs Saccades
subplot(1, 3, 2);
scatter(microsaccades, saccades, marker_size, conditions, 'filled');
lsline;
title('Microsaccades vs Saccades', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('Microsaccades [ms/s]', 'FontSize', font_size);
ylabel('Saccades [count]', 'FontSize', font_size);
ax = gca;
ax.FontSize = font_size;

% Subplot 3: Saccades vs Gaze STD
subplot(1, 3, 3);
scatter(saccades, gazeSDx + gazeSDy, marker_size, conditions, 'filled');
lsline;
title('Saccades vs Gaze STD', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('Saccades [count]', 'FontSize', font_size);
ylabel('Gaze Deviation [px]', 'FontSize', font_size);
ax = gca;
ax.FontSize = font_size;
hold off

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/relations/GCP_gaze_metrics_relations_STD.png');

%% mit gaze std machen
