%% Heatmap for GCP gaze data

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/gaze');
    load([datapath, filesep 'dataET'])

    %% Choose data 300ms after stimulus presentation to exclude evoked activity
    cfg = [];
    cfg.latency = [0.3 2];
    dataET_lc = ft_selectdata(cfg, dataET_lc);
    dataET_hc = ft_selectdata(cfg, dataET_hc);

    %% Filter data for out-of-screen data points and zeros from blinks
    condcounter=0;
    for condition = 2:2:8
        condcounter=condcounter+1;
        if condition == 2
            data=dataetL2;
            data=horzcat(dataetL2.trial{:});
        elseif condition == 4
            data=dataetL4;
            data=horzcat(dataetL4.trial{:});
        elseif condition == 6
            data=dataetL6;
            data=horzcat(dataetL6.trial{:});
        elseif condition == 8
            data=dataetL8;
            data=horzcat(dataetL8.trial{:});
        end

        % Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600; % Check that x and y positions are in boundaries of screen
        valid_data = data(:, valid_data_indices);

        % Remove data points that contain zeros (assuming your data contains x, y, and pupil size)
        window_size = 50;
        cleaned_data = remove_blink_window(data, window_size);
        data = cleaned_data;

        x_positions = data(1, :);
        y_positions = data(2, :);

for subj = 1:length(subjects)
    %% Create custom grid for heatmap in pixels
    num_bins = 100;
    x_grid_pixels = linspace(0, 800, num_bins);
    y_grid_pixels = linspace(0, 600, num_bins);

    %% Low contrast
    x_positions = cellfun(@(x) nanmean([x{:}]), num2cell(gaze_x_lc, 1));
    y_positions = cellfun(@(x) nanmean([x{:}]), num2cell(gaze_y_lc, 1));

    % Bin data
    smoothing_factor = 5;
    binned_data_pixels = histcounts2(x_positions, y_positions, x_grid_pixels, y_grid_pixels);

    % Apply gaussian smoothing
    smoothed_data_pixels(subj, 1, :, :) = imgaussfilt(binned_data_pixels, smoothing_factor);

    % Treat ET data as TFR for stats
    freq = [];
    freq.freq       = linspace(0, 600, 99);
    freq.time       = linspace(0, 800, 99);
    freq.label      = {'et'};
    freq.dimord     = 'chan_freq_time';
    tmp(1,:,:)      = squeeze(smoothed_data_pixels(subj, 1, :, :));
    freq.powspctrm  = tmp;

    lcg{subj} = freq;

    %% High contrast
    x_positions = cellfun(@(x) nanmean([x{:}]), num2cell(gaze_x_hc, 1));
    y_positions = cellfun(@(x) nanmean([x{:}]), num2cell(gaze_y_hc, 1));

    % Bin data
    binned_data_pixels = histcounts2(x_positions, y_positions, x_grid_pixels, y_grid_pixels);

    % Apply gaussian smoothing
    smoothed_data_pixels(subj, 2, :, :) = imgaussfilt(binned_data_pixels, smoothing_factor);

    % Treat ET data as TFR for stats
    tmp(1,:,:)      = squeeze(smoothed_data_pixels(subj, 2, :, :));
    freq.powspctrm  = tmp;

    hcg{subj} = freq;
end

%% Aggregate data for subjects
subject_average = squeeze(mean(smoothed_data_pixels, 1));
lc_gdensity = subject_average(1, :, :);
hc_gdensity = subject_average(2, :, :);

%% Calculate significant differences between low and high contrast
cfg                    = [];
cfg.spmversion         = 'spm12';
cfg.method             = 'analytic';
cfg.statistic          = 'ft_statfun_depsamplesT';
cfg.tail               = 0;
cfg.clustertail        = 0;
cfg.alpha              = 0.05;
cfg.numrandomization   = 1000;
cfg.neighbours         = [];

clear design
subj = length(subjects);
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat] = ft_freqstatistics(cfg, hcg{:}, lcg{:});
stat.stat(stat.mask==0) = 0; % mask out all non-significant
statstern = stat;
cohensd = 2 * ((statstern.stat) ./ sqrt(numel(design)));
statstern.stat = cohensd;

%% Common settings for plotting
centerX = 800 / 2;
centerY = 600 / 2;

%% Plot low contrast condition
freq.powspctrm(1,:,:) = squeeze(lc_gdensity)';
freq.time = x_grid_pixels(2:end);
freq.freq = y_grid_pixels(2:end);
freq.label = {'et'};
freq.dimord = 'chan_freq_time';

clf;
close all;
figure;
mycolormap = customcolormap_preset('red-white-blue');
colormap(mycolormap);
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR([], freq);
set(gcf, 'color', 'w');
set(gca, 'Fontsize', 30);
set(gca, 'YDir', 'reverse');
title('');
clim([0 650]);
hold on; plot(centerX, centerY, 'o', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
xlim([0 800]);
ylim([0 600]);

xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
cb = colorbar; % Create the colorbar
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', 32); % Label the colorbar

saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/heatmap/GCP_gaze_lc.png');

%% Plot high contrast condition
freq.powspctrm(1,:,:) = squeeze(hc_gdensity)';
freq.time = x_grid_pixels(2:end);
freq.freq = y_grid_pixels(2:end);
freq.label = {'et'};
freq.dimord = 'chan_freq_time';

clf;
close all;
figure;
addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');
mycolormap = customcolormap_preset('red-white-blue');
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR([], freq);
set(gcf, 'color', 'w');
set(gca, 'Fontsize', 30);
set(gca, 'YDir', 'reverse');
title('');
clim([0 650]);
hold on; plot(centerX, centerY, 'o', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
xlim([0 800]);
ylim([0 600]);

xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
cb = colorbar; % Create the colorbar
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', 32); % Label the colorbar

saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/heatmap/GCP_gaze_hc.png');

%% Plot difference (hc - lc)
diff = hc_gdensity - lc_gdensity;
freq.powspctrm(1,:,:) = squeeze(diff)';
freq.time = x_grid_pixels(2:end);
freq.freq = y_grid_pixels(2:end);
freq.label = {'et'};
freq.dimord = 'chan_freq_time';

clf;
close all;
figure;
addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');
mycolormap = customcolormap_preset('red-white-blue');
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR([], freq);
set(gcf, 'color', 'w');
set(gca, 'Fontsize', 30);
set(gca, 'YDir', 'reverse');
title('');
clim([-90 90]);
hold on; plot(centerX, centerY, 'o', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
xlim([0 800]);
ylim([0 600]);

xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
cb = colorbar; % Create the colorbar
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', 32); % Label the colorbar

saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/heatmap/GCP_gaze_diff.png');

%% Plot t-value stats
mycolormap = customcolormap_preset('red-white-blue');
freq.powspctrm(1,:,:)= squeeze(stat.stat)';
freq.time = x_grid_pixels(2:end);
freq.freq = y_grid_pixels(2:end);
freq.label={'et'};
freq.dimord= 'chan_freq_time';

clf;
close all;
figure;

set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
ft_singleplotTFR([],freq);
set(gcf,'color','w');
set(gca,'Fontsize',30);
set(gca, 'YDir', 'reverse')
title('');
clim([-3.45 3.45])
hold on; plot(centerX, centerY, 'o', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
xlim([0 800]);
ylim([0 600]);

xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
cb = colorbar; % Create the colorbar
ylabel(cb, 'Effect Size [Cohen''s d]', 'FontSize', 32); % Label the colorbar

saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/heatmap/SternbergSIM_gaze_stat_diff_fine.png');
clf;
close all;

%% Plot differnces between load 2 & load 8 for all subs  - INDIVIDUAL PLOTS
for subj = 1:length(subjects)
    diff = l8g{subj};
    diff = l8g{subj}.powspctrm-l2g{subj}.powspctrm;

    freq.powspctrm(1,:,:)= squeeze(diff)';
    freq.time = x_grid_pixels(2:end);
    freq.freq = y_grid_pixels(2:end);
    freq.label={'et'};
    freq.dimord= 'chan_freq_time';

    clf;
    close all;
    figure('Color', 'w');
    addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');

    % Set colormap to have white in the middle
    mycolormap = customcolormap_preset('red-white-blue');
    totalColors = size(mycolormap, 1);
    maxVal = max(freq.powspctrm(:));
    minVal = min(freq.powspctrm(:));
    climValue = max(abs(minVal), abs(maxVal));
    proportion = 2.5 / climValue;
    rangeIndices = round(totalColors * proportion);
    middleIndex = ceil(totalColors / 2);  % find the middle index
    indicesToWhite = (middleIndex - rangeIndices):(middleIndex + rangeIndices);  % find the range of indices to set to white
    mycolormap(indicesToWhite, :) = repmat([1 1 1], length(indicesToWhite), 1);  % set them to white
    colormap(mycolormap);

    set(gcf, 'Position', [0, 0, 1000, 800]);
    cfg =[];
    cfg.figure='gcf';
    ft_singleplotTFR([],freq);
    set(gcf,'color','w');
    set(gca,'Fontsize',30);
    set(gca, 'YDir', 'reverse')
    title('');
    maxVal = max(freq.powspctrm(:));
    minVal = min(freq.powspctrm(:));
    climValue = max(abs(minVal), abs(maxVal));
    clim([-climValue climValue])
    hold on; plot(centerX, centerY, 'o', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
    xlim([0 800]);
    ylim([0 600]);

    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_stat_diff_gen_subj', num2str(subj), '.png']);

end

%% MONTECARLO

% Calculate significant differences l2 and l8
stat = [];

cfg                    = [];
cfg.spmversion         = 'spm12';
% cfg.method             = 'analytic';
cfg.method             = 'montecarlo';
cfg.statistic          = 'ft_statfun_depsamplesT';
cfg.tail               = 0;
cfg.clustertail        = 0;
cfg.alpha              = 0.05;
cfg.numrandomization   = 'all';
cfg.neighbours         = [];

clear design
subj = length(subjects);
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat] = ft_freqstatistics(cfg, l8g{:}, l2g{:});
stat.stat(stat.mask==0)=0;% mask out all non significant
statstern=stat;
cohensd=2*((statstern.stat)./sqrt(numel(design)));
statstern.stat=cohensd;

% Plot t-value stats
mycolormap = customcolormap_preset('red-white-blue');
freq.powspctrm(1,:,:)= squeeze(stat.stat)';
freq.time = x_grid_pixels(2:end);
freq.freq = y_grid_pixels(2:end);
freq.label={'et'};
freq.dimord= 'chan_freq_time';

clf;
close all;
figure;

set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
ft_singleplotTFR([],freq);
set(gcf,'color','w');
set(gca,'Fontsize',30);
set(gca, 'YDir', 'reverse')
title('');
clim([-3.45 3.45])
hold on; plot(centerX, centerY, 'o', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
xlim([0 800]);
ylim([0 600]);

saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/gaze/heatmap/SternbergSIM_gaze_stat_diff_fine_montecarlo.png');