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
    for condition = {'lc', 'hc'}
        if strcmp(condition, 'lc')
            data = dataET_lc;
            data = horzcat(dataET_lc.trial{:});
        elseif strcmp(condition, 'hc')
            data = dataET_hc;
            data = horzcat(dataET_hc.trial{:});
        end

        % Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(1:3, valid_data_indices); % Excluding pupil size data

        % Remove blinks with a window of 100ms (= 50 samples/timepoints)
        win_size = 50;
        data = remove_blinks(data, win_size);

        if strcmp(condition, 'lc')
            x_positions_lc = data(1, :);
            y_positions_lc = data(2, :);
        elseif strcmp(condition, 'hc')
            x_positions_hc = data(1, :);
            y_positions_hc = data(2, :);
        end
    end

    %% Bin and smooth data
    % Create custom grid for heatmap in pixels
    num_bins = 100;
    x_grid_pixels = linspace(0, 800, num_bins);
    y_grid_pixels = linspace(0, 600, num_bins);

    % Bin data
    smoothing_factor = 5;
    binned_data_pixels_lc = histcounts2(x_positions_lc, y_positions_lc, x_grid_pixels, y_grid_pixels);
    binned_data_pixels_hc = histcounts2(x_positions_hc, y_positions_hc, x_grid_pixels, y_grid_pixels);

    % Apply gaussian smoothing
    smoothed_data_pixels_lc(subj, :, :) = imgaussfilt(binned_data_pixels_lc, smoothing_factor);
    smoothed_data_pixels_hc(subj, :, :) = imgaussfilt(binned_data_pixels_hc, smoothing_factor);

    % Treat ET data as TFR for stats
    freq = [];
    freq.freq       = linspace(0, 600, 99);
    freq.time       = linspace(0, 800, 99);
    freq.label      = {'et'};
    freq.dimord     = 'chan_freq_time';
    tmp(1,:,:)      = squeeze(smoothed_data_pixels_lc(subj, :, :));
    freq.powspctrm  = tmp;

    cfg = [];
    cfg.frequency = [0 600];
    lcg{subj} = freq;
    hcg{subj} = freq;    
end

%% Aggregate data for subjects
lc_gdensity = squeeze(mean(smoothed_data_pixels_lc, 1));
%lc_gdensity = subject_average_lc(1, :, :);
hc_gdensity = squeeze(mean(smoothed_data_pixels_hc, 1));
%hc_gdensity = subject_average_hc(2, :, :);

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
title('');
%clim([0 650]);
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
title('');
%clim([0 650]);
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