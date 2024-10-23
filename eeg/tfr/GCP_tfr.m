%% Gamma Time Frequency Analysis for GCP data
clear
clc
close all
run startup
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Compute grand average time and frequency data GATFR
% Load high contrast data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_hc
    tfrhc_0{subj} = tfr_hc_0;
    tfrhc_45{subj} = tfr_hc_45;
    tfrhc_90{subj} = tfr_hc_90;
    tfrhc_115{subj} = tfr_hc_115;
    tfrhc_all{subj} = tfr_hc_all;
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR HIGH CONTRAST loaded.'])
end

% Load low contrast data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_lc
    tfrlc_0{subj} = tfr_lc_0;
    tfrlc_45{subj} = tfr_lc_45;
    tfrlc_90{subj} = tfr_lc_90;
    tfrlc_115{subj} = tfr_lc_115;
    tfrlc_all{subj} = tfr_lc_all;
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR LOW CONTRAST loaded.'])
end

% Compute grand average
ga_tfr_hc_0 = ft_freqgrandaverage([],tfrhc_0{:});
ga_tfr_hc_45 = ft_freqgrandaverage([],tfrhc_45{:});
ga_tfr_hc_90 = ft_freqgrandaverage([],tfrhc_90{:});
ga_tfr_hc_115 = ft_freqgrandaverage([],tfrhc_115{:});
ga_tfr_hc_all = ft_freqgrandaverage([],tfrhc_all{:});

ga_tfr_lc_0 = ft_freqgrandaverage([],tfrlc_0{:});
ga_tfr_lc_45 = ft_freqgrandaverage([],tfrlc_45{:});
ga_tfr_lc_90 = ft_freqgrandaverage([],tfrlc_90{:});
ga_tfr_lc_115 = ft_freqgrandaverage([],tfrlc_115{:});
ga_tfr_lc_all = ft_freqgrandaverage([],tfrlc_all{:});

%% Define occipital channels
load('tfr_hc.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(tfr_lc_0.label)
    label = tfr_lc_0.label{i};
    if contains(label, {'O'}) 
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot TFR for HIGH and LOW CONTRAT conditions
close all

% Common configuration
cfg = [];
cfg.layout = layANThead; % your specific layout
cfg.showlabels = 'yes'; % show channel labels
cfg.channel = channels; % specify the channels to include
cfg.colorbar = 'yes'; % include color bar
cfg.zlim = 'maxabs'; % color limits
cfg.xlim = ([0 2]);
cfg.ylim = [30 90];
% clim = ([-5 5]);
color_map = flipud(cbrewer('div', 'RdBu', 64)); % 'RdBu' for blue to red diverging color map

% HIGH CONTRAST ALL
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

ft_singleplotTFR(cfg, ga_tfr_hc_all);
colormap(color_map);
% set(gca, 'CLim', clim);
set(gca, 'FontSize', 20)
colorbar;
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
title('High Contrast (ALL Orientations)');

saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/tfr/GCP_tfr_hc_all.png');

% LOW CONTRAST ALL
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

ft_singleplotTFR(cfg, ga_tfr_lc_all);
colormap(color_map);
% set(gca, 'CLim', clim);
set(gca, 'FontSize', 20)
colorbar;
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
title('Low Contrast (ALL Orientations)');

saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/tfr/GCP_tfr_lc_all.png');


%% Plot the grand averages for the DIFFERENCE 
close all

% Compute difference
diff = ga_tfr_hc_all;
diff.powspctrm = ga_tfr_hc_all.powspctrm - ga_tfr_lc_all.powspctrm;

% Define configuration for multiplot
cfg = [];
cfg.layout = layANThead; % your specific layout
cfg.channel = channels; % specify the channels to include
cfg.showlabels = 'yes'; % show channel labels
cfg.colorbar = 'yes'; % include color bar
cfg.xlim = [0 2]; 
cfg.ylim = [30 90];
color_map = flipud(cbrewer('div', 'RdBu', 64)); % 'RdBu' for blue to red diverging color map
clim = [-0.6, 0.61];
% cliom = max(diff.powspctrm(80:120, 5:10, 81))

% Plot: Difference Time-Frequency Response
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

ft_singleplotTFR(cfg, diff);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
colorbar;
colormap(color_map);
set(gca, 'FontSize', 25);
title('TFR Diff (High Contrast - Low Contrast)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'CLim', clim);
yline(30, 'LineWidth', 5, 'LineStyle', '-', 'Color', 'r');
yline(90, 'LineWidth', 5, 'LineStyle', '-', 'Color', 'r');

% Save 
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/tfr/GCP_tfr_diff.png');
