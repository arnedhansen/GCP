%% Gamma Time Frequency Analysis for GCP data

%% Setup
[subjects, path] = setup('GCP');

%% Compute grand average time and frequency data GATFR
% Load high contrast data
for subj = 1%:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load('data_tfr.mat');

    % Horizontal
    tfr_horz_lc_baselined{subj}                             = tfr_horz_lc_bl;
    tfr_horz_hc_baselined{subj}                             = tfr_horz_hc_bl;
    
    % Vertical
    tfr_vert_lc_baselined{subj}                             = tfr_vert_lc_bl;
    tfr_vert_hc_baselined{subj}                             = tfr_vert_hc_bl;
    
    % Concentric Static
    tfr_concentric_static_lc_baselined{subj}                = tfr_concentric_static_lc_bl;
    tfr_concentric_static_hc_baselined{subj}                = tfr_concentric_static_hc_bl;
   
    % Concentric Dynamic Inward
    tfr_concentric_dynamic_inward_lc_baselined{subj}        = tfr_concentric_dynamic_inward_lc_bl;
    tfr_concentric_dynamic_inward_hc_baselined{subj}        = tfr_concentric_dynamic_inward_hc_bl;
    
    % Concentric Dynamic Outward
    tfr_concentric_dynamic_outward_lc_baselined{subj}       = tfr_concentric_dynamic_outward_lc_bl;
    tfr_concentric_dynamic_outward_hc_baselined{subj}       = tfr_concentric_dynamic_outward_hc_bl;

    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR loaded.'])
end

%% Compute grand average
% Horizontal
gatfr_horz_lc_baselined = ft_freqgrandaverage([], tfr_horz_lc_baselined{:});
gatfr_horz_hc_baselined = ft_freqgrandaverage([], tfr_horz_hc_baselined{:});

% Vertical
gatfr_vert_lc_baselined= ft_freqgrandaverage([], tfr_vert_lc_baselined{:});
gatfr_vert_hc_baselined= ft_freqgrandaverage([], tfr_vert_hc_baselined{:});

% Concentric Static
gatfr_concentric_static_lc_baselined = ft_freqgrandaverage([], tfr_concentric_static_lc_baselined{:});
gatfr_concentric_static_hc_baselined = ft_freqgrandaverage([], tfr_concentric_static_hc_baselined{:});

% Concentric Dynamic Inward
gatfr_concentric_dynamic_inward_lc_baselined = ft_freqgrandaverage([], tfr_concentric_dynamic_inward_lc_baselined{:});
gatfr_concentric_dynamic_inward_hc_baselined = ft_freqgrandaverage([], tfr_concentric_dynamic_inward_hc_baselined{:});

% Concentric Dynamic Outward
gatfr_concentric_dynamic_outward_lc_baselined = ft_freqgrandaverage([], tfr_concentric_dynamic_outward_lc_baselined{:});
gatfr_concentric_dynamic_outward_hc_baselined = ft_freqgrandaverage([], tfr_concentric_dynamic_outward_hc_baselined{:});

%% Define occipital channels
% Occipital channels
occ_channels = {};
tfrlabel = tfr_horz_lc_bl.label;
for i = 1:length(tfrlabel)
    label = tfrlabel{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot TFR for HIGH and LOW CONTRAST conditions
close all
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat')
gatfr_data = {gatfr_horz_lc_baselined, gatfr_horz_hc_baselined; ...
    gatfr_vert_lc_baselined, gatfr_vert_hc_baselined; ...
    gatfr_concentric_static_lc_baselined, gatfr_concentric_static_hc_baselined; ...
    gatfr_concentric_dynamic_inward_lc_baselined, gatfr_concentric_dynamic_inward_hc_baselined; ...
    gatfr_concentric_dynamic_outward_lc_baselined, gatfr_concentric_dynamic_outward_hc_baselined};

% List of conditions
conditions_list = {'horizontal', 'vertical', 'concentric_static', 'concentric_dynamic_inward', 'concentric_dynamic_outward'};

% File path for saving plots
output_path = '/Volumes/methlab/Students/Arne/GCP/figures/eeg/tfr/';

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

% Loop through all conditions
for cond_idx = 1:length(conditions_list)
    condition = conditions_list{cond_idx};
    gatfr_lc = gatfr_data{cond_idx, 1};
    gatfr_hc = gatfr_data{cond_idx, 2};

    % LOW CONTRAST
    figure;
    set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

    ft_singleplotTFR(cfg, gatfr_lc);
    colormap(color_map);
    % set(gca, 'CLim', clim);
    set(gca, 'FontSize', 20)
    colb = colorbar;
    colb.Label.String = 'Power [dB]';
    c.Label.FontSize = 25;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title(['TFR LOW CONTRAST ', upper(condition)], 'FontSize', 30);

    % Save
    saveas(gcf, strcat(output_path, 'GCP_eeg_tfr_', condition, '_lc.png'));

    % HIGH CONTRAST
    figure;
    set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

    ft_singleplotTFR(cfg, gatfr_hc);
    colormap(color_map);
    % set(gca, 'CLim', clim);
    set(gca, 'FontSize', 20)
    colb = colorbar;
    colb.Label.String = 'Power [dB]';
    c.Label.FontSize = 25;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title(['TFR HIGH CONTRAST ', upper(condition)], 'FontSize', 30);

    % Save
    saveas(gcf, strcat(output_path, 'GCP_eeg_tfr_', condition, '_hc.png'));

    % Difference
    diff = gatfr_hc;
    diff.powspctrm = gatfr_hc.powspctrm - gatfr_lc.powspctrm;

    figure;
    set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

    ft_singleplotTFR(cfg, diff);
    colormap(color_map);
    % set(gca, 'CLim', clim);
    set(gca, 'FontSize', 20)
    colb = colorbar;
    colb.Label.String = 'Power [dB]';
    c.Label.FontSize = 25;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title(['TFR DIFFERENCE (HC-LC) ', upper(condition)], 'FontSize', 30);

    % Save
    saveas(gcf, strcat(output_path, 'GCP_eeg_tfr_', condition, '_diff.png'));

end

%% Plot the grand averages for the DIFFERENCE
close all

% Compute difference
diff = gatfr_hc_all;
diff.powspctrm = gatfr_hc_all.powspctrm - gatfr_lc_all.powspctrm;

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

% Plot: Difference Time-Frequency Response
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

ft_singleplotTFR(cfg, diff);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
colb = colorbar;
colb.Label.String = 'Power [dB]';
c.Label.FontSize = 25;
colormap(color_map);
set(gca, 'FontSize', 25);
title('Time-Frequency Representation: Difference (High Contrast - Low Contrast)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'CLim', clim);
%yline(30, 'LineWidth', 5, 'LineStyle', '-', 'Color', 'r');
%yline(90, 'LineWidth', 5, 'LineStyle', '-', 'Color', 'r');

% Save 
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/tfr/GCP_tfr_diff.png');
