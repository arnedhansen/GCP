%% Gamma Time Frequency Analysis for GCP data

%% Setup
[subjects, path] = setup('GCP');

%% Compute grand average time and frequency data GATFR
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load('data_tfr.mat');
    tfr_lc_baselined{subj}                             = tfr_lc_bl;
    tfr_hc_baselined{subj}                             = tfr_hc_bl;

    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR loaded.'])
end

%% Compute grand average
gatfr_lc_baselined = ft_freqgrandaverage([], tfr_lc_baselined{:});
gatfr_hc_baselined = ft_freqgrandaverage([], tfr_hc_baselined{:});

%% Define occipital channels
% Occipital channels
occ_channels = {};
tfrlabel = tfr_lc_bl.label;
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

% LOW CONTRAST
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

ft_singleplotTFR(cfg, gatfr_lc_baselined);
colormap(color_map);
% set(gca, 'CLim', clim);
set(gca, 'FontSize', 20)
colb = colorbar;
colb.Label.String = 'Power [dB]';
c.Label.FontSize = 25;
xlabel('Time [s]');
ylabel('Frequency [Hz]');
title('TFR LOW CONTRAST ', 'FontSize', 30);

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/tfr/GCP_eeg_tfr_lc.png');

% HIGH CONTRAST
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

ft_singleplotTFR(cfg, gatfr_hc_baselined);
colormap(color_map);
% set(gca, 'CLim', clim);
set(gca, 'FontSize', 20)
colb = colorbar;
colb.Label.String = 'Power [dB]';
c.Label.FontSize = 25;
xlabel('Time [s]');
ylabel('Frequency [Hz]');
title('TFR HIGH CONTRAST ', 'FontSize', 30);

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/tfr/GCP_eeg_tfr_hc.png');

%% Plot TFR for the DIFFERENCE
diff = gatfr_hc_baselined;
diff.powspctrm = gatfr_hc_baselined.powspctrm - gatfr_lc_baselined.powspctrm;

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
title('TFR DIFFERENCE (HC-LC)', 'FontSize', 30);

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/tfr/GCP_eeg_tfr_diff.png');