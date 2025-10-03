%% Gamma Time Frequency Analysis for GCP data

%% Setup
clear
[subjects, path, ~, layANThead] = setup('GCP');

%% Compute grand average time and frequency data GATFR
% Load high contrast data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load('data_tfr_superlets.mat');
    tfr25{subj}  = tfr_c25_fooof_bl_smooth;
    tfr50{subj}  = tfr_c50_fooof_bl_smooth;
    tfr75{subj}  = tfr_c75_fooof_bl_smooth;
    tfr100{subj} = tfr_c100_fooof_bl_smooth;

    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR loaded.'])
end

% Compute grand average
gatfr25  = ft_freqgrandaverage([], tfr25{:});
gatfr50  = ft_freqgrandaverage([], tfr50{:});
gatfr75  = ft_freqgrandaverage([], tfr75{:});
gatfr100 = ft_freqgrandaverage([], tfr100{:});

%% Define occipital channels
% Occipital channels
occ_channels = {};
tfrlabel = gatfr25.label;
for i = 1:length(tfrlabel)
    label = tfrlabel{i};
    if contains(label, {'O'}) || contains(label, {'P'}) && ~contains(label, {'T'}) ...
        && ~contains(label, {'C'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot GRAND AVERAGE TFR for each individual condition
close all

% Find maximum deviation across conditions
[~, channel_idx] = ismember(channels, gatfr25.label);
freq_idx = find(gatfr25.freq >= 30 & gatfr25.freq <= 90);
time_idx = find(gatfr25.time >= 0 & gatfr25.time <= 2);
avgscptrm25  = mean(gatfr25.powspctrm(channel_idx, freq_idx, time_idx), 1);
avgscptrm50  = mean(gatfr50.powspctrm(channel_idx, freq_idx, time_idx), 1);
avgscptrm75  = mean(gatfr75.powspctrm(channel_idx, freq_idx, time_idx), 1);
avgscptrm100 = mean(gatfr100.powspctrm(channel_idx, freq_idx, time_idx), 1);
max_spctrm = max([
    max(abs(avgscptrm25), [], 'all'), ...
    max(abs(avgscptrm50), [], 'all'), ...
    max(abs(avgscptrm75), [], 'all'), ...
    max(abs(avgscptrm100), [], 'all')]);
clim = double([-max_spctrm * 0.9, max_spctrm * 0.9]);

% Figure
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
sgtitle('Grand Average Time-Frequency Representations', 'FontSize', 30, 'FontWeight', 'bold')

% Common configuration
cfg = [];
cfg.figure = gcf;
cfg.layout = layANThead;
cfg.showlabels = 'yes';
cfg.channel = channels;
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = ([-0.5 2.5]);
cfg.ylim = [30 90];
%load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat')
color_map = flipud(cbrewer('div', 'RdBu', 64));

% TFR2
subplot(2, 2, 1);
ft_singleplotTFR(cfg, gatfr25);
colormap(color_map);
set(gca, 'CLim', clim);
colb = colorbar;
colb.Label.String = 'Power [dB]';
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0.3, 30, 1.7, 60], 'EdgeColor', 'r', 'LineWidth', 5);
set(gca, 'FontSize', 20)
title('25% Contrast');

% TFR50
subplot(2, 2, 2);
ft_singleplotTFR(cfg, gatfr50);
colormap(color_map);
set(gca, 'CLim', clim);
colb = colorbar;
colb.Label.String = 'Power [dB]';
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0.3, 30, 1.7, 60], 'EdgeColor', 'r', 'LineWidth', 5);
set(gca, 'FontSize', 20)
title('50% Contrast');

% TFR75
subplot(2, 2, 3);
ft_singleplotTFR(cfg, gatfr75);
colormap(color_map);
set(gca, 'CLim', clim);
colb = colorbar;
colb.Label.String = 'Power [dB]';
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0.3, 30, 1.7, 60], 'EdgeColor', 'r', 'LineWidth', 5);
set(gca, 'FontSize', 20)
title('75% Contrast');

% TFR100
subplot(2, 2, 4);
ft_singleplotTFR(cfg, gatfr100);
colormap(color_map);
set(gca, 'CLim', clim);
colb = colorbar;
colb.Label.String = 'Power [dB]';
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0.3, 30, 1.7, 60], 'EdgeColor', 'r', 'LineWidth', 5);
set(gca, 'FontSize', 20)
title('100% Contrast');

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/tfr/GCP_eeg_tfr_superlets.png');

%% Difference
diff = gatfr100;
diff.powspctrm = gatfr100.powspctrm - gatfr50.powspctrm;

figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
ft_singleplotTFR(cfg, diff);
colormap(color_map);
set(gca, 'CLim', clim);
colb = colorbar;
colb.Label.String = 'Power [dB]';
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0.3, 30, 1.7, 60], 'EdgeColor', 'r', 'LineWidth', 5);
set(gca, 'FontSize', 25)
title('TFR Difference: 100% Contrast - 50% Contrast');
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/tfr/GCP_eeg_tfr_diff_superlets.png');
