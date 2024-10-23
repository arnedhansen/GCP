%% GCP Gamma Peak Power and Frequency

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.2');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load data and convert TFR data to POWSCPTRM (channels x frequency)
cfg             = [];
cfg.latency     = [0 2];
cfg.frequency   = [30 120];
cfg.avgovertime = 'yes';

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath);
    load('tfr_hc.mat');
    load('tfr_lc.mat');  
    pow_hc{subj} = ft_selectdata(cfg, tfr_hc_all);
    pow_lc{subj} = ft_selectdata(cfg, tfr_lc_all);

    % Remove time dimension
    pow_lc{subj} = rmfield(pow_lc{subj}, 'time');
    pow_lc{subj}.dimord = 'chan_freq';
    pow_hc{subj} = rmfield(pow_hc{subj}, 'time');
    pow_hc{subj}.dimord = 'chan_freq';
    
    fprintf('Subject %.2s / %.3d loaded \n', num2str(subj), length(subjects))
end

% Compute grand average 
gapow_hc = ft_freqgrandaverage([], pow_hc{:});
gapow_lc = ft_freqgrandaverage([], pow_lc{:});

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, '/eeg');
cd(datapath);
% Occipital channels
occ_channels = {};
pow_label = pow_lc{1};
for i = 1:length(pow_label.label)
    label = pow_label.label{i};
    if contains(label, {'O'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot GRAND AVERAGE power spectrum
close all;
figure;
set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
colors = {'b', 'r'}; % Blue for Low Contrast, Red for High Contrast
conditions = {'Low Contrast', 'High Contrast'};
cfg = [];
cfg.channel = channels;  % Focus on occipital channels
cfg.figure = 'gcf';
cfg.linecolor = 'br';    % Line colours for both conditions
cfg.linewidth = 1;

% Plot for low and high contrast
ft_singleplotER(cfg, gapow_lc, gapow_hc);
hold on;

% Add shaded error bars
channels_seb = ismember(gapow_lc.label, cfg.channel);
lceb = shadedErrorBar(gapow_lc.freq, mean(gapow_lc.powspctrm(channels_seb, :), 1), std(gapow_lc.powspctrm(channels_seb, :))/sqrt(size(gapow_lc.powspctrm(channels_seb, :), 1)), {'b', 'markerfacecolor', 'b'});
hceb = shadedErrorBar(gapow_hc.freq, mean(gapow_hc.powspctrm(channels_seb, :), 1), std(gapow_hc.powspctrm(channels_seb, :))/sqrt(size(gapow_hc.powspctrm(channels_seb, :), 1)), {'r', 'markerfacecolor', 'r'});
transparency = 0.5;
set(lceb.patch, 'FaceAlpha', transparency);
set(hceb.patch, 'FaceAlpha', transparency);

% Adjust plot aesthetics
set(gcf,'color','w');
set(gca,'Fontsize',20);
[~, channel_idx] = ismember(channels, gapow_lc.label);
freq_idx = find(gapow_lc.freq >= 30 & gapow_lc.freq <= 90); % Adjust freq range to gamma
max_spctrm = max([mean(gapow_lc.powspctrm(channel_idx, freq_idx), 2); mean(gapow_hc.powspctrm(channel_idx, freq_idx), 2)]);
% ylim([0 max_spctrm*0.8]);
xlim([30 90])
box on;
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
title('Power Spectrum: Low vs. High Contrast', 'FontSize', 30);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_eeg_gamma_powspctrm.png');

%% Plot and save INDIVIDUAL power spectra
for subj = 1:length(subjects)
    close all;
    figure;
    set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
    
    % Extract participant data
    pow_lc_subj = pow_lc{subj};
    pow_hc_subj = pow_hc{subj};
    
    cfg = [];
    cfg.channel = channels;  % Focus on occipital channels
    cfg.figure = 'gcf';
    cfg.linecolor = 'br';    % Blue for Low Contrast, Red for High Contrast
    cfg.linewidth = 1;

    % Plot power spectrum for low and high contrast
    ft_singleplotER(cfg, pow_lc_subj, pow_hc_subj);
    hold on;
    
    % Add shaded error bars
    channels_seb = ismember(pow_lc_subj.label, cfg.channel);
    lceb = shadedErrorBar(pow_lc_subj.freq, mean(pow_lc_subj.powspctrm(channels_seb, :), 1), std(pow_lc_subj.powspctrm(channels_seb, :))/sqrt(size(pow_lc_subj.powspctrm(channels_seb, :), 1)), {'b', 'markerfacecolor', 'b'});
    hceb = shadedErrorBar(pow_hc_subj.freq, mean(pow_hc_subj.powspctrm(channels_seb, :), 1), std(pow_hc_subj.powspctrm(channels_seb, :))/sqrt(size(pow_hc_subj.powspctrm(channels_seb, :), 1)), {'r', 'markerfacecolor', 'r'});
    transparency = 0.5;
    set(lceb.patch, 'FaceAlpha', transparency);
    set(hceb.patch, 'FaceAlpha', transparency);

    % Adjust plot aesthetics
    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    [~, channel_idx] = ismember(channels, pow_lc_subj.label);
    freq_idx = find(pow_lc_subj.freq >= 30 & pow_lc_subj.freq <= 80); % Adjust freq range to gamma
    max_spctrm = max([mean(pow_lc_subj.powspctrm(channel_idx, freq_idx), 2); mean(pow_hc_subj.powspctrm(channel_idx, freq_idx), 2)]);
    % ylim([0 max_spctrm*1.25]);
    xlim([30 120]);
    box on;
    xlabel('Frequency [Hz]');
    ylabel('Power [dB]');
    legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
    title(sprintf('Subject %s: Power Spectrum', subjects{subj}), 'FontSize', 30);
    hold off;

    % Save individual plot
    saveas(gcf, sprintf('/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_eeg_gamma_powspctrm_subj%s.png', subjects{subj}));
end

%% Subplot with all INDIVIDUAL plots
figure;
set(gcf, 'Position', [0, 0, 1600, 1600], 'Color', 'w');
num_subs = length(subjects);
cols = 5;  % Number of columns
rows = 2; % Number of rows

for subj = 1:num_subs
    % Create subplot for each subject
    subplot(rows, cols, subj);
    
    % Extract participant data
    pow_lc_subj = pow_lc{subj};
    pow_hc_subj = pow_hc{subj};
    
    cfg = [];
    cfg.channel = channels;  % Focus on occipital channels
    cfg.figure = 'gcf';
    cfg.linecolor = 'br';    % Blue for Low Contrast, Red for High Contrast
    cfg.linewidth = 1;

    % Plot power spectrum for low and high contrast
    ft_singleplotER(cfg, pow_lc_subj, pow_hc_subj);
    hold on;
    
    % Add shaded error bars
    channels_seb = ismember(pow_lc_subj.label, cfg.channel);
    transparency = 0.5;

    % Adjust plot aesthetics
    set(gca, 'Fontsize', 12);
    [~, channel_idx] = ismember(channels, pow_lc_subj.label);
    freq_idx = find(pow_lc_subj.freq >= 30 & pow_lc_subj.freq <= 80); % Adjust freq range to gamma
    max_spctrm = max([mean(pow_lc_subj.powspctrm(channel_idx, freq_idx), 2); mean(pow_hc_subj.powspctrm(channel_idx, freq_idx), 2)]);
    % ylim([0 max_spctrm*1.05]);
    ylim([-4 4])
    xlim([30 120]);
    box on;
    xlabel('Freq [Hz]', 'FontSize', 10);
    ylabel('Power [dB]', 'FontSize', 10);
    title(sprintf('Subj %s', subjects{subj}), 'FontSize', 12);
    hold off;
end

% Save the combined figure with all subplots
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_eeg_gamma_powspctrm_all_subjects_subplot.png');
