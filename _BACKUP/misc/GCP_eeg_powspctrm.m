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

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, '/eeg');
cd(datapath);
load('power_lc.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload_lc_0.label)
    label = powload_lc_0.label{i};
    if contains(label, {'O'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data 
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath);
    load('power_lc.mat');
    load('power_hc.mat');
    
    pow_hc_0{subj} = powload_hc_0;
    pow_hc_45{subj} = powload_hc_45;
    pow_hc_90{subj} = powload_hc_90;
    pow_hc_115{subj} = powload_hc_115;
    
    pow_lc_0{subj} = powload_lc_0;
    pow_lc_45{subj} = powload_lc_45;
    pow_lc_90{subj} = powload_lc_90;
    pow_lc_115{subj} = powload_lc_115;
    
    fprintf('Subject %.2s / %.3d loaded \n', num2str(subj), length(subjects))
end
plot(gapow_lc.powspctrm)

% Compute grand average for low contrast condition
gapow_lc_0 = ft_freqgrandaverage([], pow_lc_0{:});
gapow_lc_45 = ft_freqgrandaverage([], pow_lc_45{:});
gapow_lc_90 = ft_freqgrandaverage([], pow_lc_90{:});
gapow_lc_115 = ft_freqgrandaverage([], pow_lc_115{:});

% Combine into a single structure
gapow_lc = ft_freqgrandaverage([], gapow_lc_0, gapow_lc_45, gapow_lc_90, gapow_lc_115);

% Compute grand average for high contrast condition
gapow_hc_0 = ft_freqgrandaverage([], pow_hc_0{:});
gapow_hc_45 = ft_freqgrandaverage([], pow_hc_45{:});
gapow_hc_90 = ft_freqgrandaverage([], pow_hc_90{:});
gapow_hc_115 = ft_freqgrandaverage([], pow_hc_115{:});

% Combine into a single structure
gapow_hc = ft_freqgrandaverage([], gapow_hc_0, gapow_hc_45, gapow_hc_90, gapow_hc_115);

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

% Log-transform
gapow_lc.powspctrm = log10(gapow_lc.powspctrm);
gapow_hc.powspctrm = log10(gapow_hc.powspctrm);

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
freq_idx = find(gapow_lc.freq >= 30 & gapow_lc.freq <= 80); % Adjust freq range to gamma
max_spctrm = max([mean(gapow_lc.powspctrm(channel_idx, freq_idx), 2); mean(gapow_hc.powspctrm(channel_idx, freq_idx), 2)]);
ylim([0 max_spctrm*1.25]);
xlim([30 90])
box on;
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
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
    pow_lc_subj = pow_lc_0{subj};
    pow_hc_subj = pow_hc_0{subj};
    
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
    ylim([0 max_spctrm*1.25]);
    xlim([30 90]);
    box on;
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
    title(sprintf('Subject %s: Power Spectrum: Low vs. High Contrast', subjects{subj}), 'FontSize', 30);
    hold off;

    % Save individual plot
    saveas(gcf, sprintf('/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_eeg_gamma_powspctrm_subj%s.png', subjects{subj}));
end
