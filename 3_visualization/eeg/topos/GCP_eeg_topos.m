%% GCP Gamma Topoplots

%% Setup
clear
[subjects, path] = setup('GCP');

%% Load power spectra data
for subj = 1:length(subjects)
    load(strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra'))
    pow25{subj}  = pow_c25_fooof_bl_smooth;
    pow50{subj}  = pow_c50_fooof_bl_smooth;
    pow75{subj}  = pow_c75_fooof_bl_smooth;
    pow100{subj} = pow_c100_fooof_bl_smooth;
    fprintf('Subject %3s loaded... \n', subjects{subj})
end

% Compute grand averages
gapow25  = ft_freqgrandaverage([], pow25{:});
gapow50  = ft_freqgrandaverage([], pow50{:});
gapow75  = ft_freqgrandaverage([], pow75{:});
gapow100 = ft_freqgrandaverage([], pow100{:});

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, '/eeg');
cd(datapath);
% Occipital channels
occ_channels = {};
pow_label = pow25{1, 1};
for i = 1:length(pow_label.label)
    label = pow_label.label{i};
    if contains(label, {'O'}) || contains(label, {'P'}) && ~contains(label, {'T'}) ...
        && ~contains(label, {'C'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;
%channels = [{'Pz'}, {'P1'}, {'P2'}, {'P3'}, {'P4'}, {'P5'}, {'P6'}, {'P7'}, {'P8'}, {'PPO1'}, {'PPO2'}, {'PPO5h'}, {'PPO6h'} {'PO3'}, {'PO4'}, {'PO7'}, {'PO8'}, {'POz'}, {'POO3h'}, {'POO4h'}]

%% Plot GRAND AVERAGE topoplots
close all

% Create figure
figure;
set(gcf, 'Position', [0, 0, 1200, 1000], 'Color', 'w');
set(gca, 'Fontsize', 25);
%sgtitle('Topographical Maps 300 ms - 2000 ms after Stimulus Presentation (30 - 90 Hz)', 'FontSize', 30, 'FontWeight', 'bold');
sgtitle('Grand Average Topographical Maps', 'FontSize', 30, 'FontWeight', 'bold')

% Common configuration
cfg = [];
cfg.figure = 'gcf';
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat')
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.highlight = 'on';
cfg.highlightchannel = channels;
cfg.highlightsymbol = '.';
cfg.highlightsize = 5;
cfg.marker = 'off';
cfg.comment = 'no';
cmap = flipud(cbrewer('div', 'RdBu', 64));
cfg.colormap = cmap;
cfg.gridscale = 300;
cfg.ylim = [30 90];
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'Power [dB]', 'FontSize', 25);

% Find max power value for frequency band
[~, channel_idx] = ismember(channels, gapow25.label);
freq_idx = find(gapow25.freq >= 30 & gapow25.freq <= 90);
avgscptrm25  = mean(gapow25.powspctrm(channel_idx, freq_idx), 1);
avgscptrm50  = mean(gapow50.powspctrm(channel_idx, freq_idx), 1);
avgscptrm75  = mean(gapow75.powspctrm(channel_idx, freq_idx), 1);
avgscptrm100 = mean(gapow100.powspctrm(channel_idx, freq_idx), 1);
max_spctrm = max([
    max(abs(avgscptrm25), [], 'all'), ...
    max(abs(avgscptrm50), [], 'all'), ...
    max(abs(avgscptrm75), [], 'all'), ...
    max(abs(avgscptrm100), [], 'all')]);
cfg.zlim = double([-max_spctrm * 0.9, max_spctrm * 0.9]);

% POW25
subplot(2, 2, 1);
ft_topoplotER(cfg, gapow25);
title('25% Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;
ylabel(cb, 'Power [dB]', 'FontSize', 25);

% POW50
subplot(2, 2, 2);
ft_topoplotER(cfg, gapow50);
title('50% Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;
ylabel(cb, 'Power [dB]', 'FontSize', 25);

% POW75
subplot(2, 2, 3);
ft_topoplotER(cfg, gapow75);
title('75% Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;
ylabel(cb, 'Power [dB]', 'FontSize', 25);

% POW100
subplot(2, 2, 4);
ft_topoplotER(cfg, gapow100);
title('100% Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;
ylabel(cb, 'Power [dB]', 'FontSize', 25);

% Save figure
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/topos/GCP_eeg_topos_ga.png');

%% Topoplots for Individual Subjects
close all;
% Common configuration
cfg = [];
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat')
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.highlight = 'on';
cfg.highlightchannel = channels;
cfg.highlightsymbol = '.';
cfg.highlightsize = 5;
cfg.marker = 'off';
cfg.comment = 'no';
cmap = flipud(cbrewer('div', 'RdBu', 64));
cfg.colormap = cmap;
cfg.gridscale = 300;
cfg.ylim = [30 90];

for subj = 1:length(subjects)
    close all
    
    % Prepare data for subj
    pow25_subj  = pow25{subj};
    pow50_subj  = pow50{subj};
    pow75_subj  = pow75{subj};
    pow100_subj = pow100{subj};

    % Create figure
    figure;
    set(gcf, 'Position', [0, 0, 2000, 800], 'Color', 'w');
    set(gca, 'Fontsize', 25);
    sgtitle(sprintf('Topographical Maps for Subject %s', subjects{subj}), 'FontSize', 30, 'FontWeight', 'bold');
    cfg.figure = 'gcf';

    % Find max power value for frequency band
    [~, channel_idx] = ismember(channels, pow25_subj.label);
    freq_idx = find(pow25_subj.freq >= 30 & pow25_subj.freq <= 90);
    avgscptrm25  = mean(pow25_subj.powspctrm(channel_idx, freq_idx), 1);
    avgscptrm50  = mean(pow50_subj.powspctrm(channel_idx, freq_idx), 1);
    avgscptrm75  = mean(pow75_subj.powspctrm(channel_idx, freq_idx), 1);
    avgscptrm100 = mean(pow100_subj.powspctrm(channel_idx, freq_idx), 1);
    max_spctrm = max([
        max(abs(avgscptrm25), [], 'all'), ...
        max(abs(avgscptrm50), [], 'all'), ...
        max(abs(avgscptrm75), [], 'all'), ...
        max(abs(avgscptrm100), [], 'all')]);
    cfg.zlim = double([-max_spctrm, max_spctrm]);

    % POW25
    subplot(2, 2, 1);
    ft_topoplotER(cfg, pow25_subj);
    title('25% Contrast', 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'Power [dB]', 'FontSize', 25);

    % POW50
    subplot(2, 2, 2);
    ft_topoplotER(cfg, pow50_subj);
    title('50% Contrast', 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'Power [dB]', 'FontSize', 25);

    % POW75
    subplot(2, 2, 3);
    ft_topoplotER(cfg, pow75_subj);
    title('75% Contrast', 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'Power [dB]', 'FontSize', 25);

    % POW100
    subplot(2, 2, 4);
    ft_topoplotER(cfg, pow100_subj);
    title('100% Contrast', 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'Power [dB]', 'FontSize', 25);

    % Save individual figure
    saveas(gcf, sprintf('/Volumes/methlab/Students/Arne/GCP/figures/eeg/topos/GCP_eeg_topos_subj%s.png', subjects{subj}));
end
