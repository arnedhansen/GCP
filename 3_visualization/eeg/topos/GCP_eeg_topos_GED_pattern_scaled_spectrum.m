%% GCP GED pattern scaled by condition spectrum
%
% This legacy visualization multiplies each subject's fixed full-window GED
% activation pattern by the condition-specific GED power spectrum. Spatial
% differences between conditions therefore reflect amplitude scaling only.

%% Setup
clear
startup
[subjects, paths] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);

ged = load(fullfile(paths.features, 'GCP_eeg_GED.mat'), ...
    'subjects', 'all_topos', 'all_topo_labels', ...
    'all_condition_powspctrm_full', 'scan_freqs');
[~, subjectIndex] = ismember(subjects, ged.subjects);

%% Reconstruct channel by frequency GED data
pow25 = cell(1, numel(subjects));
pow50 = cell(1, numel(subjects));
pow75 = cell(1, numel(subjects));
pow100 = cell(1, numel(subjects));
conditionPower = {pow25, pow50, pow75, pow100};

for subj = 1:numel(subjects)
    sourceIndex = subjectIndex(subj);
    topo = double(ged.all_topos{sourceIndex}(:));
    topo = topo ./ sqrt(mean(topo .^ 2));

    for cond = 1:4
        freqData = [];
        freqData.label = ged.all_topo_labels{sourceIndex}(:);
        freqData.freq = ged.scan_freqs(:)';
        freqData.powspctrm = topo * ...
            double(ged.all_condition_powspctrm_full{cond, sourceIndex}(:)');
        freqData.dimord = 'chan_freq';
        conditionPower{cond}{subj} = freqData;
    end
end

pow25 = conditionPower{1};
pow50 = conditionPower{2};
pow75 = conditionPower{3};
pow100 = conditionPower{4};

%% Compute grand averages
gapow25 = ft_freqgrandaverage([], pow25{:});
gapow50 = ft_freqgrandaverage([], pow50{:});
gapow75 = ft_freqgrandaverage([], pow75{:});
gapow100 = ft_freqgrandaverage([], pow100{:});

%% Define channels
occChannels = {};
for iChannel = 1:numel(gapow25.label)
    label = gapow25.label{iChannel};
    if contains(label, 'O') || ...
            (contains(label, 'P') && ~contains(label, 'T') && ...
            ~contains(label, 'C')) || contains(label, 'I')
        occChannels{end + 1} = label; %#ok<SAGROW>
    end
end

%% Plot grand average topoplots
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
set(gca, 'FontSize', 25);
sgtitle('GED Pattern Scaled by Condition Spectrum', ...
    'FontSize', 30, 'FontWeight', 'bold');

cfg = [];
load(fullfile(paths.base_students, 'toolboxes', ...
    'headmodel', 'layANThead.mat'), 'layANThead')
cfg.figure = 'gcf';
cfg.layout = layANThead;
allChannels = cfg.layout.label;
cfg.channel = allChannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.highlight = 'on';
cfg.highlightchannel = occChannels;
cfg.highlightsymbol = '.';
cfg.highlightsize = 7.5;
cfg.marker = 'off';
cfg.comment = 'no';
cfg.colormap = flipud(cbrewer('div', 'RdBu', 64));
cfg.gridscale = 300;
cfg.ylim = [30 90];

[~, channelIndex] = ismember(occChannels, gapow25.label);
frequencyIndex = gapow25.freq >= 30 & gapow25.freq <= 90;
avgSpectrum25 = mean(gapow25.powspctrm(channelIndex, frequencyIndex), 1);
avgSpectrum50 = mean(gapow50.powspctrm(channelIndex, frequencyIndex), 1);
avgSpectrum75 = mean(gapow75.powspctrm(channelIndex, frequencyIndex), 1);
avgSpectrum100 = mean(gapow100.powspctrm(channelIndex, frequencyIndex), 1);
maxSpectrum = max(abs([avgSpectrum25, avgSpectrum50, ...
    avgSpectrum75, avgSpectrum100]), [], 'all');
cfg.zlim = double([-maxSpectrum * 0.9, maxSpectrum * 0.9]);

subplot(2, 2, 1);
ft_topoplotER(cfg, gapow25);
title('25% Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;
ylabel(cb, 'Pattern-scaled power [a.u.]', 'FontSize', 25);

subplot(2, 2, 2);
ft_topoplotER(cfg, gapow50);
title('50% Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;
ylabel(cb, 'Pattern-scaled power [a.u.]', 'FontSize', 25);

subplot(2, 2, 3);
ft_topoplotER(cfg, gapow75);
title('75% Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;
ylabel(cb, 'Pattern-scaled power [a.u.]', 'FontSize', 25);

subplot(2, 2, 4);
ft_topoplotER(cfg, gapow100);
title('100% Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;
ylabel(cb, 'Pattern-scaled power [a.u.]', 'FontSize', 25);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(paths.figures, 'eeg', 'topos', ...
    'GCP_eeg_topos_GED_pattern_scaled_spectrum.png'), '-dpng', '-r600');
