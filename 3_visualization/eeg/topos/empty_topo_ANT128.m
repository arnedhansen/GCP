%% Empty ANT128 Neuro Topoplot for Channel Locations
close all
[~, paths] = setup('GCP', 0);
figure;
set(gcf, 'Position', [100, 100, 1000, 800], 'Color', 'w');
title('Electrodes ANT 128', 'FontSize', 25);
load(fullfile(paths.base_students, 'toolboxes', 'headmodel', 'layANThead.mat'));
load(fullfile(paths.features, '601', 'eeg', 'power_spectra'));
channels = pow_c25_fooof_bl_smooth.label;
pow_c25_fooof_bl_smooth.powspctrm(:, :) = 0;
cfg = [];
cfg.figure = 'gcf';
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.marker = 'labels';
cfg.comment = 'no';
cmap = flipud(cbrewer('div', 'RdBu', 64));
cfg.colormap = cmap;
cfg.gridscale = 300;
cfg.ylim = [1000 1005];
try
 ft_topoplotER(cfg, pow_c25_fooof_bl_smooth);
end
saveas(gcf, fullfile(paths.figures, 'eeg', 'topos', 'ANT128_channel_locations.png'))