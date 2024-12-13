%% GCP Gamma Topoplots

%% Setup
clear
[subjects, path] = setup('GCP');

%% Load power spectra data
for subj = 1:length(subjects)
    load(strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra'))

    % Low contrast
    power_lc{subj} = pow_lc;
    power_lc_baselined{subj} = pow_lc_baselined;
    power_lc_baseline_period{subj} = pow_lc_baseline_period;

    % High contrast
    power_hc{subj} = pow_hc;
    power_hc_baselined{subj} = pow_hc_baselined;
    power_hc_baseline_period{subj} = pow_hc_baseline_period;
end

%% Compute grand averages
gapow_lc                                               = ft_freqgrandaverage([],power_lc{:});
gapow_lc_baselined                                     = ft_freqgrandaverage([],power_lc_baselined{:});
gapow_lc_baseline_period                               = ft_freqgrandaverage([],power_lc_baseline_period{:});

gapow_hc                                               = ft_freqgrandaverage([],power_hc{:});
gapow_hc_baselined                                     = ft_freqgrandaverage([],power_hc_baselined{:});
gapow_hc_baseline_period                               = ft_freqgrandaverage([],power_hc_baseline_period{:});

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, '/eeg');
cd(datapath);
% Occipital channels
occ_channels = {};
pow_label = power_lc{1, 1};
for i = 1:length(pow_label.label)
    label = pow_label.label{i};
    if contains(label, {'O'}) && ~contains(label, {'P'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot GRAND AVERAGE topoplots for LOW, HIGH and DIFFERENCE
close all
% Define data for coditions
gapow_LOW = gapow_lc_baselined;
gapow_HIGH = gapow_hc_baselined;

% Common configuration
cfg = [];
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat')
cfg.layout = layANThead;
cfg.comment     = 'no';
cfg.gridscale   = 300;
cfg.figure      = 'gcf';
cfg.xlim = [30 90];
cfg.zlim = 'maxabs';
% cfg.zlim = [-1.1 1.1];
cfg.marker      = 'off';
cfg.colormap    = '*RdBu';
cfg.colorbartext = 'Power [dB]';

% Create figure
figure;
set(gcf, 'Position', [0, 0, 2000, 800], 'Color', 'w');
set(gca, 'Fontsize', 25);

% Set title
sgtitle('Topographical Maps 300 ms - 2000 ms after Stimulus Presentation (30 - 90 Hz)', 'FontSize', 30, 'FontWeight', 'bold');

% LOW CONTRAST
subplot(1, 3, 1);
ft_topoplotER(cfg, gapow_LOW);
title('Low Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;
ylabel(cb, 'Power [dB]', 'FontSize', 25);
title('LOW CONTRAST', 'FontSize', 25);

% HIGH CONTRAST
subplot(1, 3, 2);
ft_topoplotER(cfg, gapow_HIGH);
title('High Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;
ylabel(cb, 'Power [dB]', 'FontSize', 25);
title('HIGH CONTRAST', 'FontSize', 25);

% DIFFERENCE
gapow_DIFF = gapow_HIGH;
gapow_DIFF.powspctrm = gapow_HIGH.powspctrm - gapow_LOW.powspctrm;
subplot(1, 3, 3);
ft_topoplotER(cfg, gapow_DIFF);
title('Difference (High - Low)', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;
ylabel(cb, 'Power [dB]', 'FontSize', 25);
title('DIFFERENCE (HC-LC)', 'FontSize', 25);

% Save figure
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/topos/GCP_eeg_topos_ga.png');

%% Topoplots for Individual Subjects (HC, LC, and Difference)
close all;
cfg             = [];
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat')
cfg.layout = layANThead;
cfg.comment     = 'no';
%cfg.channels = channels;
cfg.gridscale   = 300;
cfg.figure      = 'gcf';
cfg.xlim = [30 90];
cfg.zlim = 'maxabs';
cfg.marker      = 'on';
cfg.colormap = '*RdBu';
cfg.colorbartext = 'Power [dB]';

for subj = 1:length(subjects)
    % Prepare data for HC, LC, and Difference
    pow_hc_subj = pow_hc{subj};
    pow_lc_subj = pow_lc{subj};
    pow_diff_subj = pow_hc_subj;
    pow_diff_subj.powspctrm = pow_hc_subj.powspctrm - pow_lc_subj.powspctrm;

    % Set up figure with adjusted layout and style
    figure;
    set(gcf, 'Position', [0, 0, 2000, 800], 'Color', 'w');
    sgtitle(sprintf('Topographical Maps for Subject %s', subjects{subj}), 'FontSize', 30, 'FontWeight', 'bold');

    % High Contrast
    subplot(1, 3, 1);
    ft_topoplotER(cfg, pow_hc_subj);
    title('High Contrast', 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'Power [dB]', 'FontSize', 25);

    % Low Contrast
    subplot(1, 3, 2);
    ft_topoplotER(cfg, pow_lc_subj);
    title('Low Contrast', 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'Power [dB]', 'FontSize', 25);

    % Difference (HC - LC)
    subplot(1, 3, 3);
    ft_topoplotER(cfg, pow_diff_subj);
    title('Difference (High - Low)', 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'Power [dB]', 'FontSize', 25);

    % Save individual figure
    saveas(gcf, sprintf('/Volumes/methlab/Students/Arne/GCP/figures/eeg/topos/GCP_topoplot_subj%s_HC_LC_Diff.png', subjects{subj}));
end

%% Subplot of All Subjects (HC, LC, and Difference)
close all
% Plot all subjects' topoplots in a 3-column grid: HC, LC, and HC-LC difference
figure;
set(gcf, 'Position', [0, 0, 2000, 2000], 'Color', 'w');
sgtitle('Topographical Maps for All Subjects', 'FontSize', 30, 'FontWeight', 'bold');

num_subs = length(subjects);
cols = 3;  % 3 columns for HC, LC, and Difference
rows = num_subs; % Number of rows matches the number of subjects

for subj = 1:num_subs
    pow_hc_subj = pow_hc{subj};
    pow_lc_subj = pow_lc{subj};
    pow_diff_subj = pow_hc_subj;
    pow_diff_subj.powspctrm = pow_hc_subj.powspctrm - pow_lc_subj.powspctrm;

    % High Contrast
    subplot(rows, cols, (subj-1)*cols + 1);
    ft_topoplotER(cfg, pow_hc_subj);
    title('High Contrast', 'FontSize', 10);
    cb = colorbar;
    cb.FontSize = 10;
    ylabel(cb, 'Power [dB]', 'FontSize', 10);

    % Low Contrast
    subplot(rows, cols, (subj-1)*cols + 2);
    ft_topoplotER(cfg, pow_lc_subj);
    title('Low Contrast', 'FontSize', 10);
    cb = colorbar;
    cb.FontSize = 10;
    ylabel(cb, 'Power [dB]', 'FontSize', 10);

    % Difference (HC - LC)
    subplot(rows, cols, (subj-1)*cols + 3);
    ft_topoplotER(cfg, pow_diff_subj);
    title('Difference (High - Low)', 'FontSize', 10);
    cb = colorbar;
    cb.FontSize = 10;
    ylabel(cb, 'Power [dB]', 'FontSize', 10);

end

% Save the combined subplot figure with all subjects
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/topos/GCP_topoplot_allsubs_HC_LC_Diff.png');
