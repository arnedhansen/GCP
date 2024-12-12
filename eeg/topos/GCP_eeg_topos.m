%% GCP Gamma Topoplots

%% Setup
startup
[subjects, path] = setup('GCP');

%% Load data and convert TFR data to POWSCPTRM (channels x frequency)
baseline_period = [-0.5 -0.25];
analysis_period = [0 2];
freq_range = [30 120];

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath);

    % Load data
    load('data_tfr.mat');

    %% Select analysis and baseline period data
    % (1) Analysis period data, no baseline
    % (2) Analysis period data, baselined
    % (3) Baseline period data (to compare with (1) non-baselined data for percentage change)

    pow_lc{subj}                                       = select_data(analysis_period, freq_range, tfr_lc);
    pow_lc_baselined{subj}                             = select_data(analysis_period, freq_range, tfr_lc_bl);
    pow_lc_baseline_period{subj}                       = select_data(baseline_period, freq_range, tfr_lc);

    pow_hc{subj}                                       = select_data(analysis_period, freq_range, tfr_hc);
    pow_hc_baselined{subj}                             = select_data(analysis_period, freq_range, tfr_hc_bl);
    pow_hc_baseline_period{subj}                       = select_data(baseline_period, freq_range, tfr_hc);

    %% Remove time dimension for POWSCPTRM (channels x frequency)
    pow_lc{subj}                                       = remove_time_dimension(pow_lc{subj});
    pow_lc_baselined{subj}                             = remove_time_dimension(pow_lc_baselined{subj});
    pow_lc_baseline_period{subj}                       = remove_time_dimension(pow_lc_baseline_period{subj});

    pow_hc{subj}                                       = remove_time_dimension(pow_hc{subj});
    pow_hc_baselined{subj}                             = remove_time_dimension(pow_hc_baselined{subj});
    pow_hc_baseline_period{subj}                       = remove_time_dimension(pow_hc_baseline_period{subj});

    fprintf('Subject %.3d/%.3d loaded \n', subj, length(subjects))
end

%% Compute grand averages
gapow_lc                                        = ft_freqgrandaverage([],pow_lc{subj});
gapow_lc_baselined                              = ft_freqgrandaverage([],pow_lc_baselined{subj});
gapow_lc_baseline_period                        = ft_freqgrandaverage([],pow_lc_baseline_period{subj});

gapow_hc                                        = ft_freqgrandaverage([],pow_hc{subj});
gapow_hc_baselined                              = ft_freqgrandaverage([],pow_hc_baselined{subj});
gapow_hc_baseline_period                        = ft_freqgrandaverage([],pow_hc_baseline_period{subj});

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, '/eeg');
cd(datapath);
% Occipital channels
occ_channels = {};
pow_label = pow_lc{1, 1};
for i = 1:length(pow_label.label)
    label = pow_label.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
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
