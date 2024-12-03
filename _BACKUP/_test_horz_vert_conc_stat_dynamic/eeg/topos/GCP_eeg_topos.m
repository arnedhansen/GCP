%% GCP Gamma Peak Power and Frequency

%% Setup
[subjects, path] = setup('GCP');

%% Load data and convert TFR data to POWSCPTRM (channels x frequency)
baseline_period = [-0.5 -0.25];
analysis_period = [0 2];
freq_range = [30 120];

for subj = 1%:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath);

    % Load data
    load('data_tfr.mat');

    %% Select analysis and baseline period data
    % (1) Analysis period data, no baseline
    % (2) Analysis period data, baselined
    % (3) Baseline period data (to compare with (1) non-baselined data for percentage change)

    % Horizontal
    pow_horz_lc{subj}                                       = select_data(analysis_period, freq_range, tfr_horz_lc);
    pow_horz_lc_baselined{subj}                             = select_data(analysis_period, freq_range, tfr_horz_lc_bl);
    pow_horz_lc_baseline_period{subj}                       = select_data(baseline_period, freq_range, tfr_horz_lc);

    pow_horz_hc{subj}                                       = select_data(analysis_period, freq_range, tfr_horz_hc);
    pow_horz_hc_baselined{subj}                             = select_data(analysis_period, freq_range, tfr_horz_hc_bl);
    pow_horz_hc_baseline_period{subj}                       = select_data(baseline_period, freq_range, tfr_horz_hc);

    % Vertical
    pow_vert_lc{subj}                                       = select_data(analysis_period, freq_range, tfr_vert_lc);
    pow_vert_lc_baselined{subj}                             = select_data(analysis_period, freq_range, tfr_vert_lc_bl);
    pow_vert_lc_baseline_period{subj}                       = select_data(baseline_period, freq_range, tfr_vert_lc);

    pow_vert_hc{subj}                                       = select_data(analysis_period, freq_range, tfr_vert_hc);
    pow_vert_hc_baselined{subj}                             = select_data(analysis_period, freq_range, tfr_vert_hc_bl);
    pow_vert_hc_baseline_period{subj}                       = select_data(baseline_period, freq_range, tfr_vert_hc);

    % Concentric Static
    pow_concentric_static_lc{subj}                          = select_data(analysis_period, freq_range, tfr_concentric_static_lc);
    pow_concentric_static_lc_baselined{subj}                = select_data(analysis_period, freq_range, tfr_concentric_static_lc_bl);
    pow_concentric_static_lc_baseline_period{subj}          = select_data(baseline_period, freq_range, tfr_concentric_static_lc);

    pow_concentric_static_hc{subj}                          = select_data(analysis_period, freq_range, tfr_concentric_static_hc);
    pow_concentric_static_hc_baselined{subj}                = select_data(analysis_period, freq_range, tfr_concentric_static_hc_bl);
    pow_concentric_static_hc_baseline_period{subj}          = select_data(baseline_period, freq_range, tfr_concentric_static_hc);

    % Concentric Dynamic Inward
    pow_concentric_dynamic_inward_lc{subj}                  = select_data(analysis_period, freq_range, tfr_concentric_dynamic_inward_lc);
    pow_concentric_dynamic_inward_lc_baselined{subj}        = select_data(analysis_period, freq_range, tfr_concentric_dynamic_inward_lc_bl);
    pow_concentric_dynamic_inward_lc_baseline_period{subj}  = select_data(baseline_period, freq_range, tfr_concentric_dynamic_inward_lc);

    pow_concentric_dynamic_inward_hc{subj}                  = select_data(analysis_period, freq_range, tfr_concentric_dynamic_inward_hc);
    pow_concentric_dynamic_inward_hc_baselined{subj}        = select_data(analysis_period, freq_range, tfr_concentric_dynamic_inward_hc_bl);
    pow_concentric_dynamic_inward_hc_baseline_period{subj}  = select_data(baseline_period, freq_range, tfr_concentric_dynamic_inward_hc);

    % Concentric Dynamic Outward
    pow_concentric_dynamic_outward_lc{subj}                 = select_data(analysis_period, freq_range, tfr_concentric_dynamic_outward_lc);
    pow_concentric_dynamic_outward_lc_baselined{subj}       = select_data(analysis_period, freq_range, tfr_concentric_dynamic_outward_lc_bl);
    pow_concentric_dynamic_outward_lc_baseline_period{subj} = select_data(baseline_period, freq_range, tfr_concentric_dynamic_outward_lc);

    pow_concentric_dynamic_outward_hc{subj}                 = select_data(analysis_period, freq_range, tfr_concentric_dynamic_outward_hc);
    pow_concentric_dynamic_outward_hc_baselined{subj}       = select_data(analysis_period, freq_range, tfr_concentric_dynamic_outward_hc_bl);
    pow_concentric_dynamic_outward_hc_baseline_period{subj} = select_data(baseline_period, freq_range, tfr_concentric_dynamic_outward_hc);

    %% Remove time dimension for POWSCPTRM (channels x frequency)
    % Horizontal
    pow_horz_lc{subj}                                       = remove_time_dimension(pow_horz_lc{subj});
    pow_horz_lc_baselined{subj}                             = remove_time_dimension(pow_horz_lc_baselined{subj});
    pow_horz_lc_baseline_period{subj}                       = remove_time_dimension(pow_horz_lc_baseline_period{subj});

    pow_horz_hc{subj}                                       = remove_time_dimension(pow_horz_hc{subj});
    pow_horz_hc_baselined{subj}                             = remove_time_dimension(pow_horz_hc_baselined{subj});
    pow_horz_hc_baseline_period{subj}                       = remove_time_dimension(pow_horz_hc_baseline_period{subj});

    % Vertical
    pow_vert_lc{subj}                                       = remove_time_dimension(pow_vert_lc{subj});
    pow_vert_lc_baselined{subj}                             = remove_time_dimension(pow_vert_lc_baselined{subj});
    pow_vert_lc_baseline_period{subj}                       = remove_time_dimension(pow_vert_lc_baseline_period{subj});

    pow_vert_hc{subj}                                       = remove_time_dimension(pow_vert_hc{subj});
    pow_vert_hc_baselined{subj}                             = remove_time_dimension(pow_vert_hc_baselined{subj});
    pow_vert_hc_baseline_period{subj}                       = remove_time_dimension(pow_vert_hc_baseline_period{subj});

    % Concentric Static
    pow_concentric_static_lc{subj}                          = remove_time_dimension(pow_concentric_static_lc{subj});
    pow_concentric_static_lc_baselined{subj}                = remove_time_dimension(pow_concentric_static_lc_baselined{subj});
    pow_concentric_static_lc_baseline_period{subj}          = remove_time_dimension(pow_concentric_static_lc_baseline_period{subj});

    pow_concentric_static_hc{subj}                          = remove_time_dimension(pow_concentric_static_hc{subj});
    pow_concentric_static_hc_baselined{subj}                = remove_time_dimension(pow_concentric_static_hc_baselined{subj});
    pow_concentric_static_hc_baseline_period{subj}          = remove_time_dimension(pow_concentric_static_hc_baseline_period{subj});

    % Concentric Dynamic Inward
    pow_concentric_dynamic_inward_lc{subj}                  = remove_time_dimension(pow_concentric_dynamic_inward_lc{subj});
    pow_concentric_dynamic_inward_lc_baselined{subj}        = remove_time_dimension(pow_concentric_dynamic_inward_lc_baselined{subj});
    pow_concentric_dynamic_inward_lc_baseline_period{subj}  = remove_time_dimension(pow_concentric_dynamic_inward_lc_baseline_period{subj});

    pow_concentric_dynamic_inward_hc{subj}                  = remove_time_dimension(pow_concentric_dynamic_inward_hc{subj});
    pow_concentric_dynamic_inward_hc_baselined{subj}        = remove_time_dimension(pow_concentric_dynamic_inward_hc_baselined{subj});
    pow_concentric_dynamic_inward_hc_baseline_period{subj}  = remove_time_dimension(pow_concentric_dynamic_inward_hc_baseline_period{subj});

    % Concentric Dynamic Outward
    pow_concentric_dynamic_outward_lc{subj}                 = remove_time_dimension(pow_concentric_dynamic_outward_lc{subj});
    pow_concentric_dynamic_outward_lc_baselined{subj}       = remove_time_dimension(pow_concentric_dynamic_outward_lc_baselined{subj});
    pow_concentric_dynamic_outward_lc_baseline_period{subj} = remove_time_dimension(pow_concentric_dynamic_outward_lc_baseline_period{subj});

    pow_concentric_dynamic_outward_hc{subj}                 = remove_time_dimension(pow_concentric_dynamic_outward_hc{subj});
    pow_concentric_dynamic_outward_hc_baselined{subj}       = remove_time_dimension(pow_concentric_dynamic_outward_hc_baselined{subj});
    pow_concentric_dynamic_outward_hc_baseline_period{subj} = remove_time_dimension(pow_concentric_dynamic_outward_hc_baseline_period{subj});

    fprintf('Subject %.3d/%.3d loaded \n', subj, length(subjects))
end

%% Compute grand averages
% Horizontal
gapow_horz_lc                                        = ft_freqgrandaverage([],pow_horz_lc{subj});
gapow_horz_lc_baselined                              = ft_freqgrandaverage([],pow_horz_lc_baselined{subj});
gapow_horz_lc_baseline_period                        = ft_freqgrandaverage([],pow_horz_lc_baseline_period{subj});

gapow_horz_hc                                        = ft_freqgrandaverage([],pow_horz_hc{subj});
gapow_horz_hc_baselined                              = ft_freqgrandaverage([],pow_horz_hc_baselined{subj});
gapow_horz_hc_baseline_period                        = ft_freqgrandaverage([],pow_horz_hc_baseline_period{subj});

% Vertical
gapow_vert_lc                                        = ft_freqgrandaverage([],pow_vert_lc{subj});
gapow_vert_lc_baselined                              = ft_freqgrandaverage([],pow_vert_lc_baselined{subj});
gapow_vert_lc_baseline_period                        = ft_freqgrandaverage([],pow_vert_lc_baseline_period{subj});

gapow_vert_hc                                        = ft_freqgrandaverage([],pow_vert_hc{subj});
gapow_vert_hc_baselined                              = ft_freqgrandaverage([],pow_vert_hc_baselined{subj});
gapow_vert_hc_baseline_period                        = ft_freqgrandaverage([],pow_vert_hc_baseline_period{subj});

% Concentric Static
gapow_concentric_static_lc                           = ft_freqgrandaverage([],pow_concentric_static_lc{subj});
gapow_concentric_static_lc_baselined                 = ft_freqgrandaverage([],pow_concentric_static_lc_baselined{subj});
gapow_concentric_static_lc_baseline_period           = ft_freqgrandaverage([],pow_concentric_static_lc_baseline_period{subj});

gapow_concentric_static_hc                           = ft_freqgrandaverage([],pow_concentric_static_hc{subj});
gapow_concentric_static_hc_baselined                 = ft_freqgrandaverage([],pow_concentric_static_hc_baselined{subj});
gapow_concentric_static_hc_baseline_period           = ft_freqgrandaverage([],pow_concentric_static_hc_baseline_period{subj});

% Concentric Dynamic Inward
gapow_concentric_dynamic_inward_lc                   = ft_freqgrandaverage([],pow_concentric_dynamic_inward_lc{subj});
gapow_concentric_dynamic_inward_lc_baselined         = ft_freqgrandaverage([],pow_concentric_dynamic_inward_lc_baselined{subj});
gapow_concentric_dynamic_inward_lc_baseline_period   = ft_freqgrandaverage([],pow_concentric_dynamic_inward_lc_baseline_period{subj});

gapow_concentric_dynamic_inward_hc                   = ft_freqgrandaverage([],pow_concentric_dynamic_inward_hc{subj});
gapow_concentric_dynamic_inward_hc_baselined         = ft_freqgrandaverage([],pow_concentric_dynamic_inward_hc_baselined{subj});
gapow_concentric_dynamic_inward_hc_baseline_period   = ft_freqgrandaverage([],pow_concentric_dynamic_inward_hc_baseline_period{subj});

% Concentric Dynamic Outward
gapow_concentric_dynamic_outward_lc                  = ft_freqgrandaverage([],pow_concentric_dynamic_outward_lc{subj});
gapow_concentric_dynamic_outward_lc_baselined        = ft_freqgrandaverage([],pow_concentric_dynamic_outward_lc_baselined{subj});
gapow_concentric_dynamic_outward_lc_baseline_period  = ft_freqgrandaverage([],pow_concentric_dynamic_outward_lc_baseline_period{subj});

gapow_concentric_dynamic_outward_hc                  = ft_freqgrandaverage([],pow_concentric_dynamic_outward_hc{subj});
gapow_concentric_dynamic_outward_hc_baselined        = ft_freqgrandaverage([],pow_concentric_dynamic_outward_hc_baselined{subj});
gapow_concentric_dynamic_outward_hc_baseline_period  = ft_freqgrandaverage([],pow_concentric_dynamic_outward_hc_baseline_period{subj});

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, '/eeg');
cd(datapath);
% Occipital channels
occ_channels = {};
pow_label = pow_horz_lc{1, 1};
for i = 1:length(pow_label.label)
    label = pow_label.label{i};
    if contains(label, {'O'}) || ~contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot Grand Average Topoplots for HC, LC, and Difference
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat')
gapow_data = {gapow_horz_lc_baselined, gapow_horz_hc_baselined; ...
    gapow_vert_lc_baselined, gapow_vert_hc_baselined; ...
    gapow_concentric_static_lc_baselined, gapow_concentric_static_hc_baselined; ...
    gapow_concentric_dynamic_inward_lc_baselined, gapow_concentric_dynamic_inward_hc_baselined; ...
    gapow_concentric_dynamic_outward_lc_baselined, gapow_concentric_dynamic_outward_hc_baselined};

% List of conditions
conditions_list = {'horizontal', 'vertical', 'concentric_static', 'concentric_dynamic_inward', 'concentric_dynamic_outward'};

% File path for saving plots
output_path = '/Volumes/methlab/Students/Arne/GCP/figures/eeg/topos/';

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

% Loop through all conditions
for cond_idx = 1:length(conditions_list)
    close all
    % Create figure
    figure;
    set(gcf, 'Position', [0, 0, 2000, 800], 'Color', 'w');
    set(gca, 'Fontsize', 25);

    % Set conditions and data
    condition = conditions_list{cond_idx};
    gapow_lc = gapow_data{cond_idx, 1};
    gapow_hc = gapow_data{cond_idx, 2};

    % Set title
    sgtitle(strcat('Topographical Maps - ', upper(condition)), 'FontSize', 30, 'FontWeight', 'bold');

    % LOW CONTRAST
    subplot(1, 3, 1);
    ft_topoplotER(cfg, gapow_lc);
    title('Low Contrast', 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'Power [dB]', 'FontSize', 25);
    title('LOW CONTRAST', 'FontSize', 25);

    % HIGH CONTRAST
    subplot(1, 3, 2);
    ft_topoplotER(cfg, gapow_hc);
    title('High Contrast', 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'Power [dB]', 'FontSize', 25);
    title('HIGH CONTRAST', 'FontSize', 25);

    % DIFFERENCE
    gapow_diff = gapow_hc;
    gapow_diff.powspctrm = gapow_hc.powspctrm - gapow_lc.powspctrm;
    subplot(1, 3, 3);
    ft_topoplotER(cfg, gapow_diff);
    title('Difference (High - Low)', 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'Power [dB]', 'FontSize', 25);
    title('DIFFERENCE (HC-LC)', 'FontSize', 25);

    % Save figure
    saveas(gcf, strcat(output_path, 'GCP_eeg_topos_ga_', condition, '.png'));
end

%% Topoplots for Individual Subjects (HC, LC, and Difference)
close all;
cfg             = [];
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat') 
cfg.layout = layANThead; 
cfg.comment     = 'no';
cfg.channels = channels;
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
    
    % Plot HC
    subplot(rows, cols, (subj-1)*cols + 1);
    ft_topoplotER(cfg, pow_hc_subj);
    title(sprintf('Subj %s: High Contrast', subjects{subj}), 'FontSize', 12);

    % Plot LC
    subplot(rows, cols, (subj-1)*cols + 2);
    ft_topoplotER(cfg, pow_lc_subj);
    title(sprintf('Subj %s: Low Contrast', subjects{subj}), 'FontSize', 12);
    
    % Plot Difference (HC - LC)
    subplot(rows, cols, (subj-1)*cols + 3);
    ft_topoplotER(cfg, pow_diff_subj);
    title(sprintf('Subj %s: High - Low', subjects{subj}), 'FontSize', 12);
end

% Save the combined subplot figure with all subjects
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/topos/GCP_topoplot_allsubs_HC_LC_Diff.png');
