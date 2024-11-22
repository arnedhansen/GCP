%% GCP Gamma Peak Power and Frequency

%% Setup
% setup('GCP');

projectName = 'GCP'
    baseDir = '/Volumes/methlab/Students/Arne/';
    path = fullfile(baseDir, projectName, 'data/features/');
    % List directories in the selected path
    dirs = dir(path);
    folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
    subjects = {folders.name};

    % Display the loaded subjects
    disp('Loaded subjects:');
    disp(subjects);
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
gapow_diff = gapow_hc;  
gapow_diff.powspctrm = gapow_hc.powspctrm - gapow_lc.powspctrm;

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

%% Plot Grand Average Topoplots for HC, LC, and Difference
close all;
cfg             = [];
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat') 
cfg.layout = layANThead; 
cfg.comment     = 'no';
cfg.channels = channels;
cfg.gridscale   = 300;
cfg.figure      = 'gcf';
cfg.xlim = [30 90]; 
% cfg.zlim = 'maxabs'; 
cfg.zlim = [-1.1 1.1];
cfg.marker      = 'on';
cfg.colormap = '*RdBu';
cfg.colorbartext = 'Power [dB]';

% Create figure
figure; 
set(gcf, 'Position', [0, 0, 2000, 800], 'Color', 'w');
set(gca, 'Fontsize', 25);
sgtitle('Topographical Maps', 'FontSize', 30, 'FontWeight', 'bold'); 

% High contrast
subplot(1, 3, 1);
ft_topoplotER(cfg, gapow_hc);
title('High Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;  
ylabel(cb, 'Power [dB]', 'FontSize', 25);

% Low contrast
subplot(1, 3, 2);
ft_topoplotER(cfg, gapow_lc);
title('Low Contrast', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;  
ylabel(cb, 'Power [dB]', 'FontSize', 25);

% Difference (HC - LC)
subplot(1, 3, 3);
ft_topoplotER(cfg, gapow_diff);
title('Difference (High - Low)', 'FontSize', 25);
cb = colorbar;
cb.FontSize = 20;  
ylabel(cb, 'Power [dB]', 'FontSize', 25);  

% Save figure
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/topos/GCP_topoplot_ga_hc_lc_diff.png');

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
