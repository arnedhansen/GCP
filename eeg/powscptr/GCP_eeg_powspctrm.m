%% GCP Gamma Peak Power and Frequency

%% Setup
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
    if contains(label, {'O'}) && ~contains(label, {'P'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% GA TESTS
% List of conditions
conditions_list = {'horizontal', 'vertical', 'concentric_static', 'concentric_dynamic_inward', 'concentric_dynamic_outward'};

gapow_data = {gapow_horz_lc_baselined, gapow_horz_hc_baselined; ...
    gapow_vert_lc_baselined, gapow_vert_hc_baselined; ...
    gapow_concentric_static_lc_baselined, gapow_concentric_static_hc_baselined; ...
    gapow_concentric_dynamic_inward_lc_baselined, gapow_concentric_dynamic_inward_hc_baselined; ...
    gapow_concentric_dynamic_outward_lc_baselined, gapow_concentric_dynamic_outward_hc_baselined};

% File path for saving plots
output_path = '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/';

% Loop through all conditions
for cond_idx = 1:length(conditions_list)
    condition = conditions_list{cond_idx};
    gapow_lc = gapow_data{cond_idx, 1};
    gapow_hc = gapow_data{cond_idx, 2};

    % Plotting
    close all;
    figure;
    set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
    cfg = [];
    cfg.channel = channels; % Use the same occipital channels
    cfg.figure = 'gcf';
    cfg.linecolor = 'br';
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
    % ylim([-0.45 0.45])
    % xlim([30 90])
    box on;
    xlabel('Frequency [Hz]');
    ylabel('Power [dB]');
    legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
    title(['Power Spectrum: ', upper(condition)], 'FontSize', 30);
    hold off;

    % Save the plot
    saveas(gcf, strcat(output_path, 'GCP_eeg_gamma_powspctrm_', condition, '.png'));
end

%%
%%
%%
%%
%%

%% Plot GRAND AVERAGE power spectrum
close all;
figure;
set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
colors = {'b', 'r'};
conditions = {'Low Contrast', 'High Contrast'};
cfg = [];
cfg.channel = channels;
cfg.figure = 'gcf';
cfg.linecolor = 'br';
cfg.linewidth = 1;

% Plot for low and high contrast
ft_singleplotER(cfg, gapow_horz_lc, gapow_horz_hc);
hold on;

% Add shaded error bars
channels_seb = ismember(gapow_horz_lc.label, cfg.channel);
lceb = shadedErrorBar(gapow_horz_lc.freq, mean(gapow_horz_lc.powspctrm(channels_seb, :), 1), std(gapow_horz_lc.powspctrm(channels_seb, :))/sqrt(size(gapow_horz_lc.powspctrm(channels_seb, :), 1)), {'b', 'markerfacecolor', 'b'});
hceb = shadedErrorBar(gapow_horz_hc.freq, mean(gapow_horz_hc.powspctrm(channels_seb, :), 1), std(gapow_horz_hc.powspctrm(channels_seb, :))/sqrt(size(gapow_horz_hc.powspctrm(channels_seb, :), 1)), {'r', 'markerfacecolor', 'r'});
transparency = 0.5;
set(lceb.patch, 'FaceAlpha', transparency);
set(hceb.patch, 'FaceAlpha', transparency);

% Adjust plot aesthetics
set(gcf,'color','w');
set(gca,'Fontsize',20);
[~, channel_idx] = ismember(channels, gapow_horz_lc.label);
freq_idx = find(gapow_horz_lc.freq >= 30 & gapow_horz_lc.freq <= 90); % Adjust freq range to gamma
max_spctrm = max([mean(gapow_horz_lc.powspctrm(channel_idx, freq_idx), 2); mean(gapow_horz_hc.powspctrm(channel_idx, freq_idx), 2)]);
ylim([-0.45 0.45])
xlim([30 90])
box on;
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
title('Power Spectrum: HORIZONTAL', 'FontSize', 30);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_eeg_gamma_powspctrm_horz.png');

%% Plot and save INDIVIDUAL power spectra
for subj = 1:length(subjects)
    close all;
    figure;
    set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');

    % Extract participant data
    pow_lc_subj = pow_lc{subj};
    pow_hc_subj = pow_hc{subj};

    cfg = [];
    cfg.channel = channels;
    cfg.figure = 'gcf';
    cfg.linecolor = 'br';
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

%% Plot GRAND AVERAGE power spectrum for each orientation (0, 45, 90, 115)
close all;
figure;
set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
colors = {'b', 'r', 'g', 'm', 'k'};  % Colors for each orientation
cfg = [];
cfg.figure = 'gcf';
cfg.linewidth = 1;
cfg.channel = 'all';  % Specify which channels to use (modify if needed)

% Prepare data for plotting
gapow_hc_dataplotting = {gapow_hc_0, gapow_hc_45, gapow_hc_90, gapow_hc_115, gapow_hc};
gapow_hc_baseline = {gapow_hc_0_baseline_period, gapow_hc_45_baseline_period, gapow_hc_90_baseline_period, gapow_hc_115_baseline_period, gapow_hc_baseline_period};

% Plot for each orientation
hold on;
for i = 1:5
    % Get power spectrum and baseline for this orientation
    power_spectrum = gapow_hc_dataplotting{i}.powspctrm;
    baseline = gapow_hc_baseline{i}.powspctrm;

    % Calculate percentage change relative to the baseline
    percentage_change = ((power_spectrum - baseline) ./ baseline) * 100;

    % Plot the percentage change for this orientation
    shadedErrorBar(gapow_hc_dataplotting{i}.freq, nanmean(percentage_change, 1), nanstd(percentage_change, 0, 1)/sqrt(size(percentage_change, 1)), {'color', colors{i}, 'lineWidth', 2});
end

% Adjust plot aesthetics
set(gca, 'FontSize', 20);
% xlim([30 90]);  % Adjust frequency range to gamma
% ylim([-0.45 0.45]);  % Adjust y-axis limit based on the data
xlabel('Frequency [Hz]', 'FontSize', 25);
ylabel('Percentage Change in Power Spectrum (%)', 'FontSize', 25);
title('Grand Average Power Spectrum for Each Orientation', 'FontSize', 30);
legend({'0°', '45°', '90°', '115°'}, 'FontSize', 20, 'Location', 'northeast');

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_eeg_gamma_powspctrm_orientations_hc.png');
hold off;

