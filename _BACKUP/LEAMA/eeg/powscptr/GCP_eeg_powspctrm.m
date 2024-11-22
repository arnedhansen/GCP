%% GCP Gamma Peak Power and Frequency

%% Setup
[subjects, path] = setup('GCP');

%% Load data and convert TFR data to POWSCPTRM (channels x frequency)
orientations = {'0', '45', '90', '115', 'all'};
baseline_period = [-0.5 -0.25];
analysis_period = [0 2];
freq_range = [30 120];

% Load and process data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath);

    % Load high and low contrast data
    load('tfr_hc.mat');
    load('tfr_lc.mat');

    % Initialise containers for the subject-level datasets
    for ori = 1:length(orientations)
        orientation = orientations{ori};

        % HC data and baseline data
        eval(sprintf('data_hc = tfr_hc_%s;', orientation));
        eval(sprintf('data_hc_bl = tfr_hc_%s_bl;', orientation));

        % LC data and baseline data
        eval(sprintf('data_lc = tfr_lc_%s;', orientation));
        eval(sprintf('data_lc_bl = tfr_lc_%s_bl;', orientation));

        % Select analysis and baseline periods for HC and LC
        pow_hc{subj, ori} = select_data(analysis_period, freq_range, data_hc);
        bl_hc{subj, ori}  = select_data(baseline_period, freq_range, data_hc);
        bc_hc{subj, ori}  = select_data(analysis_period, freq_range, data_hc_bl); 

        pow_lc{subj, ori} = select_data(analysis_period, freq_range, data_lc);
        bl_lc{subj, ori}  = select_data(baseline_period, freq_range, data_lc);
        bc_lc{subj, ori}  = select_data(analysis_period, freq_range, data_lc_bl); 

        % Remove time dimension for HC and LC data
        pow_hc{subj, ori} = remove_time_dimension(pow_hc{subj, ori});
        bl_hc{subj, ori}  = remove_time_dimension(bl_hc{subj, ori});
        bc_hc{subj, ori}  = remove_time_dimension(bc_hc{subj, ori});

        pow_lc{subj, ori} = remove_time_dimension(pow_lc{subj, ori});
        bl_lc{subj, ori}  = remove_time_dimension(bl_lc{subj, ori});
        bc_lc{subj, ori}  = remove_time_dimension(bc_lc{subj, ori});
    end

    fprintf('Subject %.2s / %.3d loaded \n', num2str(subj), length(subjects))
end

% Compute grand averages for HC and LC for each orientation and condition
gapow_hc_0   = ft_freqgrandaverage([], pow_hc{:, 1});  % 0 degree HC data
gapow_hc_45  = ft_freqgrandaverage([], pow_hc{:, 2});  % 45 degree HC data
gapow_hc_90  = ft_freqgrandaverage([], pow_hc{:, 3});  % 90 degree HC data
gapow_hc_115 = ft_freqgrandaverage([], pow_hc{:, 4});  % 115 degree HC data
gapow_hc     = ft_freqgrandaverage([], pow_hc{:, 5});  % 'all' HC data

gapow_hc_0_baseline_period   = ft_freqgrandaverage([], bl_hc{:, 1});  % 0 degree HC data
gapow_hc_45_baseline_period  = ft_freqgrandaverage([], bl_hc{:, 2});  % 45 degree HC data
gapow_hc_90_baseline_period  = ft_freqgrandaverage([], bl_hc{:, 3});  % 90 degree HC data
gapow_hc_115_baseline_period = ft_freqgrandaverage([], bl_hc{:, 4});  % 115 degree HC data
gapow_hc_baseline_period     = ft_freqgrandaverage([], bl_hc{:, 5});  % 'all' HC data

gapow_lc_0   = ft_freqgrandaverage([], pow_lc{:, 1});  % 0 degree LC data
gapow_lc_45  = ft_freqgrandaverage([], pow_lc{:, 2});  % 45 degree LC data
gapow_lc_90  = ft_freqgrandaverage([], pow_lc{:, 3});  % 90 degree LC data
gapow_lc_115 = ft_freqgrandaverage([], pow_lc{:, 4});  % 115 degree LC data
gapow_lc     = ft_freqgrandaverage([], pow_lc{:, 5});  % 'all' LC data

gapow_hc_0_baseline_period   = ft_freqgrandaverage([], bl_hc{:, 1});  % 0 degree HC data
gapow_hc_45_baseline_period  = ft_freqgrandaverage([], bl_hc{:, 2});  % 45 degree HC data
gapow_hc_90_baseline_period  = ft_freqgrandaverage([], bl_hc{:, 3});  % 90 degree HC data
gapow_hc_115_baseline_period = ft_freqgrandaverage([], bl_hc{:, 4});  % 115 degree HC data
gapow_hc_baseline_period     = ft_freqgrandaverage([], bl_hc{:, 5});  % 'all' HC data

%% TEST ONLY HC

orientations = {'0', '45', '90', '115', 'all'};
baseline_period = [-0.5 -0.25];
analysis_period = [0 2];
freq_range = [30 120];

% Load and process data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath);

    % Load high and low contrast data
    load('tfr_hc.mat');

    % Initialise containers for the subject-level datasets
    for ori = 1:length(orientations)
        orientation = orientations{ori};

        % HC data and baseline data
        eval(sprintf('data_hc = tfr_hc_%s;', orientation));
        eval(sprintf('data_hc_bl = tfr_hc_%s_bl;', orientation));

       
        % Select analysis and baseline periods for HC and LC
        pow_hc{subj, ori} = select_data(analysis_period, freq_range, data_hc);
        bl_hc{subj, ori}  = select_data(baseline_period, freq_range, data_hc);
        bc_hc{subj, ori}  = select_data(analysis_period, freq_range, data_hc_bl); 


        % Remove time dimension for HC and LC data
        pow_hc{subj, ori} = remove_time_dimension(pow_hc{subj, ori});
        bl_hc{subj, ori}  = remove_time_dimension(bl_hc{subj, ori});
        bc_hc{subj, ori}  = remove_time_dimension(bc_hc{subj, ori});

    end

    fprintf('Subject %.2s / %.3d loaded \n', num2str(subj), length(subjects))
end

% Compute grand averages for HC and LC for each orientation and condition
gapow_hc_0   = ft_freqgrandaverage([], pow_hc{:, 1});  % 0 degree HC data
gapow_hc_45  = ft_freqgrandaverage([], pow_hc{:, 2});  % 45 degree HC data
gapow_hc_90  = ft_freqgrandaverage([], pow_hc{:, 3});  % 90 degree HC data
gapow_hc_115 = ft_freqgrandaverage([], pow_hc{:, 4});  % 115 degree HC data
gapow_hc     = ft_freqgrandaverage([], pow_hc{:, 5});  % 'all' HC data

gapow_hc_0_baseline_period   = ft_freqgrandaverage([], bl_hc{:, 1});  % 0 degree HC data
gapow_hc_45_baseline_period  = ft_freqgrandaverage([], bl_hc{:, 2});  % 45 degree HC data
gapow_hc_90_baseline_period  = ft_freqgrandaverage([], bl_hc{:, 3});  % 90 degree HC data
gapow_hc_115_baseline_period = ft_freqgrandaverage([], bl_hc{:, 4});  % 115 degree HC data
gapow_hc_baseline_period     = ft_freqgrandaverage([], bl_hc{:, 5});  % 'all' HC data

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, '/eeg');
cd(datapath);
% Occipital channels
occ_channels = {};
pow_label = pow_hc{1, 1};
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
colors = {'b', 'r'};
conditions = {'Low Contrast', 'High Contrast'};
cfg = [];
% cfg.channel = channels;
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
ylim([-0.45 0.45])
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

