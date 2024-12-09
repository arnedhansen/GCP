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
gapow_lc                                               = ft_freqgrandaverage([],pow_lc{subj});
gapow_lc_baselined                                     = ft_freqgrandaverage([],pow_lc_baselined{subj});
gapow_lc_baseline_period                               = ft_freqgrandaverage([],pow_lc_baseline_period{subj});
       
gapow_hc                                               = ft_freqgrandaverage([],pow_hc{subj});
gapow_hc_baselined                                     = ft_freqgrandaverage([],pow_hc_baselined{subj});
gapow_hc_baseline_period                               = ft_freqgrandaverage([],pow_hc_baseline_period{subj});

%% Define channels
datapath = strcat(path, subjects{1}, '/eeg');
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

%% Extract gamma peak power and frequency
eeg_data = [];

for subj = 1:length(subjects)
    % Extract participant data
    pow_lc_subj = pow_lc_baselined{subj};
    pow_hc_subj = pow_hc_baselined{subj};

    % Find channels and frequencies of interest
    channels_idx = ismember(pow_lc_subj.label, channels);
    freq_idx = find(pow_lc_subj.freq >= 30 & pow_hc_subj.freq <= 90);

    % Find gamma peak for low contrast
    lc_gamma_power = mean(pow_lc_subj.powspctrm(channels_idx, freq_idx), 1);
    [lc_pow, lc_peak_idx] = max(lc_gamma_power);
    lc_freq = pow_lc_subj.freq(freq_idx(lc_peak_idx));

    % Find gamma peak for high contrast
    hc_gamma_power = mean(pow_hc_subj.powspctrm(channels_idx, freq_idx), 1);
    [hc_pow, hc_peak_idx] = max(hc_gamma_power);
    hc_freq = pow_hc_subj.freq(freq_idx(hc_peak_idx));

    % Create across condition structure
    subject_id = [str2num(subjects{subj}); str2num(subjects{subj})];
    subj_data_eeg = struct('ID', num2cell(subject_id(1:2)), 'Condition', num2cell([1; 2]), 'Power', num2cell([lc_pow; hc_pow]), 'Frequency', num2cell([lc_freq; hc_freq]));

    % Save data
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save eeg_matrix subj_data_eeg
    save pow lc_pow hc_pow
    save freq lc_freq hc_freq
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

    % Append to the final structure array
    eeg_data = [eeg_data; subj_data_eeg];
end
save /Volumes/methlab/Students/Arne/GCP/data/features/eeg_matrix eeg_data

%% Plot GRAND AVERAGE power spectrum BASELINED
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

% Define data for coditions
gapow_LOW = gapow_lc_baselined;
gapow_HIGH = gapow_hc_baselined;

% Plot for low and high contrast
ft_singleplotER(cfg, gapow_LOW, gapow_HIGH);
hold on;

% Add shaded error bars
channels_seb = ismember(gapow_LOW.label, cfg.channel);
lceb = shadedErrorBar(gapow_LOW.freq, mean(gapow_LOW.powspctrm(channels_seb, :), 1), std(gapow_LOW.powspctrm(channels_seb, :))/sqrt(size(gapow_LOW.powspctrm(channels_seb, :), 1)), {'b', 'markerfacecolor', 'b'});
hceb = shadedErrorBar(gapow_HIGH.freq, mean(gapow_HIGH.powspctrm(channels_seb, :), 1), std(gapow_HIGH.powspctrm(channels_seb, :))/sqrt(size(gapow_HIGH.powspctrm(channels_seb, :), 1)), {'r', 'markerfacecolor', 'r'});
transparency = 0.5;
set(lceb.patch, 'FaceAlpha', transparency);
set(hceb.patch, 'FaceAlpha', transparency);

% Adjust plot aesthetics
set(gcf,'color','w');
set(gca,'Fontsize',20);
[~, channel_idx] = ismember(channels, gapow_LOW.label);
freq_idx = find(gapow_LOW.freq >= 30 & gapow_LOW.freq <= 90); % Adjust freq range to gamma
max_spctrm = max([mean(gapow_LOW.powspctrm(channel_idx, freq_idx), 2); mean(gapow_HIGH.powspctrm(channel_idx, freq_idx), 2)]);
ylim([-1 1])
xlim([30 90])

xlabel('Frequency [Hz]');
ylabel('Power [dB]');
legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
title('Power Spectrum', 'FontSize', 30);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_eeg_gamma_powspctrm_baselined.png');

%% Plot GRAND AVERAGE power spectrum PERCENTAGE CHANGE
gapow_LOW = gapow_lc; % Non-baselined low-contrast data
gapow_HIGH = gapow_hc; % Non-baselined high-contrast data
gapow_BASELINE_PERIOD_LOW = gapow_lc_baseline_period; % Baseline-period low-contrast data
gapow_BASELINE_PERIOD_HIGH = gapow_hc_baseline_period; % Baseline-period high-contrast data

% Preallocate arrays for percentage change
percent_change_LOW = gapow_LOW;
percent_change_HIGH = gapow_HIGH;

% Percentage change formula: ((stimulus - baseline) / baseline) * 100
percent_change_LOW.powspctrm = ((gapow_LOW.powspctrm - gapow_BASELINE_PERIOD_LOW.powspctrm) ./ gapow_BASELINE_PERIOD_LOW.powspctrm) * 100;
percent_change_HIGH.powspctrm = ((gapow_HIGH.powspctrm - gapow_BASELINE_PERIOD_HIGH.powspctrm) ./ gapow_BASELINE_PERIOD_HIGH.powspctrm) * 100;

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
ft_singleplotER(cfg, percent_change_LOW, percent_change_HIGH);
hold on;

% Add shaded error bars
channels_seb = ismember(percent_change_LOW.label, cfg.channel);
lceb = shadedErrorBar(percent_change_LOW.freq, mean(percent_change_LOW.powspctrm(channels_seb, :), 1), std(percent_change_LOW.powspctrm(channels_seb, :))/sqrt(size(percent_change_LOW.powspctrm(channels_seb, :), 1)), {'b', 'markerfacecolor', 'b'});
hceb = shadedErrorBar(percent_change_HIGH.freq, mean(percent_change_HIGH.powspctrm(channels_seb, :), 1), std(percent_change_HIGH.powspctrm(channels_seb, :))/sqrt(size(percent_change_HIGH.powspctrm(channels_seb, :), 1)), {'r', 'markerfacecolor', 'r'});
transparency = 0.5;
set(lceb.patch, 'FaceAlpha', transparency);
set(hceb.patch, 'FaceAlpha', transparency);

% Extract stats for gamma peak power and frequency
% Find gamma peak for low contrast
lc_gamma_power = mean(percent_change_LOW.powspctrm(channels_seb, freq_idx), 1);
[lc_pow, lc_peak_idx] = max(lc_gamma_power);
lc_freq = percent_change_LOW.freq(freq_idx(lc_peak_idx));

% Find gamma peak for high contrast
hc_gamma_power = mean(percent_change_HIGH.powspctrm(channels_seb, freq_idx), 1);
[hc_pow, hc_peak_idx] = max(hc_gamma_power);
hc_freq = percent_change_HIGH.freq(freq_idx(hc_peak_idx));

% Adjust plot aesthetics
yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5);
plot([0 lc_freq], [lc_pow lc_pow], '--', 'Color', 'b', 'LineWidth', 0.5);
plot([lc_freq lc_freq], [-100 lc_pow], '--', 'Color', 'b', 'LineWidth', 0.5);
plot([0 hc_freq], [hc_pow hc_pow], '--', 'Color', 'r', 'LineWidth', 0.5);
plot([hc_freq hc_freq], [-100 hc_pow], '--', 'Color', 'r', 'LineWidth', 0.5);
set(gcf,'color','w');
set(gca,'Fontsize',20);
max_spctrm = max(lc_pow, hc_pow);
ylim([-max_spctrm*1.2 max_spctrm*1.2]);
xlim([30 90]) % Gamma frequency range
xlabel('Frequency [Hz]');
ylabel('Percentage Change [%]');
legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
title('Percentage Change Power Spectrum', 'FontSize', 30);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_eeg_gamma_powspctrm_percentage.png');

%% Plot and save INDIVIDUAL power spectra
output_dir = '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/';

for subj = 1:length(subjects)
    close all;
    figure;
    set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');

    % Extract participant data
    pow_lc_subj = pow_lc_baselined{subj};
    pow_hc_subj = pow_hc_baselined{subj};

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
    lceb = shadedErrorBar(pow_lc_subj.freq, mean(pow_lc_subj.powspctrm(channels_seb, :), 1), ...
        std(pow_lc_subj.powspctrm(channels_seb, :)) / sqrt(size(pow_lc_subj.powspctrm(channels_seb, :), 1)), {'b', 'markerfacecolor', 'b'});
    hceb = shadedErrorBar(pow_hc_subj.freq, mean(pow_hc_subj.powspctrm(channels_seb, :), 1), ...
        std(pow_hc_subj.powspctrm(channels_seb, :)) / sqrt(size(pow_hc_subj.powspctrm(channels_seb, :), 1)), {'r', 'markerfacecolor', 'r'});
    transparency = 0.5;
    set(lceb.patch, 'FaceAlpha', transparency);
    set(hceb.patch, 'FaceAlpha', transparency);

    % Save results
    gamma_power_frequency(subj).subject = subjects{subj};
    gamma_power_frequency(subj).low_contrast_amp = lc_pow;
    gamma_power_frequency(subj).low_contrast_freq = lc_freq;
    gamma_power_frequency(subj).high_contrast_amp = hc_pow;
    gamma_power_frequency(subj).high_contrast_freq = hc_freq;

    % Adjust plot aesthetics
    set(gca, 'FontSize', 20);
    plot([0 lc_freq], [lc_pow lc_pow], '--', 'Color', 'b', 'LineWidth', 0.5);
    plot([lc_freq lc_freq], [-100 lc_pow], '--', 'Color', 'b', 'LineWidth', 0.5);
    plot([0 hc_freq], [hc_pow hc_pow], '--', 'Color', 'r', 'LineWidth', 0.5);
    plot([hc_freq hc_freq], [-100 hc_pow], '--', 'Color', 'r', 'LineWidth', 0.5);
    max_spctrm = max(lc_pow, hc_pow);
    ylim([-max_spctrm*1.2 max_spctrm*1.2]);
    xlim([30 90]);
    xlabel('Frequency [Hz]');
    ylabel('Power [dB]');
    legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
    title(sprintf('Subject %s: Power Spectrum', subjects{subj}), 'FontSize', 30);
    hold off;

    % Save individual plot
    save_path = fullfile(output_dir, sprintf('GCP_eeg_gamma_powspctrm_subj%s.png', subjects{subj}));
    saveas(gcf, save_path);

end

%% Subplot with all INDIVIDUAL plots
output_dir = '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/';
num_subs = length(subjects);
cols = 5;  % Number of columns
rows = ceil(num_subs / cols); % Number of rows dynamically calculated

figure;
set(gcf, 'Position', [0, 0, 1600, 1200], 'Color', 'w'); % Adjust overall figure size

for subj = 1:num_subs
    % Create subplot for each subject
    subplot(rows, cols, subj);

    % Extract participant data
    pow_lc_subj = pow_lc_baselined{subj};
    pow_hc_subj = pow_hc_baselined{subj};

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
    set(gca, 'FontSize', 12);
    [~, channel_idx] = ismember(channels, pow_lc_subj.label);
    freq_idx = find(pow_lc_subj.freq >= 30 & pow_lc_subj.freq <= 80); % Adjust freq range to gamma
    max_spctrm = max([mean(pow_lc_subj.powspctrm(channel_idx, freq_idx), 2); mean(pow_hc_subj.powspctrm(channel_idx, freq_idx), 2)]);
    ylim([-max_spctrm*1.15 max_spctrm*1.15]);
    xlim([30 90]);
    xlabel('Freq [Hz]', 'FontSize', 10);
    ylabel('Power [dB]', 'FontSize', 10);
    title(sprintf('Subj %s', subjects{subj}), 'FontSize', 12);
    hold off;
end

% Save the combined figure with all subplots
save_path = fullfile(output_dir, 'GCP_eeg_gamma_powspctrm_all_subjects_subplot.png');
saveas(gcf, save_path);
