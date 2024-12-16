%% GCP Gamma Peak Power and Frequency

%% Setup
clear
[subjects, path, colors] = setup('GCP');
anal_period = 1


%% Load power spectra data
for subj = 1:length(subjects)
    if anal_period == 1
        load(strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra_300'))
    else
        load(strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra'))
    end

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
datapath = strcat(path, subjects{1}, '/eeg');
cd(datapath);
% Occipital channels
occ_channels = {};
pow_label = pow_lc;
for i = 1:length(pow_label.label)
    label = pow_label.label{i};
    if contains(label, {'O'}) && ~contains(label, {'P'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot GRAND AVERAGE power spectrum BASELINED
close all;
figure;
set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
cfg = [];
cfg.channel = channels;
cfg.figure = 'gcf';
cfg.linewidth = 1;

% Plot for low and high contrast
ft_singleplotER(cfg, gapow_lc_baselined, gapow_hc_baselined);
hold on;

% Add shaded error bars
channels_seb = ismember(gapow_lc_baselined.label, cfg.channel);
lceb = shadedErrorBar(gapow_lc_baselined.freq, mean(gapow_lc_baselined.powspctrm(channels_seb, :), 1), std(gapow_lc_baselined.powspctrm(channels_seb, :))/sqrt(size(gapow_lc_baselined.powspctrm(channels_seb, :), 1)), {'-'}, 0);
hceb = shadedErrorBar(gapow_hc_baselined.freq, mean(gapow_hc_baselined.powspctrm(channels_seb, :), 1), std(gapow_hc_baselined.powspctrm(channels_seb, :))/sqrt(size(gapow_hc_baselined.powspctrm(channels_seb, :), 1)), {'-'}, 0);
lceb.mainLine.Color = colors(1, :);
hceb.mainLine.Color = colors(2, :);
lceb.patch.FaceColor = colors(1, :);
hceb.patch.FaceColor = colors(2, :);
set(lceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
set(hceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
set(lceb.edge(1), 'Color', colors(1, :));
set(lceb.edge(2), 'Color', colors(1, :));
set(hceb.edge(1), 'Color', colors(2, :));
set(hceb.edge(2), 'Color', colors(2, :));
set(lceb.patch, 'FaceAlpha', 0.5);
set(hceb.patch, 'FaceAlpha', 0.5);

% Find GA gamma peak for LOW contrast
freq_idx = find(gapow_lc_baselined.freq >= 30 & gapow_lc_baselined.freq <= 90); % Adjust freq range to gamma
lc_gamma_power = mean(gapow_lc_baselined.powspctrm(channels_seb, freq_idx), 1);
[peaks, locs] = findpeaks(lc_gamma_power, gapow_lc_baselined.freq(freq_idx));
[lc_pow, peak_idx] = max(peaks);
lc_freq = locs(peak_idx);

% Find GA gamma peak for HIGH contrast
hc_gamma_power = mean(gapow_hc_baselined.powspctrm(channels_seb, freq_idx), 1);
[peaks, locs] = findpeaks(hc_gamma_power, gapow_hc_baselined.freq(freq_idx));
[hc_pow, peak_idx] = max(peaks);
hc_freq = locs(peak_idx);

% Adjust plot aesthetics
yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
plot([0 lc_freq], [lc_pow lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
plot([lc_freq lc_freq], [-100 lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
plot([0 hc_freq], [hc_pow hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
plot([hc_freq hc_freq], [-100 hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
set(gcf,'color','w');
set(gca,'Fontsize',20);
max_spctrm = max(lc_pow, hc_pow);
ylim([-max_spctrm*1.25 max_spctrm*1.25]);
xlim([30 90])
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
title('Grand Average Power Spectrum', 'FontSize', 30);
hold off;

% Save the plot
if anal_period == 1
    saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_powspctrm_baselined_300.png');
else
    saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_powspctrm_baselined.png');
end

%% Plot GRAND AVERAGE power spectrum PERCENTAGE CHANGE
% Percentage change formula: ((stimulus - baseline) / baseline) * 100
percent_change_LOW = gapow_lc;
percent_change_HIGH = gapow_hc;
percent_change_LOW.powspctrm = ((gapow_lc.powspctrm - gapow_lc_baseline_period.powspctrm) ./ gapow_lc_baseline_period.powspctrm) * 100;
percent_change_HIGH.powspctrm = ((gapow_hc.powspctrm - gapow_hc_baseline_period.powspctrm) ./ gapow_hc_baseline_period.powspctrm) * 100;

close all;
figure;
set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
cfg = [];
cfg.channel = channels;
cfg.figure = 'gcf';
cfg.linewidth = 1;

% Plot for low and high contrast
ft_singleplotER(cfg, percent_change_LOW, percent_change_HIGH);
hold on;

% Add shaded error bars
channels_seb = ismember(percent_change_LOW.label, cfg.channel);
lceb = shadedErrorBar(percent_change_LOW.freq, mean(percent_change_LOW.powspctrm(channels_seb, :), 1), std(percent_change_LOW.powspctrm(channels_seb, :))/sqrt(size(percent_change_LOW.powspctrm(channels_seb, :), 1)), {'-'}, 0);
hceb = shadedErrorBar(percent_change_HIGH.freq, mean(percent_change_HIGH.powspctrm(channels_seb, :), 1), std(percent_change_HIGH.powspctrm(channels_seb, :))/sqrt(size(percent_change_HIGH.powspctrm(channels_seb, :), 1)), {'-'}, 0);
lceb.mainLine.Color = colors(1, :);
hceb.mainLine.Color = colors(2, :);
lceb.patch.FaceColor = colors(1, :);
hceb.patch.FaceColor = colors(2, :);
set(lceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
set(hceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
set(lceb.edge(1), 'Color', colors(1, :));
set(lceb.edge(2), 'Color', colors(1, :));
set(hceb.edge(1), 'Color', colors(2, :));
set(hceb.edge(2), 'Color', colors(2, :));
set(lceb.patch, 'FaceAlpha', 0.5);
set(hceb.patch, 'FaceAlpha', 0.5);

% Find GA gamma peak for LOW contrast
lc_gamma_power = mean(percent_change_LOW.powspctrm(channels_seb, freq_idx), 1);
[peaks, locs] = findpeaks(lc_gamma_power, percent_change_LOW.freq(freq_idx));
[lc_pow, peak_idx] = max(peaks);
lc_freq = locs(peak_idx);

% Find GA gamma peak for HIGH contrast
hc_gamma_power = mean(percent_change_HIGH.powspctrm(channels_seb, freq_idx), 1);
[peaks, locs] = findpeaks(hc_gamma_power, percent_change_HIGH.freq(freq_idx));
[hc_pow, peak_idx] = max(peaks);
hc_freq = locs(peak_idx);

% Adjust plot aesthetics
yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
plot([0 lc_freq], [lc_pow lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
plot([lc_freq lc_freq], [-100 lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
plot([0 hc_freq], [hc_pow hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
plot([hc_freq hc_freq], [-100 hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
set(gcf,'color','w');
set(gca,'Fontsize',20);
max_spctrm = max(lc_pow, hc_pow);
ylim([-max_spctrm*1.25 max_spctrm*1.25]);
xlim([30 90]) % Gamma frequency range
xlabel('Frequency [Hz]');
ylabel('Percentage Change [%]');
legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
title('Percentage Change Power Spectrum', 'FontSize', 30);
hold off;

% Save the plot
if anal_period == 1
    saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_powspctrm_percentage_300.png');
else
    saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_powspctrm_percentage.png');
end

%% Plot INDIVIDUAL power spectra BASELINED
output_dir = '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/';

for subj = 1:length(subjects)
    close all;
    figure;
    set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');

    % Extract participant data
    pow_lc_subj = power_lc_baselined{subj};
    pow_hc_subj = power_hc_baselined{subj};

    % Figure common config
    cfg = [];
    cfg.channel = channels;
    cfg.figure = 'gcf';
    cfg.linewidth = 1;

    % Plot power spectrum for low and high contrast
    ft_singleplotER(cfg, pow_lc_subj, pow_hc_subj);
    hold on;

    % Add shaded error bars
    channels_seb = ismember(pow_lc_subj.label, cfg.channel);
    lceb = shadedErrorBar(pow_lc_subj.freq, mean(pow_lc_subj.powspctrm(channels_seb, :), 1), ...
        std(pow_lc_subj.powspctrm(channels_seb, :)) / sqrt(size(pow_lc_subj.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    hceb = shadedErrorBar(pow_hc_subj.freq, mean(pow_hc_subj.powspctrm(channels_seb, :), 1), ...
        std(pow_hc_subj.powspctrm(channels_seb, :)) / sqrt(size(pow_hc_subj.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    lceb.mainLine.Color = colors(1, :);
    hceb.mainLine.Color = colors(2, :);
    lceb.patch.FaceColor = colors(1, :);
    hceb.patch.FaceColor = colors(2, :);
    set(lceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
    set(hceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
    set(lceb.edge(1), 'Color', colors(1, :));
    set(lceb.edge(2), 'Color', colors(1, :));
    set(hceb.edge(1), 'Color', colors(2, :));
    set(hceb.edge(2), 'Color', colors(2, :));
    set(lceb.patch, 'FaceAlpha', 0.5);
    set(hceb.patch, 'FaceAlpha', 0.5);

    % Find channels and frequencies of interest
    channels_idx = ismember(pow_lc_subj.label, channels);
    freq_idx = find(pow_lc_subj.freq >= 30 & pow_hc_subj.freq <= 90);

    % Find gamma peak for LOW contrast
    lc_gamma_power = mean(pow_lc_subj.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(lc_gamma_power, pow_lc_subj.freq(freq_idx));
    [lc_pow, peak_idx] = max(peaks);
    lc_freq = locs(peak_idx);

    % Find gamma peak for HIGH contrast
    hc_gamma_power = mean(pow_hc_subj.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(hc_gamma_power, pow_hc_subj.freq(freq_idx));
    [hc_pow, peak_idx] = max(peaks);
    hc_freq = locs(peak_idx);
    if subj == 9
        hc_pow = 0.1149
        hc_freq = 86
    end

    % Adjust plot aesthetics
    yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
    plot([0 lc_freq], [lc_pow lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
    plot([lc_freq lc_freq], [-100 lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
    plot([0 hc_freq], [hc_pow hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
    plot([hc_freq hc_freq], [-100 hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
    set(gca, 'FontSize', 20);
    max_spctrm = max(lc_pow, hc_pow);
    if subj == 9
        max_spctrm = 0.4
    elseif subj == 10
        max_spctrm = 0.55
    end
    ylim([-max_spctrm*1.25 max_spctrm*1.25]);
    xlim([30 90]);
    xlabel('Frequency [Hz]');
    ylabel('Power [dB]');
    legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
    title(sprintf('Subject %s: Power Spectrum', subjects{subj}), 'FontSize', 30);
    hold off;

    % Save individual plot
    if anal_period == 1
        save_path = fullfile(output_dir, sprintf('GCP_powspctrm_subj%s_baselined_300.png', subjects{subj}));
    else
        save_path = fullfile(output_dir, sprintf('GCP_powspctrm_subj%s_baselined.png', subjects{subj}));
    end
    saveas(gcf, save_path);
end

%% Subplot with all INDIVIDUAL power spectra BASELINED
close all
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
    pow_lc_subj = power_lc_baselined{subj};
    pow_hc_subj = power_hc_baselined{subj};

    cfg = [];
    cfg.channel = channels;  % Focus on occipital channels
    cfg.figure = 'gcf';
    cfg.linewidth = 1;

    % Plot power spectrum for low and high contrast
    ft_singleplotER(cfg, pow_lc_subj, pow_hc_subj);
    hold on;

    % Add shaded error bars
    channels_seb = ismember(pow_lc_subj.label, cfg.channel);
    lceb = shadedErrorBar(pow_lc_subj.freq, mean(pow_lc_subj.powspctrm(channels_seb, :), 1), ...
        std(pow_lc_subj.powspctrm(channels_seb, :)) / sqrt(size(pow_lc_subj.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    hceb = shadedErrorBar(pow_hc_subj.freq, mean(pow_hc_subj.powspctrm(channels_seb, :), 1), ...
        std(pow_hc_subj.powspctrm(channels_seb, :)) / sqrt(size(pow_hc_subj.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    lceb.mainLine.Color = colors(1, :);
    hceb.mainLine.Color = colors(2, :);
    lceb.patch.FaceColor = colors(1, :);
    hceb.patch.FaceColor = colors(2, :);
    set(lceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
    set(hceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
    set(lceb.edge(1), 'Color', colors(1, :));
    set(lceb.edge(2), 'Color', colors(1, :));
    set(hceb.edge(1), 'Color', colors(2, :));
    set(hceb.edge(2), 'Color', colors(2, :));
    set(lceb.patch, 'FaceAlpha', 0.5);
    set(hceb.patch, 'FaceAlpha', 0.5);

     % Find channels and frequencies of interest
    channels_idx = ismember(pow_lc_subj.label, channels);
    freq_idx = find(pow_lc_subj.freq >= 30 & pow_hc_subj.freq <= 90);

    % Find gamma peak for LOW contrast
    lc_gamma_power = mean(pow_lc_subj.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(lc_gamma_power, pow_lc_subj.freq(freq_idx));
    [lc_pow, peak_idx] = max(peaks);
    lc_freq = locs(peak_idx);

    % Find gamma peak for HIGH contrast
    hc_gamma_power = mean(pow_hc_subj.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(hc_gamma_power, pow_hc_subj.freq(freq_idx));
    [hc_pow, peak_idx] = max(peaks);
    hc_freq = locs(peak_idx);
    % if subj == 9
    %     hc_pow = 0.1149
    %     hc_freq = 86
    % end

    % Adjust plot aesthetics
    set(gca, 'FontSize', 20);
    yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
    plot([0 lc_freq], [lc_pow lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
    plot([lc_freq lc_freq], [-100 lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
    plot([0 hc_freq], [hc_pow hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
    plot([hc_freq hc_freq], [-100 hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
    max_spctrm = max(lc_pow, hc_pow);
    ylim([-max_spctrm*1.25 max_spctrm*1.25]);
    % if subj == 9
    %     ylim([-0.35 0.35])
    % elseif subj == 10
    %     ylim([-0.75 0.75])
    % end
    xlim([30 90]);
    xticks(30:10:90);
    xlabel('Freq [Hz]', 'FontSize', 15);
    ylabel('Power [dB]', 'FontSize', 15);
    if mod(subj, 5) == 0
        legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
    end
    title(sprintf('Subject %s', subjects{subj}), 'FontSize', 12);
    hold off;
end

% Save the combined figure with all subplots
if anal_period == 1
    save_path = fullfile(output_dir, 'GCP_powspctrm_all_subjects_subplot_baselined_300.png');
else
    save_path = fullfile(output_dir, 'GCP_powspctrm_all_subjects_subplot_baselined.png');
end
saveas(gcf, save_path);