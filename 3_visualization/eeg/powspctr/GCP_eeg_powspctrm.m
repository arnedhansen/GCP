%% GCP Gamma Peak Power and Frequency

%% Setup
clear
[subjects, path, colors] = setup('GCP');
analysis_period = 0; % 1 = ONLY 0-300ms, otherwise 300-2000ms after stimulus presentation

%% Load power spectra data
for subj = 1:length(subjects)
    % if analysis_period == 1
    %     load(strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra_300'))
    % else
    %     load(strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra'))
    % end
    load(strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra_superlets'))

    % Raw powerspectra
    power_c25{subj}  = pow_c25;
    power_c50{subj}  = pow_c50;
    power_c75{subj}  = pow_c75;
    power_c100{subj} = pow_c100;

    % Baselined powerspectra
    power_c25_baselined{subj}  = pow_c25_baselined;
    power_c50_baselined{subj}  = pow_c50_baselined;
    power_c75_baselined{subj}  = pow_c75_baselined;
    power_c100_baselined{subj} = pow_c100_baselined;

    % Baseline period powerspectra
    power_c25_baseline_period{subj}  = pow_c25_baseline_period;
    power_c50_baseline_period{subj}  = pow_c50_baseline_period;
    power_c75_baseline_period{subj}  = pow_c75_baseline_period;
    power_c100_baseline_period{subj} = pow_c100_baseline_period;

    % FOOOFed powerspectra
    power_c25_fooof{subj}  = pow_c25_fooof;
    power_c50_fooof{subj}  = pow_c50_fooof;
    power_c75_fooof{subj}  = pow_c75_fooof;
    power_c100_fooof{subj} = pow_c100_fooof;

    % Baselined FOOOFed smoothed powerspectra
    power_c25_fooof_bl_smooth{subj}  = pow_c25_fooof_bl_smooth;
    power_c50_fooof_bl_smooth{subj}  = pow_c50_fooof_bl_smooth;
    power_c75_fooof_bl_smooth{subj}  = pow_c75_fooof_bl_smooth;
    power_c100_fooof_bl_smooth{subj} = pow_c100_fooof_bl_smooth;
end

% Compute grand averages
cfg = [];
%cfg.keepindividual = 'yes';

gapow_c25  = ft_freqgrandaverage(cfg, power_c25{:});
gapow_c50  = ft_freqgrandaverage(cfg, power_c50{:});
gapow_c75  = ft_freqgrandaverage(cfg, power_c75{:});
gapow_c100 = ft_freqgrandaverage(cfg, power_c100{:});

gapow_c25_baselined  = ft_freqgrandaverage(cfg, power_c25_baselined{:});
gapow_c50_baselined  = ft_freqgrandaverage(cfg, power_c50_baselined{:});
gapow_c75_baselined  = ft_freqgrandaverage(cfg, power_c75_baselined{:});
gapow_c100_baselined = ft_freqgrandaverage(cfg, power_c100_baselined{:});

gapow_c25_baseline_period  = ft_freqgrandaverage(cfg, power_c25_baseline_period{:});
gapow_c50_baseline_period  = ft_freqgrandaverage(cfg, power_c50_baseline_period{:});
gapow_c75_baseline_period  = ft_freqgrandaverage(cfg, power_c75_baseline_period{:});
gapow_c100_baseline_period = ft_freqgrandaverage(cfg, power_c100_baseline_period{:});

gapow_c25_fooof  = ft_freqgrandaverage(cfg, power_c25_fooof{:});
gapow_c50_fooof  = ft_freqgrandaverage(cfg, power_c50_fooof{:});
gapow_c75_fooof  = ft_freqgrandaverage(cfg, power_c75_fooof{:});
gapow_c100_fooof = ft_freqgrandaverage(cfg, power_c100_fooof{:});

gapow_c25_fooof_bl_smooth  = ft_freqgrandaverage(cfg, power_c25_fooof_bl_smooth{:});
gapow_c50_fooof_bl_smooth  = ft_freqgrandaverage(cfg, power_c50_fooof_bl_smooth{:});
gapow_c75_fooof_bl_smooth  = ft_freqgrandaverage(cfg, power_c75_fooof_bl_smooth{:});
gapow_c100_fooof_bl_smooth = ft_freqgrandaverage(cfg, power_c100_fooof_bl_smooth{:});

%% Define channels
datapath = strcat(path, subjects{1}, '/eeg');
cd(datapath);
% Occipital channels
occ_channels = {};
pow_label = pow_c25; % Assume similar structure across conditions
for i = 1:length(pow_label.label)
    label = pow_label.label{i};
    if contains(label, {'O'}) || contains(label, {'P'}) && ~contains(label, {'T'}) ...
        && ~contains(label, {'C'}) || contains(label, {'I'}) %|| contains(label, {'CPP'}) ...
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;
% Test channels
%channels = [{'Pz'}, {'P1'}, {'P2'}, {'P3'}, {'P4'}, {'P5'}, {'P6'}, {'PPO1'}, {'PPO2'}, {'PPO5h'}, {'PPO6h'} {'PO3'}, {'PO4'}, {'POz'}]

%% Plot GRAND AVERAGE power spectrum
close all;
figure;
set(gcf, 'Position', [0, 0, 1000, 2000], 'Color', 'w');
cfg = [];
cfg.channel = channels;
cfg.figure = 'gcf';
cfg.linewidth = 4;
channels_seb = ismember(gapow_c25_fooof_bl_smooth.label, cfg.channel);

% %OPTIONAL: Add extra smoothing to powerspectra
% smoothWin = 20; % Size of smoothing window in frequency bins
% gapow_c25_fooof_bl_smooth_extra  = gapow_c25_fooof_bl_smooth;
% gapow_c50_fooof_bl_smooth_extra  = gapow_c50_fooof_bl_smooth;
% gapow_c75_fooof_bl_smooth_extra  = gapow_c75_fooof_bl_smooth;
% gapow_c100_fooof_bl_smooth_extra = gapow_c100_fooof_bl_smooth;
 % Moving Window
% gapow_c25_fooof_bl_smooth_extra.powspctrm(channels_seb, :)  = movmean(gapow_c25_fooof_bl_smooth.powspctrm(channels_seb, :), smoothWin, 2);
% gapow_c50_fooof_bl_smooth_extra.powspctrm(channels_seb, :)  = movmean(gapow_c50_fooof_bl_smooth.powspctrm(channels_seb, :), smoothWin, 2);
% gapow_c75_fooof_bl_smooth_extra.powspctrm(channels_seb, :)  = movmean(gapow_c75_fooof_bl_smooth.powspctrm(channels_seb, :), smoothWin, 2);
% gapow_c100_fooof_bl_smooth_extra.powspctrm(channels_seb, :) = movmean(gapow_c100_fooof_bl_smooth.powspctrm(channels_seb, :), smoothWin, 2);
%  % Gaussian Kernel
% gapow_c25_fooof_bl_smooth_extra.powspctrm(channels_seb, :)  = smoothdata(gapow_c25_fooof_bl_smooth.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);
% gapow_c50_fooof_bl_smooth_extra.powspctrm(channels_seb, :)  = smoothdata(gapow_c50_fooof_bl_smooth.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);
% gapow_c75_fooof_bl_smooth_extra.powspctrm(channels_seb, :)  = smoothdata(gapow_c75_fooof_bl_smooth.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);
% gapow_c100_fooof_bl_smooth_extra.powspctrm(channels_seb, :) = smoothdata(gapow_c100_fooof_bl_smooth.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);

% %OPTIONAL: Add extra smoothing to powerspectra
smoothWin = 10; % Size of smoothing window in frequency bins
gapow_c25_fooof.powspctrm  = smoothdata(gapow_c25_fooof.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);
gapow_c50_fooof.powspctrm  = smoothdata(gapow_c50_fooof.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);
gapow_c75_fooof.powspctrm  = smoothdata(gapow_c75_fooof.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);
gapow_c100_fooof.powspctrm = smoothdata(gapow_c100_fooof.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);


% Set data to plot
ga25=  gapow_c25_fooof%_bl_smooth_extra ;
ga50=  gapow_c50_fooof%_bl_smooth_extra ;
ga75=  gapow_c75_fooof%_bl_smooth_extra ;
ga100=gapow_c100_fooof%_bl_smooth_extra;


% Plot for all contrasts
ft_singleplotER(cfg, ga25, ga50, ga75, ga100);
hold on;

% Add shaded error bars for each condition
conditions = {ga25, ga50, ga75, ga100};
col_indices = [1, 2, 3, 4];
for i = 1:4
    curr_cond = conditions{i};
    lp     = {'-','Color', colors(col_indices(i),:)};
    seb = shadedErrorBar(curr_cond.freq, ...
        mean(curr_cond.powspctrm(channels_seb, :), 1), ...
        std(curr_cond.powspctrm(channels_seb, :)) / sqrt(size(curr_cond.powspctrm(channels_seb, :), 1)), ...
        'lineProps', lp);
    h(i) = seb.mainLine;
    seb.mainLine.Color = colors(col_indices(i), :);
    seb.patch.FaceColor = colors(col_indices(i), :);
    set(seb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(col_indices(i), :));
    set(seb.patch, 'FaceAlpha', 0.25);
end

% Adjust plot aesthetics
yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
set(gcf, 'color', 'w');
set(gca, 'FontSize', 20);
% xlim([30 90])
%ylim([-0.04 0.04])
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Power [dB]', 'FontSize', 30);
legend(h, {' 25% Contrast',' 50% Contrast',' 75% Contrast','100% Contrast'}, 'FontName','Arial','FontSize',25);
title('Grand Average Power Spectrum', 'FontSize', 40);
hold off;

% Save the plot
if analysis_period == 1
    saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_powspctrm_fooof_bl_smooth_300.png');
else
    saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_powspctrm_fooof_bl_smooth.png');
end

%% Plot GRAND AVERAGE power spectrum PERCENTAGE CHANGE
% Percentage change formula: ((stimulus - baseline) / baseline) * 100
% percent_change_LOW = gapow_lc;
% percent_change_HIGH = gapow_hc;
% percent_change_LOW.powspctrm = ((gapow_lc.powspctrm - gapow_lc_baseline_period.powspctrm) ./ gapow_lc_baseline_period.powspctrm) * 100;
% percent_change_HIGH.powspctrm = ((gapow_hc.powspctrm - gapow_hc_baseline_period.powspctrm) ./ gapow_hc_baseline_period.powspctrm) * 100;
% 
% close all;
% figure;
% set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
% cfg = [];
% cfg.channel = channels;
% cfg.figure = 'gcf';
% cfg.linewidth = 1;
% 
% % Plot for low and high contrast
% ft_singleplotER(cfg, percent_change_LOW, percent_change_HIGH);
% hold on;
% 
% % Add shaded error bars
% channels_seb = ismember(percent_change_LOW.label, cfg.channel);
% lceb = shadedErrorBar(percent_change_LOW.freq, mean(percent_change_LOW.powspctrm(channels_seb, :), 1), std(percent_change_LOW.powspctrm(channels_seb, :))/sqrt(size(percent_change_LOW.powspctrm(channels_seb, :), 1)), {'-'}, 0);
% hceb = shadedErrorBar(percent_change_HIGH.freq, mean(percent_change_HIGH.powspctrm(channels_seb, :), 1), std(percent_change_HIGH.powspctrm(channels_seb, :))/sqrt(size(percent_change_HIGH.powspctrm(channels_seb, :), 1)), {'-'}, 0);
% lceb.mainLine.Color = colors(1, :);
% hceb.mainLine.Color = colors(2, :);
% lceb.patch.FaceColor = colors(1, :);
% hceb.patch.FaceColor = colors(2, :);
% set(lceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
% set(hceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
% set(lceb.edge(1), 'Color', colors(1, :));
% set(lceb.edge(2), 'Color', colors(1, :));
% set(hceb.edge(1), 'Color', colors(2, :));
% set(hceb.edge(2), 'Color', colors(2, :));
% set(lceb.patch, 'FaceAlpha', 0.5);
% set(hceb.patch, 'FaceAlpha', 0.5);
% 
% % Find GA gamma peak for LOW contrast
% lc_gamma_power = mean(percent_change_LOW.powspctrm(channels_seb, freq_idx), 1);
% [peaks, locs] = findpeaks(lc_gamma_power, percent_change_LOW.freq(freq_idx));
% [lc_pow, peak_idx] = max(peaks);
% lc_freq = locs(peak_idx);
% 
% % Find GA gamma peak for HIGH contrast
% hc_gamma_power = mean(percent_change_HIGH.powspctrm(channels_seb, freq_idx), 1);
% [peaks, locs] = findpeaks(hc_gamma_power, percent_change_HIGH.freq(freq_idx));
% [hc_pow, peak_idx] = max(peaks);
% hc_freq = locs(peak_idx);
% 
% % Adjust plot aesthetics
% yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
% plot([0 lc_freq], [lc_pow lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
% plot([lc_freq lc_freq], [-100 lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
% plot([0 hc_freq], [hc_pow hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
% plot([hc_freq hc_freq], [-100 hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
% set(gcf,'color','w');
% set(gca,'Fontsize',20);
% max_spctrm = max(lc_pow, hc_pow);
% ylim([-max_spctrm*1.25 max_spctrm*1.25]);
% xlim([30 90]) % Gamma frequency range
% xlabel('Frequency [Hz]');
% ylabel('Percentage Change [%]');
% legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
% title('Percentage Change Power Spectrum', 'FontSize', 30);
% hold off;
% 
% % Save the plot
% if analysis_period == 1
%     saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/powspctrm_analysis_period_300/GCP_powspctrm_percentage_300.png');
% else
%     saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_powspctrm_percentage.png');
% end

%% Plot INDIVIDUAL powerspectra
% output_dir = '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/';
% for subj = 1:length(subjects)
%     close all;
%     figure;
%     set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
% 
%     % Extract participant data
%     pow_lc_subj = power_lc_baselined{subj};
%     pow_hc_subj = power_hc_baselined{subj};
% 
%     % Figure common config
%     cfg = [];
%     cfg.channel = channels;
%     cfg.figure = 'gcf';
%     cfg.linewidth = 1;
% 
%     % Plot power spectrum for low and high contrast
%     ft_singleplotER(cfg, pow_lc_subj, pow_hc_subj);
%     hold on;
% 
%     % Add shaded error bars
%     channels_seb = ismember(pow_lc_subj.label, cfg.channel);
%     lceb = shadedErrorBar(pow_lc_subj.freq, mean(pow_lc_subj.powspctrm(channels_seb, :), 1), ...
%         std(pow_lc_subj.powspctrm(channels_seb, :)) / sqrt(size(pow_lc_subj.powspctrm(channels_seb, :), 1)), {'-'}, 0);
%     hceb = shadedErrorBar(pow_hc_subj.freq, mean(pow_hc_subj.powspctrm(channels_seb, :), 1), ...
%         std(pow_hc_subj.powspctrm(channels_seb, :)) / sqrt(size(pow_hc_subj.powspctrm(channels_seb, :), 1)), {'-'}, 0);
%     lceb.mainLine.Color = colors(1, :);
%     hceb.mainLine.Color = colors(2, :);
%     lceb.patch.FaceColor = colors(1, :);
%     hceb.patch.FaceColor = colors(2, :);
%     set(lceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
%     set(hceb.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
%     set(lceb.edge(1), 'Color', colors(1, :));
%     set(lceb.edge(2), 'Color', colors(1, :));
%     set(hceb.edge(1), 'Color', colors(2, :));
%     set(hceb.edge(2), 'Color', colors(2, :));
%     set(lceb.patch, 'FaceAlpha', 0.5);
%     set(hceb.patch, 'FaceAlpha', 0.5);
% 
%     % Find channels and frequencies of interest
%     channels_idx = ismember(pow_lc_subj.label, channels);
%     freq_idx = find(pow_lc_subj.freq >= 30 & pow_hc_subj.freq <= 90);
% 
%     % Find gamma peak for LOW contrast
%     lc_gamma_power = mean(pow_lc_subj.powspctrm(channels_idx, freq_idx), 1);
%     [peaks, locs] = findpeaks(lc_gamma_power, pow_lc_subj.freq(freq_idx));
%     [lc_pow, peak_idx] = max(peaks);
%     lc_freq = locs(peak_idx);
% 
%     % Find gamma peak for HIGH contrast
%     hc_gamma_power = mean(pow_hc_subj.powspctrm(channels_idx, freq_idx), 1);
%     [peaks, locs] = findpeaks(hc_gamma_power, pow_hc_subj.freq(freq_idx));
%     [hc_pow, peak_idx] = max(peaks);
%     hc_freq = locs(peak_idx);
%     % if subj == 9
%     %     hc_pow = 0.1149
%     %     hc_freq = 86
%     % end
% 
%     % Adjust plot aesthetics
%     yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
%     plot([0 lc_freq], [lc_pow lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
%     plot([lc_freq lc_freq], [-100 lc_pow], '--', 'Color', colors(1, :), 'LineWidth', 2);
%     plot([0 hc_freq], [hc_pow hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
%     plot([hc_freq hc_freq], [-100 hc_pow], '--', 'Color', colors(2, :), 'LineWidth', 2);
%     set(gca, 'FontSize', 20);
%     max_spctrm = max(lc_pow, hc_pow);
%     % if subj == 9
%     %     max_spctrm = 0.4
%     % elseif subj == 10
%     %     max_spctrm = 0.55
%     % end
%     ylim([-max_spctrm*1.25 max_spctrm*1.25]);
%     xlim([30 90]);
%     xlabel('Frequency [Hz]');
%     ylabel('Power [dB]');
%     legend([lceb.mainLine, hceb.mainLine], {'Low Contrast', 'High Contrast'}, 'FontName', 'Arial', 'FontSize', 20);
%     title(sprintf('Subject %s: Power Spectrum', subjects{subj}), 'FontSize', 30);
%     hold off;
% 
%     % Save individual plot
%     if analysis_period == 1
%         output_dir = '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/powspctrm_analysis_period_300';
%         save_path = fullfile(output_dir, sprintf('GCP_powspctrm_subj%s_baselined_300.png', subjects{subj}));
%     else
%         save_path = fullfile(output_dir, sprintf('GCP_powspctrm_subj%s_baselined.png', subjects{subj}));
%     end
%     saveas(gcf, save_path);
% end

%% Subplot with all INDIVIDUAL power spectra
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

    % Extract participant data for all four conditions
    pow_c25_subj  = power_c25_fooof_bl_smooth{subj};
    pow_c50_subj  = power_c50_fooof_bl_smooth{subj};
    pow_c75_subj  = power_c75_fooof_bl_smooth{subj};
    pow_c100_subj = power_c100_fooof_bl_smooth{subj};

    % OPTIONAL: Add extra smoothing with moving window
    channels_seb = ismember(pow_c25_subj.label, cfg.channel);
    smoothWin = 10; % Size of smoothing window in frequency bins    smoothWin = 10
    pow_c25_subj.powspctrm(channels_seb, :)  = movmean(pow_c25_subj.powspctrm(channels_seb, :), smoothWin, 2);
    pow_c50_subj.powspctrm(channels_seb, :)  = movmean(pow_c50_subj.powspctrm(channels_seb, :), smoothWin, 2);
    pow_c75_subj.powspctrm(channels_seb, :)  = movmean(pow_c75_subj.powspctrm(channels_seb, :), smoothWin, 2);
    pow_c100_subj.powspctrm(channels_seb, :) = movmean(pow_c100_subj.powspctrm(channels_seb, :), smoothWin, 2);

    cfg = [];
    cfg.channel = channels;
    cfg.figure = 'gcf';
    cfg.linewidth = 3;

    % Plot power spectrum for all four conditions
    ft_singleplotER(cfg, pow_c25_subj, pow_c50_subj, pow_c75_subj, pow_c100_subj);
    hold on;

    % Define colors for each condition
    condition_colors = colors(1:4, :);

    % Add shaded error bars for each condition
    conditions = {pow_c25_subj, pow_c50_subj, pow_c75_subj, pow_c100_subj};
    shadedEBs = cell(1, 4);
    for i = 1:4
        channels_seb = ismember(conditions{i}.label, cfg.channel);
        shadedEBs{i} = shadedErrorBar(conditions{i}.freq, mean(conditions{i}.powspctrm(channels_seb, :), 1), ...
            std(conditions{i}.powspctrm(channels_seb, :)) / sqrt(size(conditions{i}.powspctrm(channels_seb, :), 1)));
        shadedEBs{i}.mainLine.Color = condition_colors(i, :);
        shadedEBs{i}.patch.FaceColor = condition_colors(i, :);
        set(shadedEBs{i}.mainLine, 'LineWidth', cfg.linewidth);
        set(shadedEBs{i}.patch, 'FaceAlpha', 0.5);
    end

    % Find channels and frequencies of interest
    channels_idx = ismember(pow_c25_subj.label, channels);
    freq_idx = find(pow_c25_subj.freq >= 30 & pow_c100_subj.freq <= 90);

    peak_freqs = zeros(1, 4);
    peak_powers = zeros(1, 4);
    for i = 1:4
        gamma_power = mean(conditions{i}.powspctrm(channels_idx, freq_idx), 1);
        [peaks, locs] = findpeaks(gamma_power, conditions{i}.freq(freq_idx));
        if ~isempty(peaks)
            [peak_powers(i), peak_idx] = max(peaks);
            peak_freqs(i) = locs(peak_idx);
        end
    end

    % Adjust plot aesthetics
    set(gca, 'FontSize', 20);
    yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
    for i = 1:4
        if peak_powers(i) > 0
            plot([0 peak_freqs(i)], [peak_powers(i) peak_powers(i)], '--', 'Color', condition_colors(i, :), 'LineWidth', 2);
            plot([peak_freqs(i) peak_freqs(i)], [-100 peak_powers(i)], '--', 'Color', condition_colors(i, :), 'LineWidth', 2);
        end
    end

    max_spctrm = max(peak_powers);
    set(gca, "YLim", [-max_spctrm*1.25 max_spctrm*1.25])
    xlim([30 90]);
    xticks(30:10:90);
    xlabel('Freq [Hz]', 'FontSize', 20);
    ylabel('Power [dB]', 'FontSize', 20);
    if subj == 5
        legend([shadedEBs{1}.mainLine, shadedEBs{2}.mainLine, shadedEBs{3}.mainLine, shadedEBs{4}.mainLine], ...
            {'C25', 'C50', 'C75', 'C100'}, 'FontName', 'Arial', 'FontSize', 20, 'Location', 'best');
    end
    title(sprintf('Subject %d', subj), 'FontSize', 20);
    hold off;
end

% Save the combined figure with all subplots
if analysis_period == 1
    save_path = fullfile(output_dir, 'GCP_powspctrm_all_subjects_subplot_300.png');
else
    save_path = fullfile(output_dir, 'GCP_powspctrm_all_subjects_subplot.png');
end
saveas(gcf, save_path);
