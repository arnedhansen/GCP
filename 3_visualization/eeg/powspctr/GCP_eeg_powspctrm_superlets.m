%% GCP Gamma Peak Power and Frequency SUPERLETS (Percentage Change Version)

%% Setup
clear
[subjects, path, colors] = setup('GCP');
analysis_period = 0; % 1 = ONLY 0–300 ms, otherwise 300–2000 ms after stimulus presentation

%% Load power spectra data
for subj = 1:length(subjects)
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

    % Re-adjust frequency smoothing
    orig_freq = 30:5:90; % 13 values
    new_freq  = 30:1:90; % 61 values
    power_c25_baseline_period{subj}  = smooth_tfr(power_c25_baseline_period{subj}, orig_freq, new_freq);
    power_c50_baseline_period{subj}  = smooth_tfr(power_c50_baseline_period{subj}, orig_freq, new_freq);
    power_c75_baseline_period{subj}  = smooth_tfr(power_c75_baseline_period{subj}, orig_freq, new_freq);
    power_c100_baseline_period{subj} = smooth_tfr(power_c100_baseline_period{subj}, orig_freq, new_freq);

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
gapow_c25_fooof_bl_smooth  = ft_freqgrandaverage(cfg, power_c25_fooof_bl_smooth{:});
gapow_c50_fooof_bl_smooth  = ft_freqgrandaverage(cfg, power_c50_fooof_bl_smooth{:});
gapow_c75_fooof_bl_smooth  = ft_freqgrandaverage(cfg, power_c75_fooof_bl_smooth{:});
gapow_c100_fooof_bl_smooth = ft_freqgrandaverage(cfg, power_c100_fooof_bl_smooth{:});

gapow_c25_baseline_period  = ft_freqgrandaverage(cfg, power_c25_baseline_period{:});
gapow_c50_baseline_period  = ft_freqgrandaverage(cfg, power_c50_baseline_period{:});
gapow_c75_baseline_period  = ft_freqgrandaverage(cfg, power_c75_baseline_period{:});
gapow_c100_baseline_period = ft_freqgrandaverage(cfg, power_c100_baseline_period{:});

%% Define channels
datapath = strcat(path, subjects{1}, '/eeg');
cd(datapath);
occ_channels = {};
pow_label = pow_c25; 
for i = 1:length(pow_label.label)
    label = pow_label.label{i};
    if (contains(label, 'O') || contains(label, 'P')) && ...
       ~contains(label, 'T') && ~contains(label, 'C') || ...
       contains(label, 'I')
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Figure 1: Grand average percentage change spectrum
% Compute percentage change grand averages
pct_c25 = gapow_c25_fooof_bl_smooth;
pct_c50 = gapow_c50_fooof_bl_smooth;
pct_c75 = gapow_c75_fooof_bl_smooth;
pct_c100 = gapow_c100_fooof_bl_smooth;

pct_c25.powspctrm  = (gapow_c25_fooof_bl_smooth.powspctrm  - gapow_c25_baseline_period.powspctrm ) ...
                     ./ gapow_c25_baseline_period.powspctrm  * 100;
pct_c50.powspctrm  = (gapow_c50_fooof_bl_smooth.powspctrm  - gapow_c50_baseline_period.powspctrm ) ...
                     ./ gapow_c50_baseline_period.powspctrm  * 100;
pct_c75.powspctrm  = (gapow_c75_fooof_bl_smooth.powspctrm  - gapow_c75_baseline_period.powspctrm ) ...
                     ./ gapow_c75_baseline_period.powspctrm  * 100;
pct_c100.powspctrm = (gapow_c100_fooof_bl_smooth.powspctrm - gapow_c100_baseline_period.powspctrm) ...
                     ./ gapow_c100_baseline_period.powspctrm * 100;

% %OPTIONAL: Add extra smoothing to powerspectra
smoothWin = 20; % Size of smoothing window in frequency bins
pct_c25_smooth_extra  = pct_c25;
pct_c50_smooth_extra  = pct_c50;
pct_c75_smooth_extra  = pct_c75;
pct_c100_smooth_extra = pct_c100;
 % Moving Window
% pct_c25_smooth_extra.powspctrm  = movmean(pct_c25_smooth_extra.powspctrm, smoothWin, 2);
% pct_c50_smooth_extra.powspctrm  = movmean(pct_c50_smooth_extra.powspctrm, smoothWin, 2);
% pct_c75_smooth_extra.powspctrm  = movmean(pct_c75_smooth_extra.powspctrm, smoothWin, 2);
% pct_c100_smooth_extra.powspctrm = movmean(pct_c100_smooth_extra.powspctrm, smoothWin, 2);
%  % Gaussian Kernel
% pct_c25_smooth_extra.powspctrm  = smoothdata(gapow_c25_fooof_bl_smooth.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);
% pct_c50_smooth_extra.powspctrm  = smoothdata(gapow_c50_fooof_bl_smooth.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);
% pct_c75_smooth_extra.powspctrm  = smoothdata(gapow_c75_fooof_bl_smooth.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);
% pct_c100_smooth_extra.powspctrm = smoothdata(gapow_c100_fooof_bl_smooth.powspctrm(channels_seb, :), 2, 'gaussian', smoothWin);

% Plot grand average percentage change
figure;
set(gcf, 'Position', [0, 0, 1000, 2000], 'Color', 'w');
cfg = [];
cfg.channel  = channels;
cfg.figure  = 'gcf';
cfg.linewidth = 4;
channels_seb = ismember(pct_c25.label, cfg.channel);

ft_singleplotER(cfg, pct_c25_smooth_extra, pct_c50_smooth_extra, pct_c75_smooth_extra, pct_c100_smooth_extra);
hold on;

conditions = {pct_c25_smooth_extra, pct_c50_smooth_extra, pct_c75_smooth_extra, pct_c100_smooth_extra};
col_indices = [1, 2, 3, 4];
h = gobjects(1,4);
for i = 1:4
    curr = conditions{i};
    lp = {'-','Color', colors(col_indices(i),:)};
    seb = shadedErrorBar(curr.freq, ...
        mean(curr.powspctrm(channels_seb,:),1), ...
        std(curr.powspctrm(channels_seb,:),[],1) ./ sqrt(sum(channels_seb)), ...
        'lineProps', lp);
    h(i) = seb.mainLine;
    seb.mainLine.Color      = colors(col_indices(i), :);
    seb.patch.FaceColor     = colors(col_indices(i), :);
    set(seb.mainLine, 'LineWidth', cfg.linewidth);
    set(seb.patch, 'FaceAlpha', 0.25);
end

yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
set(gca, 'FontSize', 20);
xlabel('Frequency [Hz]', 'FontSize', 30);
ylabel('Percentage Change [%]', 'FontSize', 30);
legend(h, {' 25% Contrast',' 50% Contrast',' 75% Contrast','100% Contrast'}, ...
       'FontName','Arial','FontSize',25);
title('Grand Average Percentage Change Power Spectrum', 'FontSize', 40);
hold off;

% Save the plot
if analysis_period == 1
    saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_pct_powspctrm_fooof_bl_smooth_300.png');
else
    saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/GCP_pct_powspctrm_fooof_bl_smooth.png');
end

%% Figure 2: Subplots of individual percentage change spectra
output_dir = '/Volumes/methlab/Students/Arne/GCP/figures/eeg/powspctrm/';
num_subs = length(subjects);
cols = 5;
rows = ceil(num_subs / cols);

figure;
set(gcf, 'Position', [0, 0, 1600, 1200], 'Color', 'w');

for subj = 1:num_subs
    subplot(rows, cols, subj);

    % Load subject data
    fooof25 = power_c25_fooof_bl_smooth{subj};
    fooof50 = power_c50_fooof_bl_smooth{subj};
    fooof75 = power_c75_fooof_bl_smooth{subj};
    fooof100= power_c100_fooof_bl_smooth{subj};

    base25 = power_c25_baseline_period{subj}.powspctrm;
    base50 = power_c50_baseline_period{subj}.powspctrm;
    base75 = power_c75_baseline_period{subj}.powspctrm;
    base100= power_c100_baseline_period{subj}.powspctrm;

    % Compute percentage change for this subject
    fooof25.powspctrm  = (fooof25.powspctrm  - base25)  ./ base25  * 100;
    fooof50.powspctrm  = (fooof50.powspctrm  - base50)  ./ base50  * 100;
    fooof75.powspctrm  = (fooof75.powspctrm  - base75)  ./ base75  * 100;
    fooof100.powspctrm = (fooof100.powspctrm - base100) ./ base100 * 100;

    cfg = [];
    cfg.channel   = channels;
    cfg.figure   = 'gcf';
    cfg.linewidth = 3;

    ft_singleplotER(cfg, fooof25, fooof50, fooof75, fooof100);
    hold on;

    conds = {fooof25, fooof50, fooof75, fooof100};
    for i = 1:4
        seb = shadedErrorBar(conds{i}.freq, ...
            mean(conds{i}.powspctrm(ismember(conds{i}.label,channels),:),1), ...
            std(conds{i}.powspctrm(ismember(conds{i}.label,channels),:),[],1) ...
            ./ sqrt(sum(ismember(conds{i}.label,channels))));
        seb.mainLine.Color  = colors(i,:);
        seb.patch.FaceColor = colors(i,:);
        set(seb.mainLine, 'LineWidth', cfg.linewidth);
        set(seb.patch, 'FaceAlpha', 0.5);
    end

    set(gca, 'FontSize', 16);
    yline(0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
    xlim([30 90]);
    xticks(30:10:90);
    xlabel('Freq [Hz]', 'FontSize', 12);
    ylabel('Pct Change [%]', 'FontSize', 12);
    if subj == 5
        legend('C25','C50','C75','C100','Location','best');
    end
    title(sprintf('Subject %d', subj), 'FontSize', 14);
    hold off;
end

% Save combined figure
if analysis_period == 1
    saveas(gcf, fullfile(output_dir, 'GCP_pct_powspctrm_all_subjects_subplot_300.png'));
else
    saveas(gcf, fullfile(output_dir, 'GCP_pct_powspctrm_all_subjects_subplot.png'));
end
