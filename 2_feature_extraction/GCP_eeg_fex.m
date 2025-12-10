%% GCP EEG Feature Extraction
%
% Extracted features:
%   TFR
%       Raw
%       FOOOF
%       Baselined
%   Power Spectrum
%   Gamma Peak Power
%   Gamma Peak Frequency

%% Setup
startup
[subjects, path, ~, ~] = setup('GCP');

%% Extract TFR
% Read data, segment and convert to FieldTrip data structure
for subj = 1 : length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    %if ~isfile(strcat([datapath, '/data_tfr.mat'])) % only new data
    cd(datapath)
    close all
    load dataEEG
    load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions
    ind61 = find(dataEEG_c25.trialinfo  == 61);
    ind62 = find(dataEEG_c50.trialinfo  == 62);
    ind63 = find(dataEEG_c75.trialinfo  == 63);
    ind64 = find(dataEEG_c100.trialinfo == 64);

    %% Time frequency analysis
    cfg             = [];
    cfg.output      = 'pow';
    cfg.method      = 'mtmconvol';
    cfg.taper       = 'dpss';
    cfg.foi         = 30:5:120;                         % Analysis 30 to 120 Hz in steps of 5 Hz
    cfg.tapsmofrq   = 5;                                % Analysis of FOI +/- 5 Hz
    cfg.t_ftimwin   = ones(length(cfg.foi),1).*0.5;     % Length of time window = 0.5 sec
    cfg.toi          = -1.75:0.05:2;                 % Time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
    cfg.keeptrials  = 'no';

    % 25% contrast concentric dynamic inward
    cfg.trials = ind61;
    tfr_c25 = ft_freqanalysis(cfg,dataEEG_c25);

    % 50% contrast concentric dynamic inward
    cfg.trials = ind62;
    tfr_c50 = ft_freqanalysis(cfg,dataEEG_c50);

    % 75% contrast concentric dynamic inward
    cfg.trials = ind63;
    tfr_c75 = ft_freqanalysis(cfg,dataEEG_c75);

    % 100% contrast concentric dynamic inward
    cfg.trials = ind64;
    tfr_c100 = ft_freqanalysis(cfg,dataEEG_c100);

    %% FOOOF
    orig_freq = 30:5:120;
    tfrs = {tfr_c25, tfr_c50, tfr_c75, tfr_c100};
    for tfr_contrast = 1:4
        clc
        disp('FOOOFing...')
        clear fspctrm
        tfr = tfrs{1, tfr_contrast};
        for t = 1 :length(tfr.time)
            cfg = [];
            cfg.latency = tfr.time(t);
            tmp = ft_selectdata(cfg,tfr);

            for chan = 1:length(tmp.label)

                % Transpose, to make inputs row vectors
                freqs = tmp.freq';
                psd = tmp.powspctrm(chan,:)';

                % FOOOF settings
                settings = struct();  % Use defaults
                settings.verbose = false; % Suppress warnings about too low peak_width_limits
                f_range = [tfr.freq(1), tfr.freq(end)];

                % Run FOOOF
                freqs = orig_freq'; % Equidistant freq distribution
                fooof_results = fooof(freqs, psd, f_range, settings, true);
                powspctrmff(chan,:) = fooof_results.fooofed_spectrum-fooof_results.ap_fit;
            end
            fspctrm(:,:,t) = powspctrmff;
        end
        fooofedtrl(:,:,:) = fspctrm;
        if tfr_contrast == 1
            tfr_c25_fooof = tfr;
            tfr_c25_fooof.powspctrm = fooofedtrl;
        elseif tfr_contrast == 2
            tfr_c50_fooof = tfr;
            tfr_c50_fooof.powspctrm = fooofedtrl;
        elseif tfr_contrast == 3
            tfr_c75_fooof = tfr;
            tfr_c75_fooof.powspctrm = fooofedtrl;
        elseif tfr_contrast == 4
            tfr_c100_fooof = tfr;
            tfr_c100_fooof.powspctrm = fooofedtrl;
        end
    end
    disp(upper('FOOOF done...'))

    %% Baseline
    % Raw powspctrm baselined
    cfg              = [];
    cfg.baseline     = [-1.5 -.25];
    cfg.baselinetype = 'db';
    tfr_c25_bl                       = ft_freqbaseline(cfg, tfr_c25);
    tfr_c50_bl                       = ft_freqbaseline(cfg, tfr_c50);
    tfr_c75_bl                       = ft_freqbaseline(cfg, tfr_c75);
    tfr_c100_bl                      = ft_freqbaseline(cfg, tfr_c100);

    % FOOOFed powspctrm baselined
    cfg              = [];
    cfg.baseline     = [-1.5 -.25];
    cfg.baselinetype = 'absolute';   % FOOOF already sets log scale, so no 'dB' here
    tfr_c25_fooof_bl                 = ft_freqbaseline(cfg, tfr_c25_fooof);
    tfr_c50_fooof_bl                 = ft_freqbaseline(cfg, tfr_c50_fooof);
    tfr_c75_fooof_bl                 = ft_freqbaseline(cfg, tfr_c75_fooof);
    tfr_c100_fooof_bl                = ft_freqbaseline(cfg, tfr_c100_fooof);
    disp(upper('Baseline done...'))

    %% Smooth powerspectra
    orig_freq = 30:5:120; % 19 values
    new_freq  = 30:1:120; % 91 values
    tfr_c25_fooof_bl_smooth  = smooth_tfr(tfr_c25_fooof_bl, orig_freq, new_freq);
    tfr_c50_fooof_bl_smooth  = smooth_tfr(tfr_c50_fooof_bl, orig_freq, new_freq);
    tfr_c75_fooof_bl_smooth  = smooth_tfr(tfr_c75_fooof_bl, orig_freq, new_freq);
    tfr_c100_fooof_bl_smooth = smooth_tfr(tfr_c100_fooof_bl, orig_freq, new_freq);
    disp(upper('Smoothing done...'))

    %% Save data
    cd(datapath)
    save data_tfr tfr_c25 tfr_c50 tfr_c75 tfr_c100 ...
        tfr_c25_fooof tfr_c50_fooof tfr_c75_fooof tfr_c100_fooof ...
        tfr_c25_bl tfr_c50_bl tfr_c75_bl tfr_c100_bl ...
        tfr_c25_fooof_bl tfr_c50_fooof_bl tfr_c75_fooof_bl tfr_c100_fooof_bl ...
        tfr_c25_fooof_bl_smooth tfr_c50_fooof_bl_smooth tfr_c75_fooof_bl_smooth tfr_c100_fooof_bl_smooth
    clc
    fprintf('Subject GCP %s (%.3d/%.3d) TFR DATA computed... \n', num2str(subjects{subj}), subj, length(subjects))
    %end
end

%% Convert TFR data to POWSCPTRM (channels x frequency)
clc
disp('POWER ANALYSIS')
baseline_period = [-1.5 -0.25];
analysis_period = [0.3 2]; % only start from 300ms after stimulus presentation
freq_range = [30 90];
[subjects, path] = setup('GCP');

for subj = 1 : length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    %if ~isfile(strcat([datapath, '/power_spectra.mat'])) % only new data
    % Load data
    cd(datapath);
    load('data_tfr.mat');

    %% Select analysis and baseline period data
    % (1) Analysis period data, no baseline
    % (2) Analysis period data, baselined
    % (3) Analysis period data, FOOOFed, baselined, and smoothed
    % (4) Baseline period data (to compare with (1) non-baselined data for percentage change)

    %% Select data
    pow_c25                                       = select_data(analysis_period, freq_range, tfr_c25);
    pow_c25_baselined                             = select_data(analysis_period, freq_range, tfr_c25_bl);
    pow_c25_fooof_bl_smooth                       = select_data(analysis_period, freq_range, tfr_c25_fooof_bl_smooth);
    pow_c25_baseline_period                       = select_data(baseline_period, freq_range, tfr_c25);

    pow_c50                                       = select_data(analysis_period, freq_range, tfr_c50);
    pow_c50_baselined                             = select_data(analysis_period, freq_range, tfr_c50_bl);
    pow_c50_fooof_bl_smooth                       = select_data(analysis_period, freq_range, tfr_c50_fooof_bl_smooth);
    pow_c50_baseline_period                       = select_data(baseline_period, freq_range, tfr_c50);

    pow_c75                                       = select_data(analysis_period, freq_range, tfr_c75);
    pow_c75_baselined                             = select_data(analysis_period, freq_range, tfr_c75_bl);
    pow_c75_fooof_bl_smooth                       = select_data(analysis_period, freq_range, tfr_c75_fooof_bl_smooth);
    pow_c75_baseline_period                       = select_data(baseline_period, freq_range, tfr_c75);

    pow_c100                                       = select_data(analysis_period, freq_range, tfr_c100);
    pow_c100_baselined                             = select_data(analysis_period, freq_range, tfr_c100_bl);
    pow_c100_fooof_bl_smooth                       = select_data(analysis_period, freq_range, tfr_c100_fooof_bl_smooth);
    pow_c100_baseline_period                       = select_data(baseline_period, freq_range, tfr_c100);

    %% Remove time dimension for POWSCPTRM (channels x frequency)
    pow_c25                                       = remove_time_dimension(pow_c25);
    pow_c25_baselined                             = remove_time_dimension(pow_c25_baselined);
    pow_c25_fooof_bl_smooth                       = remove_time_dimension(pow_c25_fooof_bl_smooth);
    pow_c25_baseline_period                       = remove_time_dimension(pow_c25_baseline_period);

    pow_c50                                       = remove_time_dimension(pow_c50);
    pow_c50_baselined                             = remove_time_dimension(pow_c50_baselined);
    pow_c50_fooof_bl_smooth                       = remove_time_dimension(pow_c50_fooof_bl_smooth);
    pow_c50_baseline_period                       = remove_time_dimension(pow_c50_baseline_period);

    pow_c75                                       = remove_time_dimension(pow_c75);
    pow_c75_baselined                             = remove_time_dimension(pow_c75_baselined);
    pow_c75_fooof_bl_smooth                       = remove_time_dimension(pow_c75_fooof_bl_smooth);
    pow_c75_baseline_period                       = remove_time_dimension(pow_c75_baseline_period);

    pow_c100                                       = remove_time_dimension(pow_c100);
    pow_c100_baselined                             = remove_time_dimension(pow_c100_baselined);
    pow_c100_fooof_bl_smooth                       = remove_time_dimension(pow_c100_fooof_bl_smooth);
    pow_c100_baseline_period                       = remove_time_dimension(pow_c100_baseline_period);

    %% Save data
    savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save power_spectra pow_c25 pow_c25_baselined pow_c25_fooof_bl_smooth pow_c25_baseline_period ...
        pow_c50 pow_c50_baselined pow_c50_fooof_bl_smooth pow_c50_baseline_period ...
        pow_c75 pow_c75_baselined pow_c75_fooof_bl_smooth pow_c75_baseline_period ...
        pow_c100 pow_c100_baselined pow_c100_fooof_bl_smooth pow_c100_baseline_period
end
%end

%% Define channels
datapath = strcat(path, subjects{1}, '/eeg');
cd(datapath);
% Occipital channels
occ_channels = {};
pow_label = pow_c25;
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
    % Load power spectra data
    load(strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra'))

    % Find channels and frequencies of interest
    channels_idx = ismember(pow_c25_fooof_bl_smooth.label, channels);
    freq_idx = find(pow_c25_fooof_bl_smooth.freq >= 30 & pow_c25_fooof_bl_smooth.freq <= 90);

    % Extract gamma peak power and freq for analysis time window
    % -------- 25% contrast --------
    c25_gamma_power = mean(pow_c25_fooof_bl_smooth.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c25_gamma_power, pow_c25_fooof_bl_smooth.freq(freq_idx));

    if isempty(peaks)
        c25_pow  = NaN;
        c25_freq = NaN;
    else
        [c25_pow, peak_idx] = max(peaks);
        c25_freq = locs(peak_idx);

        % ±5 Hz range
        c25_freq_range = c25_freq + [-5 5];
        c25_freq_idx_range = find(pow_c25_fooof_bl_smooth.freq >= c25_freq_range(1) & ...
            pow_c25_fooof_bl_smooth.freq <= c25_freq_range(2));

        c25_pow = mean(c25_gamma_power(c25_freq_idx_range));
    end


    % -------- 50% contrast --------
    c50_gamma_power = mean(pow_c50_fooof_bl_smooth.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c50_gamma_power, pow_c50_fooof_bl_smooth.freq(freq_idx));

    if isempty(peaks)
        c50_pow  = NaN;
        c50_freq = NaN;
    else
        [c50_pow, peak_idx] = max(peaks);
        c50_freq = locs(peak_idx);

        % ±5 Hz range
        c50_freq_range = c50_freq + [-5 5];
        c50_freq_idx_range = find(pow_c50_fooof_bl_smooth.freq >= c50_freq_range(1) & ...
            pow_c50_fooof_bl_smooth.freq <= c50_freq_range(2));

        c50_pow = mean(c50_gamma_power(c50_freq_idx_range));
    end


    % -------- 75% contrast --------
    c75_gamma_power = mean(pow_c75_fooof_bl_smooth.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c75_gamma_power, pow_c75_fooof_bl_smooth.freq(freq_idx));

    if isempty(peaks)
        c75_pow  = NaN;
        c75_freq = NaN;
    else
        [c75_pow, peak_idx] = max(peaks);
        c75_freq = locs(peak_idx);

        % ±5 Hz range
        c75_freq_range = c75_freq + [-5 5];
        c75_freq_idx_range = find(pow_c75_fooof_bl_smooth.freq >= c75_freq_range(1) & ...
            pow_c75_fooof_bl_smooth.freq <= c75_freq_range(2));

        c75_pow = mean(c75_gamma_power(c75_freq_idx_range));
    end


    % -------- 100% contrast --------
    c100_gamma_power = mean(pow_c100_fooof_bl_smooth.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c100_gamma_power, pow_c100_fooof_bl_smooth.freq(freq_idx));

    if isempty(peaks)
        c100_pow  = NaN;
        c100_freq = NaN;
    else
        [c100_pow, peak_idx] = max(peaks);
        c100_freq = locs(peak_idx);

        % ±5 Hz range
        c100_freq_range = c100_freq + [-5 5];
        c100_freq_idx_range = find(pow_c100_fooof_bl_smooth.freq >= c100_freq_range(1) & ...
            pow_c100_fooof_bl_smooth.freq <= c100_freq_range(2));

        c100_pow = mean(c100_gamma_power(c100_freq_idx_range));
    end

    % Control figure
    colors = color_def('GCP');
    close all
    figure('Position',  [0 0 1512 982], 'Color', 'W');
    plot(pow_c25_fooof_bl_smooth.freq(freq_idx), c25_gamma_power, 'Color', colors(1, :), 'LineWidth', 5)
    xlabel('Frequency')
    ylabel('Power [dB]')
    hold on
    plot(pow_c50_fooof_bl_smooth.freq(freq_idx), c50_gamma_power, 'Color', colors(2, :), 'LineWidth', 5)
    plot(pow_c75_fooof_bl_smooth.freq(freq_idx), c75_gamma_power, 'Color', colors(3, :), 'LineWidth', 5)
    plot(pow_c100_fooof_bl_smooth.freq(freq_idx), c100_gamma_power, 'Color', colors(4, :), 'LineWidth', 5)
    legend({'25%', '50%', '75%', '100%'})
    set(gca, 'FontSize', 25)
    xline(c25_freq, 'LineStyle', '--', 'LineWidth', 5, 'Color', colors(1, :))
    xline(c50_freq, 'LineStyle', '--', 'LineWidth', 5, 'Color', colors(2, :))
    xline(c75_freq, 'LineStyle', '--',  'LineWidth', 5, 'Color', colors(3, :))
    xline(c100_freq, 'LineStyle', '--',  'LineWidth', 5, 'Color', colors(4, :))
    yline(0, '--');
    title(['Subj ' num2str(subjects{subj}) ': FOOOFed baselined powerspectrum'])
    maxylim = max([abs(c25_gamma_power), abs(c50_gamma_power), abs(c75_gamma_power), abs(c100_gamma_power)]);
    ylim([-maxylim*1.2 maxylim*1.2])
    text(30, -maxylim*0.5, [num2str(c25_freq) 'Hz'], 'FontSize', 25, 'Color', colors(1, :))
    text(38, -maxylim*0.5, [num2str(c50_freq) 'Hz'], 'FontSize', 25, 'Color', colors(2, :))
    text(46, -maxylim*0.5, [num2str(c75_freq) 'Hz'], 'FontSize', 25, 'Color', colors(3, :))
    text(54, -maxylim*0.5, [num2str(c100_freq) 'Hz'], 'FontSize', 25, 'Color', colors(4, :))
    text(30, -maxylim*0.65, [num2str(round(c25_pow, 4)) ], 'FontSize', 25, 'Color', colors(1, :))
    text(38, -maxylim*0.65, [num2str(round(c50_pow, 4)) ], 'FontSize', 25, 'Color', colors(2, :))
    text(46, -maxylim*0.65, [num2str(round(c75_pow, 4)) ], 'FontSize', 25, 'Color', colors(3, :))
    text(54, -maxylim*0.65, [num2str(round(c100_pow, 4))], 'FontSize', 25, 'Color', colors(4, :))
    saveas(gcf, ['/Users/Arne/Documents/GitHub/GCP/controls/powspctrm/GCP_controls_powscptrm_subj' num2str(subjects{subj}), '.png'])

    % Create across condition structure
    subject_id = repmat(str2num(subjects{subj}), 4, 1);
    conditions = [1; 2; 3; 4];
    powers = [c25_pow; c50_pow; c75_pow; c100_pow];
    frequencies = [c25_freq; c50_freq; c75_freq; c100_freq];

    subj_data_eeg = struct('ID', num2cell(subject_id), ...
        'Condition',     num2cell(conditions), ...
        'Power',         num2cell(powers), ...
        'Frequency',     num2cell(frequencies));

    % Save data
    savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save eeg_matrix_subj subj_data_eeg
    save pow c25_pow c50_pow c75_pow c100_pow
    save freq c25_freq c50_freq c75_freq c100_freq

    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' gamma peak POWER and FREQUENCY extracted.'])

    % Append to the final structure array
    eeg_data = [eeg_data; subj_data_eeg];
end
save /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/eeg_matrix eeg_data
clc
disp('EEG Feature Matrix created')
