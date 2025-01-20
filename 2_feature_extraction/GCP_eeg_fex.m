%% GCP EEG Feature Extraction
%
% Extracted features:
%   TFR
%   Power Spectrum
%   Gamma Peak Power
%   Gamma Peak Frequency

%% Setup
startup
[subjects, path] = setup('GCP');

%% Extract TFR HIGH CONTRAST
% Read data, segment and convert to FieldTrip data structure
for subj = 1 : length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    %if ~isfile(strcat([datapath, '/data_tfr.mat'])) % only new data
        cd(datapath)
        close all
        load dataEEG
        load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

        %% Identify indices of trials belonging to conditions
        ind61 = find(dataEEG_c25.trialinfo == 61);
        ind62 = find(dataEEG_c50.trialinfo == 62);
        ind63 = find(dataEEG_c75.trialinfo == 63);
        ind64 = find(dataEEG_c100.trialinfo == 64);

        %% Time frequency analysis
        cfg             = [];
        cfg.output      = 'pow';
        cfg.method      = 'mtmconvol';
        cfg.taper       = 'dpss';
        cfg.foi         = 30:4:120;                         % Analysis 30 to 90 Hz in steps of 4 Hz
        cfg.tapsmofrq   = 11;                               % 10 tapers
        cfg.t_ftimwin   = ones(length(cfg.foi),1).*0.5;     % Length of time window = 0.5 sec
        cfg.toi         = -1:0.05:3;
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

        %% Baseline
        cfg             = [];
        cfg.baseline    = [-.5 -.25];
        cfg.baselinetype = 'db';
        tfr_c25_bl                       = ft_freqbaseline(cfg, tfr_c25);
        tfr_c50_bl                       = ft_freqbaseline(cfg, tfr_c50);
        tfr_c75_bl                       = ft_freqbaseline(cfg, tfr_c75);
        tfr_c100_bl                      = ft_freqbaseline(cfg, tfr_c100);

        %% Save data
        cd(datapath)
        save data_tfr tfr_c25 tfr_c50 tfr_c75 tfr_c100 tfr_c25_bl tfr_c50_bl tfr_c75_bl tfr_c100_bl

        clc
        fprintf('Subject %.3d/%.3d TFR DATA computed \n', subj, length(subjects))
    %end
end

%% Load data and convert TFR data to POWSCPTRM (channels x frequency)
clc
disp('POWER ANALYSIS')
% Set analysis to 0-300ms or 300ms-2000ms after stimulus presentation
anal_period = 2;
baseline_period = [-0.5 -0.25];
if anal_period == 1
    analysis_period = [0 0.3]; % only starting stimulus activity
else
    analysis_period = [0.3 2]; % only start from 300ms after stimulus presentation
end
freq_range = [30 90];
[subjects, path] = setup('GCP');

for subj = 1 : length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    if ~isfile(strcat([datapath, '/power_spectra.mat'])) % only new data
        cd(datapath);

        % Load data
        load('data_tfr.mat');

        %% Select analysis and baseline period data
        % (1) Analysis period data, no baseline
        % (2) Analysis period data, baselined
        % (3) Baseline period data (to compare with (1) non-baselined data for percentage change)

        %% Power spectra calculations
        pow_c25                                       = select_data(analysis_period, freq_range, tfr_c25);
        pow_c25_baselined                             = select_data(analysis_period, freq_range, tfr_c25_bl);
        pow_c25_baseline_period                       = select_data(baseline_period, freq_range, tfr_c25);

        pow_c50                                       = select_data(analysis_period, freq_range, tfr_c50);
        pow_c50_baselined                             = select_data(analysis_period, freq_range, tfr_c50_bl);
        pow_c50_baseline_period                       = select_data(baseline_period, freq_range, tfr_c50);

        pow_c75                                       = select_data(analysis_period, freq_range, tfr_c75);
        pow_c75_baselined                             = select_data(analysis_period, freq_range, tfr_c75_bl);
        pow_c75_baseline_period                       = select_data(baseline_period, freq_range, tfr_c75);

        pow_c100                                       = select_data(analysis_period, freq_range, tfr_c100);
        pow_c100_baselined                             = select_data(analysis_period, freq_range, tfr_c100_bl);
        pow_c100_baseline_period                       = select_data(baseline_period, freq_range, tfr_c100);

        %% Remove time dimension for POWSCPTRM (channels x frequency)
        pow_c25                                       = remove_time_dimension(pow_c25);
        pow_c25_baselined                             = remove_time_dimension(pow_c25_baselined);
        pow_c25_baseline_period                       = remove_time_dimension(pow_c25_baseline_period);

        pow_c50                                       = remove_time_dimension(pow_c50);
        pow_c50_baselined                             = remove_time_dimension(pow_c50_baselined);
        pow_c50_baseline_period                       = remove_time_dimension(pow_c50_baseline_period);

        pow_c75                                       = remove_time_dimension(pow_c75);
        pow_c75_baselined                             = remove_time_dimension(pow_c75_baselined);
        pow_c75_baseline_period                       = remove_time_dimension(pow_c75_baseline_period);

        pow_c100                                       = remove_time_dimension(pow_c100);
        pow_c100_baselined                             = remove_time_dimension(pow_c100_baselined);
        pow_c100_baseline_period                       = remove_time_dimension(pow_c100_baseline_period);

        %% Save data
        savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/');
        mkdir(savepath)
        cd(savepath)
        if anal_period == 1
            save power_spectra_300 pow_c25 pow_c25_baselined pow_c25_baseline_period ...
                pow_c50 pow_c50_baselined pow_c50_baseline_period ...
                pow_c75 pow_c75_baselined pow_c75_baseline_period ...
                pow_c100 pow_c100_baselined pow_c100_baseline_period
        else
            save power_spectra pow_c25 pow_c25_baselined pow_c25_baseline_period ...
                pow_c50 pow_c50_baselined pow_c50_baseline_period ...
                pow_c75 pow_c75_baselined pow_c75_baseline_period ...
                pow_c100 pow_c100_baselined pow_c100_baseline_period
        end
    end
end

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
    if anal_period == 1
        load(strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra_300'))
    else
        load(strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra'))
    end

    % Find channels and frequencies of interest
    channels_idx = ismember(pow_c25_baselined.label, channels);
    freq_idx = find(pow_c25_baselined.freq >= 30 & pow_c100_baselined.freq <= 90);

    % Find gamma peak for 25% contrast
    c25_gamma_power = mean(pow_c25_baselined.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c25_gamma_power, pow_c25_baselined.freq(freq_idx));
    [c25_pow, peak_idx] = max(peaks);
    c25_freq = locs(peak_idx);

    % % Caclculate the range of frequencies for ±5 Hz around the peak frequency
    % c25_freq_range = c25_freq + [-5, 5];  % ±5 Hz range
    % c25_freq_idx_range = find(pow_c25_baselined.freq >= c25_freq_range(1) & pow_c25_baselined.freq <= c25_freq_range(2));
    %
    % % Compute the peak gamma power as the mean power in the ±5 Hz range
    % c25_pow = mean(c25_gamma_power(c25_freq_idx_range));

    % Find gamma peak for 50% contrast
    c50_gamma_power = mean(pow_c50_baselined.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c50_gamma_power, pow_c50_baselined.freq(freq_idx));
    [c50_pow, peak_idx] = max(peaks);
    c50_freq = locs(peak_idx);

    % Find gamma peak for 75% contrast
    c75_gamma_power = mean(pow_c75_baselined.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c75_gamma_power, pow_c75_baselined.freq(freq_idx));
    [c75_pow, peak_idx] = max(peaks);
    c75_freq = locs(peak_idx);

    % Find gamma peak for 100% contrast
    c100_gamma_power = mean(pow_c100_baselined.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c100_gamma_power, pow_c100_baselined.freq(freq_idx));
    [c100_pow, peak_idx] = max(peaks);
    c100_freq = locs(peak_idx);

    % Create across condition structure
    subject_id = repmat(str2num(subjects{subj}), 4, 1);
    conditions = [1; 2; 3; 4];
    powers = [c25_pow; c50_pow; c75_pow; c100_pow];
    frequencies = [c25_freq; c50_freq; c75_freq; c100_freq];

    subj_data_eeg = struct('ID', num2cell(subject_id), ...
        'Condition', num2cell(conditions), ...
        'Power', num2cell(powers), ...
        'Frequency', num2cell(frequencies));

    % Save data
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    if anal_period == 1
        save eeg_matrix_subj300 subj_data_eeg
        save pow300 c25_pow c50_pow c75_pow c100_pow
        save freq300 c25_freq c50_freq c75_freq c100_freq
    else
        save eeg_matrix_subj subj_data_eeg
        save pow c25_pow c50_pow c75_pow c100_pow
        save freq c25_freq c50_freq c75_freq c100_freq
    end

    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' gamma peak POWER and FREQUENCY extracted.'])

    % Append to the final structure array
    eeg_data = [eeg_data; subj_data_eeg];
end
save /Volumes/methlab/Students/Arne/GCP/data/features/eeg_matrix eeg_data
disp('EEG Feature Matrix created')