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
    if ~isfile(strcat([datapath, '/data_tfr_LCHC.mat'])) % only new data
        cd(datapath)
        close all
        load dataEEG_LCHC
        load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

        %% Identify indices of trials belonging to conditions
        %ind_LC = find(dataEEG_lc.trialinfo == [61, 62]);
        %ind_HC = find(dataEEG_hc.trialinfo == [63, 64]);

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

        % LOW contrast concentric dynamic inward
        %cfg.trials = ind_LC;
        tfr_lc = ft_freqanalysis(cfg,dataEEG_lc);

        % HIGH contrast concentric dynamic inward
        %cfg.trials = ind_HC;
        tfr_hc = ft_freqanalysis(cfg,dataEEG_hc);

        %% Baseline
        cfg             = [];
        cfg.baseline    = [-.5 -.25];
        cfg.baselinetype = 'db';
        tfr_lc_bl                       = ft_freqbaseline(cfg, tfr_lc);
        tfr_hc_bl                       = ft_freqbaseline(cfg, tfr_hc);

        %% Save data
        cd(datapath)
        save data_tfr_LCHC tfr_lc  tfr_lc_bl tfr_hc tfr_hc_bl

        clc
        fprintf('Subject %.3d/%.3d TFR DATA computed \n', subj, length(subjects))
    end
end

%% Load data and convert TFR data to POWSCPTRM (channels x frequency)
clc
disp('POWER ANALYSIS')
% Set analysis to 0-300ms or 300ms-2000ms after stimulus presentation
anal_period = 0;
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
    %if ~isfile(strcat([datapath, '/power_spectra.mat'])) % only new data
    cd(datapath);

    % Load data
    load('data_tfr_LCHC.mat');

    %% Select analysis and baseline period data
    % (1) Analysis period data, no baseline
    % (2) Analysis period data, baselined
    % (3) Baseline period data (to compare with (1) non-baselined data for percentage change)

    pow_lc                                       = select_data(analysis_period, freq_range, tfr_lc);
    pow_lc_baselined                             = select_data(analysis_period, freq_range, tfr_lc_bl);
    pow_lc_baseline_period                       = select_data(baseline_period, freq_range, tfr_lc);

    pow_hc                                       = select_data(analysis_period, freq_range, tfr_hc);
    pow_hc_baselined                             = select_data(analysis_period, freq_range, tfr_hc_bl);
    pow_hc_baseline_period                       = select_data(baseline_period, freq_range, tfr_hc);

    %% Remove time dimension for POWSCPTRM (channels x frequency)
    pow_lc                                       = remove_time_dimension(pow_lc);
    pow_lc_baselined                             = remove_time_dimension(pow_lc_baselined);
    pow_lc_baseline_period                       = remove_time_dimension(pow_lc_baseline_period);

    pow_hc                                       = remove_time_dimension(pow_hc);
    pow_hc_baselined                             = remove_time_dimension(pow_hc_baselined);
    pow_hc_baseline_period                       = remove_time_dimension(pow_hc_baseline_period);

    % Save data
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    if anal_period == 1
            save power_spectra_300 pow_lc pow_lc_baselined pow_lc_baseline_period pow_hc pow_hc_baselined pow_hc_baseline_period
    else
    save power_spectra_LCHC pow_lc pow_lc_baselined pow_lc_baseline_period pow_hc pow_hc_baselined pow_hc_baseline_period
    end
    %end
end

%% Define channels
datapath = strcat(path, subjects{1}, '/eeg');
cd(datapath);
% Occipital channels
occ_channels = {};
pow_label = pow_lc;
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
        load(strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra_LCHC'))
    end

    % Find channels and frequencies of interest
    channels_idx = ismember(pow_lc_baselined.label, channels);
    freq_idx = find(pow_lc_baselined.freq >= 30 & pow_hc_baselined.freq <= 90);

    % Find gamma peak for LOW contrast
    lc_gamma_power = mean(pow_lc_baselined.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(lc_gamma_power, pow_lc_baselined.freq(freq_idx));
    [lc_pow, peak_idx] = max(peaks);
    lc_freq = locs(peak_idx);

    % % Calculate the range of frequencies for ±5 Hz around the peak frequency
    % lc_freq_range = lc_freq + [-5, 5];  % ±5 Hz range
    % lc_freq_idx_range = find(pow_lc_baselined.freq >= lc_freq_range(1) & pow_lc_baselined.freq <= lc_freq_range(2));
    % 
    % % Compute the peak gamma power as the mean power in the ±5 Hz range
    % lc_pow = mean(lc_gamma_power(lc_freq_idx_range));

    % Find gamma peak for HIGH contrast
    hc_gamma_power = mean(pow_hc_baselined.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(hc_gamma_power, pow_hc_baselined.freq(freq_idx));
    [hc_pow, peak_idx] = max(peaks);
    hc_freq = locs(peak_idx);

    % % Calculate the range of frequencies for ±5 Hz around the peak frequency
    % hc_freq_range = hc_freq + [-5, 5];  % ±5 Hz range
    % hc_freq_idx_range = find(pow_hc_baselined.freq >= hc_freq_range(1) & pow_hc_baselined.freq <= hc_freq_range(2));
    % 
    % % Compute the peak gamma power as the mean power in the ±5 Hz range
    % hc_pow = mean(hc_gamma_power(hc_freq_idx_range));
    if subj == 9
        hc_pow = 0.1149;
        hc_freq = 86;
    end

    % Create across condition structure
    subject_id = [str2num(subjects{subj}); str2num(subjects{subj})];
    subj_data_eeg = struct('ID', num2cell(subject_id(1:2)), 'Condition', num2cell([1; 2]), 'Power', num2cell([lc_pow; hc_pow]), 'Frequency', num2cell([lc_freq; hc_freq]));

    % Save data
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    if anal_period == 1
        save eeg_matrix_subj300 subj_data_eeg
        save pow300 lc_pow hc_pow
        save freq300 lc_freq hc_freq
    else
        save eeg_matrix_subj_LCHC subj_data_eeg
        save pow_LCHC lc_pow hc_pow
        save freq_LCHC lc_freq hc_freq
    end

    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' gamma peak POWER and FREQUENCY extracted.'])

    % Append to the final structure array
    eeg_data = [eeg_data; subj_data_eeg];
end
save /Volumes/methlab/Students/Arne/GCP/data/features/eeg_matrix_LCHC eeg_data
disp('EEG Feature Matrix created')