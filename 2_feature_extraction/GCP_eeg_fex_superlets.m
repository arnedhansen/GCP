%% GCP EEG Superlets
%
%  Computes the adaptive superresolution wavelet (superlet) transform on
%  input data to produce a time-frequency representation. For each
%  frequency of interest, the closest integer order from the order
%  interval will be chosen to produce each superlet. A superlet is a set
%  of wavelets with the same center frequency but different number of
%  cycles.
%
%  Reference: Time-frequency super-resolution with superlets
%  https://www.nature.com/articles/s41467-020-20539-9

%% Setup
startup
[subjects, path] = setup('GCP');

%% Define channels
datapath = strcat(path, subjects{1}, '/eeg');
cd(datapath);
load('power_spectra.mat')
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

%% Extract TFR WITH SUPERLETS
% Read data, segment and convert to FieldTrip data structure
for subj = 1 : length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath)
    close all
    fprintf('Loading data for Subject GCP %s (%.2d/%.2d)... \n', num2str(subjects{subj}), subj, length(subjects))
    load dataEEG
    load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions
    ind61 = find(dataEEG_c25.trialinfo  == 61);
    ind62 = find(dataEEG_c50.trialinfo  == 62);
    ind63 = find(dataEEG_c75.trialinfo  == 63);
    ind64 = find(dataEEG_c100.trialinfo == 64);

    %% Define index of channels for superlet TFR
    chn = find(ismember(dataEEG_c25.label, occ_channels));

    %% SUPERLET Time Frequency Analysis
    % Parameters for superlets
    Fs     = dataEEG_c25.fsample;         % sampling rate
    F      = 30:5:120;                    % frequencies of interest
    Ncyc   = 3;                           % base number of cycles
    ord    = [1 30];                      % superresolution orders 1→7
    mult   = 1;                           % multiplicative superlets

    % Preallocate output arrays: [nChan × nFreq × nTime]
    nCh = numel(chn);
    nF  = numel(F);
    nT  = numel(dataEEG_c25.time{1});
    suplet25_pow  = zeros(nCh, nF, nT);
    suplet50_pow  = zeros(nCh, nF, nT);
    suplet75_pow  = zeros(nCh, nF, nT);
    suplet100_pow = zeros(nCh, nF, nT);

    % Loop over channels
    for c = 1:nCh
        % Build [buffers×samples] for each condition and channel c
        X25 = cell2mat( cellfun(@(trl) trl(chn(c),:), dataEEG_c25.trial(ind61),  'uni',0) );
        X25 = reshape(X25, length(ind61), []);
        X50 = cell2mat( cellfun(@(trl) trl(chn(c),:), dataEEG_c50.trial(ind62),  'uni',0) );
        X50 = reshape(X50, length(ind62), []);
        X75 = cell2mat( cellfun(@(trl) trl(chn(c),:), dataEEG_c75.trial(ind63),  'uni',0) );
        X75 = reshape(X75, length(ind63), []);
        X100= cell2mat( cellfun(@(trl) trl(chn(c),:), dataEEG_c100.trial(ind64),'uni',0) );
        X100= reshape(X100, length(ind64), []);

        % Compute Superlet TFR for this channel
        suplet25  = aslt(X25,  Fs, F, Ncyc, ord, mult); % [nFreq × nTime]
        suplet50  = aslt(X50,  Fs, F, Ncyc, ord, mult);
        suplet75  = aslt(X75,  Fs, F, Ncyc, ord, mult);
        suplet100 = aslt(X100, Fs, F, Ncyc, ord, mult);

        % Store into 3D arrays
        suplet25_pow(c,:,:)  = suplet25;
        suplet50_pow(c,:,:)  = suplet50;
        suplet75_pow(c,:,:)  = suplet75;
        suplet100_pow(c,:,:) = suplet100;
    end

    % Convert back to FieldTrip structure
    % 25% contrast concentric dynamic inward
    tfr_c25            = [];
    tfr_c25.powspctrm  = suplet25_pow; % [nChan × nFreq × nTime]
    tfr_c25.freq       = F;
    tfr_c25.time       = dataEEG_c25.time{1};
    tfr_c25.label      = dataEEG_c25.label(chn);
    tfr_c25.dimord     = 'chan_freq_time';

    % 50% contrast concentric dynamic inward
    tfr_c50            = [];
    tfr_c50.powspctrm  = suplet50_pow;
    tfr_c50.freq       = F;
    tfr_c50.time       = dataEEG_c50.time{1};
    tfr_c50.label      = dataEEG_c50.label(chn);
    tfr_c50.dimord     = 'chan_freq_time';

    % 75% contrast concentric dynamic inward
    tfr_c75            = [];
    tfr_c75.powspctrm  = suplet75_pow;
    tfr_c75.freq       = F;
    tfr_c75.time       = dataEEG_c75.time{1};
    tfr_c75.label      = dataEEG_c75.label(chn);
    tfr_c75.dimord     = 'chan_freq_time';

    % 100% contrast concentric dynamic inward
    tfr_c100           = [];
    tfr_c100.powspctrm = suplet100_pow;
    tfr_c100.freq      = F;
    tfr_c100.time      = dataEEG_c100.time{1};
    tfr_c100.label     = dataEEG_c100.label(chn);
    tfr_c100.dimord     = 'chan_freq_time';

    %% Baseline
    % Raw powspctrm baselined
    cfg              = [];
    cfg.baseline     = [-1.5 -.25];
    cfg.baselinetype = 'db';
    tfr_c25_bl       = ft_freqbaseline(cfg, tfr_c25);
    tfr_c50_bl       = ft_freqbaseline(cfg, tfr_c50);
    tfr_c75_bl       = ft_freqbaseline(cfg, tfr_c75);
    tfr_c100_bl      = ft_freqbaseline(cfg, tfr_c100);

    %% Save data
    cd(datapath)
    save data_tfr_superlets tfr_c25 tfr_c50 tfr_c75 tfr_c100 ...
        tfr_c25_bl tfr_c50_bl tfr_c75_bl tfr_c100_bl
    save superlets suplet25_pow suplet50_pow suplet75_pow suplet100_pow
    clc
    fprintf('Subject GCP %s (%.3d/%.3d) TFR DATA computed... \n', num2str(subjects{subj}), subj, length(subjects))
end

%% Convert TFR data to POWSCPTRM (channels x frequency)
startup
clc
disp('POWER ANALYSIS')
analysis_period = [0.3 2];
baseline_period = [-1.5 -0.25];
freq_range = [30 90];
[subjects, path] = setup('GCP');

for subj = 1 : length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    % Load data
    cd(datapath);
    load('data_tfr.mat');

    %% Select data
    pow_c25                                       = select_data(analysis_period, freq_range, tfr_c25);
    pow_c25_baselined                             = select_data(analysis_period, freq_range, tfr_c25_bl);

    pow_c50                                       = select_data(analysis_period, freq_range, tfr_c50);
    pow_c50_baselined                             = select_data(analysis_period, freq_range, tfr_c50_bl);

    pow_c75                                       = select_data(analysis_period, freq_range, tfr_c75);
    pow_c75_baselined                             = select_data(analysis_period, freq_range, tfr_c75_bl);

    pow_c100                                      = select_data(analysis_period, freq_range, tfr_c100);
    pow_c100_baselined                            = select_data(analysis_period, freq_range, tfr_c100_bl);

    %% Remove time dimension for POWSCPTRM (channels x frequency)
    pow_c25                                       = remove_time_dimension(pow_c25);
    pow_c25_baselined                             = remove_time_dimension(pow_c25_baselined);

    pow_c50                                       = remove_time_dimension(pow_c50);
    pow_c50_baselined                             = remove_time_dimension(pow_c50_baselined);

    pow_c75                                       = remove_time_dimension(pow_c75);
    pow_c75_baselined                             = remove_time_dimension(pow_c75_baselined);

    pow_c100                                      = remove_time_dimension(pow_c100);
    pow_c100_baselined                            = remove_time_dimension(pow_c100_baselined);

    %% Save data
    savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save power_spectra_superlets pow_c25 pow_c25_baselined    ...
        pow_c50 pow_c50_baselined    ...
        pow_c75 pow_c75_baselined    ...
        pow_c100 pow_c100_baselined
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
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' gamma peak POWER and FREQUENCY extracting...'])

    % Load power spectra data
    load(strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/power_spectra_superlets'))

    % Find channels and frequencies of interest
    channels_idx = ismember(pow_c25.label, channels);
    freq_idx = find(pow_c25.freq >= 30 & pow_c25.freq <= 90);

    % Find gamma peak for 25% contrast
    c25_gamma_power = mean(pow_c25.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c25_gamma_power, pow_c25.freq(freq_idx));
    [c25_pow, peak_idx] = max(peaks);
    c25_freq = locs(peak_idx);
    % Compute the peak gamma power as the mean power in the ±5 Hz range
    c25_freq_range = c25_freq + [-5, 5];  % ±5 Hz range
    c25_freq_idx_range = find(pow_c25_baselined.freq >= c25_freq_range(1) & pow_c25_baselined.freq <= c25_freq_range(2));
    c25_pow = mean(c25_gamma_power(c25_freq_idx_range));

    % Find gamma peak for 50% contrast
    c50_gamma_power = mean(pow_c50.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c50_gamma_power, pow_c50.freq(freq_idx));
    [c50_pow, peak_idx] = max(peaks);
    c50_freq = locs(peak_idx);
    % Compute the peak gamma power as the mean power in the ±5 Hz range
    c50_freq_range = c50_freq + [-5, 5];  % ±5 Hz range
    c50_freq_idx_range = find(pow_c50_baselined.freq >= c50_freq_range(1) & pow_c50_baselined.freq <= c50_freq_range(2));
    c50_pow = mean(c50_gamma_power(c50_freq_idx_range));

    % Find gamma peak for 75% contrast
    c75_gamma_power = mean(pow_c75.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c75_gamma_power, pow_c75.freq(freq_idx));
    if isempty(peaks)
        c75_pow = NaN;
        c75_freq = NaN;
    else
        [c75_pow, peak_idx] = max(peaks);
        c75_freq = locs(peak_idx);
        % Compute the peak gamma power as the mean power in the ±5 Hz range
        c75_freq_range = c75_freq + [-5, 5];  % ±5 Hz range
        c75_freq_idx_range = find(pow_c75_baselined.freq >= c75_freq_range(1) & pow_c75_baselined.freq <= c75_freq_range(2));
        c75_pow = mean(c75_gamma_power(c75_freq_idx_range));
    end

    % Find gamma peak for 100% contrast
    c100_gamma_power = mean(pow_c100.powspctrm(channels_idx, freq_idx), 1);
    [peaks, locs] = findpeaks(c100_gamma_power, pow_c100.freq(freq_idx));
    [c100_pow, peak_idx] = max(peaks);
    c100_freq = locs(peak_idx);
    % Compute the peak gamma power as the mean power in the ±5 Hz range
    c100_freq_range = c100_freq + [-5, 5];  % ±5 Hz range
    c100_freq_idx_range = find(pow_c100_baselined.freq >= c100_freq_range(1) & pow_c100_baselined.freq <= c100_freq_range(2));
    c100_pow = mean(c100_gamma_power(c100_freq_idx_range));

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
    savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/', subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save eeg_matrix_subj_superlets subj_data_eeg
    save pow_superlets c25_pow c50_pow c75_pow c100_pow
    save freq_superlets c25_freq c50_freq c75_freq c100_freq

    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' gamma peak POWER and FREQUENCY extracted.'])

    % Append to the final structure array
    eeg_data = [eeg_data; subj_data_eeg];
end
save /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/eeg_matrix_superlets eeg_data
disp('EEG Feature Matrix created')