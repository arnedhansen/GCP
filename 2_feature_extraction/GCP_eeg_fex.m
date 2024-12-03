%% GCP EEG Feature Extraction
%
% Extracted features:
%   TFR

%% Setup
[subjects, path] = setup('GCP');

%% Extract TFR HIGH CONTRAST
% Read data, segment and convert to FieldTrip data structure
for subj = 1 : length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    close all
    load dataEEG
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions 
    ind61 = find(dataEEG_lc.trialinfo == 61);
    ind62 = find(dataEEG_hc.trialinfo == 62);

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
    cfg.trials = ind61;
    tfr_lc = ft_freqanalysis(cfg,dataEEG_lc);

    % HIGH contrast concentric dynamic inward
    cfg.trials = ind62;
    tfr_hc = ft_freqanalysis(cfg,dataEEG_hc);

    %% Baseline
    cfg             = [];
    cfg.baseline    = [-.5 -.25];
    cfg.baselinetype = 'db';
    tfr_lc_bl                       = ft_freqbaseline(cfg, tfr_lc);
    tfr_hc_bl                       = ft_freqbaseline(cfg, tfr_hc);

    %% Save data
    cd(datapath)
    save data_tfr tfr_lc  tfr_lc_bl tfr_hc tfr_hc_bl

    fprintf('Subject %.3d/%.3d TFR DATA computed \n', subj, length(subjects))
end