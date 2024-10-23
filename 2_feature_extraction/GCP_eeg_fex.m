%% GCP EEG Feature Extraction
%
% Extracted features:
%   Power Spectrum
%   TFR    

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.2');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Extract TFR HIGH CONTRAST 
% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    close all
    load dataEEG_TFRhc
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to orientation conditions
    ind1=find(dataTFRhc.trialinfo==0);   % HC 0°
    ind2=find(dataTFRhc.trialinfo==45);  % HC 45°
    ind3=find(dataTFRhc.trialinfo==90);  % HC 90°
    ind4=find(dataTFRhc.trialinfo==115); % HC 115°
    indall = find(ismember(dataTFRhc.trialinfo, [0 45 90 115])); % HC all orientations

    %% Time frequency analysis
    cfg             = [];
    cfg.output      = 'pow';
    cfg.method      = 'mtmconvol';
    cfg.taper       = 'dpss';
    cfg.foi         = 30:4:120;                          % Analysis 30 to 90 Hz in steps of 1 Hz
    cfg.tapsmofrq   = 11;                               % 10 tapers
    cfg.t_ftimwin   = ones(length(cfg.foi),1).*0.5;     % Length of time window = 0.5 sec
    cfg.toi         = -1:0.05:3;                       % What time window we are interested in, rsp. to stim presentation
    cfg.keeptrials  = 'no';
    cfg.trials = ind1;
    tfr_hc_0 = ft_freqanalysis(cfg,dataTFRhc);
    cfg.trials = ind2;
    tfr_hc_45 = ft_freqanalysis(cfg,dataTFRhc);
    cfg.trials = ind3;
    tfr_hc_90 = ft_freqanalysis(cfg,dataTFRhc);
    cfg.trials = ind4;
    tfr_hc_115 = ft_freqanalysis(cfg,dataTFRhc);
    cfg.trials = indall;
    tfr_hc_all = ft_freqanalysis(cfg,dataTFRhc);

    %% Baseline
    cfg             = [];
    cfg.baseline    = [-Inf -.25];
    cfg.baselinetype = 'db';
    tfr_hc_0 = ft_freqbaseline(cfg,tfr_hc_0);
    tfr_hc_45 = ft_freqbaseline(cfg,tfr_hc_45);
    tfr_hc_90 = ft_freqbaseline(cfg,tfr_hc_90);
    tfr_hc_115 = ft_freqbaseline(cfg,tfr_hc_115);
    tfr_hc_all = ft_freqbaseline(cfg,tfr_hc_all);

    %% Save data
    cd(datapath)
    save tfr_hc tfr_hc_0 tfr_hc_45 tfr_hc_90 tfr_hc_115 tfr_hc_all
    fprintf('Subject %.2s / %.3d TFR HC computed \n', num2str(subj), length(subjects))
end

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.2');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Extract TFR LOW CONTRAST
% Read data, segment and convert to FieldTrip data structure
for subj = 6:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    close all
    load dataEEG_TFRlc
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to orientation conditions
    ind1=find(dataTFRlc.trialinfo==1000);   % LC 0°
    ind2=find(dataTFRlc.trialinfo==1045);   % LC 45°
    ind3=find(dataTFRlc.trialinfo==1090);   % LC 90°
    ind4=find(dataTFRlc.trialinfo==1115);   % LC 115°
    indall = find(ismember(dataTFRlc.trialinfo, [1000 1045 1090 1115])); % LC all orientations

    %% Time frequency analysis
    cfg             = [];
    cfg.output      = 'pow';
    cfg.method      = 'mtmconvol';
    cfg.taper       = 'dpss';
    cfg.foi         = 30:4:120;                          % Analysis 30 to 90 Hz in steps of 1 Hz
    cfg.tapsmofrq   = 11;                               % 10 tapers
    cfg.t_ftimwin   = ones(length(cfg.foi),1).*0.5;     % Length of time window = 0.5 sec
    cfg.toi         = -1:0.05:3;                        % What time window we are interested in, rsp. to stim presentation
    cfg.keeptrials  = 'no';
    cfg.trials = ind1;
    tfr_lc_0 = ft_freqanalysis(cfg,dataTFRlc);
    cfg.trials = ind2;
    tfr_lc_45 = ft_freqanalysis(cfg,dataTFRlc);
    cfg.trials = ind3;
    tfr_lc_90 = ft_freqanalysis(cfg,dataTFRlc);
    cfg.trials = ind4;
    tfr_lc_115 = ft_freqanalysis(cfg,dataTFRlc);
    cfg.trials = indall;
    tfr_lc_all = ft_freqanalysis(cfg,dataTFRlc);

    %% Baseline
    cfg             = [];
    cfg.baseline    = [-Inf -.25];
    cfg.baselinetype = 'db';
    tfr_lc_0 = ft_freqbaseline(cfg,tfr_lc_0);
    tfr_lc_45 = ft_freqbaseline(cfg,tfr_lc_45);
    tfr_lc_90 = ft_freqbaseline(cfg,tfr_lc_90);
    tfr_lc_115 = ft_freqbaseline(cfg,tfr_lc_115);
    tfr_lc_all = ft_freqbaseline(cfg,tfr_lc_all);

    %% Save data
    cd(datapath)
    save tfr_lc tfr_lc_0 tfr_lc_45 tfr_lc_90 tfr_lc_115 tfr_lc_all
    fprintf('Subject %.2s / %.3d TFR LC computed \n', num2str(subj), length(subjects))
end