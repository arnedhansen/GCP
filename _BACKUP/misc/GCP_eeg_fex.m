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

%% Extract POWER HIGH CONTRAST
% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    close all
    load dataEEGhc
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to orientation conditions
    ind1=find(dataHC.trialinfo==0);   % HC 0°
    ind2=find(dataHC.trialinfo==45);  % HC 45°
    ind3=find(dataHC.trialinfo==90);  % HC 90°
    ind4=find(dataHC.trialinfo==115); % HC 115°

    %% Frequency analysis
    cfg = [];
    cfg.latency = [0 2];
    dat = ft_selectdata(cfg,dataHC);
    cfg = [];% empty config
    cfg.output = 'pow';% estimates power only
    cfg.method = 'mtmfft';% multi taper fft method
    cfg.taper = 'dpss';% multiple tapers
    cfg.tapsmofrq = 4;% smoothening frequency around foi
    cfg.foilim = [4 120];% frequencies of interest (foi)
    cfg.keeptrials = 'no';
    cfg.pad = 10;
    cfg.trials = ind1;
    powload_hc_0 = ft_freqanalysis(cfg,dat);
    cfg.trials = ind2;
    powload_hc_45 = ft_freqanalysis(cfg,dat);
    cfg.trials = ind3;
    powload_hc_90 = ft_freqanalysis(cfg,dat);
    cfg.trials = ind4;
    powload_hc_115 = ft_freqanalysis(cfg,dat);

    % Add baseline
    cfg             = [];
    cfg.baseline    = [-Inf -.25];
    cfg.baselinetype = 'db';
    tfr_hc = ft_freqbaseline(cfg, tfr_hc_mul);
    tfr_lc = ft_freqbaseline(cfg, tfr_lc_mul);

    %% Save data
    cd(datapath)
    save power_hc powload_hc_0 powload_hc_45 powload_hc_90 powload_hc_115
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

%% Extract POWER LOW CONTRAST
% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    close all
    load dataEEGlc
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions
    ind5=find(dataLC.trialinfo==1000);
    ind6=find(dataLC.trialinfo==1045);
    ind7=find(dataLC.trialinfo==1090);
    ind8=find(dataLC.trialinfo==1115);

    %% Frequency analysis
    cfg = [];
    cfg.latency = [0 2];
    dat = ft_selectdata(cfg,dataLC);
    cfg = [];% empty config
    cfg.output = 'pow';% estimates power only
    cfg.method = 'mtmfft';% multi taper fft method
    cfg.taper = 'dpss';% multiple tapers
    cfg.tapsmofrq = 1;% smoothening frequency around foi
    cfg.foilim = [4 120];% frequencies of interest (foi)
    cfg.keeptrials = 'no';
    cfg.pad = 10;
    cfg.trials = ind5;
    powload_lc_0 = ft_freqanalysis(cfg,dat);
    cfg.trials = ind6;
    powload_lc_45 = ft_freqanalysis(cfg,dat);
    cfg.trials = ind7;
    powload_lc_90 = ft_freqanalysis(cfg,dat);
    cfg.trials = ind8;
    powload_lc_115 = ft_freqanalysis(cfg,dat);

    %% Save data
    cd(datapath)
    save power_lc powload_lc_0 powload_lc_45 powload_lc_90 powload_lc_115
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
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 4:1:120;                         % analysis 4 to 120 Hz in steps of 1 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -1:0.05:3;
    cfg.keeptrials = 'no';
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
    tfr_hc_0 = ft_freqanalysis(cfg,tfr_hc_0);
    tfr_hc_45 = ft_freqanalysis(cfg,tfr_hc_45);
    tfr_hc_90 = ft_freqanalysis(cfg,tfr_hc_90);
    tfr_hc_115 = ft_freqanalysis(cfg,tfr_hc_115);
    tfr_hc_all = ft_freqanalysis(cfg,tfr_hc_all);

    %% Save data
    cd(datapath)
    save tfr_hc tfr_hc_0 tfr_hc_45 tfr_hc_90 tfr_hc_115 tfr_hc_all
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
for subj = 1:length(subjects)
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
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 4:1:120;                         % analysis 4 to 120 Hz in steps of 1 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -1:0.05:3;
    cfg.keeptrials = 'no';

    cfg.trials = ind1;
    tfr_lc_0 = ft_freqanalysis(cfg,dataTFRlc);
    cfg.trials = ind2;
    tfr_lc_45 = ft_freqanalysis(cfg,dataTFRlc);v
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
    tfr_hc_0 = ft_freqanalysis(cfg,tfr_hc_0);
    tfr_hc_45 = ft_freqanalysis(cfg,tfr_hc_45);
    tfr_hc_90 = ft_freqanalysis(cfg,tfr_hc_90);
    tfr_hc_115 = ft_freqanalysis(cfg,tfr_hc_115);
    tfr_hc_all = ft_freqanalysis(cfg,tfr_hc_all);

    %% Save data
    cd(datapath)
    save tfr_lc tfr_lc_0 tfr_lc_45 tfr_lc_90 tfr_lc_115 tfr_lc_all
end