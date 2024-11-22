%% GCP EEG Feature Extraction
%
% Extracted features:
%   TFR

%% Setup
[subjects, path] = setup('GCP');

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
    cfg.foi         = 30:4:120;                         % Analysis 30 to 90 Hz in steps of 4 Hz
    cfg.tapsmofrq   = 11;                               % 10 tapers
    cfg.t_ftimwin   = ones(length(cfg.foi),1).*0.5;     % Length of time window = 0.5 sec
    cfg.toi         = -1:0.05:3;                       
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
    cfg.baseline    = [-.5 -.25];
    cfg.baselinetype = 'db';
    tfr_hc_0_bl = ft_freqbaseline(cfg,tfr_hc_0);
    tfr_hc_45_bl = ft_freqbaseline(cfg,tfr_hc_45);
    tfr_hc_90_bl = ft_freqbaseline(cfg,tfr_hc_90);
    tfr_hc_115_bl = ft_freqbaseline(cfg,tfr_hc_115);
    tfr_hc_all_bl = ft_freqbaseline(cfg,tfr_hc_all);

    %% Save data
    cd(datapath)
    save data_tfr ...
        tfr_hc_0 tfr_hc_0_bl ...
        tfr_hc_45 tfr_hc_45_bl ...
        tfr_hc_90 tfr_hc_90_bl ...
        tfr_hc_115 tfr_hc_115_bl ...
        tfr_hc_all tfr_hc_all_bl
    fprintf('Subject %.2s / %.3d TFR DATA computed \n', num2str(subj), length(subjects))
end