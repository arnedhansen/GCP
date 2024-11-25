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
    load dataEEG
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions 
    % Horizontal
    ind21 = find(dataEEG_horz_lc.trialinfo == 21);
    ind22 = find(dataEEG_horz_hc.trialinfo == 22);

    % Vertical
    ind31 = find(dataEEG_vert_lc.trialinfo == 31);
    ind32 = find(dataEEG_vert_hc.trialinfo == 32);

    % Concentric Static
    ind41 = find(dataEEG_concentric_static_lc.trialinfo == 41);
    ind42 = find(dataEEG_concentric_static_hc.trialinfo == 42);

    % Concentric Dynamic Inward
    ind51 = find(dataEEG_concentric_dynamic_inward_lc.trialinfo == 51);
    ind52 = find(dataEEG_concentric_dynamic_inward_hc.trialinfo == 52);

    % Concentric Dynamic Outward
    ind61 = find(dataEEG_concentric_dynamic_outward_lc.trialinfo == 61);
    ind62 = find(dataEEG_concentric_dynamic_outward_hc.trialinfo == 62);

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

    % low contrast horizontal
    cfg.trials = ind21;
    tfr_horz_lc = ft_freqanalysis(cfg,dataEEG_horz_lc);

    % high contrast horizontal
    cfg.trials = ind22;
    tfr_horz_hc = ft_freqanalysis(cfg,dataEEG_horz_hc);

    % low contrast vertical
    cfg.trials = ind31;
    tfr_vert_lc = ft_freqanalysis(cfg,dataEEG_vert_lc);

    % high contrast vertical
    cfg.trials = ind32;
    tfr_vert_hc = ft_freqanalysis(cfg,dataEEG_vert_hc);

    % low contrast concentric static
    cfg.trials = ind41;
    tfr_concentric_static_lc = ft_freqanalysis(cfg,dataEEG_concentric_static_lc);

    % high contrast concentric static
    cfg.trials = ind42;
    tfr_concentric_static_hc = ft_freqanalysis(cfg,dataEEG_concentric_static_hc);

    % low contrast concentric dynamic inward
    cfg.trials = ind51;
    tfr_concentric_dynamic_inward_lc = ft_freqanalysis(cfg,dataEEG_concentric_dynamic_inward_lc);

    % high contrast concentric dynamic inward
    cfg.trials = ind52;
    tfr_concentric_dynamic_inward_hc = ft_freqanalysis(cfg,dataEEG_concentric_dynamic_inward_hc);

    % low contrast concentric dynamic outward
    cfg.trials = ind61;
    tfr_concentric_dynamic_outward_lc = ft_freqanalysis(cfg,dataEEG_concentric_dynamic_outward_lc);

    % high contrast concentric dynamic outward
    cfg.trials = ind62;
    tfr_concentric_dynamic_outward_hc = ft_freqanalysis(cfg,dataEEG_concentric_dynamic_outward_hc);

    %% Baseline
    cfg             = [];
    cfg.baseline    = [-.5 -.25];
    cfg.baselinetype = 'db';
    tfr_horz_lc_bl                       = ft_freqbaseline(cfg, tfr_horz_lc);
    tfr_horz_hc_bl                       = ft_freqbaseline(cfg, tfr_horz_hc);
    tfr_vert_lc_bl                       = ft_freqbaseline(cfg, tfr_vert_lc);
    tfr_vert_hc_bl                       = ft_freqbaseline(cfg, tfr_vert_hc);
    tfr_concentric_static_lc_bl          = ft_freqbaseline(cfg, tfr_concentric_static_lc);
    tfr_concentric_static_hc_bl          = ft_freqbaseline(cfg, tfr_concentric_static_hc);
    tfr_concentric_dynamic_inward_lc_bl  = ft_freqbaseline(cfg, tfr_concentric_dynamic_inward_lc);
    tfr_concentric_dynamic_inward_hc_bl  = ft_freqbaseline(cfg, tfr_concentric_dynamic_inward_hc);
    tfr_concentric_dynamic_outward_lc_bl = ft_freqbaseline(cfg, tfr_concentric_dynamic_outward_lc);
    tfr_concentric_dynamic_outward_hc_bl = ft_freqbaseline(cfg, tfr_concentric_dynamic_outward_hc);

    %% Save data
    cd(datapath)
    save data_tfr ...
        tfr_horz_lc  tfr_horz_lc_bl ...
        tfr_horz_hc  tfr_horz_hc_bl ...
        tfr_vert_lc  tfr_vert_lc_bl ...
        tfr_vert_hc  tfr_vert_hc_bl ...
        tfr_concentric_static_lc  tfr_concentric_static_lc_bl ...
        tfr_concentric_static_hc  tfr_concentric_static_hc_bl ...
        tfr_concentric_dynamic_inward_lc  tfr_concentric_dynamic_inward_lc_bl ...
        tfr_concentric_dynamic_inward_hc  tfr_concentric_dynamic_inward_hc_bl ...
        tfr_concentric_dynamic_outward_lc  tfr_concentric_dynamic_outward_lc_bl ...
        tfr_concentric_dynamic_outward_hc  tfr_concentric_dynamic_outward_hc_bl

    fprintf('Subject %.3d/%.3d TFR DATA computed \n', subj, length(subjects))
end