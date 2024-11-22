%% Preprocessing GCP

% Description
% This script merges the blocks to a single data set for all participants, 
% and adds trial information corresponding to the different contrast and 
% orientation conditions. Furthermore, preprocessing is completed by 
% changing the reference to the average across electrodes.

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.2');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/merged/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    clearvars -except subjects subj path
    datapath = strcat(path,subjects{subj});
    cd(datapath)

    %% Read blocks
    for block = 1:4
        try % Do not load emtpy blocks
            load(strcat(subjects{subj}, '_EEG_ET_GCP_block',num2str(block),'_merged.mat'))
            alleeg{block} = EEG;
            clear EEG
            fprintf('Subject GCP %.3s (%.3d/%.3d): Block %.1d loaded \n', subjects{subj}, subj, length(subjects), block)
        end
    end

    %% Segment data into epochs -2s before and 3.5s after stim onset and 
    %  convert to Fieldtrip data structure
    for block = 1:4
        try % Do not segment empty blocks

            % High Contrast: 21-24
            % 21 = horizontal (0 deg), 22 = vertical (90  deg), 23 = 45 degrees, 24 = 115 degrees
            EEG_hc_0deg = pop_epoch(alleeg{block},{'21'},[-2 3.5]);
            data0_hc{block} = eeglab2fieldtrip(EEG_hc_0deg, 'raw');
            EEG_hc_90deg = pop_epoch(alleeg{block},{'22'},[-2 3.5]);
            data90_hc{block} = eeglab2fieldtrip(EEG_hc_90deg, 'raw');
            EEG_hc_45deg = pop_epoch(alleeg{block},{'23'},[-2 3.5]);
            data45_hc{block} = eeglab2fieldtrip(EEG_hc_45deg, 'raw');
            EEG_hc_115deg = pop_epoch(alleeg{block},{'24'},[-2 3.5]);
            data115_hc{block} = eeglab2fieldtrip(EEG_hc_115deg, 'raw');

            % Low Contrast: 25-28
            % 25 = horizontal (0 deg), 26 = vertical (90  deg), 27 = 45 degrees, 28 = 115 degrees
            EEG_lc_0deg = pop_epoch(alleeg{block}, {'25'}, [-2 3.5]);
            data0_lc{block} = eeglab2fieldtrip(EEG_lc_0deg, 'raw');
            EEG_lc_90deg = pop_epoch(alleeg{block}, {'26'}, [-2 3.5]);
            data90_lc{block} = eeglab2fieldtrip(EEG_lc_90deg, 'raw');
            EEG_lc_45deg = pop_epoch(alleeg{block}, {'27'}, [-2 3.5]);
            data45_lc{block} = eeglab2fieldtrip(EEG_lc_45deg, 'raw');
            EEG_lc_115deg = pop_epoch(alleeg{block}, {'28'}, [-2 3.5]);
            data115_lc{block} = eeglab2fieldtrip(EEG_lc_115deg, 'raw');
        end
    end

    %% Remove empty blocks
    data0_hc = data0_hc(~cellfun('isempty', data0_hc));
    data45_hc = data45_hc(~cellfun('isempty', data45_hc));
    data90_hc = data90_hc(~cellfun('isempty', data90_hc));
    data115_hc = data115_hc(~cellfun('isempty', data115_hc));
    
    data0_lc = data0_lc(~cellfun('isempty', data0_lc));
    data45_lc = data45_lc(~cellfun('isempty', data45_lc));
    data90_lc = data90_lc(~cellfun('isempty', data90_lc));
    data115_lc = data115_lc(~cellfun('isempty', data115_lc));

    %% Equalize labels
    update_labels(data0_hc);
    update_labels(data90_hc);
    update_labels(data45_hc);
    update_labels(data115_hc);
    update_labels(data0_lc);
    update_labels(data45_lc);
    update_labels(data90_lc);
    update_labels(data115_lc);

    %% Append data for orientation conditions
    cfg = [];
    cfg.keepsampleinfo = 'no';   
    data0_hc = ft_appenddata(cfg, data0_hc{:});
    data45_hc = ft_appenddata(cfg, data45_hc{:});
    data90_hc = ft_appenddata(cfg, data90_hc{:});
    data115_hc = ft_appenddata(cfg, data115_hc{:});
    data0_lc = ft_appenddata(cfg, data0_lc{:});
    data45_lc = ft_appenddata(cfg, data45_lc{:});
    data90_lc = ft_appenddata(cfg, data90_lc{:});
    data115_lc = ft_appenddata(cfg, data115_lc{:});

    %% Add trialinfo 
    % High contrast
    data0_hc.trialinfo = zeros(numel(data0_hc.trial),1);
    data45_hc.trialinfo = zeros(numel(data45_hc.trial),1) + 45;
    data90_hc.trialinfo = zeros(numel(data90_hc.trial),1) + 90;
    data115_hc.trialinfo = zeros(numel(data115_hc.trial),1) + 115;

    % Low contrast (trialinfo + 1000 to avoid mixing the conditions)
    data0_lc.trialinfo = zeros(numel(data0_lc.trial), 1) + 1000;
    data45_lc.trialinfo = zeros(numel(data45_lc.trial), 1) + 1045;
    data90_lc.trialinfo = zeros(numel(data90_lc.trial), 1) + 1090;
    data115_lc.trialinfo = zeros(numel(data115_lc.trial) ,1) + 1115;

    %% Append data for contrast conditions
    cfg = [];
    cfg.keepsampleinfo = 'no';    
    dataHC = ft_appenddata([], data0_hc, data45_hc, data90_hc, data115_hc);
    dataLC = ft_appenddata([], data0_lc, data45_lc, data90_lc, data115_lc);

    %% Get EyeTracking data
    cfg = [];
    cfg.channel = {'L-GAZE-X'  'L-GAZE-Y' 'L-AREA'};
    dataethc = ft_selectdata(cfg,dataHC);
    dataetlc = ft_selectdata(cfg,dataLC);

    %% Get EEG data (excl. ET and EOG data)
    cfg = [];
    cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' , '-L-AREA'};
    dataHC = ft_selectdata(cfg,dataHC);
    dataLC = ft_selectdata(cfg,dataLC);

    %% Re-segment data to avoid filter ringing
    cfg = [];
    dataTFRhc = ft_selectdata(cfg,dataHC); % TRF data long HC
    dataTFRlc = ft_selectdata(cfg,dataLC); % TRF data long LC
    cfg = [];
    cfg.latency = [0 2]; % Time window for grating analysis
    dataHC = ft_selectdata(cfg,dataHC); % EEG data HC
    dataethc = ft_selectdata(cfg,dataethc); % ET data HC
    dataLC = ft_selectdata(cfg,dataLC); % EEG data LC
    dataetlc = ft_selectdata(cfg,dataetlc); % ET data LC

    %% Re-reference data to average or common reference
    cfg = [];
    cfg.reref   = 'yes';
    cfg.refchannel = 'all';
    dataHC = ft_preprocessing(cfg,dataHC);
    dataLC = ft_preprocessing(cfg,dataLC);

    %% Save data
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save dataEEGhc dataHC
    save dataEEGlc dataLC
    save dataEEG_TFRhc dataTFRhc
    save dataEEG_TFRlc dataTFRlc
    savepathET = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/gaze/');
    mkdir(savepathET)
    cd(savepathET)
    save dataEThc dataethc
    save dataETlc dataetlc
    clc
    if subj == length(subjects)
        disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. PREPROCESSING FINALIZED.'])
    else
        disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. Loading next subject...'])
    end
end

%% Function to update labels
function update_labels(data)
blocks = size(data);
for block = 1:blocks
    if isempty(data{block})
        break;
    else
        try
            for i = 1:blocks
                if ~isempty(data{i}.label)
                    data{block}.label = data{i}.label;
                    break;
                end
            end
        catch
            warning('Error occurred while processing block %d in data structure.', block);
        end
    end
end
end