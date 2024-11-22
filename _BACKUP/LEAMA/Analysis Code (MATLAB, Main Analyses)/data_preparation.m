%% Preprocessing Script for OCC Gamma

%% Description %%
% This script merges the blocks to a single data set for all participants, and adds trial information
% corresponding to the different contrast and orientation conditions. Furthermore, preprocessing is completed
% by changing the reference to the average across electrodes.

%% Script Structure
% 0. Preparation of paths and participant IDs, as well as EEGLAB
% 1. Appending all blocks with good or okay EEG data quality, and convert to a fieldtrip structure
% 2. Saving eye tracker and EEG data separately
% 3. Re-referencing into common average
% 4. Saving data in participant folder



%% 0. Preparation

%% 0.1 Path and subjects
subjects = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11'};
path = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results\';

%% 0.2 Trigger information

% High Contrast: 21-24
% 21 = horizontal (0 deg), 22 = vertical (90  deg), 23 = 45 degrees, 24 = 115 degrees

% Low Contrast: 25-28
% 25 = horizontal (0 deg), 26 = vertical (90  deg), 27 = 45 degrees, 28 = 115 degrees



for subj = 1:length(subjects)
    keep subj subjects path


    %% 0.3  Prepare eeglab (has to be done for each participant because of interfering FieldTrip installation below)
    addpath('\\psyger-stor02.d.uzh.ch\methlab\Students\Arne\MA\eeglab2022.1');
    eeglab
    ft_defaults;
    close all hidden


    %% 1. Append blocks

    % Prepare paths
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    disp(datapath) % For information about current participant

    %% 1.1 Find merged data, sort out bad blocks, and load
    filter = strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results\', num2str(subj), '\*ip_EEGmerged*.mat');
    MergedData = dir(filter(:));

    % Prepare and make sure the variables are empty
    counter = 0;
    xx = {};

    % Create structure containing all blocks
    for block = 1:length(MergedData)

        counter = counter + 1;

        % Select only blocks that are not 'bip' (bad quality)
        if startsWith(MergedData(block).name, 'bip') && block < 4
            block = block + 1;
        elseif startsWith(MergedData(block).name, 'bip') && block == 4
            break
        end

        load(MergedData(block).name)

        % Gather all blocks in one structure
        AllBlocks{block} = EEG;
        AllBlocks(cellfun('isempty', AllBlocks)) = [];

    end

    % Merge all blocks
    EEGcomplete = AllBlocks{1};
    for eeg = 1 : length(AllBlocks)
        EEGcomplete = pop_mergeset(EEGcomplete, AllBlocks{eeg}, 0);
    end

    %% 1.2 Segment data into epochs -2 before and 3.5 after stim onset

    % Triggers for high contrast: 21 hor, 22 vert, 23 45deg, 24 115deg
    EEG_hc_0deg = pop_epoch(EEGcomplete,{'21'},[-2 3.5]);
    EEG_hc_90deg = pop_epoch(EEGcomplete,{'22'},[-2 3.5]);
    EEG_hc_45deg = pop_epoch(EEGcomplete,{'23'},[-2 3.5]);
    EEG_hc_115deg = pop_epoch(EEGcomplete,{'24'},[-2 3.5]);

    % Triggers for low contrast: 25 hor, 26 vert, 27 45deg, 28 115deg
    EEG_lc_0deg = pop_epoch(EEGcomplete, {'25'}, [-2 3.5]);
    EEG_lc_90deg = pop_epoch(EEGcomplete, {'26'}, [-2 3.5]);
    EEG_lc_45deg = pop_epoch(EEGcomplete, {'27'}, [-2 3.5]);
    EEG_lc_115deg = pop_epoch(EEGcomplete, {'28'}, [-2 3.5]);

    %% 1.3 Convert into fieldtrip

    % Convert high contrast
    data0 = eeglab2fieldtrip(EEG_hc_0deg, 'raw');
    data0.trialinfo = zeros(numel(data0.trial),1);
    data45 = eeglab2fieldtrip(EEG_hc_45deg, 'raw');
    data45.trialinfo = zeros(numel(data45.trial),1) + 45;
    data90 = eeglab2fieldtrip(EEG_hc_90deg, 'raw');
    data90.trialinfo = zeros(numel(data90.trial),1) + 90;
    data115 = eeglab2fieldtrip(EEG_hc_115deg, 'raw');
    data115.trialinfo = zeros(numel(data115.trial),1) + 115;

    % Convert low contrast (trialinfo + 1 to avoid mixing the conditions)
    data1 = eeglab2fieldtrip(EEG_lc_0deg, 'raw');
    data1.trialinfo = ones(numel(data1.trial), 1); % Ones here to differentiate between contrasts
    data46 = eeglab2fieldtrip(EEG_lc_45deg, 'raw');
    data46.trialinfo = ones(numel(data46.trial), 1) + 45;
    data91 = eeglab2fieldtrip(EEG_lc_90deg, 'raw');
    data91.trialinfo = ones(numel(data91.trial), 1) + 90;
    data116 = eeglab2fieldtrip(EEG_lc_115deg, 'raw');
    data116.trialinfo = ones(numel(data116.trial) ,1) + 115;


    %% 1.4 Append

    % Add clean version of Fieldtrip (so it does not interfere with EEGLAB)
    cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin')
    run startup

    % Append all orientation conditions
    datahigh = ft_appenddata([], data0, data45, data90, data115);
    datalow = ft_appenddata([], data1, data46, data91, data116);

    % Append both condition data sets to one data set
    data = ft_appenddata([], datahigh, datalow);


    %% 2. Save eye tracker and EEG data separately

    % Save the eye tracker data
    cfg             = [];
    cfg.channel     = {'L-GAZE-X', 'L-GAZE-Y',  'L-AREA'}; % Extract "electrodes" that contain ET movement data
    dataet = ft_selectdata(cfg, data);

    dataet.event = EEGcomplete.event; % Add eye tracker events

    % Remove unnecessary electrodes and save EEG data separately
    cfg             = [];
    cfg.channel     = {'all'  '-B*' '-HEOGR', '-HEOGL', '-VEOGU', '-VEOGL', '-L-GAZE-X','-L-GAZE-Y','-L-AREA' }; % Ocular electrodes and ET movement data
    data = ft_selectdata(cfg, data);

    %% 3. Re-reference EEG data into common average

    % Select data set without the original reference CPz
    cfg             = [];
    cfg.channel     = {'all' '-CPz'};
    datarv = ft_selectdata(cfg, data);

    % Re-reference
    cfg             = [];
    cfg.reref       = 'yes';
    cfg.refchannel  = 'all';
    cfg.implicitref = 'CPz';
    data = ft_preprocessing(cfg, datarv);

    data.event = EEGcomplete.event; % Add eye tracker events


    %% 4. Save the data

    save(strcat(datapath, '\data.mat'), 'data', '-v7.3') % Save data in participant folder (-v7.3 is for making sure that it works for large volume)
    save(strcat(datapath, '\dataet.mat'), 'dataet') % Save dataet in participant folder

end % End participant loop

