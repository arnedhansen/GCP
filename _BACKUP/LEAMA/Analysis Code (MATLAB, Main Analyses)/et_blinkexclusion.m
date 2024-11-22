%% ET: Exclusion of Trials With Closed Eye Periods %%

%% Description %%
% This script produces data sets with eye tracker data excluding the trials that contain a) blinks
% (as defined by the eye trracking device) and b) longer periods of closed eyes.

%% Script Structure
% 0. Preparation
% 1. Blink exclusion
% 2. Closed Eye Periods Exclusion (lenient blink exclusion)



%% 0. Preparation

clear all
subjects = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11'};
path = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results\';


%% 1. Exclusion of all trials with blinks

%% 1.1 Find indices of blink trials

for subj = 1:length(subjects)

    keep subj subjects path amount_blktrials_low amount_blktrials_high amount_trials_low amount_trials_high

    datapath = strcat(path, subjects{subj});
    cd(datapath)

    %% 1.1.1  Prepare eeglab and FieldTrip (has to be done for each subject because of fieldtrip installation below)

    addpath('\\psyger-stor02.d.uzh.ch\methlab\Students\Arne\MA\eeglab2022.1');
    eeglab
    ft_defaults;
    close all hidden


    %% 1.1.2 Find merged data, sort out bad blocks, and load

    filter = strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results\', num2str(subj), '\*ip_EEGmerged*.mat');
    MergedData = dir(filter(:));

    % Prepare and make sure variables are zero
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

    % Merge all blocks together
    EEGcomplete = AllBlocks{1};

    for eeg = 1 : length(AllBlocks)
        EEGcomplete = pop_mergeset(EEGcomplete, AllBlocks{eeg}, 0);
    end



    %% 1.1.3 Segment data into epochs -2 before and 3.5 after stim onset

    % Triggers for high contrast: 21 hor, 22 vert, 23 45deg , 24 115deg
    EEG_hc = pop_epoch(EEGcomplete,{'21', '22', '23', '24'}, [-2 3.5]);
    EEG_lc = pop_epoch(EEGcomplete,{'25', '26', '27', '28'}, [-2 3.5]);


    %% 1.2 Find trials that contain blinks registered by eye tracker

    % Blink indices

    blinkind_high =  find(strcmp({EEG_hc.event.type},'L_fixation'));
    blinktrials_high = [EEG_hc.event(blinkind_high).epoch];
    blinktrials_high = unique(blinktrials_high);

    blinkind_low =  find(strcmp({EEG_lc.event.type},'L_fixation'));
    blinktrials_low = [EEG_lc.event(blinkind_low).epoch];
    blinktrials_low = unique(blinktrials_low);

    % Count for later percentage calculation
    amount_blktrials_low{subj} = length(blinktrials_low);
    amount_trials_low{subj} = EEG_lc.trials;
    amount_blktrials_high{subj} = length(blinktrials_high);
    amount_trials_high{subj} = EEG_hc.trials;

    % Deselect blink trials (low contrast)
    if length(blinktrials_low) ~= length(EEG_lc.epoch) % If not all the trials are blink trials
        EEG_lc = pop_select(EEG_lc, 'notrial', blinktrials_high);
        data_low = eeglab2fieldtrip(EEG_lc, 'raw');
        data_low.trialinfo = ones(numel(data_low.trial),1);
    end

    % Deselect blink trials (high contrast)
    if length(blinktrials_high) ~= length(EEG_hc.epoch) % If not all the trials are blink trials
        EEG_hc = pop_select(EEG_hc, 'notrial', blinktrials_high);
        data_high = eeglab2fieldtrip(EEG_hc, 'raw');
        data_high.trialinfo = zeros(numel(data_high.trial),1);
    end


    %% 1.3 Append data

    cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin')
    run startup

    % If data exists, append to one data set
    if exist('data_high') && exist('data_low')
        data_noblinkET = ft_appenddata([], data_high, data_low);
    elseif exist('data_high')
        data_noblinkET = ft_appenddata([], data_high);
    elseif exist('data_low')
        data_noblinkET = ft_appenddata([], data_low);
    end

    % Re-reference and save
    if exist('data_noblinkET')

        cfg             = [];
        cfg.channel     = {'all' '-CPz'};
        datarv = ft_selectdata(cfg, data_noblinkET);

        cfg             = [];
        cfg.reref       = 'yes';
        cfg.refchannel  = 'all';
        cfg.implicitref = 'CPz';
        data_noblinkET = ft_preprocessing(cfg, datarv);

        %% Save
        save(strcat(datapath, '\data_noblinkET.mat'), 'data_noblinkET', '-v7.3')

    end
end


%% 1.4 Calculate percentages

% Percentages
for n = 1:length(amount_trials_high)
    perctge_low(n) = (cell2mat(amount_blktrials_low(n))/cell2mat(amount_trials_low(n)))*100;
    perctge_high(n) = (cell2mat(amount_blktrials_high(n))/cell2mat(amount_trials_high(n)))*100;
end

% Mean percentage of data lost and standard deviation
mean(perctge_low)
mean(perctge_high)
std(perctge_low)
std(perctge_high)



%% 1.5 Fast Fourier Transform D

subjects = {'1' '2' '3' '6' '8' '9' '11'};
path = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results\';

% FieldTrip
cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin')
run startup

% FFT
for subj = 1:length(subjects)

    %% 1.5.1 Load data
    keep subj subjects path
    datapath = strcat(path, subjects{subj});
    cd(datapath)

    load data_noblinkET

    indhc = find(data_noblinkET.trialinfo == 0); % High contrast
    indlc = find(data_noblinkET.trialinfo == 1); % Low contrast


    %% 1.5.2 Time frequency analysis (single hanning taper)

    % Preparation
    cfg             = [];
    cfg.output      = 'pow';
    cfg.method      = 'mtmconvol';
    cfg.taper       = 'hanning';
    cfg.foi         = 4:1:120;                        % analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin   = ones(length(cfg.foi), 1).*0.5;  % length of time window = 0.5 sec (not adapted to the respective frequency!)
    cfg.toi         = -1:0.05:3;                      % what time window we are interested in, rsp. to stim presentation
    cfg.keeptrials  = 'no';

    % FFT
    cfg.trials = [indhc];
    tfr_hc_nbET = ft_freqanalysis(cfg, data_noblinkET);
    cfg.trials = [indlc];
    tfr_lc_nbET = ft_freqanalysis(cfg, data_noblinkET);

    tfr_lc_nbET.cfg = [];
    tfr_hc_nbET.cfg = [];

    % Save
    save tfr_lc_nbET tfr_lc_nbET
    save tfr_hc_nbET tfr_hc_nbET

end



%% 2. Exclude trials with 10 or more timepoints with no eye tracker information

%% 2.1 Get indices for trials without blinks from all subjects

for subj = 1:length(subjects)

    datapath = strcat(path, subjects{subj});
    cd(datapath)
    disp(datapath)
    load dataet

    counter  = 0;

    for tr = 1:length(dataet.trial)

        % Find indices for all timepoints in trial that are not 0 (= blink)
        i = find(dataet.trial{1,tr}(1,:) == 0);

        % Remove all indices before stimuli presentation
        i(i < 1001) =  [];

        % Save remaining indic
        if length(i) < 10
            counter = counter + 1;
            trind(subj).ind(counter) = tr;
        end

        % Total trial number
        total(subj) = length(dataet.trial);
    end
end

%% 2.2 Calculate percentage lost

for subj = 1:length(subjects)
    left(subj) = (length(trind(subj).ind)/total(subj));
end

lost = 1-left(:);

mean(lost)
std(lost)


%% 1.3. Get rid of indexed trials and save as data_noblink

keep subjects path trind
cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\')
run startup

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    disp(datapath)
    load data

    % Select correct index and use it for selecting trials
    index = trind(subj);
    cfg.trials = struct2array(index);
    data_noblink = ft_preprocessing(cfg, data);

    % Save as new data
    save data_noblink data_noblink
end


