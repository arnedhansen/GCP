%% Preprocessing GCP

%% Setup
startup
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
for subj = 6:length(subjects)
    clearvars -except subjects subj path
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    %if isempty(dir(['/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/eeg/dataEEG.mat']))
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
    %  convert to Fieldtrip data structure AND extract gaze metrics from raw EEG data
    epoch_window = [-2 3.5];
    analysis_window = [0.3 2]; % Analysis window for eye metric extraction
    for block = 1:4
        % 51 = Trigger for presentation of low contrast concentric dynamic inward grating WITH button press response (PRESENTATION_LC_TASK)
        % 52 = Trigger for presentation of high contrast concentric dynamic inward grating WITH button press response (PRESENTATION_HC_TASK)
        % 61 = Trigger for presentation of low contrast concentric dynamic inward grating WITHOUT button press response (PRESENTATION_LC_NOTASK)
        % 62 = Trigger for presentation of high contrast concentric dynamic inward grating WITHOUT button press response (PRESENTATION_HC_NOTASK)
        try
            % Low contrast EEG data (trigger = 61)
            EEG_lc = pop_epoch(alleeg{block}, {'61'}, epoch_window);
            data_lc{block} = eeglab2fieldtrip(EEG_lc, 'raw');

            % Low contrast gaze metrics extraction
            lc_gaze_metrics = pop_epoch(alleeg{block}, {'61'}, analysis_window);
            lc_trl(block) = lc_gaze_metrics.trials;

            % Extract blink timepoints
            blink_times = [lc_gaze_metrics.event(strcmp({lc_gaze_metrics.event.type}, 'L_blink') | strcmp({lc_gaze_metrics.event.type}, 'R_blink')).latency];

            % Extract saccades timepoints
            saccade_events = lc_gaze_metrics.event(strcmp({lc_gaze_metrics.event.type}, 'L_saccade') | strcmp({lc_gaze_metrics.event.type}, 'R_saccade'));

            % Exclude saccades around blinks
            valid_saccades = 0;
            for s = 1:length(saccade_events)
                saccade_time = saccade_events(s).latency;
                % Check if this saccade is within 100 ms of any blink
                near_blink = any(abs(saccade_time - blink_times) <= 50); % 50 samples = 100 ms
                if ~near_blink
                    valid_saccades = valid_saccades + 1;
                end
            end

            % Count low contrast gaze metrics
            lc_sacc(block) = valid_saccades;
            lc_fix(block) = sum(ismember({lc_gaze_metrics.event.type}, {'L_fixation', 'R_fixation'}));
            lc_blink(block) = numel(blink_times);

            % High contrast EEG data (trigger = 62)
            EEG_hc = pop_epoch(alleeg{block}, {'62'}, epoch_window);
            data_hc{block} = eeglab2fieldtrip(EEG_hc, 'raw');

            % High contrast gaze metrics extraction
            hc_gaze_metrics = pop_epoch(alleeg{block}, {'62'}, analysis_window);
            hc_trl(block) = hc_gaze_metrics.trials;

            % Extract blink timepoints
            blink_times = [hc_gaze_metrics.event(strcmp({hc_gaze_metrics.event.type}, 'L_blink') | strcmp({hc_gaze_metrics.event.type}, 'R_blink')).latency];

            % Extract saccades timepoints
            saccade_events = hc_gaze_metrics.event(strcmp({hc_gaze_metrics.event.type}, 'L_saccade') | strcmp({hc_gaze_metrics.event.type}, 'R_saccade'));

            % Exclude saccades around blinks
            valid_saccades = 0;
            for s = 1:length(saccade_events)
                saccade_time = saccade_events(s).latency;
                % Check if this saccade is within 100 ms of any blink
                near_blink = any(abs(saccade_time - blink_times) <= 50); % 50 samples = 100 ms
                if ~near_blink
                    valid_saccades = valid_saccades + 1;
                end
            end

            % Count high contrast gaze metrics
            hc_sacc(block) = valid_saccades;
            hc_fix(block) = sum(ismember({hc_gaze_metrics.event.type}, {'L_fixation', 'R_fixation'}));
            hc_blink(block) = numel(blink_times);
        end
    end

    %% Remove empty blocks
    data_lc = data_lc(~cellfun('isempty', data_lc));
    data_hc = data_hc(~cellfun('isempty', data_hc));

    %% Equalize labels
    update_labels(data_lc);
    update_labels(data_hc);

    %% Add trialinfo
    for block = 1:4
        try
            data_lc{block}.trialinfo = zeros(numel(data_lc{block}.trial), 1) + 61;
            data_hc{block}.trialinfo = zeros(numel(data_hc{block}.trial), 1) + 62;
        end
    end

    %% Append data for conditions
    cfg = [];
    cfg.keepsampleinfo = 'no';
    data_lc = ft_appenddata(cfg, data_lc{:});
    data_hc = ft_appenddata(cfg, data_hc{:});

    %% Get EyeTracking data
    cfg = [];
    cfg.channel = {'L-GAZE-X'  'L-GAZE-Y' 'L-AREA', 'R-GAZE-X'  'R-GAZE-Y' 'R-AREA'};
    dataET_lc = ft_selectdata(cfg, data_lc);
    dataET_hc = ft_selectdata(cfg, data_hc);

    %% Get EEG data (excl. ET and EOG data)
    cfg = [];
    cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' , '-L-AREA', '-R-GAZE-X'  '-R-GAZE-Y' '-R-AREA'};
    dataEEG_lc = ft_selectdata(cfg, data_lc);
    dataEEG_hc = ft_selectdata(cfg, data_hc);

    %% Re-reference data to average/common reference
    cfg = [];
    cfg.reref   = 'yes';
    cfg.refchannel = 'all';
    dataEEG_lc = ft_preprocessing(cfg, dataEEG_lc);
    dataEEG_hc = ft_preprocessing(cfg, dataEEG_hc);

    %% Compute gaze metric data
    % Low contrast gaze metrics average across trials
    lc_saccades = sum(lc_sacc(:))/sum(lc_trl(:));
    lc_fixations = sum(lc_fix(:))/sum(lc_trl(:));
    lc_blinks = sum(lc_blink(:))/sum(lc_trl(:));

    % High contrast gaze metrics average across trials
    hc_saccades = sum(hc_sacc(:))/sum(hc_trl(:));
    hc_fixations = sum(hc_fix(:))/sum(hc_trl(:));
    hc_blinks = sum(hc_blink(:))/sum(hc_trl(:));

    %% Save data
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save dataEEG dataEEG_lc dataEEG_hc
    savepathET = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/gaze/');
    mkdir(savepathET)
    cd(savepathET)
    save dataET dataET_lc dataET_hc
    save gaze_metrics lc_saccades lc_fixations lc_blinks hc_saccades hc_fixations hc_blinks
    clc
    if subj == length(subjects)
        disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. PREPROCESSING FINALIZED.'])
    else
        disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. Loading next subject...'])
    end
    %end
end