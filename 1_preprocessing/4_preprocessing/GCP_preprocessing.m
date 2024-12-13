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
for subj = 1:length(subjects)
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
                EEG_lc = pop_epoch(alleeg{block}, {'61'}, epoch_window);
                data_lc{block} = eeglab2fieldtrip(EEG_lc, 'raw');
                % Low contrast gaze metrics extraction
                lc_gaze_metrics = pop_epoch(alleeg{block}, {'61'}, analysis_window); % Stimulus presentation time
                lc_trl(block) = lc_gaze_metrics.trials;
                lc_sacc(block) = sum(ismember({lc_gaze_metrics.event.type}, {'L_saccade', 'R_saccade'}));
                lc_fix(block) = sum(ismember({lc_gaze_metrics.event.type}, {'L_fixation', 'R_fixation'}));
                lc_blink(block) = sum(ismember({lc_gaze_metrics.event.type}, {'L_blink', 'R_blink'}));

                EEG_hc = pop_epoch(alleeg{block}, {'62'}, epoch_window);
                data_hc{block} = eeglab2fieldtrip(EEG_hc, 'raw');
                % High contrast gaze metrics extraction
                hc_gaze_metrics = pop_epoch(alleeg{block}, {'62'}, analysis_window); % Stimulus presentation time
                hc_trl(block) = hc_gaze_metrics.trials;
                hc_sacc(block) = sum(ismember({hc_gaze_metrics.event.type}, {'L_saccade', 'R_saccade'}));
                hc_fix(block) = sum(ismember({hc_gaze_metrics.event.type}, {'L_fixation', 'R_fixation'}));
                hc_blink(block) = sum(ismember({hc_gaze_metrics.event.type}, {'L_blink', 'R_blink'}));
            end
        end

        % %% Remove empty blocks
        % data_lc = data_lc(~cellfun('isempty', data_lc));
        % data_hc = data_hc(~cellfun('isempty', data_hc));
        % 
        % %% Equalize labels
        % update_labels(data_lc);
        % update_labels(data_hc);
        % 
        % %% Append data for conditions
        % cfg = [];
        % cfg.keepsampleinfo = 'no';
        % data_lc = ft_appenddata(cfg, data_lc{:});
        % data_hc = ft_appenddata(cfg, data_hc{:});
        % 
        % %% Add trialinfo
        % data_lc.trialinfo = zeros(numel(data_lc.trial), 1) + 61;
        % data_hc.trialinfo = zeros(numel(data_hc.trial), 1) + 62;
        % 
        % %% Get EyeTracking data
        % cfg = [];
        % cfg.channel = {'L-GAZE-X'  'L-GAZE-Y' 'L-AREA', 'R-GAZE-X'  'R-GAZE-Y' 'R-AREA'};
        % dataET_lc = ft_selectdata(cfg, data_lc);
        % dataET_hc = ft_selectdata(cfg, data_hc);
        % 
        % %% Get EEG data (excl. ET and EOG data)
        % cfg = [];
        % cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' , '-L-AREA', '-R-GAZE-X'  '-R-GAZE-Y' '-R-AREA'};
        % dataEEG_lc = ft_selectdata(cfg, data_lc);
        % dataEEG_hc = ft_selectdata(cfg, data_hc);
        % 
        % %% Re-reference data to average/common reference
        % cfg = [];
        % cfg.reref   = 'yes';
        % cfg.refchannel = 'all';
        % dataEEG_lc = ft_preprocessing(cfg, dataEEG_lc);
        % dataEEG_hc = ft_preprocessing(cfg, dataEEG_hc);

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
        % savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/eeg/');
        % mkdir(savepath)
        % cd(savepath)
        % save dataEEG dataEEG_lc dataEEG_hc
        savepathET = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/gaze/');
        mkdir(savepathET)
        cd(savepathET)
        %save dataET dataET_lc dataET_hc
        save gaze_metrics lc_saccades lc_fixations lc_blinks hc_saccades hc_fixations hc_blinks
        clc
        if subj == length(subjects)
            disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. PREPROCESSING FINALIZED.'])
        else
            disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. Loading next subject...'])
        end
    %end
end