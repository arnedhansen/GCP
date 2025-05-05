%% Preprocessing GCP

%% Setup
startup
clear
addEEGLab
path = '/Volumes/methlab/Students/Arne/GCP/data/merged/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    clearvars -except subjects subj path
    datapath = strcat(path,subjects{subj});
    cd(datapath)

    if isempty(dir(['/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/eeg/dataEEG.mat']))
        %% Read blocks
        for block = 1:4
            try % Do not load emtpy blocks
                load(strcat(subjects{subj}, '_EEG_ET_GCP_block',num2str(block),'_merged.mat'))
                alleeg{block} = EEG;
                clear EEG
                fprintf('Subject GCP %.3s (%.3d/%.3d): Block %.1d loaded \n', subjects{subj}, subj, length(subjects), block)
            catch ME
                ME.message
                disp(['ERROR loading Block ' num2str(block) '!'])
            end
        end

        %% Segment data into epochs -2s before and 3.5s after stim onset and
        %  convert to Fieldtrip data structure AND extract gaze metrics from raw EEG data
        epoch_window = [-2 3.5];
        analysis_window = [0.3 2]; % Analysis window for eye metric extraction
        for block = 1:4
            % 51 = PRESENTATION_C25_TASK    (Trigger for presentation of 25% contrast concentric dynamic inward grating WITH button press response)
            % 52 = PRESENTATION_C50_TASK    (Trigger for presentation of 50% contrast concentric dynamic inward grating WITH button press response)
            % 53 = PRESENTATION_C75_TASK    (Trigger for presentation of 75% contrast concentric dynamic inward grating WITH button press response)
            % 54 = PRESENTATION_C100_TASK   (Trigger for presentation of 100% contrast concentric dynamic inward grating WITH button press response)
            % 61 = PRESENTATION_C25_NOTASK  (Trigger for presentation of 25% contrast concentric dynamic inward grating WITHOUT button press response)
            % 62 = PRESENTATION_C50_NOTASK  (Trigger for presentation of 50% contrast concentric dynamic inward grating WITHOUT button press response)
            % 63 = PRESENTATION_C75_NOTASK  (Trigger for presentation of 75% contrast concentric dynamic inward grating WITHOUT button press response)
            % 64 = PRESENTATION_C100_NOTASK (Trigger for presentation of 100% contrast concentric dynamic inward grating WITHOUT button press response)
            try
                %% Segment 25% contrast data
                % 25% contrast EEG data (trigger = 61)
                EEG_c25 = pop_epoch(alleeg{block}, {'61'}, epoch_window);
                data_c25{block} = eeglab2fieldtrip(EEG_c25, 'raw');

                % 25% contrast gaze metrics extraction
                c25_gaze_metrics = pop_epoch(alleeg{block}, {'61'}, analysis_window);
                c25_trl(block) = c25_gaze_metrics.trials;
                % Extract blink timepoints
                blink_times = [c25_gaze_metrics.event(strcmp({c25_gaze_metrics.event.type}, 'L_blink') | strcmp({c25_gaze_metrics.event.type}, 'R_blink')).latency];
                % Extract saccades timepoints
                saccade_events = c25_gaze_metrics.event(strcmp({c25_gaze_metrics.event.type}, 'L_saccade') | strcmp({c25_gaze_metrics.event.type}, 'R_saccade'));

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

                % Count 25% contrast gaze metrics
                c25_sacc(block) = valid_saccades;
                c25_fix(block) = sum(ismember({c25_gaze_metrics.event.type}, {'L_fixation', 'R_fixation'}));
                c25_blink(block) = numel(blink_times);

                %% Segment 50% contrast data
                % 50% contrast EEG data (trigger = 62)
                EEG_c50 = pop_epoch(alleeg{block}, {'62'}, epoch_window);
                data_c50{block} = eeglab2fieldtrip(EEG_c50, 'raw');

                % 50% contrast gaze metrics extraction
                c50_gaze_metrics = pop_epoch(alleeg{block}, {'62'}, analysis_window);
                c50_trl(block) = c50_gaze_metrics.trials;

                % Extract blink timepoints
                blink_times = [c50_gaze_metrics.event(strcmp({c50_gaze_metrics.event.type}, 'L_blink') | strcmp({c50_gaze_metrics.event.type}, 'R_blink')).latency];

                % Extract saccades timepoints
                saccade_events = c50_gaze_metrics.event(strcmp({c50_gaze_metrics.event.type}, 'L_saccade') | strcmp({c50_gaze_metrics.event.type}, 'R_saccade'));

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

                % Count 50% contrast gaze metrics
                c50_sacc(block) = valid_saccades;
                c50_fix(block) = sum(ismember({c50_gaze_metrics.event.type}, {'L_fixation', 'R_fixation'}));
                c50_blink(block) = numel(blink_times);

                %% Segment 75% contrast data
                % 75% contrast EEG data (trigger = 63)
                EEG_c75 = pop_epoch(alleeg{block}, {'63'}, epoch_window);
                data_c75{block} = eeglab2fieldtrip(EEG_c75, 'raw');

                % 75% contrast gaze metrics extraction
                c75_gaze_metrics = pop_epoch(alleeg{block}, {'63'}, analysis_window);
                c75_trl(block) = c75_gaze_metrics.trials;

                % Extract blink timepoints
                blink_times = [c75_gaze_metrics.event(strcmp({c75_gaze_metrics.event.type}, 'L_blink') | strcmp({c75_gaze_metrics.event.type}, 'R_blink')).latency];

                % Extract saccades timepoints
                saccade_events = c75_gaze_metrics.event(strcmp({c75_gaze_metrics.event.type}, 'L_saccade') | strcmp({c75_gaze_metrics.event.type}, 'R_saccade'));

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

                % Count 75% contrast gaze metrics
                c75_sacc(block) = valid_saccades;
                c75_fix(block) = sum(ismember({c75_gaze_metrics.event.type}, {'L_fixation', 'R_fixation'}));
                c75_blink(block) = numel(blink_times);

                %% Segment 100% contrast data
                % 100% contrast EEG data (trigger = 64)
                EEG_c100 = pop_epoch(alleeg{block}, {'64'}, epoch_window);
                data_c100{block} = eeglab2fieldtrip(EEG_c100, 'raw');

                % 100% contrast gaze metrics extraction
                c100_gaze_metrics = pop_epoch(alleeg{block}, {'64'}, analysis_window);
                c100_trl(block) = c100_gaze_metrics.trials;

                % Extract blink timepoints
                blink_times = [c100_gaze_metrics.event(strcmp({c100_gaze_metrics.event.type}, 'L_blink') | strcmp({c100_gaze_metrics.event.type}, 'R_blink')).latency];

                % Extract saccades timepoints
                saccade_events = c100_gaze_metrics.event(strcmp({c100_gaze_metrics.event.type}, 'L_saccade') | strcmp({c100_gaze_metrics.event.type}, 'R_saccade'));

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
                c100_sacc(block) = valid_saccades;
                c100_fix(block) = sum(ismember({c100_gaze_metrics.event.type}, {'L_fixation', 'R_fixation'}));
                c100_blink(block) = numel(blink_times);
            catch ME
                ME.message
                disp(['ERROR segmenting Block ' num2str(block) '!'])
            end
        end

        %% Remove empty blocks
        data_c25 = data_c25(~cellfun('isempty', data_c25));
        data_c50 = data_c50(~cellfun('isempty', data_c50));
        data_c75 = data_c75(~cellfun('isempty', data_c75));
        data_c100 = data_c100(~cellfun('isempty', data_c100));

        %% Equalize labels
        update_labels(data_c25);
        update_labels(data_c50);
        update_labels(data_c75);
        update_labels(data_c100);

        %% Add trialinfo
        for block = 1:4
            try
                data_c25{block}.trialinfo = zeros(numel(data_c25{block}.trial), 1) + 61;
                data_c50{block}.trialinfo = zeros(numel(data_c50{block}.trial), 1) + 62;
                data_c75{block}.trialinfo = zeros(numel(data_c75{block}.trial), 1) + 63;
                data_c100{block}.trialinfo = zeros(numel(data_c100{block}.trial), 1) + 64;
            catch ME
                ME.message
                disp(['ERROR adding trialinfo in Block ' num2str(block) '!'])
            end
        end

        %% Append data for conditions
        cfg = [];
        cfg.keepsampleinfo = 'no';
        data_c25 = ft_appenddata(cfg, data_c25{:});
        data_c50 = ft_appenddata(cfg, data_c50{:});
        data_c75 = ft_appenddata(cfg, data_c75{:});
        data_c100 = ft_appenddata(cfg, data_c100{:});

        %% Select EyeTracking data
        cfg = [];
        cfg.channel = {'L-GAZE-X'  'L-GAZE-Y' 'L-AREA', 'R-GAZE-X'  'R-GAZE-Y' 'R-AREA'};
        dataET_c25 = ft_selectdata(cfg, data_c25);
        dataET_c50 = ft_selectdata(cfg, data_c50);
        dataET_c75 = ft_selectdata(cfg, data_c75);
        dataET_c100 = ft_selectdata(cfg, data_c100);

        %% Select EEG data (excl. ET and EOG data)
        cfg = [];
        cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' , '-L-AREA', '-R-GAZE-X'  '-R-GAZE-Y' '-R-AREA'};
        dataEEG_c25 = ft_selectdata(cfg, data_c25);
        dataEEG_c50 = ft_selectdata(cfg, data_c50);
        dataEEG_c75 = ft_selectdata(cfg, data_c75);
        dataEEG_c100 = ft_selectdata(cfg, data_c100);

        %% Re-reference data to average/common reference
        cfg = [];
        cfg.reref   = 'yes';
        cfg.refchannel = 'all';
        dataEEG_c25 = ft_preprocessing(cfg, dataEEG_c25);
        dataEEG_c50 = ft_preprocessing(cfg, dataEEG_c50);
        dataEEG_c75 = ft_preprocessing(cfg, dataEEG_c75);
        dataEEG_c100 = ft_preprocessing(cfg, dataEEG_c100);

        %% Compute gaze metric data
        % 25% contrast gaze metrics average across trials
        c25_saccades = sum(c25_sacc(:)) / sum(c25_trl(:));
        c25_fixations = sum(c25_fix(:)) / sum(c25_trl(:));
        c25_blinks = sum(c25_blink(:)) / sum(c25_trl(:));

        % 50% contrast gaze metrics average across trials
        c50_saccades = sum(c50_sacc(:)) / sum(c50_trl(:));
        c50_fixations = sum(c50_fix(:)) / sum(c50_trl(:));
        c50_blinks = sum(c50_blink(:)) / sum(c50_trl(:));

        % 75% contrast gaze metrics average across trials
        c75_saccades = sum(c75_sacc(:)) / sum(c75_trl(:));
        c75_fixations = sum(c75_fix(:)) / sum(c75_trl(:));
        c75_blinks = sum(c75_blink(:)) / sum(c75_trl(:));

        % 100% contrast gaze metrics average across trials
        c100_saccades = sum(c100_sacc(:)) / sum(c100_trl(:));
        c100_fixations = sum(c100_fix(:)) / sum(c100_trl(:));
        c100_blinks = sum(c100_blink(:)) / sum(c100_trl(:));

        %% Save data
        savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/eeg/');
        mkdir(savepath)
        cd(savepath)
        save dataEEG dataEEG_c25 dataEEG_c50 dataEEG_c75 dataEEG_c100
        savepathET = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/gaze/');
        mkdir(savepathET)
        cd(savepathET)
        save dataET dataET_c25 dataET_c50 dataET_c75 dataET_c100
        save gaze_metrics c25_saccades c25_fixations c25_blinks ...
            c50_saccades c50_fixations c50_blinks ...
            c75_saccades c75_fixations c75_blinks ...
            c100_saccades c100_fixations c100_blinks
        clc
        if subj == length(subjects)
            disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. PREPROCESSING FINALIZED.'])
        else
            disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. Loading next subject...'])
        end
    end
end