%% GCP Gaze Feature Extraction Sternberg
%
% Extracted features:
%   Gaze deviation (Euclidean distances)
%   Gaze standard deviation
%   Pupil size
%   Microsaccades
%
% Gaze metrics labelled by eye-tracker (saccades, blinks and
% fixations) are extracted already in GCP_preprocessing.m

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
gaze_data = struct('ID', {}, 'Condition', {}, 'GazeDeviation', {}, ...
    'GazeStdX', {}, 'GazeStdY', {}, 'PupilSize', {}, 'MSRate', {}, ...
    'Blinks', {}, 'Fixations', {}, 'Saccades', {});

%% Load all eye movements
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    load([datapath, filesep, 'dataET'])

    for conds = {'lc', 'hc'}
        if strcmp(conds, 'lc')
            dataET = dataET_lc;
        elseif strcmp(conds, 'hc')
            dataET = dataET_hc;
        end

        %% Initialize arrays
        subject_id = [];
        trial_num = [];
        num_trials = length(dataET.trialinfo);
        condition = [];
        gazeDev = [];
        gazeSDx = [];
        gazeSDy = [];
        pupilSize = [];
        microsaccadeRate = [];

        %% Get trial-by-trial gaze data
        for trl = 1:length(dataET.trialinfo)
            close all
            data = dataET.trial{trl};

            %% Choose data 300ms after stimulus presentation to exclude evoked activity
            analysis_period = [0.3 2];
            time_vector = dataET.time{trl};
            analysis_idx = (time_vector >= analysis_period(1)) & (time_vector <= analysis_period(2));
            data = data(:, analysis_idx);

            %% Filter out data points outside the screen boundaries
            valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
            valid_data = data(1:3, valid_data_indices);
            data = valid_data;

            %% Remove blinks with a window of 100ms (= 50 samples/timepoints)
            win_size = 50;
            data = remove_blinks(data, win_size);

            %% Extract gaze data and pupil size
            gaze_x{subj, trl} = data(1, :);
            gaze_y{subj, trl} = data(2, :);
            pupil_size{subj, trl} = mean(data(3, :), 'omitnan') / 1000;
            pups = pupil_size{subj, trl};

            % Save raw gaze data
            if strcmp(conds, 'lc')
                gaze_x_lc{subj, trl} = gaze_x{subj, trl};
                gaze_y_lc{subj, trl} = gaze_y{subj, trl};
            elseif strcmp(conds, 'hc')
                gaze_x_hc{subj, trl} = gaze_x{subj, trl};
                gaze_y_hc{subj, trl} = gaze_y{subj, trl};
            end

            %% Compute gaze deviation as euclidean distances from the center
            x_coords = gaze_x{subj, trl};
            y_coords = gaze_y{subj, trl};
            gaze_euclidean_dev = zeros(1, length(x_coords) - 1);
            % Calculate Euclidean distances
            for samples = 1:length(x_coords)
                dx = x_coords(samples) - 400; % Distance from middle of x-axis (total 800 px)
                dy = y_coords(samples) - 300; % Distance from middle of y-axis (total 600 px)
                gaze_euclidean_dev(samples) = sqrt(dx^2 + dy^2);
            end
            % Calculate the mean Euclidean distance
            mean_euclidean_distance = mean(gaze_euclidean_dev, 'omitnan');
            gaze_standard_deviation_x = nanstd(x_coords);
            gaze_standard_deviation_y = nanstd(y_coords);

            %% Compute microsaccades
            fsample = 500; % Sample rate of 500 Hz
            velData = [gaze_x{subj, trl}; gaze_y{subj, trl}]; % Concatenate x and y gaze coordinates to compute the velocity of eye movements in a 2D space
            trlLength = length(dataET.time{trl});
            [microsaccade_rate, microsaccade_details] = detect_microsaccades(fsample, velData, trlLength);
            if strcmp(conds, 'lc')
                ms_data_lc{trl} = microsaccade_details; % Save trial by trial ms data
            elseif strcmp(conds, 'hc')
                ms_data_hc{trl} = microsaccade_details; % Save trial by trial ms data
            end

            %% Append data for this trial
            subject_id = [subject_id; str2num(subjects{subj})];
            trial_num = [trial_num; trl];
            condition = [condition; dataET.trialinfo(trl)-60];
            gazeDev = [gazeDev; mean_euclidean_distance];
            gazeSDx = [gazeSDx; gaze_standard_deviation_x];
            gazeSDy = [gazeSDy; gaze_standard_deviation_y];
            pupilSize = [pupilSize; pups];
            microsaccadeRate = [microsaccadeRate; microsaccade_rate];
        end

        %% Create a trial-by-trial structure array for this subject
        subj_data_gaze_trial = struct('ID', num2cell(subject_id), 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), 'GazeDeviation', num2cell(gazeDev), 'GazeStdX', num2cell(gazeSDx), 'GazeStdY', num2cell(gazeSDy), 'PupilSize', num2cell(pupilSize), 'MSRate', num2cell(microsaccadeRate));

        %% Calculate subject-specific data by condition (GazeDev, PupilSize, MSRate)
        if strcmp(conds, 'lc')
            subj_data_gaze_trial_lc = subj_data_gaze_trial;
            lc_gdev = mean([subj_data_gaze_trial_lc.GazeDeviation], 'omitnan');
            lc_gSDx = mean([subj_data_gaze_trial_lc.GazeStdX], 'omitnan');
            lc_gSDy = mean([subj_data_gaze_trial_lc.GazeStdY], 'omitnan');
            lc_pups = mean([subj_data_gaze_trial_lc.PupilSize], 'omitnan');
            lc_msrate = mean([subj_data_gaze_trial_lc.MSRate], 'omitnan');
        elseif strcmp(conds, 'hc')
            subj_data_gaze_trial_hc = subj_data_gaze_trial;
            hc_gdev = mean([subj_data_gaze_trial_hc.GazeDeviation], 'omitnan');
            hc_gSDx = mean([subj_data_gaze_trial_hc.GazeStdX], 'omitnan');
            hc_gSDy = mean([subj_data_gaze_trial_hc.GazeStdY], 'omitnan');
            hc_pups = mean([subj_data_gaze_trial_hc.PupilSize], 'omitnan');
            hc_msrate = mean([subj_data_gaze_trial_hc.MSRate], 'omitnan');
        end
    end

    %% Load gaze metrics (extracted in GCP_preprocessing.m)
    load([datapath, filesep, 'gaze_metrics'])

    %% Create across condition structure
    subj_data_gaze = struct('ID', num2cell(subject_id(1:2)), ...
        'Condition', num2cell([1; 2]), ...
        'GazeDeviation', num2cell([lc_gdev; hc_gdev]), ...
        'GazeStdX', num2cell([lc_gSDx; hc_gSDx]), ...
        'GazeStdY', num2cell([lc_gSDy; hc_gSDy]), ...
        'PupilSize', num2cell([lc_pups; hc_pups]), ...
        'MSRate', num2cell([lc_msrate; hc_msrate]), ...
        'Blinks', num2cell([lc_blinks; hc_blinks]),...
        'Fixations', num2cell([lc_fixations; hc_fixations]),...
        'Saccades', num2cell([lc_saccades; hc_saccades]));

    %% Save data
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/gaze/');
    if isempty(dir(savepath))
        mkdir(savepath)
    end
    cd(savepath)
    save gaze_matrix_trial subj_data_gaze_trial_lc subj_data_gaze_trial_hc
    save gaze_matrix_subj subj_data_gaze
    save gaze_dev lc_gdev hc_gdev
    save gaze_std lc_gSDx lc_gSDy hc_gSDx hc_gSDy
    save pupil_size lc_pups hc_pups
    save ms_rate lc_msrate hc_msrate
    save ms_data ms_data_lc ms_data_hc
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

    % Append to the final structure array
    gaze_data = [gaze_data; subj_data_gaze];
end
save /Volumes/methlab/Students/Arne/GCP/data/features/gaze_raw gaze_x_lc gaze_y_lc gaze_x_hc gaze_y_hc
save /Volumes/methlab/Students/Arne/GCP/data/features/gaze_matrix gaze_data 
