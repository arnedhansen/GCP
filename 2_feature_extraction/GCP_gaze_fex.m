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

    for conds = {'c25', 'c50', 'c75', 'c100'}
        if strcmp(conds, 'c25')
            dataET = dataET_c25;
        elseif strcmp(conds, 'c50')
            dataET = dataET_c50;
        elseif strcmp(conds, 'c75')
            dataET = dataET_c75;
        elseif strcmp(conds, 'c100')
            dataET = dataET_c100;
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
            data(2, :) = 600 - data(2, :); % Invert Y-axis

            %% Remove blinks with a window of 100ms (= 50 samples/timepoints)
            win_size = 50;
            data = remove_blinks(data, win_size);

            %% Extract gaze data and pupil size
            gaze_x{subj, trl} = data(1, :);
            gaze_y{subj, trl} = data(2, :);
            pupil_size{subj, trl} = mean(data(3, :), 'omitnan') / 1000;
            pups = pupil_size{subj, trl};

            % Save raw gaze data
            if strcmp(conds, 'c25')
                gaze_x_c25{subj, trl} = gaze_x{subj, trl};
                gaze_y_c25{subj, trl} = gaze_y{subj, trl};
            elseif strcmp(conds, 'c50')
                gaze_x_c50{subj, trl} = gaze_x{subj, trl};
                gaze_y_c50{subj, trl} = gaze_y{subj, trl};
            elseif strcmp(conds, 'c75')
                gaze_x_c75{subj, trl} = gaze_x{subj, trl};
                gaze_y_c75{subj, trl} = gaze_y{subj, trl};
            elseif strcmp(conds, 'c100')
                gaze_x_c100{subj, trl} = gaze_x{subj, trl};
                gaze_y_c100{subj, trl} = gaze_y{subj, trl};
            end

            %% Compute gaze deviation as euclidean distances from the center and gaze standard deviation
            x_coords = gaze_x{subj, trl};
            y_coords = gaze_y{subj, trl};

            % Calculate Euclidean distances
            dx = x_coords - 400 - nanmean(x_coords - 400); % Distance from mean (to get rid of impact of gaze shifts)
            dy = y_coords - 300 - nanmean(y_coords - 300); % Distance from mean (to get rid of impact of gaze shifts)
            gaze_euclidean_dev = sqrt(dx.^2 + dy.^2);

            % Calculate the mean Euclidean distance
            mean_euclidean_distance = nanmean(gaze_euclidean_dev);
            gaze_standard_deviation_x = nanstd(x_coords);
            gaze_standard_deviation_y = nanstd(y_coords);

            %% Compute microsaccades
            fsample = 500; % Sample rate of 500 Hz
            velData = [gaze_x{subj, trl}; gaze_y{subj, trl}]; % Concatenate x and y gaze coordinates to compute the velocity of eye movements in a 2D space
            trlLength = length(dataET.time{trl});
            [microsaccade_rate, microsaccade_details] = detect_microsaccades(fsample, velData, trlLength);
            if strcmp(conds, 'c25')
                ms_data_c25{trl} = microsaccade_details; % Save trial by trial ms data
            elseif strcmp(conds, 'c50')
                ms_data_c50{trl} = microsaccade_details; % Save trial by trial ms data
            elseif strcmp(conds, 'c75')
                ms_data_c75{trl} = microsaccade_details; % Save trial by trial ms data
            elseif strcmp(conds, 'c100')
                ms_data_c100{trl} = microsaccade_details; % Save trial by trial ms data
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

        %% Check data by visualizing raw gaze data
        close all
        % Preallocate arrays for averaged gaze data
        num_samples = 850; % Assuming each trial has ~850 samples
        mean_gaze_x = nan(num_samples, 1);
        mean_gaze_y = nan(num_samples, 1);

        % Stack trials into matrices for averaging
        gaze_x_matrix = cell2mat(cellfun(@(x) x(:), gaze_x(subj, :), 'UniformOutput', false)');
        gaze_y_matrix = cell2mat(cellfun(@(y) y(:), gaze_y(subj, :), 'UniformOutput', false)');

        % Calculate the mean over trials for each sample
        mean_gaze_x = nanmean(gaze_x_matrix, 2);
        mean_gaze_y = nanmean(gaze_y_matrix, 2);

        %% Plot the averaged gaze data
        % figure;
        % set(gcf, "Position", [200, 200, 1000, 600]);
        % plot(mean_gaze_x, mean_gaze_y, 'o');
        % hold on;
        % plot(400, 300, 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Centre point
        % title('Averaged Gaze Data Distribution Across Samples');
        % xlabel('X Coordinates');
        % ylabel('Y Coordinates');
        % xlim([0 800]);
        % ylim([0 600]);

        %% Create a trial-by-trial structure array for this subject
        subj_data_gaze_trial = struct('ID', num2cell(subject_id), 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), 'GazeDeviation', num2cell(gazeDev), 'GazeStdX', num2cell(gazeSDx), 'GazeStdY', num2cell(gazeSDy), 'PupilSize', num2cell(pupilSize), 'MSRate', num2cell(microsaccadeRate));

        %% Calculate subject-specific data by condition (GazeDev, PupilSize, MSRate)
        if strcmp(conds, 'c25')
            subj_data_gaze_trial_c25 = subj_data_gaze_trial;
            c25_gdev = mean([subj_data_gaze_trial_c25.GazeDeviation], 'omitnan');
            c25_gSDx = mean([subj_data_gaze_trial_c25.GazeStdX], 'omitnan');
            c25_gSDy = mean([subj_data_gaze_trial_c25.GazeStdY], 'omitnan');
            c25_pups = mean([subj_data_gaze_trial_c25.PupilSize], 'omitnan');
            c25_msrate = mean([subj_data_gaze_trial_c25.MSRate], 'omitnan');
        elseif strcmp(conds, 'c50')
            subj_data_gaze_trial_c50 = subj_data_gaze_trial;
            c50_gdev = mean([subj_data_gaze_trial_c50.GazeDeviation], 'omitnan');
            c50_gSDx = mean([subj_data_gaze_trial_c50.GazeStdX], 'omitnan');
            c50_gSDy = mean([subj_data_gaze_trial_c50.GazeStdY], 'omitnan');
            c50_gSDy = mean([subj_data_gaze_trial_c50.GazeStdY], 'omitnan');
            c50_pups = mean([subj_data_gaze_trial_c50.PupilSize], 'omitnan');
            c50_msrate = mean([subj_data_gaze_trial_c50.MSRate], 'omitnan');
        elseif strcmp(conds, 'c75')
            subj_data_gaze_trial_c75 = subj_data_gaze_trial;
            c75_gdev = mean([subj_data_gaze_trial_c75.GazeDeviation], 'omitnan');
            c75_gSDx = mean([subj_data_gaze_trial_c75.GazeStdX], 'omitnan');
            c75_gSDy = mean([subj_data_gaze_trial_c75.GazeStdY], 'omitnan');
            c75_gSDy = mean([subj_data_gaze_trial_c75.GazeStdY], 'omitnan');
            c75_pups = mean([subj_data_gaze_trial_c75.PupilSize], 'omitnan');
            c75_msrate = mean([subj_data_gaze_trial_c75.MSRate], 'omitnan');
        elseif strcmp(conds, 'c100')
            subj_data_gaze_trial_c100 = subj_data_gaze_trial;
            c100_gdev = mean([subj_data_gaze_trial_c100.GazeDeviation], 'omitnan');
            c100_gSDx = mean([subj_data_gaze_trial_c100.GazeStdX], 'omitnan');
            c100_gSDy = mean([subj_data_gaze_trial_c100.GazeStdY], 'omitnan');
            c100_gSDy = mean([subj_data_gaze_trial_c100.GazeStdY], 'omitnan');
            c100_pups = mean([subj_data_gaze_trial_c100.PupilSize], 'omitnan');
            c100_msrate = mean([subj_data_gaze_trial_c100.MSRate], 'omitnan');
        end
    end

    %% Load gaze metrics (extracted in GCP_preprocessing.m)
    load([datapath, filesep, 'gaze_metrics'])

    %% Create across condition structure
    subj_data_gaze = struct('ID', num2cell(subject_id(1:4)), ...
        'Condition', num2cell([1; 2; 3; 4]), ...
        'GazeDeviation', num2cell([c25_gdev; c50_gdev; c75_gdev; c100_gdev]), ...
        'GazeStdX', num2cell([c25_gSDx; c50_gSDx; c75_gSDx; c100_gSDx]), ...
        'GazeStdY', num2cell([c25_gSDy; c50_gSDy; c75_gSDy; c100_gSDy]), ...
        'PupilSize', num2cell([c25_pups; c50_pups; c75_pups; c100_pups]), ...
        'MSRate', num2cell([c25_msrate; c50_msrate; c75_msrate; c100_msrate]), ...
        'Blinks', num2cell([c25_blinks; c50_blinks; c75_blinks; c100_blinks]), ...
        'Fixations', num2cell([c25_fixations; c50_fixations; c75_fixations; c100_fixations]), ...
        'Saccades', num2cell([c25_saccades; c50_saccades; c75_saccades; c100_saccades]));

    %% Save data
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/', subjects{subj}, '/gaze/');
    if isempty(dir(savepath))
        mkdir(savepath)
    end
    cd(savepath)
    save gaze_matrix_trial subj_data_gaze_trial_c25 subj_data_gaze_trial_c50 subj_data_gaze_trial_c75 subj_data_gaze_trial_c100
    save gaze_matrix_subj subj_data_gaze
    save gaze_dev c25_gdev c50_gdev c75_gdev c100_gdev
    save gaze_std c25_gSDx c25_gSDy c50_gSDx c50_gSDy c75_gSDx c75_gSDy c100_gSDx c100_gSDy
    save pupil_size c25_pups c50_pups c75_pups c100_pups
    save ms_rate c25_msrate c50_msrate c75_msrate c100_msrate
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

    % Append to the final structure array
    gaze_data = [gaze_data; subj_data_gaze];
end
save /Volumes/methlab/Students/Arne/GCP/data/features/gaze_raw ...
    gaze_x_c25 gaze_y_c25 gaze_x_c50 gaze_y_c50 gaze_x_c75 gaze_y_c75 gaze_x_c100 gaze_y_c100
save /Volumes/methlab/Students/Arne/GCP/data/features/gaze_matrix gaze_data
