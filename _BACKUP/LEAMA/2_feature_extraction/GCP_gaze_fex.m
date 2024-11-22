%% GCP Gaze Feature Extraction Sternberg
%
% Extracted features:
%   Gaze deviation
%   Pupil size

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
gaze_data_hc = struct('ID', {}, 'Trial', {}, 'Condition', {}, 'GazeDeviation', {}, 'PupilSize', {});

%% Load all eye movements HIGH CONTRAST
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    load([datapath, filesep, 'dataEThc.mat'])

    %% Initialize arrays
    subject_id = [];
    trial_num = [];
    num_trials = length(dataethc.trialinfo);
    condition = [];
    gazeDev = [];
    pupilSize = [];

    %% Get trial-by-trial gaze data
    for trl = 1:num_trials
        close all
        data = dataethc.trial{trl};

        %% Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(1:3, valid_data_indices); % Excluding pupil size data

        %% Remove blinks
        data = remove_blink_window(data, 50);

        %% Extract gaze data and pupil size
        gaze_x_hc{subj, trl} = data(1, :);
        gaze_y_hc{subj, trl} = data(2, :);
        pupil_size = data(3, :);

        %% Compute gaze deviation with euclidean distances
        x_coords = gaze_x_hc{subj, trl};
        y_coords = gaze_y_hc{subj, trl};
        num_points = length(x_coords);
        gaze_euclidean_dev = zeros(1, num_points - 1);
        % Calculate the Euclidean distances between successive points
        for i = 1:num_points - 1
            dx = x_coords(i + 1) - x_coords(i);
            dy = y_coords(i + 1) - y_coords(i);
            gaze_euclidean_dev(i) = sqrt(dx^2 + dy^2);
        end
        % Calculate the mean Euclidean distance
        mean_euclidean_distance = mean(gaze_euclidean_dev);

        % Sanity check
        % plot(gaze_euclidean_dev)

        %% Append data for this trial
        subject_id = [subject_id; str2num(subjects{subj})];
        trial_num = [trial_num; trl];
        condition = [condition; dataethc.trialinfo(trl)-50];
        gazeDev = [gazeDev; mean_euclidean_distance];
        pupilSize = [pupilSize; mean(pupil_size)/1000];
    end
    %% Create a structure array for this subject
    subj_data_gaze = struct('ID', num2cell(subject_id), 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), 'GazeDeviation', num2cell(gazeDev), 'PupilSize', num2cell(pupilSize));

    %% Calculate subject-specific GazeDev by condition
    l1 = subj_data_gaze([subj_data_gaze.Condition] == 1);
    l1gdev_hc = mean([l1.GazeDeviation], 'omitnan');
    l2 = subj_data_gaze([subj_data_gaze.Condition] == 2);
    l2gdev_hc = mean([l2.GazeDeviation], 'omitnan');
    l3 = subj_data_gaze([subj_data_gaze.Condition] == 3);
    l3gdev_hc = mean([l3.GazeDeviation], 'omitnan');
    l4 = subj_data_gaze([subj_data_gaze.Condition] == 4);
    l4gdev_hc = mean([l4.GazeDeviation], 'omitnan');
    l5 = subj_data_gaze([subj_data_gaze.Condition] == 5);
    l5gdev_hc = mean([l5.GazeDeviation], 'omitnan');
    l6 = subj_data_gaze([subj_data_gaze.Condition] == 6);
    l6gdev_hc = mean([l6.GazeDeviation], 'omitnan');
    l7 = subj_data_gaze([subj_data_gaze.Condition] == 7);
    l7gdev_hc = mean([l7.GazeDeviation], 'omitnan');
    l8 = subj_data_gaze([subj_data_gaze.Condition] == 8);
    l8gdev_hc = mean([l8.GazeDeviation], 'omitnan');

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/gaze/');
    mkdir(savepath)
    cd(savepath)
    subj_data_gaze_hc = subj_data_gaze;
    save gaze_matrix_subj_hc subj_data_gaze_hc
    save gaze_dev_hc l1gdev_hc l2gdev_hc l3gdev_hc l4gdev_hc l5gdev_hc l6gdev_hc l7gdev_hc l8gdev_hc
    clc
    disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done.'])

    % Append to the final structure array
    gaze_data_hc = [gaze_data_hc; subj_data_gaze_hc];
end
save /Volumes/methlab/Students/Arne/GCP/data/features/gaze_hc gaze_x_hc gaze_y_hc
save /Volumes/methlab/Students/Arne/GCP/data/features/gaze_matrix_hc gaze_data_hc

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
gaze_data_lc = struct('ID', {}, 'Trial', {}, 'Condition', {}, 'GazeDeviation', {}, 'PupilSize', {});

%% Load all eye movements LOW CONTRAST
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    load([datapath, filesep, 'dataETlc'])

    %% Initialize arrays
    subject_id = [];
    trial_num = [];
    num_trials = length(dataetlc.trialinfo);
    condition = [];
    gazeDev = [];
    pupilSize = [];

    %% Get trial-by-trial gaze data
    for trl = 1:length(dataetlc.trialinfo)
        close all
        data = dataetlc.trial{trl};

        %% Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(1:3, valid_data_indices); % Excluding pupil size data

        %% Remove blinks
        data = remove_blink_window(data, 50);

        %% Extract gaze data and pupil size
        gaze_x_lc{subj, trl} = data(1, :);
        gaze_y_lc{subj, trl} = data(2, :);
        pupil_size = data(3, :);

        %% Compute gaze deviation with euclidean distances
        x_coords = gaze_x_lc{subj, trl};
        y_coords = gaze_y_lc{subj, trl};
        num_points = length(x_coords);
        gaze_euclidean_dev = zeros(1, num_points - 1);
        % Calculate the Euclidean distances between successive points
        for i = 1:num_points - 1
            dx = x_coords(i + 1) - x_coords(i);
            dy = y_coords(i + 1) - y_coords(i);
            gaze_euclidean_dev(i) = sqrt(dx^2 + dy^2);
        end
        % Calculate the mean Euclidean distance
        mean_euclidean_distance = mean(gaze_euclidean_dev);

        % Sanity check
        % plot(gaze_euclidean_dev)

        %% Append data for this trial
        subject_id = [subject_id; str2num(subjects{subj})];
        trial_num = [trial_num; trl];
        condition = [condition; dataetlc.trialinfo(trl)-90];
        gazeDev = [gazeDev; mean_euclidean_distance];
        pupilSize = [pupilSize; mean(pupil_size)/1000];
    end
    %% Create a structure array for this subject
    subj_data_gaze = struct('ID', num2cell(subject_id), 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), 'GazeDeviation', num2cell(gazeDev), 'PupilSize', num2cell(pupilSize));

    %% Calculate subject-specific GazeDev by condition
    l1 = subj_data_gaze([subj_data_gaze.Condition] == 1);
    l1gdev_lc = mean([l1.GazeDeviation], 'omitnan');
    l2 = subj_data_gaze([subj_data_gaze.Condition] == 2);
    l2gdev_lc = mean([l2.GazeDeviation], 'omitnan');
    l3 = subj_data_gaze([subj_data_gaze.Condition] == 3);
    l3gdev_lc = mean([l3.GazeDeviation], 'omitnan');
    l4 = subj_data_gaze([subj_data_gaze.Condition] == 4);
    l4gdev_lc = mean([l4.GazeDeviation], 'omitnan');
    l5 = subj_data_gaze([subj_data_gaze.Condition] == 5);
    l5gdev_lc = mean([l5.GazeDeviation], 'omitnan');
    l6 = subj_data_gaze([subj_data_gaze.Condition] == 6);
    l6gdev_lc = mean([l6.GazeDeviation], 'omitnan');
    l7 = subj_data_gaze([subj_data_gaze.Condition] == 7);
    l7gdev_lc = mean([l7.GazeDeviation], 'omitnan');
    l8 = subj_data_gaze([subj_data_gaze.Condition] == 8);
    l8gdev_lc = mean([l8.GazeDeviation], 'omitnan');

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/gaze/');
    mkdir(savepath)
    cd(savepath)
    subj_data_gaze_lc = subj_data_gaze;
    save gaze_matrix_subj_lc subj_data_gaze_lc
    save gaze_dev_lc l1gdev_lc l2gdev_lc l3gdev_lc l4gdev_lc l5gdev_lc l6gdev_lc l7gdev_lc l8gdev_lc
    clc
    disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done.'])

    % Append to the final structure array
    gaze_data_lc = [gaze_data_lc; subj_data_gaze_lc];
end
save /Volumes/methlab/Students/Arne/GCP/data/features/gaze_lc gaze_x_lc gaze_y_lc
save /Volumes/methlab/Students/Arne/GCP/data/features/gaze_matrix_lc gaze_data_lc

%% Define function for blink removal
function cleaned_data = remove_blink_window(data, window_size)
blink_indices = find(all(data(1:2, :) == 0, 1));
removal_indices = [];
for i = 1:length(blink_indices)
    start_idx = max(1, blink_indices(i) - window_size);
    end_idx = min(size(data, 2), blink_indices(i) + window_size);
    removal_indices = [removal_indices, start_idx:end_idx];
end
data(:, removal_indices) = [];
cleaned_data = data;
end
