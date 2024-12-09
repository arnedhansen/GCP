%% GCP Behavioral Feature Extraction
%
% Extracted features:
%   Accuracy
%   Reaction Times

%% Setup
clear
clc
close all
path = '/Volumes/methlab_data/OCC/GCP/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
behav_data = struct('ID', {}, 'Condition', {}, 'Accuracy', {}, 'ReactionTime', {});

%% Read data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)

    % Initialize subject-specific arrays
    subject_id = [];
    trial_num = [];
    condition = [];
    accuracy = [];
    reaction_time = [];

    %% Read blocks
    trial_counter = 1;
    for block = 1:2 %% %% %% 4
        load(strcat(subjects{subj}, '_GCP_block', num2str(block), '.mat'))
        num_trials = length(saves.data.correct);

        % Append data for this block
        subject_id = [subject_id; repmat({saves.subjectID}, num_trials, 1)];
        trial_num = [trial_num; (trial_counter:(trial_counter + num_trials - 1))'];
        condition = [condition; saves.data.grating'];
        accuracy = [accuracy; saves.data.correct'];
        reaction_time = [reaction_time; saves.data.reactionTime'];
        trial_counter = trial_counter + num_trials;
    end

    %% Create a trial-by-trial structure array for this subject
    subj_data_behav_trial = struct('ID', subject_id, 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), ...
        'Accuracy', num2cell(accuracy), 'ReactionTime', num2cell(reaction_time));

    %% Calculate subject-specific data by condition (GazeDev, PupilSize, MSRate)
    lc = subj_data_behav_trial(ismember([subj_data_behav_trial.Condition], 1));
    lc_acc = sum([lc.Accuracy])/length(lc)*100;
    lc_rt = mean([lc.ReactionTime], 'omitnan');
    hc = subj_data_behav_trial(ismember([subj_data_behav_trial.Condition], 2));
    hc_acc = sum([hc.Accuracy])/length(hc)*100;
    hc_rt = mean([hc.ReactionTime], 'omitnan');

    %% Create across condition structure
    subject_id = str2num(subjects{subj});
    subj_data_behav = struct('ID', subject_id, 'Condition', num2cell([1; 2]), 'Accuracy', num2cell([lc_acc; hc_acc]), 'ReactionTime', num2cell([lc_rt; hc_rt]));

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/behavioral/');
    mkdir(savepath)
    cd(savepath)
    save behavioral_matrix_trial subj_data_behav_trial
    save gaze_matrix_subj subj_data_behav
    save acc lc_acc hc_acc
    save rt lc_rt hc_rt
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' behavioral feature extraction done.'])

    % Append to the final structure array
    behav_data = [behav_data; subj_data_behav];
end
save /Volumes/methlab/Students/Arne/GCP/data/features/behavioral_matrix behav_data