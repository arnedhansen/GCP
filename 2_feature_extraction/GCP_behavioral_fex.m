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
    try
        for block = 1:4
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
    end

    %% Create a trial-by-trial structure array for this subject
    subj_data_behav_trial = struct('ID', subject_id, 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), ...
        'Accuracy', num2cell(accuracy), 'ReactionTime', num2cell(reaction_time));

    %% Calculate subject-specific data by condition (GazeDev, PupilSize, MSRate)
    c25 = subj_data_behav_trial(ismember([subj_data_behav_trial.Condition], 1));
    c25_acc = sum([c25.Accuracy])/length(c25)*100;
    % valid_reactions = ~isnan([c25.ReactionTime]); % Logical array of non-NaN ReactionTime
    % correct_responses = [c25.Accuracy] == 1; % Logical array of Accuracy == 1
    % c25_acc = sum(valid_reactions & correct_responses) / sum(valid_reactions) * 100; % Percentage of valid and correct instances
    c25_rt = mean([c25.ReactionTime], 'omitnan'); 

    c50 = subj_data_behav_trial(ismember([subj_data_behav_trial.Condition], 2));
    c50_acc = sum([c50.Accuracy])/length(c50)*100;
    % valid_reactions = ~isnan([c50.ReactionTime]); % Logical array of non-NaN ReactionTime
    % correct_responses = [c50.Accuracy] == 1; % Logical array of Accuracy == 1
    % c50_acc = sum(valid_reactions & correct_responses) / sum(valid_reactions) * 100; % Percentage of valid and correct instances
    c50_rt = mean([c50.ReactionTime], 'omitnan');

    c75 = subj_data_behav_trial(ismember([subj_data_behav_trial.Condition], 3));
    c75_acc = sum([c75.Accuracy])/length(c75)*100;
    c75_rt = mean([c75.ReactionTime], 'omitnan');

    c100 = subj_data_behav_trial(ismember([subj_data_behav_trial.Condition], 4));
    c100_acc = sum([c100.Accuracy])/length(c100)*100;
    c100_rt = mean([c100.ReactionTime], 'omitnan');

    %% Create across condition structure
    subject_id = str2num(subjects{subj});
    subj_data_behav = struct('ID', subject_id, 'Condition', num2cell([1; 2; 3 ; 4]), ...
        'Accuracy', num2cell([c25_acc; c50_acc; c75_acc; c100_acc]), 'ReactionTime', num2cell([c25_rt; c50_rt; c75_rt; c100_rt]));

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/behavioral/');
    mkdir(savepath)
    cd(savepath)
    save behavioral_matrix_trial subj_data_behav_trial
    save gaze_matrix_subj subj_data_behav
    save acc c25_acc c50_acc c75_acc c100_acc
    save rt c25_rt c50_rt c75_rt c100_rt
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' behavioral feature extraction done.'])

    % Append to the final structure array
    behav_data = [behav_data; subj_data_behav];
end
save /Volumes/methlab/Students/Arne/GCP/data/features/behavioral_matrix behav_data