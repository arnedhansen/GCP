%% GCP Behavioral Feature Extraction
%
% Extracted features:
%   Accuracy

%% Setup
clear
clc
close all
path = '/Volumes/methlab_data/OCC/GCP/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
behav_data = struct('ID', {}, 'Trial', {}, 'Condition', {}, 'Accuracy', {}, 'ReactionTime', {});

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
    for block = 1:4
        load(strcat(subjects{subj}, '_GCP_block', num2str(block), '.mat'))
        num_trials = length(saves.data.correct);

        % Append data for this block
        subject_id = [subject_id; repmat({saves.subjectID}, num_trials, 1)];
        trial_num = [trial_num; (trial_counter:(trial_counter + num_trials - 1))'];
        condition = [condition; saves.data.GratingSequence'];
        accuracy = [accuracy; saves.data.allCorrect'];
        reaction_time = [reaction_time; saves.data.reactionTime'];
        trial_counter = trial_counter + num_trials;
    end

    %% Create a structure array for this subject
    sub_data_behav = struct('ID', subject_id, 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), ...
        'Accuracy', num2cell(accuracy), 'ReactionTime', num2cell(reaction_time));

    %% Calculate subject-specific Acc and RT by grating orientation
    l1 = sub_data_behav([sub_data_behav.Condition] == 1);
    l1acc = sum([l1.Accuracy])/length(l1)*100;
    l1rt = mean([l1.ReactionTime], 'omitnan');
    l2 = sub_data_behav([sub_data_behav.Condition] == 2);
    l2acc = sum([l2.Accuracy])/length(l2)*100;
    l2rt = mean([l2.ReactionTime], 'omitnan');
    l3 = sub_data_behav([sub_data_behav.Condition] == 3);
    l3acc = sum([l3.Accuracy])/length(l3)*100;
    l3rt = mean([l3.ReactionTime], 'omitnan');
    l4 = sub_data_behav([sub_data_behav.Condition] == 4);
    l4acc = sum([l4.Accuracy])/length(l4)*100;
    l4rt = mean([l4.ReactionTime], 'omitnan');
    l5 = sub_data_behav([sub_data_behav.Condition] == 5);
    l5acc = sum([l5.Accuracy])/length(l5)*100;
    l5rt = mean([l5.ReactionTime], 'omitnan');
    l6 = sub_data_behav([sub_data_behav.Condition] == 6);
    l6acc = sum([l6.Accuracy])/length(l6)*100;
    l6rt = mean([l6.ReactionTime], 'omitnan');
    l7 = sub_data_behav([sub_data_behav.Condition] == 7);
    l7acc = sum([l7.Accuracy])/length(l7)*100;
    l7rt = mean([l7.ReactionTime], 'omitnan');
    l8 = sub_data_behav([sub_data_behav.Condition] == 8);
    l8acc = sum([l8.Accuracy])/length(l8)*100;
    l8rt = mean([l8.ReactionTime], 'omitnan');

    %% Calculate subject-specific Acc and RT by grating contrast
    hc = sub_data_behav(ismember([sub_data_behav.Condition], [1, 2, 3, 4]));
    hc_acc = sum([hc.Accuracy])/length(hc)*100;
    hc_rt = mean([hc.ReactionTime], 'omitnan');
    lc = sub_data_behav(ismember([sub_data_behav.Condition], [5, 6, 7, 8]));
    lc_acc = sum([lc.Accuracy])/length(lc)*100;
    lc_rt = mean([lc.ReactionTime], 'omitnan');

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/behavioral/');
    mkdir(savepath)
    cd(savepath)
    save behavioral_matrix_subj sub_data_behav
    save acc hc_acc lc_acc l1acc l2acc l3acc l4acc l5acc l6acc l7acc l8acc
    save rt hc_rt lc_rt l1rt l2rt l3rt l4rt l5rt l6rt l7rt l8rt
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

    % Append to the final structure array
    behav_data = [behav_data; sub_data_behav];
end
save /Volumes/methlab/Students/Arne/GCP/data/features/behavioral_matrix behav_data