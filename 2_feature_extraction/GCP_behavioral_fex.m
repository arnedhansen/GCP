%% GCP Behavioral Feature Extraction
%
% Extracted features:
%   Accuracy (catch-hit rate only: space on white-fixation trials /
%             all white-fixation trials)
%   Reaction Times
%   WhiteCross (trial table only)

%% Setup
startup
[subjects, paths] = setup('GCP');
rawPath = paths.raw_occ;
behav_data = struct('ID', {}, 'Condition', {}, 'Accuracy', {}, 'ReactionTime', {});

%% Read data
for subj = 1:length(subjects)
    clc; fprintf('[BEHAV FEX] Subject %d / %d (%s)\n', subj, length(subjects), subjects{subj});
    datapath = fullfile(rawPath, subjects{subj});
    cd(datapath)

    % Initialize subject-specific arrays
    subject_id = [];
    trial_num = [];
    condition = [];
    accuracy = [];
    reaction_time = [];
    white_cross = [];

    %% Read blocks
    trial_counter = 1;
    for block = 1:4
        load(sprintf('%s_GCP_block%d.mat', subjects{subj}, block))
        num_trials = length(saves.data.correct);

        % Append data for this block
        subject_id = [subject_id; repmat({saves.subjectID}, num_trials, 1)];
        trial_num = [trial_num; (trial_counter:(trial_counter + num_trials - 1))'];
        condition = [condition; saves.data.grating'];
        accuracy = [accuracy; saves.data.correct(:)];
        reaction_time = [reaction_time; saves.data.reactionTime(:)];
        if ~isfield(saves.data, 'whiteCross')
            if ~isfield(saves.data, 'redCross')
                error('GCP_behavioral_fex:MissingWhiteCross', ...
                    'whiteCross missing in %s_GCP_block%d.mat', subjects{subj}, block);
            end
            saves.data.whiteCross = saves.data.redCross;
        end
        white_cross = [white_cross; saves.data.whiteCross(:)];
        trial_counter = trial_counter + num_trials;
    end

    %% Create a trial-by-trial structure array for this subject
    subj_data_behav_trial = struct('ID', subject_id, 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), ...
        'Accuracy', num2cell(accuracy), 'ReactionTime', num2cell(reaction_time), ...
        'WhiteCross', num2cell(white_cross));

    %% Catch-hit accuracy by condition (white-fixation trials only)
    c25 = subj_data_behav_trial(ismember([subj_data_behav_trial.Condition], 1));
    c25_acc = catch_hit_rate(c25);
    c25_rt = mean([c25.ReactionTime], 'omitnan');

    c50 = subj_data_behav_trial(ismember([subj_data_behav_trial.Condition], 2));
    c50_acc = catch_hit_rate(c50);
    c50_rt = mean([c50.ReactionTime], 'omitnan');

    c75 = subj_data_behav_trial(ismember([subj_data_behav_trial.Condition], 3));
    c75_acc = catch_hit_rate(c75);
    c75_rt = mean([c75.ReactionTime], 'omitnan');

    c100 = subj_data_behav_trial(ismember([subj_data_behav_trial.Condition], 4));
    c100_acc = catch_hit_rate(c100);
    c100_rt = mean([c100.ReactionTime], 'omitnan');

    %% Create across condition structure
    subject_id = str2num(subjects{subj});
    subj_data_behav = struct('ID', subject_id, 'Condition', num2cell([1; 2; 3 ; 4]), ...
        'Accuracy', num2cell([c25_acc; c50_acc; c75_acc; c100_acc]), 'ReactionTime', num2cell([c25_rt; c50_rt; c75_rt; c100_rt]));

    %% Save
    savepath = fullfile(paths.features, subjects{subj}, 'behavioral');
    mkdir(savepath)
    cd(savepath)
    save behavioral_matrix_trial subj_data_behav_trial
    save behavioral_matrix_subj subj_data_behav
    save acc c25_acc c50_acc c75_acc c100_acc
    save rt c25_rt c50_rt c75_rt c100_rt
    clc
    fprintf('[BEHAV FEX] Subject %d / %d (%s) done\n', subj, length(subjects), subjects{subj})

    % Append to the final structure array
    behav_data = [behav_data; subj_data_behav];
end
save(fullfile(paths.features, 'GCP_behavioral_matrix.mat'), 'behav_data')
fprintf('[BEHAV FEX] Done. %d/%d subjects\n', length(subjects), length(subjects))

function acc = catch_hit_rate(trials)
% Space bar on white-fixation trials / all white-fixation trials, in percent.
if isempty(trials)
    acc = NaN;
    return
end
wc = [trials.WhiteCross];
wc = wc ~= 0;
n = sum(wc);
if n == 0
    acc = NaN;
    return
end
acc = sum([trials(wc).Accuracy]) / n * 100;
end
