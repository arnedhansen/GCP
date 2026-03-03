%% GCP Master Matrix (Trial-Level)
%
% Builds one trial-level merged table across behavioral, gaze, and GED
% metrics. Legacy subject-level script remains untouched:
%   - GCP_master_matrix.m
%
% Output:
%   - /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_merged_data_trials.mat
%   - /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_merged_data_trials.csv
%
% Keys:
%   ID, Condition, Trial

%% Setup
clear
clc
close all

features_root = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features';
dirs = dir(features_root);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
subjects = subjects(cellfun(@(s) ~isnan(str2double(s)), subjects));

fprintf('Building trial-level merged matrix for %d subjects.\n', numel(subjects));

%% Load trial-level behavioral data
T_behav = table();
for s = 1:numel(subjects)
    subj = subjects{s};
    fpath = fullfile(features_root, subj, 'behavioral', 'behavioral_matrix_trial.mat');
    if ~isfile(fpath)
        continue
    end
    S = load(fpath);
    if ~isfield(S, 'subj_data_behav_trial')
        continue
    end
    B = struct2table(S.subj_data_behav_trial);
    B = standardize_id_condition_trial(B);

    keep = {'ID','Condition','Trial','Accuracy','ReactionTime'};
    keep = keep(ismember(keep, B.Properties.VariableNames));
    B = B(:, keep);

    T_behav = [T_behav; B]; %#ok<AGROW>
end

%% Load trial-level gaze data
T_gaze = table();
for s = 1:numel(subjects)
    subj = subjects{s};
    fpath = fullfile(features_root, subj, 'gaze', 'gaze_matrix_trial.mat');
    if ~isfile(fpath)
        continue
    end
    S = load(fpath);

    cond_vars = {'subj_data_gaze_trial_c25', 'subj_data_gaze_trial_c50', ...
                 'subj_data_gaze_trial_c75', 'subj_data_gaze_trial_c100'};

    for ci = 1:numel(cond_vars)
        if ~isfield(S, cond_vars{ci})
            continue
        end
        G = S.(cond_vars{ci});
        if ~isstruct(G)
            continue
        end

        GT = struct_of_vectors_to_table(G);
        if isempty(GT)
            continue
        end
        GT = standardize_id_condition_trial(GT);

        keep = {'ID','Condition','Trial', ...
                'GazeDeviation','GazeStdX','GazeStdY','BCEA','PupilSize','MSRate', ...
                'VelH','VelV','Vel2D', ...
                'PctGazeDeviation','PctGazeStdX','PctGazeStdY','PctBCEA', ...
                'PctPupilSize','PctMSRate','PctVelH','PctVelV','PctVel2D'};
        keep = keep(ismember(keep, GT.Properties.VariableNames));
        GT = GT(:, keep);

        T_gaze = [T_gaze; GT]; %#ok<AGROW>
    end
end

%% Load GED trial-level table
T_ged = table();
ged_mat = fullfile(features_root, 'GCP_eeg_GED_gamma_metrics_trials.mat');
ged_csv = fullfile(features_root, 'GCP_eeg_GED_gamma_metrics_trials.csv');

if isfile(ged_mat)
    S = load(ged_mat);
    if isfield(S, 'T') && istable(S.T)
        T_ged = S.T;
    end
elseif isfile(ged_csv)
    T_ged = readtable(ged_csv);
end

if ~isempty(T_ged)
    T_ged = standardize_id_condition_trial(T_ged);
end

%% Prefix non-key variables for unambiguous merged names
T_behav = prefix_nonkeys(T_behav, {'ID','Condition','Trial'}, 'Behavior_');
T_gaze  = prefix_nonkeys(T_gaze,  {'ID','Condition','Trial'}, 'Gaze_');
T_ged   = prefix_nonkeys(T_ged,   {'ID','Condition','Trial'}, 'GED_');

%% Merge (outer over behavior and gaze; GED joined afterwards)
keys = {'ID','Condition','Trial'};

if isempty(T_behav) && isempty(T_gaze)
    error('No trial-level behavioral or gaze data found.');
elseif isempty(T_behav)
    T_bg = T_gaze;
elseif isempty(T_gaze)
    T_bg = T_behav;
else
    T_bg = outerjoin(T_behav, T_gaze, 'Keys', keys, 'MergeKeys', true, 'Type', 'full');
end

if isempty(T_ged)
    T_merge = T_bg;
else
    T_merge = outerjoin(T_bg, T_ged, 'Keys', keys, 'MergeKeys', true, 'Type', 'full');
end

T_merge = sortrows(T_merge, {'ID','Condition','Trial'});

%% Diagnostics
n_behav = height(T_behav);
n_gaze = height(T_gaze);
n_ged = height(T_ged);
n_merge = height(T_merge);

fprintf('Rows loaded: behavior=%d, gaze=%d, ged=%d\n', n_behav, n_gaze, n_ged);
fprintf('Rows merged: %d\n', n_merge);

%% Save outputs with GCP_ prefix
GCP_merged_table_trials = T_merge; %#ok<NASGU>
GCP_merged_data_trials = table2struct(T_merge); %#ok<NASGU>

save(fullfile(features_root, 'GCP_merged_data_trials.mat'), ...
    'GCP_merged_data_trials', 'GCP_merged_table_trials');
writetable(T_merge, fullfile(features_root, 'GCP_merged_data_trials.csv'));

fprintf('Saved:\n');
fprintf('  %s\n', fullfile(features_root, 'GCP_merged_data_trials.mat'));
fprintf('  %s\n', fullfile(features_root, 'GCP_merged_data_trials.csv'));

%% Local helper functions
function T = standardize_id_condition_trial(T)
if ~ismember('ID', T.Properties.VariableNames)
    if ismember('Subject', T.Properties.VariableNames)
        T.ID = T.Subject;
    end
end

if ismember('ID', T.Properties.VariableNames)
    T.ID = to_numeric_col(T.ID);
end
if ismember('Condition', T.Properties.VariableNames)
    T.Condition = to_numeric_col(T.Condition);
end
if ismember('Trial', T.Properties.VariableNames)
    T.Trial = to_numeric_col(T.Trial);
end
end

function c = to_numeric_col(c)
if isnumeric(c)
    c = double(c);
elseif iscell(c)
    if isempty(c)
        c = [];
    elseif all(cellfun(@isnumeric, c))
        c = cellfun(@double, c);
    else
        c = str2double(string(c));
    end
elseif isstring(c) || ischar(c) || iscategorical(c)
    c = str2double(string(c));
else
    c = str2double(string(c));
end
end

function T = prefix_nonkeys(T, keys, prefix)
if isempty(T)
    return
end
vars = T.Properties.VariableNames;
for i = 1:numel(vars)
    v = vars{i};
    if ~ismember(v, keys) && ~startsWith(v, prefix)
        T.Properties.VariableNames{v} = [prefix v];
    end
end
end

function T = struct_of_vectors_to_table(S)
fn = fieldnames(S);
if isempty(fn)
    T = table();
    return
end
len = numel(S.(fn{1}));
for i = 2:numel(fn)
    if numel(S.(fn{i})) ~= len
        T = table();
        return
    end
end

T = table();
for i = 1:numel(fn)
    T.(fn{i}) = S.(fn{i})(:);
end
end
