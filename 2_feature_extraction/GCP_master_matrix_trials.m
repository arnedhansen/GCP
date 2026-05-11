%% GCP Master Matrix (Trial-Level)
%
% Builds trial-level merged table across behavioral, gaze, and GED metrics
%
% Output:
%   /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_merged_data_trials.mat
%   /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_merged_data_trials.csv

%% Setup
[subjects, paths] = setup('GCP', 0);
features_root = paths.features;
dirs = dir(features_root);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
fprintf('Building trial-level merged matrix for %d subjects.\n', numel(subjects));

%% Load trial-level behavioral data
tbl_behav = table();
for s = 1:numel(subjects)
    subj = subjects{s};
    fpath = fullfile(features_root, subj, 'behavioral', 'behavioral_matrix_trial.mat');
    if ~isfile(fpath)
        continue
    end
    dat = load(fpath);
    if ~isfield(dat, 'subj_data_behav_trial')
        continue
    end
    B = struct2table(dat.subj_data_behav_trial);
    B = standardize_id_condition_trial(B);

    keep = {'ID','Condition','Trial','Accuracy','ReactionTime'};
    keep = keep(ismember(keep, B.Properties.VariableNames));
    B = B(:, keep);

    tbl_behav = [tbl_behav; B]; %#ok<AGROW>
end

%% Load trial-level gaze data
tbl_gaze = table();
for s = 1:numel(subjects)
    subj = subjects{s};
    fpath = fullfile(features_root, subj, 'gaze', 'gaze_matrix_trial.mat');
    if ~isfile(fpath)
        continue
    end
    dat = load(fpath);

    cond_vars = {'subj_data_gaze_trial_c25', 'subj_data_gaze_trial_c50', ...
                 'subj_data_gaze_trial_c75', 'subj_data_gaze_trial_c100'};

    for ci = 1:numel(cond_vars)
        if ~isfield(dat, cond_vars{ci})
            continue
        end
        G = dat.(cond_vars{ci});
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

        tbl_gaze = [tbl_gaze; GT]; %#ok<AGROW>
    end
end

%% Load GED trial-level table
tbl_ged = table();
ged_mat = fullfile(features_root, 'GCP_eeg_GED_gamma_metrics_trials.mat');
dat = load(ged_mat);
tbl_ged = dat.T;
if ~isempty(tbl_ged)
    tbl_ged = standardize_id_condition_trial(tbl_ged);
end

%% Prefix non-key variables for unambiguous merged names
tbl_behav = prefix_nonkeys(tbl_behav, {'ID','Condition','Trial'}, 'Behavior_');
tbl_gaze  = prefix_nonkeys(tbl_gaze,  {'ID','Condition','Trial'}, 'Gaze_');
tbl_ged   = prefix_nonkeys(tbl_ged,   {'ID','Condition','Trial'}, 'GED_');

%% Merge (outer over behavior and gaze; GED joined afterwards)
keys = {'ID','Condition','Trial'};

if isempty(tbl_behav) && isempty(tbl_gaze)
    error('No trial-level behavioral or gaze data found.');
elseif isempty(tbl_behav)
    tbl_bg = tbl_gaze;
elseif isempty(tbl_gaze)
    tbl_bg = tbl_behav;
else
    tbl_bg = outerjoin(tbl_behav, tbl_gaze, 'Keys', keys, 'MergeKeys', true, 'Type', 'full');
end

if isempty(tbl_ged)
    tbl_merge = tbl_bg;
else
    tbl_merge = outerjoin(tbl_bg, tbl_ged, 'Keys', keys, 'MergeKeys', true, 'Type', 'full');
end

tbl_merge = sortrows(tbl_merge, {'ID','Condition','Trial'});

%% Diagnostics
n_behav = height(tbl_behav);
n_gaze = height(tbl_gaze);
n_ged = height(tbl_ged);
n_merge = height(tbl_merge);

fprintf('Rows loaded: behavior=%d, gaze=%d, ged=%d\n', n_behav, n_gaze, n_ged);
fprintf('Rows merged: %d\n', n_merge);

%% Save outputs with GCP_ prefix
GCP_merged_table_trials = tbl_merge; %#ok<NASGU>
GCP_merged_data_trials = table2struct(tbl_merge); %#ok<NASGU>
merged_table_trials = tbl_merge; %#ok<NASGU>
merged_data_trials = GCP_merged_data_trials; %#ok<NASGU>

save(fullfile(features_root, 'GCP_merged_data_trials.mat'), ...
    'GCP_merged_data_trials', 'GCP_merged_table_trials', ...
    'merged_data_trials', 'merged_table_trials');
writetable(tbl_merge, fullfile(features_root, 'GCP_merged_data_trials.csv'));

fprintf('Saved:\n');
fprintf('  %s\n', fullfile(features_root, 'GCP_merged_data_trials.mat'));
fprintf('  %s\n', fullfile(features_root, 'GCP_merged_data_trials.csv'));

%% Local helper functions
function tbl = standardize_id_condition_trial(tbl)
if ~ismember('ID', tbl.Properties.VariableNames)
    if ismember('Subject', tbl.Properties.VariableNames)
        tbl.ID = tbl.Subject;
    end
end

if ismember('ID', tbl.Properties.VariableNames)
    tbl.ID = to_numeric_col(tbl.ID);
end
if ismember('Condition', tbl.Properties.VariableNames)
    tbl.Condition = to_numeric_col(tbl.Condition);
end
if ismember('Trial', tbl.Properties.VariableNames)
    tbl.Trial = to_numeric_col(tbl.Trial);
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

function tbl = prefix_nonkeys(tbl, keys, prefix)
if isempty(tbl)
    return
end
vars = tbl.Properties.VariableNames;
for i = 1:numel(vars)
    v = vars{i};
    if ~ismember(v, keys) && ~startsWith(v, prefix)
        tbl.Properties.VariableNames{v} = [prefix v];
    end
end
end

function tbl = struct_of_vectors_to_table(strct)
fn = fieldnames(strct);
if isempty(fn)
    tbl = table();
    return
end
len = numel(strct.(fn{1}));
for i = 2:numel(fn)
    if numel(strct.(fn{i})) ~= len
        tbl = table();
        return
    end
end

tbl = table();
for i = 1:numel(fn)
    tbl.(fn{i}) = strct.(fn{i})(:);
end
end
