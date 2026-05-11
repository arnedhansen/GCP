%% GCP Master Matrix (Subject-Level)
%
% Builds subject-level merged table across behavioral, gaze, and GED metrics
%
% Output:
%   /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_merged_data.mat
%   /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_merged_data.csv

%% Setup
[subjects, paths] = setup('GCP', 0);
features_root = paths.features;
fprintf('Building subject-level merged matrix for %d subjects.\n', numel(subjects));

%% Load subject-level behavioral data
T_behav = table();
for s = 1:numel(subjects)
    subj = subjects{s};
    fpath = fullfile(features_root, subj, 'behavioral', 'behavioral_matrix_subj.mat');
    if ~isfile(fpath)
        continue
    end

    S = load(fpath);
    if ~isfield(S, 'subj_data_behav')
        continue
    end

    B = struct2table(S.subj_data_behav);
    B = standardize_id_condition(B);

    keep = {'ID','Condition','Accuracy','ReactionTime'};
    keep = keep(ismember(keep, B.Properties.VariableNames));
    B = B(:, keep);

    T_behav = [T_behav; B]; %#ok<AGROW>
end

%% Load subject-level gaze data
T_gaze = load_subject_level_gaze_table(features_root, subjects);

%% Load subject-level GED data
T_ged = load_subject_level_ged_table(features_root, subjects);

%% Merge (outer over behavior and gaze; GED joined afterwards)
keys = {'ID','Condition'};

if isempty(T_behav) && isempty(T_gaze)
    error('No subject-level behavioral or gaze data found.');
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

T_merge = sortrows(T_merge, {'ID','Condition'});

%% Diagnostics
n_behav = height(T_behav);
n_gaze = height(T_gaze);
n_ged = height(T_ged);
n_merge = height(T_merge);

fprintf('Rows loaded: behavior=%d, gaze=%d, ged=%d\n', n_behav, n_gaze, n_ged);
fprintf('Rows merged: %d\n', n_merge);

%% Save outputs
merged_table = T_merge; %#ok<NASGU>
merged_data = table2struct(T_merge); %#ok<NASGU>
GCP_merged_table = T_merge; %#ok<NASGU>
GCP_merged_data = merged_data; %#ok<NASGU>

save(fullfile(features_root, 'GCP_merged_data.mat'), ...
    'merged_data', 'merged_table', 'GCP_merged_data', 'GCP_merged_table');
writetable(T_merge, fullfile(features_root, 'GCP_merged_data.csv'));

fprintf('Saved:\n');
fprintf('  %s\n', fullfile(features_root, 'GCP_merged_data.mat'));
fprintf('  %s\n', fullfile(features_root, 'GCP_merged_data.csv'));

%% Local helper functions
function T = standardize_id_condition(T)
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

function T_gaze = load_subject_level_gaze_table(features_root, subjects)
T_gaze = table();

for s = 1:numel(subjects)
    subj = subjects{s};
    fpath = fullfile(features_root, subj, 'gaze', 'gaze_matrix_subj.mat');
    if ~isfile(fpath)
        continue
    end

    S = load(fpath);
    if ~isfield(S, 'subj_data_gaze')
        continue
    end
    G = struct2table(S.subj_data_gaze);
    G = standardize_id_condition(G);
    T_gaze = [T_gaze; G]; %#ok<AGROW>
end
end

function T_ged = load_subject_level_ged_table(features_root, subjects)
T_ged = table();

% Preferred GED summary file
ged_summary_path = fullfile(features_root, 'GCP_eeg_GED_gamma_metrics.mat');
if isfile(ged_summary_path)
    S = load(ged_summary_path);

    if isfield(S, 'T') && istable(S.T)
        T_ged = standardize_id_condition(S.T);
    elseif isfield(S, 'GCP_eeg_GED_gamma_metrics')
        G = S.GCP_eeg_GED_gamma_metrics;
        if istable(G)
            T_ged = standardize_id_condition(G);
        elseif isstruct(G)
            T_ged = standardize_id_condition(struct2table(G));
        end
    elseif isfield(S, 'eeg_data') && isstruct(S.eeg_data)
        T_ged = standardize_id_condition(struct2table(S.eeg_data));
    else
        T_ged = build_ged_table_from_arrays(S, subjects);
    end
end

% Fallback to legacy EEG matrix if GED summary is unavailable
if isempty(T_ged)
    legacy_path = fullfile(features_root, 'GCP_eeg_matrix.mat');
    if isfile(legacy_path)
        L = load(legacy_path);
        if isfield(L, 'eeg_data') && isstruct(L.eeg_data)
            T_ged = standardize_id_condition(struct2table(L.eeg_data));
        end
    end
end

% Keep only key GED fields expected in downstream stats
if ~isempty(T_ged)
    keep = {'ID','Condition','Power','Frequency'};
    keep = keep(ismember(keep, T_ged.Properties.VariableNames));
    if numel(keep) >= 2
        T_ged = T_ged(:, keep);
    end
end
end

function T = build_ged_table_from_arrays(S, subjects)
T = table();

freq_names = {'all_trial_median_single', 'all_trial_mean_single'};
pow_names = {'all_trial_gamma_power_plotstat', 'all_trial_gamma_power'};

freq_mat = pick_first_numeric_matrix(S, freq_names);
pow_mat = pick_first_numeric_matrix(S, pow_names);

if isempty(freq_mat) && isempty(pow_mat)
    return
end

nCond = 4;
nSubj = numel(subjects);
if ~isempty(freq_mat)
    nSubj = min(nSubj, size(freq_mat, 2));
end
if ~isempty(pow_mat)
    nSubj = min(nSubj, size(pow_mat, 2));
end

ID = nan(nCond * nSubj, 1);
Condition = nan(nCond * nSubj, 1);
Frequency = nan(nCond * nSubj, 1);
Power = nan(nCond * nSubj, 1);

row = 0;
for s = 1:nSubj
    sid = str2double(subjects{s});
    for c = 1:nCond
        row = row + 1;
        ID(row) = sid;
        Condition(row) = c;
        if ~isempty(freq_mat) && size(freq_mat, 1) >= c
            Frequency(row) = freq_mat(c, s);
        end
        if ~isempty(pow_mat) && size(pow_mat, 1) >= c
            Power(row) = pow_mat(c, s);
        end
    end
end

T = table(ID, Condition, Power, Frequency);
end

function M = pick_first_numeric_matrix(S, candidate_names)
M = [];
for i = 1:numel(candidate_names)
    name = candidate_names{i};
    if isfield(S, name)
        v = S.(name);
        if isnumeric(v) && ~isempty(v)
            M = double(v);
            return
        end
    end
end
end
