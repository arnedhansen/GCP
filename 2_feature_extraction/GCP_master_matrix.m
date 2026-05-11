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
tbl_behav = table();
for s = 1:numel(subjects)
    subj = subjects{s};
    fpath = fullfile(features_root, subj, 'behavioral', 'behavioral_matrix_subj.mat');
    if ~isfile(fpath)
        continue
    end

    dat = load(fpath);
    if ~isfield(dat, 'subj_data_behav')
        continue
    end

    B = struct2table(dat.subj_data_behav);
    B = standardize_id_condition(B);

    keep = {'ID','Condition','Accuracy','ReactionTime'};
    keep = keep(ismember(keep, B.Properties.VariableNames));
    B = B(:, keep);

    tbl_behav = [tbl_behav; B]; %#ok<AGROW>
end

%% Load subject-level gaze data
tbl_gaze = load_subject_level_gaze_table(features_root, subjects);

%% Load subject-level GED data
tbl_ged = load_subject_level_ged_table(features_root, subjects);

%% Merge (outer over behavior and gaze; GED joined afterwards)
keys = {'ID','Condition'};

if isempty(tbl_behav) && isempty(tbl_gaze)
    error('No subject-level behavioral or gaze data found.');
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

tbl_merge = sortrows(tbl_merge, {'ID','Condition'});

%% Diagnostics
n_behav = height(tbl_behav);
n_gaze = height(tbl_gaze);
n_ged = height(tbl_ged);
n_merge = height(tbl_merge);

fprintf('Rows loaded: behavior=%d, gaze=%d, ged=%d\n', n_behav, n_gaze, n_ged);
fprintf('Rows merged: %d\n', n_merge);

%% Save outputs
merged_table = tbl_merge; %#ok<NASGU>
merged_data = table2struct(tbl_merge); %#ok<NASGU>
GCP_merged_table = tbl_merge; %#ok<NASGU>
GCP_merged_data = merged_data; %#ok<NASGU>

save(fullfile(features_root, 'GCP_merged_data.mat'), ...
    'merged_data', 'merged_table', 'GCP_merged_data', 'GCP_merged_table');
writetable(tbl_merge, fullfile(features_root, 'GCP_merged_data.csv'));

fprintf('Saved:\n');
fprintf('  %s\n', fullfile(features_root, 'GCP_merged_data.mat'));
fprintf('  %s\n', fullfile(features_root, 'GCP_merged_data.csv'));

%% Local helper functions
function tbl = standardize_id_condition(tbl)
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

function tbl_gaze = load_subject_level_gaze_table(features_root, subjects)
tbl_gaze = table();

for s = 1:numel(subjects)
    subj = subjects{s};
    fpath = fullfile(features_root, subj, 'gaze', 'gaze_matrix_subj.mat');
    if ~isfile(fpath)
        continue
    end

    dat = load(fpath);
    if ~isfield(dat, 'subj_data_gaze')
        continue
    end
    G = struct2table(dat.subj_data_gaze);
    G = standardize_id_condition(G);
    tbl_gaze = [tbl_gaze; G]; %#ok<AGROW>
end
end

function tbl_ged = load_subject_level_ged_table(features_root, subjects)
tbl_ged = table();

% Preferred GED summary file
ged_summary_path = fullfile(features_root, 'GCP_eeg_GED_gamma_metrics.mat');
if isfile(ged_summary_path)
    dat = load(ged_summary_path);

    if isfield(dat, 'T') && istable(dat.T)
        tbl_ged = standardize_id_condition(dat.T);
    elseif isfield(dat, 'GCP_eeg_GED_gamma_metrics')
        G = dat.GCP_eeg_GED_gamma_metrics;
        if istable(G)
            tbl_ged = standardize_id_condition(G);
        elseif isstruct(G)
            tbl_ged = standardize_id_condition(struct2table(G));
        end
    elseif isfield(dat, 'eeg_data') && isstruct(dat.eeg_data)
        tbl_ged = standardize_id_condition(struct2table(dat.eeg_data));
    else
        tbl_ged = build_ged_table_from_arrays(dat, subjects);
    end
end

% Fallback to legacy EEG matrix if GED summary is unavailable
if isempty(tbl_ged)
    legacy_path = fullfile(features_root, 'GCP_eeg_matrix.mat');
    if isfile(legacy_path)
        dat = load(legacy_path);
        if isfield(dat, 'eeg_data') && isstruct(dat.eeg_data)
            tbl_ged = standardize_id_condition(struct2table(dat.eeg_data));
        end
    end
end

% Keep only key GED fields expected in downstream stats
if ~isempty(tbl_ged)
    keep = {'ID','Condition','Power','Frequency'};
    keep = keep(ismember(keep, tbl_ged.Properties.VariableNames));
    if numel(keep) >= 2
        tbl_ged = tbl_ged(:, keep);
    end
end
end

function tbl = build_ged_table_from_arrays(dat, subjects)
tbl = table();

freq_names = {'all_trial_median_single', 'all_trial_mean_single'};
pow_names = {'all_trial_gamma_power_plotstat', 'all_trial_gamma_power'};

freq_mat = pick_first_numeric_matrix(dat, freq_names);
pow_mat = pick_first_numeric_matrix(dat, pow_names);

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

tbl = table(ID, Condition, Power, Frequency);
end

function M = pick_first_numeric_matrix(dat, candidate_names)
M = [];
for i = 1:numel(candidate_names)
    name = candidate_names{i};
    if isfield(dat, name)
        v = dat.(name);
        if isnumeric(v) && ~isempty(v)
            M = double(v);
            return
        end
    end
end
end
