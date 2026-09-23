%% GCP Master Matrix (Subject-Level)
%
% Builds subject-level merged table across behavioral, gaze, GED metrics,
% GED inclusion, and VP demographics (Gender, Age, Handedness).
%
% Output:
%   /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_merged_data.mat
%   /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_merged_data.csv

%% Setup
[subjects, paths] = setup('GCP', 0);
features_root = paths.features;
fprintf('[MASTER MATRIX] Building subject-level merged matrix for %d subjects.\n', numel(subjects));

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
inc = load(fullfile(paths.controls, 'GCP_subject_inclusion.mat'), 'subject_inclusion');
inc = inc.subject_inclusion;
inc_ids = to_numeric_col(inc.SubjID);
merge_ids = to_numeric_col(tbl_merge.ID);
[is_match, loc] = ismember(merge_ids, inc_ids);
tbl_merge.Include = nan(height(tbl_merge), 1);
tbl_merge.Include(is_match) = double(inc.Include(loc(is_match)));

%% Demographics from VP table
demog = readtable(paths.vp_table);
demog = gcp_vp_demographics(demog);
demoIDs = to_numeric_col(demog.ID);
[has_vp, vp_loc] = ismember(merge_ids, demoIDs);
tbl_merge.Gender = repmat({''}, height(tbl_merge), 1);
tbl_merge.Age = nan(height(tbl_merge), 1);
tbl_merge.Handedness = repmat({''}, height(tbl_merge), 1);
if ~isempty(demog)
    tbl_merge.Gender(has_vp) = demog.Gender(vp_loc(has_vp));
    tbl_merge.Age(has_vp) = demog.Age(vp_loc(has_vp));
    tbl_merge.Handedness(has_vp) = demog.Handedness(vp_loc(has_vp));
    n_missing_vp = sum(~has_vp);
    if n_missing_vp > 0
        warning('GCP_master_matrix:MissingVP', ...
            '%d merged rows have no VP match in %s.', n_missing_vp, paths.vp_table);
    end
end

front = {'ID', 'Condition', 'Gender', 'Age', 'Handedness'};
front = front(ismember(front, tbl_merge.Properties.VariableNames));
rest = setdiff(tbl_merge.Properties.VariableNames, front, 'stable');
tbl_merge = tbl_merge(:, [front, rest]);

%% Diagnostics
n_behav = height(tbl_behav);
n_gaze = height(tbl_gaze);
n_ged = height(tbl_ged);
n_merge = height(tbl_merge);

fprintf('[MASTER MATRIX] Rows loaded: behavior=%d, gaze=%d, ged=%d\n', n_behav, n_gaze, n_ged);
fprintf('[MASTER MATRIX] Rows merged: %d\n', n_merge);

%% Save outputs
merged_table = tbl_merge;
merged_data = table2struct(tbl_merge);
GCP_merged_table = tbl_merge;
GCP_merged_data = merged_data;

save(fullfile(features_root, 'GCP_merged_data.mat'), ...
    'merged_data', 'merged_table', 'GCP_merged_data', 'GCP_merged_table');
writetable(tbl_merge, fullfile(features_root, 'GCP_merged_data.csv'));

fprintf('[MASTER MATRIX] Saved:\n');
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

    keep = {'ID','Condition', ...
            'BCEA','BCEA_Direction','BCEA_Eccentricity', ...
            'PupilSize','MSRate', ...
            'VelH','VelV','Vel2D', ...
            'Blinks','Fixations','Saccades', ...
            'BCEA_bl','PupilSize_bl','MSRate_bl', ...
            'VelH_bl','VelV_bl','Vel2D_bl', ...
            'Blinks_bl','Fixations_bl','Saccades_bl'};
    keep = keep(ismember(keep, G.Properties.VariableNames));
    G = G(:, keep);

    tbl_gaze = [tbl_gaze; G]; %#ok<AGROW>
end
end

function tbl_ged = load_subject_level_ged_table(features_root, subjects)
tbl_ged = table();

% Subject-level GED gamma metrics come from peaks of condition-averaged power spectra
ged_path = fullfile(features_root, 'GCP_eeg_GED.mat');
if ~isfile(ged_path)
    warning('GCP_master_matrix:NoGED', ...
        'Current GED file not found: %s. GED columns will be empty.', ged_path);
    return
end

dat = load(ged_path, ...
    'all_condition_peak_freq', 'all_condition_peak_power', ...
    'subjects');

freq_mat  = pick_first_numeric_matrix(dat, {'all_condition_peak_freq'});
pow_mat   = pick_first_numeric_matrix(dat, {'all_condition_peak_power'});

if isfield(dat, 'subjects') && ~isempty(dat.subjects)
    ged_subjects = dat.subjects;
else
    ged_subjects = subjects;
end

tbl_ged = build_ged_table_from_arrays(pow_mat, freq_mat, ged_subjects);

% Keep only key GED fields expected in downstream stats
if ~isempty(tbl_ged)
    keep = {'ID','Condition', 'Power','Frequency'};
    keep = keep(ismember(keep, tbl_ged.Properties.VariableNames));
    if numel(keep) >= 2
        tbl_ged = tbl_ged(:, keep);
    end
end
end

function tbl = build_ged_table_from_arrays(pow_mat, freq_mat, subjects)
tbl = table();

mats = {pow_mat, freq_mat};
if all(cellfun(@isempty, mats))
    return
end

nCond = 4;
nSubj = numel(subjects);
for i = 1:numel(mats)
    if ~isempty(mats{i})
        nSubj = min(nSubj, size(mats{i}, 2));
    end
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
        Frequency(row) = read_cond_subj(freq_mat, c, s);
        Power(row) = read_cond_subj(pow_mat, c, s);
    end
end

tbl = table(ID, Condition, Power, Frequency);
end

function v = read_cond_subj(M, c, s)
v = NaN;
if ~isempty(M) && size(M, 1) >= c && size(M, 2) >= s
    v = M(c, s);
end
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

function demog = gcp_vp_demographics(vp)
ID = to_numeric_col(vp.ID);
keep = isfinite(ID);
vp = vp(keep, :);
ID = ID(keep);

Gender = cellstr(strtrim(string(vp_column(vp, {'Gender', 'Geschlecht'}))));
Handedness = cellstr(strtrim(string(vp_column(vp, {'Handedness', 'H_ndigkeit'}))));

if vp_has_column(vp, {'Alter', 'Age'})
    Age = to_numeric_col(vp_column(vp, {'Alter', 'Age'}));
else
    testdate = to_datetime_col(vp_column(vp, {'Datum'}));
    dob = to_datetime_col(vp_column(vp, {'Geburtsdatum'}));
    Age = nan(size(ID));
    has_dates = ~isnat(testdate) & ~isnat(dob);
    Age(has_dates) = days(testdate(has_dates) - dob(has_dates)) / 365.25;
end

demog = table(ID, Gender, Age, Handedness);
[~, uniq_idx] = unique(demog.ID, 'stable');
demog = demog(uniq_idx, :);
end

function tf = vp_has_column(vp, candidates)
tf = ~isempty(find_vp_column(vp.Properties.VariableNames, candidates));
end

function c = vp_column(vp, candidates)
name = find_vp_column(vp.Properties.VariableNames, candidates);
if isempty(name)
    c = strings(height(vp), 1);
    return
end
c = vp.(name);
end

function name = find_vp_column(names, candidates)
name = '';
norm = lower(regexprep(string(names), '[^a-zA-Z0-9]+', ''));
cand = lower(regexprep(string(candidates), '[^a-zA-Z0-9]+', ''));
for i = 1:numel(cand)
    hit = find(norm == cand(i), 1);
    if ~isempty(hit)
        name = names{hit};
        return
    end
end
end

function dt = to_datetime_col(c)
if isdatetime(c)
    dt = c;
elseif isnumeric(c)
    dt = datetime(c, 'ConvertFrom', 'excel');
else
    dt = datetime(c);
end
end
