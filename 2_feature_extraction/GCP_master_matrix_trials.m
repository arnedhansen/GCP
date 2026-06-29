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
                'dBGazeDeviation','dBGazeStdX','dBGazeStdY','dBBCEA', ...
                'dBPupilSize','dBMSRate','dBVelH','dBVelV','dBVel2D'};
        keep = keep(ismember(keep, GT.Properties.VariableNames));
        GT = GT(:, keep);

        tbl_gaze = [tbl_gaze; GT]; %#ok<AGROW>
    end
end

%% Load GED trial-level table (current pipeline: GCP_eeg_GED.mat)
% Per-trial peak gamma frequency [Hz] (trials_peaks) and peak gamma power [dB]
% (mean powratio within peak +/- 5 Hz, with the same power-outlier mask used
% by GCP_eeg_fex_GED.m / GCP_stats_rainclouds.py) so the trial-level merged
% table matches the GED spectra and the subject-level boxplots.
tbl_ged = build_ged_trial_table(features_root, subjects);
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
inc = load(fullfile(paths.controls, 'GCP_subject_inclusion.mat'), 'subject_inclusion');
inc = inc.subject_inclusion;
[~, loc] = ismember(tbl_merge.ID, inc.SubjID);
tbl_merge.Include = inc.Include(loc);

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

function tbl_ged = build_ged_trial_table(features_root, subjects)
tbl_ged = table();

ged_path = fullfile(features_root, 'GCP_eeg_GED.mat');
if ~isfile(ged_path)
    warning('GCP_master_matrix_trials:NoGED', ...
        'Current GED file not found: %s. GED trial columns will be empty.', ged_path);
    return
end

dat = load(ged_path, 'trials_peaks', 'trials_powratio_fullscan', ...
    'trials_outlier_mask_power_full', 'scan_freqs', 'subjects');
if ~isfield(dat, 'trials_peaks') || ~isfield(dat, 'trials_powratio_fullscan')
    warning('GCP_master_matrix_trials:BadGED', ...
        'GED file missing trial-level fields; GED trial columns will be empty.');
    return
end

scan_freqs = dat.scan_freqs(:)';
peak_power_halfwidth_hz = 5;  % matches GCP_eeg_fex_GED.m / GCP_stats_rainclouds.py

if isfield(dat, 'subjects') && ~isempty(dat.subjects)
    ged_subjects = dat.subjects;
else
    ged_subjects = subjects;
end

nCond = size(dat.trials_peaks, 1);
nSubj = size(dat.trials_peaks, 2);

ID = [];
Condition = [];
Trial = [];
GammaFrequency = [];
GammaPower = [];

for c = 1:nCond
    for s = 1:nSubj
        pf = dat.trials_peaks{c, s};
        pr = dat.trials_powratio_fullscan{c, s};
        if isempty(pf) || isempty(pr)
            continue
        end

        pf = pf(:);
        nTrl = numel(pf);
        if size(pr, 1) ~= nTrl
            continue
        end

        pp = reconstruct_trial_peak_power(pf, pr, scan_freqs, peak_power_halfwidth_hz);

        if isfield(dat, 'trials_outlier_mask_power_full') && ...
                ~isempty(dat.trials_outlier_mask_power_full) && ...
                c <= size(dat.trials_outlier_mask_power_full, 1) && ...
                s <= size(dat.trials_outlier_mask_power_full, 2)
            mask = dat.trials_outlier_mask_power_full{c, s};
            if ~isempty(mask) && numel(mask) == nTrl
                pp(logical(mask(:))) = NaN;
            end
        end

        sid = str2double(ged_subjects{s});
        ID = [ID; repmat(sid, nTrl, 1)]; %#ok<AGROW>
        Condition = [Condition; repmat(c, nTrl, 1)]; %#ok<AGROW>
        Trial = [Trial; (1:nTrl)']; %#ok<AGROW>
        GammaFrequency = [GammaFrequency; pf]; %#ok<AGROW>
        GammaPower = [GammaPower; pp(:)]; %#ok<AGROW>
    end
end

if isempty(ID)
    return
end

tbl_ged = table(ID, Condition, Trial, GammaFrequency, GammaPower);
end

function pp = reconstruct_trial_peak_power(pf, pr, scan_freqs, halfwidth_hz)
nTrl = numel(pf);
pp = nan(nTrl, 1);
scan_freqs = scan_freqs(:)';

for t = 1:nTrl
    if ~isfinite(pf(t))
        continue
    end
    band = abs(scan_freqs - pf(t)) <= halfwidth_hz;
    if ~any(band)
        continue
    end
    pp(t) = mean(pr(t, band), 'omitnan');
end
end
