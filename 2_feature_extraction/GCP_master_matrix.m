%% GCP MASTER Matrix

%% Setup
clear
clc
close all
path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load data
% demog_data = readtable('/Volumes/g_psyplafor_methlab$/VP/OCC/GCP/GCP_VPs.xlsx');
% demog_data = demog_data(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
% demog_data = table2struct(demog_data(1:120, :));

load('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/behavioral_matrix.mat');
load('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/eeg_matrix.mat');
load('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/gaze_matrix.mat');
ged_file = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_eeg_GED.mat';

%% Merge structures
% Convert structs to tables
T_behav = struct2table(behav_data);
T_eeg   = struct2table(eeg_data);
T_gaze  = struct2table(gaze_data);

% Keep only the key variables for merging
keys = {'ID','Condition'};

% Merge step by step
T_merge = innerjoin(T_behav, T_eeg,  'Keys', keys);
T_merge = innerjoin(T_merge, T_gaze, 'Keys', keys);

% Append GED gamma metrics without replacing legacy EEG values.
if isfile(ged_file)
    ged = load(ged_file);

    if isfield(ged, 'subjects') && ...
            isfield(ged, 'all_trial_gamma_power') && ...
            isfield(ged, 'all_trial_median_single')

        n_subj_ged = numel(ged.subjects);
        n_rows_ged = n_subj_ged * 4;

        ID_ged = nan(n_rows_ged, 1);
        Condition_ged = nan(n_rows_ged, 1);
        Power_GED = nan(n_rows_ged, 1);
        Frequency_GED = nan(n_rows_ged, 1);

        row_idx = 1;
        for subj = 1:n_subj_ged
            subj_id = str2double(ged.subjects{subj});
            for cond = 1:4
                ID_ged(row_idx) = subj_id;
                Condition_ged(row_idx) = cond;
                Power_GED(row_idx) = ged.all_trial_gamma_power(cond, subj);
                Frequency_GED(row_idx) = ged.all_trial_median_single(cond, subj);
                row_idx = row_idx + 1;
            end
        end

        T_ged = table(ID_ged, Condition_ged, Power_GED, Frequency_GED, ...
            'VariableNames', {'ID', 'Condition', 'Power_GED', 'Frequency_GED'});

        T_merge = leftjoin(T_merge, T_ged, 'Keys', keys);
    else
        warning('GED file found but required fields are missing. GED columns were not added.');
    end
else
    warning('GED file not found: %s. GED columns were not added.', ged_file);
end

% Convert back to struct array
merged_data = table2struct(T_merge);

%% Save as .mat
save /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/merged_data.mat merged_data

%% Save as .csv
merged_table = struct2table(merged_data)
csv_filename = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/merged_data.csv';
writetable(merged_table, csv_filename);