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

% Convert back to struct array
merged_data = table2struct(T_merge);

%% Save as .mat
save /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/merged_data.mat merged_data

%% Save as .csv
merged_table = struct2table(merged_data);
csv_filename = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/merged_data.csv';
writetable(merged_table, csv_filename);