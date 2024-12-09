%% GCP MASTER Matrix Sternberg

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load data
load('/Volumes/methlab/Students/Arne/GCP/data/features/behavioral_matrix.mat');
load('/Volumes/methlab/Students/Arne/GCP/data/features/eeg_matrix.mat');
load('/Volumes/methlab/Students/Arne/GCP/data/features/gaze_matrix.mat');

%% Sort behavioral structure
conds = [behav_data.Condition];
[~, sortedIndices] = sort(conds);
behav = behav_data(sortedIndices);

%% Merge structures
% Initialize the merged structure array with the same size as the original structures
merged_data = struct('ID', {behav.ID}, ...
                               'Condition', {behav.Condition}, ...
                               'Accuracy', {behav.Accuracy}, ...
                               'ReactionTime', {behav.ReactionTime}, ...
                               'GazeDeviation', {gaze_data.GazeDeviation}, ...
                               'PupilSize', {gaze_data.PupilSize}, ...
                               'AlphaPower', {eeg_data.AlphaPower}, ...
                               'IAF', {eeg_data.IAF});

%% Save as .mat
save /Volumes/methlab/Students/Arne/GCP/data/features/merged_data.mat merged_data

%% Save as .csv
merged_table = struct2table(merged_data);
csv_filename = '/Volumes/methlab/Students/Arne/GCP/data/features/merged_data.csv';
writetable(merged_table, csv_filename);