%% AOC MASTER Matrix Sternberg

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load data
load('/Volumes/methlab/Students/Arne/AOC/data/features/behavioral_matrix_sternberg.mat');
load('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg.mat');
load('/Volumes/methlab/Students/Arne/AOC/data/features/gaze_matrix_sternberg.mat');

%% Sort behavioral structure
conds = [behav_data_sternberg.Condition];
[~, sortedIndices] = sort(conds);
behav = behav_data_sternberg(sortedIndices);

%% Merge structures
% Initialize the merged structure array with the same size as the original structures
merged_data_sternberg = struct('ID', {behav.ID}, ...
                               'Condition', {behav.Condition}, ...
                               'Accuracy', {behav.Accuracy}, ...
                               'ReactionTime', {behav.ReactionTime}, ...
                               'GazeDeviation', {gaze_data_sternberg.GazeDeviation}, ...
                               'PupilSize', {gaze_data_sternberg.PupilSize}, ...
                               'AlphaPower', {eeg_data_sternberg.AlphaPower}, ...
                               'IAF', {eeg_data_sternberg.IAF});

%% Save as .mat
save /Volumes/methlab/Students/Arne/AOC/data/features/merged_data_sternberg.mat merged_data_sternberg

%% Save as .csv
merged_table = struct2table(merged_data_sternberg);
csv_filename = '/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_sternberg.csv';
writetable(merged_table, csv_filename);