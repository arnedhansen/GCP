%% GCP MASTER Matrix

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
merged_data = struct('ID', {behav_data.ID}, ...
                     'Condition', {behav_data.Condition}, ...
                     'Accuracy', {behav_data.Accuracy}, ...
                     'ReactionTime', {behav_data.ReactionTime}, ...
                     'GazeDeviation', {gaze_data.GazeDeviation}, ...
                     'PupilSize', {gaze_data.PupilSize}, ...
                     'MSRate', {gaze_data.MSRate}, ...
                     'Blinks', {gaze_data.Blinks}, ...
                     'Fixations', {gaze_data.Fixations}, ...
                     'Saccades', {gaze_data.Saccades}, ...
                     'GammaPower', {eeg_data.Power}, ...
                     'GammaFreq', {eeg_data.Frequency});

%% Save as .mat
save /Volumes/methlab/Students/Arne/GCP/data/features/merged_data.mat merged_data

%% Save as .csv
merged_table = struct2table(merged_data);
csv_filename = '/Volumes/methlab/Students/Arne/GCP/data/features/merged_data.csv';
writetable(merged_table, csv_filename);