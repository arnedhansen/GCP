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
load('/Volumes/methlab/Students/Arne/GCP/data/features/gaze_matrix_bl.mat');

%% Merge structures
merged_data = struct('ID', {behav_data.ID}, ...
                     'Condition', {behav_data.Condition}, ...
                     'Accuracy', {behav_data.Accuracy}, ...
                     'ReactionTime', {behav_data.ReactionTime}, ...
                     'GazeDeviation', {gaze_data_bl.PctGazeDeviation}, ...
                     'GazeSTD', {gaze_data_bl.PctGazeStdX}, ...
                     'PupilSize', {gaze_data_bl.PctPupilSize}, ...
                     'MSRate', {gaze_data_bl.PctMSRate}, ...
                     'GammaPower', {eeg_data.Power}, ...
                     'GammaFreq', {eeg_data.Frequency});

%% Save as .mat
save /Volumes/methlab/Students/Arne/GCP/data/features/merged_data.mat merged_data

%% Save as .csv
merged_table = struct2table(merged_data);
csv_filename = '/Volumes/methlab/Students/Arne/GCP/data/features/merged_data.csv';
writetable(merged_table, csv_filename);