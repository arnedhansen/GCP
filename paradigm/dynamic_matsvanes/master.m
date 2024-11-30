%% Master script for the GCP (Gamma Contrast Perception) Study
%
% - Resting EEG
% - Gratings Training (10 trials)
% - Gratings (200 trials)
%   - 100x low contrast concentric dynamic inward
%   - 100x high contrast concentric dynamic inward

%% General settings, screens and paths

% Set up MATLAB workspace
clear all;
close all;
clc;
rootFilepath = pwd; % Retrieve the present working directory

% Define paths
PPDEV_PATH = '/home/methlab/Documents/MATLAB/ppdev-mex-master'; % For sending EEG triggers
DATA_PATH = '/home/methlab/Desktop/GCP_data'; % Folder to save data
FUNS_PATH = '/home/methlab/Desktop/GCP' ; % Folder with all functions

addpath(FUNS_PATH) % Add path to folder with functions
screenSettings % Manage screens

%% Collect ID and Age
dialogID;

%% Protect Matlab code from participant keyboard input
ListenChar(2);

%% Check for existing files and start tasks

if ~isfile([DATA_PATH, '/', num2str(subject.ID), '/', [num2str(subject.ID), '_Resting.mat']])
    restingEEG
else
    disp('RESTING EEG DATA ALREADY EXISTS');
end

if ~isfile([DATA_PATH, '/', num2str(subject.ID), '/', [num2str(subject.ID), '_training.mat']])
    TRAINING = 1;
    start = 1;
    TASK = 'GCP';
    BLOCK = 1;
    GCP_gratingsTask;
else
    disp('TRAINING BLOCK DATA ALREADY EXISTS');
end

TRAINING = 0;
TASK = 'GCP';
if isfile([DATA_PATH, '/', num2str(subject.ID), '/', [num2str(subject.ID), '_GCP_block4.mat']])
    disp('BLOCK 4 DATA ALREADY EXISTS');
    return
elseif isfile([DATA_PATH, '/', num2str(subject.ID), '/', [num2str(subject.ID), '_GCP_block3.mat']])
    disp('BLOCK 3 DATA ALREADY EXISTS');
    start = 4;
elseif isfile([DATA_PATH, '/', num2str(subject.ID), '/', [num2str(subject.ID), '_GCP_block2.mat']])
    disp('BLOCK 2 DATA ALREADY EXISTS');
    start = 3;
elseif isfile([DATA_PATH, '/', num2str(subject.ID), '/', [num2str(subject.ID), '_GCP_block1.mat']])
    disp('BLOCK 1 DATA ALREADY EXISTS');
    start = 2;
else
    start = 1;
end

for BLOCK = start:4
    disp([TASK, ' STARTING...']);
    GCP_gratingsTask; % Run the task
end

%% Allow keyboard input into Matlab code
ListenChar(0);
disp('GCP RECORDING FINISHED')