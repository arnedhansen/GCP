%% Merge GCP ET and EEG data
% Use AFTER AUTOMAGIC

%% Setup
clear
addpath /Volumes/methlab/4marius_bdf/eeglab % for pop_importeyetracker (EYE-EEG)
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/automagic/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjectIDs = {folders.name};

%% Merge data
tic;
for subjects = 1 : length(subjectIDs)
    subjectID = subjectIDs(subjects);
    % Check if subject files have already been merged
    if isempty(dir(['/Volumes/methlab/Students/Arne/GCP/data/merged/', char(subjectID), filesep, char(subjectID), '*_merged.mat']))
        % Set up data paths
        filePathET = ['/Volumes/methlab_data/OCC/GCP/data/', char(subjectID)];
        filePathEEG = ['/Volumes/methlab/Students/Arne/GCP/data/automagic/',  char(subjectID)];
        resultFolder = ['/Volumes/methlab/Students/Arne/GCP/data/merged/', char(subjectID)];
        mkdir(resultFolder)
        dEEG = dir([filePathEEG, filesep, '*ip*EEG.mat']);
        dET = dir([filePathET, filesep, '*ET.mat']);

        for files = 1 : size(dEEG, 1)
                ETnameShort = dET(files).name(1:end-7);
                ETname = dET(files).name;

                idxEEG = contains({dEEG.name}, ETnameShort);

                EEGname = dEEG(idxEEG).name;

                load(fullfile(dEEG(idxEEG).folder, EEGname));
                ETfile = fullfile(dET(1).folder, ETname);

                fileTaskName = strsplit(EEGname, '_');
                task = sprintf('%s', char(fileTaskName(3)));
                if ~strcmp(task, 'Resting_EEG.mat')
                    block = sprintf('%s', char(fileTaskName(4)));
                    block = block(6:end);
                end

                %% Define start and end triggers
                % Resting
                if strcmp(task, 'Resting_EEG.mat')
                    startTrigger = 10;
                    endTrigger = 90;
                % Grating Task    
                else
                    startTriggers = 11:14;
                    endTriggers = 71:74;
                    startTriggersCell = arrayfun(@num2str, 11:14, 'UniformOutput', 0);

                    startTrigger = startTriggers(ismember(startTriggersCell, {EEG.event.type}));
                    endTrigger = endTriggers(ismember(startTriggersCell, {EEG.event.type}));
                end
                % End trigger
                endTrigger = startTrigger + 60;

                %% Merge files
                EEG = pop_importeyetracker(EEG, ETfile,[startTrigger endTrigger],[2 3 4],{'L_GAZE_X', 'L_GAZE_Y', 'L_AREA', 'R_GAZE_X', 'R_GAZE_Y', 'R_AREA'},1,1,1,0);

                %% Save to disk
                if strcmp(task, 'Resting_EEG.mat') == 1
                    fileName = [char(subjectID) '_EEG_ET_RestingEO_merged'];
                elseif strcmp(task, 'GCP') == 1
                    fileName = [char(subjectID) '_EEG_ET_GCP_block' num2str(block) '_merged'];
                end
                save(fullfile(resultFolder, fileName), 'EEG', '-v7.3')
                if strcmp(task, 'Resting_EEG.mat') == 1
                    disp(['GCP' char(subjectID) ': Resting done' ])
                else
                    step = sprintf('%s', char(fileTaskName(3)), '_', char(fileTaskName(4)));
                    disp(['GCP' char(subjectID) ': ' step ' done' ])
                end
        end
    end
end
disp('SYNCHRONIZATION COMPLETE')
toc