%% Merge GCP ET and EEG data
% Use AFTER AUTOMAGIC

%% Setup
startup
clear
[~, paths] = setup('GCP', 0);
addpath(paths.eeglab) % for pop_importeyetracker (EYE-EEG)
eeglab
clc
close all
path = paths.automagic;
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjectIDs = {folders.name};

%% Merge data
tic;
for subjects = 1 : length(subjectIDs)
    subjectID = subjectIDs(subjects);
    % Check if subject files have already been merged
    mergedPattern = fullfile(paths.merged, char(subjectID), [char(subjectID) '*_merged.mat']);
    if isempty(dir(mergedPattern))
        % Set up data paths
        filePathET = fullfile(paths.raw_occ, char(subjectID));
        filePathEEG = fullfile(paths.automagic, char(subjectID));
        resultFolder = fullfile(paths.merged, char(subjectID));
        mkdir(resultFolder)
        dEEG = dir(fullfile(filePathEEG, '*ip*EEG.mat'));
        dET = dir(fullfile(filePathET, '*ET.mat'));
        try
            for files = 1 : size(dEEG, 1)
                ETnameShort = dET(files).name(1:end-7);
                if files == size(dEEG, 1)
                    ETnameShort = [dET(files).name(1:3), '_GCP_', dET(files).name(5:end-7)];
                end
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
        catch ME
            fprintf('ERROR processing file %d: %s\n', files, ME.message);
        end
    end
end
disp('SYNCHRONIZATION COMPLETE')
toc