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
                EEGname = dEEG(files).name;
                fileTaskName = strsplit(EEGname, '_');
                isResting = contains(EEGname, 'Resting');

                if isResting
                    ETname = [char(subjectID) '_Resting_ET.mat'];
                    block = '';
                else
                    block = '';
                    for k = 1:numel(fileTaskName)
                        tok = char(fileTaskName{k});
                        if startsWith(tok, 'block')
                            block = tok(6:end);
                            break
                        end
                    end
                    if isempty(block)
                        error('No block token in %s', EEGname);
                    end
                    ETname = [char(subjectID) '_GCP_block' block '_task_ET.mat'];
                end

                load(fullfile(dEEG(files).folder, EEGname));
                ETfile = fullfile(filePathET, ETname);

                %% Define start and end triggers
                if isResting
                    startTrigger = 10;
                    endTrigger = 90;
                else
                    startTriggers = 11:14;
                    startTriggersCell = arrayfun(@num2str, 11:14, 'UniformOutput', 0);
                    startTrigger = startTriggers(ismember(startTriggersCell, {EEG.event.type}));
                    endTrigger = startTrigger + 60;
                end

                %% Merge files
                EEG = pop_importeyetracker(EEG, ETfile,[startTrigger endTrigger],[2 3 4],{'L_GAZE_X', 'L_GAZE_Y', 'L_AREA', 'R_GAZE_X', 'R_GAZE_Y', 'R_AREA'},1,1,1,0);

                %% Save to disk
                if isResting
                    fileName = [char(subjectID) '_EEG_ET_RestingEO_merged'];
                else
                    fileName = [char(subjectID) '_EEG_ET_GCP_block' block '_merged'];
                end
                save(fullfile(resultFolder, fileName), 'EEG', '-v7.3')
                if isResting
                    fprintf('[MERGE] Subject %s: Resting done\n', char(subjectID))
                else
                    fprintf('[MERGE] Subject %s: block %s done\n', char(subjectID), block)
                end
            end
        catch ME
            fprintf('[MERGE] ERROR processing file %d: %s\n', files, ME.message);
        end
    end
end
fprintf('[MERGE] SYNCHRONIZATION COMPLETE\n')
toc