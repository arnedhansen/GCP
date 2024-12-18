%% GCP ANT Neuro Impedances Visualization
% This script extracts impedances before and after tasks from ANT Neuro .cnt files
% and visualises the results for each subject and block.
clear
[subjects, path, colors] = setup('GCP');
addpath('/Volumes/methlab/Utility/scriptsForEEG')

%% Define paths
baseDir = '/Volumes/methlab_data/OCC/GCP/data/'; % Base directory
saveDir = '/Volumes/methlab/Students/Arne/GCP/figures/controls'; % Directory to save figures
if ~exist(saveDir, 'dir'), mkdir(saveDir); end % Create save directory if it doesn't exist

% Find all subjects (subfolders in baseDir)
subjectDirs = dir(fullfile(baseDir, '*/archive')); % Subfolders with "archive"
cntFiles = {}; % Store all .cnt file paths

% Collect all .cnt files from archive folders
for i = 1:length(subjectDirs)
    subjPath = fullfile(subjectDirs(i).folder, subjectDirs(i).name);
    files = dir(fullfile(subjPath, '*.cnt')); % Find .cnt files in archive
    for f = 1:length(files)
        cntFiles{end+1} = fullfile(files(f).folder, files(f).name);
    end
end

%% Define task start and end triggers
% Update these triggers based on your experiment
taskStartTrig = 10; % Example task start trigger
taskEndTrig = 90;   % Example task end trigger

%% Preallocate storage
allSubjects = {}; % To store subject IDs
impData = {};     % To store impedance data

%% Loop through all .cnt files
fprintf('Starting impedance extraction and visualisation...\n');
for i = 1%:length(cntFiles)
    % Get current file path
    filePath = cntFiles{i};
    fprintf('Processing file: %s\n', filePath);
    
    % Extract subject ID from path
    parts = strsplit(filePath, filesep);
    subjectID = parts{end-2}; % Extract subject ID from the folder structure
    allSubjects{i} = subjectID;

    % Call the impedance extraction function
    try
        [imps, offsetBefore, offsetAfter, imps_bipolar] = ...
            getImpsTEST(filePath, taskStartTrig, taskEndTrig);
    catch ME
        fprintf('Error processing file %s: %s\n', filePath, ME.message);
        continue;
    end

    % Store impedance data
    impData{i}.subjectID = subjectID;
    impData{i}.filePath = filePath;
    impData{i}.imps = imps;
    impData{i}.offsetBefore = offsetBefore;
    impData{i}.offsetAfter = offsetAfter;
    impData{i}.imps_bipolar = imps_bipolar;

    %% Visualise impedances for the current subject/block
    figure('Name', ['Impedances - ' subjectID], 'NumberTitle', 'off', 'Color', 'w');
    
    % Plot impedances before task
    subplot(1, 2, 1);
    bar(imps(1, :));
    title('Impedances Before Task');
    xlabel('Electrodes');
    ylabel('Impedance (Ohms)');
    grid on;
    set(gca, 'FontSize', 12);

    % Plot impedances after task
    subplot(1, 2, 2);
    bar(imps(2, :));
    title('Impedances After Task');
    xlabel('Electrodes');
    ylabel('Impedance (Ohms)');
    grid on;
    set(gca, 'FontSize', 12);

    % Save figure
    saveFileName = fullfile(saveDir, [subjectID '_impedances.png']);
    saveas(gcf, saveFileName);
    close(gcf);
end

%% Summary Output
fprintf('Finished processing %d files.\n', length(cntFiles));


function [imps, offsetBefore, offsetAfter, imps_bipolar] = getImpsTEST(fullFileName, taskStartTrig, taskEndTrig)
% Extract impedances before and after task start/stop triggers from ANT Neuro .cnt files.

% Initialise outputs
imps = nan(2, 129);
offsetBefore = nan;
offsetAfter = nan;
imps_bipolar = nan(2, 4);

try
    % Read file information and EEG data
    r.v4_info = eepv4_read_info(fullFileName);
    r.sample1 = 1;
    r.sample2 = r.v4_info.sample_count;
    r.v4_data = eepv4_read(fullFileName, r.sample1, r.sample2);
    
    % Find triggers for task start and stop
    startTriggerIdx = find(strcmp({r.v4_data.triggers.label}, taskStartTrig));
    stopTriggerIdx = find(strcmp({r.v4_data.triggers.label}, taskEndTrig));
    
    if isempty(startTriggerIdx) || isempty(stopTriggerIdx)
        error('Task start or stop trigger not found in file: %s', fullFileName);
    end
    
    restStart = r.v4_data.triggers(startTriggerIdx).offset_in_file;
    restStop = r.v4_data.triggers(stopTriggerIdx).offset_in_file;
    
    % Find impedance events
    imp_events = r.v4_data.triggers(strcmp({r.v4_data.triggers.description}, 'Impedance'));
    if isempty(imp_events)
        error('No impedance events found in file: %s', fullFileName);
    end
    
    % Closest impedance before task start
    [offsetBefore, closestToRestStart_ind] = min(abs([imp_events.offset_in_file] - restStart));
    imps_before = str2double(strsplit(imp_events(closestToRestStart_ind).impedances, ' '));
    
    % Closest impedance after task stop
    [offsetAfter, closestToRestEnd_ind] = min(abs([imp_events.offset_in_file] - restStop));
    imps_after = str2double(strsplit(imp_events(closestToRestEnd_ind).impedances, ' '));
    
    % Check if valid impedance data is extracted
    if length(imps_before) < 129 || length(imps_after) < 129
        error('Insufficient impedance data in file: %s', fullFileName);
    end
    
    % Assign impedance values (first 129 channels)
    imps(1, :) = imps_before(1:129);
    imps(2, :) = imps_after(1:129);
    
    % Handle bipolar channels (if available)
    if length(imps_after) > 129
        imps_bipolar(1, :) = imps_before(130:end); % Skip the photodiode
        imps_bipolar(2, :) = imps_after(130:end);
    end
catch ME
    % Display error message and continue execution
    fprintf('Error processing file %s: %s\n', fullFileName, ME.message);
end

end


