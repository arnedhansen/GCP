function  cutData(filePath)

% ANT EEG data comes in 1 file containig all tasks (e.g. Resting & Grating).
% This function cuts the data into distinct tasks, saves the tasks
% to .mat files and moves the original file to filePath/archive

% EEGlab with 'pop_loadeep_v4' required

%% Load the data
[EEGorig, command] = pop_loadeep_v4(filePath);

%% Extract path and filename
splitpath = strsplit(filePath, filesep);
subjectID = splitpath{7};
filePath = fullfile(filesep, splitpath{1:end-1});
fileName = splitpath{end};
splitpath_str = strsplit(fileName, '_');

%% Remove photodiode data and save to a file
diode = pop_select(EEGorig, 'channel', 129);
savename_pd = [subjectID, '_Photodiode.mat'];
save(fullfile(filePath, savename_pd), 'diode', '-v7.3')

% Exclude photodiode data
if EEGorig.nbchan > 128
    EEGorig = pop_select(EEGorig, 'nochannel', 129:EEGorig.nbchan);
end

%% Add Ref channel and data, and load channel location file
EEGorig.data(129, :) = 0;
EEGorig.nbchan = 129;
EEGorig.chanlocs(129).labels = 'CPz';
locspath = 'standard_1005.elc';
EEGorig = pop_chanedit(EEGorig, 'lookup', locspath);

%% Find start and end triggers of the data recording:

% Block 1
i11 = find(ismember({EEGorig.event.type}, '11'));
i71 = find(ismember({EEGorig.event.type}, '71'));

% Block 2
i12 = find(ismember({EEGorig.event.type}, '12'));
i72 = find(ismember({EEGorig.event.type}, '72'));

%% Cut

% Block 1
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i11(end)).latency, EEGorig.event(i71(end)).latency]);
    task = [subjectID, '_GCP_block1_task_EEG.mat'];
    % Save to a file
    save(fullfile(filePath, task), 'EEG', '-v7.3')
catch ME
    ME.message
    warning('BLOCK 1 is missing...')
end

% Block 2
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i12(end)).latency, EEGorig.event(i72(end)).latency]);
    task = [subjectID, '_GCP_block2_task_EEG.mat'];
    % Save to a file
    save(fullfile(filePath, task), 'EEG', '-v7.3')
catch ME
    ME.message
    warning('BLOCK 2 is missing...')
end

%% mkdir archive and move the orig files there
    source = fullfile(filePath, fileName);
    destination = fullfile(fullfile(filePath, 'archive'));
    mkdir(destination)
    movefile(source,destination) % .cnt file
    source = fullfile(filePath, [fileName(1:end-4), '.evt']);
    movefile(source,destination) % .evt file
    source = fullfile(filePath, [fileName(1:end-4), '.seg']);
    movefile(source,destination) % .seg file
end
