%% Cutting of GCP EEG data for the GCP study
% Finds .cnt files that are still in the subject folder.
% Subjects that already have cut EEG files are left unchanged.

%% EEGlab
clear
clc
[~, paths] = setup('GCP');
p = pwd;
cd(paths.eeglab)
eeglab
close()
cd(p)
clc

%% Define path, cut data and convert asc to mat files.
dRaw = dir(paths.raw_occ);
rawSubjects = sort({dRaw([dRaw.isdir] & ~ismember({dRaw.name}, {'.', '..'})).name});
clc; disp('Loaded subjects:'); disp(str2double(string(rawSubjects)));
d = dir(fullfile(paths.raw_occ, '*', '*.cnt'));
d = d(~[d.isdir]);
nFiles = size(d, 1);
ids = {};
for f = 1 : nFiles
    filePath = fullfile(d(f).folder, d(f).name);
    id = strsplit(d(f).name, '_');
    subjectID = id{1};
    clc
    fprintf('[CUTTING] Subject GCP%s (%d/%d)\n', subjectID, f, nFiles)
    if subject_already_cut(d(f).folder, subjectID)
        fprintf('[CUTTING] Subject GCP%s (%d/%d) already cut. Skipping.\n', subjectID, f, nFiles)
        hasAsc = ~isempty(dir(fullfile(d(f).folder, '*.asc')));
        hasEdf = ~isempty(dir(fullfile(d(f).folder, '*.edf')));
        if hasAsc || hasEdf
            ids{end + 1} = subjectID; %#ok<SAGROW>
            fprintf('[CUTTING] Subject GCP%s (%d/%d) still has eye tracking files. Converting those.\n', ...
                subjectID, f, nFiles)
        end
        continue
    end
    GCP_cutData(filePath) % Perform cutting
    ids{end + 1} = subjectID; %#ok<SAGROW>
    fprintf('[CUTTING] Subject GCP%s (%d/%d) done\n', subjectID, f, nFiles)
end
ids = unique(ids, 'stable');

%% Load and synchronize EEG & Eyelink
nSync = length(ids);
for id = 1 : nSync
    ID = ids{id};
    clc
    fprintf('[CUTTING] Subject GCP%s (%d/%d)\n', ID, id, nSync)

    filePath = fullfile(paths.raw_occ, ID);
    d = dir(fullfile(filePath, '*.asc'));
    
    if not(isempty(d))
        for f = 1 : size(d, 1)
            % Convert ET asc to mat
            inputFile = fullfile(d(f).folder, d(f).name);
            
            % Rename the files (EyeLink can't deal with long filenames, edf filenames has to be short)
            x = strsplit(d(f).name, '_');
            name = x{2};
            id = x{end-1};
            block = str2double(name(1));

            if isnan(block)
                if strcmp(name, 'Res.asc')
                    task = 'Resting';
                    eegfile = [id, '_Resting_EEG.mat'];
                    etfile = [id, '_Resting_ET.mat'];
                else
                    task = 'Training';
                    continue
                end
            else
                task = 'GCP';
                etfile = [id, '_' task, '_block', num2str(block), '_task_ET.mat'];
            end
            outputFile = fullfile(d(f).folder, etfile);
            ET = parseeyelink(inputFile, outputFile);
        end
    end

    % Move asc files to archive
    d = dir(fullfile(filePath, '*.asc'));
    for f = 1 : size(d, 1)
        source = fullfile(d(f).folder, d(f).name);
        destination = fullfile(fullfile(filePath, 'archive'));
        movefile(source,destination)
    end

    % Move edf files to archive
    d = dir(fullfile(filePath, '*.edf'));
    for f = 1 : size(d, 1)
        source = fullfile(d(f).folder, d(f).name);
        destination = fullfile(fullfile(filePath, 'archive'));
        movefile(source,destination)
    end
end
fprintf('[CUTTING] SYNCHRONIZATION COMPLETE\n')

function tf = subject_already_cut(subjectDir, subjectID)
names = { ...
    [subjectID '_GCP_Resting_EEG.mat'] ...
    [subjectID '_GCP_block1_task_EEG.mat'] ...
    [subjectID '_GCP_block2_task_EEG.mat'] ...
    [subjectID '_GCP_block3_task_EEG.mat'] ...
    [subjectID '_GCP_block4_task_EEG.mat']};
tf = false;
for k = 1:numel(names)
    inFolder = isfile(fullfile(subjectDir, names{k}));
    inArchive = isfile(fullfile(subjectDir, 'archive', names{k}));
    if inFolder || inArchive
        tf = true;
        return
    end
end
end
