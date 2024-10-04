%% Computation of GCP EyeLink CALIBRATION and VALIDATION data

%% Setup
clear
clc
close all
path = '/Volumes/methlab_data/GCP/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load gaze data
% Define the number of files for each pattern (N and S)
numFiles = 4;

% Loop through each subject
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/archive');
    cd(datapath);

    % Initialize
    fileCounter = 1;
    VAL{subj} = cell(1, numFiles);

    for fileNum = 1:numFiles
        % Create filename
        filePath = sprintf('%s_%d.asc', subjects{subj}, fileNum);

        % Load data file
        [lastValidationOffsets] = parseASCFile(filePath);

        % Save positions
        VAL{subj}{fileCounter} = lastValidationOffsets;

        % Increment file counter
        fileCounter = fileCounter + 1;
    end
    disp(['Subject ' num2str(subjects{subj}) ' loaded.'])
end

%% VISUALIZE for each subject
numSubjects = numel(VAL);
for subjIdx = 1:numSubjects
    try
        close all
        figure;
        set(gcf, 'Position', [300, 200, 1000, 600], 'Color', 'w');

        % Get data
        validationData = VAL{subjIdx};

        % Plot
        hold on;
        for fileIdx = 1:numel(validationData)
            % Extract validation values for the current file
            values = validationData{fileIdx};

            % Plot the values
            plot(fileIdx * ones(size(values)), values, 'bo');
        end

        % Add a red line at y=1
        yline(1, 'r', 'LineWidth', 2);
        yline(1.25, 'k', 'LineWidth', 0.5, 'LineStyle', '--');

        % Set the x-axis and y-axis labels
        xlabel('Sternberg Blocks');
        ylabel('Validation Values [dva]');

        % Set the x-axis limits
        xlim([0, 7]);

        % Set the x-axis limits, ticks and labels
        xticks(1:6);
        xticklabels({'1', '2', '3', '4', '5', '6'});

        % Set the y-axis limits and ticks
        numericArray = [validationData{:}];
        maxValueVAL = max(numericArray);
        if maxValueVAL < 1.25
            maxValueVAL = 1.25;
        end
        ylim([0, maxValueVAL*1.1])
        yticks(0:0.1:10);

        % Add a title
        title(['Validation Data for Subject ' num2str(subjects{subjIdx})], 'FontSize', 20);
        hold off;

        % Save
        savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/controls/',subjects{subjIdx});
        warning('off', 'MATLAB:MKDIR:DirectoryExists')
        mkdir(savepath)
        cd(savepath)
        saveName = [savepath, filesep, num2str(subjects{subjIdx}) '_Validations.png'];
        saveas(gcf, saveName)
    catch
        disp(['Could not create VALIDATION overview for subject ' num2str(subjects{subjIdx})])
    end
end

%% FUNCTION
function [lastValidationOffsets] = parseASCFile(filePath)
% Initialize cell arrays to store calibration and validation data
calibrationData = {};
validationData = {};

% Open the file for reading
fid = fopen(filePath, 'rt');

if fid == -1
    error('Cannot open file: %s', filePath);
end

% Read the file line by line
line = fgetl(fid);

currentCalibration = [];
currentValidation = [];

while ~contains(line, '!MODE RECORD CR')

    % Check for validation offset lines
    if contains(line, 'VALIDATE') && contains(line, 'POINT')
        tokens = regexp(line, 'OFFSET (\d+\.\d+) deg.', 'tokens');
        if ~isempty(tokens)
            offset = str2double(tokens{1}{1});
            currentValidation = [currentValidation, offset];
        end
    end

    % Read the next line
    line = fgetl(fid);
end

fclose(fid);
try
    lastValidationOffsets = currentValidation(end-8:end);
catch
    lastValidationOffsets = currentValidation;
end
end