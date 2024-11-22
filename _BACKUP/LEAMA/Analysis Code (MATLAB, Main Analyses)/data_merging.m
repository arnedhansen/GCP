%% EEG and Eye Tracker Data Merging %%

%% Description %%
% This script merges the eye tracking and EEG data for these blocks. The
% output are EEGmerged.mat files with indicators of the subject and EEG data quality.


%% Script Structure
% 0. Preparation of Paths
% 1. Extraction of EEG Data Files
% 2. Merging of EEG and eye tracker data



%% 0. Preparation

%% 0.1 Set Current Directory and Clear Figures
p = pwd;
close()
cd(p)

%% 0.2 Add Path for Eye-EEG-Master Toolbox and initialize folder for results
addpath(genpath('\\psyger-stor02.d.uzh.ch\methlab\Students\Arne\MA/toolboxes\eye-eeg-master'));
resultFolder = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results';


%% 1. Define Path and Extract EEG Data Files
% Define the path for EEG data files and extract relevant information
dEEG = dir(strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results\*\', '*ip*EEG.mat'));

% Initialize cell arrays for storing file information
ids = {};
idfile = {};

% Extract unique participant IDs from file names
for f = 1 : size(dEEG, 1)
    filePath = fullfile(dEEG(f).folder, dEEG(f).name);
    idfile{f} = strsplit(dEEG(f).name, '_');
    ids{f} = idfile{f}{2};
end

ids_unique = unique(ids);



%% 2. Load and merge EEG & EyeLink Data

% Define start and end markers of the blocks
START = [11, 12, 13, 14, 10];
END = [91 92 93 94 90];

% For each participant...
for id = 1 : length(ids_unique) %-(length(ids_unique)-1)
    ID = ids_unique{id};

    % Define path for EyeLink data
    filePath = fullfile('\\psyger-stor02.d.uzh.ch\methlab_data\OCC\LEAMA\data', ID);
    dET = dir([filePath, filesep, '*ET.mat']);

    % Set the mask for finding the blocks and prepare the numbers array
    dEEGsub = dir((strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results\', num2str(ID), '\', '*ip*EEG.mat')));
    numbers  = [];

    % Extract block names for existing blocks
    for block = 1:length(dEEGsub)
        all_num = regexp(dEEGsub(block).name, '\d+', 'match');
        numbers(block) = str2num(all_num{2});
    end

    % Process EEG and EyeLink data for each block
    for f = numbers
        try
            % Find the EEG block index
            ind = find(numbers == f);

            % Load EyeLink data
            ETFile =  load(fullfile(dET(f).folder, dET(f).name));

            % Load EEG and merge with EyeLink data
            load(fullfile(dEEGsub(ind).folder, dEEGsub(ind).name));
            EEG = pop_importeyetracker(EEG, ETFile, [START(f) END(f)], [2 3 4],{'L_GAZE_X', 'L_GAZE_Y', 'L_AREA'}, 1, 1, 1, 0);
            mkdir(fullfile(resultFolder, ID))

            % Save the merged EEG data
            if f < 5
                save(fullfile(resultFolder, ID, strcat(dEEGsub(ind).name(1:3), '_EEGmerged', num2str(f))), 'EEG', '-v7.3')
            else % Block 5 is the Resting State block
                save(fullfile(resultFolder, ID, strcat(dEEGsub(ind).name(1:3), 'EEGmergedRes')), 'EEG', '-v7.3')
            end

        catch ME
            % Display an error message if an exception occurs
            ME.message
        end
    end
end % End participant loop


