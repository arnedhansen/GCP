%% GCP Eye Velocity

%% Setup
startup
[subjects, path, colors, headmodel] = setup(projectName, 0);

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/gaze');
    load([datapath, filesep 'dataET'])
