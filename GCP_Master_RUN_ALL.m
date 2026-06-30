%% GCP Master Analysis
%
% Runs the GCP pipeline: preprocessing, feature extraction, visualizations,
% and creation of CSVs for stats. Change the paths to your file locations first.
%
% IMPORTANT: Scripts call startup/clear, which wipe the workspace. Use
% absolute paths for the 'run' commands.
%
% For raw data, DO THIS FIRST: in 1_preprocessing, run in order
%   (1) GCP_doCutting.m (1_cut) to cut the raw ANT .cnt and EyeLink .asc data,
%   (2) Automagic (refer to README for settings),
%   (3) GCP_mergeData.m (3_merge) to merge Automagic EEG with the ET files

%% Setup
startup

%% 1_preprocessing/4_preprocessing FieldTrip
run('/Users/Arne/Documents/GitHub/GCP/1_preprocessing/4_preprocessing/GCP_preprocessing.m');

%% 2 Subject-level feature extraction
% Behavioral
run('/Users/Arne/Documents/GitHub/GCP/2_feature_extraction/GCP_behavioral_fex.m');

% Gaze
run('/Users/Arne/Documents/GitHub/GCP/2_feature_extraction/GCP_gaze_fex.m');

% EEG (GED)
run('/Users/Arne/Documents/GitHub/GCP/2_feature_extraction/GCP_eeg_fex_GED.m');
run('/Users/Arne/Documents/GitHub/GCP/2_feature_extraction/GCP_eeg_fex_GED_TFR.m');

% Master Matrices (CSVs)
run('/Users/Arne/Documents/GitHub/GCP/2_feature_extraction/GCP_master_matrix.m');
run('/Users/Arne/Documents/GitHub/GCP/2_feature_extraction/GCP_master_matrix_trials.m');

%% 3 Visualization
% Behavioral
run('/Users/Arne/Documents/GitHub/GCP/3_visualization/behavioral/GCP_behav.m');

% EEG (GED): power spectra, TFRs, trial-level boxplots
run('/Users/Arne/Documents/GitHub/GCP/3_visualization/eeg/powspctr/GCP_eeg_powspctrm_GED.m');
run('/Users/Arne/Documents/GitHub/GCP/3_visualization/eeg/tfr/GCP_TFR_GED.m');
run('/Users/Arne/Documents/GitHub/GCP/3_visualization/eeg/GCP_eeg_GED_trial_boxplots.m');

% Gaze time courses
run('/Users/Arne/Documents/GitHub/GCP/3_visualization/gaze/microsaccades/GCP_gaze_microsaccades_TC.m');
run('/Users/Arne/Documents/GitHub/GCP/3_visualization/gaze/pupilSize/GCP_gaze_pupilSize_TC.m');
run('/Users/Arne/Documents/GitHub/GCP/3_visualization/gaze/velocity/GCP_gaze_velocity_TC.m');
run('/Users/Arne/Documents/GitHub/GCP/3_visualization/gaze/fixations/GCP_gaze_fixations_TC.m');

% Hypotheses schematic
run('/Users/Arne/Documents/GitHub/GCP/3_visualization/hypotheses/GCP_hypotheses_plot.m');

%% 4 Stats
% Subject-level overview and boxplots (MATLAB)
run('/Users/Arne/Documents/GitHub/GCP/4_stats/GCP_stats_overview.m');
run('/Users/Arne/Documents/GitHub/GCP/4_stats/GCP_stats_boxplots.m');

% Hypothesis testing on trial-level data (MATLAB)
run('/Users/Arne/Documents/GitHub/GCP/4_stats/GCP_hypotheses_trials.m');

% Rainclouds: run the Python script GCP_stats_rainclouds.py separately.

%%
disp(datestr(now))
