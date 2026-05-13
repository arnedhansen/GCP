%% GCP Test suite for GED microsaccade analyses
% Runs all three exploratory GED strategies.

clear; clc; close all;

fprintf('Running GED microsaccade covariance regression test...\n');
run(fullfile(fileparts(mfilename('fullpath')), 'GCP_test_GED_microsaccade_spoc.m'));

fprintf('Running GED contrast modulation test...\n');
run(fullfile(fileparts(mfilename('fullpath')), 'GCP_test_GED_microsaccade_contrast_modulation.m'));

fprintf('Running GED event locked test...\n');
run(fullfile(fileparts(mfilename('fullpath')), 'GCP_test_GED_microsaccade_eventlocked.m'));

fprintf('GED microsaccade test suite finished.\n');
