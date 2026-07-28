%% GCP Test suite for GED microsaccade analyses
% Runs exploratory GED strategies plus confirmatory Engbert-locked test.

clear; clc; close all;

fprintf('Running GED microsaccade covariance regression test...\n');
run(fullfile(fileparts(mfilename('fullpath')), 'GCP_test_GED_microsaccade_spoc.m'));

fprintf('Running GED contrast modulation test...\n');
run(fullfile(fileparts(mfilename('fullpath')), 'GCP_test_GED_microsaccade_contrast_modulation.m'));

fprintf('Running GED event locked test (exploratory MS-vs-control GED)...\n');
run(fullfile(fileparts(mfilename('fullpath')), 'GCP_test_GED_microsaccade_eventlocked.m'));

fprintf('Running confirmatory GED event locked test (Engbert + fixed GED)...\n');
run(fullfile(fileparts(mfilename('fullpath')), 'GCP_test_GED_microsaccade_eventlocked_confirmatory.m'));

fprintf('GED microsaccade test suite finished.\n');
