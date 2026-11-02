%% GCP exclusion list
% Prints IDs in controls/GCP_exclusion_participants.rtf.
% Later scripts drop these IDs. This script does not write a table.
startup
[~, paths] = setup('GCP', 0);
ids = gcp_manual_exclusions(paths);
disp(ids);
