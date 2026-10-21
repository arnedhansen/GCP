%% GCP subject inclusion list (manual seed / override)
startup
[subjects, paths] = setup('GCP', 0);
if ispc
    addpath('C:\Users\Administrator\Documents\GitHub\GCP\_BACKUP\_PILOT');
else
    addpath('/Users/Arne/Documents/GitHub/GCP/_BACKUP/_PILOT');
end
[paths, subjects] = gcp_pilot_paths(paths);

SubjID = str2double(string(subjects));
Include = ~ismember(SubjID, [602, 604, 608]);
subject_inclusion = table(SubjID, Include, 'VariableNames', {'SubjID', 'Include'});
save(fullfile(paths.controls, 'GCP_subject_inclusion.mat'), 'subject_inclusion', '-v7.3');
disp(subject_inclusion);
