%% GCP subject inclusion list (manual seed / override)
startup
[subjects, paths] = setup('GCP', 0);

SubjID = str2double(string(subjects));
Include = ~ismember(SubjID, [602, 604, 608]);
subject_inclusion = table(SubjID, Include, 'VariableNames', {'SubjID', 'Include'});
save(fullfile(paths.controls, 'GCP_subject_inclusion.mat'), 'subject_inclusion', '-v7.3');
disp(subject_inclusion);
