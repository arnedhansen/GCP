%% GCP subject inclusion list (manual seed / override)
% Every discovered subject starts included. A subject already excluded
% because fewer than two blocks remained keeps that exclusion and its reason.
startup
[subjects, paths] = setup('GCP', 0);

reason = gcp_inclusion_table('reason');
incFile = fullfile(paths.controls, 'GCP_subject_inclusion.mat');
if isfile(incFile)
    S = load(incFile, 'subject_inclusion');
    old = S.subject_inclusion;
else
    old = table();
end
if ~ismember('Reason', old.Properties.VariableNames)
    old.Reason = repmat("", height(old), 1);
end
if isempty(old)
    old.SubjID = zeros(0, 1);
    old.Include = false(0, 1);
    old.Reason = strings(0, 1);
end
old.SubjID = double(old.SubjID);
old.Reason = string(old.Reason);

SubjID = str2double(string(subjects(:)));
Include = true(numel(SubjID), 1);
Reason = repmat("", numel(SubjID), 1);
for i = 1:numel(SubjID)
    row = find(old.SubjID == SubjID(i), 1);
    if ~isempty(row) && old.Reason(row) == reason
        Include(i) = false;
        Reason(i) = reason;
    end
end
subject_inclusion = table(SubjID, Include, Reason);
save(incFile, 'subject_inclusion', '-v7.3');
disp(subject_inclusion);
