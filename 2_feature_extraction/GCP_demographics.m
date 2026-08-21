%% GCP Demographics
% Participant summary statistics from the merged subject-level table.
% Output: console summary of N, age mean and SD, female percent, and
% right-handed percent for the full sample and the GED-included subsample.

startup
[~, paths] = setup('GCP', 0);
load(fullfile(paths.features, 'GCP_merged_data.mat'));
if exist('merged_table', 'var')
    T = merged_table;
else
    T = struct2table(merged_data);
end

[~, idx] = unique(T.ID, 'stable');
T = T(idx, :);

print_demographics(T, 'Full sample');
if ismember('Include', T.Properties.VariableNames)
    print_demographics(T(T.Include == 1, :), 'GED-included sample');
end

function print_demographics(dat, label)
if isempty(dat)
    fprintf('%s: N = 0\n', label);
    return
end

ages = dat.Age;
genders = dat.Gender;
if isstring(genders) || ischar(genders)
    genders = cellstr(genders);
end

mean_age = mean(ages, 'omitnan');
sd_age = std(ages, 'omitnan');
n_female = sum(strcmp(genders, 'W'));
perc_female = (n_female / numel(genders)) * 100;

fprintf('%s\n', label);
fprintf('N: %.f participants\n', height(dat));
fprintf('Mean age: %.2f years\n', mean_age);
fprintf('SD age: %.2f years\n', sd_age);
fprintf('Female: %.1f%%\n', perc_female);

if ismember('Handedness', dat.Properties.VariableNames)
    handedness = dat.Handedness;
    if isstring(handedness) || ischar(handedness)
        handedness = cellstr(handedness);
    end
    n_right = sum(strcmp(handedness, 'R'));
    perc_right = (n_right / numel(handedness)) * 100;
    fprintf('Right-handed: %.1f%%\n', perc_right);
end
fprintf('\n');
end
