function [tfr_cond_trials, tfr_cond_avg, ged_filter_meta] = compute_ged_tfr( ...
    subjects, paths, all_combined_filter_full, all_topo_labels, ...
    all_component_selection_stats_full, ...
    baseline_window, ...
    tfr_foi, tfr_toi, tfr_win_sec, tfr_tapsmofrq, ...
    condNames, condCodes)
% Project trials through each subject GED filter and compute condition TFRs in dB.
% Prefer compute_ged_tfr_subject from the main loop to avoid reloading EEG.
nSubj = numel(subjects);
nCond = numel(condNames);
tfr_cond_trials = cell(nCond, nSubj);
tfr_cond_avg = cell(nCond, nSubj);
ged_filter_meta = cell(1, nSubj);

for subj = 1:nSubj
    clc; fprintf('[GED TFR] Subject %s (%d/%d)\n', subjects{subj}, subj, nSubj);

    subj_eeg_path = fullfile(paths.features, subjects{subj}, 'eeg', 'dataEEG.mat');
    if ~isfile(subj_eeg_path)
        warning('Missing EEG file for %s: %s', subjects{subj}, subj_eeg_path);
        continue;
    end
    E = load(subj_eeg_path, 'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    dat_by_cond = {E.dataEEG_c25, E.dataEEG_c50, E.dataEEG_c75, E.dataEEG_c100};

    combined_filter = all_combined_filter_full{subj};
    if isempty(combined_filter)
        warning('No valid combined GED filter for %s. Skipping subject.', subjects{subj});
        continue;
    end
    saved_labels = all_topo_labels{subj};
    if isempty(saved_labels) || isempty(dat_by_cond{1}) || ~isfield(dat_by_cond{1}, 'label')
        warning('Missing channel labels for %s. Skipping subject.', subjects{subj});
        continue;
    end
    current_labels = dat_by_cond{1}.label;
    [labels_found, filter_order] = ismember(current_labels, saved_labels);
    if ~all(labels_found) || numel(combined_filter) ~= numel(saved_labels)
        warning('GED filter channels do not match the EEG channels for %s. Skipping subject.', subjects{subj});
        continue;
    end
    combined_filter = combined_filter(filter_order);

    stat_full = all_component_selection_stats_full{subj};
    if isempty(stat_full)
        stat_full = struct();
    end

    [tfr_cond_trials(:, subj), tfr_cond_avg(:, subj), ged_filter_meta{subj}] = ...
        compute_ged_tfr_subject( ...
        dat_by_cond, combined_filter, current_labels, subjects{subj}, subj, ...
        stat_full, baseline_window, tfr_foi, tfr_toi, tfr_win_sec, tfr_tapsmofrq, condCodes);
end
end
