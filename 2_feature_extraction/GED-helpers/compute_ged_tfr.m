function [tfr_cond_trials, tfr_cond_avg, ged_filter_meta] = compute_ged_tfr( ...
    subjects, paths, all_combined_filter_full, all_topo_labels, ...
    all_component_selection_stats_full, ...
    baseline_window, ...
    tfr_foi, tfr_toi, tfr_win_sec, tfr_tapsmofrq, ...
    condNames, condCodes)
% Project trials through each subject GED filter and compute condition TFRs in dB.
nSubj = numel(subjects);
nCond = numel(condNames);
tfr_cond_trials = cell(nCond, nSubj);
tfr_cond_avg = cell(nCond, nSubj);
ged_filter_meta = cell(1, nSubj);

tfr_baseline_window = [baseline_window(1) + tfr_win_sec / 2, ...
    baseline_window(2) - tfr_win_sec / 2];
if tfr_baseline_window(1) > tfr_baseline_window(2)
    error('The TFR window is too long for the requested baseline interval.');
end

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
    ged_filter_meta{subj} = struct( ...
        'subject', subjects{subj}, ...
        'source_subject_index', subj, ...
        'channel_labels', {current_labels}, ...
        'combined_filter', combined_filter, ...
        'component_selection', stat_full);

    for c = 1:nCond
        dat = dat_by_cond{c};
        if isempty(dat) || ~isfield(dat, 'trial') || isempty(dat.trial)
            continue;
        end

        trl_idx = find(dat.trialinfo == condCodes(c));
        if isempty(trl_idx)
            continue;
        end
        cfg_sel = [];
        cfg_sel.trials = trl_idx;
        dat = ft_selectdata(cfg_sel, dat);

        ged_dat = dat;
        ged_dat.label = {'GED'};
        for tr = 1:numel(dat.trial)
            x = double(dat.trial{tr});
            ged_dat.trial{tr} = combined_filter(:)' * x;
        end

        cfg_tfr = [];
        cfg_tfr.method = 'mtmconvol';
        cfg_tfr.output = 'pow';
        cfg_tfr.taper = 'dpss';
        cfg_tfr.foi = tfr_foi;
        cfg_tfr.toi = tfr_toi;
        cfg_tfr.t_ftimwin = tfr_win_sec * ones(size(tfr_foi));
        cfg_tfr.tapsmofrq = tfr_tapsmofrq * ones(size(tfr_foi));
        cfg_tfr.keeptrials = 'yes';
        tfr_trials = ft_freqanalysis(cfg_tfr, ged_dat);

        cfg_baseline = [];
        cfg_baseline.baseline = tfr_baseline_window;
        cfg_baseline.baselinetype = 'db';
        cfg_baseline.parameter = 'powspctrm';
        tfr_trials = ft_freqbaseline(cfg_baseline, tfr_trials);

        cfg_avg = [];
        cfg_avg.keeptrials = 'no';
        tfr_avg = ft_freqdescriptives(cfg_avg, tfr_trials);

        tfr_cond_trials{c, subj} = tfr_trials;
        tfr_cond_avg{c, subj} = tfr_avg;
    end
end
end
