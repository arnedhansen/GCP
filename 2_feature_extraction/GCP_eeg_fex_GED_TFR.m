%% GCP GED-based Time-Frequency Representation (TFR)
%
% This script applies the exact subject-specific combined GED spatial filter
% saved by the GED pipeline and computes trial-resolved TFRs per condition.
%
% Output:
%   data/features/GCP_eeg_GED_TFR.mat
%
% Notes:
%   - Signed, noise-normalized combined GED filters are loaded from
%     GCP_eeg_GED.mat and matched to the current EEG channel order.
%   - Baseline correction is applied to linear power before dB conversion:
%       power_change(t) = 10*log10(power(t) / mean(power(baseline))).

%% Setup
startup
[subjects, paths, ~] = setup('GCP');

%% Parameters
baseline_window = [-1.5, -0.5];
stim_window = [0, 2.0];

% TFR settings
tfr_foi = 30:1:90;
tfr_toi = -1.75:0.05:2.00;
tfr_win_sec = 0.50;
tfr_tapsmofrq = 5;
tfr_baseline_window = [baseline_window(1) + tfr_win_sec / 2, ...
    baseline_window(2) - tfr_win_sec / 2];
if tfr_baseline_window(1) > tfr_baseline_window(2)
    error('The TFR window is too long for the requested baseline interval.');
end

condNames = {'c25', 'c50', 'c75', 'c100'};
condCodes = [61, 62, 63, 64];
condLabels = {'25%', '50%', '75%', '100%'};

%% Paths
ged_path = fullfile(paths.features, 'GCP_eeg_GED.mat');
if ~isfile(ged_path)
    error('Missing GED features file: %s', ged_path);
end

Sged = load(ged_path, ...
    'all_combined_filter_full', ...
    'all_topo_labels', ...
    'all_component_selection_stats_full', ...
    'subjects');
if ~isfield(Sged, 'all_combined_filter_full') || isempty(Sged.all_combined_filter_full)
    error(['Combined GED filters are missing in %s. Rerun ' ...
        'GCP_eeg_fex_GED.m before computing GED TFRs.'], ged_path);
end

out_path = fullfile(paths.features, 'GCP_eeg_GED_TFR.mat');

%% Storage
nSubj = numel(subjects);
nCond = numel(condNames);

tfr_cond_trials = cell(nCond, nSubj);
tfr_cond_avg = cell(nCond, nSubj);
ged_filter_meta = cell(1, nSubj);

%% Subject loop
for subj = 1:nSubj
    fprintf('GED-TFR Subject %s (%d/%d)\n', subjects{subj}, subj, nSubj);

    % Load subject EEG
    subj_eeg_path = fullfile(paths.features, subjects{subj}, 'eeg', 'dataEEG.mat');
    if ~isfile(subj_eeg_path)
        warning('Missing EEG file for %s: %s', subjects{subj}, subj_eeg_path);
        continue;
    end
    E = load(subj_eeg_path, 'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    dat_by_cond = {E.dataEEG_c25, E.dataEEG_c50, E.dataEEG_c75, E.dataEEG_c100};

    source_idx = find(strcmp(Sged.subjects, subjects{subj}), 1);
    if isempty(source_idx)
        warning('Subject %s is absent from the GED feature file. Skipping subject.', subjects{subj});
        continue;
    end
    combined_filter = Sged.all_combined_filter_full{source_idx};
    if isempty(combined_filter)
        warning('No valid combined GED filter for %s. Skipping subject.', subjects{subj});
        continue;
    end
    saved_labels = Sged.all_topo_labels{source_idx};
    current_labels = dat_by_cond{1}.label;
    [labels_found, filter_order] = ismember(current_labels, saved_labels);
    if ~all(labels_found) || numel(combined_filter) ~= numel(saved_labels)
        warning('GED filter channels do not match the EEG channels for %s. Skipping subject.', subjects{subj});
        continue;
    end
    combined_filter = combined_filter(filter_order);

    stat_full = Sged.all_component_selection_stats_full{source_idx};
    if isempty(stat_full)
        stat_full = struct();
    end
    ged_filter_meta{subj} = struct( ...
        'subject', subjects{subj}, ...
        'source_subject_index', source_idx, ...
        'channel_labels', {current_labels}, ...
        'combined_filter', combined_filter, ...
        'component_selection', stat_full);

    % Condition loop
    for c = 1:nCond
        dat = dat_by_cond{c};
        if isempty(dat) || ~isfield(dat, 'trial') || isempty(dat.trial)
            continue;
        end

        % Keep only task trials matching condition code
        trl_idx = find(dat.trialinfo == condCodes(c));
        if isempty(trl_idx)
            continue;
        end
        cfg_sel = [];
        cfg_sel.trials = trl_idx;
        dat = ft_selectdata(cfg_sel, dat);

        % Project each trial through the exact combined GED filter.
        ged_dat = dat;
        ged_dat.label = {'GED'};
        for tr = 1:numel(dat.trial)
            x = double(dat.trial{tr}); % chan x time
            ged_dat.trial{tr} = combined_filter(:)' * x;
        end

        % Trial-resolved TFR
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

        % Subject-condition average
        cfg_avg = [];
        cfg_avg.keeptrials = 'no';
        tfr_avg = ft_freqdescriptives(cfg_avg, tfr_trials);

        tfr_cond_trials{c, subj} = tfr_trials;
        tfr_cond_avg{c, subj} = tfr_avg;
    end
end

%% Save
save(out_path, ...
    'tfr_cond_trials', 'tfr_cond_avg', ...
    'ged_filter_meta', ...
    'subjects', 'condNames', 'condLabels', ...
    'baseline_window', 'tfr_baseline_window', 'stim_window', ...
    'tfr_foi', 'tfr_toi', 'tfr_win_sec', 'tfr_tapsmofrq', ...
    '-v7.3');

fprintf('Saved GED-TFR data to: %s\n', out_path);
