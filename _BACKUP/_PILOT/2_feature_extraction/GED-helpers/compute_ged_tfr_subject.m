function [tfr_cond_trials, tfr_cond_avg, ged_filter_meta] = compute_ged_tfr_subject( ...
    dat_per_cond, combined_filter, channel_labels, subject_id, source_subject_index, ...
    component_selection_stats, ...
    baseline_window, tfr_foi, tfr_toi, tfr_win_sec, tfr_tapsmofrq, condCodes)
% GED-projected condition TFRs for one subject using in-memory trial data.
nCond = numel(dat_per_cond);
tfr_cond_trials = cell(nCond, 1);
tfr_cond_avg = cell(nCond, 1);
ged_filter_meta = struct( ...
    'subject', subject_id, ...
    'source_subject_index', source_subject_index, ...
    'channel_labels', {channel_labels}, ...
    'combined_filter', combined_filter(:), ...
    'component_selection', component_selection_stats);

tfr_baseline_window = [baseline_window(1) + tfr_win_sec / 2, ...
    baseline_window(2) - tfr_win_sec / 2];
if tfr_baseline_window(1) > tfr_baseline_window(2)
    error('The TFR window is too long for the requested baseline interval.');
end
if isempty(combined_filter) || isempty(channel_labels)
    return;
end
combined_filter = combined_filter(:);
if numel(combined_filter) ~= numel(channel_labels)
    warning('GED filter channels do not match the EEG channels for %s. Skipping TFR.', subject_id);
    return;
end

for c = 1:nCond
    dat = dat_per_cond{c};
    if isempty(dat) || ~isfield(dat, 'trial') || isempty(dat.trial)
        continue;
    end
    if nargin >= 12 && ~isempty(condCodes) && isfield(dat, 'trialinfo')
        trl_idx = find(dat.trialinfo == condCodes(c));
        if isempty(trl_idx)
            continue;
        end
        if ~isequal(trl_idx(:).', 1:numel(dat.trial))
            dat = select_trials_light(dat, trl_idx);
        end
    end

    ged_dat = dat;
    ged_dat.label = {'GED'};
    for tr = 1:numel(dat.trial)
        x = double(dat.trial{tr});
        ged_dat.trial{tr} = combined_filter' * x;
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
    cfg_tfr.feedback = 'none';
    tfr_trials = ft_freqanalysis(cfg_tfr, ged_dat);

    cfg_baseline = [];
    cfg_baseline.baseline = tfr_baseline_window;
    cfg_baseline.baselinetype = 'db';
    cfg_baseline.parameter = 'powspctrm';
    tfr_trials = ft_freqbaseline(cfg_baseline, tfr_trials);

    cfg_avg = [];
    cfg_avg.keeptrials = 'no';
    tfr_avg = ft_freqdescriptives(cfg_avg, tfr_trials);

    tfr_cond_trials{c} = tfr_trials;
    tfr_cond_avg{c} = tfr_avg;
end
end

function dat_out = select_trials_light(dat_in, trl_idx)
dat_out = dat_in;
dat_out.trial = dat_in.trial(trl_idx);
dat_out.time = dat_in.time(trl_idx);
if isfield(dat_in, 'trialinfo')
    dat_out.trialinfo = dat_in.trialinfo(trl_idx, :);
end
if isfield(dat_in, 'sampleinfo')
    dat_out.sampleinfo = dat_in.sampleinfo(trl_idx, :);
end
end
