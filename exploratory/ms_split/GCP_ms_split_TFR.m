%% GCP high minus low microsaccade rate difference TFR
%
% Same median split as GCP_ms_split_gamma.m. Within each participant and
% each contrast, trials are split at the median stimulus window
% microsaccade rate (0 to 2 s, Hz). Trials equal to that median are
% dropped. A contrast is omitted when either side has fewer than
% min_trials_per_bin trials with a finite GED spectrum.
%
% EEG trials are projected through the saved GED filter. The filter is
% not refit. The multitaper and dB baseline match the confirmatory TFR
% in GCP_eeg_fex_GED.m. Each trial is baseline corrected, then the trials
% in a bin are averaged. The map is high minus low.
%
% This is one GED virtual channel. It is not a sensor level TFR.
%
% Figures (figures/controls/ms_rate_median_split):
%   GCP_split_tfr_diff.png
%   GCP_split_tfr_diff_pooled.png
%   GCP_ms_split_tfr.mat
%
% Inputs:
%   paths.features/GCP_eeg_GED.mat
%       trials_powratio, scan_freqs, subjects,
%       all_combined_filter, all_topo_labels
%   paths.features/<subj>/eeg/dataEEG.mat
%   paths.features/<subj>/gaze/gaze_matrix_trial.mat (MSRate, Trial)

%% Setup
startup
[subjects, paths, ~, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
nSubj = numel(subjects);
nCond = 4;
condPct = [25 50 75 100];
condCodes = [61 62 63 64];
eegFields = {'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100'};
gazeFields = {'subj_data_gaze_trial_c25', 'subj_data_gaze_trial_c50', ...
    'subj_data_gaze_trial_c75', 'subj_data_gaze_trial_c100'};
condTitles = {'25% Contrast', '50% Contrast', '75% Contrast', '100% Contrast'};
min_trials_per_bin = 8;

% Matched to GCP_eeg_fex_GED.m
baseline_window = [-1.5 -0.5];
tfr_foi = 30:1:90;
tfr_toi = -1.75:0.05:2.00;
tfr_win_sec = 0.50;
tfr_tapsmofrq = 5;
tfr_baseline_window = [baseline_window(1) + tfr_win_sec / 2, ...
    baseline_window(2) - tfr_win_sec / 2];
if tfr_baseline_window(1) > tfr_baseline_window(2)
    error('GCP_ms_split_TFR:Baseline', ...
        'The TFR window is too long for the requested baseline interval.');
end
plot_latency = [-0.5 2];
plot_frequency = [30 90];

fig_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/controls/ms_rate_median_split';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

fprintf('[MS SPLIT TFR] Included cohort N = %d (%s)\n', nSubj, strjoin(subjects, ', '));

%% Load saved GED filters and the trial spectra used for the split gate
ged = load(fullfile(paths.features, 'GCP_eeg_GED.mat'), ...
    'trials_powratio', 'scan_freqs', 'subjects', ...
    'all_combined_filter', 'all_topo_labels');
scan_freqs = ged.scan_freqs(:)';
nFreq = numel(scan_freqs);
ged_idx = match_subjects(ged.subjects, subjects);

tfr_low = cell(nCond, nSubj);
tfr_high = cell(nCond, nSubj);
n_low = nan(nCond, nSubj);
n_high = nan(nCond, nSubj);
n_tie = nan(nCond, nSubj);
n_zero = nan(nCond, nSubj);
n_usable = nan(nCond, nSubj);
ms_median = nan(nCond, nSubj);

%% Median split and bin TFRs
for s = 1:nSubj
    gi = ged_idx(s);
    filt = ged.all_combined_filter{gi};
    saved_labels = ged.all_topo_labels{gi};
    if isempty(filt) || isempty(saved_labels)
        fprintf('[MS SPLIT TFR] %s: no saved GED filter, skipped\n', subjects{s});
        continue
    end
    filt = filt(:);

    eeg_file = fullfile(paths.features, subjects{s}, 'eeg', 'dataEEG.mat');
    gaze_file = fullfile(paths.features, subjects{s}, 'gaze', 'gaze_matrix_trial.mat');
    if ~isfile(eeg_file)
        error('GCP_ms_split_TFR:MissingEEG', ...
            'EEG file missing for %s: %s', subjects{s}, eeg_file);
    end
    if ~isfile(gaze_file)
        error('GCP_ms_split_TFR:MissingGaze', ...
            'Trial gaze file missing for %s: %s', subjects{s}, gaze_file);
    end
    eeg = load(eeg_file, eegFields{:});
    gaze = load(gaze_file);

    current_labels = {};
    for c = 1:nCond
        if isfield(eeg, eegFields{c}) && isstruct(eeg.(eegFields{c})) && ...
                isfield(eeg.(eegFields{c}), 'label')
            current_labels = eeg.(eegFields{c}).label;
            break
        end
    end
    if isempty(current_labels)
        error('GCP_ms_split_TFR:MissingLabels', ...
            'No EEG channel labels for %s.', subjects{s});
    end
    [labels_found, filter_order] = ismember(current_labels, saved_labels);
    if ~all(labels_found) || numel(filt) ~= numel(saved_labels)
        error('GCP_ms_split_TFR:FilterMismatch', ...
            'GED filter channels do not match the EEG channels for %s.', subjects{s});
    end
    filt = filt(filter_order);

    for c = 1:nCond
        pr = ged.trials_powratio{c, gi};
        if isempty(pr)
            fprintf('[MS SPLIT TFR] %s %d%%: empty GED spectra, skipped\n', ...
                subjects{s}, condPct(c));
            continue
        end
        pr = orient_powratio(pr, nFreq);
        if ~isfield(eeg, eegFields{c}) || isempty(eeg.(eegFields{c}))
            error('GCP_ms_split_TFR:MissingCond', ...
                '%s has no %s', subjects{s}, eegFields{c});
        end
        dat = eeg.(eegFields{c});
        if ~isfield(dat, 'trialinfo')
            error('GCP_ms_split_TFR:MissingTrialinfo', ...
                '%s %s has no trialinfo.', subjects{s}, eegFields{c});
        end
        trlIdx = find(dat.trialinfo == condCodes(c));
        if isempty(trlIdx)
            fprintf('[MS SPLIT TFR] %s %d%%: no trials with code %d, skipped\n', ...
                subjects{s}, condPct(c), condCodes(c));
            continue
        end
        cfg_trl = [];
        cfg_trl.trials = trlIdx(:)';
        dat = ft_selectdata(cfg_trl, dat);
        nTrl = numel(dat.trial);
        if size(pr, 1) ~= nTrl
            error('GCP_ms_split_TFR:TrialMismatch', ...
                ['%s %d%%: GED spectra have %d trials and EEG has %d ' ...
                'after the condition code selection.'], ...
                subjects{s}, condPct(c), size(pr, 1), nTrl);
        end
        if ~isfield(gaze, gazeFields{c})
            error('GCP_ms_split_TFR:MissingGazeCond', ...
                '%s has no %s', subjects{s}, gazeFields{c});
        end
        [ms, aligned] = align_gaze_to_trials(gaze.(gazeFields{c}), nTrl);
        if ~aligned
            error('GCP_ms_split_TFR:TrialMismatch', ...
                ['%s %d%%: gaze trials do not align with EEG trials ' ...
                '(nEEG = %d).'], subjects{s}, condPct(c), nTrl);
        end

        spec_ok = any(isfinite(pr), 2);
        usable = isfinite(ms) & spec_ok;
        n_usable(c, s) = nnz(usable);
        n_zero(c, s) = nnz(usable & ms == 0);
        if n_usable(c, s) < 1
            fprintf('[MS SPLIT TFR] %s %d%%: no usable trials, skipped\n', ...
                subjects{s}, condPct(c));
            continue
        end

        med = median(ms(usable));
        ms_median(c, s) = med;
        low = usable & ms < med;
        high = usable & ms > med;
        n_low(c, s) = nnz(low);
        n_high(c, s) = nnz(high);
        n_tie(c, s) = nnz(usable & ms == med);
        pct_zero = 100 * n_zero(c, s) / n_usable(c, s);
        if n_low(c, s) < min_trials_per_bin || n_high(c, s) < min_trials_per_bin
            fprintf('[MS SPLIT TFR] %s %d%%: 0 MS %.0f%%, low %d, high %d, ties %d, skipped\n', ...
                subjects{s}, condPct(c), pct_zero, n_low(c, s), n_high(c, s), n_tie(c, s));
            continue
        end

        tfr_low{c, s} = drop_cfg(ged_bin_tfr(dat, find(low), filt, ...
            tfr_foi, tfr_toi, tfr_win_sec, tfr_tapsmofrq, tfr_baseline_window));
        tfr_high{c, s} = drop_cfg(ged_bin_tfr(dat, find(high), filt, ...
            tfr_foi, tfr_toi, tfr_win_sec, tfr_tapsmofrq, tfr_baseline_window));
        fprintf('[MS SPLIT TFR] %s %d%%: 0 MS %.0f%%, low %d, high %d, ties dropped %d\n', ...
            subjects{s}, condPct(c), pct_zero, n_low(c, s), n_high(c, s), n_tie(c, s));
    end
    clear eeg gaze dat
end

%% Grand average of high minus low, per contrast
grand_diff = cell(1, nCond);
n_grand = zeros(1, nCond);
for c = 1:nCond
    pair_low = tfr_low(c, :);
    pair_high = tfr_high(c, :);
    keep = ~cellfun(@isempty, pair_low) & ~cellfun(@isempty, pair_high);
    n_grand(c) = nnz(keep);
    if n_grand(c) < 1
        fprintf('[MS SPLIT TFR] %d%%: no subjects in the grand average\n', condPct(c));
        continue
    end
    grand_diff{c} = high_minus_low(pair_high(keep), pair_low(keep));
    fprintf('[MS SPLIT TFR] %d%% grand average N = %d\n', condPct(c), n_grand(c));
end

%% One map: mean of the within contrast differences, equal weight per subject
pooled_subj = cell(1, nSubj);
for s = 1:nSubj
    dsub = cell(1, nCond);
    for c = 1:nCond
        if isempty(tfr_low{c, s}) || isempty(tfr_high{c, s})
            continue
        end
        dsub{c} = high_minus_low(tfr_high(c, s), tfr_low(c, s));
    end
    pooled_subj{s} = mean_tfr(dsub);
end
pooled_subj = pooled_subj(~cellfun(@isempty, pooled_subj));
n_pooled = numel(pooled_subj);
if n_pooled < 1
    error('GCP_ms_split_TFR:NoData', ...
        'No subject had a contrast with enough low and high trials.');
end
grand_pooled = mean_tfr(pooled_subj);
fprintf('[MS SPLIT TFR] Pooled grand average N = %d\n', n_pooled);

save(fullfile(fig_dir, 'GCP_ms_split_tfr.mat'), ...
    'subjects', 'condPct', 'min_trials_per_bin', ...
    'tfr_foi', 'tfr_toi', 'tfr_win_sec', 'tfr_tapsmofrq', 'tfr_baseline_window', ...
    'tfr_low', 'tfr_high', 'grand_diff', 'grand_pooled', 'n_grand', 'n_pooled', ...
    'n_low', 'n_high', 'n_tie', 'n_zero', 'n_usable', 'ms_median', ...
    '-v7.3');
fprintf('[MS SPLIT TFR] Saved %s\n', fullfile(fig_dir, 'GCP_ms_split_tfr.mat'));

%% Figures
clim = diff_clim([grand_diff, {grand_pooled}], plot_latency, plot_frequency);
plot_split_diff_tfr(grand_diff, condTitles, n_grand, clim, plot_latency, plot_frequency, ...
    fullfile(fig_dir, 'GCP_split_tfr_diff.png'));
plot_pooled_diff_tfr(grand_pooled, n_pooled, clim, plot_latency, plot_frequency, ...
    fullfile(fig_dir, 'GCP_split_tfr_diff_pooled.png'));

fprintf('[MS SPLIT TFR] Done %s\n', datestr(now));

%% Local functions
function idx = match_subjects(allSubjects, keepSubjects)
allSubjects = cellstr(string(allSubjects(:)));
keepSubjects = cellstr(string(keepSubjects(:)));
idx = nan(1, numel(keepSubjects));
for s = 1:numel(keepSubjects)
    hit = find(strcmp(allSubjects, keepSubjects{s}), 1);
    if isempty(hit)
        error('GCP_ms_split_TFR:SubjectMismatch', ...
            'Included subject %s is missing from GCP_eeg_GED.mat.', keepSubjects{s});
    end
    idx(s) = hit;
end
end

function pr = orient_powratio(pr, nFreq)
pr = double(pr);
if size(pr, 2) ~= nFreq && size(pr, 1) == nFreq
    pr = pr.';
end
if size(pr, 2) ~= nFreq
    error('GCP_ms_split_TFR:FreqMismatch', ...
        'Trial spectra have %d frequencies, scan grid has %d.', size(pr, 2), nFreq);
end
end

function [ms, aligned] = align_gaze_to_trials(G, nTrl)
ms = nan(nTrl, 1);
aligned = false;
if ~isstruct(G) || ~isfield(G, 'MSRate') || ~isfield(G, 'Trial')
    return
end
trial = double(G.Trial(:));
rate = double(G.MSRate(:));
if numel(rate) ~= numel(trial)
    return
end
if numel(trial) == nTrl && isequal(trial(:).', 1:nTrl)
    ms = rate;
    aligned = true;
    return
end
trial = round(trial);
if numel(unique(trial)) == numel(trial) && all(isfinite(trial)) && ...
        all(trial >= 1) && all(trial <= nTrl)
    ms(trial) = rate;
    aligned = true;
end
end

function tfr_avg = ged_bin_tfr(dat, trl_idx, filt, foi, toi, win_sec, tapsmofrq, baseline_window)
cfg_trl = [];
cfg_trl.trials = trl_idx(:)';
dat_bin = ft_selectdata(cfg_trl, dat);

ged_dat = dat_bin;
ged_dat.label = {'GED'};
if isfield(ged_dat, 'elec')
    ged_dat = rmfield(ged_dat, 'elec');
end
if isfield(ged_dat, 'grad')
    ged_dat = rmfield(ged_dat, 'grad');
end
for tr = 1:numel(dat_bin.trial)
    ged_dat.trial{tr} = filt(:).' * double(dat_bin.trial{tr});
end

cfg_tfr = [];
cfg_tfr.method = 'mtmconvol';
cfg_tfr.output = 'pow';
cfg_tfr.taper = 'dpss';
cfg_tfr.foi = foi;
cfg_tfr.toi = toi;
cfg_tfr.t_ftimwin = win_sec * ones(size(foi));
cfg_tfr.tapsmofrq = tapsmofrq * ones(size(foi));
cfg_tfr.keeptrials = 'yes';
cfg_tfr.feedback = 'none';
tfr_trials = ft_freqanalysis(cfg_tfr, ged_dat);

cfg_bl = [];
cfg_bl.baseline = baseline_window;
cfg_bl.baselinetype = 'db';
cfg_bl.parameter = 'powspctrm';
tfr_trials = ft_freqbaseline(cfg_bl, tfr_trials);

cfg_avg = [];
cfg_avg.keeptrials = 'no';
tfr_avg = ft_freqdescriptives(cfg_avg, tfr_trials);
end

function tfr = drop_cfg(tfr)
if isstruct(tfr) && isfield(tfr, 'cfg')
    tfr = rmfield(tfr, 'cfg');
end
end

function tfr_diff = high_minus_low(tfr_high, tfr_low)
ga_high = mean_tfr(tfr_high);
ga_low = mean_tfr(tfr_low);
tfr_diff = ga_high;
tfr_diff.powspctrm = ga_high.powspctrm - ga_low.powspctrm;
end

function tfr = mean_tfr(tfrs)
if iscell(tfrs)
    tfrs = tfrs(~cellfun(@isempty, tfrs));
else
    tfrs = {tfrs};
end
if isempty(tfrs)
    tfr = [];
elseif isscalar(tfrs)
    tfr = tfrs{1};
else
    tfr = ft_freqgrandaverage([], tfrs{:});
    tfr = drop_cfg(tfr);
end
end

function clim = diff_clim(grand, latency, frequency)
max_abs = 0;
for c = 1:numel(grand)
    if isempty(grand{c})
        continue
    end
    freq_idx = grand{c}.freq >= frequency(1) & grand{c}.freq <= frequency(2);
    time_idx = grand{c}.time >= latency(1) & grand{c}.time <= latency(2);
    vals = grand{c}.powspctrm(:, freq_idx, time_idx);
    max_abs = max(max_abs, max(abs(vals(:))));
end
if ~isfinite(max_abs) || max_abs <= 0
    max_abs = 1;
end
clim = double([-max_abs, max_abs]);
end

function color_map = diverging_map()
if exist('cbrewer', 'file') == 2
    color_map = flipud(cbrewer('div', 'RdBu', 64));
else
    color_map = flipud(parula(64));
end
end

function plot_one_diff(tfr, clim, latency, frequency, fontSize)
cfg = [];
cfg.channel = 'GED';
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = latency;
cfg.ylim = frequency;
cfg.shading = 'interp';
cfg.figure = gcf;
ft_singleplotTFR(cfg, tfr);
colormap(diverging_map());
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power difference [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
yticks(frequency(1):10:frequency(2));
xticks([-0.5 0 1 2]);
xline(0, '--');
set(gca, 'FontSize', fontSize);
end

function plot_split_diff_tfr(grand, condTitles, nUsed, clim, latency, frequency, outPath)
fontSize = 22;
fig = figure('Position', [0 0 1512 982], 'Color', 'w');
set(fig, 'PaperPositionMode', 'auto');
for c = 1:numel(grand)
    subplot(2, 2, c);
    if isempty(grand{c})
        axis off
        title(sprintf('%s (no data)', condTitles{c}));
        continue
    end
    plot_one_diff(grand{c}, clim, latency, frequency, fontSize);
    title(sprintf('%s (N = %d)', condTitles{c}, nUsed(c)));
end
sgtitle('High minus low microsaccade rate', 'FontSize', fontSize + 4, 'FontWeight', 'bold');
set(fig, 'Renderer', 'painters');
drawnow;
print(fig, outPath, '-dpng', '-r600');
fprintf('[MS SPLIT TFR] Saved %s\n', outPath);
end

function plot_pooled_diff_tfr(tfr, nUsed, clim, latency, frequency, outPath)
fontSize = 25;
fig = figure('Position', [0 0 1512 982], 'Color', 'w');
set(fig, 'PaperPositionMode', 'auto');
plot_one_diff(tfr, clim, latency, frequency, fontSize);
title(sprintf('High minus low microsaccade rate, contrasts averaged (N = %d)', nUsed));
set(fig, 'Renderer', 'painters');
drawnow;
print(fig, outPath, '-dpng', '-r600');
fprintf('[MS SPLIT TFR] Saved %s\n', outPath);
end
