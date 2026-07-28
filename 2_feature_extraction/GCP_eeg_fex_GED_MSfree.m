%% GCP GED gamma after peri-microsaccade removal (neural vs artifactual)
%
% Keeps the confirmatory GED spatial filter from GCP_eeg_GED.mat fixed.
% Recomputes trial-level peak frequency / peak power after excluding samples
% within +/- ms_guard_s of each Engbert microsaccade onset.
%
% Requires:
%   data/features/GCP_eeg_GED.mat
%   <subj>/eeg/dataEEG.mat
%   <subj>/gaze/gaze_microsaccade_events.mat  (from GCP_gaze_fex.m)
%   Fallback if events missing: redetect Engbert MS from dataET.mat
%
% Outputs:
%   data/features/GCP_eeg_GED_MSfree.mat

%% Setup
startup
[subjects, paths] = setup('GCP', 0);
nSubj = numel(subjects);

ged_path = fullfile(paths.features, 'GCP_eeg_GED.mat');
if ~isfile(ged_path)
    error('GCP_eeg_fex_GED_MSfree:NoGED', 'Missing %s. Run GCP_eeg_fex_GED.m first.', ged_path);
end
GED = load(ged_path, 'all_combined_filter_full', 'scan_freqs', 'subjects', ...
    'baseline_window', 'full_window');

if isfield(GED, 'subjects') && ~isempty(GED.subjects)
    subjects = GED.subjects;
    nSubj = numel(subjects);
end

scan_freqs = GED.scan_freqs(:)';
baseline_window = GED.baseline_window;
full_window = GED.full_window;
if isempty(baseline_window), baseline_window = [-1.5, -0.5]; end
if isempty(full_window), full_window = [0, 2.0]; end

% Analysis parameters (match confirmatory GED where relevant)
ms_guard_s = 0.10;                 % exclude +/- 100 ms around MS onset
min_keep_frac = 0.50;              % drop trial if fewer samples retained
min_seg_samp = 8;                  % minimum samples for a kept chunk
mtmfft_tapsmofrq_hz = 3;
powratio_smooth_bins = 5;
peak_power_halfwidth_hz = 5;
analysis_freq_range = [30 90];
condNames = {'c25', 'c50', 'c75', 'c100'};

trials_peaks_MSfree = cell(4, nSubj);
trials_gamma_power_MSfree = cell(4, nSubj);
trials_powratio_MSfree = cell(4, nSubj);
trials_pct_kept_MSfree = cell(4, nSubj);
trials_n_ms_MSfree = cell(4, nSubj);

fprintf('[GED-MSfree] guard=+/-%.0f ms | min keep=%.0f%% | subjects=%d\n', ...
    1000 * ms_guard_s, 100 * min_keep_frac, nSubj);

for si = 1:nSubj
    subj = subjects{si};
    fprintf('[GED-MSfree] Subject %s (%d/%d)\n', subj, si, nSubj);

    w = [];
    if isfield(GED, 'all_combined_filter_full') && numel(GED.all_combined_filter_full) >= si
        w = GED.all_combined_filter_full{si};
    end
    if isempty(w) || ~all(isfinite(w(:)))
        warning('GCP_eeg_fex_GED_MSfree:NoFilter', ...
            'No valid combined GED filter for subject %s. Skipping.', subj);
        continue
    end
    w = w(:);

    eeg_path = fullfile(paths.features, subj, 'eeg', 'dataEEG.mat');
    if ~isfile(eeg_path)
        warning('GCP_eeg_fex_GED_MSfree:NoEEG', 'Missing %s', eeg_path);
        continue
    end
    E = load(eeg_path, 'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    data_eeg = {E.dataEEG_c25, E.dataEEG_c50, E.dataEEG_c75, E.dataEEG_c100};

    ms_events = load_ms_events(paths.features, subj, condNames);

    for ci = 1:4
        dat = data_eeg{ci};
        if isempty(dat) || ~isfield(dat, 'trial')
            continue
        end
        nTrl = numel(dat.trial);
        fs = dat.fsample;
        if size(dat.trial{1}, 1) ~= numel(w)
            warning('GCP_eeg_fex_GED_MSfree:ChanMismatch', ...
                'Subject %s cond %s: filter length %d vs channels %d. Skipping.', ...
                subj, condNames{ci}, numel(w), size(dat.trial{1}, 1));
            continue
        end

        peaks = nan(nTrl, 1);
        powers = nan(nTrl, 1);
        pct_kept = nan(nTrl, 1);
        n_ms = zeros(nTrl, 1);
        powratio = nan(nTrl, numel(scan_freqs));

        onsets_cond = {};
        if ~isempty(ms_events) && numel(ms_events) >= ci && isfield(ms_events{ci}, 'Onset')
            onsets_cond = ms_events{ci}.Onset;
        end

        for trl = 1:nTrl
            x = double(dat.trial{trl});
            t = dat.time{trl};
            z = (w' * x);
            z = z(:).';
            t = t(:).';

            idx_base = t >= baseline_window(1) & t <= baseline_window(2);
            idx_stim = t >= full_window(1) & t <= full_window(2);
            if ~any(idx_base) || ~any(idx_stim)
                continue
            end

            onsets = [];
            if iscell(onsets_cond) && trl <= numel(onsets_cond) && ~isempty(onsets_cond{trl})
                onsets = onsets_cond{trl}(:);
            end
            onsets = onsets(isfinite(onsets));
            n_ms(trl) = numel(onsets);

            keep = true(size(t));
            for oi = 1:numel(onsets)
                keep = keep & ~(t >= onsets(oi) - ms_guard_s & t <= onsets(oi) + ms_guard_s);
            end

            keep_stim = keep & idx_stim;
            keep_base = keep & idx_base;
            pct_kept(trl) = sum(keep_stim) / max(sum(idx_stim), 1);
            if pct_kept(trl) < min_keep_frac
                continue
            end

            sig_stim = concatenate_kept_chunks(z(keep_stim), min_seg_samp);
            sig_base = concatenate_kept_chunks(z(keep_base), min_seg_samp);
            if isempty(sig_stim) || isempty(sig_base)
                continue
            end

            [p_stim, p_base] = compute_scan_power_pair(sig_stim, sig_base, fs, scan_freqs, mtmfft_tapsmofrq_hz);
            if isempty(p_stim) || isempty(p_base)
                continue
            end
            valid = isfinite(p_stim) & isfinite(p_base) & (p_stim > 0) & (p_base > 0);
            pr = nan(1, numel(scan_freqs));
            pr(valid) = 10 * log10(p_stim(valid) ./ p_base(valid));
            powratio(trl, :) = pr;

            band = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
            [peak_hz, peak_pow] = pick_tallest_peak(pr(band), scan_freqs(band), ...
                powratio_smooth_bins, peak_power_halfwidth_hz);
            peaks(trl) = peak_hz;
            powers(trl) = peak_pow;
        end

        trials_peaks_MSfree{ci, si} = peaks;
        trials_gamma_power_MSfree{ci, si} = powers;
        trials_powratio_MSfree{ci, si} = powratio;
        trials_pct_kept_MSfree{ci, si} = pct_kept;
        trials_n_ms_MSfree{ci, si} = n_ms;
    end
end

% Subject-level condition means (for quick H5/H6 checks)
subj_mean_peak_MSfree = nan(4, nSubj);
subj_mean_power_MSfree = nan(4, nSubj);
subj_mean_pct_kept = nan(4, nSubj);
for si = 1:nSubj
    for ci = 1:4
        pf = trials_peaks_MSfree{ci, si};
        pp = trials_gamma_power_MSfree{ci, si};
        pk = trials_pct_kept_MSfree{ci, si};
        if ~isempty(pf)
            subj_mean_peak_MSfree(ci, si) = mean(pf, 'omitnan');
        end
        if ~isempty(pp)
            subj_mean_power_MSfree(ci, si) = mean(pp, 'omitnan');
        end
        if ~isempty(pk)
            subj_mean_pct_kept(ci, si) = mean(pk, 'omitnan');
        end
    end
end

save_path = fullfile(paths.features, 'GCP_eeg_GED_MSfree.mat');
save(save_path, ...
    'trials_peaks_MSfree', 'trials_gamma_power_MSfree', 'trials_powratio_MSfree', ...
    'trials_pct_kept_MSfree', 'trials_n_ms_MSfree', ...
    'subj_mean_peak_MSfree', 'subj_mean_power_MSfree', 'subj_mean_pct_kept', ...
    'scan_freqs', 'subjects', 'condNames', ...
    'ms_guard_s', 'min_keep_frac', 'baseline_window', 'full_window', ...
    'mtmfft_tapsmofrq_hz', 'peak_power_halfwidth_hz', '-v7.3');
fprintf('[GED-MSfree] Saved %s\n', save_path);

%% Local helpers
function ms_events = load_ms_events(features_root, subj, condNames)
ms_events = cell(1, 4);
ev_path = fullfile(features_root, subj, 'gaze', 'gaze_microsaccade_events.mat');
if isfile(ev_path)
    S = load(ev_path);
    for ci = 1:4
        fn = sprintf('ms_events_%s', condNames{ci});
        if isfield(S, fn)
            ms_events{ci} = S.(fn);
        end
    end
    return
end

% Fallback: Engbert redetect from dataET (same detector as gaze_fex)
et_path = fullfile(features_root, subj, 'gaze', 'dataET.mat');
if ~isfile(et_path)
    warning('GCP_eeg_fex_GED_MSfree:NoMS', ...
        'No MS events or dataET for %s; MS-free masks will be empty.', subj);
    return
end
fprintf('  [fallback] Detecting Engbert MS from dataET for %s\n', subj);
G = load(et_path);
data_et = {G.dataET_c25, G.dataET_c50, G.dataET_c75, G.dataET_c100};
for ci = 1:4
    gaze = data_et{ci};
    nTrl = numel(gaze.trial);
    onsets = cell(1, nTrl);
    offsets = cell(1, nTrl);
    for trl = 1:nTrl
        [on_s, off_s] = detect_engbert_ms_times(gaze.trial{trl}, gaze.time{trl}, gaze.fsample);
        onsets{trl} = on_s;
        offsets{trl} = off_s;
    end
    ms_events{ci} = struct();
    ms_events{ci}.Onset = onsets;
    ms_events{ci}.Offset = offsets;
end
end

function [onset_s, offset_s] = detect_engbert_ms_times(raw, tVec, fsample)
onset_s = [];
offset_s = [];
if isempty(raw) || isempty(tVec)
    return
end
win = [-1.5, 2.0];
full_idx = tVec >= win(1) & tVec <= win(2);
t_full = tVec(full_idx);
x = raw(1, full_idx);
y = 600 - raw(2, full_idx);
valid = x >= 0 & x <= 800 & y >= 0 & y <= 600 & isfinite(x) & isfinite(y);
if sum(valid) < round(0.5 * fsample)
    return
end
x_val = x(valid);
y_val = y(valid);
try
    tmp = [x_val; y_val; nan(1, numel(x_val))];
    tmp = remove_blinks(tmp, round(0.2 * fsample));
    x_val = tmp(1, :);
    y_val = tmp(2, :);
catch
end
valid_clean = isfinite(x_val) & isfinite(y_val);
if sum(valid_clean) < round(0.5 * fsample)
    return
end
x_clean = x_val(valid_clean);
y_clean = y_val(valid_clean);
[~, ms_det] = detect_microsaccades(fsample, [x_clean; y_clean], numel(x_clean));
if isempty(ms_det.Onset)
    return
end
idx_full_valid = find(valid);
idx_clean_in_full = idx_full_valid(valid_clean);
nEv = min(numel(ms_det.Onset), numel(ms_det.Offset));
on_idx = ms_det.Onset(1:nEv);
off_idx = ms_det.Offset(1:nEv);
keep = on_idx >= 1 & on_idx <= numel(idx_clean_in_full) & ...
       off_idx >= 1 & off_idx <= numel(idx_clean_in_full);
onset_s = t_full(idx_clean_in_full(on_idx(keep)))';
offset_s = t_full(idx_clean_in_full(off_idx(keep)))';
end

function sig = concatenate_kept_chunks(z, min_seg_samp)
sig = [];
z = z(:)';
if isempty(z)
    return
end
% kept samples are already contiguous in the logical mask extract; still
% split if zeros-length gaps were removed (they are already contiguous)
if numel(z) < min_seg_samp
    return
end
sig = z.';
end

function [p_stim, p_base] = compute_scan_power_pair(sig_stim, sig_base, fs, scan_freqs, tapsmofrq_hz)
p_stim = nan(1, numel(scan_freqs));
p_base = nan(1, numel(scan_freqs));
p_all = compute_scan_power_mtmfft({sig_stim, sig_base}, fs, scan_freqs, tapsmofrq_hz);
if size(p_all, 1) < 2
    return
end
p_stim = p_all(1, :);
p_base = p_all(2, :);
end

function p_scan = compute_scan_power_mtmfft(sig_cell, fs, scan_freqs, tapsmofrq_hz)
nSig = numel(sig_cell);
p_scan = nan(nSig, numel(scan_freqs));
dat = [];
dat.label = {'GED'};
dat.fsample = fs;
dat.trial = {};
dat.time = {};
row_map = [];
for si = 1:nSig
    x = double(sig_cell{si});
    if iscolumn(x), x = x.'; end
    if numel(x) < 8 || any(~isfinite(x))
        continue
    end
    x = x - mean(x);
    dat.trial{end+1} = x;
    dat.time{end+1} = (0:(numel(x)-1)) / fs;
    row_map(end+1) = si; %#ok<AGROW>
end
if isempty(dat.trial)
    return
end
nTrials_valid = numel(dat.trial);
dat.sampleinfo = zeros(nTrials_valid, 2);
sample_start = 1;
for ti = 1:nTrials_valid
    nSamp = numel(dat.trial{ti});
    dat.sampleinfo(ti, :) = [sample_start, sample_start + nSamp - 1];
    sample_start = sample_start + nSamp;
end

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'dpss';
cfg.foi = scan_freqs;
cfg.tapsmofrq = tapsmofrq_hz;
cfg.pad = 'nextpow2';
cfg.keeptrials = 'yes';
cfg.feedback = 'none';
try
    freq = ft_freqanalysis(cfg, dat);
catch
    return
end
pow = freq.powspctrm;
if ndims(pow) == 3
    pow = squeeze(pow(:, 1, :));
elseif isvector(pow)
    pow = pow(:).';
end
if size(pow, 1) ~= numel(row_map)
    return
end
for ri = 1:numel(row_map)
    p_scan(row_map(ri), :) = pow(ri, :);
end
end

function [peak_hz, peak_power] = pick_tallest_peak(y, x, smooth_n, peak_power_halfwidth_hz)
peak_hz = NaN;
peak_power = NaN;
y = y(:);
x = x(:);
if ~isfinite(smooth_n) || smooth_n < 1
    smooth_n = 1;
end
if isempty(y) || numel(y) ~= numel(x)
    return
end
y = movmean(y, max(1, round(smooth_n)), 'omitnan');
valid = isfinite(y) & isfinite(x);
if sum(valid) < 3
    return
end
x_use = x(valid);
y_use = y(valid);
dx = diff(x_use);
dx = dx(isfinite(dx) & dx > 0);
if isempty(dx)
    return
end
freq_step = median(dx);
min_peak_width_hz = max(2 * freq_step, 2.0);
local_spread = iqr(y_use);
if ~isfinite(local_spread) || local_spread <= 0
    local_spread = std(y_use, 'omitnan');
end
if ~isfinite(local_spread) || local_spread <= 0
    local_spread = 1;
end
min_peak_prom_db = max(0.15, 0.25 * local_spread);
try
    [pks, locs] = findpeaks(y_use, x_use, ...
        'MinPeakProminence', min_peak_prom_db, ...
        'MinPeakWidth', min_peak_width_hz, ...
        'SortStr', 'descend', ...
        'NPeaks', 1);
catch
    return
end
if isempty(pks) || isempty(locs)
    return
end
peak_hz = locs(1);
peak_power = pks(1);
if peak_power_halfwidth_hz > 0
    band_mask = abs(x_use - peak_hz) <= peak_power_halfwidth_hz;
    band_power = y_use(band_mask);
    band_power = band_power(isfinite(band_power));
    if ~isempty(band_power)
        peak_power = mean(band_power);
    end
end
end
