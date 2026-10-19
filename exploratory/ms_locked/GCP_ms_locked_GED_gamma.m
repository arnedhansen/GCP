%% GCP microsaccade GED versus the existing gamma filter
%
% Exploratory. Not confirmatory. Does not overwrite GCP_eeg_GED.mat.
%
% Question
%   Which spatial pattern separates the moment of a microsaccade from a
%   window in the same trial with no microsaccade? The pattern may be
%   broadband or spiky. It is not required to be occipital gamma.
%   The saved gamma filter is still scored on the same snippets.
%
% One filter per subject, contrasts pooled
%   Signal:    EEG from -0.15 to +0.15 s around an Engbert onset
%   Reference: a window of the same length in the same trial, inside 0 to 2 s,
%              containing no microsaccade, and as close in time to that onset
%              as possible. Time matching stops the filter from picking up an
%              early versus late stimulus effect.
%   Events closer than the window length are dropped, so one snippet does not
%   contain two microsaccades. Snippets from all four contrasts are
%   concatenated, so each event has equal weight. The fit needs at least
%   min_pairs pairs in that pool.
%   Spectra and the gamma filter ratio are then scored inside each contrast
%   with that one filter.
%
% The GED is fit on 30 to 90 Hz covariance, with the same regularisation
% and search depth as GCP_eeg_fex_GED.m. The confirmatory occipital and
% peak form gates are not used. The leading component is kept when its
% eigenvalue is at least 1.1 and no single electrode exceeds 25% of the
% map. Polarity makes the strongest electrode positive. EMG (rising high
% frequency slope) is plotted and kept out of the group spectrum.
%
% Comparison, same peri and control snippets
%   Rayleigh quotient of the saved gamma filter (power ratio in dB).
%   Topography correlation between a non EMG leading component and the
%   saved gamma topography.
%   Peri over control spectrum of the gamma filter and of the non EMG
%   leading components.
%
% The window length is a choice. Shorter windows favour the ocular spike.
% Longer windows mix in grating gamma that is not locked to the movement.

%% Setup
startup
prev_vis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
cleanup_vis = onCleanup(@() set(0, 'DefaultFigureVisible', prev_vis)); %#ok<NASGU>

[subjects, paths, colors, headmodel] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
if ispc
    addpath('C:\Users\Administrator\Documents\GitHub\GCP\2_feature_extraction\GED-helpers');
else
    addpath('/Users/Arne/Documents/GitHub/GCP/2_feature_extraction/GED-helpers');
end
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar');

nSubj = numel(subjects);
nCond = 4;
condPct = [25 50 75 100];
condTitles = {'25% Contrast', '50% Contrast', '75% Contrast', '100% Contrast'};
eegFields = {'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100'};
evFields = {'ms_events_c25', 'ms_events_c50', 'ms_events_c75', 'ms_events_c100'};
if size(colors, 1) < nCond
    colors = color_def('GCP');
end

% Same GED settings as GCP_eeg_fex_GED.m
analysis_freq_range = [30 90];
scan_freqs = analysis_freq_range(1):1:analysis_freq_range(2);
nFreq = numel(scan_freqs);
lambda = 0.05;
ged_search_n = 10;
min_eigval = 1.1;
mtmfft_tapsmofrq_hz = 3;
stim_window = [0 2];
peri_half_s = 0.15;
control_step_s = 0.02;
min_pairs = 30;

params = struct( ...
    'lambda', lambda, ...
    'ged_search_n', ged_search_n, ...
    'min_eigval', min_eigval, ...
    'tapsmo', mtmfft_tapsmofrq_hz, ...
    'analysis_freq_range', analysis_freq_range, ...
    'scan_freqs', scan_freqs);

fig_dir = fullfile(paths.figures, 'exploratory', 'ms_locked');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

fprintf('[MS GED] Included cohort N = %d (%s)\n', nSubj, strjoin(subjects, ', '));

%% Existing gamma filter (not refit, not replaced)
ged = load(fullfile(paths.features, 'GCP_eeg_GED.mat'), ...
    'all_combined_filter', 'all_topos', 'all_topo_labels', 'subjects');

n_events = nan(nCond, nSubj);
n_pairs = nan(nCond, nSubj);
median_abs_dt = nan(nCond, nSubj);
n_pairs_pooled = nan(1, nSubj);
passed = false(1, nSubj);
topo_r = nan(1, nSubj);
rayleigh_db = nan(nCond, nSubj);
eig_selected = nan(1, nSubj);
top_eig = nan(1, nSubj);
top_pf = nan(1, nSubj);
top_post_front = nan(1, nSubj);
top_post_temp = nan(1, nSubj);
top_hf_slope = nan(1, nSubj);
top_peak_frac = nan(1, nSubj);
top_class = repmat({''}, 1, nSubj);
spec_ms = nan(nCond, nSubj, nFreq);
spec_gamma = nan(nCond, nSubj, nFreq);
topos_ms = cell(1, nSubj);
topos_gamma = cell(1, nSubj);
topo_labels = cell(1, nSubj);

%% Subject loop
for s = 1:nSubj
    gi = find(strcmp(ged.subjects, subjects{s}), 1);
    if isempty(gi)
        fprintf('[MS GED] %s: not in GCP_eeg_GED.mat, skipped\n', subjects{s});
        continue
    end
    w_gamma = ged.all_combined_filter{gi};
    if isempty(w_gamma) || ~any(isfinite(w_gamma))
        fprintf('[MS GED] %s: empty gamma filter, skipped\n', subjects{s});
        continue
    end
    w_gamma = w_gamma(:);
    ged_labels = ged.all_topo_labels{gi};
    if numel(w_gamma) ~= numel(ged_labels)
        fprintf('[MS GED] %s: filter length does not match labels, skipped\n', subjects{s});
        continue
    end
    if ~isempty(ged.all_topos{gi})
        topos_gamma{s} = ged.all_topos{gi}(:);
    end
    topo_labels{s} = ged_labels;

    eeg_file = fullfile(paths.features, subjects{s}, 'eeg', 'dataEEG.mat');
    ev_file = fullfile(paths.features, subjects{s}, 'gaze', 'gaze_microsaccade_events.mat');
    if ~isfile(eeg_file) || ~isfile(ev_file)
        fprintf('[MS GED] %s: missing EEG or microsaccade events, skipped\n', subjects{s});
        continue
    end
    eeg = load(eeg_file, eegFields{:});
    ev = load(ev_file, evFields{:});
    pair_store = cell(1, nCond);
    fs_here = NaN;
    for c = 1:nCond
        if ~isfield(eeg, eegFields{c}) || isempty(eeg.(eegFields{c}))
            fprintf('[MS GED] %s %d%%: empty EEG, skipped\n', subjects{s}, condPct(c));
            continue
        end
        dat = order_channels(eeg.(eegFields{c}), ged_labels);
        fs_here = dat.fsample;
        if ~isfield(ev, evFields{c}) || ~isfield(ev.(evFields{c}), 'Onset')
            fprintf('[MS GED] %s %d%%: no onset field, skipped\n', subjects{s}, condPct(c));
            continue
        end
        onsets = ev.(evFields{c}).Onset;
        if ~iscell(onsets) || numel(onsets) ~= numel(dat.trial)
            fprintf('[MS GED] %s %d%%: onset count %d does not match EEG trials %d, skipped\n', ...
                subjects{s}, condPct(c), numel(onsets), numel(dat.trial));
            continue
        end

        cfg_filt = [];
        cfg_filt.bpfilter = 'yes';
        cfg_filt.bpfreq = analysis_freq_range;
        cfg_filt.bpfilttype = 'fir';
        cfg_filt.bpfiltord = round(3 * dat.fsample / analysis_freq_range(1));
        cfg_filt.feedback = 'none';
        dat_gamma = ft_preprocessing(cfg_filt, dat);

        pairs = pair_windows(dat, dat_gamma, onsets, stim_window, peri_half_s, control_step_s);
        n_events(c, s) = pairs.n_events;
        n_pairs(c, s) = pairs.n_pairs;
        median_abs_dt(c, s) = pairs.median_abs_dt;
        if pairs.n_pairs < 1
            fprintf('[MS GED] %s %d%%: no pairs, events in 0 to 2 s = %d\n', ...
                subjects{s}, condPct(c), pairs.n_events);
            continue
        end
        pair_store{c} = pairs;

        t_snip = pairs.time_snip;
        datS_c = snippets_to_ft(pairs.peri_gam, ged_labels, dat.fsample, t_snip);
        datR_c = snippets_to_ft(pairs.ctrl_gam, ged_labels, dat.fsample, t_snip);
        win_c = [t_snip(1) - 1e-6, t_snip(end) + 1e-6];
        covS_c = compute_pooled_covariance_window(datS_c, win_c);
        covR_c = compute_pooled_covariance_window(datR_c, win_c);
        covS_c = (covS_c + covS_c') / 2;
        covR_c = (covR_c + covR_c') / 2;
        num = real(w_gamma' * covS_c * w_gamma);
        den = real(w_gamma' * covR_c * w_gamma);
        if isfinite(num) && isfinite(den) && den > 0
            rayleigh_db(c, s) = 10 * log10(num / den);
        end
    end

    use_c = find(~cellfun(@isempty, pair_store));
    n_pairs_pooled(s) = sum(n_pairs(use_c, s), 'omitnan');
    if isempty(use_c) || n_pairs_pooled(s) < min_pairs
        fprintf('[MS GED] %s: %d pooled pairs (need %d), no microsaccade filter\n', ...
            subjects{s}, n_pairs_pooled(s), min_pairs);
        for c = use_c
            fprintf('[MS GED] %s %d%%: pairs=%d gamma ratio=%.2f dB\n', ...
                subjects{s}, condPct(c), n_pairs(c, s), rayleigh_db(c, s));
        end
        continue
    end

    datS_parts = cell(1, numel(use_c));
    datR_parts = cell(1, numel(use_c));
    paired_parts = cell(1, numel(use_c));
    for ic = 1:numel(use_c)
        pairs = pair_store{use_c(ic)};
        t_snip = pairs.time_snip;
        datS_parts{ic} = snippets_to_ft(pairs.peri_gam, ged_labels, fs_here, t_snip);
        datR_parts{ic} = snippets_to_ft(pairs.ctrl_gam, ged_labels, fs_here, t_snip);
        paired_parts{ic} = paired_ft(pairs.ctrl_raw, pairs.peri_raw, ged_labels, ...
            fs_here, pairs.time_ctrl, pairs.time_peri);
    end
    datS = append_ft(datS_parts);
    datR = append_ft(datR_parts);
    paired_pool = append_ft(paired_parts);
    t_snip = pair_store{use_c(1)}.time_snip;
    win = [t_snip(1) - 1e-6, t_snip(end) + 1e-6];
    covS = compute_pooled_covariance_window(datS, win);
    covR = compute_pooled_covariance_window(datR, win);
    covS = (covS + covS') / 2;
    covR = (covR + covR') / 2;
    peri_win = [pair_store{use_c(1)}.time_peri(1) - 1e-6, pair_store{use_c(1)}.time_peri(end) + 1e-6];
    ctrl_win = [pair_store{use_c(1)}.time_ctrl(1) - 1e-6, pair_store{use_c(1)}.time_ctrl(end) + 1e-6];

    sel = select_ms_components(covS, covR, paired_pool, ged_labels, ...
        fs_here, peri_win, ctrl_win, params);
    top_eig(s) = sel.top_eig;
    top_pf(s) = sel.top_pf;
    top_post_front(s) = sel.top_post_front;
    top_post_temp(s) = sel.top_post_temp;
    top_hf_slope(s) = sel.top_hf_slope;
    top_peak_frac(s) = sel.top_peak_frac;
    top_class{s} = sel.top_class;
    passed(s) = sel.in_spectrum;
    if sel.in_spectrum
        eig_selected(s) = sel.eig;
    end
    if sel.show && ~isempty(sel.w)
        topos_ms{s} = sel.topo(:);
        if ~isempty(topos_gamma{s}) && numel(topos_gamma{s}) == numel(sel.topo) && sel.in_spectrum
            topo_r(s) = corr(sel.topo(:), topos_gamma{s}(:), 'rows', 'complete');
        end
    end
    dt_all = [];
    for ic = 1:numel(use_c)
        dt_all = [dt_all, pair_store{use_c(ic)}.abs_dt]; %#ok<AGROW>
    end
    if isempty(dt_all)
        dt_pool = NaN;
    else
        dt_pool = median(dt_all);
    end
    fprintf(['[MS GED] %s: pooled pairs=%d median |dt|=%.0f ms | ', ...
        'top eig=%.2f (%s, PF=%.2f, peak frac=%.2f) map=%d spectrum=%d | topo r=%.2f\n'], ...
        subjects{s}, n_pairs_pooled(s), 1000 * dt_pool, sel.top_eig, sel.top_class, ...
        sel.top_pf, sel.top_peak_frac, sel.show, sel.in_spectrum, topo_r(s));

    for c = use_c
        pairs = pair_store{c};
        paired = paired_ft(pairs.ctrl_raw, pairs.peri_raw, ged_labels, fs_here, ...
            pairs.time_ctrl, pairs.time_peri);
        peri_c = [pairs.time_peri(1) - 1e-6, pairs.time_peri(end) + 1e-6];
        ctrl_c = [pairs.time_ctrl(1) - 1e-6, pairs.time_ctrl(end) + 1e-6];
        Wplot = w_gamma;
        plot_which = {'gamma'};
        if sel.in_spectrum && ~isempty(sel.w)
            Wplot = [sel.w(:), w_gamma]; %#ok<AGROW>
            plot_which = {'ms', 'gamma'};
        end
        prox = estimate_components_artifact_proxies( ...
            Wplot, {paired}, peri_c, ctrl_c, fs_here, scan_freqs, ...
            mtmfft_tapsmofrq_hz, '');
        for ip = 1:numel(plot_which)
            sp = prox(ip).mean_pr_spectrum(:)';
            if numel(sp) ~= nFreq
                continue
            end
            if strcmp(plot_which{ip}, 'ms')
                spec_ms(c, s, :) = sp;
            else
                spec_gamma(c, s, :) = sp;
            end
        end
        fprintf('[MS GED] %s %d%%: pairs=%d gamma ratio=%.2f dB\n', ...
            subjects{s}, condPct(c), n_pairs(c, s), rayleigh_db(c, s));
    end
    clear pair_store datS datR paired_pool
end

out_mat = fullfile(paths.features, 'GCP_ms_locked_GED_gamma.mat');
save(out_mat, ...
    'subjects', 'scan_freqs', 'condPct', 'peri_half_s', 'min_pairs', ...
    'n_events', 'n_pairs', 'n_pairs_pooled', 'median_abs_dt', 'passed', ...
    'topo_r', 'rayleigh_db', 'eig_selected', ...
    'top_eig', 'top_pf', 'top_post_front', 'top_post_temp', ...
    'top_hf_slope', 'top_peak_frac', 'top_class', ...
    'spec_ms', 'spec_gamma', 'topos_ms', 'topos_gamma', 'topo_labels', '-v7.3');
fprintf('[MS GED] Saved %s\n', out_mat);

%% Figures
plot_ms_ged_spectra(spec_ms, spec_gamma, scan_freqs, colors, condTitles, ...
    fullfile(fig_dir, 'GCP_ms_ged_spectra.png'));
plot_ms_ged_metrics(topo_r, rayleigh_db, colors, condPct, ...
    fullfile(fig_dir, 'GCP_ms_ged_metrics.png'));
if isstruct(headmodel) && isfield(headmodel, 'layANThead')
    plot_ms_ged_topos(topos_ms, topos_gamma, topo_labels, topo_r, top_class, ...
        subjects, headmodel.layANThead, fig_dir);
end
fprintf('[MS GED] Figures in %s\n', fig_dir);

%% Local functions

function dat = order_channels(dat, labels)
cfg = [];
cfg.channel = labels;
cfg.feedback = 'no';
dat = ft_selectdata(cfg, dat);
if isequal(dat.label(:), labels(:))
    return
end
[ok, idx] = ismember(labels, dat.label);
if ~all(ok)
    missing = labels(~ok);
    error('GCP_ms_locked_GED_gamma:Channel', ...
        'EEG is missing channels needed by the gamma filter: %s', strjoin(missing, ', '));
end
for tr = 1:numel(dat.trial)
    dat.trial{tr} = dat.trial{tr}(idx, :);
end
dat.label = labels;
end

function pairs = pair_windows(dat_raw, dat_gam, onsets, stim_window, peri_half, step_s)
% One peri snippet and one time matched control snippet per isolated event.
fs = dat_raw.fsample;
nSamp = 2 * round(peri_half * fs) + 1;
peri_raw = {};
ctrl_raw = {};
peri_gam = {};
ctrl_gam = {};
abs_dt = [];
n_events = 0;
nTr = numel(dat_raw.trial);
for tr = 1:nTr
    on = double(onsets{tr}(:));
    on = on(isfinite(on));
    n_events = n_events + sum(on >= stim_window(1) & on <= stim_window(2));
    if isempty(on) || numel(dat_raw.time{tr}) < nSamp
        continue
    end
    t = dat_raw.time{tr}(:)';
    xraw = double(dat_raw.trial{tr});
    xgam = double(dat_gam.trial{tr});
    if size(xgam, 2) ~= size(xraw, 2)
        continue
    end
    keep = false(numel(on), 1);
    for iEv = 1:numel(on)
        lo = on(iEv) - peri_half;
        hi = on(iEv) + peri_half;
        in_stim = lo >= stim_window(1) && hi <= stim_window(2);
        n_inside = sum(abs(on - on(iEv)) <= peri_half);
        keep(iEv) = in_stim && n_inside == 1;
    end
    use = on(keep);
    if isempty(use)
        continue
    end
    use = sort(use);
    centers = (stim_window(1) + peri_half):step_s:(stim_window(2) - peri_half);
    cand = centers(:)';
    keep_c = false(size(cand));
    for k = 1:numel(cand)
        clo = cand(k) - peri_half;
        chi = cand(k) + peri_half;
        hit = false;
        for iOn = 1:numel(on)
            if clo < (on(iOn) + peri_half) && (on(iOn) - peri_half) < chi
                hit = true;
                break
            end
        end
        keep_c(k) = ~hit;
    end
    cand = cand(keep_c);
    % Closest match first, then drop every candidate that overlaps that control.
    remaining = use(:)';
    while ~isempty(remaining) && ~isempty(cand)
        best_dt = inf;
        best_ev = NaN;
        best_c = NaN;
        for iEv = 1:numel(remaining)
            [dt, ix] = min(abs(cand - remaining(iEv)));
            if dt < best_dt
                best_dt = dt;
                best_ev = iEv;
                best_c = cand(ix);
            end
        end
        if ~isfinite(best_c)
            break
        end
        idx_p = center_idx(t, remaining(best_ev), nSamp);
        idx_c = center_idx(t, best_c, nSamp);
        clo = cand - peri_half;
        chi = cand + peri_half;
        overlap = clo < (best_c + peri_half) & (best_c - peri_half) < chi;
        cand = cand(~overlap);
        if isempty(idx_p)
            remaining(best_ev) = [];
            continue
        end
        if isempty(idx_c)
            continue
        end
        remaining(best_ev) = [];
        peri_raw{end + 1} = xraw(:, idx_p); %#ok<AGROW>
        peri_gam{end + 1} = xgam(:, idx_p); %#ok<AGROW>
        ctrl_raw{end + 1} = xraw(:, idx_c); %#ok<AGROW>
        ctrl_gam{end + 1} = xgam(:, idx_c); %#ok<AGROW>
        abs_dt(end + 1) = best_dt; %#ok<AGROW>
    end
end
pairs.peri_raw = peri_raw;
pairs.ctrl_raw = ctrl_raw;
pairs.peri_gam = peri_gam;
pairs.ctrl_gam = ctrl_gam;
pairs.n_events = n_events;
pairs.n_pairs = numel(peri_raw);
pairs.abs_dt = abs_dt;
if isempty(abs_dt)
    pairs.median_abs_dt = NaN;
else
    pairs.median_abs_dt = median(abs_dt);
end
pairs.time_snip = (0:nSamp - 1) / fs;
pairs.time_ctrl = pairs.time_snip;
pairs.time_peri = 1 + pairs.time_snip;
end

function idx = center_idx(t, center, nSamp)
idx = [];
dt = median(diff(t));
if ~isfinite(dt) || dt <= 0 || ~isfinite(center)
    return
end
[dmin, ic] = min(abs(t - center));
if dmin > max(dt * 1.5, 0)
    return
end
halfL = floor((nSamp - 1) / 2);
i1 = ic - halfL;
i2 = i1 + nSamp - 1;
if i1 < 1 || i2 > numel(t)
    return
end
idx = i1:i2;
end

function dat = append_ft(parts)
parts = parts(~cellfun(@isempty, parts));
if isempty(parts)
    dat = [];
    return
end
if numel(parts) == 1
    dat = parts{1};
    return
end
cfg = [];
cfg.feedback = 'no';
dat = ft_appenddata(cfg, parts{:});
end

function dat = snippets_to_ft(snips, labels, fs, time_vec)
dat = [];
dat.label = labels;
dat.fsample = fs;
n = numel(snips);
dat.trial = cell(1, n);
dat.time = cell(1, n);
for i = 1:n
    dat.trial{i} = snips{i};
    dat.time{i} = time_vec;
end
end

function dat = paired_ft(ctrl, peri, labels, fs, time_ctrl, time_peri)
dat = [];
dat.label = labels;
dat.fsample = fs;
n = numel(ctrl);
dat.trial = cell(1, n);
dat.time = cell(1, n);
tcat = [time_ctrl, time_peri];
for i = 1:n
    dat.trial{i} = [ctrl{i}, peri{i}];
    dat.time{i} = tcat;
end
end

function sel = select_ms_components(covS, covR, paired, labels, fs, peri_win, ctrl_win, params)
% Leading component only. Eigenvalue and single electrode fraction.
% Peak form and occipital location are recorded and not used as gates.
% EMG is returned for plotting and marked out of the group spectrum.
nChan = size(covS, 1);
sel = empty_selection(numel(params.scan_freqs));
if nChan < 2 || ~all(isfinite(covS(:))) || ~all(isfinite(covR(:)))
    return
end
if mean(diag(covS)) <= 0 || mean(diag(covR)) <= 0
    return
end
lam = params.lambda;
eyeN = eye(nChan);
covS_reg = (1 - lam) * covS + lam * mean(diag(covS)) * eyeN;
covR_reg = (1 - lam) * covR + lam * mean(diag(covR)) * eyeN;
try
    [W, D] = eig(covS_reg, covR_reg);
catch
    return
end
evals = real(diag(D));
[evals, ord] = sort(evals, 'descend');
W = real(W(:, ord));
nSearch = min(params.ged_search_n, size(W, 2));
if nSearch < 1
    return
end

front_idx = find(mask_labels(labels, '^(Fp|AF|F)'));
post_idx = find(mask_labels(labels, '^(O|I|PO|PPO|P)'));
temp_idx = find(mask_labels(labels, '^(T|TP|FT)'));

searchFilters = nan(nChan, nSearch);
searchTopos = nan(nChan, nSearch);
post_front = nan(nSearch, 1);
post_temp = nan(nSearch, 1);
topo_peak_frac = nan(nSearch, 1);
for ci = 1:nSearch
    w = W(:, ci);
    topo = covS_reg * w;
    if ~isempty(post_idx)
        post_strength = mean(abs(topo(post_idx)));
    else
        post_strength = 0;
    end
    if ~isempty(front_idx)
        front_strength = mean(abs(topo(front_idx)));
        post_front(ci) = post_strength / max(front_strength, eps);
    else
        post_front(ci) = Inf;
    end
    if ~isempty(temp_idx)
        temp_strength = mean(abs(topo(temp_idx)));
        post_temp(ci) = post_strength / max(temp_strength, eps);
    else
        post_temp(ci) = Inf;
    end
    topo_abs = abs(topo(isfinite(topo)));
    if ~isempty(topo_abs) && sum(topo_abs) > 0
        topo_peak_frac(ci) = max(topo_abs) / sum(topo_abs);
    end
    searchFilters(:, ci) = w;
    searchTopos(:, ci) = topo;
end

prox = estimate_components_artifact_proxies( ...
    searchFilters, {paired}, peri_win, ctrl_win, fs, params.scan_freqs, params.tapsmo, '');
spec = nan(nSearch, numel(params.scan_freqs));
hf_slope = nan(nSearch, 1);
for ci = 1:nSearch
    hf_slope(ci) = prox(ci).hf_slope;
    sp = prox(ci).mean_pr_spectrum(:)';
    if numel(sp) == size(spec, 2)
        spec(ci, :) = sp;
    end
end
hf_slope(~isfinite(hf_slope)) = 0;
[pf, ~] = compute_powspctrm_form_laplacian_score_from_spectra( ...
    spec, params.scan_freqs, params.analysis_freq_range);
eval_raw = evals(1:nSearch);

[~, top_i] = max(eval_raw);
if isempty(top_i) || ~isfinite(eval_raw(top_i))
    top_i = 1;
end
sel.top_eig = eval_raw(top_i);
sel.top_pf = pf(top_i);
sel.top_post_front = post_front(top_i);
sel.top_post_temp = post_temp(top_i);
sel.top_hf_slope = hf_slope(top_i);
sel.top_peak_frac = topo_peak_frac(top_i);
sel.top_class = component_class(top_i, hf_slope, topo_peak_frac, post_front, post_temp);

lead_ok = isfinite(eval_raw(top_i)) && eval_raw(top_i) >= params.min_eigval && ...
    isfinite(topo_peak_frac(top_i)) && topo_peak_frac(top_i) <= 0.25;
if ~lead_ok
    return
end
w = searchFilters(:, top_i);
topo = searchTopos(:, top_i);
[~, imax] = max(abs(topo));
if ~isempty(imax) && isfinite(topo(imax)) && topo(imax) < 0
    topo = -topo;
    w = -w;
end
w = normalize_filters_to_noise_metric(w, covR);
sel.show = true;
sel.w = w;
sel.topo = topo;
sel.eig = eval_raw(top_i);
sel.spec = spec(top_i, :);
sel.in_spectrum = ~(hf_slope(top_i) > 0);
sel.passed = sel.in_spectrum;
end

function sel = empty_selection(nFreq)
sel.passed = false;
sel.show = false;
sel.in_spectrum = false;
sel.w = [];
sel.topo = [];
sel.eig = NaN;
sel.spec = nan(1, nFreq);
sel.top_eig = NaN;
sel.top_pf = NaN;
sel.top_post_front = NaN;
sel.top_post_temp = NaN;
sel.top_hf_slope = NaN;
sel.top_peak_frac = NaN;
sel.top_class = 'none';
end

function cls = component_class(ci, hf_slope, topo_peak_frac, post_front, post_temp)
if topo_peak_frac(ci) > 0.25
    cls = 'focal';
elseif hf_slope(ci) > 0
    cls = 'EMG';
elseif ~(post_front(ci) > 1)
    cls = 'frontal';
elseif ~(post_temp(ci) > 1)
    cls = 'temporal';
else
    cls = 'posterior';
end
end

function mask = mask_labels(labels, expr)
mask = cellfun(@(l) ~isempty(regexp(l, expr, 'once')), labels);
end

function plot_ms_ged_spectra(spec_ms, spec_gamma, scan_freqs, colors, condTitles, outPath)
nCond = size(spec_ms, 1);
nFreq = numel(scan_freqs);
ym = nan(nCond, 2, nFreq);
ye = nan(nCond, 2, nFreq);
n_line = zeros(nCond, 2);
specs = {spec_gamma, spec_ms};
for c = 1:nCond
    for b = 1:2
        curves = reshape(specs{b}(c, :, :), size(specs{b}, 2), nFreq);
        [ym(c, b, :), ye(c, b, :), n_line(c, b)] = grand_mean_sem_ft(curves, scan_freqs);
    end
end
yShow = [ym(:) - ye(:); ym(:) + ye(:); ym(:); 0];
yShow = yShow(isfinite(yShow));
if isempty(yShow)
    yLim = [-1 1];
else
    yPad = 0.08 * range(yShow);
    if yPad == 0
        yPad = 0.5;
    end
    yLim = [min(yShow) - yPad, max(yShow) + yPad];
end

fig = new_hidden_figure();
closeFig = onCleanup(@() close_if_open(fig)); %#ok<NASGU>
for c = 1:nCond
    subplot(2, 2, c);
    hold on
    leg_p = gobjects(0);
    leg_lbl = {};
    series_color = {[0.25 0.25 0.25], colors(c, :)};
    series_name = {'Gamma filter', 'MS filter'};
    series_style = {'-', '--'};
    for b = [1 2]
        mu = squeeze(ym(c, b, :))';
        sem = squeeze(ye(c, b, :))';
        if ~any(isfinite(mu))
            continue
        end
        if all(isfinite(mu)) && all(isfinite(sem))
            eb = shadedErrorBar(scan_freqs(:)', mu, sem, ...
                'lineProps', {series_style{b}}, 'transparent', true);
            set(eb.mainLine, 'Color', series_color{b}, 'LineWidth', 2.5);
            set(eb.patch, 'FaceColor', series_color{b}, 'FaceAlpha', 0.15);
            set(eb.edge(1), 'Color', 'none');
            set(eb.edge(2), 'Color', 'none');
            h = eb.mainLine;
        else
            h = plot(scan_freqs, mu, series_style{b}, ...
                'Color', series_color{b}, 'LineWidth', 2.5);
        end
        leg_p(end + 1) = h; %#ok<AGROW>
        leg_lbl{end + 1} = sprintf('%s (N=%d)', series_name{b}, n_line(c, b)); %#ok<AGROW>
    end
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xlim([scan_freqs(1), scan_freqs(end)]);
    ylim(yLim);
    xticks(30:10:90);
    title(condTitles{c}, 'FontSize', 28, 'FontWeight', 'bold');
    xlabel('Frequency [Hz]', 'FontSize', 22);
    ylabel('Power ratio [dB]', 'FontSize', 22);
    set(gca, 'FontSize', 18, 'Box', 'off');
    if ~isempty(leg_p)
        legend(leg_p, leg_lbl, 'Location', 'northeast', 'Box', 'off', 'FontSize', 16);
    else
        text(0.5, 0.5, 'Too few pairs', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 18);
    end
    box off
    hold off
end
save_png(fig, outPath);
fprintf('[MS GED] Saved %s\n', outPath);
end

function plot_ms_ged_metrics(topo_r, rayleigh_db, colors, condPct, outPath)
fig = new_hidden_figure();
closeFig = onCleanup(@() close_if_open(fig)); %#ok<NASGU>
subplot(1, 2, 1);
hold on
y = topo_r(isfinite(topo_r));
if isempty(y)
    axis off
    text(0.5, 0.5, 'No estimate', 'HorizontalAlignment', 'center', 'FontSize', 18);
else
    boxplot(y, 'Positions', 1, 'Colors', [0.2 0.2 0.2], 'Symbol', '', 'Widths', 0.55);
    styleCurrentBoxplot([0.35 0.35 0.35]);
    rng(1, 'twister');
    x = 1 + (rand(numel(y), 1) - 0.5) * 0.25;
    scatter(x, y, 160, [0.2 0.2 0.2], 'filled', ...
        'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'none');
    yline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1, 'HandleVisibility', 'off');
    xlim([0.4, 1.6]);
    xticks(1);
    xticklabels({'Pooled'});
    ylabel('Topography correlation', 'FontSize', 22);
    xlabel('Microsaccade filter', 'FontSize', 22);
    set(gca, 'FontSize', 18, 'Box', 'off');
    box off
end
hold off

subplot(1, 2, 2);
hold on
Y = rayleigh_db;
nCond = size(Y, 1);
if ~any(isfinite(Y(:)))
    axis off
    text(0.5, 0.5, 'No estimate', 'HorizontalAlignment', 'center', 'FontSize', 18);
else
    boxplot(Y', 'Positions', 1:nCond, 'Colors', [0.2 0.2 0.2], 'Symbol', '', 'Widths', 0.55);
    styleCurrentBoxplot(colors(1:nCond, :));
    rng(1, 'twister');
    for c = 1:nCond
        yc = Y(c, :);
        yc = yc(isfinite(yc));
        if isempty(yc)
            continue
        end
        x = c + (rand(numel(yc), 1) - 0.5) * 0.25;
        scatter(x, yc, 160, colors(c, :), 'filled', ...
            'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'none');
    end
    yline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1, 'HandleVisibility', 'off');
    xlim([0.4, nCond + 0.6]);
    xticks(1:nCond);
    xticklabels(arrayfun(@(p) sprintf('%d%%', p), condPct, 'UniformOutput', false));
    ylabel('Gamma filter ratio [dB]', 'FontSize', 22);
    xlabel('Contrast', 'FontSize', 22);
    set(gca, 'FontSize', 18, 'Box', 'off');
    box off
end
hold off
save_png(fig, outPath);
fprintf('[MS GED] Saved %s\n', outPath);
end

function plot_ms_ged_topos(topos_ms, topos_gamma, topo_labels, topo_r, top_class, ...
    subjects, layout, fig_dir)
cfg0 = [];
cfg0.layout = layout;
cfg0.comment = 'no';
cfg0.marker = 'off';
cfg0.style = 'straight';
cfg0.gridscale = 150;
cfg0.colormap = '*RdBu';
cfg0.interactive = 'no';
cfg0.figure = 'gcf';
for s = 1:numel(subjects)
    if isempty(topo_labels{s})
        continue
    end
    fig = new_hidden_figure(491);
    closeFig = onCleanup(@() close_if_open(fig)); %#ok<NASGU>
    panels = cell(1, 2);
    panel_titles = cell(1, 2);
    panels{1} = topos_gamma{s};
    panel_titles{1} = 'Gamma filter';
    panels{2} = topos_ms{s};
    cls = top_class{s};
    if isempty(cls)
        cls = 'none';
    end
    if isempty(panels{2}) || ~any(isfinite(panels{2}))
        panel_titles{2} = 'Microsaccade filter   none';
    elseif isfinite(topo_r(s))
        panel_titles{2} = sprintf('Microsaccade filter   %s   r = %.2f', cls, topo_r(s));
    else
        panel_titles{2} = sprintf('Microsaccade filter   %s', cls);
    end
    for p = 1:2
        set(0, 'CurrentFigure', fig);
        subplot(1, 2, p);
        topo = panels{p};
        if isempty(topo) || ~any(isfinite(topo))
            axis off
            text(0.5, 0.5, 'None', 'HorizontalAlignment', 'center', 'FontSize', 16);
            title(panel_titles{p}, 'FontSize', 16, 'FontWeight', 'bold');
            continue
        end
        topo_data = [];
        topo_data.label = topo_labels{s};
        topo_data.avg = topo(:);
        topo_data.dimord = 'chan';
        clim = max(abs(topo(isfinite(topo))));
        if ~isfinite(clim) || clim <= 0
            clim = 1;
        end
        cfg = cfg0;
        cfg.zlim = [-clim clim];
        try
            ft_topoplotER(cfg, topo_data);
        catch
            axis off
            text(0.5, 0.5, 'Topo failed', 'HorizontalAlignment', 'center', 'FontSize', 14);
        end
        title(panel_titles{p}, 'FontSize', 16, 'FontWeight', 'bold');
    end
    outPath = fullfile(fig_dir, sprintf('GCP_ms_ged_topo_%s.png', subjects{s}));
    save_png(fig, outPath);
    fprintf('[MS GED] Saved %s\n', outPath);
    clear closeFig
end
end

function [mu, sem, n_used] = grand_mean_sem_ft(curves, scan_freqs)
nFreq = numel(scan_freqs);
mu = nan(1, nFreq);
sem = nan(1, nFreq);
n_used = 0;
pow_per_subj = {};
for i = 1:size(curves, 1)
    row = curves(i, :);
    if ~any(isfinite(row))
        continue
    end
    pow_per_subj{end + 1} = ged_powcurve_to_freq_ft(row, scan_freqs); %#ok<AGROW>
end
n_used = numel(pow_per_subj);
if n_used < 1
    return
end
cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter = 'powspctrm';
ga = ft_freqgrandaverage(cfg, pow_per_subj{:});
nf = numel(ga.freq);
pow = reshape(double(ga.powspctrm), size(ga.powspctrm, 1), nf);
n_fin = sum(isfinite(pow), 1);
mu = mean(pow, 1, 'omitnan');
sem = std(pow, 0, 1, 'omitnan') ./ sqrt(max(1, n_fin));
sem(n_fin < 2) = NaN;
end

function styleCurrentBoxplot(boxColors)
hBoxes = findobj(gca, 'Tag', 'Box');
if isempty(hBoxes)
    return
end
nB = numel(hBoxes);
mx = zeros(nB, 1);
for ii = 1:nB
    xd = get(hBoxes(ii), 'XData');
    mx(ii) = mean(xd(:), 'omitnan');
end
[~, ord] = sort(mx);
hBoxes = hBoxes(ord);
for bi = 1:min(numel(hBoxes), size(boxColors, 1))
    patch(get(hBoxes(bi), 'XData'), get(hBoxes(bi), 'YData'), ...
        boxColors(bi, :), 'FaceAlpha', 0.25, ...
        'EdgeColor', boxColors(bi, :), 'LineWidth', 1.5);
end
delete(findobj(gca, 'Tag', 'Box'));
end

function fig = new_hidden_figure(fig_height)
if nargin < 1 || isempty(fig_height)
    fig_height = 982;
end
fig = figure('Visible', 'off', 'Position', [0 0 1512 fig_height], 'Color', 'w', ...
    'MenuBar', 'none', 'ToolBar', 'none', 'DockControls', 'off');
set(fig, 'DefaultAxesToolbarVisible', 'off');
end

function close_if_open(fig)
if isgraphics(fig)
    close(fig);
end
end

function save_png(fig, outPath)
set(fig, 'Visible', 'off');
axs = findall(fig, 'Type', 'axes');
for iAx = 1:numel(axs)
    if isprop(axs(iAx), 'Toolbar') && ~isempty(axs(iAx).Toolbar)
        axs(iAx).Toolbar.Visible = 'off';
    end
end
set(fig, 'PaperPositionMode', 'auto');
print(fig, outPath, '-dpng', '-r600');
end
