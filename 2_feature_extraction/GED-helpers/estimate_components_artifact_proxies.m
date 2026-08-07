function proxies = estimate_components_artifact_proxies(filters_w, dat_per_cond, stim_window, base_window, fs, scan_freqs, tapsmofrq_hz, cw_prefix)
% Batched component spectrum proxies for one or more GED filters.
% Projects all trials onto all filters, then one multi-channel mtmfft pair
% per window with a shared nextpow2 pad.
if nargin < 8
    cw_prefix = '';
end

nFreq = numel(scan_freqs);
empty_proxy = struct( ...
    'lineharm_ratio', NaN, ...
    'stationarity_cv', NaN, ...
    'burst_ratio', NaN, ...
    'hf_slope', NaN, ...
    'mean_pr_spectrum', nan(1, nFreq), ...
    'n_trials_used', 0);

if isempty(filters_w) || isempty(dat_per_cond)
    proxies = repmat(empty_proxy, 0, 1);
    return;
end
if isvector(filters_w)
    filters_w = filters_w(:);
end
nComp = size(filters_w, 2);
proxies = repmat(empty_proxy, nComp, 1);

band_mask = scan_freqs >= 30 & scan_freqs <= 90;
hf_mask = scan_freqs >= 70 & scan_freqs <= min(110, max(scan_freqs));
harm_mask = ismember(round(scan_freqs), [50 60 100]);
if ~any(harm_mask)
    [~, idx50] = min(abs(scan_freqs - 50));
    harm_mask(idx50) = true;
end

n_per_cond = zeros(1, numel(dat_per_cond));
n_guess = 0;
for cond = 1:numel(dat_per_cond)
    dat = dat_per_cond{cond};
    if ~isempty(dat) && isfield(dat, 'trial')
        n_guess = n_guess + numel(dat.trial);
    end
end
stim_seg = cell(n_guess, 1);
base_seg = cell(n_guess, 1);
nTrials = 0;
for cond = 1:numel(dat_per_cond)
    dat = dat_per_cond{cond};
    if isempty(dat) || ~isfield(dat, 'trial')
        continue;
    end
    for trl = 1:numel(dat.trial)
        x = double(dat.trial{trl});
        t = dat.time{trl};
        if isempty(x) || isempty(t)
            continue;
        end
        idx_base = t >= base_window(1) & t <= base_window(2);
        idx_stim = t >= stim_window(1) & t <= stim_window(2);
        if sum(idx_base) < 5 || sum(idx_stim) < 5
            continue;
        end
        nTrials = nTrials + 1;
        base_seg{nTrials} = x(:, idx_base);
        stim_seg{nTrials} = x(:, idx_stim);
        n_per_cond(cond) = n_per_cond(cond) + 1;
    end
end
stim_seg = stim_seg(1:nTrials);
base_seg = base_seg(1:nTrials);
for ci = 1:nComp
    proxies(ci).n_trials_used = nTrials;
end
if nTrials < 4
    return;
end
if ~isfinite(fs) || fs <= 0
    return;
end

if ~isempty(cw_prefix)
    cond_parts = arrayfun(@(n) sprintf('%d', n), n_per_cond, 'UniformOutput', false);
    fprintf('%s %s = %d trials × %d components × 2 (stim + base)\n', ...
        cw_prefix, strjoin(cond_parts, ' + '), nTrials, nComp);
end

stim_trials = cell(nTrials, 1);
base_trials = cell(nTrials, 1);
for trl = 1:nTrials
    stim_trials{trl} = filters_w' * stim_seg{trl};
    base_trials{trl} = filters_w' * base_seg{trl};
end
[p_stim, p_base] = compute_scan_power_mtmfft_ft_pair_trialchans( ...
    stim_trials, base_trials, fs, scan_freqs, tapsmofrq_hz);
if isempty(p_stim) || isempty(p_base)
    return;
end

pr_all = nan(nComp * nTrials, nFreq);
for trl = 1:nTrials
    for ci = 1:nComp
        p_stim_row = double(squeeze(p_stim(trl, ci, :))).';
        p_base_row = double(squeeze(p_base(trl, ci, :))).';
        valid_base = isfinite(p_base_row) & (p_base_row > 0);
        if ~any(valid_base)
            continue;
        end
        base_anchor = prctile(p_base_row(valid_base), 20);
        base_median = median(p_base_row(valid_base), 'omitnan');
        if ~isfinite(base_anchor) || base_anchor <= 0
            base_anchor = base_median;
        end
        if ~isfinite(base_median) || base_median <= 0
            base_median = base_anchor;
        end
        base_floor = max(0.25 * base_anchor, eps);
        valid = isfinite(p_stim_row) & isfinite(p_base_row) & (p_base_row > 0);
        ratio_row = nan(1, nFreq);
        ratio_row(valid) = 10 * log10((p_stim_row(valid) + base_floor) ./ (p_base_row(valid) + base_floor));
        pr_all((ci - 1) * nTrials + trl, :) = ratio_row;
    end
end

for ci = 1:nComp
    trial_pr = pr_all((ci - 1) * nTrials + (1:nTrials), :);
    proxies(ci) = summarize_component_proxy_from_trial_pr( ...
        trial_pr, nTrials, scan_freqs, band_mask, hf_mask, harm_mask, empty_proxy);
end
end

function proxy = summarize_component_proxy_from_trial_pr(trial_pr, nTrials, scan_freqs, band_mask, hf_mask, harm_mask, empty_proxy)
proxy = empty_proxy;
proxy.n_trials_used = nTrials;
if isempty(trial_pr)
    return;
end

trial_gamma = [];
lineharm_acc = [];
for ri = 1:size(trial_pr, 1)
    pr_full = trial_pr(ri, :);
    if all(~isfinite(pr_full)) || ~any(band_mask)
        continue;
    end
    pr_band = pr_full(band_mask);
    if all(~isfinite(pr_band))
        continue;
    end
    nonharm_val = nanmean(pr_full(band_mask & ~harm_mask));
    harm_val = nanmean(pr_full(band_mask & harm_mask));
    if ~isfinite(nonharm_val) || abs(nonharm_val) <= eps || ~isfinite(harm_val)
        continue;
    end
    trial_gamma(end+1, 1) = nanmean(pr_band); %#ok<AGROW>
    lineharm_acc(end+1, 1) = max(harm_val, 0) / max(abs(nonharm_val), eps); %#ok<AGROW>
end
if isempty(trial_gamma) || isempty(lineharm_acc)
    return;
end
proxy.mean_pr_spectrum = nanmean(trial_pr, 1);
proxy.lineharm_ratio = nanmean(lineharm_acc);
proxy.stationarity_cv = nanstd(trial_gamma) / max(abs(nanmean(trial_gamma)), eps);
burst_thr = prctile(abs(trial_gamma), 85);
if ~isfinite(burst_thr) || burst_thr <= 0
    burst_thr = median(abs(trial_gamma)) + mad(abs(trial_gamma), 1);
end
proxy.burst_ratio = mean(abs(trial_gamma) >= max(burst_thr, eps));

if sum(hf_mask) >= 3
    f_hf = scan_freqs(hf_mask);
    spec_hf = proxy.mean_pr_spectrum(hf_mask);
    xh_full = log(f_hf(:));
    yh_full = log(abs(spec_hf(:)) + eps);
    valid_hf = isfinite(xh_full) & isfinite(yh_full) & ...
        isfinite(spec_hf(:)) & (spec_hf(:) > -inf);
    n_hf = nnz(valid_hf);
    if n_hf >= 3
        xh = xh_full(valid_hf);
        yh = yh_full(valid_hf);
        p = polyfit(xh, yh, 1);
        if isfinite(p(1))
            proxy.hf_slope = p(1);
        end
    end
end
end
