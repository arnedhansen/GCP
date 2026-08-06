function proxy = estimate_component_artifact_proxies(filter_w, dat_per_cond, stim_window, base_window, fs, scan_freqs, tapsmofrq_hz)
% Component spectrum proxies: mean power ratio, HF slope, and line-harmonic ratio.
proxy = struct( ...
    'lineharm_ratio', NaN, ...
    'stationarity_cv', NaN, ...
    'burst_ratio', NaN, ...
    'hf_slope', NaN, ...
    'mean_pr_spectrum', nan(1, numel(scan_freqs)), ...
    'n_trials_used', 0);
if isempty(filter_w) || isempty(dat_per_cond)
    return;
end

band_mask = scan_freqs >= 30 & scan_freqs <= 90;
hf_mask = scan_freqs >= 70 & scan_freqs <= min(110, max(scan_freqs));
harm_mask = ismember(round(scan_freqs), [50 60 100]);
if ~any(harm_mask)
    [~, idx50] = min(abs(scan_freqs - 50));
    harm_mask(idx50) = true;
end

trial_gamma = [];
trial_pr = [];
lineharm_acc = [];
trial_count = 0;
stim_sig_cells = {};
base_sig_cells = {};
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
        x_proj = filter_w' * x;
        idx_base = t >= base_window(1) & t <= base_window(2);
        idx_stim = t >= stim_window(1) & t <= stim_window(2);
        if sum(idx_base) < 5 || sum(idx_stim) < 5
            continue;
        end
        base_sig_cells{end+1, 1} = x_proj(idx_base);
        stim_sig_cells{end+1, 1} = x_proj(idx_stim);
        trial_count = trial_count + 1;
    end
end
proxy.n_trials_used = trial_count;
if trial_count < 4
    return;
end

if ~isfinite(fs) || fs <= 0
    return;
end
trial_pr = compute_simple_power_ratio_scan_batch(stim_sig_cells, base_sig_cells, fs, scan_freqs, tapsmofrq_hz);
if isempty(trial_pr)
    return;
end
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
    trial_gamma(end+1, 1) = nanmean(pr_band);
    lineharm_acc(end+1, 1) = max(harm_val, 0) / max(abs(nonharm_val), eps);
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
