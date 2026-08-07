function [ratio_trials, near_floor_freq_count_per_trial, valid_freq_count_per_trial] = ...
    compute_scan_ratio_for_combined_filter_batch(trial_cache, filter_vec, stim_field, trial_mask, fs, scan_freqs, tapsmofrq_hz, base_floor, near_floor_mult)
% Stim/baseline dB power-ratio spectra for one combined GED filter across trials.
nTrl = numel(trial_cache);
ratio_trials = nan(nTrl, numel(scan_freqs));
near_floor_freq_count_per_trial = zeros(nTrl, 1);
valid_freq_count_per_trial = zeros(nTrl, 1);
if nTrl == 0 || isempty(filter_vec) || ~any(trial_mask)
    return;
end
filter_vec = filter_vec(:);
if ~all(isfinite(filter_vec))
    return;
end

valid_trls = find(trial_mask(:));
nValid = numel(valid_trls);
stim_trials = cell(nValid, 1);
base_trials = cell(nValid, 1);
keep = false(nValid, 1);
for vi = 1:nValid
    trl = valid_trls(vi);
    tc = trial_cache{trl};
    x_base = tc.x_base;
    x_stim = tc.(stim_field);
    if isempty(x_base) || isempty(x_stim)
        continue;
    end
    stim_trials{vi} = filter_vec' * x_stim;
    base_trials{vi} = filter_vec' * x_base;
    keep(vi) = true;
end
if ~any(keep)
    return;
end
stim_trials = stim_trials(keep);
base_trials = base_trials(keep);
valid_trls = valid_trls(keep);

[p_stim, p_base] = compute_scan_power_mtmfft_ft_pair_trialchans( ...
    stim_trials, base_trials, fs, scan_freqs, tapsmofrq_hz);
if isempty(p_stim) || isempty(p_base)
    return;
end

if ~isfinite(near_floor_mult) || near_floor_mult <= 0
    near_floor_mult = 1.5;
end
floor_fallback = max(base_floor, eps);
for vi = 1:numel(valid_trls)
    trl = valid_trls(vi);
    p_stim_row = double(squeeze(p_stim(vi, 1, :))).';
    p_base_row = double(squeeze(p_base(vi, 1, :))).';
    valid_base = isfinite(p_base_row) & (p_base_row > 0);
    if ~any(valid_base)
        continue;
    end
    base_anchor = prctile(p_base_row(valid_base), 20);
    base_median = median(p_base_row(valid_base), 'omitnan');
    if ~isfinite(base_anchor) || base_anchor <= 0
        base_anchor = base_median;
    end
    if ~isfinite(base_anchor) || base_anchor <= 0
        base_anchor = floor_fallback;
    end
    floor_row = max(0.25 * base_anchor, eps);
    valid = isfinite(p_stim_row) & isfinite(p_base_row) & ...
        (p_stim_row > 0) & (p_base_row > 0);
    ratio_row = nan(1, numel(scan_freqs));
    ratio_row(valid) = 10 * log10(p_stim_row(valid) ./ p_base_row(valid));
    ratio_trials(trl, :) = ratio_row;
    valid_freq_count_per_trial(trl) = sum(isfinite(ratio_row));
    near_floor_freq_count_per_trial(trl) = sum(p_base_row(valid) <= near_floor_mult * floor_row);
end
end
