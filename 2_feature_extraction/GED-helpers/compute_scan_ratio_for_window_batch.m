function [ratio_cube, near_floor_freq_count_per_trial] = compute_scan_ratio_for_window_batch(trial_cache, search_filters, stim_field, trial_mask, fs, scan_freqs, tapsmofrq_hz, base_floor, near_floor_mult)
% Stim/baseline dB power-ratio spectra for each search component across trials.
nTrl = numel(trial_cache);
nComp = size(search_filters, 2);
nFreq = numel(scan_freqs);
ratio_cube = nan(nComp, nTrl, nFreq);
near_floor_freq_count_per_trial = zeros(nTrl, 1);
if nComp == 0 || nTrl == 0 || ~any(trial_mask)
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
    stim_trials{vi} = search_filters' * x_stim;
    base_trials{vi} = search_filters' * x_base;
    keep(vi) = true;
end
if ~any(keep)
    return;
end
stim_trials = stim_trials(keep);
base_trials = base_trials(keep);
valid_trls = valid_trls(keep);
nValid = numel(valid_trls);

[p_stim, p_base] = compute_scan_power_mtmfft_ft_pair_trialchans( ...
    stim_trials, base_trials, fs, scan_freqs, tapsmofrq_hz);
if isempty(p_stim) || isempty(p_base)
    return;
end

if ~isfinite(near_floor_mult) || near_floor_mult <= 0
    near_floor_mult = 1.5;
end
floor_fallback = max(base_floor, eps);
near_floor_comp_rows = false(nValid, nComp, nFreq);
for vi = 1:nValid
    trl = valid_trls(vi);
    for ci = 1:nComp
        p_stim_row = double(squeeze(p_stim(vi, ci, :))).';
        p_base_row = double(squeeze(p_base(vi, ci, :))).';
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
        ratio_row = nan(1, nFreq);
        ratio_row(valid) = 10 * log10(p_stim_row(valid) ./ p_base_row(valid));
        ratio_cube(ci, trl, :) = ratio_row;
        near_floor_comp_rows(vi, ci, valid) = p_base_row(valid) <= near_floor_mult * floor_row;
    end
    nf = reshape(near_floor_comp_rows(vi, :, :), nComp, nFreq);
    near_floor_freq_count_per_trial(trl) = sum(mean(nf, 1) >= 0.5);
end
end
