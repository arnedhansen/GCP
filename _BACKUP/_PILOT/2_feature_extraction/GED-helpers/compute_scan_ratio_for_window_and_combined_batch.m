function [ratio_cube, ratio_combined, near_floor_comp_count, near_floor_comb_count, valid_freq_comb_count] = ...
    compute_scan_ratio_for_window_and_combined_batch(trial_cache, search_filters, w_components, stim_field, trial_mask, fs, scan_freqs, tapsmofrq_hz, base_floor, near_floor_mult)
% Component + combined-filter stim/baseline dB ratios in one spectral batch.
% Combined series is formed in time as w'*(W'*x) before FFT (not after).
nTrl = numel(trial_cache);
nComp = size(search_filters, 2);
nFreq = numel(scan_freqs);
ratio_cube = nan(nComp, nTrl, nFreq);
ratio_combined = nan(nTrl, nFreq);
near_floor_comp_count = zeros(nTrl, 1);
near_floor_comb_count = zeros(nTrl, 1);
valid_freq_comb_count = zeros(nTrl, 1);
if nComp == 0 || nTrl == 0 || ~any(trial_mask)
    return;
end

filter_vec = build_combined_filter_vector(search_filters, w_components);
include_combined = ~isempty(w_components) && ~isempty(filter_vec) && all(isfinite(filter_vec));
nSeries = nComp + double(include_combined);

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
    comp_base = search_filters' * x_base;
    comp_stim = search_filters' * x_stim;
    if include_combined
        series_base = [comp_base; filter_vec' * x_base];
        series_stim = [comp_stim; filter_vec' * x_stim];
    else
        series_base = comp_base;
        series_stim = comp_stim;
    end
    stim_trials{vi} = series_stim;
    base_trials{vi} = series_base;
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
    near_floor_rows_trl = false(nSeries, nFreq);
    for si = 1:nSeries
        p_stim_row = double(squeeze(p_stim(vi, si, :))).';
        p_base_row = double(squeeze(p_base(vi, si, :))).';
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
        near_floor_rows_trl(si, valid) = p_base_row(valid) <= near_floor_mult * floor_row;
        if si <= nComp
            ratio_cube(si, trl, :) = ratio_row;
            near_floor_comp_rows(vi, si, :) = near_floor_rows_trl(si, :);
        else
            ratio_combined(trl, :) = ratio_row;
            valid_freq_comb_count(trl) = sum(isfinite(ratio_row));
            near_floor_comb_count(trl) = sum(near_floor_rows_trl(si, :));
        end
    end
    if nComp > 0
        nf = reshape(near_floor_comp_rows(vi, :, :), nComp, nFreq);
        near_floor_comp_count(trl) = sum(mean(nf, 1) >= 0.5);
    end
end
end
