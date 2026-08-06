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

sig_stim_cells = cell(0, 1);
sig_base_cells = cell(0, 1);
row_trial_idx = zeros(0, 1);
for trl = 1:nTrl
    if ~trial_mask(trl)
        continue;
    end
    tc = trial_cache{trl};
    x_base = tc.x_base;
    x_stim = tc.(stim_field);
    if isempty(x_base) || isempty(x_stim)
        continue;
    end
    sig_base_cells{end+1, 1} = (filter_vec' * x_base);
    sig_stim_cells{end+1, 1} = (filter_vec' * x_stim);
    row_trial_idx(end+1, 1) = trl;
end
if isempty(sig_stim_cells)
    return;
end

[ratio_rows, ~, near_floor_row_mask] = compute_scan_ratio_from_timeseries( ...
    sig_stim_cells, sig_base_cells, fs, scan_freqs, tapsmofrq_hz, base_floor, near_floor_mult);
for ri = 1:size(ratio_rows, 1)
    trl = row_trial_idx(ri);
    ratio_trials(trl, :) = ratio_rows(ri, :);
    valid_freq_count_per_trial(trl) = sum(isfinite(ratio_rows(ri, :)));
    near_floor_freq_count_per_trial(trl) = sum(near_floor_row_mask(ri, :));
end
end
