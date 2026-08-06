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

sig_stim_cells = cell(0, 1);
sig_base_cells = cell(0, 1);
row_comp_idx = zeros(0, 1);
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
    comp_base = search_filters' * x_base;
    comp_stim = search_filters' * x_stim;
    for ci = 1:nComp
        sig_base_cells{end+1, 1} = comp_base(ci, :);
        sig_stim_cells{end+1, 1} = comp_stim(ci, :);
        row_comp_idx(end+1, 1) = ci;
        row_trial_idx(end+1, 1) = trl;
    end
end
if isempty(sig_stim_cells)
    return;
end

[ratio_rows, ~, near_floor_row_mask] = compute_scan_ratio_from_timeseries( ...
    sig_stim_cells, sig_base_cells, fs, scan_freqs, tapsmofrq_hz, base_floor, near_floor_mult);
for ri = 1:size(ratio_rows, 1)
    ratio_cube(row_comp_idx(ri), row_trial_idx(ri), :) = ratio_rows(ri, :);
end
for trl = 1:nTrl
    rows = row_trial_idx == trl;
    if ~any(rows)
        continue;
    end
    near_floor_trial_mask = mean(near_floor_row_mask(rows, :), 1) >= 0.5;
    near_floor_freq_count_per_trial(trl) = sum(near_floor_trial_mask);
end
end
