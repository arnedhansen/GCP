function pr_mat = compute_simple_power_ratio_scan_batch(sig_stim_cells, sig_base_cells, fs, scan_freqs, tapsmofrq_hz)
% dB stim/baseline power ratios for component-selection proxy spectra.
nRows = numel(sig_stim_cells);
pr_mat = nan(nRows, numel(scan_freqs));
if nRows == 0 || numel(sig_base_cells) ~= nRows || fs <= 0
    return;
end
[p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(sig_stim_cells, sig_base_cells, fs, scan_freqs, tapsmofrq_hz);
if isempty(p_stim_scan) || isempty(p_base_scan)
    return;
end
for ri = 1:nRows
    p_stim_row = p_stim_scan(ri, :);
    p_base_row = p_base_scan(ri, :);
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
    ratio_row = nan(1, numel(scan_freqs));
    ratio_row(valid) = 10 * log10((p_stim_row(valid) + base_floor) ./ (p_base_row(valid) + base_floor));
    pr_mat(ri, :) = ratio_row;
end
end
