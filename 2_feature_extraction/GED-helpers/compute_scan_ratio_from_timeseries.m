function [ratio_db, near_floor_freq_mask, near_floor_row_mask] = compute_scan_ratio_from_timeseries(sig_stim, sig_base, fs, scan_freqs, tapsmofrq_hz, base_floor, near_floor_mult)
% Stim/baseline dB power-ratio spectra and near-floor frequency masks.
if ~iscell(sig_stim) && isvector(sig_stim)
    sig_stim = sig_stim(:)';
end
if ~iscell(sig_base) && isvector(sig_base)
    sig_base = sig_base(:)';
end
if iscell(sig_stim)
    nSig = numel(sig_stim);
else
    nSig = size(sig_stim, 1);
end
ratio_db = nan(nSig, numel(scan_freqs));
near_floor_freq_mask = false(1, numel(scan_freqs));
near_floor_row_mask = false(nSig, numel(scan_freqs));
if isempty(sig_stim) || isempty(sig_base) || fs <= 0
    return;
end
if iscell(sig_stim) ~= iscell(sig_base)
    return;
end
if iscell(sig_stim)
    if numel(sig_stim) ~= numel(sig_base)
        return;
    end
else
    if size(sig_stim, 1) ~= size(sig_base, 1)
        return;
    end
end
if nSig == 0
    return;
end
if ~isfinite(near_floor_mult) || near_floor_mult <= 0
    near_floor_mult = 1.5;
end
% Each signal is represented as one virtual channel (GED component time series).
[p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(sig_stim, sig_base, fs, scan_freqs, tapsmofrq_hz);
if isempty(p_stim_scan) || isempty(p_base_scan)
    return;
end
% Derive a baseline floor for instability detection only. Valid power ratios
% use linear stimulus and baseline power directly before conversion to dB.
floor_fallback = max(base_floor, eps);
for ri = 1:nSig
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
    if ~isfinite(base_anchor) || base_anchor <= 0
        base_anchor = floor_fallback;
    end
    floor_row = max(0.25 * base_anchor, eps);

    valid = isfinite(p_stim_row) & isfinite(p_base_row) & ...
        (p_stim_row > 0) & (p_base_row > 0);
    ratio_row = nan(1, numel(scan_freqs));
    ratio_row(valid) = 10 * log10(p_stim_row(valid) ./ p_base_row(valid));
    ratio_db(ri, :) = ratio_row;
    near_floor_row_mask(ri, valid) = p_base_row(valid) <= near_floor_mult * floor_row;
end

for fi = 1:numel(scan_freqs)
    valid_fi = isfinite(ratio_db(:, fi));
    if any(valid_fi)
        near_floor_freq_mask(fi) = mean(near_floor_row_mask(valid_fi, fi)) >= 0.5;
    end
end
end
