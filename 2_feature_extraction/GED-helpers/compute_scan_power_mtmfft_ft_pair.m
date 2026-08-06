function [p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(sig_stim, sig_base, fs, scan_freqs, tapsmofrq_hz)
% Joint stim/baseline multitaper power spectra for matched signal sets.
if ~iscell(sig_stim) && isvector(sig_stim)
    sig_stim = sig_stim(:)';
end
if ~iscell(sig_base) && isvector(sig_base)
    sig_base = sig_base(:)';
end
nFreq = numel(scan_freqs);
if iscell(sig_stim)
    nStim = numel(sig_stim);
else
    nStim = size(sig_stim, 1);
end
p_stim_scan = nan(nStim, nFreq);
p_base_scan = nan(nStim, nFreq);
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
    sig_all = [sig_stim(:); sig_base(:)];
else
    if size(sig_stim, 1) ~= size(sig_base, 1)
        return;
    end
    sig_all = [sig_stim; sig_base];
end
if isempty(sig_all)
    return;
end
p_all = compute_scan_power_mtmfft_ft(sig_all, fs, scan_freqs, tapsmofrq_hz);
if isempty(p_all) || size(p_all, 1) < 2 * nStim
    return;
end
p_stim_scan = p_all(1:nStim, :);
p_base_scan = p_all((nStim + 1):(2 * nStim), :);
end
