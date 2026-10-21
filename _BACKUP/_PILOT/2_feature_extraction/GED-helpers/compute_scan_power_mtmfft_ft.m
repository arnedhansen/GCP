function p_scan = compute_scan_power_mtmfft_ft(sig, fs, scan_freqs, tapsmofrq_hz)
% Multitaper power spectra on scan_freqs using FieldTrip mtmfft.
%
% Input sig:
%   - cell of row/column vectors: one spectrum row per cell
%   - numeric matrix [nSig x nTime]: one spectrum row per matrix row
%
% Equal-length series are packed as channels of shared FieldTrip trials with an
% explicit common nextpow2 pad so results match the previous 1-channel packing.
if ~iscell(sig) && isvector(sig)
    sig = sig(:)';
end
if iscell(sig)
    nSig = numel(sig);
else
    nSig = size(sig, 1);
end
p_scan = nan(nSig, numel(scan_freqs));
if isempty(sig) || fs <= 0 || isempty(scan_freqs) || ~isfinite(tapsmofrq_hz) || tapsmofrq_hz <= 0
    return;
end
if ~iscell(sig) && size(sig, 2) < 8
    return;
end

valid_rows = false(nSig, 1);
row_data = cell(nSig, 1);
row_len = zeros(nSig, 1);
for si = 1:nSig
    if iscell(sig)
        x = double(sig{si});
    else
        x = double(sig(si, :));
    end
    if iscolumn(x)
        x = x.';
    end
    if any(~isfinite(x))
        continue;
    end
    if numel(x) < 8
        continue;
    end
    x = x - mean(x);
    valid_rows(si) = true;
    row_data{si} = x;
    row_len(si) = numel(x);
end
if ~any(valid_rows)
    return;
end

max_len = max(row_len(valid_rows));
pad_samp = 2^ceil(log2(max_len));
pad_sec = pad_samp / fs;
valid_idx = find(valid_rows);
lengths = unique(row_len(valid_idx));
pow_valid = nan(numel(valid_idx), numel(scan_freqs));

for li = 1:numel(lengths)
    L = lengths(li);
    group_local = find(row_len(valid_idx) == L);
    group_sig = valid_idx(group_local);
    nGroup = numel(group_sig);
    % Pack independent equal-length series as channels (mtmfft is per-channel).
    max_chan = 256;
    nBlock = ceil(nGroup / max_chan);
    for bi = 1:nBlock
        i0 = (bi - 1) * max_chan + 1;
        i1 = min(bi * max_chan, nGroup);
        block_sig = group_sig(i0:i1);
        nChan = numel(block_sig);
        dat = [];
        dat.label = arrayfun(@(k) sprintf('GED%03d', k), 1:nChan, 'UniformOutput', false);
        dat.fsample = fs;
        trial_mat = zeros(nChan, L);
        for ci = 1:nChan
            trial_mat(ci, :) = row_data{block_sig(ci)};
        end
        dat.trial = {trial_mat};
        dat.time = {(0:(L - 1)) / fs};
        dat.sampleinfo = [1, L];

        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'pow';
        cfg.taper = 'dpss';
        cfg.foi = scan_freqs;
        cfg.tapsmofrq = tapsmofrq_hz;
        cfg.pad = pad_sec;
        cfg.keeptrials = 'yes';
        cfg.feedback = 'none';
        try
            freq = ft_freqanalysis(cfg, dat);
        catch
            continue;
        end
        pow = double(freq.powspctrm);
        if ndims(pow) == 3
            % [1 x nChan x nFreq] or [nChan x 1 x nFreq] after squeeze quirks
            pow = reshape(pow, [], size(pow, 3));
            if size(pow, 1) ~= nChan && size(pow, 2) == nChan
                pow = pow.';
            end
        elseif isvector(pow)
            pow = reshape(pow, 1, []);
        end
        if isfield(freq, 'freq') && ~isempty(freq.freq)
            freq_axis = freq.freq(:)';
        else
            freq_axis = scan_freqs(:);
        end
        if size(pow, 2) ~= numel(scan_freqs)
            pow_interp = nan(size(pow, 1), numel(scan_freqs));
            for ri = 1:size(pow, 1)
                if numel(pow(ri, :)) ~= numel(freq_axis)
                    continue;
                end
                pow_interp(ri, :) = interp1(freq_axis, pow(ri, :), scan_freqs, 'linear', NaN);
            end
            pow = pow_interp;
        end
        if size(pow, 1) ~= nChan
            continue;
        end
        pow_valid(group_local(i0:i1), :) = pow;
    end
end

p_scan(valid_rows, :) = pow_valid;
end
