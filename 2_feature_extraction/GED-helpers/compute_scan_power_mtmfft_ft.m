function p_scan = compute_scan_power_mtmfft_ft(sig, fs, scan_freqs, tapsmofrq_hz)
% Multitaper power spectra on scan_freqs using FieldTrip mtmfft.
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
dat = [];
dat.label = {'GED'};
dat.fsample = fs;
dat.trial = {};
dat.time = {};
for si = 1:nSig
    if iscell(sig)
        x = double(sig{si});
    else
        x = double(sig(si, :));
    end
    if iscolumn(x)
        x = x';
    end
    if any(~isfinite(x))
        continue;
    end
    if numel(x) < 8
        continue;
    end
    x = x - mean(x);
    valid_rows(si) = true;
    dat.trial{end+1} = x;
    dat.time{end+1} = (0:(numel(x)-1)) / fs;
end
if ~any(valid_rows)
    return;
end

% Provide explicit sampleinfo to avoid FieldTrip reconstruction warnings.
nTrials_valid = numel(dat.trial);
dat.sampleinfo = zeros(nTrials_valid, 2);
sample_start = 1;
for ti = 1:nTrials_valid
    nSamp = numel(dat.trial{ti});
    dat.sampleinfo(ti, :) = [sample_start, sample_start + nSamp - 1];
    sample_start = sample_start + nSamp;
end

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'dpss';
cfg.foi = scan_freqs;
cfg.tapsmofrq = tapsmofrq_hz;
cfg.pad = 'nextpow2';
cfg.keeptrials = 'yes';
cfg.feedback = 'none';
try
    freq = ft_freqanalysis(cfg, dat);
catch
    return;
end

pow = freq.powspctrm;
if ndims(pow) == 3
    pow = squeeze(pow(:, 1, :));
elseif isvector(pow)
    pow = pow(:)';
end
if isfield(freq, 'freq') && ~isempty(freq.freq)
    n_freq = numel(freq.freq);
else
    n_freq = numel(scan_freqs);
end
if isvector(pow) && numel(pow) == n_freq
    % Keep one-trial output as row [1 x nFreq].
    pow = reshape(pow, 1, []);
elseif size(pow, 1) == n_freq && size(pow, 2) ~= n_freq
    % Guard against squeeze() producing [nFreq x nTrials].
    pow = pow.';
end
if size(pow, 1) == 1 && sum(valid_rows) > 1
    pow = repmat(pow, sum(valid_rows), 1);
end

if size(pow, 2) ~= numel(scan_freqs) && isfield(freq, 'freq') && ~isempty(freq.freq)
    freq_axis = freq.freq(:)';
    pow_interp = nan(size(pow, 1), numel(scan_freqs));
    for ri = 1:size(pow, 1)
        row_pow = pow(ri, :);
        if numel(row_pow) ~= numel(freq_axis) && size(pow, 1) == numel(freq_axis) && size(pow, 2) == size(pow_interp, 1)
            row_pow = pow(:, ri).';
        end
        if numel(row_pow) ~= numel(freq_axis)
            continue;
        end
        pow_interp(ri, :) = interp1(freq_axis, row_pow, scan_freqs, 'linear', NaN);
    end
    pow = pow_interp;
end

p_scan(valid_rows, :) = pow;
end
