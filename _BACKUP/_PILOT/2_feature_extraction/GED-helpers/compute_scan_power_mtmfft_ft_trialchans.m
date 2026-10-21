function pow = compute_scan_power_mtmfft_ft_trialchans(trial_cell, fs, scan_freqs, tapsmofrq_hz, pad_sec)
% Multitaper power for multi-channel trials: trial_cell{i} is [nChan x nTime].
% Returns pow [nTrials x nChan x nFreq]. Each channel is demeaned before FFT.
nTrials = numel(trial_cell);
nFreq = numel(scan_freqs);
pow = nan(nTrials, 0, nFreq);
if nTrials == 0 || fs <= 0 || isempty(scan_freqs) || ~isfinite(tapsmofrq_hz) || tapsmofrq_hz <= 0
    return;
end

nChan = size(trial_cell{1}, 1);
pow = nan(nTrials, nChan, nFreq);
valid = false(nTrials, 1);
dat = [];
dat.label = arrayfun(@(k) sprintf('GED%03d', k), 1:nChan, 'UniformOutput', false);
dat.fsample = fs;
dat.trial = {};
dat.time = {};
map_out = zeros(0, 1);
for tr = 1:nTrials
    x = double(trial_cell{tr});
    if isempty(x) || size(x, 1) ~= nChan || size(x, 2) < 8
        continue;
    end
    if any(~isfinite(x(:)))
        continue;
    end
    x = x - mean(x, 2);
    valid(tr) = true;
    dat.trial{end+1} = x; %#ok<AGROW>
    dat.time{end+1} = (0:(size(x, 2) - 1)) / fs; %#ok<AGROW>
    map_out(end+1, 1) = tr; %#ok<AGROW>
end
if ~any(valid)
    return;
end

nValid = numel(dat.trial);
dat.sampleinfo = zeros(nValid, 2);
sample_start = 1;
for ti = 1:nValid
    nSamp = size(dat.trial{ti}, 2);
    dat.sampleinfo(ti, :) = [sample_start, sample_start + nSamp - 1];
    sample_start = sample_start + nSamp;
end

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'dpss';
cfg.foi = scan_freqs;
cfg.tapsmofrq = tapsmofrq_hz;
if nargin >= 5 && isfinite(pad_sec) && pad_sec > 0
    cfg.pad = pad_sec;
else
    cfg.pad = 'nextpow2';
end
cfg.keeptrials = 'yes';
cfg.feedback = 'none';
try
    freq = ft_freqanalysis(cfg, dat);
catch
    return;
end

P = double(freq.powspctrm);
if ndims(P) < 3
    if isvector(P)
        P = reshape(P, 1, 1, []);
    else
        P = reshape(P, size(P, 1), size(P, 2), 1);
    end
end
% Expect [nValid x nChan x nFreq]
if size(P, 1) ~= nValid && size(P, 2) == nValid
    P = permute(P, [2 1 3]);
end
if isfield(freq, 'freq') && ~isempty(freq.freq) && numel(freq.freq) ~= nFreq
    freq_axis = freq.freq(:)';
    P_interp = nan(size(P, 1), size(P, 2), nFreq);
    for ti = 1:size(P, 1)
        for ci = 1:size(P, 2)
            row = squeeze(P(ti, ci, :));
            if numel(row) ~= numel(freq_axis)
                continue;
            end
            P_interp(ti, ci, :) = interp1(freq_axis, row, scan_freqs, 'linear', NaN);
        end
    end
    P = P_interp;
end
if size(P, 2) ~= nChan
    return;
end
for vi = 1:numel(map_out)
    pow(map_out(vi), :, :) = P(vi, :, :);
end
end
