function [p_stim, p_base] = compute_scan_power_mtmfft_ft_pair_trialchans(stim_trials, base_trials, fs, scan_freqs, tapsmofrq_hz)
% Stim/baseline mtmfft for multi-channel trials with a shared pad.
% stim_trials{i} and base_trials{i} are [nChan x nTime].
% Returns p_stim/p_base as [nTrials x nChan x nFreq].
% Pad is nextpow2 while that bin is at least as fine as the requested step.
% Shorter data (for example a 300 ms snippet) is padded so the FFT bin
% equals that step, which keeps 1 Hz bins from collapsing onto each other.
nTrials = numel(stim_trials);
nFreq = numel(scan_freqs);
p_stim = nan(nTrials, 0, nFreq);
p_base = nan(nTrials, 0, nFreq);
if nTrials == 0 || numel(base_trials) ~= nTrials || fs <= 0
    return;
end
nChan = size(stim_trials{1}, 1);
p_stim = nan(nTrials, nChan, nFreq);
p_base = nan(nTrials, nChan, nFreq);

max_len = 0;
for tr = 1:nTrials
    max_len = max([max_len, size(stim_trials{tr}, 2), size(base_trials{tr}, 2)]);
end
if max_len < 8
    return;
end
nfft = 2 ^ ceil(log2(max_len));
if numel(scan_freqs) >= 2
    foi_step = min(diff(sort(scan_freqs(:))));
    if isfinite(foi_step) && foi_step > 0 && fs / nfft > foi_step
        nfft_step = round(fs / foi_step);
        if nfft_step >= max_len && abs(fs / nfft_step - foi_step) <= 1e-6 * foi_step
            nfft = nfft_step;
        else
            while fs / nfft > foi_step
                nfft = nfft + 1;
            end
        end
    end
end
pad_sec = nfft / fs;
p_stim = compute_scan_power_mtmfft_ft_trialchans(stim_trials, fs, scan_freqs, tapsmofrq_hz, pad_sec);
p_base = compute_scan_power_mtmfft_ft_trialchans(base_trials, fs, scan_freqs, tapsmofrq_hz, pad_sec);
end
