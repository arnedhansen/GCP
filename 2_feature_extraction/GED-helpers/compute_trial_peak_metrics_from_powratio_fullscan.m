function [trl_peaks, trl_peak_power, trl_centroid] = ...
    compute_trial_peak_metrics_from_powratio_fullscan(powratio_trials_fullscan, scan_freqs_full, analysis_mask, ...
    smooth_n, peak_power_halfwidth_hz)
% Per-trial peak frequency, peak power, and spectral centroid from power-ratio scans.
nTrl = size(powratio_trials_fullscan, 1);
trl_peaks = nan(nTrl, 1);
trl_peak_power = nan(nTrl, 1);
trl_centroid = nan(nTrl, 1);
scan_freqs_analysis = scan_freqs_full(analysis_mask);
centroid_band_mask = scan_freqs_full >= 30 & scan_freqs_full <= 90;
freq_band = scan_freqs_full(centroid_band_mask);
if ~isfinite(peak_power_halfwidth_hz) || peak_power_halfwidth_hz < 0
    peak_power_halfwidth_hz = 0;
end
for trl = 1:nTrl
    pr_full = powratio_trials_fullscan(trl, :);
    if all(~isfinite(pr_full))
        continue;
    end
    pr_proc = pr_full(analysis_mask);
    pr_proc = pr_proc(:)';
    valid = isfinite(pr_proc) & isfinite(scan_freqs_analysis);
    pr_proc = pr_proc(valid);
    x_use = scan_freqs_analysis(valid);
    if isempty(pr_proc)
        continue;
    end

    % Peak search over the full 30-90 Hz gamma band.
    [peak_hz, peak_power] = pick_tallest_peak(pr_proc, x_use, smooth_n, peak_power_halfwidth_hz);
    if isfinite(peak_hz)
        trl_peaks(trl) = peak_hz;
        trl_peak_power(trl) = peak_power;
    end

    pr_proc_full = movmean(pr_full, max(1, round(smooth_n)), 'omitnan');
    pr_proc_band = pr_proc_full(centroid_band_mask);
    w_pos = max(pr_proc_band, 0);
    pos_mass = sum(w_pos);
    if pos_mass > 0
        trl_centroid(trl) = sum(freq_band .* w_pos) / pos_mass;
    end
end
end
