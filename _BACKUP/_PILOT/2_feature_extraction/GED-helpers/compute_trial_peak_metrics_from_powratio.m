function [trl_peaks, trl_peak_power, trl_centroid] = ...
    compute_trial_peak_metrics_from_powratio(powratio_trials, scan_freqs, analysis_mask, ...
    smooth_n, peak_power_halfwidth_hz)
% Per-trial peak frequency, peak power, and spectral centroid from power-ratio scans.
nTrl = size(powratio_trials, 1);
trl_peaks = nan(nTrl, 1);
trl_peak_power = nan(nTrl, 1);
trl_centroid = nan(nTrl, 1);
scan_freqs_analysis = scan_freqs(analysis_mask);
centroid_band_mask = scan_freqs >= 30 & scan_freqs <= 90;
freq_band = scan_freqs(centroid_band_mask);
if ~isfinite(peak_power_halfwidth_hz) || peak_power_halfwidth_hz < 0
    peak_power_halfwidth_hz = 0;
end
for trl = 1:nTrl
    pr_scan = powratio_trials(trl, :);
    if all(~isfinite(pr_scan))
        continue;
    end
    pr_proc = pr_scan(analysis_mask);
    pr_proc = pr_proc(:)';
    valid = isfinite(pr_proc) & isfinite(scan_freqs_analysis);
    pr_proc = pr_proc(valid);
    x_use = scan_freqs_analysis(valid);
    if isempty(pr_proc)
        continue;
    end

    % Peak search over the 30-90 Hz gamma band.
    [peak_hz, peak_power] = pick_tallest_peak(pr_proc, x_use, smooth_n, peak_power_halfwidth_hz);
    if isfinite(peak_hz)
        trl_peaks(trl) = peak_hz;
        trl_peak_power(trl) = peak_power;
    end

    pr_smooth = movmean(pr_scan, max(1, round(smooth_n)), 'omitnan');
    pr_band = pr_smooth(centroid_band_mask);
    w_pos = max(pr_band, 0);
    pos_mass = sum(w_pos);
    if pos_mass > 0
        trl_centroid(trl) = sum(freq_band .* w_pos) / pos_mass;
    end
end
end
