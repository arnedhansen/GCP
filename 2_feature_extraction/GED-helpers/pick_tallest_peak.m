function [peak_hz, peak_power] = pick_tallest_peak(y, x, smooth_n, peak_power_halfwidth_hz)
% Peak frequency and local mean power around the tallest spectral peak.
peak_hz = NaN;
peak_power = NaN;
y = y(:)';
x = x(:)';
if ~isfinite(smooth_n) || smooth_n < 1
    smooth_n = 1;
end
if ~isfinite(peak_power_halfwidth_hz) || peak_power_halfwidth_hz < 0
    peak_power_halfwidth_hz = 0;
end
if isempty(y) || numel(y) ~= numel(x)
    return;
end
y = movmean(y, max(1, round(smooth_n)), 'omitnan');
core_mask = x >= 30 & x <= 90 & isfinite(y) & isfinite(x);
if any(core_mask)
    x_use = x(core_mask);
    y_use = y(core_mask);
else
    valid = isfinite(y) & isfinite(x);
    if ~any(valid)
        return;
    end
    x_use = x(valid);
    y_use = y(valid);
end
if numel(x_use) < 3
    return;
end

% Use local maxima only; if none pass criteria, keep NaN.
dx = diff(x_use);
dx = dx(isfinite(dx) & dx > 0);
if isempty(dx)
    return;
end
freq_step = median(dx);
if ~isfinite(freq_step) || freq_step <= 0
    return;
end
min_peak_width_hz = max(2 * freq_step, 2.0);
local_spread = iqr(y_use);
if ~isfinite(local_spread) || local_spread <= 0
    local_spread = std(y_use, 'omitnan');
end
if ~isfinite(local_spread) || local_spread <= 0
    local_spread = max(abs(y_use), [], 'omitnan');
end
if ~isfinite(local_spread) || local_spread <= 0
    local_spread = 1;
end
min_peak_prom_db = max(0.15, 0.25 * local_spread);

[pks, locs] = findpeaks(y_use, x_use, ...
    'MinPeakProminence', min_peak_prom_db, ...
    'MinPeakWidth', min_peak_width_hz, ...
    'SortStr', 'descend', ...
    'NPeaks', 1);
if isempty(pks) || isempty(locs) || ~isfinite(pks(1)) || ~isfinite(locs(1))
    return;
end

peak_hz = locs(1);
peak_power = pks(1);
if peak_power_halfwidth_hz > 0
    band_mask = abs(x_use - peak_hz) <= peak_power_halfwidth_hz;
    band_power = y_use(band_mask);
    band_power = band_power(isfinite(band_power));
    if ~isempty(band_power)
        peak_power = mean(band_power);
    end
end

end
