function add_peak_point_overlay(freq_axis, spec_vec, peak_hz)
% Mark the spectrum peak frequency on the current axes.
if isempty(freq_axis) || isempty(spec_vec) || numel(freq_axis) ~= numel(spec_vec)
    return;
end
valid = isfinite(freq_axis) & isfinite(spec_vec);
x = freq_axis(valid);
y = spec_vec(valid);
if numel(x) < 2 || numel(y) < 2
    return;
end
[x, sort_idx] = sort(x(:)');
y = y(sort_idx);
baseline_y = min(y);
if ~isfinite(baseline_y)
    return;
end

if ~isfinite(peak_hz)
    [~, idx_peak] = max(y);
    peak_hz = x(idx_peak);
end

peak_power = interp1(x, y, peak_hz, 'linear', NaN);
if ~isfinite(peak_power)
    [~, idx_near] = min(abs(x - peak_hz));
    peak_hz = x(idx_near);
    peak_power = y(idx_near);
end
if isfinite(peak_hz) && isfinite(peak_power)
    plot([peak_hz peak_hz], [baseline_y peak_power], '--', ...
        'Color', [0.55 0.55 0.55], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    plot(peak_hz, peak_power, 'o', ...
        'Color', [0.55 0.55 0.55], 'MarkerFaceColor', [0.55 0.55 0.55], ...
        'MarkerSize', 4, 'HandleVisibility', 'off');
end
end
