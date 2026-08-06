function [base_floor, base_median] = compute_baseline_floor_stats(base_power_vals, floor_prctile, floor_frac)
% Estimate baseline power floor and median for near-floor instability checks.
base_floor = eps;
base_median = eps;
if ~isfinite(floor_prctile)
    floor_prctile = 20;
end
if ~isfinite(floor_frac)
    floor_frac = 0.25;
end
if isempty(base_power_vals)
    return;
end
valid = isfinite(base_power_vals) & (base_power_vals > 0);
if ~any(valid)
    return;
end
vals = base_power_vals(valid);
base_median = median(vals, 'omitnan');
anchor = prctile(vals, floor_prctile);
if ~isfinite(anchor) || anchor <= 0
    anchor = base_median;
end
if ~isfinite(base_median) || base_median <= 0
    base_median = anchor;
end
base_floor = max(anchor * max(floor_frac, 0), eps);
base_median = max(base_median, base_floor);
end
