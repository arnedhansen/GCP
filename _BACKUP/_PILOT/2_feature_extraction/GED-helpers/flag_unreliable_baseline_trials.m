function bad_mask = flag_unreliable_baseline_trials(base_power_vals, mad_mult)
% Flag trials with abnormally low baseline power (MAD criterion).
bad_mask = false(size(base_power_vals));
if isempty(base_power_vals)
    return;
end
valid = isfinite(base_power_vals) & (base_power_vals > 0);
if ~any(valid)
    return;
end
logp = log10(base_power_vals(valid));
center = median(logp, 'omitnan');
spread = robust_mad(logp);
if ~isfinite(spread) || spread <= eps
    spread = iqr(logp) / 1.349;
end
if ~isfinite(spread) || spread <= eps
    return;
end
zabs = abs((logp - center) / spread);
tmp = false(size(logp));
tmp(zabs > mad_mult) = true;
bad_mask(valid) = tmp;
end
