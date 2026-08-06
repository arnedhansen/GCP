function ylims = compute_symmetric_ylim(vals, pad_frac, min_half_range)
% Symmetric y-limits around zero from data range and padding.
if ~isfinite(pad_frac) || pad_frac < 0
    pad_frac = 0.15;
end
if ~isfinite(min_half_range) || min_half_range <= 0
    min_half_range = 1;
end
vals = vals(:);
vals = vals(isfinite(vals));
if isempty(vals)
    half_range = min_half_range;
else
    max_abs = max(abs(vals));
    if ~isfinite(max_abs)
        max_abs = min_half_range;
    end
    half_range = max(max_abs * (1 + pad_frac), min_half_range);
end
ylims = [-half_range, half_range];
end
