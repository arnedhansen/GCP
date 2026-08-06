function m = robust_trial_mean(x)
% Trimmed mean across trials; falls back to the ordinary mean for small n.
x = x(:);
x = x(isfinite(x));
if isempty(x)
    m = NaN;
    return;
end
if numel(x) < 4
    m = mean(x);
    return;
end
q = quantile(x, [0.25 0.75]);
iqr_val = q(2) - q(1);
if ~isfinite(iqr_val) || iqr_val <= eps
    m = median(x);
    return;
end
lo_thr = q(1) - 1.5 * iqr_val;
hi_thr = q(2) + 1.5 * iqr_val;
x_trim = x(x >= lo_thr & x <= hi_thr);
if isempty(x_trim)
    m = median(x);
else
    m = mean(x_trim);
end
end
