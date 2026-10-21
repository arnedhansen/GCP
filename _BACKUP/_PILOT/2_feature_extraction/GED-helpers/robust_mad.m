function m = robust_mad(x)
% Median absolute deviation (normalized) of finite samples.
x = x(:);
x = x(isfinite(x));
if isempty(x)
    m = 0;
    return;
end
m = mad(x, 1);
if ~isfinite(m) || m <= eps
    m = std(x);
end
if ~isfinite(m) || m <= eps
    m = 0;
end
end
