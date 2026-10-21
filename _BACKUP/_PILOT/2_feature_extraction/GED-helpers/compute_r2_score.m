function r2 = compute_r2_score(y, yhat)
% Coefficient of determination between observed and predicted values.
if numel(y) ~= numel(yhat)
    r2 = 0;
    return;
end
valid = isfinite(y) & isfinite(yhat);
y = y(valid);
yhat = yhat(valid);
if numel(y) < 3
    r2 = 0;
    return;
end
sse = sum((y - yhat).^2);
sst = sum((y - mean(y)).^2);
if ~isfinite(sst) || sst <= eps
    r2 = 0;
else
    r2 = 1 - (sse / sst);
end
r2 = max(0, min(1, r2));
end
