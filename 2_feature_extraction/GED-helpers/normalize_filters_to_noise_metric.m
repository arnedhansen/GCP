function Wn = normalize_filters_to_noise_metric(W, covBase)
% Scale GED filters so w'*covBase*w = 1.
Wn = W;
if isempty(W) || isempty(covBase)
    return;
end
for ci = 1:size(W, 2)
    w = W(:, ci);
    if ~all(isfinite(w))
        continue;
    end
    denom = real(w' * covBase * w);
    if ~isfinite(denom) || denom <= eps
        continue;
    end
    Wn(:, ci) = w / sqrt(denom);
end
end
