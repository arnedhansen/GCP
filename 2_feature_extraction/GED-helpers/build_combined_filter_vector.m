function w_combined = build_combined_filter_vector(W_combined, w_components)
% Form one channel filter as the lambda-weighted sum of selected GED columns.
w_combined = [];
if isempty(W_combined)
    return;
end
nComp = size(W_combined, 2);
if nComp < 1
    return;
end
if isempty(w_components) || numel(w_components) ~= nComp
    w_components = ones(nComp, 1);
else
    w_components = w_components(:);
end
w_components(~isfinite(w_components) | w_components <= 0) = 0;
if sum(w_components) <= 0
    w_components = ones(nComp, 1);
end
w_components = w_components / sum(w_components);
w_combined = W_combined * w_components;
if ~all(isfinite(w_combined))
    w_combined = [];
end
end
