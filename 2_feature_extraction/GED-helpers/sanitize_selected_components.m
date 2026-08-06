function [sel_idx, sel_w] = sanitize_selected_components(sel_idx, sel_w, nComp)
% Validate selected component indices and renormalize combination weights.
if isempty(nComp) || ~isfinite(nComp)
    nComp = 0;
end
sel_idx = sel_idx(:);
sel_idx = sel_idx(isfinite(sel_idx));
sel_idx = round(sel_idx);
sel_idx = sel_idx(sel_idx >= 1 & sel_idx <= nComp);
if isempty(sel_idx)
    sel_w = [];
    return;
end
if isempty(sel_w) || numel(sel_w) ~= numel(sel_idx)
    sel_w = ones(numel(sel_idx), 1);
else
    sel_w = sel_w(:);
end
sel_w(~isfinite(sel_w) | sel_w <= 0) = 0;
if sum(sel_w) <= 0
    sel_w = ones(numel(sel_idx), 1);
end
sel_w = sel_w / sum(sel_w);
end
