function avg_curve = compute_condition_average_powratio_ft(pr_mat, scan_freqs)
% Robust mean power-ratio spectrum across trials for one condition.
avg_curve = nan(1, numel(scan_freqs));
if isempty(pr_mat) || isempty(scan_freqs)
    return;
end
if size(pr_mat, 2) ~= numel(scan_freqs)
    return;
end
valid_trials = any(isfinite(pr_mat), 2);
if ~any(valid_trials)
    return;
end
avg_curve = mean(pr_mat(valid_trials, :), 1, 'omitnan');
end
