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
freq_dat = [];
freq_dat.label = {'GED'};
freq_dat.freq = scan_freqs(:)';
freq_dat.dimord = 'rpt_chan_freq';
freq_dat.powspctrm = nan(sum(valid_trials), 1, numel(scan_freqs));
trial_rows = find(valid_trials);
for ti = 1:numel(trial_rows)
    freq_dat.powspctrm(ti, 1, :) = pr_mat(trial_rows(ti), :);
end
cfg = [];
cfg.avgoverrpt = 'yes';
try
    freq_avg = ft_selectdata(cfg, freq_dat);
catch
    avg_curve = mean(pr_mat(valid_trials, :), 1, 'omitnan');
    return;
end
pow_avg = freq_avg.powspctrm;
if ndims(pow_avg) == 3
    pow_avg = squeeze(pow_avg(1, 1, :));
elseif ismatrix(pow_avg)
    pow_avg = squeeze(pow_avg);
end
avg_curve = pow_avg(:)';
if numel(avg_curve) ~= numel(scan_freqs)
    avg_curve = mean(pr_mat(valid_trials, :), 1, 'omitnan');
end
if all(~isfinite(avg_curve))
    avg_curve = mean(pr_mat(valid_trials, :), 1, 'omitnan');
end
end
