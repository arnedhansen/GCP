function C = compute_pooled_covariance_window(dat, time_window)
% Trial-pooled covariance over a latency window (demean per trial, then X*X').
% Matches the intent of ft_timelockanalysis covariance + removemean='yes'.
nChan = numel(dat.label);
C = zeros(nChan);
nSamp_total = 0;
if isempty(dat) || ~isfield(dat, 'trial') || isempty(dat.trial)
    return;
end
for trl = 1:numel(dat.trial)
    x = double(dat.trial{trl});
    t = dat.time{trl};
    if isempty(x) || isempty(t)
        continue;
    end
    idx = t >= time_window(1) & t <= time_window(2);
    if sum(idx) < 2
        continue;
    end
    xw = x(:, idx);
    xw = xw - mean(xw, 2);
    C = C + (xw * xw');
    nSamp_total = nSamp_total + size(xw, 2);
end
if nSamp_total > 1
    C = C / (nSamp_total - 1);
end
end
