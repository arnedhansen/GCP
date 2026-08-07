function C = compute_pooled_covariance_window(dat, time_window)
% Trial-pooled covariance over a latency window.
% Matches FieldTrip ft_timelockanalysis with covariance='yes', removemean='yes',
% keeptrials='no': per-trial demean, sum(X*X'), divide by sum(nTime-1).
nChan = numel(dat.label);
C = zeros(nChan);
dof = 0;
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
    nTime = sum(idx);
    if nTime < 2
        continue;
    end
    xw = x(:, idx);
    xw = xw - mean(xw, 2);
    C = C + (xw * xw');
    dof = dof + max(nTime - 1, 1);
end
if dof > 0
    C = C / dof;
end
end
