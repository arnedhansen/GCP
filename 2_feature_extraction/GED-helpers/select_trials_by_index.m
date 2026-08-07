function dat_out = select_trials_by_index(dat_in, trl_idx)
% Lightweight trial subset without ft_selectdata copying overhead.
dat_out = dat_in;
dat_out.trial = dat_in.trial(trl_idx);
dat_out.time = dat_in.time(trl_idx);
if isfield(dat_in, 'trialinfo')
    dat_out.trialinfo = dat_in.trialinfo(trl_idx, :);
end
if isfield(dat_in, 'sampleinfo')
    dat_out.sampleinfo = dat_in.sampleinfo(trl_idx, :);
end
end
