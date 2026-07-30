%% GCP Oculo-Cortical Coherence (OCC), Contrasts Pooled
%
% Magnitude-squared coherence between eye velocity (velOCC)
% and EEG, with all contrast conditions (25/50/75/100%) pooled per subject.
% Purpose: test whether gamma-band coherence between cortical activity and
% oculomotor velocity is present, irrespective of contrast.
%
% Method:
%   1. Load per-subject velOCC, dataET (for sampleinfo), and EEG.
%   2. Match EEG↔ET trials by sampleinfo within each contrast; pool contrasts.
%   3. Multitaper Fourier (dpss, 3 Hz smoothing); magnitude-squared coherence.
%   4. Trial-shuffle Monte Carlo null (Vel2D trial order), kept as a full
%      surrogate distribution (not averaged before testing).
%   5. Subject-level Monte Carlo p-values; group tests by resampling one
%      surrogate draw per subject.
%
% Outputs:
%   - Grand-average occipital coherence spectrum (observed vs. null mean±SEM)
%   - Scalp topography of gamma-band OCC
%   - Subject-level observed vs. null with Monte Carlo p-values
%   - Group Monte Carlo / cluster test on occipital spectra

%% Setup
startup
[subjects, paths, ~, ~] = setup('GCP');
nSubj = numel(subjects);

condNames = {'c25','c50','c75','c100'};
nCond = numel(condNames);

foi_range   = [30 90];
tapsmofrq   = 3;
coh_latency = [0.3 2.0];
occ_pattern = '^(O|PO|I)';
nSurrogate  = 99;    % subject-level MC resolution: min p = 1/(nSurrogate+1)
nGroupPerm  = 1000;  % group MC by resampling stored subject surrogates

fig_dir = fullfile(paths.figures, 'coherence');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Preallocate
coh_spectra      = cell(1, nSubj);           % observed: chan x freq
surr_occ_spectra = cell(1, nSubj);           % nSurr x nFreq (occipital mean)
coh_gamma_obs    = nan(1, nSubj);
coh_gamma_surr   = nan(nSurrogate, nSubj);   % full null distribution
coh_gamma_peak   = nan(1, nSubj);
p_subj_mc        = nan(1, nSubj);

%% Subject loop
for subj = 1:nSubj
    fprintf('Subject %d/%d: %s\n', subj, nSubj, subjects{subj});

    eeg_path = fullfile(paths.features, subjects{subj}, 'eeg', 'dataEEG.mat');
    if ~isfile(eeg_path)
        warning('No EEG data for %s, skipping.', subjects{subj});
        continue
    end
    eeg_data = load(eeg_path, 'dataEEG_c25','dataEEG_c50','dataEEG_c75','dataEEG_c100');

    gaze_path = fullfile(paths.features, subjects{subj}, 'gaze', 'gaze_velocity_timeseries.mat');
    if ~isfile(gaze_path)
        warning('No gaze velocity data for %s, skipping.', subjects{subj});
        continue
    end
    vel_data = load(gaze_path, 'velOCC_trials_c25','velOCC_trials_c50', ...
        'velOCC_trials_c75','velOCC_trials_c100');

    et_path = fullfile(paths.features, subjects{subj}, 'gaze', 'dataET.mat');
    if ~isfile(et_path)
        warning('No dataET.mat for %s (needed for sampleinfo matching), skipping.', subjects{subj});
        continue
    end
    et_data = load(et_path, 'dataET_c25','dataET_c50','dataET_c75','dataET_c100');

    eeg_structs = {eeg_data.dataEEG_c25, eeg_data.dataEEG_c50, ...
                   eeg_data.dataEEG_c75, eeg_data.dataEEG_c100};
    vel_structs = {vel_data.velOCC_trials_c25, vel_data.velOCC_trials_c50, ...
                   vel_data.velOCC_trials_c75, vel_data.velOCC_trials_c100};
    et_structs  = {et_data.dataET_c25, et_data.dataET_c50, ...
                   et_data.dataET_c75, et_data.dataET_c100};

    pooledEEG = [];
    pooledVel = [];
    for ci = 1:nCond
        dataEEG = eeg_structs{ci};
        dataVel = vel_structs{ci};
        dataET  = et_structs{ci};

        cfg_sel = [];
        cfg_sel.channel = {'Vel2D'};
        dataVel = ft_selectdata(cfg_sel, dataVel);

        % Match EEG to ET/velOCC by sampleinfo (trialinfo is condition-constant)
        [idxEEG, idxVel] = local_match_by_sampleinfo(dataEEG, dataET, dataVel, subjects{subj}, condNames{ci});
        if numel(idxEEG) < 1
            continue
        end

        cfg_eeg = [];
        cfg_eeg.trials = idxEEG;
        dataEEG = ft_selectdata(cfg_eeg, dataEEG);

        cfg_vel = [];
        cfg_vel.trials = idxVel;
        dataVel = ft_selectdata(cfg_vel, dataVel);

        cfg_lat = [];
        cfg_lat.latency = coh_latency;
        dataEEG_cut = ft_selectdata(cfg_lat, dataEEG);
        dataVel_cut = ft_selectdata(cfg_lat, dataVel);

        if isempty(pooledEEG)
            pooledEEG = dataEEG_cut;
            pooledVel = dataVel_cut;
        else
            cfg_app = [];
            cfg_app.keepsampleinfo = 'no';
            pooledEEG = ft_appenddata(cfg_app, pooledEEG, dataEEG_cut);
            pooledVel = ft_appenddata(cfg_app, pooledVel, dataVel_cut);
        end
    end

    if isempty(pooledEEG)
        warning('%s: no pooled trials, skipping.', subjects{subj});
        continue
    end

    nTrl = numel(pooledEEG.trial);
    if nTrl ~= numel(pooledVel.trial)
        warning('%s: pooled EEG/Vel trial count mismatch (%d vs %d), skipping.', ...
            subjects{subj}, nTrl, numel(pooledVel.trial));
        continue
    end
    if nTrl < 20
        warning('%s: only %d pooled trials, skipping.', subjects{subj}, nTrl);
        continue
    end

    % Append Vel2D as extra channel
    dataCombined = pooledEEG;
    vel_label = pooledVel.label{1};
    dataCombined.label{end+1} = vel_label;
    for trl = 1:nTrl
        nSampEEG = size(dataCombined.trial{trl}, 2);
        nSampVel = numel(pooledVel.trial{trl});
        nSamp = min(nSampEEG, nSampVel);
        dataCombined.trial{trl} = [dataCombined.trial{trl}(:,1:nSamp); ...
                                   pooledVel.trial{trl}(1:nSamp)];
        dataCombined.time{trl} = dataCombined.time{trl}(1:nSamp);
    end
    vel_chan = numel(dataCombined.label);

    % Observed coherence
    coh = local_compute_coh(dataCombined, foi_range, tapsmofrq, vel_label);
    if isempty(coh)
        warning('%s: coherence failed, skipping.', subjects{subj});
        continue
    end
    coh_spectra{subj} = coh;

    occ_mask = cellfun(@(l) ~isempty(regexp(l, occ_pattern, 'once')), coh.label);
    if ~any(occ_mask)
        warning('%s: no occipital channels, skipping.', subjects{subj});
        continue
    end
    occ_obs = mean(coh.cohspctrm(occ_mask, :), 1, 'omitnan');
    coh_gamma_obs(subj) = mean(occ_obs, 'omitnan');
    [~, pk_idx] = max(occ_obs);
    coh_gamma_peak(subj) = coh.freq(pk_idx);

    % Surrogate distribution: shuffle Vel2D trial order (keep each draw)
    surr_occ = nan(nSurrogate, numel(coh.freq));
    fprintf('  Surrogates (%d)...\n', nSurrogate);
    for s = 1:nSurrogate
        dataSurr = dataCombined;
        shuf = randperm(nTrl);
        for trl = 1:nTrl
            nSamp = size(dataSurr.trial{trl}, 2);
            vel_src = pooledVel.trial{shuf(trl)};
            nSampVel = numel(vel_src);
            nUse = min(nSamp, nSampVel);
            row = dataSurr.trial{trl};
            row(vel_chan, 1:nUse) = vel_src(1:nUse);
            if nUse < nSamp
                row(vel_chan, nUse+1:nSamp) = 0;
            end
            dataSurr.trial{trl} = row;
        end
        coh_s = local_compute_coh(dataSurr, foi_range, tapsmofrq, vel_label);
        if isempty(coh_s), continue; end
        occ_s = mean(coh_s.cohspctrm(occ_mask, :), 1, 'omitnan');
        surr_occ(s,:) = occ_s;
        coh_gamma_surr(s, subj) = mean(occ_s, 'omitnan');
    end
    surr_occ_spectra{subj} = surr_occ;

    % Subject-level Monte Carlo: P(null >= observed)
    nOK = sum(isfinite(coh_gamma_surr(:, subj)));
    if nOK > 0 && isfinite(coh_gamma_obs(subj))
        p_subj_mc(subj) = (1 + sum(coh_gamma_surr(:, subj) >= coh_gamma_obs(subj), 'omitnan')) ...
            / (nOK + 1);
    end
end

fprintf('Pooled coherence computation complete.\n');

%% Grand-average occipital spectra
valid_subj = find(~cellfun(@isempty, coh_spectra) & ~cellfun(@isempty, surr_occ_spectra));
if isempty(valid_subj)
    error('No valid subjects for pooled coherence.');
end

freq_axis = coh_spectra{valid_subj(1)}.freq;
nFreq = numel(freq_axis);
nValid = numel(valid_subj);

all_obs = nan(nValid, nFreq);
all_surr_draws = nan(nValid, nSurrogate, nFreq); % keep full null
for si = 1:nValid
    s = valid_subj(si);
    occ_mask = cellfun(@(l) ~isempty(regexp(l, occ_pattern, 'once')), ...
        coh_spectra{s}.label);
    all_obs(si,:) = mean(coh_spectra{s}.cohspctrm(occ_mask,:), 1, 'omitnan');
    all_surr_draws(si,:,:) = surr_occ_spectra{s};
end

% Null display: mean across surrogate draws, then across subjects
surr_subj_mean = squeeze(mean(all_surr_draws, 2, 'omitnan')); % nValid x nFreq
ga_obs   = mean(all_obs, 1, 'omitnan');
sem_obs  = std(all_obs, [], 1, 'omitnan') ./ sqrt(nValid);
ga_surr  = mean(surr_subj_mean, 1, 'omitnan');
sem_surr = std(surr_subj_mean, [], 1, 'omitnan') ./ sqrt(nValid);

%% Figure 1: Grand-average spectrum (observed vs. null)
fig1 = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
cObs  = [0.10 0.35 0.75];
cSurr = [0.55 0.55 0.55];

patch([freq_axis fliplr(freq_axis)], [ga_surr-sem_surr fliplr(ga_surr+sem_surr)], ...
    cSurr, 'FaceColor', 0.8*cSurr+0.2, 'EdgeColor', 'none', 'FaceAlpha', 0.35);
hS = plot(freq_axis, ga_surr, 'Color', cSurr, 'LineWidth', 2.5, 'LineStyle', '--');

patch([freq_axis fliplr(freq_axis)], [ga_obs-sem_obs fliplr(ga_obs+sem_obs)], ...
    cObs, 'FaceColor', 0.8*cObs+0.2, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hO = plot(freq_axis, ga_obs, 'Color', cObs, 'LineWidth', 2.5);

xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('Coherence (Vel2D \leftrightarrow Occipital EEG)', 'FontSize', 16);
title(sprintf('Pooled OCC (all contrasts): N = %d', nValid), 'FontSize', 20);
legend([hO hS], {'Observed', sprintf('Shuffle null (mean of %d)', nSurrogate)}, ...
    'Location', 'northeast', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim(foi_range); grid on; box on;
saveas(fig1, fullfile(fig_dir, 'GCP_OCC_pooled_spectra.png'));

%% Figure 2: Topography of mean gamma-band coherence
fig2 = figure('Position', [0 0 1512 982], 'Color', 'w');
ref_coh = coh_spectra{valid_subj(1)};
topo_data = [];
topo_data.label = ref_coh.label;
topo_data.freq = mean(freq_axis);
topo_data.dimord = 'chan_freq';

all_coh_topo = nan(nValid, numel(ref_coh.label));
for si = 1:nValid
    s = valid_subj(si);
    all_coh_topo(si,:) = mean(coh_spectra{s}.cohspctrm, 2, 'omitnan')';
end
topo_data.powspctrm = mean(all_coh_topo, 1, 'omitnan')';

cfg_topo = [];
cfg_topo.layout = 'EEG1010.lay';
cfg_topo.parameter = 'powspctrm';
cfg_topo.comment = 'no';
cfg_topo.colorbar = 'yes';
cfg_topo.figure = 'gca';
cfg_topo.style = 'straight';
cfg_topo.marker = 'off';
ft_topoplotER(cfg_topo, topo_data);
title('Scalp Topography of Pooled Gamma OCC', 'FontSize', 20);
cb = colorbar; cb.Label.String = 'Coherence'; cb.FontSize = 14;
saveas(fig2, fullfile(fig_dir, 'GCP_OCC_pooled_topography.png'));

%% Group Monte Carlo on mean occipital gamma (resample one surrogate per subject)
obs_v = coh_gamma_obs(valid_subj)';
surr_mat = coh_gamma_surr(:, valid_subj); % nSurr x nValid

T_obs = mean(obs_v, 'omitnan');
T_null = nan(nGroupPerm, 1);
for p = 1:nGroupPerm
    picked = nan(nValid, 1);
    for si = 1:nValid
        row = surr_mat(:, si);
        ok = find(isfinite(row));
        if isempty(ok), continue; end
        picked(si) = row(ok(randi(numel(ok))));
    end
    T_null(p) = mean(picked, 'omitnan');
end
p_group_mc = (1 + sum(T_null >= T_obs, 'omitnan')) / (sum(isfinite(T_null)) + 1);

% Subject MC summary (Stouffer combination of one-sided subject p-values)
p_valid = p_subj_mc(valid_subj);
p_valid = max(p_valid, 1/(nSurrogate+1)); % numerical floor
z_stouffer = sum(norminv(1 - p_valid)) / sqrt(numel(p_valid));
p_stouffer = 1 - normcdf(z_stouffer);

%% Figure 3: Subject observed vs. null distribution mean, with MC p
fig3 = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

surr_mean_subj = mean(surr_mat, 1, 'omitnan')';
valid_pair = isfinite(obs_v) & isfinite(surr_mean_subj);
obs_plot  = obs_v(valid_pair);
surr_plot = surr_mean_subj(valid_pair);
nPair = numel(obs_plot);

bar(1, mean(surr_plot), 0.5, 'FaceColor', 0.8*cSurr+0.2, 'EdgeColor', cSurr, 'LineWidth', 1.5, 'FaceAlpha', 0.35);
bar(2, mean(obs_plot),  0.5, 'FaceColor', 0.8*cObs+0.2,  'EdgeColor', cObs,  'LineWidth', 1.5, 'FaceAlpha', 0.35);
errorbar(1, mean(surr_plot), std(surr_plot)/sqrt(nPair), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
errorbar(2, mean(obs_plot),  std(obs_plot)/sqrt(nPair),  'k', 'LineStyle', 'none', 'LineWidth', 1.5);

jit = (rand(nPair,1)-0.5)*0.15;
for i = 1:nPair
    plot([1 2]+jit(i), [surr_plot(i) obs_plot(i)], '-', 'Color', [0 0 0 0.25], 'LineWidth', 0.8);
end
scatter(1+jit, surr_plot, 40, cSurr, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'w');
scatter(2+jit, obs_plot,  40, cObs,  'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'w');

set(gca, 'XTick', [1 2], 'XTickLabel', {'Null (subj mean)', 'Observed'}, 'FontSize', 16);
ylabel('Mean Occipital Coherence (30-90 Hz)', 'FontSize', 16);
title(sprintf(['Pooled OCC vs. Null (group MC p = %.4g; Stouffer p = %.4g; N = %d)'], ...
    p_group_mc, p_stouffer, nPair), 'FontSize', 18);
xlim([0.4 2.6]); grid on; box on;
saveas(fig3, fullfile(fig_dir, 'GCP_OCC_pooled_obs_vs_null.png'));

%% Group Monte Carlo cluster test on occipital spectra
% Observed group-mean spectrum vs. null group-mean spectra formed by
% drawing one surrogate spectrum per subject on each permutation.
T_obs_spec = ga_obs;
T_null_spec = nan(nGroupPerm, nFreq);
for p = 1:nGroupPerm
    picked = nan(nValid, nFreq);
    for si = 1:nValid
        ok = find(all(isfinite(squeeze(all_surr_draws(si,:,:))), 2));
        if isempty(ok), continue; end
        r = ok(randi(numel(ok)));
        picked(si,:) = squeeze(all_surr_draws(si, r, :))';
    end
    T_null_spec(p,:) = mean(picked, 1, 'omitnan');
end

% Excess over null mean; cluster on frequencies where excess > 0
null_mean_spec = mean(T_null_spec, 1, 'omitnan');
null_std_spec  = std(T_null_spec, [], 1, 'omitnan');
null_std_spec(null_std_spec == 0) = inf;
obs_z = (T_obs_spec - null_mean_spec) ./ null_std_spec;
null_z = (T_null_spec - null_mean_spec) ./ null_std_spec;

z_thresh = norminv(0.95); % one-sided cluster-forming threshold
[obs_clusters, obs_mass] = local_freq_clusters(obs_z, z_thresh);
max_null_mass = zeros(nGroupPerm, 1);
for p = 1:nGroupPerm
    [~, masses] = local_freq_clusters(null_z(p,:), z_thresh);
    if isempty(masses)
        max_null_mass(p) = 0;
    else
        max_null_mass(p) = max(masses);
    end
end

cluster_p = nan(numel(obs_mass), 1);
for c = 1:numel(obs_mass)
    cluster_p(c) = (1 + sum(max_null_mass >= obs_mass(c))) / (nGroupPerm + 1);
end
mask_sig = false(1, nFreq);
for c = 1:numel(obs_clusters)
    if cluster_p(c) < 0.05
        mask_sig(obs_clusters{c}) = true;
    end
end

%% Figure 4: Difference spectrum with significant clusters
fig4 = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
diff_spec = all_obs - surr_subj_mean;
d_mean = mean(diff_spec, 1, 'omitnan');
d_sem  = std(diff_spec, [], 1, 'omitnan') ./ sqrt(nValid);
patch([freq_axis fliplr(freq_axis)], [d_mean-d_sem fliplr(d_mean+d_sem)], ...
    cObs, 'FaceColor', 0.8*cObs+0.2, 'EdgeColor', 'none', 'FaceAlpha', 0.35);
plot(freq_axis, d_mean, 'Color', cObs, 'LineWidth', 2.5);
yline(0, 'k--', 'LineWidth', 1);

if any(mask_sig)
    yl = ylim;
    sig_runs = diff([false mask_sig false]);
    starts = find(sig_runs == 1);
    stops  = find(sig_runs == -1) - 1;
    for r = 1:numel(starts)
        patch([freq_axis(starts(r)) freq_axis(stops(r)) ...
               freq_axis(stops(r)) freq_axis(starts(r))], ...
            [yl(1) yl(1) yl(2) yl(2)], [0 0 0], ...
            'FaceAlpha', 0.12, 'EdgeColor', 'none');
    end
end

xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('\Delta Coherence (Observed − Null mean)', 'FontSize', 16);
if any(isfinite(cluster_p))
    title(sprintf('Pooled OCC excess over null (min cluster p = %.4g)', ...
        min(cluster_p)), 'FontSize', 20);
else
    title('Pooled OCC excess over null (no clusters)', 'FontSize', 20);
end
set(gca, 'FontSize', 14); xlim(foi_range); grid on; box on;
saveas(fig4, fullfile(fig_dir, 'GCP_OCC_pooled_diff_cluster.png'));

%% Save
stat_coh = struct();
stat_coh.freq = freq_axis;
stat_coh.obs_z = obs_z;
stat_coh.mask = mask_sig;
stat_coh.cluster_p = cluster_p;
stat_coh.obs_clusters = obs_clusters;
stat_coh.p_group_mc = p_group_mc;
stat_coh.p_stouffer = p_stouffer;
stat_coh.p_subj_mc = p_subj_mc;

save(fullfile(paths.data, 'features', 'GCP_OCC_coherence_pooled.mat'), ...
    'coh_spectra', 'surr_occ_spectra', 'coh_gamma_obs', 'coh_gamma_surr', ...
    'coh_gamma_peak', 'p_subj_mc', 'ga_obs', 'sem_obs', 'ga_surr', 'sem_surr', ...
    'freq_axis', 'stat_coh', 'p_group_mc', 'p_stouffer', ...
    'coh_latency', 'foi_range', 'tapsmofrq', 'nSurrogate', 'nGroupPerm', '-v7.3');

fprintf('Group MC (mean gamma OCC): p = %.4g\n', p_group_mc);
fprintf('Stouffer combined subject MC: p = %.4g\n', p_stouffer);
fprintf('All figures saved to: %s\n', fig_dir);
fprintf('Data saved to: %s\n', fullfile(paths.data, 'features', 'GCP_OCC_coherence_pooled.mat'));

%% Local helpers
function [idxEEG, idxVel] = local_match_by_sampleinfo(dataEEG, dataET, dataVel, subjLabel, condLabel)
% Match EEG trials to velOCC via dataET sampleinfo.
% velOCC is built 1:1 with dataET trial order in GCP_gaze_fex.
    idxEEG = [];
    idxVel = [];

    nVel = numel(dataVel.trial);
    nET  = numel(dataET.trial);
    if nVel ~= nET
        warning('%s %s: velOCC (%d) and dataET (%d) trial counts differ.', ...
            subjLabel, condLabel, nVel, nET);
        nCommon = min(nVel, nET);
    else
        nCommon = nVel;
    end

    if ~isfield(dataEEG, 'sampleinfo') || isempty(dataEEG.sampleinfo) ...
            || ~isfield(dataET, 'sampleinfo') || isempty(dataET.sampleinfo)
        warning(['%s %s: sampleinfo missing; falling back to ordinal match ' ...
            'after verifying equal trial counts.'], subjLabel, condLabel);
        nEEG = numel(dataEEG.trial);
        if nEEG ~= nCommon
            warning('%s %s: cannot match without sampleinfo (%d EEG vs %d ET).', ...
                subjLabel, condLabel, nEEG, nCommon);
            return
        end
        idxEEG = (1:nEEG)';
        idxVel = (1:nCommon)';
        return
    end

    eegStarts = dataEEG.sampleinfo(:,1);
    etStarts  = dataET.sampleinfo(1:nCommon, 1);
    [~, idxEEG, idxET] = intersect(eegStarts, etStarts, 'stable');
    idxVel = idxET;

    if isempty(idxEEG)
        warning('%s %s: no sampleinfo overlap between EEG and ET.', subjLabel, condLabel);
        return
    end
    if numel(idxEEG) < numel(eegStarts) || numel(idxVel) < nCommon
        fprintf('  %s %s: kept %d/%d EEG and %d/%d ET trials via sampleinfo.\n', ...
            subjLabel, condLabel, numel(idxEEG), numel(eegStarts), numel(idxVel), nCommon);
    end
end

function coh = local_compute_coh(dataCombined, foi_range, tapsmofrq, vel_label)
    coh = [];
    cfg_freq = [];
    cfg_freq.method     = 'mtmfft';
    cfg_freq.output     = 'fourier';
    cfg_freq.foilim     = foi_range;
    cfg_freq.tapsmofrq  = tapsmofrq;
    cfg_freq.keeptrials = 'yes';
    cfg_freq.channel    = 'all';
    cfg_freq.pad        = 'nextpow2';
    freq = ft_freqanalysis(cfg_freq, dataCombined);

    cfg_coh = [];
    cfg_coh.method = 'coh';
    coh_full = ft_connectivityanalysis(cfg_coh, freq);

    vel_idx = find(strcmp(coh_full.label, vel_label));
    eeg_idx = setdiff(1:numel(coh_full.label), vel_idx);
    if isempty(vel_idx) || isempty(eeg_idx)
        return
    end

    coh = [];
    coh.label  = coh_full.label(eeg_idx);
    coh.freq   = coh_full.freq;
    coh.dimord = 'chan_freq';
    coh.cohspctrm = squeeze(coh_full.cohspctrm(eeg_idx, vel_idx, :));
    if size(coh.cohspctrm, 2) ~= numel(coh.freq)
        coh.cohspctrm = coh.cohspctrm';
    end
end

function [clusters, masses] = local_freq_clusters(z, z_thresh)
% Positive one-sided contiguous clusters above z_thresh; mass = sum(z).
    above = z > z_thresh & isfinite(z);
    clusters = {};
    masses = [];
    if ~any(above), return; end
    d = diff([false, above, false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;
    for i = 1:numel(starts)
        idx = starts(i):stops(i);
        clusters{end+1} = idx; %#ok<AGROW>
        masses(end+1) = sum(z(idx)); %#ok<AGROW>
    end
end
