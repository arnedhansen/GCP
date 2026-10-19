%% GED eligibility with and without the Laplace peak-form gate
% Same 30 to 90 Hz covariance and the same other gates as GCP_eeg_fex_GED.m.
% Does not overwrite GCP_eeg_GED.mat.

startup
[subjects, paths, colors, headmodel] = setup('GCP'); %#ok<ASGLU>
addpath('/Users/Arne/Documents/GitHub/GCP/2_feature_extraction/GED-helpers');
global ft_default
ft_default.showcallinfo = 'no';

out_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/controls';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
set(0, 'DefaultFigureVisible', 'off');

baseline_window = [-1.5, -0.5];
stim_window = [0, 2.0];
band = [30 90];
scan_freqs = 30:1:90;
lambda = 0.05;
min_eigval = 1.1;
min_pf = 0.80;
condCodes = [61 62 63 64];
condNames = {'c25', 'c50', 'c75', 'c100'};

nSubj = numel(subjects);
nSearch = 10;
nFreq = numel(scan_freqs);
elig_pf = false(nSubj, nSearch);
elig_nopf = false(nSubj, nSearch);
only_pf = false(nSubj, nSearch);
lam = nan(nSubj, nSearch);
pf = nan(nSubj, nSearch);
post_front = nan(nSubj, nSearch);
post_temp = nan(nSubj, nSearch);
topo_frac = nan(nSubj, nSearch);
hf = nan(nSubj, nSearch);
spec = nan(nSubj, nSearch, nFreq);
spec_pf = nan(nSubj, nFreq);
spec_nopf = nan(nSubj, nFreq);
peak_pf = nan(nSubj, 1);
peak_nopf = nan(nSubj, 1);
ex50_pf = nan(nSubj, 1);
ex50_nopf = nan(nSubj, 1);
n_pf = nan(nSubj, 1);
n_nopf = nan(nSubj, 1);

for s = 1:nSubj
    t_subj = tic;
    fprintf('[NO PF] Loading %s\n', subjects{s});
    eeg = load(fullfile(paths.features, subjects{s}, 'eeg', 'dataEEG.mat'), ...
        'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    dataStructs = {eeg.dataEEG_c25, eeg.dataEEG_c50, eeg.dataEEG_c75, eeg.dataEEG_c100};
    clear eeg
    dat_per_cond = cell(1, 4);
    for c = 1:4
        dat = dataStructs{c};
        if isempty(dat) || ~isfield(dat, 'trialinfo')
            continue
        end
        trlIdx = find(dat.trialinfo == condCodes(c));
        if isempty(trlIdx)
            continue
        end
        dat_per_cond{c} = select_trials_by_index(dat, trlIdx);
    end
    clear dataStructs
    out = ged_gates(dat_per_cond, band, baseline_window, stim_window, ...
        lambda, headmodel, scan_freqs, min_eigval, min_pf);
    nHere = size(out.elig_pf, 1);
    elig_pf(s, 1:nHere) = out.elig_pf';
    elig_nopf(s, 1:nHere) = out.elig_nopf';
    only_pf(s, 1:nHere) = out.only_pf';
    lam(s, 1:nHere) = out.lambda';
    pf(s, 1:nHere) = out.pf';
    post_front(s, 1:nHere) = out.post_front';
    post_temp(s, 1:nHere) = out.post_temp';
    topo_frac(s, 1:nHere) = out.topo_frac';
    hf(s, 1:nHere) = out.hf';
    spec(s, 1:nHere, :) = out.spectra;
    spec_pf(s, :) = out.spec_pf;
    spec_nopf(s, :) = out.spec_nopf;
    peak_pf(s) = out.peak_pf;
    peak_nopf(s) = out.peak_nopf;
    ex50_pf(s) = out.ex50_pf;
    ex50_nopf(s) = out.ex50_nopf;
    n_pf(s) = out.n_pf;
    n_nopf(s) = out.n_nopf;
    clear dat_per_cond
    fprintf('[NO PF] %s  with PF %d, without PF %d, added by dropping PF %d, %.1f min\n', ...
        subjects{s}, out.n_pf, out.n_nopf, sum(out.only_pf), toc(t_subj) / 60);
    added = find(out.only_pf);
    for k = 1:numel(added)
        ci = added(k);
        fprintf('  C%d  lambda %.2f  PF %.2f  post/front %.2f  post/temp %.2f  topo %.2f  HF %.2f\n', ...
            ci, out.lambda(ci), out.pf(ci), out.post_front(ci), out.post_temp(ci), ...
            out.topo_frac(ci), out.hf(ci));
    end
end

fprintf('\n[NO PF] Median eligible with PF %.1f, without %.1f\n', ...
    median(n_pf, 'omitnan'), median(n_nopf, 'omitnan'));
fprintf('[NO PF] Median |peak shift| %.2f Hz\n', median(abs(peak_nopf - peak_pf), 'omitnan'));
fprintf('[NO PF] Median 50 Hz excess with PF %.2f dB, without %.2f dB\n', ...
    median(ex50_pf, 'omitnan'), median(ex50_nopf, 'omitnan'));

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
for s = 1:nSubj
    subplot(2, 5, s); hold on;
    plot(scan_freqs, spec_pf(s, :), 'Color', [0.15 0.15 0.15], 'LineWidth', 1.8);
    plot(scan_freqs, spec_nopf(s, :), 'Color', [0.85 0.35 0.1], 'LineWidth', 1.6);
    xline(50, ':', 'Color', [0.6 0.6 0.6]);
    xlim([30 90]);
    title(sprintf('%s  %d to %d', subjects{s}, n_pf(s), n_nopf(s)));
    box off;
    if s == 1
        ylabel('Stim / base [dB]');
        legend({'With PF', 'Without PF'}, 'Box', 'off', 'Location', 'northeast');
    end
    if s > 5
        xlabel('Frequency [Hz]');
    end
end
print(fig, fullfile(out_dir, 'GCP_control_ged_no_pf.png'), '-dpng', '-r300');
close(fig);

save(fullfile(out_dir, 'GCP_control_ged_no_pf.mat'), ...
    'subjects', 'elig_pf', 'elig_nopf', 'only_pf', 'lam', 'pf', ...
    'post_front', 'post_temp', 'topo_frac', 'hf', 'spec', ...
    'spec_pf', 'spec_nopf', 'peak_pf', 'peak_nopf', 'ex50_pf', 'ex50_nopf', ...
    'n_pf', 'n_nopf', 'scan_freqs', '-v7.3');
fprintf('[NO PF] Saved %s\n', out_dir);

function out = ged_gates(dat_per_cond, band, baseline_window, stim_window, ...
    lambda, headmodel, scan_freqs, min_eigval, min_pf)
out = struct();
ref = [];
for c = 1:numel(dat_per_cond)
    if ~isempty(dat_per_cond{c})
        ref = dat_per_cond{c};
        break
    end
end
labels = ref.label;
nChans = numel(labels);
fs = ref.fsample;
cfg_filt = [];
cfg_filt.bpfilter = 'yes';
cfg_filt.bpfreq = band;
cfg_filt.bpfilttype = 'fir';
cfg_filt.bpfiltord = round(3 * fs / band(1));
cfg_filt.feedback = 'none';
covStim = zeros(nChans);
covBase = zeros(nChans);
nTot = 0;
for c = 1:numel(dat_per_cond)
    dat = dat_per_cond{c};
    if isempty(dat) || ~isfield(dat, 'trial') || isempty(dat.trial)
        continue
    end
    dat_f = ft_preprocessing(cfg_filt, dat);
    nTrl = numel(dat_f.trial);
    if nTrl < 1
        continue
    end
    covBase = covBase + compute_pooled_covariance_window(dat_f, baseline_window) * nTrl;
    covStim = covStim + compute_pooled_covariance_window(dat_f, stim_window) * nTrl;
    nTot = nTot + nTrl;
    clear dat_f
end
covStim = covStim / nTot;
covBase = covBase / nTot;
covStim_reg = (1 - lambda) * covStim + lambda * mean(diag(covStim)) * eye(nChans);
covBase_reg = (1 - lambda) * covBase + lambda * mean(diag(covBase)) * eye(nChans);
[W, D] = eig(covStim_reg, covBase_reg);
[evals, ord] = sort(real(diag(D)), 'descend');
W = W(:, ord);
nSearch = min(10, size(W, 2));
[post_idx, front_idx, temp_idx] = region_indices(labels);
sim_template = occipital_template(labels, headmodel, post_idx, front_idx);
filters = nan(nChans, nSearch);
topos = nan(nChans, nSearch);
post_front = nan(nSearch, 1);
post_temp = nan(nSearch, 1);
topo_frac = nan(nSearch, 1);
for ci = 1:nSearch
    w = W(:, ci);
    topo = covStim_reg * w;
    r = corr(topo, sim_template, 'rows', 'complete');
    if ~isnan(r) && r < 0
        w = -w;
        topo = -topo;
    end
    filters(:, ci) = w;
    topos(:, ci) = topo;
    post_s = mean(abs(topo(post_idx)));
    front_s = mean(abs(topo(front_idx)));
    temp_s = mean(abs(topo(temp_idx)));
    post_front(ci) = post_s / max(front_s, eps);
    post_temp(ci) = post_s / max(temp_s, eps);
    topo_abs = abs(topo);
    topo_frac(ci) = max(topo_abs) / max(sum(topo_abs), eps);
end
proxies = estimate_components_artifact_proxies( ...
    filters, dat_per_cond, stim_window, baseline_window, ...
    fs, scan_freqs, 3, '[NO PF]');
spectra = nan(nSearch, numel(scan_freqs));
hf_slope = zeros(nSearch, 1);
for ci = 1:nSearch
    spectra(ci, :) = proxies(ci).mean_pr_spectrum(:)';
    if isfinite(proxies(ci).hf_slope)
        hf_slope(ci) = proxies(ci).hf_slope;
    end
end
[pf_score, ~] = compute_powspctrm_form_laplacian_score_from_spectra(spectra, scan_freqs, band);
finite_core = isfinite(evals(1:nSearch)) & isfinite(post_front) & ...
    isfinite(post_temp) & isfinite(topo_frac);
pass_eig = finite_core & (evals(1:nSearch) >= min_eigval);
pass_space = post_front > 1 & post_temp > 1 & topo_frac <= 0.25 & hf_slope <= 0;
pass_pf = isfinite(pf_score) & (pf_score >= min_pf);
elig_pf = pass_eig & pass_space & pass_pf;
elig_nopf = pass_eig & pass_space;
only_pf = elig_nopf & ~elig_pf;

out.elig_pf = elig_pf;
out.elig_nopf = elig_nopf;
out.only_pf = only_pf;
out.lambda = evals(1:nSearch);
out.pf = pf_score;
out.post_front = post_front;
out.post_temp = post_temp;
out.topo_frac = topo_frac;
out.hf = hf_slope;
out.spectra = spectra;
out.n_pf = sum(elig_pf);
out.n_nopf = sum(elig_nopf);
out.spec_pf = combine_spec(elig_pf, evals(1:nSearch), spectra);
out.spec_nopf = combine_spec(elig_nopf, evals(1:nSearch), spectra);
out.peak_pf = peak_hz_of(scan_freqs, out.spec_pf);
out.peak_nopf = peak_hz_of(scan_freqs, out.spec_nopf);
out.ex50_pf = band_excess_db(scan_freqs, out.spec_pf, 50);
out.ex50_nopf = band_excess_db(scan_freqs, out.spec_nopf, 50);
end

function spec = combine_spec(mask, evals, spectra)
if ~any(mask)
    spec = nan(1, size(spectra, 2));
    return
end
wts = evals(mask);
wts = wts(:)';
wts(~isfinite(wts) | wts <= 0) = 0;
if sum(wts) <= 0
    wts = ones(size(wts));
end
wts = wts / sum(wts);
spec = wts * spectra(mask, :);
end

function pk = peak_hz_of(freq, spec)
pk = NaN;
if isempty(spec) || all(~isfinite(spec))
    return
end
[pk, ~] = pick_tallest_peak(spec, freq, 0, 0);
if ~isfinite(pk)
    mask = freq >= 30 & freq <= 90 & isfinite(spec);
    yy = spec(mask);
    xx = freq(mask);
    [~, ix] = max(yy);
    if ~isempty(ix)
        pk = xx(ix);
    end
end
end

function [post_idx, front_idx, temp_idx] = region_indices(labels)
post_idx = find(cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P)', 'once')), labels));
front_idx = find(cellfun(@(l) ~isempty(regexp(l, '^(Fp|AF|F)', 'once')), labels));
temp_idx = find(cellfun(@(l) ~isempty(regexp(l, '^(T|TP|FT)', 'once')), labels));
end

function tmpl = occipital_template(labels, headmodel, post_idx, front_idx)
nChans = numel(labels);
tmpl = zeros(nChans, 1);
use_lay = isstruct(headmodel) && isfield(headmodel, 'layANThead') && ...
    isfield(headmodel.layANThead, 'label') && isfield(headmodel.layANThead, 'pos');
if use_lay
    lay_labels = headmodel.layANThead.label;
    lay_pos = headmodel.layANThead.pos;
    pos = nan(nChans, 2);
    for ch = 1:nChans
        li = find(strcmp(lay_labels, labels{ch}), 1, 'first');
        if ~isempty(li)
            pos(ch, :) = lay_pos(li, :);
        end
    end
    has = ~any(isnan(pos), 2);
    occ_idx = find(cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P10|P9)', 'once')), labels));
    if any(has(occ_idx)) && any(has(front_idx))
        occ_ctr = mean(pos(occ_idx(has(occ_idx)), :), 1);
        front_ctr = mean(pos(front_idx(has(front_idx)), :), 1);
        d_occ = sum((pos - occ_ctr).^2, 2);
        d_front = sum((pos - front_ctr).^2, 2);
        tmpl = exp(-d_occ / (2 * 0.12^2)) - 0.75 * exp(-d_front / (2 * 0.25^2));
        tmpl(~has) = 0;
    end
end
if all(tmpl == 0)
    tmpl(post_idx) = 1;
    tmpl(front_idx) = -0.75;
end
if std(tmpl) > 0
    tmpl = (tmpl - mean(tmpl)) / std(tmpl);
end
end

function ex = band_excess_db(freq, spec_db, f0)
ex = NaN;
if isempty(spec_db) || all(~isfinite(spec_db))
    return
end
i0 = find(abs(freq - f0) < 0.51, 1, 'first');
neigh = (abs(freq - f0) >= 4 & abs(freq - f0) <= 8);
if isempty(i0) || nnz(neigh) < 4 || ~isfinite(spec_db(i0))
    return
end
ex = spec_db(i0) - mean(spec_db(neigh), 'omitnan');
end
