%% GCP standalone GED + PF comparison
%
% Scope:
% - Full-window GED only (same setup as main GED script)
% - Component eigenvalue ranking
% - PF metrics (0-1):
%   (1) Current PF (full current implementation)
%   (2) Kurtosis-based PF
%   (3) Laplace fit PF (R^2)
%   (4) Truncated Gaussian fit PF (R^2)
% - Visualizations:
%   (a) Scree + PF curves
%   (b) Component topographies + spectra with all PF values printed

%% Setup
startup
[subjects, paths, ~, headmodel] = setup('GCP');
nSubj = numel(subjects);

%% Parameters
baseline_window = [-1.5, -0.25];
full_window = [0, 2.0];
analysis_freq_range = [30 90];
scan_freq_step_hz = 1;
scan_freqs = analysis_freq_range(1):scan_freq_step_hz:analysis_freq_range(2);
scan_width = 2.0;

lambda = 0.05;      % regularization (same as main script)
ged_search_n = 10;  % search first N GED components
random_seed = 123;

condNames = {'c25', 'c50', 'c75', 'c100'};
trialCodes = [61, 62, 63, 64];

fig_save_dir = fullfile(paths.figures, 'eeg', 'ged', 'pf_computation');
if ~exist(fig_save_dir, 'dir')
    mkdir(fig_save_dir);
end

gcp_feature_data_path = paths.features;
if ~exist(gcp_feature_data_path, 'dir')
    gcp_feature_data_path = paths.root;
end

%% Subject loop
for subj = 1:nSubj
    subject_id = subjects{subj};
    subj_tag = ['subj', subject_id];
    clc
    fprintf('\n[Powspec Form Computation] Subject %s (%d/%d)\n', subject_id, subj, nSubj);
    rng(random_seed + subj, 'twister');

    datapath = fullfile(gcp_feature_data_path, subject_id, 'eeg');
    data_file = fullfile(datapath, 'dataEEG.mat');
    if ~exist(data_file, 'file')
        warning('Data file not found: %s', data_file);
        continue;
    end

    eeg_data = load(data_file, 'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    dataStructs = {eeg_data.dataEEG_c25, eeg_data.dataEEG_c50, eeg_data.dataEEG_c75, eeg_data.dataEEG_c100};
    fsample = dataStructs{1}.fsample;
    nChans = numel(dataStructs{1}.label);
    nSearch = min(ged_search_n, nChans);

    trialIndices = cell(1, 4);
    for cond = 1:4
        trialIndices{cond} = find(dataStructs{cond}.trialinfo == trialCodes(cond));
    end

    % Channel groups for signed template orientation
    occ_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P10|P9)', 'once')), dataStructs{1}.label);
    occ_idx  = find(occ_mask);
    front_mask = cellfun(@(l) ~isempty(regexp(l, '^(Fp|AF|F)', 'once')), dataStructs{1}.label);
    front_idx  = find(front_mask);

    %% Build pooled covariance (full window only) exactly as in main GED setup
    covStim_full = zeros(nChans);
    covBase_full = zeros(nChans);
    nTrials_total = 0;
    dat_per_cond = cell(1, 4);

    for cond = 1:4
        dat = dataStructs{cond};
        trlIdx = trialIndices{cond};
        if isempty(trlIdx)
            continue;
        end

        cfg = [];
        cfg.trials = trlIdx;
        dat = ft_selectdata(cfg, dat);
        dat_per_cond{cond} = dat;

        cfg_filt = [];
        cfg_filt.bpfilter   = 'yes';
        cfg_filt.bpfreq     = analysis_freq_range;
        cfg_filt.bpfilttype = 'fir';
        cfg_filt.bpfiltord  = round(3 * fsample / analysis_freq_range(1));
        dat_gamma = ft_preprocessing(cfg_filt, dat);

        cfg_t = [];
        cfg_t.latency = baseline_window;
        dat_base = ft_selectdata(cfg_t, dat_gamma);

        cfg_t.latency = full_window;
        dat_stim_full = ft_selectdata(cfg_t, dat_gamma);

        nTrl = numel(dat_stim_full.trial);
        if nTrl > 0
            cfg_cov = [];
            cfg_cov.covariance = 'yes';
            cfg_cov.covariancewindow = 'all';
            cfg_cov.removemean = 'yes';

            tl_base = ft_timelockanalysis(cfg_cov, dat_base);
            tl_stim_full = ft_timelockanalysis(cfg_cov, dat_stim_full);

            covBase_full = covBase_full + double(tl_base.cov) * nTrl;
            covStim_full = covStim_full + double(tl_stim_full.cov) * nTrl;
        end
        nTrials_total = nTrials_total + nTrl;
    end

    if nTrials_total < 1
        warning('No valid trials for subject %s after trial selection.', subject_id);
        continue;
    end
    covStim_full = covStim_full / nTrials_total;
    covBase_full = covBase_full / nTrials_total;

    %% GED (full window)
    covStim_reg = (1 - lambda) * covStim_full + lambda * mean(diag(covStim_full)) * eye(nChans);
    covBase_reg = (1 - lambda) * covBase_full + lambda * mean(diag(covBase_full)) * eye(nChans);

    [W_full, D_full] = eig(covStim_reg, covBase_reg);
    [evals_sorted, sortIdx] = sort(real(diag(D_full)), 'descend');
    W_full = W_full(:, sortIdx);

    searchFilters = nan(nChans, nSearch);
    searchTopos = nan(nChans, nSearch);
    searchMeanPrSpectrum = nan(nSearch, numel(scan_freqs));

    % Signed occipital template for stable topography orientation
    sim_template = build_signed_occipital_template(dataStructs{1}.label, headmodel, occ_idx, front_idx);

    for ci = 1:nSearch
        w_ci = W_full(:, ci);
        topo_ci = covStim_reg * w_ci;

        r_ci = corr(topo_ci, sim_template, 'rows', 'complete');
        if ~isnan(r_ci) && r_ci < 0
            w_ci = -w_ci;
            topo_ci = -topo_ci;
        end

        proxy_ci = estimate_component_artifact_proxies( ...
            w_ci, dat_per_cond, full_window, baseline_window, fsample, scan_freqs, scan_width);

        searchFilters(:, ci) = w_ci;
        searchTopos(:, ci) = topo_ci;
        searchMeanPrSpectrum(ci, :) = proxy_ci.mean_pr_spectrum(:)';
    end

    %% PF metrics (all mapped to 0-1)
    pf_current = compute_peak_form_template_score_from_spectra( ...
        searchMeanPrSpectrum, scan_freqs, analysis_freq_range, 5, 0.8, 0.35);

    [pf_kurtosis, kurtosis_raw] = compute_pf_kurtosis_from_spectra( ...
        searchMeanPrSpectrum, scan_freqs, analysis_freq_range);

    [pf_laplace, laplace_r2] = compute_pf_laplace_fit_from_spectra( ...
        searchMeanPrSpectrum, scan_freqs, analysis_freq_range);

    [pf_tgauss, tgauss_r2] = compute_pf_trunc_gauss_fit_from_spectra( ...
        searchMeanPrSpectrum, scan_freqs, analysis_freq_range);

    %% Figure 1: Scree + PF curves
    fig1 = figure('Position', [0 0 1512 982], 'Color', 'w');

    subplot(2, 2, 1);
    nEigPlot = min(30, numel(evals_sorted));
    plot(1:nEigPlot, evals_sorted(1:nEigPlot), 'k-o', 'LineWidth', 1.8, 'MarkerSize', 5);
    xlabel('Component rank');
    ylabel('\lambda');
    title('GED eigenvalue scree');
    box off;

    subplot(2, 2, 2); hold on;
    comp_idx = 1:nSearch;
    plot(comp_idx, pf_current(:)', '-o', 'LineWidth', 1.6, 'MarkerSize', 5, 'Color', [0.10 0.10 0.10], 'DisplayName', 'Current PF');
    plot(comp_idx, pf_kurtosis(:)', '-o', 'LineWidth', 1.6, 'MarkerSize', 5, 'Color', [0.16 0.45 0.76], 'DisplayName', 'Kurtosis PF');
    plot(comp_idx, pf_laplace(:)', '-o', 'LineWidth', 1.6, 'MarkerSize', 5, 'Color', [0.88 0.49 0.10], 'DisplayName', 'Laplace PF (R^2)');
    plot(comp_idx, pf_tgauss(:)', '-o', 'LineWidth', 1.6, 'MarkerSize', 5, 'Color', [0.22 0.66 0.33], 'DisplayName', 'TruncGaussian PF (R^2)');
    ylim([0 1]);
    xlim([1 nSearch]);
    xlabel('Component');
    ylabel('PF (0-1)');
    title('PF methods (primary scores)');
    legend('Location', 'southoutside', 'NumColumns', 2);
    box off;

    subplot(2, 2, 3); hold on;
    plot(comp_idx, laplace_r2(:)', '-o', 'LineWidth', 1.6, 'MarkerSize', 5, 'Color', [0.88 0.49 0.10], 'DisplayName', 'Laplace R^2');
    plot(comp_idx, tgauss_r2(:)', '-o', 'LineWidth', 1.6, 'MarkerSize', 5, 'Color', [0.22 0.66 0.33], 'DisplayName', 'TruncGaussian R^2');
    ylim([0 1]);
    xlim([1 nSearch]);
    xlabel('Component');
    ylabel('R^2 (0-1)');
    title('Fit R^2 diagnostics');
    legend('Location', 'best');
    box off;

    subplot(2, 2, 4);
    axis off;
    text(0.01, 0.95, sprintf('Subject: %s', subject_id), 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
    text(0.01, 0.84, sprintf('nTrials pooled: %d', nTrials_total), 'FontSize', 11, 'Interpreter', 'none');
    text(0.01, 0.74, sprintf('nComponents shown: %d', nSearch), 'FontSize', 11, 'Interpreter', 'none');
    text(0.01, 0.62, 'PF settings:', 'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none');
    text(0.03, 0.53, '(1) Current PF = full existing method', 'FontSize', 10, 'Interpreter', 'none');
    text(0.03, 0.45, '(2) Kurtosis PF = spectral-shape kurtosis mapping', 'FontSize', 10, 'Interpreter', 'none');
    text(0.03, 0.37, '(3) Laplace PF = R^2', 'FontSize', 10, 'Interpreter', 'none');
    text(0.03, 0.29, '(4) Truncated Gaussian PF = R^2', 'FontSize', 10, 'Interpreter', 'none');

    sgtitle(sprintf('GED PF comparison (full window): %s', subject_id), 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
    save_figure_png(fig1, fullfile(fig_save_dir, sprintf('GCP_eeg_GED_pf_scree_metrics_%s.png', subj_tag)));
    close(fig1);

    %% Figure 2: Topographies + spectra with all outcome values printed
    nShow = min(10, nSearch);
    nColsComp = 5;
    fig2 = figure('Position', [0 0 1512 982], 'Color', 'w');

    cfg_topo = [];
    cfg_topo.layout = headmodel.layANThead;
    cfg_topo.comment = 'no';
    cfg_topo.marker = 'off';
    cfg_topo.style = 'straight';
    cfg_topo.gridscale = 300;
    cfg_topo.colormap = '*RdBu';
    cfg_topo.figure = 'gcf';

    for k = 1:nShow
        ci = k;
        if ci <= nColsComp
            topo_pos = ci;
            spec_pos = nColsComp + ci;
        else
            topo_pos = (2 * nColsComp) + (ci - nColsComp);
            spec_pos = (3 * nColsComp) + (ci - nColsComp);
        end

        subplot(4, nColsComp, topo_pos);
        topo_data = [];
        topo_data.label = dataStructs{1}.label;
        topo_data.avg = searchTopos(:, ci);
        topo_data.dimord = 'chan';
        topo_vals = topo_data.avg(isfinite(topo_data.avg));
        topo_clim = max(abs(topo_vals));
        if ~isfinite(topo_clim) || topo_clim <= 0
            topo_clim = 1;
        end
        cfg_ci = cfg_topo;
        cfg_ci.zlim = [-topo_clim topo_clim];
        try
            ft_topoplotER(cfg_ci, topo_data);
        catch
            imagesc(topo_data.avg(:));
            axis tight;
            caxis([-topo_clim topo_clim]);
            colorbar;
        end
        title(sprintf('C%d, \\lambda=%.2f', ci, evals_sorted(ci)), 'FontSize', 9, 'Interpreter', 'none');

        subplot(4, nColsComp, spec_pos); hold on;
        spec = searchMeanPrSpectrum(ci, :);
        plot(scan_freqs, spec, 'k-', 'LineWidth', 1.5);
        yline(0, 'k--', 'LineWidth', 0.8);
        xlim(analysis_freq_range);
        spec_finite = spec(isfinite(spec));
        if ~isempty(spec_finite)
            ymin = min(spec_finite);
            ymax = max(spec_finite);
            if isfinite(ymin) && isfinite(ymax) && ymin < ymax
                yr = ymax - ymin;
                ylim([ymin - 0.12 * yr, ymax + 0.30 * yr]);
            end
        end
        xlabel('Hz', 'FontSize', 8);
        ylabel('Power [dB]', 'FontSize', 8);
        box on;
        set(gca, 'FontSize', 8);
        format_power_change_db_axis(gca);

        txt = sprintf(['Current=%.2f | Kurt=%.2f\n' ...
            'Laplace=%.2f (R^2=%.2f)\n' ...
            'TGauss=%.2f (R^2=%.2f)'], ...
            pf_current(ci), pf_kurtosis(ci), ...
            pf_laplace(ci), laplace_r2(ci), ...
            pf_tgauss(ci), tgauss_r2(ci));
        text(0.02, 0.98, txt, 'Units', 'normalized', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontSize', 7, 'Interpreter', 'none');
        title(sprintf('Component %d spectrum', ci), 'FontSize', 9, 'Interpreter', 'none');
    end

    sgtitle(sprintf('Top 10 components with PF values: %s', subject_id), 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
    save_figure_png(fig2, fullfile(fig_save_dir, sprintf('GCP_eeg_GED_pf_components_%s.png', subj_tag)));
    close(fig2);
end

%% ------------------------------------------------------------------------
%% PF methods

function [pf_score_vec, kurtosis_raw_vec] = compute_pf_kurtosis_from_spectra(mean_pr_spectrum, scan_freqs, analysis_freq_range)
nComp = size(mean_pr_spectrum, 1);
pf_score_vec = nan(nComp, 1);
kurtosis_raw_vec = nan(nComp, 1);
if isempty(mean_pr_spectrum) || isempty(scan_freqs)
    return;
end
freq_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
for ci = 1:nComp
    y = mean_pr_spectrum(ci, :);
    yb = y(freq_mask);
    valid = isfinite(yb);
    yb = yb(valid);
    if numel(yb) < 7
        continue;
    end
    yb = movmean(yb, 5, 'omitnan');
    yshape = max(yb - median(yb, 'omitnan'), 0);
    if all(~isfinite(yshape)) || max(yshape) <= eps
        pf_score_vec(ci) = 0;
        kurtosis_raw_vec(ci) = NaN;
        continue;
    end
    k = kurtosis(yshape);
    kurtosis_raw_vec(ci) = k;
    % Anchored absolute mapping to [0,1]: excess kurtosis normalized by a fixed cap.
    excess = max(0, k - 3);
    pf_score_vec(ci) = max(0, min(1, excess / 12));
end
pf_score_vec(~isfinite(pf_score_vec)) = 0;
end

function [pf_score_vec, r2_vec] = compute_pf_laplace_fit_from_spectra(mean_pr_spectrum, scan_freqs, analysis_freq_range)
nComp = size(mean_pr_spectrum, 1);
pf_score_vec = nan(nComp, 1);
r2_vec = nan(nComp, 1);
freq_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
for ci = 1:nComp
    y = mean_pr_spectrum(ci, :);
    x = scan_freqs(freq_mask);
    y = y(freq_mask);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    if numel(y) < 7
        continue;
    end
    [x, sidx] = sort(x(:)');
    y = y(sidx);
    y = movmean(y, 5, 'omitnan');
    yshape = max(y - prctile(y, 20), 0);
    yscale = max(yshape);
    if ~isfinite(yscale) || yscale <= eps
        r2_vec(ci) = 0;
        pf_score_vec(ci) = 0;
        continue;
    end
    yfit_target = yshape / yscale;
    [~, idx_pk] = max(yfit_target);
    p0 = [1.0, x(idx_pk), 6.0];
    lb = [0.20, analysis_freq_range(1), 1.0];
    ub = [1.50, analysis_freq_range(2), 20.0];

    obj = @(p) mean((yfit_target - laplace_model(x, p)).^2, 'omitnan');
    options = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 3000);
    p_est = fminsearch(@(pp) bounded_obj(pp, lb, ub, obj), p0, options);
    p_est = min(max(p_est, lb), ub);
    yhat = laplace_model(x, p_est);

    r2 = compute_r2_score(yfit_target, yhat);
    r2_vec(ci) = r2;
    pf_score_vec(ci) = r2; % user-selected primary metric
end
pf_score_vec(~isfinite(pf_score_vec)) = 0;
r2_vec(~isfinite(r2_vec)) = 0;
end

function [pf_score_vec, r2_vec] = compute_pf_trunc_gauss_fit_from_spectra(mean_pr_spectrum, scan_freqs, analysis_freq_range)
nComp = size(mean_pr_spectrum, 1);
pf_score_vec = nan(nComp, 1);
r2_vec = nan(nComp, 1);
freq_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
for ci = 1:nComp
    y = mean_pr_spectrum(ci, :);
    x = scan_freqs(freq_mask);
    y = y(freq_mask);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    if numel(y) < 7
        continue;
    end
    [x, sidx] = sort(x(:)');
    y = y(sidx);
    y = movmean(y, 5, 'omitnan');
    yshape = max(y - prctile(y, 20), 0);
    yscale = max(yshape);
    if ~isfinite(yscale) || yscale <= eps
        r2_vec(ci) = 0;
        pf_score_vec(ci) = 0;
        continue;
    end
    yfit_target = yshape / yscale;
    [~, idx_pk] = max(yfit_target);
    p0 = [1.0, x(idx_pk), 6.0];
    lb = [0.20, analysis_freq_range(1), 1.0];
    ub = [1.50, analysis_freq_range(2), 20.0];

    obj = @(p) mean((yfit_target - trunc_gauss_model(x, p, analysis_freq_range)).^2, 'omitnan');
    options = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 3000);
    p_est = fminsearch(@(pp) bounded_obj(pp, lb, ub, obj), p0, options);
    p_est = min(max(p_est, lb), ub);
    yhat = trunc_gauss_model(x, p_est, analysis_freq_range);

    r2 = compute_r2_score(yfit_target, yhat);
    r2_vec(ci) = r2;
    pf_score_vec(ci) = r2; % user-selected primary metric
end
pf_score_vec(~isfinite(pf_score_vec)) = 0;
r2_vec(~isfinite(r2_vec)) = 0;
end

function v = bounded_obj(p, lb, ub, obj_fun)
p_clip = min(max(p, lb), ub);
v = obj_fun(p_clip);
if ~isfinite(v)
    v = 1e6;
end
end

function yhat = laplace_model(x, p)
A = p(1);
mu = p(2);
b = p(3);
yhat = A .* exp(-abs(x - mu) ./ max(b, eps));
yhat = max(yhat, 0);
if max(yhat) > 0
    yhat = yhat ./ max(yhat);
end
end

function yhat = trunc_gauss_model(x, p, frange)
A = p(1);
mu = p(2);
sigma = p(3);
yhat = A .* exp(-0.5 * ((x - mu) ./ max(sigma, eps)).^2);
yhat(x < frange(1) | x > frange(2)) = 0;
yhat = max(yhat, 0);
if max(yhat) > 0
    yhat = yhat ./ max(yhat);
end
end

function r2 = compute_r2_score(y, yhat)
if numel(y) ~= numel(yhat)
    r2 = 0;
    return;
end
valid = isfinite(y) & isfinite(yhat);
y = y(valid);
yhat = yhat(valid);
if numel(y) < 3
    r2 = 0;
    return;
end
sse = sum((y - yhat).^2);
sst = sum((y - mean(y)).^2);
if ~isfinite(sst) || sst <= eps
    r2 = 0;
else
    r2 = 1 - (sse / sst);
end
r2 = max(0, min(1, r2));
end

%% ------------------------------------------------------------------------
%% Current PF function (kept equivalent to current implementation)

function [peak_form_score_vec, peak_form_mode_vec, diag] = compute_peak_form_template_score_from_spectra( ...
    mean_pr_spectrum, scan_freqs, analysis_freq_range, multi_peak_sep_min_hz, multi_peak_height_ratio_min, multi_peak_penalty_mult)
if nargin < 4 || isempty(multi_peak_sep_min_hz) || ~isfinite(multi_peak_sep_min_hz)
    multi_peak_sep_min_hz = 15;
end
if nargin < 5 || isempty(multi_peak_height_ratio_min) || ~isfinite(multi_peak_height_ratio_min)
    multi_peak_height_ratio_min = 0.80;
end
if nargin < 6 || isempty(multi_peak_penalty_mult) || ~isfinite(multi_peak_penalty_mult)
    multi_peak_penalty_mult = 0.62;
end
multi_peak_penalty_mult = max(0, min(1, multi_peak_penalty_mult));
shift_max_hz = 10;
single_widths_hz = [5 7 9 12];
double_separations_hz = [8 12 16 20];
double_widths_hz = [3 4 5 6];
min_trough_depth = 0.10;
min_similarity = 0.40;
smooth_n = 5;
peak_width_min_hz = 2.0;
peak_width_max_hz = 20.0;
edge_ratio_soft = 1.25;
edge_ratio_limit = 1.75;
edge_run_soft = 0.025;
edge_run_limit = 0.060;
peak_form_prom_abs_floor = 0.02;
nComp = size(mean_pr_spectrum, 1);
peak_form_score_vec = zeros(nComp, 1);
peak_form_mode_vec = repmat({'none'}, nComp, 1);
diag = struct( ...
    'best_single_similarity', nan(nComp, 1), ...
    'best_double_similarity', nan(nComp, 1), ...
    'best_similarity_raw', nan(nComp, 1), ...
    'pre_penalty_similarity', nan(nComp, 1), ...
    'total_penalty_raw', nan(nComp, 1), ...
    'total_penalty_used', nan(nComp, 1), ...
    'edge_penalty', nan(nComp, 1), ...
    'dominance_penalty', nan(nComp, 1), ...
    'dominance_rival_count_major', nan(nComp, 1), ...
    'dominance_rival_count_minor', nan(nComp, 1), ...
    'dominance_rival_amp_ratio_sum', nan(nComp, 1), ...
    'hf_rise_penalty', nan(nComp, 1), ...
    'best_shift_hz', nan(nComp, 1), ...
    'best_center_hz', nan(nComp, 1), ...
    'best_width_hz', nan(nComp, 1), ...
    'best_separation_hz', nan(nComp, 1), ...
    'best_trough_depth', nan(nComp, 1), ...
    'edge_ratio', nan(nComp, 1), ...
    'edge_run_score', nan(nComp, 1), ...
    'edge_artifact_flag', false(nComp, 1), ...
    'roughness_ratio', nan(nComp, 1), ...
    'roughness_penalty', nan(nComp, 1), ...
    'multi_peak_penalty', nan(nComp, 1), ...
    'multi_peak_flag', false(nComp, 1), ...
    'dominant_penalty_tag', {repmat({'none'}, nComp, 1)});
if isempty(mean_pr_spectrum) || isempty(scan_freqs)
    return;
end
freq_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
if ~any(freq_mask)
    return;
end
smooth_n = max(1, round(smooth_n));
for ci = 1:nComp
    y = mean_pr_spectrum(ci, :);
    if all(~isfinite(y))
        continue;
    end
    x_band = scan_freqs(freq_mask);
    y_proc = y(freq_mask);
    y_proc = movmean(y_proc, smooth_n, 'omitnan');
    y_band = y_proc;
    valid = isfinite(x_band) & isfinite(y_band);
    x_band = x_band(valid);
    y_band = y_band(valid);
    y_resid = y_band;
    if numel(y_band) < 7
        continue;
    end
    y_shape = max(y_band - median(y_band(isfinite(y_band))), 0);
    if ~any(isfinite(y_shape))
        y_shape = zeros(size(y_band));
    end
    [single_best, single_meta] = evaluate_single_template_bank(y_shape, x_band, single_widths_hz, shift_max_hz);
    [double_best, double_meta] = evaluate_double_template_bank( ...
        y_shape, y_band, x_band, double_widths_hz, double_separations_hz, shift_max_hz, min_trough_depth);
    diag.best_single_similarity(ci) = single_best;
    diag.best_double_similarity(ci) = double_best;

    y_floor = prctile(y_band(isfinite(y_band)), 20);
    if ~isfinite(y_floor)
        y_floor = median(y_band(isfinite(y_band)));
    end
    if ~isfinite(y_floor)
        y_floor = 0;
    end
    y_pos = max(y_band - y_floor, 0);
    robust_scale = robust_mad(y_pos);
    if ~isfinite(robust_scale) || robust_scale <= eps
        robust_scale = iqr(y_pos);
    end
    if ~isfinite(robust_scale) || robust_scale <= eps
        robust_scale = std(y_pos(isfinite(y_pos)));
    end
    if ~isfinite(robust_scale) || robust_scale <= eps
        robust_scale = 1;
    end
    rel_prom = 0.08 * max(y_pos);
    min_prom = max([0, rel_prom, peak_form_prom_abs_floor, 0.15 * robust_scale]);
    [pks, locs, widths] = findpeaks(y_pos, x_band, 'MinPeakProminence', min_prom, 'MinPeakDistance', 5);
    if isempty(pks)
        continue;
    end
    pks = pks(:);
    locs = locs(:);
    widths = widths(:);
    [~, dom_idx] = max(pks);
    dom_amp = pks(dom_idx);
    dom_width = widths(dom_idx);
    dom_loc = locs(dom_idx);

    multi_peak_pen = 1;
    multi_peak_hit = false;
    if numel(pks) >= 2 && isfinite(dom_amp) && dom_amp > eps
        sep_ok = abs(locs - dom_loc) >= multi_peak_sep_min_hz;
        height_ok = pks >= multi_peak_height_ratio_min * dom_amp;
        rival_mask = ((1:numel(pks))' ~= dom_idx) & sep_ok & height_ok;
        if any(rival_mask)
            multi_peak_hit = true;
            n_rivals = sum(rival_mask);
            rival_amp_ratio = pks(rival_mask) ./ max(dom_amp, eps);
            amp_excess = max(0, mean(rival_amp_ratio) - multi_peak_height_ratio_min);
            amp_scale = amp_excess / max(1 - multi_peak_height_ratio_min, eps);
            extra_loss = 0.18 * max(0, n_rivals - 1) + 0.15 * min(1, amp_scale);
            multi_peak_pen = max(0.15, multi_peak_penalty_mult - extra_loss);
        end
    end

    dominance_pen = 1;
    n_rival_major = 0;
    n_rival_minor = 0;
    rival_amp_ratio_sum = 0;
    if numel(pks) >= 2 && isfinite(dom_amp) && dom_amp > eps
        sep_from_dom = abs(locs - dom_loc);
        is_rival = ((1:numel(pks))' ~= dom_idx);
        rival_ratios = pks ./ max(dom_amp, eps);

        major_ratio_min = 0.55;
        minor_ratio_min = 0.30;
        major_sep_min_hz = 4.0;
        minor_sep_min_hz = 2.5;

        major_mask = is_rival & (sep_from_dom >= major_sep_min_hz) & (rival_ratios >= major_ratio_min);
        minor_mask = is_rival & (sep_from_dom >= minor_sep_min_hz) & (rival_ratios >= minor_ratio_min);

        n_rival_major = sum(major_mask);
        n_rival_minor = sum(minor_mask);
        if any(major_mask)
            rival_amp_ratio_sum = sum(rival_ratios(major_mask));
        end

        clutter_load = 0;
        clutter_load = clutter_load + 0.35 * max(0, n_rival_major - 1);
        clutter_load = clutter_load + 0.12 * max(0, n_rival_minor - 2);
        clutter_load = clutter_load + 0.45 * max(0, rival_amp_ratio_sum - 0.85);

        max_loss = 0.55;
        dominance_pen = max(0.45, 1 - min(max_loss, clutter_load));
    end

    width_score = 1;
    if isfinite(dom_width)
        if dom_width < peak_width_min_hz
            width_score = max(0.15, dom_width / max(peak_width_min_hz, eps));
        elseif dom_width > peak_width_max_hz
            width_score = max(0.20, peak_width_max_hz / max(dom_width, eps));
        end
    end
    snr_score = dom_amp / (dom_amp + robust_scale);
    amp_score = dom_amp / (dom_amp + 0.70 * robust_scale);
    shape_score = max(single_best, double_best);
    peak_core_hz = 6;
    peak_core_mask = abs(x_band - dom_loc) <= peak_core_hz;
    conc_ratio = sum(y_pos(peak_core_mask)) / max(sum(y_pos), eps);
    concentration_pen = 1;
    if isfinite(conc_ratio) && conc_ratio < 0.32
        concentration_pen = max(0.35, conc_ratio / 0.32);
    end
    dominant_quality = max(0, min(1, (0.40 * amp_score + 0.30 * snr_score + 0.30 * shape_score) * width_score * concentration_pen));

    best_pre_penalty = dominant_quality;
    dual_candidate = isfinite(double_best) && isfinite(single_best) && ...
        (double_best > (single_best + 0.05)) && isfinite(double_meta.trough_depth) && ...
        (double_meta.trough_depth >= max(0.05, 0.75 * min_trough_depth));
    if dual_candidate
        mode_raw = 'dual';
        best_shift = double_meta.shift_hz;
        best_width = double_meta.width_hz;
        best_sep = double_meta.sep_hz;
        best_trough = double_meta.trough_depth;
        if ~isempty(double_meta.centers_hz)
            best_centers = double_meta.centers_hz;
        else
            best_centers = dom_loc;
        end
    else
        mode_raw = 'single';
        best_shift = single_meta.shift_hz;
        best_width = dom_width;
        best_sep = double_meta.sep_hz;
        best_trough = double_meta.trough_depth;
        best_centers = dom_loc;
    end

    hf_pen = 1;
    hf_mask = x_band >= max(70, analysis_freq_range(2) - 15);
    if sum(hf_mask) >= 5
        hf_idx = find(hf_mask);
        hf_rho = corr((1:numel(hf_idx))', y_band(hf_idx)', 'rows', 'complete', 'type', 'Spearman');
        if isfinite(hf_rho) && hf_rho > 0.70
            hf_pen = max(0.65, 1 - 0.30 * (hf_rho - 0.70) / 0.30);
        end
    end

    [edge_ratio, edge_run_score, edge_flag] = compute_edge_artifact_indicators( ...
        y_resid, x_band, analysis_freq_range, edge_ratio_limit, edge_run_limit);
    diag.edge_ratio(ci) = edge_ratio;
    diag.edge_run_score(ci) = edge_run_score;
    diag.edge_artifact_flag(ci) = edge_flag;
    edge_pen = 1;
    if isfinite(edge_ratio) && edge_ratio > edge_ratio_soft
        edge_pen = edge_pen * max(0.80, 1 - 0.20 * min(1, (edge_ratio - edge_ratio_soft) / max(edge_ratio_limit - edge_ratio_soft, eps)));
    end
    if isfinite(edge_run_score) && edge_run_score > edge_run_soft
        edge_pen = edge_pen * max(0.80, 1 - 0.20 * min(1, (edge_run_score - edge_run_soft) / max(edge_run_limit - edge_run_soft, eps)));
    end
    if edge_flag
        edge_pen = edge_pen * 0.90;
    end

    roughness_ratio = NaN;
    roughness_pen = 1;
    y_finite = y_band(isfinite(y_band));
    if numel(y_finite) >= 3
        y_range = max(y_finite) - min(y_finite);
        if y_range > eps
            y_scaled = y_band / y_range;
        else
            y_scaled = y_band;
        end
        dy = diff(y_scaled);
        amp_scale = prctile(y_scaled, 75) - prctile(y_scaled, 25);
        if ~isfinite(amp_scale) || amp_scale <= eps
            amp_scale = std(y_scaled(isfinite(y_scaled)));
        end
        if ~isfinite(amp_scale) || amp_scale <= eps
            amp_scale = 1;
        end
        roughness_ratio = robust_mad(dy) / amp_scale;
        rough_ref = 0.40;
        rough_span = 0.70;
        rough_pen_floor = 0.55;
        rough_pen_max_loss = 0.45;
        if isfinite(roughness_ratio) && roughness_ratio > rough_ref
            loss_frac = min(1, (roughness_ratio - rough_ref) / rough_span);
            roughness_pen = max(rough_pen_floor, 1 - rough_pen_max_loss * loss_frac);
        end
    end

    penalty_raw = edge_pen * dominance_pen * hf_pen * roughness_pen * multi_peak_pen;
    penalty_floor = 0.20;
    penalty_used = max(penalty_floor, penalty_raw);
    best_raw = best_pre_penalty * penalty_used;

    penalty_names = {'edge', 'dominance', 'hf_rise', 'roughness', 'multi_peak'};
    penalty_vals = [edge_pen, dominance_pen, hf_pen, roughness_pen, multi_peak_pen];
    [worst_pen, worst_idx] = min(penalty_vals);
    if isfinite(worst_pen) && (worst_pen < 0.999)
        dominant_penalty_tag = penalty_names{worst_idx};
    else
        dominant_penalty_tag = 'none';
    end

    diag.best_similarity_raw(ci) = best_raw;
    diag.pre_penalty_similarity(ci) = best_pre_penalty;
    diag.total_penalty_raw(ci) = penalty_raw;
    diag.total_penalty_used(ci) = penalty_used;
    diag.edge_penalty(ci) = edge_pen;
    diag.dominance_penalty(ci) = dominance_pen;
    diag.dominance_rival_count_major(ci) = n_rival_major;
    diag.dominance_rival_count_minor(ci) = n_rival_minor;
    diag.dominance_rival_amp_ratio_sum(ci) = rival_amp_ratio_sum;
    diag.hf_rise_penalty(ci) = hf_pen;
    diag.best_shift_hz(ci) = best_shift;
    if ~isempty(best_centers) && all(isfinite(best_centers))
        diag.best_center_hz(ci) = mean(best_centers);
    end
    diag.best_width_hz(ci) = best_width;
    diag.best_separation_hz(ci) = best_sep;
    diag.best_trough_depth(ci) = best_trough;
    diag.roughness_ratio(ci) = roughness_ratio;
    diag.roughness_penalty(ci) = roughness_pen;
    diag.multi_peak_penalty(ci) = multi_peak_pen;
    diag.multi_peak_flag(ci) = multi_peak_hit;
    diag.dominant_penalty_tag{ci} = dominant_penalty_tag;

    if ~isfinite(best_raw)
        peak_form_score_vec(ci) = 0;
        peak_form_mode_vec{ci} = 'none';
    else
        score_mapped = best_raw;
        if isfinite(min_similarity) && (min_similarity > 0) && (best_raw < min_similarity)
            score_mapped = 0.75 * score_mapped;
        end
        peak_form_score_vec(ci) = max(0, min(1, score_mapped));
        peak_form_mode_vec{ci} = mode_raw;
    end
end
end

function [best_sim, meta] = evaluate_single_template_bank(y_shape, x_band, widths_hz, shift_max_hz)
best_sim = 0;
meta = struct('shift_hz', NaN, 'width_hz', NaN, 'centers_hz', []);
if isempty(widths_hz)
    return;
end
x_mid = mean(x_band);
x_min = min(x_band);
x_max = max(x_band);
for wi = 1:numel(widths_hz)
    w = widths_hz(wi);
    if ~isfinite(w) || w <= 0
        continue;
    end
    shift_vals = candidate_shift_values_hz(x_band, shift_max_hz);
    for si = 1:numel(shift_vals)
        shift_hz = shift_vals(si);
        center_hz = x_mid + shift_hz;
        if ~isfinite(center_hz) || center_hz < x_min || center_hz > x_max
            continue;
        end
        tpl = gaussian_template(x_band, center_hz, w);
        sim = safe_template_similarity(y_shape, tpl);
        if sim > best_sim
            best_sim = sim;
            meta.shift_hz = shift_hz;
            meta.width_hz = w;
            meta.centers_hz = center_hz;
        end
    end
end
end

function [best_sim, meta] = evaluate_double_template_bank(y_shape, y_band, x_band, widths_hz, separations_hz, shift_max_hz, min_trough_depth)
best_sim = 0;
meta = struct('shift_hz', NaN, 'width_hz', NaN, 'sep_hz', NaN, 'trough_depth', NaN, 'centers_hz', []);
if isempty(widths_hz) || isempty(separations_hz)
    return;
end
x_mid = mean(x_band);
x_min = min(x_band);
x_max = max(x_band);
shift_vals = candidate_shift_values_hz(x_band, shift_max_hz);
for wi = 1:numel(widths_hz)
    w = widths_hz(wi);
    if ~isfinite(w) || w <= 0
        continue;
    end
    for di = 1:numel(separations_hz)
        sep = separations_hz(di);
        if ~isfinite(sep) || sep <= 0
            continue;
        end
        for si = 1:numel(shift_vals)
            shift_hz = shift_vals(si);
            c1 = x_mid + shift_hz - sep / 2;
            c2 = x_mid + shift_hz + sep / 2;
            if ~isfinite(c1) || ~isfinite(c2) || c1 < x_min || c2 > x_max
                continue;
            end
            tpl = gaussian_template(x_band, c1, w) + gaussian_template(x_band, c2, w);
            sim = safe_template_similarity(y_shape, tpl);
            trough_depth = estimate_trough_depth(y_band, x_band, c1, c2);
            trough_scale = min(1, max(0, trough_depth) / max(min_trough_depth, eps));
            sim_adj = sim * trough_scale;
            if sim_adj > best_sim
                best_sim = sim_adj;
                meta.shift_hz = shift_hz;
                meta.width_hz = w;
                meta.sep_hz = sep;
                meta.trough_depth = trough_depth;
                meta.centers_hz = [c1 c2];
            end
        end
    end
end
end

function shifts = candidate_shift_values_hz(x_band, shift_max_hz)
if numel(x_band) >= 2
    df = median(diff(x_band));
else
    df = 1;
end
if ~isfinite(df) || df <= 0
    df = 1;
end
shift_max_hz = max(0, shift_max_hz);
max_possible = max(abs(x_band - mean(x_band)));
if ~isfinite(max_possible) || max_possible <= 0
    max_possible = shift_max_hz;
end
shift_limit = min(shift_max_hz, max_possible);
shifts = -shift_limit:df:shift_limit;
if isempty(shifts)
    shifts = 0;
end
if ~any(abs(shifts) < 1e-9)
    shifts = unique(sort([shifts, 0]));
end
end

function depth = estimate_trough_depth(y_band, x_band, c1, c2)
depth = 0;
if ~(isfinite(c1) && isfinite(c2))
    return;
end
if c1 > c2
    tmp = c1;
    c1 = c2;
    c2 = tmp;
end
[~, i1] = min(abs(x_band - c1));
[~, i2] = min(abs(x_band - c2));
if i1 == i2
    return;
end
idx_lo = min(i1, i2);
idx_hi = max(i1, i2);
y_pos = max(y_band(:), 0);
p1 = y_pos(i1);
p2 = y_pos(i2);
if idx_hi - idx_lo < 2 || ~isfinite(p1) || ~isfinite(p2)
    return;
end
valley = min(y_pos(idx_lo:idx_hi));
peak_ref = max(min(p1, p2), eps);
depth = max(0, min(1, 1 - valley / peak_ref));
end

function sim = safe_template_similarity(y_shape, tpl)
sim = 0;
if isempty(y_shape) || isempty(tpl)
    return;
end
ys = y_shape(:);
ys(~isfinite(ys)) = 0;
ys = max(ys, 0);
tpl = tpl(:);
tpl = max(tpl, 0);
if numel(tpl) ~= numel(ys)
    return;
end
mass = sum(ys);
if ~isfinite(mass) || mass <= eps
    return;
end
sim = sum(ys .* tpl) / mass;
if ~isfinite(sim)
    sim = 0;
end
sim = max(0, min(1, sim));
end

function g = gaussian_template(x, mu, sigma)
if ~isfinite(mu) || ~isfinite(sigma) || sigma <= 0
    g = zeros(size(x));
    return;
end
g = exp(-0.5 * ((x - mu) ./ sigma).^2);
end

function [edge_ratio, edge_run_score, edge_artifact_flag] = compute_edge_artifact_indicators( ...
    y_resid, x_band, analysis_freq_range, edge_ratio_limit, edge_run_limit)
edge_ratio = NaN;
edge_run_score = NaN;
edge_artifact_flag = false;
if isempty(y_resid) || isempty(x_band) || numel(y_resid) ~= numel(x_band)
    return;
end
x_lo = analysis_freq_range(1);
x_hi = analysis_freq_range(2);
edge_mask = (x_band <= (x_lo + 5)) | (x_band >= (x_hi - 5));
interior_mask = x_band >= (x_lo + 10) & x_band <= (x_hi - 10);
if sum(edge_mask) < 3 || sum(interior_mask) < 3
    return;
end
edge_amp = prctile(abs(y_resid(edge_mask)), 90);
interior_amp = prctile(abs(y_resid(interior_mask)), 90);
edge_ratio = edge_amp / max(interior_amp, eps);
left_mask = x_band <= (x_lo + 6);
right_mask = x_band >= (x_hi - 6);
left_run = 0;
right_run = 0;
if sum(left_mask) >= 4
    left_run = abs(mean(diff(y_resid(left_mask))));
end
if sum(right_mask) >= 4
    right_run = abs(mean(diff(y_resid(right_mask))));
end
edge_run_score = max(left_run, right_run);
edge_artifact_flag = (isfinite(edge_ratio) && edge_ratio > edge_ratio_limit) || ...
    (isfinite(edge_run_score) && edge_run_score > edge_run_limit);
end

%% ------------------------------------------------------------------------
%% GED component spectrum proxy functions

function proxy = estimate_component_artifact_proxies(filter_w, dat_per_cond, stim_window, base_window, fs, scan_freqs, scan_width)
proxy = struct( ...
    'lineharm_ratio', NaN, ...
    'stationarity_cv', NaN, ...
    'burst_ratio', NaN, ...
    'hf_slope', NaN, ...
    'mean_pr_spectrum', nan(1, numel(scan_freqs)), ...
    'n_trials_used', 0);
if isempty(filter_w) || isempty(dat_per_cond)
    return;
end

band_mask = scan_freqs >= 30 & scan_freqs <= 90;
hf_mask = scan_freqs >= 70 & scan_freqs <= min(110, max(scan_freqs));
harm_mask = ismember(round(scan_freqs), [50 60 100]);
if ~any(harm_mask)
    [~, idx50] = min(abs(scan_freqs - 50));
    harm_mask(idx50) = true;
end

trial_gamma = [];
lineharm_acc = [];
trial_count = 0;
stim_sig_cells = {};
base_sig_cells = {};
for cond = 1:numel(dat_per_cond)
    dat = dat_per_cond{cond};
    if isempty(dat) || ~isfield(dat, 'trial')
        continue;
    end
    for trl = 1:numel(dat.trial)
        x = double(dat.trial{trl});
        t = dat.time{trl};
        if isempty(x) || isempty(t)
            continue;
        end
        x_proj = filter_w' * x;
        idx_base = t >= base_window(1) & t <= base_window(2);
        idx_stim = t >= stim_window(1) & t <= stim_window(2);
        if sum(idx_base) < 5 || sum(idx_stim) < 5
            continue;
        end
        base_sig_cells{end+1, 1} = x_proj(idx_base);
        stim_sig_cells{end+1, 1} = x_proj(idx_stim);
        trial_count = trial_count + 1;
    end
end
proxy.n_trials_used = trial_count;
if trial_count < 4
    return;
end
if ~isfinite(fs) || fs <= 0
    return;
end
trial_pr = compute_simple_power_ratio_scan_batch(stim_sig_cells, base_sig_cells, fs, scan_freqs, scan_width);
if isempty(trial_pr)
    return;
end
for ri = 1:size(trial_pr, 1)
    pr_full = trial_pr(ri, :);
    if all(~isfinite(pr_full)) || ~any(band_mask)
        continue;
    end
    pr_band = pr_full(band_mask);
    if all(~isfinite(pr_band))
        continue;
    end
    nonharm_val = nanmean(pr_full(band_mask & ~harm_mask));
    harm_val = nanmean(pr_full(band_mask & harm_mask));
    if ~isfinite(nonharm_val) || abs(nonharm_val) <= eps || ~isfinite(harm_val)
        continue;
    end
    trial_gamma(end+1, 1) = nanmean(pr_band); %#ok<AGROW>
    lineharm_acc(end+1, 1) = max(harm_val, 0) / max(abs(nonharm_val), eps); %#ok<AGROW>
end
if isempty(trial_gamma) || isempty(lineharm_acc)
    return;
end
proxy.mean_pr_spectrum = nanmean(trial_pr, 1);
proxy.lineharm_ratio = nanmean(lineharm_acc);
proxy.stationarity_cv = nanstd(trial_gamma) / max(abs(nanmean(trial_gamma)), eps);
burst_thr = prctile(abs(trial_gamma), 85);
if ~isfinite(burst_thr) || burst_thr <= 0
    burst_thr = median(abs(trial_gamma)) + mad(abs(trial_gamma), 1);
end
proxy.burst_ratio = mean(abs(trial_gamma) >= max(burst_thr, eps));

if sum(hf_mask) >= 3
    f_hf = scan_freqs(hf_mask);
    spec_hf = proxy.mean_pr_spectrum(hf_mask);
    xh_full = log(f_hf(:));
    yh_full = log(abs(spec_hf(:)) + eps);
    valid_hf = isfinite(xh_full) & isfinite(yh_full) & isfinite(spec_hf(:));
    if nnz(valid_hf) >= 3
        p = polyfit(xh_full(valid_hf), yh_full(valid_hf), 1);
        if isfinite(p(1))
            proxy.hf_slope = p(1);
        end
    end
end
end

function pr_mat = compute_simple_power_ratio_scan_batch(sig_stim_cells, sig_base_cells, fs, scan_freqs, scan_width)
nRows = numel(sig_stim_cells);
pr_mat = nan(nRows, numel(scan_freqs));
if nRows == 0 || numel(sig_base_cells) ~= nRows || fs <= 0
    return;
end
[p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(sig_stim_cells, sig_base_cells, fs, scan_freqs, scan_width);
if isempty(p_stim_scan) || isempty(p_base_scan)
    return;
end
for ri = 1:nRows
    p_stim_row = p_stim_scan(ri, :);
    p_base_row = p_base_scan(ri, :);
    valid_base = isfinite(p_base_row) & (p_base_row > 0);
    if ~any(valid_base)
        continue;
    end
    base_anchor = prctile(p_base_row(valid_base), 20);
    base_median = median(p_base_row(valid_base), 'omitnan');
    if ~isfinite(base_anchor) || base_anchor <= 0
        base_anchor = base_median;
    end
    if ~isfinite(base_median) || base_median <= 0
        base_median = base_anchor;
    end
    base_floor = max(0.25 * base_anchor, eps);
    valid = isfinite(p_stim_row) & isfinite(p_base_row) & (p_base_row > 0);
    ratio_row = nan(1, numel(scan_freqs));
    ratio_row(valid) = 10 * log10((p_stim_row(valid) + base_floor) ./ (p_base_row(valid) + base_floor));
    pr_mat(ri, :) = ratio_row;
end
end

function [p_stim_scan, p_base_scan] = compute_scan_power_mtmfft_ft_pair(sig_stim, sig_base, fs, scan_freqs, scan_width_hz)
if ~iscell(sig_stim) && isvector(sig_stim)
    sig_stim = sig_stim(:)';
end
if ~iscell(sig_base) && isvector(sig_base)
    sig_base = sig_base(:)';
end
nFreq = numel(scan_freqs);
if iscell(sig_stim)
    nStim = numel(sig_stim);
else
    nStim = size(sig_stim, 1);
end
p_stim_scan = nan(nStim, nFreq);
p_base_scan = nan(nStim, nFreq);
if isempty(sig_stim) || isempty(sig_base) || fs <= 0
    return;
end
if iscell(sig_stim) ~= iscell(sig_base)
    return;
end
if iscell(sig_stim)
    if numel(sig_stim) ~= numel(sig_base)
        return;
    end
    sig_all = [sig_stim(:); sig_base(:)];
else
    if size(sig_stim, 1) ~= size(sig_base, 1)
        return;
    end
    sig_all = [sig_stim; sig_base];
end
if isempty(sig_all)
    return;
end
p_all = compute_scan_power_mtmfft_ft(sig_all, fs, scan_freqs, scan_width_hz);
if isempty(p_all) || size(p_all, 1) < 2 * nStim
    return;
end
p_stim_scan = p_all(1:nStim, :);
p_base_scan = p_all((nStim + 1):(2 * nStim), :);
end

function p_scan = compute_scan_power_mtmfft_ft(sig, fs, scan_freqs, scan_width_hz)
if ~iscell(sig) && isvector(sig)
    sig = sig(:)';
end
if iscell(sig)
    nSig = numel(sig);
else
    nSig = size(sig, 1);
end
p_scan = nan(nSig, numel(scan_freqs));
if isempty(sig) || fs <= 0 || isempty(scan_freqs) || ~isfinite(scan_width_hz) || scan_width_hz <= 0
    return;
end
if ~iscell(sig) && size(sig, 2) < 8
    return;
end

valid_rows = false(nSig, 1);
dat = [];
dat.label = {'comp'};
dat.fsample = fs;
dat.trial = {};
dat.time = {};
for si = 1:nSig
    if iscell(sig)
        x = double(sig{si});
    else
        x = double(sig(si, :));
    end
    if iscolumn(x)
        x = x';
    end
    if any(~isfinite(x)) || numel(x) < 8
        continue;
    end
    x = x - mean(x);
    valid_rows(si) = true;
    dat.trial{end+1} = x;
    dat.time{end+1} = (0:(numel(x)-1)) / fs;
end
if ~any(valid_rows)
    return;
end

nTrials_valid = numel(dat.trial);
dat.sampleinfo = zeros(nTrials_valid, 2);
sample_start = 1;
for ti = 1:nTrials_valid
    nSamp = numel(dat.trial{ti});
    dat.sampleinfo(ti, :) = [sample_start, sample_start + nSamp - 1];
    sample_start = sample_start + nSamp;
end

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'dpss';
cfg.foi = scan_freqs;
cfg.tapsmofrq = scan_width_hz;
cfg.pad = 'nextpow2';
cfg.keeptrials = 'yes';
cfg.feedback = 'none';
try
    freq = ft_freqanalysis(cfg, dat);
catch
    return;
end

pow = freq.powspctrm;
if ndims(pow) == 3
    pow = squeeze(pow(:, 1, :));
elseif isvector(pow)
    pow = pow(:)';
end
if isfield(freq, 'freq') && ~isempty(freq.freq)
    n_freq = numel(freq.freq);
else
    n_freq = numel(scan_freqs);
end
if isvector(pow) && numel(pow) == n_freq
    pow = reshape(pow, 1, []);
elseif size(pow, 1) == n_freq && size(pow, 2) ~= n_freq
    pow = pow.';
end
if size(pow, 1) == 1 && sum(valid_rows) > 1
    pow = repmat(pow, sum(valid_rows), 1);
end

if size(pow, 2) ~= numel(scan_freqs) && isfield(freq, 'freq') && ~isempty(freq.freq)
    freq_axis = freq.freq(:)';
    pow_interp = nan(size(pow, 1), numel(scan_freqs));
    for ri = 1:size(pow, 1)
        row_pow = pow(ri, :);
        if numel(row_pow) ~= numel(freq_axis) && size(pow, 1) == numel(freq_axis) && size(pow, 2) == size(pow_interp, 1)
            row_pow = pow(:, ri).';
        end
        if numel(row_pow) ~= numel(freq_axis)
            continue;
        end
        pow_interp(ri, :) = interp1(freq_axis, row_pow, scan_freqs, 'linear', NaN);
    end
    pow = pow_interp;
end

p_scan(valid_rows, :) = pow;
end

%% ------------------------------------------------------------------------
%% Utilities

function sim_template = build_signed_occipital_template(chan_labels, headmodel, occ_idx, front_idx)
nChans = numel(chan_labels);
template_front_weight = 0.75;
template_sigma_occ = 0.12;
template_sigma_front = 0.25;
sim_template = zeros(nChans, 1);

lay_labels = headmodel.layANThead.label;
lay_pos = headmodel.layANThead.pos;
chan_pos = nan(nChans, 2);
for ch = 1:nChans
    li = find(strcmp(lay_labels, chan_labels{ch}), 1, 'first');
    if ~isempty(li)
        chan_pos(ch, :) = lay_pos(li, :);
    end
end
has_pos = ~any(isnan(chan_pos), 2);
occ_pos_idx = intersect(occ_idx, find(has_pos));
front_pos_idx = intersect(front_idx, find(has_pos));
if ~isempty(occ_pos_idx) && ~isempty(front_pos_idx)
    occ_ctr = mean(chan_pos(occ_pos_idx, :), 1);
    front_ctr = mean(chan_pos(front_pos_idx, :), 1);
    occ_d2 = sum((chan_pos - occ_ctr).^2, 2);
    front_d2 = sum((chan_pos - front_ctr).^2, 2);
    w_occ = exp(-occ_d2 / (2 * template_sigma_occ^2));
    w_front = exp(-front_d2 / (2 * template_sigma_front^2));
    sim_template = w_occ - template_front_weight * w_front;
    sim_template(~has_pos) = 0;
else
    sim_template(occ_idx) = 1;
    sim_template(front_idx) = -template_front_weight;
end
if std(sim_template) > 0
    sim_template = (sim_template - mean(sim_template)) / std(sim_template);
end
end

function m = robust_mad(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    m = 0;
    return;
end
m = mad(x, 1);
if ~isfinite(m) || m <= eps
    m = std(x);
end
if ~isfinite(m) || m <= eps
    m = 0;
end
end

function format_power_change_db_axis(axh)
if nargin < 1 || isempty(axh) || ~isgraphics(axh, 'axes')
    axh = gca;
end
yt = get(axh, 'YTick');
if isempty(yt)
    return;
end
yt_lbl = arrayfun(@(v) sprintf('%g', v), yt, 'UniformOutput', false);
set(axh, 'YTickLabel', yt_lbl);
end

function save_figure_png(fig_handle, out_path)
if nargin < 1 || isempty(fig_handle) || ~ishandle(fig_handle)
    return;
end
if nargin < 2 || isempty(out_path)
    return;
end
out_dir = fileparts(out_path);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
drawnow;
pause(0.05);
exportgraphics(fig_handle, out_path, 'Resolution', 300);
end
