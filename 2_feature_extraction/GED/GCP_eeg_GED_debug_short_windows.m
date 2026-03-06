%% GED Debug: short-window component stability
%
% This script isolates Phase-1 GED behavior for full/early/late windows.
% It compares:
%   1) shared baseline regularization (current pipeline behavior)
%   2) window-matched baseline regularization (diagnostic alternative)
%
% Output:
%   - command window tables with conditioning and component diagnostics
%   - topography grids for the first N components per window/mode
%   - eigenvalue and degeneracy summary plots

%% Setup
startup
[subjects, path, ~, headmodel] = setup('GCP');

%% User settings
subject_selector = '601';   % subject ID string (e.g., '601') or numeric subject index
trial_codes = [61, 62, 63, 64];

baseline_window = [-1.5, -0.25];
full_window = [0.0, 2.0];
early_window = [0.0, 0.6];
late_window = [1.0, 2.0];
win_defs = {full_window, early_window, late_window};
win_names = {'full', 'early', 'late'};

gamma_range = [30, 90];
lambda_shared = 0.05;
lambda_per_win = [0.05, 0.10, 0.08];
n_components_debug = 8; % components visualized and scored per window

save_figures = true;
if ispc
    gcp_root_path = 'W:\Students\Arne\GCP';
else
    gcp_root_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP';
end
debug_fig_dir = fullfile(gcp_root_path, 'figures', 'eeg', 'ged', 'debug_short_windows');
if save_figures && ~exist(debug_fig_dir, 'dir')
    mkdir(debug_fig_dir);
end

%% Subject resolution
if isnumeric(subject_selector)
    subj_idx = round(subject_selector);
    if subj_idx < 1 || subj_idx > numel(subjects)
        error('subject_selector index %d out of range [1..%d].', subj_idx, numel(subjects));
    end
else
    subj_idx = find(strcmp(subjects, subject_selector), 1, 'first');
    if isempty(subj_idx)
        error('Subject "%s" not found in setup("GCP") list.', subject_selector);
    end
end
subject_id = subjects{subj_idx};
fprintf('\n=== GED short-window debug: subject %s (%d/%d) ===\n', subject_id, subj_idx, numel(subjects));

%% Load subject data
datapath = fullfile(path, subject_id, 'eeg');
eeg_data = load(fullfile(datapath, 'dataEEG.mat'), 'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
dataEEG_c25 = eeg_data.dataEEG_c25;
dataEEG_c50 = eeg_data.dataEEG_c50;
dataEEG_c75 = eeg_data.dataEEG_c75;
dataEEG_c100 = eeg_data.dataEEG_c100;
data_structs = {dataEEG_c25, dataEEG_c50, dataEEG_c75, dataEEG_c100};

fsample = dataEEG_c25.fsample;
n_chans = numel(dataEEG_c25.label);

% Channel groups for spatial diagnostics
occ_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO)', 'once')), dataEEG_c25.label);
post_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|P)', 'once')), dataEEG_c25.label);
if ~any(post_mask)
    post_mask = occ_mask;
end
if ~any(post_mask)
    post_mask(:) = true;
end
nonpost_mask = ~post_mask;
if ~any(nonpost_mask)
    nonpost_mask = true(size(post_mask));
end
center_mask = compute_center_channel_mask(dataEEG_c25.label, headmodel.layANThead);

%% Build pooled covariances per window
cov_stim_per_win = {zeros(n_chans), zeros(n_chans), zeros(n_chans)};
cov_base = zeros(n_chans);
n_trials_total = 0;
n_trials_per_cond = zeros(1, numel(data_structs));
n_samples_base = 0;
n_samples_stim_per_win = zeros(1, 3);

for cond = 1:numel(data_structs)
    dat = data_structs{cond};
    trl_idx = find(dat.trialinfo == trial_codes(cond));
    if isempty(trl_idx)
        continue;
    end
    cfg = [];
    cfg.trials = trl_idx;
    dat = ft_selectdata(cfg, dat);

    cfg_f = [];
    cfg_f.bpfilter = 'yes';
    cfg_f.bpfreq = gamma_range;
    cfg_f.bpfilttype = 'fir';
    cfg_f.bpfiltord = round(3 * fsample / gamma_range(1));
    dat_gamma = ft_preprocessing(cfg_f, dat);

    cfg_t = [];
    cfg_t.latency = baseline_window;
    dat_base = ft_selectdata(cfg_t, dat_gamma);

    cfg_t.latency = full_window;
    dat_full = ft_selectdata(cfg_t, dat_gamma);
    cfg_t.latency = early_window;
    dat_early = ft_selectdata(cfg_t, dat_gamma);
    cfg_t.latency = late_window;
    dat_late = ft_selectdata(cfg_t, dat_gamma);

    n_trl = numel(dat_full.trial);
    n_trials_per_cond(cond) = n_trl;

    for trl = 1:n_trl
        xb = double(dat_base.trial{trl});
        xb = xb - mean(xb, 2);
        cov_base = cov_base + (xb * xb') / size(xb, 2);
        n_samples_base = n_samples_base + size(xb, 2);

        xf = double(dat_full.trial{trl});
        xf = xf - mean(xf, 2);
        cov_stim_per_win{1} = cov_stim_per_win{1} + (xf * xf') / size(xf, 2);
        n_samples_stim_per_win(1) = n_samples_stim_per_win(1) + size(xf, 2);

        xe = double(dat_early.trial{trl});
        xe = xe - mean(xe, 2);
        cov_stim_per_win{2} = cov_stim_per_win{2} + (xe * xe') / size(xe, 2);
        n_samples_stim_per_win(2) = n_samples_stim_per_win(2) + size(xe, 2);

        xl = double(dat_late.trial{trl});
        xl = xl - mean(xl, 2);
        cov_stim_per_win{3} = cov_stim_per_win{3} + (xl * xl') / size(xl, 2);
        n_samples_stim_per_win(3) = n_samples_stim_per_win(3) + size(xl, 2);
    end

    n_trials_total = n_trials_total + n_trl;
end

if n_trials_total < 1
    error('No valid trials found for subject %s.', subject_id);
end

for w = 1:3
    cov_stim_per_win{w} = cov_stim_per_win{w} / n_trials_total;
end
cov_base = cov_base / n_trials_total;

fprintf('Trials per condition [c25 c50 c75 c100] = [%d %d %d %d], total=%d\n', ...
    n_trials_per_cond(1), n_trials_per_cond(2), n_trials_per_cond(3), n_trials_per_cond(4), n_trials_total);
fprintf('Samples in windows (total over trials): base=%d, full=%d, early=%d, late=%d\n', ...
    n_samples_base, n_samples_stim_per_win(1), n_samples_stim_per_win(2), n_samples_stim_per_win(3));

%% GED debug modes
modes = {'shared_baseline', 'window_matched_baseline'};
results = struct();

for mi = 1:numel(modes)
    mode_name = modes{mi};
    results.(mode_name) = struct();
    fprintf('\n--- Mode: %s ---\n', mode_name);

    if strcmp(mode_name, 'shared_baseline')
        cov_base_reg_shared = regularize_cov(cov_base, lambda_shared);
    end

    for w = 1:3
        lam_w = lambda_per_win(w);
        cov_stim_w = cov_stim_per_win{w};
        cov_stim_reg = regularize_cov(cov_stim_w, lam_w);

        if strcmp(mode_name, 'shared_baseline')
            cov_base_reg = cov_base_reg_shared;
        else
            cov_base_reg = regularize_cov(cov_base, lam_w);
        end

        [W, D] = eig(symmetrize_cov(cov_stim_reg), symmetrize_cov(cov_base_reg));
        [evals_sorted, ord] = sort(real(diag(D)), 'descend');
        W = real(W(:, ord));

        n_keep = min(n_components_debug, numel(evals_sorted));
        comp_idx = (1:n_keep)';
        eigvals = evals_sorted(1:n_keep);

        topo_mat = nan(n_chans, n_keep);
        residual = nan(n_keep, 1);
        center_ratio = nan(n_keep, 1);
        post_nonpost_ratio = nan(n_keep, 1);
        topo_std = nan(n_keep, 1);
        topo_span = nan(n_keep, 1);

        for ci = 1:n_keep
            w_ci = W(:, ci);
            topo_ci = cov_stim_reg * w_ci;
            topo_mat(:, ci) = topo_ci;

            num = norm(cov_stim_reg * w_ci - eigvals(ci) * cov_base_reg * w_ci, 2);
            den = max(norm(cov_stim_reg * w_ci, 2), eps);
            residual(ci) = num / den;

            [center_ratio(ci), post_nonpost_ratio(ci), topo_std(ci), topo_span(ci)] = ...
                score_topo_degeneracy(topo_ci, center_mask, post_mask, nonpost_mask);
        end

        stats_table = table(comp_idx, eigvals, residual, center_ratio, post_nonpost_ratio, topo_std, topo_span, ...
            'VariableNames', {'component', 'eigenvalue', 'eig_residual', 'center_ratio', ...
            'post_nonpost_ratio', 'topo_std', 'topo_span'});
        disp(stats_table);

        diag_struct = struct();
        diag_struct.window = win_names{w};
        diag_struct.lambda = lam_w;
        diag_struct.evals_sorted = evals_sorted;
        diag_struct.W = W;
        diag_struct.topos = topo_mat;
        diag_struct.stats_table = stats_table;
        diag_struct.n_eig_gt1 = sum(evals_sorted > 1);
        diag_struct.cond_stim_raw = safe_cond(cov_stim_w);
        diag_struct.cond_stim_reg = safe_cond(cov_stim_reg);
        diag_struct.cond_base_reg = safe_cond(cov_base_reg);
        diag_struct.rank_stim_raw = rank(cov_stim_w);
        diag_struct.rank_stim_reg = rank(cov_stim_reg);
        diag_struct.rank_base_reg = rank(cov_base_reg);

        results.(mode_name).(win_names{w}) = diag_struct;

        fprintf('Window=%s | lambda=%.3f | eig>1=%d | cond(stimRaw)=%.3e | cond(stimReg)=%.3e | cond(baseReg)=%.3e\n', ...
            win_names{w}, lam_w, diag_struct.n_eig_gt1, diag_struct.cond_stim_raw, ...
            diag_struct.cond_stim_reg, diag_struct.cond_base_reg);
    end
end

%% Plot: topography grids by mode
for mi = 1:numel(modes)
    mode_name = modes{mi};
    fig_topo = figure('Position', [0 0 1512 982]);
    tiledlayout(3, n_components_debug, 'Padding', 'compact', 'TileSpacing', 'compact');
    for w = 1:3
        ws = results.(mode_name).(win_names{w});
        n_show = min(n_components_debug, size(ws.topos, 2));
        for ci = 1:n_components_debug
            nexttile;
            if ci <= n_show
                topo_data = [];
                topo_data.label = dataEEG_c25.label;
                topo_data.avg = ws.topos(:, ci);
                topo_data.dimord = 'chan';
                cfg_topo = [];
                cfg_topo.layout = headmodel.layANThead;
                cfg_topo.comment = 'no';
                cfg_topo.marker = 'off';
                cfg_topo.style = 'straight';
                cfg_topo.gridscale = 300;
                cfg_topo.colormap = '*RdBu';
                cfg_topo.figure = 'gcf';
                topo_lim = max(abs(topo_data.avg(isfinite(topo_data.avg))));
                if ~isfinite(topo_lim) || topo_lim <= 0
                    topo_lim = 1;
                end
                cfg_topo.zlim = [-topo_lim topo_lim];
                try
                    ft_topoplotER(cfg_topo, topo_data);
                catch
                    imagesc(topo_data.avg(:)); axis tight;
                    caxis([-topo_lim topo_lim]);
                end
                if ci == 1
                    title(sprintf('%s C%d', win_names{w}, ci), 'FontSize', 8, 'Interpreter', 'none');
                else
                    title(sprintf('C%d', ci), 'FontSize', 8, 'Interpreter', 'none');
                end
            else
                axis off;
            end
        end
    end
    sgtitle(sprintf('GED topographies (%s) | subject %s', mode_name, subject_id), 'FontWeight', 'bold');
    if save_figures
        saveas(fig_topo, fullfile(debug_fig_dir, sprintf('GED_debug_topos_subj%s_%s.png', subject_id, mode_name)));
    end
end

%% Plot: eigen spectra and degeneracy metrics
fig_metrics = figure('Position', [0 0 1512 982]);
tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for w = 1:3
    ws_shared = results.shared_baseline.(win_names{w});
    ws_match = results.window_matched_baseline.(win_names{w});

    nexttile((w - 1) * 3 + 1);
    k = min(20, numel(ws_shared.evals_sorted));
    plot(1:k, ws_shared.evals_sorted(1:k), 'o-', 'LineWidth', 1.2); hold on;
    km = min(20, numel(ws_match.evals_sorted));
    plot(1:km, ws_match.evals_sorted(1:km), 's-', 'LineWidth', 1.2);
    yline(1, 'k--');
    xlabel('Component');
    ylabel('\lambda');
    title(sprintf('%s eigenvalues', win_names{w}), 'Interpreter', 'none');
    legend({'shared baseline', 'matched baseline', '\lambda=1'}, 'Location', 'best', 'Box', 'off');
    box on;

    nexttile((w - 1) * 3 + 2);
    bar([ ...
        median(ws_shared.stats_table.center_ratio, 'omitnan'), ...
        median(ws_match.stats_table.center_ratio, 'omitnan'); ...
        median(ws_shared.stats_table.post_nonpost_ratio, 'omitnan'), ...
        median(ws_match.stats_table.post_nonpost_ratio, 'omitnan')]);
    set(gca, 'XTickLabel', {'center_ratio', 'post_nonpost_ratio'});
    ylabel('Median across displayed components');
    title(sprintf('%s spatial metrics', win_names{w}), 'Interpreter', 'none');
    legend({'shared', 'matched'}, 'Location', 'best', 'Box', 'off');
    box on;

    nexttile((w - 1) * 3 + 3);
    bar([ ...
        ws_shared.cond_stim_raw, ws_shared.cond_stim_reg, ws_shared.cond_base_reg; ...
        ws_match.cond_stim_raw, ws_match.cond_stim_reg, ws_match.cond_base_reg]);
    set(gca, 'XTickLabel', {'shared', 'matched'});
    set(gca, 'YScale', 'log');
    ylabel('Condition number (log scale)');
    title(sprintf('%s conditioning', win_names{w}), 'Interpreter', 'none');
    legend({'stim raw', 'stim reg', 'base reg'}, 'Location', 'best', 'Box', 'off');
    box on;
end

sgtitle(sprintf('GED short-window diagnostics | subject %s', subject_id), 'FontWeight', 'bold');
if save_figures
    saveas(fig_metrics, fullfile(debug_fig_dir, sprintf('GED_debug_metrics_subj%s.png', subject_id)));
end

%% Save debug structure
debug_out = struct();
debug_out.subject_id = subject_id;
debug_out.windows = win_names;
debug_out.lambda_shared = lambda_shared;
debug_out.lambda_per_win = lambda_per_win;
debug_out.n_trials_total = n_trials_total;
debug_out.n_trials_per_cond = n_trials_per_cond;
debug_out.n_samples_base = n_samples_base;
debug_out.n_samples_stim_per_win = n_samples_stim_per_win;
debug_out.results = results;

if save_figures
    save(fullfile(debug_fig_dir, sprintf('GED_debug_struct_subj%s.mat', subject_id)), 'debug_out');
end

fprintf('\nDebug finished for subject %s.\n', subject_id);
fprintf('Check center_ratio and post_nonpost_ratio in early/late for degeneration signatures.\n');

%% Local helpers
function Creg = regularize_cov(C, lam)
C = symmetrize_cov(C);
mu = mean(diag(C), 'omitnan');
if ~isfinite(mu)
    mu = 1;
end
Creg = (1 - lam) * C + lam * mu * eye(size(C));
Creg = symmetrize_cov(Creg);
end

function C = symmetrize_cov(C)
C = real((C + C') / 2);
end

function k = safe_cond(C)
C = symmetrize_cov(C);
try
    k = cond(C);
catch
    s = svd(C);
    s = s(isfinite(s));
    if isempty(s) || max(s) <= 0
        k = Inf;
    else
        k = max(s) / max(min(s), eps);
    end
end
if ~isfinite(k)
    k = Inf;
end
end

function center_mask = compute_center_channel_mask(chan_labels, layout)
n_chans = numel(chan_labels);
center_mask = false(n_chans, 1);
if isempty(layout) || ~isfield(layout, 'label') || ~isfield(layout, 'pos')
    center_mask(:) = true;
    return;
end
layout_labels = layout.label;
layout_pos = layout.pos;
chan_pos = nan(n_chans, 2);
for ch = 1:n_chans
    idx = find(strcmp(layout_labels, chan_labels{ch}), 1, 'first');
    if ~isempty(idx)
        chan_pos(ch, :) = layout_pos(idx, :);
    end
end
valid = ~any(isnan(chan_pos), 2);
if ~any(valid)
    center_mask(:) = true;
    return;
end
d = sqrt(sum(chan_pos(valid, :).^2, 2));
thr = prctile(d, 35);
if ~isfinite(thr) || thr <= 0
    thr = median(d);
end
center_mask(valid) = d <= thr;
if ~any(center_mask)
    [~, min_idx] = min(d);
    valid_idx = find(valid);
    center_mask(valid_idx(min_idx)) = true;
end
end

function [center_ratio, post_nonpost_ratio, topo_std, topo_span] = score_topo_degeneracy(topo_vec, center_mask, post_mask, nonpost_mask)
topo_vec = topo_vec(:);
valid = isfinite(topo_vec);
if ~any(valid)
    center_ratio = NaN;
    post_nonpost_ratio = NaN;
    topo_std = NaN;
    topo_span = NaN;
    return;
end
abs_topo = abs(topo_vec);
global_mean = mean(abs_topo(valid));
if ~isfinite(global_mean) || global_mean <= eps
    global_mean = eps;
end

center_idx = valid & center_mask;
if ~any(center_idx)
    center_idx = valid;
end
center_ratio = mean(abs_topo(center_idx)) / global_mean;

post_idx = valid & post_mask;
nonpost_idx = valid & nonpost_mask;
if ~any(post_idx)
    post_idx = valid;
end
if ~any(nonpost_idx)
    nonpost_idx = valid;
end
post_nonpost_ratio = mean(abs_topo(post_idx)) / max(mean(abs_topo(nonpost_idx)), eps);

vals = topo_vec(valid);
topo_std = std(vals);
topo_span = max(vals) - min(vals);
end
