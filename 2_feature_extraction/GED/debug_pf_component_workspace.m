%% Debug PF computation for one component from current workspace
% Usage:
% 1) Run GCP_eeg_GED until data for subject/window are in workspace.
% 2) Set the parameters below.
% 3) Run this script.

%% User settings
window_name = 'full';   % 'full' | 'early' | 'late'
component_idx = 3;      % example: C3 from figure

%% Resolve window-specific variables
spec_var = ['searchMeanPrSpectrum_' window_name];
table_var = ['candidate_table_' window_name];

if exist(spec_var, 'var')
    spec_mat = eval(spec_var);
elseif exist('searchMeanPrSpectrum', 'var')
    spec_mat = searchMeanPrSpectrum;
else
    error('No spectrum matrix found in workspace for window "%s".', window_name);
end

if exist(table_var, 'var')
    cand = eval(table_var);
elseif exist('candidate_table', 'var')
    cand = candidate_table;
else
    error('No candidate table found in workspace for window "%s".', window_name);
end

if ~exist('scan_freqs', 'var')
    error('scan_freqs is missing in workspace.');
end
if ~exist('analysis_freq_range', 'var')
    analysis_freq_range = [30 90];
end

assert(component_idx >= 1 && component_idx <= size(spec_mat,1), 'component_idx out of range');
y_raw = spec_mat(component_idx, :);
x = scan_freqs(:)';

%% Pull PF settings (fallback to defaults if absent)
poly_order_local = get_ws_or_default('poly_order', 2);
detrend_edge_exclude_n_local = get_ws_or_default('detrend_edge_exclude_n', 5);
detrend_in_log_local = get_ws_or_default('detrend_in_log', false);
detrend_flat_edges_local = get_ws_or_default('detrend_flat_edges', true);
peak_form_smooth_n_local = max(1, round(get_ws_or_default('peak_form_smooth_n', 3)));
peak_form_shift_max_hz_local = get_ws_or_default('peak_form_shift_max_hz', 10);
peak_form_single_widths_local = get_ws_or_default('peak_form_single_widths', [5 7 9 12]);
peak_form_double_widths_local = get_ws_or_default('peak_form_double_widths', [3 4 5 6]);
peak_form_double_separations_local = get_ws_or_default('peak_form_double_separations', [8 12 16 20]);
peak_form_min_trough_depth_local = get_ws_or_default('peak_form_min_trough_depth', 0.10);
peak_form_min_similarity_local = get_ws_or_default('peak_form_min_similarity', 0.50);

%% Recompute processed trace used by PF scoring
y_dt = detrend_power_ratio_dbg(y_raw, x, poly_order_local, detrend_edge_exclude_n_local, detrend_in_log_local, detrend_flat_edges_local);
band_mask = x >= analysis_freq_range(1) & x <= analysis_freq_range(2);
x_band = x(band_mask);
y_band = y_dt(band_mask);
valid = isfinite(x_band) & isfinite(y_band);
x_band = x_band(valid);
y_band = y_band(valid);

if numel(y_band) < 7
    error('Too few valid bins (%d) after preprocessing.', numel(y_band));
end

y_band_smooth = movmean(y_band, peak_form_smooth_n_local);
y_band_smooth = apply_soft_edge_attenuation_dbg(y_band_smooth, x_band, [30 90], [40 80], 0.01, 4.0, 3.0);
y_norm = normalize_positive_shape_dbg(y_band_smooth);

%% Compute template similarities and penalties for this component
[single_best, single_meta] = evaluate_single_template_bank_dbg( ...
    y_norm, x_band, peak_form_single_widths_local, peak_form_shift_max_hz_local);
[double_best, double_meta] = evaluate_double_template_bank_dbg( ...
    y_norm, y_band_smooth, x_band, peak_form_double_widths_local, ...
    peak_form_double_separations_local, peak_form_shift_max_hz_local, peak_form_min_trough_depth_local);

mode_raw = 'single';
best_pre_penalty = single_best;
best_centers = single_meta.centers_hz;
if double_best > best_pre_penalty
    mode_raw = 'double';
    best_pre_penalty = double_best;
    best_centers = double_meta.centers_hz;
end

edge_pen = 1;
edge_margin_hz = 2;
if ~isempty(best_centers) && any(best_centers <= (analysis_freq_range(1) + edge_margin_hz) | ...
        best_centers >= (analysis_freq_range(2) - edge_margin_hz))
    edge_pen = 0.85;
end

max_y_band = max(y_band_smooth);
if ~isfinite(max_y_band), max_y_band = 0; end
[pks, locs] = findpeaks(y_band_smooth, x_band, 'MinPeakProminence', max(0, 0.20 * max_y_band));
dominance_pen = 1;
dominance_ratio = NaN;
if numel(pks) > 1
    pks_sorted = sort(pks, 'descend');
    dominance_ratio = pks_sorted(1) / max(pks_sorted(2), eps);
    if dominance_ratio < 2.0
        dominance_pen = dominance_pen * (0.55 + 0.45 * min(1, (dominance_ratio - 1) / 1.0));
    end
    if numel(pks) >= 3
        n_secondary = numel(pks) - 1;
        dominance_pen = dominance_pen * max(0.65, 1 - 0.08 * (n_secondary - 1));
    end
elseif isempty(pks)
    dominance_pen = dominance_pen * 0.40;
end

hf_pen = 1;
hf_mask = x_band >= max(70, analysis_freq_range(2) - 15);
if sum(hf_mask) >= 5
    hf_idx = find(hf_mask);
    hf_rho = corr((1:numel(hf_idx))', y_band_smooth(hf_idx)', 'rows', 'complete', 'type', 'Spearman');
    if isfinite(hf_rho) && hf_rho > 0.70
        hf_pen = max(0.65, 1 - 0.30 * (hf_rho - 0.70) / 0.30);
    end
end

roughness_pen = 1;
roughness_ratio = NaN;
y_finite = y_band_smooth(isfinite(y_band_smooth));
if numel(y_finite) >= 3
    y_range = max(y_finite) - min(y_finite);
    if y_range > eps
        y_scaled = y_band_smooth / y_range;
    else
        y_scaled = y_band_smooth;
    end
    dy = diff(y_scaled);
    amp_scale = prctile(y_scaled, 75) - prctile(y_scaled, 25);
    if ~isfinite(amp_scale) || amp_scale <= eps
        amp_scale = std(y_scaled(isfinite(y_scaled)));
    end
    if ~isfinite(amp_scale) || amp_scale <= eps
        amp_scale = 1;
    end
    roughness_ratio = robust_mad_dbg(dy) / amp_scale;
    if isfinite(roughness_ratio) && roughness_ratio > 0.50
        loss_frac = min(1, (roughness_ratio - 0.50) / 0.80);
        roughness_pen = max(0.70, 1 - 0.30 * loss_frac);
    end
end

penalty_raw = edge_pen * dominance_pen * hf_pen * roughness_pen;
penalty_used = max(0.35, penalty_raw);
best_raw = best_pre_penalty * penalty_used;

if best_raw <= peak_form_min_similarity_local
    below = peak_form_min_similarity_local - best_raw;
    pf_score_dbg = 0.20 * exp(-4 * below / max(peak_form_min_similarity_local, eps));
else
    score_high = (best_raw - peak_form_min_similarity_local) / max(1 - peak_form_min_similarity_local, eps);
    pf_score_dbg = 0.20 + 0.80 * score_high;
end
pf_score_dbg = max(0, min(1, pf_score_dbg));

%% Read pipeline-recorded diagnostics for comparison
pf_tbl = cand.peak_form_score(component_idx);
mode_tbl = string(cand.peak_form_mode(component_idx));
pre_tbl = cand.peak_form_similarity_pre_penalty(component_idx);
sim_tbl = cand.peak_form_similarity_raw(component_idx);
pen_tbl = cand.peak_form_total_penalty_used(component_idx);
dom_tag_tbl = string(cand.peak_form_dominant_penalty(component_idx));

%% Print concise diagnosis
fprintf('\n=== PF DEBUG: window=%s, component=C%d ===\n', window_name, component_idx);
fprintf('Pipeline PF score:      %.4f\n', pf_tbl);
fprintf('Recomputed PF score:    %.4f\n', pf_score_dbg);
fprintf('Mode (table/recompute): %s / %s\n', mode_tbl, mode_raw);
fprintf('Pre-penalty sim (tbl):  %.4f\n', pre_tbl);
fprintf('Post-penalty sim (tbl): %.4f\n', sim_tbl);
fprintf('Penalty used (tbl):     %.4f\n', pen_tbl);
fprintf('Dominant penalty tag:   %s\n', dom_tag_tbl);
fprintf('--- Recomputed penalties ---\n');
fprintf('edge_pen=%.3f, dominance_pen=%.3f, hf_pen=%.3f, roughness_pen=%.3f\n', ...
    edge_pen, dominance_pen, hf_pen, roughness_pen);
fprintf('dominance peaks=%d, dominance_ratio=%.3f\n', numel(pks), dominance_ratio);
fprintf('single_best=%.4f, double_best=%.4f, best_pre=%.4f, penalty_used=%.4f\n', ...
    single_best, double_best, best_pre_penalty, penalty_used);

%% Plot traces and detected prominent peaks
figure('Position', [0 0 1512 982]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
plot(x, y_raw, 'k-', 'LineWidth', 1.3);
xline(40, 'k:'); xline(80, 'k:');
xline(30, 'k--'); xline(90, 'k--');
xlim([20 100]); xlabel('Hz'); ylabel('Raw PR');
title(sprintf('C%d raw spectrum', component_idx), 'Interpreter', 'none');
box on;

nexttile; hold on;
plot(x_band, y_band_smooth, '-', 'Color', [0.85 0.20 0.20], 'LineWidth', 1.2);
plot(x_band, y_norm, '-', 'Color', [0.15 0.45 0.90], 'LineWidth', 1.2);
if ~isempty(pks)
    scatter(locs, pks, 45, 'm', 'filled');
end
xline(40, 'k:'); xline(80, 'k:');
xline(30, 'k--'); xline(90, 'k--');
xlim([20 100]); xlabel('Hz'); ylabel('Processed');
legend({'detrended+tapered+smoothed', 'normalized', 'dominance peaks'}, 'Location', 'best');
title(sprintf('PF internals | mode=%s | PF(tbl)=%.3f | PF(recomp)=%.3f', ...
    mode_raw, pf_tbl, pf_score_dbg), 'Interpreter', 'none');
box on;

%% ---------- local helpers ----------
function v = get_ws_or_default(varname, fallback)
if evalin('base', sprintf('exist(''%s'',''var'')', varname))
    v = evalin('base', varname);
else
    v = fallback;
end
end

function y_dt = detrend_power_ratio_dbg(y, x, ord, edge_exclude_n, in_log, flat_edges)
y_fit = y;
if in_log
    y_fit = log(max(y, eps));
end
y_dt = detrend_poly_stable_dbg(y_fit, x, ord, edge_exclude_n, flat_edges);
y_dt = apply_soft_edge_attenuation_dbg(y_dt, x, [30 90], [40 80], 0.01, 4.0, 3.0);
end

function y_dt = detrend_poly_stable_dbg(y, x, ord, edge_exclude_n, flat_edges)
y_dt = y;
edge_exclude_n = max(0, round(edge_exclude_n));
n = numel(x);
fit_idx = (1 + edge_exclude_n):(n - edge_exclude_n);
if numel(fit_idx) < (ord + 1), fit_idx = 1:n; end
fit_idx = fit_idx(fit_idx <= min(numel(y), numel(x)) & fit_idx >= 1);
yv = y(fit_idx); xv = x(fit_idx);
valid_mask = isfinite(yv(:)) & isfinite(xv(:));
valid_fit = fit_idx(valid_mask);
if numel(valid_fit) < (ord + 1), valid_fit = find(isfinite(y(:)) & isfinite(x(:))); end
if numel(valid_fit) < (ord + 1), return; end
p = polyfit(x(valid_fit), y(valid_fit), ord);
baseline = polyval(p, x);
baseline = reshape(baseline, size(y));
y_dt = y - baseline;
if flat_edges && n >= 2
    x1 = x(1); xn = x(n);
    if isfinite(x1) && isfinite(xn) && xn > x1
        d1 = y_dt(1); dn = y_dt(n);
        if isfinite(d1) && isfinite(dn)
            ramp = d1 + (dn - d1) * (x(:) - x1) / (xn - x1);
            ramp = reshape(ramp, size(y_dt));
            y_dt = y_dt - ramp;
        end
    end
end
end

function y_out = apply_soft_edge_attenuation_dbg(y_in, x, outer_band, core_band, edge_floor, shoulder_pow, outside_decay_hz)
y_out = y_in;
if isempty(y_in) || isempty(x) || numel(y_in) ~= numel(x), return; end
edge_floor = max(0, min(1, edge_floor));
shoulder_pow = max(1, shoulder_pow);
outside_decay_hz = max(eps, outside_decay_hz);
outer_lo = min(outer_band); outer_hi = max(outer_band);
core_lo = min(core_band); core_hi = max(core_band);
w = zeros(size(x));
inside_outer = (x >= outer_lo) & (x <= outer_hi);
inside_core = (x >= core_lo) & (x <= core_hi);
w(inside_core) = 1;
shoulder_mask = inside_outer & ~inside_core;
if any(shoulder_mask)
    d_to_outer = min(abs(x(shoulder_mask) - outer_lo), abs(x(shoulder_mask) - outer_hi));
    shoulder_half_width = max(eps, min(core_lo - outer_lo, outer_hi - core_hi));
    t = max(0, min(1, d_to_outer / shoulder_half_width));
    w(shoulder_mask) = edge_floor + (1 - edge_floor) * (t .^ shoulder_pow);
end
outside_mask = ~inside_outer;
if any(outside_mask)
    d_out = min(abs(x(outside_mask) - outer_lo), abs(x(outside_mask) - outer_hi));
    w(outside_mask) = edge_floor * exp(-(d_out / outside_decay_hz) .^ 2);
end
y_out = y_in .* w;
end

function [best_sim, meta] = evaluate_single_template_bank_dbg(y_norm, x_band, widths_hz, shift_max_hz)
best_sim = 0; meta = struct('shift_hz', NaN, 'width_hz', NaN, 'centers_hz', []);
x_mid = mean(x_band); x_min = min(x_band); x_max = max(x_band);
shift_vals = candidate_shift_values_hz_dbg(x_band, shift_max_hz);
for wi = 1:numel(widths_hz)
    w = widths_hz(wi);
    if ~isfinite(w) || w <= 0, continue; end
    for si = 1:numel(shift_vals)
        center_hz = x_mid + shift_vals(si);
        if ~isfinite(center_hz) || center_hz < x_min || center_hz > x_max, continue; end
        sim = safe_template_similarity_dbg(y_norm, gaussian_template_dbg(x_band, center_hz, w));
        if sim > best_sim
            best_sim = sim;
            meta.shift_hz = shift_vals(si);
            meta.width_hz = w;
            meta.centers_hz = center_hz;
        end
    end
end
end

function [best_sim, meta] = evaluate_double_template_bank_dbg(y_norm, y_band, x_band, widths_hz, separations_hz, shift_max_hz, min_trough_depth)
best_sim = 0;
meta = struct('shift_hz', NaN, 'width_hz', NaN, 'sep_hz', NaN, 'trough_depth', NaN, 'centers_hz', []);
x_mid = mean(x_band); x_min = min(x_band); x_max = max(x_band);
shift_vals = candidate_shift_values_hz_dbg(x_band, shift_max_hz);
for wi = 1:numel(widths_hz)
    w = widths_hz(wi);
    if ~isfinite(w) || w <= 0, continue; end
    for di = 1:numel(separations_hz)
        sep = separations_hz(di);
        if ~isfinite(sep) || sep <= 0, continue; end
        for si = 1:numel(shift_vals)
            c1 = x_mid + shift_vals(si) - sep/2;
            c2 = x_mid + shift_vals(si) + sep/2;
            if ~isfinite(c1) || ~isfinite(c2) || c1 < x_min || c2 > x_max, continue; end
            tpl = gaussian_template_dbg(x_band, c1, w) + gaussian_template_dbg(x_band, c2, w);
            sim = safe_template_similarity_dbg(y_norm, tpl);
            trough_depth = estimate_trough_depth_dbg(y_band, x_band, c1, c2);
            sim_adj = sim * min(1, max(0, trough_depth) / max(min_trough_depth, eps));
            if sim_adj > best_sim
                best_sim = sim_adj;
                meta.shift_hz = shift_vals(si); meta.width_hz = w; meta.sep_hz = sep;
                meta.trough_depth = trough_depth; meta.centers_hz = [c1 c2];
            end
        end
    end
end
end

function shifts = candidate_shift_values_hz_dbg(x_band, shift_max_hz)
if numel(x_band) >= 2, df = median(diff(x_band)); else, df = 1; end
if ~isfinite(df) || df <= 0, df = 1; end
shift_max_hz = max(0, shift_max_hz);
shift_max_hz = max(shift_max_hz, max(abs(x_band - mean(x_band))));
shifts = -shift_max_hz:df:shift_max_hz;
if isempty(shifts), shifts = 0; end
end

function depth = estimate_trough_depth_dbg(y_band, x_band, c1, c2)
depth = 0;
if c1 > c2, tmp = c1; c1 = c2; c2 = tmp; end
[~, i1] = min(abs(x_band - c1)); [~, i2] = min(abs(x_band - c2));
if i1 == i2, return; end
idx_lo = min(i1, i2); idx_hi = max(i1, i2);
y_pos = max(y_band(:), 0);
p1 = y_pos(i1); p2 = y_pos(i2);
if idx_hi - idx_lo < 2 || ~isfinite(p1) || ~isfinite(p2), return; end
valley = min(y_pos(idx_lo:idx_hi));
peak_ref = max(min(p1, p2), eps);
depth = max(0, min(1, 1 - valley / peak_ref));
end

function y_norm = normalize_positive_shape_dbg(y)
y_norm = [];
if isempty(y), return; end
y = y(:);
y = y - median(y(isfinite(y)));
y(~isfinite(y)) = 0;
y = max(y, 0);
y_mag = norm(y);
if ~isfinite(y_mag) || y_mag <= eps, return; end
y_norm = y / y_mag;
end

function sim = safe_template_similarity_dbg(y_norm, tpl)
sim = 0;
if isempty(y_norm) || isempty(tpl), return; end
tpl = tpl(:); tpl = max(tpl, 0);
tpl_mag = norm(tpl);
if ~isfinite(tpl_mag) || tpl_mag <= eps, return; end
tpl = tpl / tpl_mag;
if numel(tpl) ~= numel(y_norm), return; end
sim = y_norm(:)' * tpl(:);
if ~isfinite(sim), sim = 0; end
sim = max(0, min(1, sim));
end

function g = gaussian_template_dbg(x, mu, sigma)
if ~isfinite(mu) || ~isfinite(sigma) || sigma <= 0
    g = zeros(size(x)); return;
end
g = exp(-0.5 * ((x - mu) ./ sigma).^2);
end

function m = robust_mad_dbg(x)
x = x(:); x = x(isfinite(x));
if isempty(x), m = 0; return; end
m = mad(x, 1);
if ~isfinite(m) || m <= eps, m = std(x); end
if ~isfinite(m) || m <= eps, m = 0; end
end
