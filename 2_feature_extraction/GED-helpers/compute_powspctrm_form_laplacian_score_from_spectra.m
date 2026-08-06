function [pf_score_vec, pf_deduction_vec, diag] = compute_powspctrm_form_laplacian_score_from_spectra( ...
    mean_pr_spectrum, scan_freqs, analysis_freq_range)
% Spectral peak-form (PF) scores via Laplace-template fit, with multi-peak deduction.
rival_sep_min_hz = 10;
rival_height_ratio_min = 0.90;
deduct_per_rival = 0.05;
deduct_per_rival = max(0, min(1, deduct_per_rival));
nComp = size(mean_pr_spectrum, 1);
pf_score_vec = zeros(nComp, 1);
pf_deduction_vec = zeros(nComp, 1);
diag = struct( ...
    'laplace_r2', zeros(nComp, 1), ...
    'dominant_peak_hz', nan(nComp, 1), ...
    'n_rival_peaks', zeros(nComp, 1));
if isempty(mean_pr_spectrum) || isempty(scan_freqs)
    return;
end
freq_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
if ~any(freq_mask)
    return;
end
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
    y_smooth = movmean(y, 5, 'omitnan');
    y_floor = prctile(y_smooth, 20);
    if ~isfinite(y_floor)
        y_floor = median(y_smooth, 'omitnan');
    end
    if ~isfinite(y_floor)
        y_floor = 0;
    end
    y_pos = max(y_smooth - y_floor, 0);
    y_scale = max(y_pos);
    if ~isfinite(y_scale) || y_scale <= eps
        continue;
    end
    y_norm = y_pos / y_scale;
    [dom_amp, dom_idx] = max(y_pos);
    if ~isfinite(dom_amp) || dom_amp <= eps || isempty(dom_idx)
        continue;
    end
    dom_hz = x(dom_idx);
    diag.dominant_peak_hz(ci) = dom_hz;

    p0 = [1.0, dom_hz, 6.0];
    lb = [0.20, min(x), 1.0];
    ub = [1.50, max(x), 20.0];
    obj = @(p) mean((y_norm - laplace_template_model(x, p)).^2, 'omitnan');
    opts = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 3000);
    p_est = fminsearch(@(pp) bounded_obj(pp, lb, ub, obj), p0, opts);
    p_est = min(max(p_est, lb), ub);
    y_hat = laplace_template_model(x, p_est);
    r2 = compute_r2_score(y_norm, y_hat);
    diag.laplace_r2(ci) = r2;

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
    min_prom = max([0, rel_prom, 0.02, 0.15 * robust_scale]);
    [pks, locs] = findpeaks(y_pos, x, 'MinPeakProminence', min_prom, 'MinPeakDistance', 5);
    n_rivals = 0;
    if ~isempty(pks) && numel(pks) >= 2
        pks = pks(:);
        locs = locs(:);
        dom_mask = abs(locs - dom_hz) <= 1e-9;
        if ~any(dom_mask)
            [~, nearest_idx] = min(abs(locs - dom_hz));
            dom_mask = false(size(locs));
            dom_mask(nearest_idx) = true;
        end
        rival_mask = ~dom_mask & ...
            (abs(locs - dom_hz) >= rival_sep_min_hz) & ...
            (pks >= rival_height_ratio_min * dom_amp);
        n_rivals = sum(rival_mask);
    end
    diag.n_rival_peaks(ci) = n_rivals;
    deduct_val = min(r2, n_rivals * deduct_per_rival);
    pf_deduction_vec(ci) = deduct_val;
    pf_score_vec(ci) = max(0, min(1, r2 - deduct_val));
end
pf_score_vec(~isfinite(pf_score_vec)) = 0;
pf_deduction_vec(~isfinite(pf_deduction_vec)) = 0;
end
