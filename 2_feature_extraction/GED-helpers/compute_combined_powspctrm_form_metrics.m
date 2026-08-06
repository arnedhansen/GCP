function [pf_score, pf_peak_hz] = compute_combined_powspctrm_form_metrics( ...
    spec_vec, scan_freqs, analysis_freq_range)
% PF score and peak frequency for a single combined GED spectrum.
pf_score = NaN;
pf_peak_hz = NaN;
if isempty(spec_vec) || isempty(scan_freqs) || numel(spec_vec) ~= numel(scan_freqs)
    return;
end
[pf_vec, ~] = compute_powspctrm_form_laplacian_score_from_spectra( ...
    spec_vec(:)', scan_freqs, analysis_freq_range);
if ~isempty(pf_vec) && isfinite(pf_vec(1))
    pf_score = pf_vec(1);
end
freq_mask = scan_freqs >= analysis_freq_range(1) & scan_freqs <= analysis_freq_range(2);
if any(freq_mask)
    x = scan_freqs(freq_mask);
    y = spec_vec(freq_mask);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    if numel(y) >= 3
        peak_lo = max(30, analysis_freq_range(1));
        peak_hi = min(90, analysis_freq_range(2));
        inner_mask = x >= peak_lo & x <= peak_hi;
        if sum(inner_mask) >= 3
            x_peak = x(inner_mask);
            y_peak = y(inner_mask);
        else
            x_peak = x;
            y_peak = y;
        end
        [~, idx_max] = max(y_peak);
        pf_peak_hz = x_peak(idx_max);
    end
end
if ~isfinite(pf_peak_hz)
    if sum(freq_mask) >= 3
        if numel(y) >= 3
            peak_lo = max(30, analysis_freq_range(1));
            peak_hi = min(90, analysis_freq_range(2));
            inner_mask = x >= peak_lo & x <= peak_hi;
            if sum(inner_mask) >= 3
                x_peak = x(inner_mask);
                y_peak = y(inner_mask);
            else
                x_peak = x;
                y_peak = y;
            end
            [~, idx_max] = max(y_peak);
            pf_peak_hz = x_peak(idx_max);
        end
    end
end
end
