function [pf_score, pf_peak_hz] = compute_combined_powspctrm_form_metrics( ...
    spec_vec, scan_freqs, analysis_freq_range)
% PF score and peak frequency for a single combined GED spectrum.
% Peak frequency uses pick_tallest_peak, the same rule as condition means.
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
[pf_peak_hz, ~] = pick_tallest_peak(spec_vec, scan_freqs, 0, 5);
end
