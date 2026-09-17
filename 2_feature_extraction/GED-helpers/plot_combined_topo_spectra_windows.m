function plot_combined_topo_spectra_windows(save_dir, subject_id, scan_freqs, cfg_topo, topo_labels, ...
    searchTopos_full, searchMeanPrSpectrum_full, selected_idx_full, w_combined_full, ...
    analysis_freq_range, file_tag)
% Subject figure of combined GED topo and spectrum for the full stimulus window.
if nargin < 11 || isempty(file_tag)
    file_tag = 'full';
end
fig = figure('Position', [0 0 1512 982], 'Color', 'w');
[sel_idx, sel_w] = sanitize_selected_components( ...
    selected_idx_full, w_combined_full, size(searchMeanPrSpectrum_full, 1));

subplot(1, 2, 1);
if isempty(sel_idx) || isempty(searchTopos_full)
    axis off;
    text(0.5, 0.5, 'No combined components (FULL)', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 11, 'Color', [0.7 0.1 0.1], 'Interpreter', 'none');
else
    topo_vec = searchTopos_full(:, sel_idx) * sel_w(:);
    topo_data = [];
    topo_data.label = topo_labels;
    topo_data.avg = topo_vec;
    topo_data.dimord = 'chan';
    topo_vals = topo_vec(isfinite(topo_vec));
    topo_clim = max(abs(topo_vals));
    if ~isfinite(topo_clim) || topo_clim <= 0
        topo_clim = 1;
    end
    cfg_ci = cfg_topo;
    cfg_ci.zlim = [-topo_clim topo_clim];
    try
        ft_topoplotER(cfg_ci, topo_data);
    catch
        imagesc(topo_vec(:)); axis tight;
        caxis([-topo_clim topo_clim]); colorbar;
    end
    title(sprintf('Combined Topography (FULL, n=%d)', numel(sel_idx)), ...
        'FontSize', 11, 'Interpreter', 'none');
end
set(gca, 'FontSize', 10);

subplot(1, 2, 2); hold on;
if isempty(sel_idx) || isempty(searchMeanPrSpectrum_full)
    axis off;
    text(0.5, 0.5, 'No combined spectrum (FULL)', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 11, 'Color', [0.7 0.1 0.1], 'Interpreter', 'none');
else
    spec_vec = sel_w(:)' * searchMeanPrSpectrum_full(sel_idx, :);
    [pf_score, pf_peak_hz] = compute_combined_powspctrm_form_metrics( ...
        spec_vec, scan_freqs, analysis_freq_range);
    plot(scan_freqs, spec_vec, '-', 'Color', [0 0 0], 'LineWidth', 2.0);
    add_peak_point_overlay(scan_freqs, spec_vec, pf_peak_hz);
    yline(0, 'k--', 'LineWidth', 0.8);
    xlim([analysis_freq_range(1) analysis_freq_range(2)]);
    spec_finite = spec_vec(isfinite(spec_vec));
    if ~isempty(spec_finite)
        sp_min = min(spec_finite);
        sp_max = max(spec_finite);
        if isfinite(sp_min) && isfinite(sp_max) && sp_min < sp_max
            sp_range = sp_max - sp_min;
            ylim([sp_min - 0.12 * sp_range, sp_max + 0.20 * sp_range]);
        end
    end
    format_power_change_db_axis(gca);
    xlabel('Hz'); ylabel('Power [dB]');
    box on;
    text(0.02, 0.98, sprintf('PF = %.2f | peak = %.1f Hz', ...
        pf_score, pf_peak_hz), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 9, 'Interpreter', 'none', 'Color', [0.1 0.1 0.1]);
    title('Combined Spectrum (FULL)', ...
        'FontSize', 11, 'Interpreter', 'none');
    set(gca, 'FontSize', 10);
end

sgtitle(sprintf('Combined GED Components: %s (%s)', subject_id, file_tag), ...
    'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
save_figure_png(fig, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_components_combined_%s.png', subject_id, file_tag)));
close(fig);
end
