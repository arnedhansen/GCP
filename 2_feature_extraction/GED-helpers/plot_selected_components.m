function plot_selected_components(save_dir, subject_id, win_name, scan_freqs, searchTopos, ...
    searchMeanPrSpectrum, eigval_vec, emg_class, ...
    eligible, ...
    post_front_vec, post_temp_vec, topo_peak_frac_vec, emg_hf_slope_vec, ...
    cfg_topo, topo_labels, ...
    powspctrm_form_score, ...
    selected_idx)
% Selected vs rejected GED components with topography, spectrum, and gate metrics.
nComp = numel(eigval_vec);
if nComp < 1
    return;
end
selected_mask = false(nComp, 1);
if ~isempty(selected_idx)
    selected_idx = unique(selected_idx(:));
    selected_idx = selected_idx(isfinite(selected_idx));
    selected_idx = selected_idx(selected_idx >= 1 & selected_idx <= nComp);
    selected_idx = round(selected_idx);
    selected_mask(selected_idx) = true;
end
final_eligible_mask = logical(eligible(:));
final_rejected_mask = ~final_eligible_mask;
sel_idx = find(selected_mask);
rej_idx = find(final_rejected_mask & ~selected_mask);
[~, so] = sort(eigval_vec(sel_idx), 'descend');
sel_idx = sel_idx(so);
[~, ro] = sort(eigval_vec(rej_idx), 'descend');
rej_idx = rej_idx(ro);

% Top selected/rejected components in one figure (max 5 per row):
% row 1 = selected topographies, row 2 = selected spectra,
% row 3 = rejected topographies, row 4 = rejected spectra.
nCols = 5;
nShowSel = min(nCols, numel(sel_idx));
nShowRej = min(nCols, numel(rej_idx));
figSel = figure('Position', [0 0 1512 982], 'Color', 'w');
for k = 1:nCols
    % Row 1: selected topographies
    subplot(4, nCols, k);
    if k <= nShowSel
        ci = sel_idx(k);
        topo_data = [];
        topo_data.label = topo_labels;
        topo_data.avg = searchTopos(:, ci);
        topo_data.dimord = 'chan';
        cfg_ci = cfg_topo;
        topo_vals = topo_data.avg(isfinite(topo_data.avg));
        topo_clim_ci = max(abs(topo_vals));
        if ~isfinite(topo_clim_ci) || topo_clim_ci <= 0
            topo_clim_ci = 1;
        end
        cfg_ci.zlim = [-topo_clim_ci topo_clim_ci];
        try
            ft_topoplotER(cfg_ci, topo_data);
        catch
            imagesc(topo_data.avg(:)); axis tight;
            caxis([-topo_clim_ci topo_clim_ci]);
        end
        ttl = sprintf('C%d (%s)', ci, emg_class{ci});
        title(ttl, 'FontSize', 8, 'Interpreter', 'none');
    else
        axis off;
    end

    % Row 2: selected spectra
    subplot(4, nCols, nCols + k); hold on;
    if k <= nShowSel
        ci = sel_idx(k);
        spec_data = searchMeanPrSpectrum(ci, :);
        plot(scan_freqs, spec_data, '-', 'Color', [0 0 0], 'LineWidth', 1.4);
        yline(0, 'k--');
        spec_min = min(spec_data(isfinite(spec_data)));
        spec_max = max(spec_data(isfinite(spec_data)));
        if isfinite(spec_min) && isfinite(spec_max) && spec_min ~= spec_max
            spec_range = spec_max - spec_min;
            ylim([spec_min - 0.10 * spec_range, spec_max + 0.10 * spec_range]);
        end
        [info_lines, info_viol] = build_selection_gate_info_columns( ...
            ci, post_front_vec, post_temp_vec, topo_peak_frac_vec, emg_hf_slope_vec);
        plot_rejection_info_text_columns(info_lines, info_viol);
        format_power_change_db_axis(gca);
        xlabel('Hz'); ylabel('Power [dB]');
        title(sprintf('\\lambda=%.2f, PF=%.2f', ...
            eigval_vec(ci), powspctrm_form_score(ci)), 'FontSize', 7);
        box on;
    else
        axis off;
    end

    % Row 3: rejected topographies
    subplot(4, nCols, 2 * nCols + k);
    if k <= nShowRej
        ci = rej_idx(k);
        topo_data = [];
        topo_data.label = topo_labels;
        topo_data.avg = searchTopos(:, ci);
        topo_data.dimord = 'chan';
        cfg_ci = cfg_topo;
        topo_vals = topo_data.avg(isfinite(topo_data.avg));
        topo_clim_ci = max(abs(topo_vals));
        if ~isfinite(topo_clim_ci) || topo_clim_ci <= 0
            topo_clim_ci = 1;
        end
        cfg_ci.zlim = [-topo_clim_ci topo_clim_ci];
        try
            ft_topoplotER(cfg_ci, topo_data);
        catch
            imagesc(topo_data.avg(:)); axis tight;
            caxis([-topo_clim_ci topo_clim_ci]);
        end
        ttl = sprintf('C%d (%s)', ci, emg_class{ci});
        title(ttl, 'FontSize', 8, 'Interpreter', 'none');
    else
        axis off;
    end

    % Row 4: rejected spectra
    subplot(4, nCols, 3 * nCols + k); hold on;
    if k <= nShowRej
        ci = rej_idx(k);
        spec_data = searchMeanPrSpectrum(ci, :);
        plot(scan_freqs, spec_data, '-', 'Color', [0 0 0], 'LineWidth', 1.4);
        yline(0, 'k--');
        spec_min = min(spec_data(isfinite(spec_data)));
        spec_max = max(spec_data(isfinite(spec_data)));
        if isfinite(spec_min) && isfinite(spec_max) && spec_min ~= spec_max
            spec_range = spec_max - spec_min;
            ylim([spec_min - 0.10 * spec_range, spec_max + 0.10 * spec_range]);
        end
        [info_lines, info_viol] = build_selection_gate_info_columns( ...
            ci, post_front_vec, post_temp_vec, topo_peak_frac_vec, emg_hf_slope_vec);
        plot_rejection_info_text_columns(info_lines, info_viol);
        format_power_change_db_axis(gca);
        xlabel('Hz'); ylabel('Power [dB]');
        title(sprintf('\\lambda=%.2f, PF=%.2f', ...
            eigval_vec(ci), powspctrm_form_score(ci)), 'FontSize', 7);
        box on;
    else
        axis off;
    end
end
if nShowSel == 0
    annotation(figSel, 'textbox', [0.25 0.77 0.50 0.06], ...
        'String', 'NO ELIGIBLE COMPONENTS', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 16, 'Color', [0.85 0.05 0.05], ...
        'Interpreter', 'none');
end
annotation(figSel, 'textbox', [0.010 0.56 0.040 0.34], ...
    'String', 'Selected Components', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontWeight', 'bold', 'FontSize', 14, 'Rotation', 90, 'FitBoxToText', 'off');
annotation(figSel, 'textbox', [0.010 0.09 0.040 0.34], ...
    'String', 'Rejected Components', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontWeight', 'bold', 'FontSize', 14, 'Rotation', 90, 'FitBoxToText', 'off');
sgtitle(sprintf('Selected and Rejected Components: %s (%s)', subject_id, win_name), ...
    'FontSize', 14, 'FontWeight', 'bold');
drawnow;
pause(0.05);
save_figure_png(figSel, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_component_selection_%s.png', subject_id, win_name)));
close(figSel);

end
