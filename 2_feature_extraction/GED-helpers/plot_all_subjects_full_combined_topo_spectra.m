function plot_all_subjects_full_combined_topo_spectra( ...
    save_dir, subjects, scan_freqs, analysis_freq_range, cfg_topo, ...
    all_topo_labels, all_topos, all_combined_spectrum_full, all_combined_eigenvalue_full)
% Group overview of combined GED topographies and spectra for the full window.
% Electrodes used in further analyses (occipital region) are marked on the topoplots.
nSubj = numel(subjects);
n_rows = 2;
n_cols = 5;
n_slots = n_rows * n_cols;
if nSubj > n_slots
    warning('GED:AllSubjectsFullOverview', ...
        'More than %d subjects; only the first %d are shown in the FULL overview figure.', ...
        n_slots, n_slots);
end
n_plot = min(nSubj, n_slots);

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
left_m = 0.04;
right_m = 0.99;
bottom_m = 0.05;
top_m = 0.97;
gap_x = 0.025;
gap_y = 0.06;
inner_gap = 0.02;
topo_frac = 0.55;
cell_w = (right_m - left_m - (n_cols - 1) * gap_x) / n_cols;
cell_h = (top_m - bottom_m - (n_rows - 1) * gap_y) / n_rows;

for subj = 1:n_plot
    row = ceil(subj / n_cols);
    col = mod(subj - 1, n_cols) + 1;
    x0 = left_m + (col - 1) * (cell_w + gap_x);
    y0 = bottom_m + (n_rows - row) * (cell_h + gap_y);
    h_spec = cell_h * (1 - topo_frac) - inner_gap / 2;
    h_topo = cell_h * topo_frac - inner_gap / 2;
    has_component = ~(isempty(all_topos{subj}) || isempty(all_topo_labels{subj}) || ...
        all(~isfinite(all_topos{subj}(:))) || isempty(all_combined_spectrum_full{subj}) || ...
        all(~isfinite(all_combined_spectrum_full{subj}(:))));
    subj_title = sprintf('Participant %s', subjects{subj});

    axes('Position', [x0, y0 + h_spec + inner_gap, cell_w, h_topo]);
    if ~has_component
        axis off;
        text(0.5, 0.5, 'NO ELIGIBLE GED COMPONENT', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.7 0.1 0.1], 'Interpreter', 'none');
    else
        topo_vec = all_topos{subj};
        topo_labels = all_topo_labels{subj};
        topo_data = [];
        topo_data.label = topo_labels;
        topo_data.avg = topo_vec(:);
        topo_data.dimord = 'chan';
        topo_vals = topo_vec(isfinite(topo_vec));
        topo_clim = max(abs(topo_vals));
        if ~isfinite(topo_clim) || topo_clim <= 0
            topo_clim = 1;
        end
        cfg_ci = cfg_topo;
        cfg_ci.zlim = [-topo_clim topo_clim];

        % Highlight occipital electrodes used in further analyses
        occ_mask = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P10|P9)', 'once')), topo_labels);
        occ_highlight = topo_labels(occ_mask);
        if ~isempty(occ_highlight)
            cfg_ci.highlight          = {'on'};
            cfg_ci.highlightchannel   = {occ_highlight};
            cfg_ci.highlightsymbol    = {'.'};
            cfg_ci.highlightsize      = {10};
            cfg_ci.highlightcolor     = {[0 0 0]};
        end

        try
            ft_topoplotER(cfg_ci, topo_data);
        catch
            imagesc(topo_vec(:)); axis tight;
            caxis([-topo_clim topo_clim]);
        end
    end
    set(gca, 'FontSize', 8);
    title(subj_title, 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'none');

    axes('Position', [x0, y0, cell_w, h_spec]);
    hold on;
    if ~has_component
        axis off;
    else
        spec_vec = all_combined_spectrum_full{subj};
        plot(scan_freqs, spec_vec, '-', 'Color', [0 0 0], 'LineWidth', 1.5);
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
        if col == 1
            ylabel('Power [dB]', 'FontSize', 8);
        end
        if row == n_rows
            xlabel('Frequency [Hz]', 'FontSize', 8);
        end
        box on;
        set(gca, 'FontSize', 8);
    end
end

save_figure_png(fig, fullfile(save_dir, 'GCP_eeg_GED_components_full_allsubjects.png'));
close(fig);
end
