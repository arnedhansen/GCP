function plot_subjects_combined_topo_spectra_windows_grid( ...
    save_dir, subjects, subject_idx, scan_freqs, analysis_freq_range, cfg_topo, ...
    all_topo_labels, all_topos_full, all_combined_spectrum_full, file_tag)
% Grid overview of combined GED topographies and spectra (full window).
% Rows: subjects in subject_idx.
if nargin < 10 || isempty(file_tag)
    file_tag = 'batch';
end

n_rows = numel(subject_idx);
n_cols = 1;
fig = figure('Position', [0 0 1512 982], 'Color', 'w');
left_m = 0.04;
right_m = 0.99;
bottom_m = 0.04;
top_m = 0.96;
gap_x = 0.03;
gap_y = 0.035;
inner_gap = 0.012;
topo_frac = 0.55;
cell_w = (right_m - left_m - (n_cols - 1) * gap_x) / n_cols;
cell_h = (top_m - bottom_m - (n_rows - 1) * gap_y) / n_rows;

for ri = 1:n_rows
    subj = subject_idx(ri);
    if subj < 1 || subj > numel(subjects)
        continue;
    end
    topo_labels = {};
    if numel(all_topo_labels) >= subj
        topo_labels = all_topo_labels{subj};
    end
    subj_title = sprintf('Participant %s', subjects{subj});
    y0 = bottom_m + (n_rows - ri) * (cell_h + gap_y);
    x0 = left_m;
    h_spec = cell_h * (1 - topo_frac) - inner_gap / 2;
    h_topo = cell_h * topo_frac - inner_gap / 2;

    topo_vec = [];
    spec_vec = [];
    if numel(all_topos_full) >= subj
        topo_vec = all_topos_full{subj};
    end
    if numel(all_combined_spectrum_full) >= subj
        spec_vec = all_combined_spectrum_full{subj};
    end
    has_component = ~(isempty(topo_vec) || isempty(topo_labels) || ...
        all(~isfinite(topo_vec(:))) || isempty(spec_vec) || ...
        all(~isfinite(spec_vec(:))));

    axes('Position', [x0, y0 + h_spec + inner_gap, cell_w, h_topo]);
    if ~has_component
        axis off;
        text(0.5, 0.5, 'NO ELIGIBLE GED COMPONENT', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.7 0.1 0.1], 'Interpreter', 'none');
    else
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
        try
            ft_topoplotER(cfg_ci, topo_data);
        catch
            imagesc(topo_vec(:)); axis tight;
            caxis([-topo_clim topo_clim]);
        end
    end
    set(gca, 'FontSize', 7);
    title(sprintf('%s | FULL', subj_title), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none');

    axes('Position', [x0, y0, cell_w, h_spec]);
    hold on;
    if ~has_component
        axis off;
    else
        plot(scan_freqs, spec_vec, '-', 'Color', [0 0 0], 'LineWidth', 1.25);
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
        ylabel('Power [dB]', 'FontSize', 7);
        if ri == n_rows
            xlabel('Frequency [Hz]', 'FontSize', 7);
        end
        box on;
        set(gca, 'FontSize', 7);
    end
end

outName = sprintf('GCP_eeg_GED_components_windows_subjects_%s.png', file_tag);
save_figure_png(fig, fullfile(save_dir, outName));
close(fig);
end
