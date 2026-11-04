function plot_subjects_combined_topo_spectra_grid( ...
    save_dir, subjects, subject_idx, scan_freqs, analysis_freq_range, cfg_topo, ...
    all_topo_labels, all_topos, all_combined_spectrum, file_tag)
% Grid overview of combined GED topographies and spectra.
% Same layout as the all-subjects overview: topography above spectrum, 2 x 5.
% subject_idx selects which subjects fill the grid (at most 10).
if nargin < 10 || isempty(file_tag)
    file_tag = 'batch';
end

n_rows = 2;
n_cols = 5;
n_slots = n_rows * n_cols;
subject_idx = subject_idx(:)';
subject_idx = subject_idx(subject_idx >= 1 & subject_idx <= numel(subjects));
if numel(subject_idx) > n_slots
    warning('GED:SubjectGridOverview', ...
        'More than %d subjects in one grid; only the first %d are shown.', ...
        n_slots, n_slots);
end
n_plot = min(numel(subject_idx), n_slots);
if n_plot < 1
    return;
end

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

for slot = 1:n_plot
    subj = subject_idx(slot);
    row = ceil(slot / n_cols);
    col = mod(slot - 1, n_cols) + 1;
    x0 = left_m + (col - 1) * (cell_w + gap_x);
    y0 = bottom_m + (n_rows - row) * (cell_h + gap_y);
    h_spec = cell_h * (1 - topo_frac) - inner_gap / 2;
    h_topo = cell_h * topo_frac - inner_gap / 2;

    topo_labels = {};
    topo_vec = [];
    spec_vec = [];
    if numel(all_topo_labels) >= subj
        topo_labels = all_topo_labels{subj};
    end
    if numel(all_topos) >= subj
        topo_vec = all_topos{subj};
    end
    if numel(all_combined_spectrum) >= subj
        spec_vec = all_combined_spectrum{subj};
    end
    has_component = ~(isempty(topo_vec) || isempty(topo_labels) || ...
        all(~isfinite(topo_vec(:))) || isempty(spec_vec) || ...
        all(~isfinite(spec_vec(:))));
    subj_title = sprintf('Participant %s', subjects{subj});

    axes('Position', [x0, y0 + h_spec + inner_gap, cell_w, h_topo]);
    if ~has_component
        axis off;
        text(0.5, 0.5, 'NO ELIGIBLE GED COMPONENT', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.7 0.1 0.1], 'Interpreter', 'none');
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
    set(gca, 'FontSize', 8);
    title(subj_title, 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'none');

    axes('Position', [x0, y0, cell_w, h_spec]);
    hold on;
    if ~has_component
        axis off;
    else
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

outName = sprintf('GCP_eeg_GED_components_subjects_%s.png', file_tag);
save_figure_png(fig, fullfile(save_dir, outName));
close(fig);
end
