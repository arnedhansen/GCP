function plot_subjects_combined_topo_spectra_grid( ...
    save_dir, subjects, subject_idx, scan_freqs, analysis_freq_range, cfg_topo, ...
    all_topo_labels, all_topos, all_combined_spectrum, file_tag)
% Grid overview of combined GED topographies and spectra.
% Rows: subjects in subject_idx. Topography left, spectrum right.
if nargin < 10 || isempty(file_tag)
    file_tag = 'batch';
end

n_rows = numel(subject_idx);
fig = figure('Position', [0 0 1512 982], 'Color', 'w');
left_m = 0.04;
right_m = 0.98;
bottom_m = 0.05;
top_m = 0.97;
gap_y = 0.035;
inner_gap = 0.012;
topo_w = 0.13;
cell_h = (top_m - bottom_m - (n_rows - 1) * gap_y) / n_rows;
spec_x = left_m + topo_w + inner_gap;
spec_w = right_m - spec_x;

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

    topo_vec = [];
    spec_vec = [];
    if numel(all_topos) >= subj
        topo_vec = all_topos{subj};
    end
    if numel(all_combined_spectrum) >= subj
        spec_vec = all_combined_spectrum{subj};
    end
    has_component = ~(isempty(topo_vec) || isempty(topo_labels) || ...
        all(~isfinite(topo_vec(:))) || isempty(spec_vec) || ...
        all(~isfinite(spec_vec(:))));

    axes('Position', [left_m, y0, topo_w, cell_h]);
    if ~has_component
        axis off;
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
    title(subj_title, 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');

    axes('Position', [spec_x, y0, spec_w, cell_h]);
    hold on;
    if ~has_component
        axis off;
        text(0.5, 0.5, 'NO ELIGIBLE GED COMPONENT', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.7 0.1 0.1], 'Interpreter', 'none');
    else
        plot(scan_freqs, spec_vec, '-', 'Color', [0 0 0], 'LineWidth', 1.6);
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
        ylabel('Power [dB]', 'FontSize', 9);
        if ri == n_rows
            xlabel('Frequency [Hz]', 'FontSize', 9);
        end
        box on;
        set(gca, 'FontSize', 9);
    end
end

outName = sprintf('GCP_eeg_GED_components_subjects_%s.png', file_tag);
save_figure_png(fig, fullfile(save_dir, outName));
close(fig);
end
