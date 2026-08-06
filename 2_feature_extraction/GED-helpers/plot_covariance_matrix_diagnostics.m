function plot_covariance_matrix_diagnostics(save_dir, subject_id, chan_labels, covBase_full, covStim_per_win, win_names_cap, lambdas)
% Diagnostic figure of baseline and stimulus covariance matrices per window.
if isempty(covStim_per_win) || isempty(covBase_full)
    return;
end
nWins = min(3, numel(covStim_per_win));
if nWins < 1
    return;
end
if isempty(win_names_cap)
    win_names_cap = {'Full', 'Early', 'Late'};
end
if isempty(lambdas)
    lambdas = repmat(0.05, 1, nWins);
end

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
nMatCols = 4;
nPlotCols = 5;
tiledlayout(nWins, nPlotCols, 'Padding', 'compact', 'TileSpacing', 'compact');

nChans = size(covBase_full, 1);
if size(covBase_full, 2) ~= nChans
    close(fig);
    return;
end

% Build matrices first so column-wise shared color limits are used.
all_mats = cell(nWins, nMatCols);
qc_lines_by_win = cell(nWins, 1);
for wi = 1:nWins
    covStim_w = covStim_per_win{wi};
    if isempty(covStim_w) || any(size(covStim_w) ~= [nChans nChans])
        covStim_w = nan(nChans);
    end
    lam_w = lambdas(min(wi, numel(lambdas)));

    d_stim = diag(covStim_w);
    d_stim = d_stim(isfinite(d_stim));
    stim_diag_mean = mean(d_stim);
    if ~isfinite(stim_diag_mean)
        stim_diag_mean = 1;
    end
    d_base = diag(covBase_full);
    d_base = d_base(isfinite(d_base));
    base_diag_mean = mean(d_base);
    if ~isfinite(base_diag_mean)
        base_diag_mean = 1;
    end

    covStim_reg = (1 - lam_w) * covStim_w + lam_w * stim_diag_mean * eye(nChans);
    covBase_reg = (1 - lam_w) * covBase_full + lam_w * base_diag_mean * eye(nChans);
    covStim_reg = real(0.5 * (covStim_reg + covStim_reg'));
    covBase_reg = real(0.5 * (covBase_reg + covBase_reg'));
    covDiff = covStim_reg - covBase_reg;

    rc_base = rcond(covBase_reg);
    if ~isfinite(rc_base) || rc_base < 1e-12
        gedOp = pinv(covBase_reg) * covStim_reg;
    else
        gedOp = covBase_reg \ covStim_reg;
    end
    gedOp = real(gedOp);

    all_mats{wi, 1} = covStim_reg;
    all_mats{wi, 2} = covBase_reg;
    all_mats{wi, 3} = covDiff;
    all_mats{wi, 4} = gedOp;

    sym_stim = norm(covStim_reg - covStim_reg', 'fro') / max(norm(covStim_reg, 'fro'), eps);
    sym_base = norm(covBase_reg - covBase_reg', 'fro') / max(norm(covBase_reg, 'fro'), eps);
    eig_base = eig(covBase_reg);
    eig_base = real(eig_base(isfinite(eig_base)));
    if isempty(eig_base)
        min_eig_base = NaN;
    else
        min_eig_base = min(eig_base);
    end
    cond_base = cond(covBase_reg);
    try
        eigvals_ged = eig(covStim_reg, covBase_reg);
    catch
        eigvals_ged = eig(gedOp);
    end
    eigvals_ged = real(eigvals_ged(isfinite(eigvals_ged)));
    if isempty(eigvals_ged)
        top_ged = nan(1, 5);
    else
        eigvals_ged = sort(eigvals_ged, 'descend');
        top_ged = nan(1, 5);
        nTop = min(5, numel(eigvals_ged));
        top_ged(1:nTop) = eigvals_ged(1:nTop);
    end
    qc_lines_by_win{wi} = sprintf([...
        'sym(S): %.2e\n', ...
        'sym(R): %.2e\n', ...
        'minEig(R): %.3g\n', ...
        'cond(R): %.2e\n', ...
        'GED top5:\n', ...
        '1) %.3f\n', ...
        '2) %.3f\n', ...
        '3) %.3f\n', ...
        '4) %.3f\n', ...
        '5) %.3f'], ...
        sym_stim, sym_base, min_eig_base, cond_base, ...
        top_ged(1), top_ged(2), top_ged(3), top_ged(4), top_ged(5));
end

col_clims = ones(1, nMatCols);
for pi = 1:nMatCols
    col_vals = [];
    for wi = 1:nWins
        mat_vals = abs(all_mats{wi, pi}(:));
        mat_vals = mat_vals(isfinite(mat_vals));
        col_vals = [col_vals; mat_vals];
    end
    if ~isempty(col_vals)
        clim = prctile(col_vals, 99);
        if ~isfinite(clim) || clim <= 0
            clim = max(col_vals);
        end
        if isfinite(clim) && clim > 0
            col_clims(pi) = clim;
        end
    end
end

% Blue-white-red diverging map for signed matrix interpretation.
nC = 256;
nHalf = nC / 2;
t = linspace(0, 1, nHalf)';
low_col = [0.10 0.25 0.80];
mid_col = [1.00 1.00 1.00];
high_col = [0.80 0.15 0.15];
cmap_low = [ ...
    low_col(1) + (mid_col(1) - low_col(1)) * t, ...
    low_col(2) + (mid_col(2) - low_col(2)) * t, ...
    low_col(3) + (mid_col(3) - low_col(3)) * t];
cmap_high = [ ...
    mid_col(1) + (high_col(1) - mid_col(1)) * t, ...
    mid_col(2) + (high_col(2) - mid_col(2)) * t, ...
    mid_col(3) + (high_col(3) - mid_col(3)) * t];
colormap(fig, [cmap_low; cmap_high]);

for wi = 1:nWins
    mats = all_mats(wi, :);
    panel_titles = { ...
        sprintf('%s Stim (reg)', win_names_cap{wi}), ...
        sprintf('%s Base (reg)', win_names_cap{wi}), ...
        sprintf('%s Stim-Base', win_names_cap{wi}), ...
        sprintf('%s R^{-1}S', win_names_cap{wi})};

    for pi = 1:nMatCols
        nexttile;
        mat = mats{pi};
        if all(~isfinite(mat(:)))
            mat = zeros(nChans);
        end
        imagesc(mat);
        axis image;
        caxis([-col_clims(pi) col_clims(pi)]);
        colorbar;
        title(panel_titles{pi}, 'FontSize', 10, 'Interpreter', 'none');
        set(gca, 'FontSize', 8, 'YDir', 'normal');
        if nChans <= 80 && ~isempty(chan_labels) && numel(chan_labels) == nChans
            xticks(1:nChans); yticks(1:nChans);
            xticklabels(chan_labels); yticklabels(chan_labels);
            xtickangle(90);
        else
            xlabel('Channels');
            ylabel('Channels');
        end
    end

    nexttile;
    axis off;
    title(sprintf('%s Window Values', win_names_cap{wi}), 'FontSize', 10, 'Interpreter', 'none');
    text(0.01, 0.99, qc_lines_by_win{wi}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 8, 'Color', [0.05 0.05 0.05], 'Interpreter', 'none');
end

sgtitle(sprintf('GED Covariance Diagnostics: Subject %s', subject_id), ...
    'FontSize', 14, 'FontWeight', 'bold');
save_figure_png(fig, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_covariance_matrix.png', subject_id)));
close(fig);
end
