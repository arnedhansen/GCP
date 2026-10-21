function plot_covariance_matrix_diagnostics(save_dir, subject_id, chan_labels, covBase, covStim, lambda)
% Diagnostic figure of baseline and stimulus covariance matrices.
if isempty(covStim) || isempty(covBase)
    return;
end
if nargin < 6 || isempty(lambda)
    lambda = 0.05;
end

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
nChans = size(covBase, 1);
if size(covBase, 2) ~= nChans || any(size(covStim) ~= [nChans nChans])
    close(fig);
    return;
end

d_stim = diag(covStim);
d_stim = d_stim(isfinite(d_stim));
stim_diag_mean = mean(d_stim);
if ~isfinite(stim_diag_mean)
    stim_diag_mean = 1;
end
d_base = diag(covBase);
d_base = d_base(isfinite(d_base));
base_diag_mean = mean(d_base);
if ~isfinite(base_diag_mean)
    base_diag_mean = 1;
end

covStim_reg = (1 - lambda) * covStim + lambda * stim_diag_mean * eye(nChans);
covBase_reg = (1 - lambda) * covBase + lambda * base_diag_mean * eye(nChans);
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

all_mats = {covStim_reg, covBase_reg, covDiff, gedOp};
nMatCols = numel(all_mats);

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
qc_lines = sprintf([...
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

col_clims = ones(1, nMatCols);
for pi = 1:nMatCols
    mat_vals = abs(all_mats{pi}(:));
    mat_vals = mat_vals(isfinite(mat_vals));
    if ~isempty(mat_vals)
        clim = prctile(mat_vals, 99);
        if ~isfinite(clim) || clim <= 0
            clim = max(mat_vals);
        end
        if isfinite(clim) && clim > 0
            col_clims(pi) = clim;
        end
    end
end

nC = 256;
nHalf = nC / 2;
tmap = linspace(0, 1, nHalf)';
low_col = [0.10 0.25 0.80];
mid_col = [1.00 1.00 1.00];
high_col = [0.80 0.15 0.15];
cmap_low = [ ...
    low_col(1) + (mid_col(1) - low_col(1)) * tmap, ...
    low_col(2) + (mid_col(2) - low_col(2)) * tmap, ...
    low_col(3) + (mid_col(3) - low_col(3)) * tmap];
cmap_high = [ ...
    mid_col(1) + (high_col(1) - mid_col(1)) * tmap, ...
    mid_col(2) + (high_col(2) - mid_col(2)) * tmap, ...
    mid_col(3) + (high_col(3) - mid_col(3)) * tmap];
colormap(fig, [cmap_low; cmap_high]);

% One compact row of square matrices plus QC text, vertically centered.
t = tiledlayout(1, 5, 'Padding', 'compact', 'TileSpacing', 'compact');
t.OuterPosition = [0.04 0.32 0.92 0.48];
panel_titles = {'Stim (reg)', 'Base (reg)', 'Stim-Base', 'R^{-1}S'};

for pi = 1:nMatCols
    nexttile;
    mat = all_mats{pi};
    if all(~isfinite(mat(:)))
        mat = zeros(nChans);
    end
    imagesc(mat);
    axis image;
    caxis([-col_clims(pi) col_clims(pi)]);
    colorbar;
    title(panel_titles{pi}, 'FontSize', 11, 'Interpreter', 'tex');
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
title('Values', 'FontSize', 11, 'Interpreter', 'none');
text(0.02, 0.98, qc_lines, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'Color', [0.05 0.05 0.05], 'Interpreter', 'none');

sgtitle(sprintf('GED Covariance Diagnostics: Subject %s', subject_id), ...
    'FontSize', 14, 'FontWeight', 'bold');
save_figure_png(fig, fullfile(save_dir, sprintf('GCP_eeg_GED_subj%s_covariance_matrix.png', subject_id)));
close(fig);
end
