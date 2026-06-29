%% GCP GED Trial-Level Boxplots
%
% Pooled trial-by-trial boxplots for gamma peak frequency and peak power
% after GED (full analysis window). One figure per metric.
%
% Data source: GCP_eeg_GED.mat
%   - trials_peaks{cond, subj}: per-trial peak frequency [Hz]
%   - trials_powratio_fullscan{cond, subj}: per-trial power spectra [dB]
%   - Peak power is reconstructed as mean power within peak frequency +/- 5 Hz
%     (same definition as GCP_eeg_fex_GED.m).

clear; close all;
startup
[subjects, paths, colors, ~] = setup('GCP', 0);

ged_path = fullfile(paths.features, 'GCP_eeg_GED.mat');
if ~isfile(ged_path)
    error('GED results not found: %s', ged_path);
end

ged = load(ged_path);
condLabels = ged.condLabels;
if isfield(ged, 'subjects')
    subjects = ged.subjects;
end
nCond = numel(condLabels);
nSubj = size(ged.trials_peaks, 2);
peak_power_halfwidth_hz = 5;

fig_dir = fullfile(paths.figures, 'eeg', 'ged');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

if isfield(ged, 'trials_outlier_mask_power_full')
    power_outlier_masks = ged.trials_outlier_mask_power_full;
else
    power_outlier_masks = cell(nCond, nSubj);
end

[y_freq, g_freq, s_freq, y_power, g_power, s_power] = pool_ged_trial_metrics( ...
    ged.trials_peaks, ged.trials_powratio_fullscan, ged.scan_freqs, ...
    power_outlier_masks, peak_power_halfwidth_hz, nCond, nSubj, subjects);

fprintf('Trial-level GED data pooled: freq n=%d, power n=%d\n', ...
    numel(y_freq), numel(y_power));

run_ged_trial_glmm(y_freq, g_freq, s_freq, condLabels, 'GammaFrequency');
run_ged_trial_glmm(y_power, g_power, s_power, condLabels, 'GammaPower');

plot_ged_trial_boxplot(y_freq, g_freq, condLabels, colors, nCond, ...
    'Peak Gamma Frequency [Hz]', 'GED Trial Gamma Frequency', [30 90]);
exportgraphics(gcf, fullfile(fig_dir, 'GCP_eeg_GED_boxplot_trial_gamma_freq.png'), ...
    'Resolution', 600);

plot_ged_trial_boxplot(y_power, g_power, condLabels, colors, nCond, ...
    'Peak Gamma Power [dB]', 'GED Trial Gamma Power', []);
exportgraphics(gcf, fullfile(fig_dir, 'GCP_eeg_GED_boxplot_trial_gamma_power.png'), ...
    'Resolution', 600);

fprintf('Saved figures to %s\n', fig_dir);

function [y_freq, g_freq, s_freq, y_power, g_power, s_power] = pool_ged_trial_metrics( ...
    trials_peaks, trials_powratio_fullscan, scan_freqs, ...
    trials_outlier_mask_power_full, peak_power_halfwidth_hz, nCond, nSubj, subjects)

y_freq = [];
g_freq = [];
s_freq = {};
y_power = [];
g_power = [];
s_power = {};

for c = 1:nCond
    for s = 1:nSubj
        pf = trials_peaks{c, s};
        pr = trials_powratio_fullscan{c, s};
        if isempty(pf) || isempty(pr)
            continue
        end

        pf = pf(:);
        nTrl = numel(pf);
        if size(pr, 1) ~= nTrl
            warning('GCP_eeg_GED_trial_boxplots:SizeMismatch', ...
                'Condition %d, subject %d: %d peak freqs but %d spectra; skipping.', ...
                c, s, nTrl, size(pr, 1));
            continue
        end

        pp = reconstruct_trial_peak_power(pf, pr, scan_freqs, peak_power_halfwidth_hz);

        if ~isempty(trials_outlier_mask_power_full) && ...
                ~isempty(trials_outlier_mask_power_full{c, s})
            mask = trials_outlier_mask_power_full{c, s}(:);
            if numel(mask) == nTrl
                pp(mask) = NaN;
            end
        end

        valid_f = isfinite(pf);
        valid_p = isfinite(pp);
        subj_id = subjects{s};

        y_freq = [y_freq; pf(valid_f)]; %#ok<AGROW>
        g_freq = [g_freq; repmat(c, sum(valid_f), 1)]; %#ok<AGROW>
        s_freq = [s_freq; repmat({subj_id}, sum(valid_f), 1)]; %#ok<AGROW>

        y_power = [y_power; pp(valid_p)]; %#ok<AGROW>
        g_power = [g_power; repmat(c, sum(valid_p), 1)]; %#ok<AGROW>
        s_power = [s_power; repmat({subj_id}, sum(valid_p), 1)]; %#ok<AGROW>
    end
end
end

function pp = reconstruct_trial_peak_power(pf, pr, scan_freqs, halfwidth_hz)
nTrl = numel(pf);
pp = nan(nTrl, 1);
scan_freqs = scan_freqs(:)';

for t = 1:nTrl
    if ~isfinite(pf(t))
        continue
    end
    band = abs(scan_freqs - pf(t)) <= halfwidth_hz;
    if ~any(band)
        continue
    end
    pp(t) = mean(pr(t, band), 'omitnan');
end
end

function plot_ged_trial_boxplot(y, g, condLabels, colors, nCond, yLabel, figTitle, yLimFixed)
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

dens_width = 0.35;
for c = 1:nCond
    vals = y(g == c);
    vals = vals(isfinite(vals));
    if numel(vals) < 3
        continue
    end
    [f_dens, xi] = ksdensity(vals);
    f_dens = f_dens / max(f_dens) * dens_width;
    patch(c - f_dens - 0.05, xi, colors(c, :), ...
        'FaceAlpha', 0.3, 'EdgeColor', colors(c, :), 'LineWidth', 1, ...
        'HandleVisibility', 'off');
end

if ~isempty(y)
    boxplot(y, g, 'Colors', 'k', 'Symbol', '', 'Widths', 0.35);
end

for c = 1:nCond
    vals = y(g == c);
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue
    end
    xJit = c + (rand(numel(vals), 1) - 0.5) * 0.12;
    scatter(xJit, vals, 18, colors(c, :), 'filled', ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
end

xlim([0.5 nCond + 0.5]);
if ~isempty(yLimFixed)
    ylim(yLimFixed);
elseif ~isempty(y)
    yMax = max(y, [], 'omitnan');
    if isempty(yMax) || ~isfinite(yMax) || yMax == 0
        ylim([0 1]);
    else
        ylim([0 yMax * 1.1]);
    end
end

set(gca, 'XTick', 1:nCond, 'XTickLabel', condLabels, ...
    'FontSize', 20, 'Box', 'off');
ylabel(yLabel);
title(figTitle, 'FontSize', 30, 'FontWeight', 'bold');
hold off;
end

function run_ged_trial_glmm(y, g, s, condLabels, outcomeName)
fprintf('\n============================================================\n');
fprintf('Trial-Level GLMM: %s ~ Condition + (1|SubjectID)\n', outcomeName);
fprintf('============================================================\n');

valid = isfinite(y) & g >= 1 & g <= numel(condLabels);
y = y(valid);
g = g(valid);
s = s(valid);

if isempty(y)
    fprintf('GLMM skipped: no valid trial-level %s values were found.\n', outcomeName);
    return
end

cond_cat = categorical(g, 1:numel(condLabels), condLabels, 'Ordinal', true);
tbl = table(y, cond_cat, categorical(s), ...
    'VariableNames', {outcomeName, 'Condition', 'SubjectID'});

fprintf('Observations: %d\n', height(tbl));
fprintf('Subjects: %d\n', numel(unique(tbl.SubjectID)));
fprintf('Condition counts:\n');
cond_levels = categories(tbl.Condition);
cond_counts = zeros(numel(cond_levels), 1);
for ci = 1:numel(cond_levels)
    cond_counts(ci) = sum(tbl.Condition == cond_levels{ci});
end
disp(table(categorical(cond_levels, cond_levels, 'Ordinal', true), cond_counts, ...
    'VariableNames', {'Condition', 'nTrials'}));

formula = sprintf('%s ~ Condition + (1|SubjectID)', outcomeName);
try
    mdl = fitglme(tbl, formula, ...
        'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace');

    fprintf('\nModel fit summary:\n');
    disp(mdl);
    fprintf('\nFixed-effects coefficients:\n');
    disp(mdl.Coefficients);

    coef_tbl = mdl.Coefficients;
    coef_names = cellstr(string(coef_tbl.Name));
    is_cond_coef = startsWith(coef_tbl.Name, 'Condition_');
    if any(is_cond_coef)
        cond_names = coef_names(is_cond_coef);
        p_raw = coef_tbl.pValue(is_cond_coef);
        p_fdr = bh_fdr_adjust(p_raw);
        sig_fdr = p_fdr < 0.05;

        fprintf('\nFDR-corrected p-values (Benjamini-Hochberg, q=0.05):\n');
        disp(table(cond_names, p_raw, p_fdr, sig_fdr, ...
            'VariableNames', {'Contrast', 'pRaw', 'pFDR', 'isSignificantFDR'}));
    else
        fprintf('\nNo condition contrasts found for FDR correction.\n');
    end

    report_pairwise_ttests(tbl, outcomeName, condLabels);
catch ME
    fprintf('GLMM fit failed: %s\n', ME.message);
    report_pairwise_ttests(tbl, outcomeName, condLabels);
end
end

function report_pairwise_ttests(tbl, outcomeName, condLabels)
fprintf('\nPost-hoc pairwise paired t-tests (subject means, Benjamini-Hochberg FDR):\n');

cond_levels = categories(tbl.Condition);
nCond = numel(cond_levels);
subj_ids = categories(tbl.SubjectID);
nSubj = numel(subj_ids);
subj_means = nan(nSubj, nCond);

for si = 1:nSubj
    for ci = 1:nCond
        idx = tbl.SubjectID == subj_ids(si) & tbl.Condition == cond_levels{ci};
        vals = tbl.(outcomeName)(idx);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            subj_means(si, ci) = mean(vals);
        end
    end
end

pair_idx = nchoosek(1:nCond, 2);
nPairs = size(pair_idx, 1);
comparison = cell(nPairs, 1);
n_used = zeros(nPairs, 1);
mean_diff = nan(nPairs, 1);
t_stat = nan(nPairs, 1);
df_val = nan(nPairs, 1);
p_raw = nan(nPairs, 1);

for pi = 1:nPairs
    c1 = pair_idx(pi, 1);
    c2 = pair_idx(pi, 2);
    x1 = subj_means(:, c1);
    x2 = subj_means(:, c2);
    valid = isfinite(x1) & isfinite(x2);
    n_used(pi) = sum(valid);
    comparison{pi} = sprintf('%s vs %s', cond_levels{c1}, cond_levels{c2});

    if n_used(pi) < 2
        continue
    end

    [~, p_raw(pi), ~, stats] = ttest(x1(valid), x2(valid));
    mean_diff(pi) = mean(x1(valid) - x2(valid));
    t_stat(pi) = stats.tstat;
    df_val(pi) = stats.df;
end

p_fdr = bh_fdr_adjust(p_raw);
sig_fdr = p_fdr < 0.05;

result_tbl = table(comparison, n_used, mean_diff, t_stat, df_val, p_raw, p_fdr, sig_fdr, ...
    'VariableNames', {'Comparison', 'nSubjects', 'MeanDiff', 't', 'df', 'pRaw', 'pFDR', 'isSignificantFDR'});
disp(result_tbl);
end

function p_adj = bh_fdr_adjust(p_raw)
p = p_raw(:);
n = numel(p);
p_adj = nan(size(p));
if n == 0
    return
end

[p_sorted, sort_idx] = sort(p);
rank_idx = (1:n)';
q_sorted = p_sorted .* n ./ rank_idx;
q_sorted = min(q_sorted, 1);

for k = n - 1:-1:1
    q_sorted(k) = min(q_sorted(k), q_sorted(k + 1));
end

p_adj(sort_idx) = q_sorted;
p_adj = reshape(p_adj, size(p_raw));
end
