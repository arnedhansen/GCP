%% GCP neural vs artifactual gamma checks
%
% Uses the trial-level merged table after GCP_eeg_fex_GED_MSfree.m:
%   1) H5/H6 on MS-free peak frequency / power
%   2) Contrast effects with microsaccade rate as covariate
%   3) Survival summary: standard vs MS-free contrast slopes
%
% Input:
%   data/features/GCP_merged_data_trials.csv
%
% Output:
%   figures/stats/neural_vs_artifact/
%   data/features/GCP_stats_neural_vs_artifact.mat

%% Setup
startup
[subjects, paths, colors] = setup('GCP', 0);

merged_csv = fullfile(paths.features, 'GCP_merged_data_trials.csv');
if ~isfile(merged_csv)
    error('GCP_stats_neural_vs_artifact:NoMerge', ...
        'Missing %s. Run GCP_master_matrix_trials.m after GED MSfree.', merged_csv);
end
T = readtable(merged_csv);

if ismember('Include', T.Properties.VariableNames)
    T = T(T.Include == 1 | isnan(T.Include), :);
end

id_var = pick_var(T, {'ID'});
cond_var = pick_var(T, {'Condition'});
freq_var = pick_var(T, {'GED_GammaFrequency', 'GammaFrequency'});
pow_var = pick_var(T, {'GED_GammaPower', 'GammaPower'});
freq_ms_var = pick_var(T, {'GED_GammaFrequency_MSfree', 'GammaFrequency_MSfree'});
pow_ms_var = pick_var(T, {'GED_GammaPower_MSfree', 'GammaPower_MSfree'});
pct_var = pick_var(T, {'GED_MSfree_pctKept', 'MSfree_pctKept'});
ms_var = pick_var(T, {'Gaze_MSRate_bl', 'Gaze_MSRate', 'MSRate_bl', 'MSRate'});

need = {id_var, cond_var, freq_var, pow_var};
if any(cellfun(@isempty, need))
    error('GCP_stats_neural_vs_artifact:MissingCols', ...
        'Merged table missing required GED/ID/Condition columns.');
end
if isempty(freq_ms_var) || isempty(pow_ms_var)
    error('GCP_stats_neural_vs_artifact:NoMSfree', ...
        ['No GED_*_MSfree columns. Run GCP_eeg_fex_GED_MSfree.m and ', ...
         'GCP_master_matrix_trials.m first.']);
end

contrast_pct = [25 50 75 100];
condLabels = {'25%', '50%', '75%', '100%'};
fig_dir = fullfile(paths.figures, 'stats', 'neural_vs_artifact');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Contrast code -> percent
T.Contrast = nan(height(T), 1);
for c = 1:4
    T.Contrast(T.(cond_var) == c) = contrast_pct(c);
end

%% Subject-condition means
uids = unique(T.(id_var)(isfinite(T.(id_var))));
nSubj = numel(uids);
mean_freq = nan(4, nSubj);
mean_pow = nan(4, nSubj);
mean_freq_ms = nan(4, nSubj);
mean_pow_ms = nan(4, nSubj);
mean_msrate = nan(4, nSubj);
mean_pct = nan(4, nSubj);

for si = 1:nSubj
    for c = 1:4
        idx = T.(id_var) == uids(si) & T.(cond_var) == c;
        mean_freq(c, si) = mean(T.(freq_var)(idx), 'omitnan');
        mean_pow(c, si) = mean(T.(pow_var)(idx), 'omitnan');
        mean_freq_ms(c, si) = mean(T.(freq_ms_var)(idx), 'omitnan');
        mean_pow_ms(c, si) = mean(T.(pow_ms_var)(idx), 'omitnan');
        if ~isempty(ms_var)
            mean_msrate(c, si) = mean(T.(ms_var)(idx), 'omitnan');
        end
        if ~isempty(pct_var)
            mean_pct(c, si) = mean(T.(pct_var)(idx), 'omitnan');
        end
    end
end

slope_freq = contrast_slope(mean_freq, contrast_pct);
slope_freq_ms = contrast_slope(mean_freq_ms, contrast_pct);
slope_pow = contrast_slope(mean_pow, contrast_pct);
slope_pow_ms = contrast_slope(mean_pow_ms, contrast_pct);

%% Trial-level LMEs
lme_freq = fit_contrast_ms_lme(T, id_var, freq_var, ms_var, 'PeakFreq');
lme_freq_msfree = fit_contrast_ms_lme(T, id_var, freq_ms_var, ms_var, 'PeakFreq_MSfree');
lme_pow = fit_contrast_ms_lme(T, id_var, pow_var, ms_var, 'PeakPower');
lme_pow_msfree = fit_contrast_ms_lme(T, id_var, pow_ms_var, ms_var, 'PeakPower_MSfree');

%% Figure: H5/H6 standard vs MS-free
fig = figure('Position', [0 0 1512 982]);
set(fig, 'Color', 'w');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ms_guard_ms = 100;
msfree_meta = fullfile(paths.features, 'GCP_eeg_GED_MSfree.mat');
if isfile(msfree_meta)
    tmp = load(msfree_meta, 'ms_guard_s');
    if isfield(tmp, 'ms_guard_s') && isfinite(tmp.ms_guard_s)
        ms_guard_ms = 1000 * tmp.ms_guard_s;
    end
end

nexttile; plot_crf(mean_freq, contrast_pct, colors, nSubj);
ylabel('Peak frequency [Hz]'); title('H5 standard GED');

nexttile; plot_crf(mean_freq_ms, contrast_pct, colors, nSubj);
ylabel('Peak frequency [Hz]'); title(sprintf('H5 MS-free (+/-%.0f ms)', ms_guard_ms));

nexttile; plot_crf(mean_pow, contrast_pct, colors, nSubj);
ylabel('Peak power [dB]'); title('H6 standard GED');

nexttile; plot_crf(mean_pow_ms, contrast_pct, colors, nSubj);
ylabel('Peak power [dB]'); title('H6 MS-free');

exportgraphics(fig, fullfile(fig_dir, 'GCP_H5H6_standard_vs_MSfree.png'), 'Resolution', 300);

%% Figure: slope survival
fig2 = figure('Position', [0 0 1512 982]);
set(fig2, 'Color', 'w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold on;
scatter(slope_freq, slope_freq_ms, 40, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.7);
lims = [min([slope_freq; slope_freq_ms], [], 'omitnan'), max([slope_freq; slope_freq_ms], [], 'omitnan')];
if all(isfinite(lims)) && lims(1) < lims(2)
    plot(lims, lims, 'k:');
end
xlabel('Contrast slope (standard)');
ylabel('Contrast slope (MS-free)');
title('H5 frequency slope survival');
hold off;

nexttile;
hold on;
scatter(slope_pow, slope_pow_ms, 40, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.7);
lims = [min([slope_pow; slope_pow_ms], [], 'omitnan'), max([slope_pow; slope_pow_ms], [], 'omitnan')];
if all(isfinite(lims)) && lims(1) < lims(2)
    plot(lims, lims, 'k:');
end
xlabel('Contrast slope (standard)');
ylabel('Contrast slope (MS-free)');
title('H6 power slope survival');
hold off;

exportgraphics(fig2, fullfile(fig_dir, 'GCP_contrast_slope_survival_MSfree.png'), 'Resolution', 300);

%% Console summary
fprintf('\n=== Neural vs artifactual summary ===\n');
fprintf('Mean %% samples kept (MS-free): %.1f%%\n', 100 * mean(mean_pct(:), 'omitnan'));
print_lme('PeakFreq ~ Contrast + MSRate', lme_freq);
print_lme('PeakFreq_MSfree ~ Contrast + MSRate', lme_freq_msfree);
print_lme('PeakPower ~ Contrast + MSRate', lme_pow);
print_lme('PeakPower_MSfree ~ Contrast + MSRate', lme_pow_msfree);
fprintf('Mean H5 slope standard: %.4f | MS-free: %.4f\n', ...
    mean(slope_freq, 'omitnan'), mean(slope_freq_ms, 'omitnan'));
fprintf('Mean H6 slope standard: %.4f | MS-free: %.4f\n', ...
    mean(slope_pow, 'omitnan'), mean(slope_pow_ms, 'omitnan'));

out_mat = fullfile(paths.features, 'GCP_stats_neural_vs_artifact.mat');
save(out_mat, 'uids', 'mean_freq', 'mean_pow', 'mean_freq_ms', 'mean_pow_ms', ...
    'mean_msrate', 'mean_pct', 'slope_freq', 'slope_freq_ms', 'slope_pow', 'slope_pow_ms', ...
    'lme_freq', 'lme_freq_msfree', 'lme_pow', 'lme_pow_msfree', ...
    'contrast_pct', 'condLabels', 'subjects');
fprintf('Saved %s\n', out_mat);

%% Helpers
function name = pick_var(T, candidates)
name = '';
for i = 1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        name = candidates{i};
        return
    end
end
end

function slopes = contrast_slope(mat, contrast_pct)
% mat: conditions x subjects
nSubj = size(mat, 2);
slopes = nan(nSubj, 1);
x = contrast_pct(:);
for si = 1:nSubj
    y = mat(:, si);
    ok = isfinite(x) & isfinite(y);
    if sum(ok) >= 3
        p = polyfit(x(ok), y(ok), 1);
        slopes(si) = p(1);
    end
end
end

function plot_crf(mat, contrast_pct, colors, nSubj)
hold on;
mu = mean(mat, 2, 'omitnan');
se = std(mat, 0, 2, 'omitnan') ./ sqrt(sum(isfinite(mat), 2));
for c = 1:4
    jitter = 0.04 * randn(1, nSubj);
    scatter(contrast_pct(c) + 8 * jitter, mat(c, :), 18, colors(c, :), ...
        'filled', 'MarkerFaceAlpha', 0.35, 'HandleVisibility', 'off');
end
errorbar(contrast_pct, mu, se, '-o', 'Color', [0.1 0.1 0.1], ...
    'LineWidth', 1.5, 'MarkerFaceColor', [0.1 0.1 0.1]);
xlim([15 110]);
set(gca, 'XTick', contrast_pct, 'XTickLabel', {'25', '50', '75', '100'});
xlabel('Contrast [%]');
hold off;
end

function S = fit_contrast_ms_lme(T, id_var, y_var, ms_var, label)
S = struct('label', label, 'ok', false, 'coef', [], 'p', [], 'formula', '');
y = T.(y_var);
ok = isfinite(y) & isfinite(T.Contrast) & isfinite(T.(id_var));
if ~isempty(ms_var)
    ok = ok & isfinite(T.(ms_var));
end
if sum(ok) < 50
    return
end
tbl = table(y(ok), T.Contrast(ok), T.(id_var)(ok), ...
    'VariableNames', {'Y', 'Contrast', 'ID'});
if ~isempty(ms_var)
    tbl.MSRate = T.(ms_var)(ok);
    % z-score MS rate for comparable coefficients
    mu = mean(tbl.MSRate, 'omitnan');
    sd = std(tbl.MSRate, 0, 'omitnan');
    if isfinite(sd) && sd > 0
        tbl.MSRate = (tbl.MSRate - mu) / sd;
    end
    formula = 'Y ~ Contrast + MSRate + (1|ID)';
else
    formula = 'Y ~ Contrast + (1|ID)';
end
try
    lme = fitlme(tbl, formula);
catch
    return
end
S.ok = true;
S.formula = formula;
S.coef = lme.Coefficients;
cn = S.coef.Name;
est = S.coef.Estimate;
pv = S.coef.pValue;
S.p = containers.Map;
S.est = containers.Map;
for i = 1:numel(cn)
    S.p(cn{i}) = pv(i);
    S.est(cn{i}) = est(i);
end
end

function print_lme(title_str, S)
fprintf('\n%s\n', title_str);
if ~S.ok
    fprintf('  (not fit)\n');
    return
end
fprintf('  formula: %s\n', S.formula);
keys = S.est.keys;
for i = 1:numel(keys)
    k = keys{i};
    if strcmp(k, '(Intercept)'), continue; end
    fprintf('  %s: b=%.4g, p=%.4g\n', k, S.est(k), S.p(k));
end
end
