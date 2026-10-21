%% GCP microsaccade rate regression of GED gamma
%
% Within each participant and each contrast, stimulus window microsaccade
% rate (0 to 2 s, Hz) is the predictor. Rate is centered within that
% participant and contrast, so each slope is a within participant effect.
% A participant contrast is omitted when fewer than min_trials trials are
% usable or the rate does not vary.
%
% Power is the mean of the trial spectrum in a +/- 5 Hz band around the
% peak of that participant contrast mean spectrum (same peak rule as the
% median split). Frequency is the trial peak from
% compute_trial_peak_metrics_from_powratio, smoothed as in the trial level
% GED summary. Pupil size is PupilSize_bl on the same usable trials.
% Spectral slopes are fit independently at each frequency.
%
% Repeated measures summaries, fit separately per contrast:
%   Participant slopes, then a one sample t test against zero
%   fitlme on the trials: Y ~ MSRate + (1|ID), REML, random intercept
% The GED spatial filter is not refit.
%
% Figures (figures/controls/ms_rate_regression):
%   Trial scatters with the within contrast line: one panel per participant
%   Slope spectra: one panel per contrast
%   Microsaccade rate check: one panel per participant
%
% Inputs:
%   paths.features/GCP_eeg_GED.mat (trials_powratio, scan_freqs)
%   paths.features/<subj>/gaze/gaze_matrix_trial.mat (MSRate, PupilSize_bl)

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');
if ispc
    addpath('C:\Users\Administrator\Documents\GitHub\GCP\_BACKUP\_PILOT');
else
    addpath('/Users/Arne/Documents/GitHub/GCP/_BACKUP/_PILOT');
end
[paths, subjects] = gcp_pilot_paths(paths);
subjects = gcp_subject_inclusion(subjects, paths);
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar');
addpath('/Users/Arne/Documents/GitHub/GCP/_BACKUP/_PILOT/2_feature_extraction/GED-helpers');
nSubj = numel(subjects);
nCond = 4;
condPct = [25 50 75 100];
condFields = {'subj_data_gaze_trial_c25', 'subj_data_gaze_trial_c50', ...
    'subj_data_gaze_trial_c75', 'subj_data_gaze_trial_c100'};
xtickLabs = {'25%', '50%', '75%', '100%'};
condTitles = {'25% Contrast', '50% Contrast', '75% Contrast', '100% Contrast'};
min_trials = 10;
peak_halfwidth_hz = 5;
peak_smooth_n = 5;

fig_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/_PILOT/figures/controls/ms_rate_regression';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

fprintf('[MS REG] Included cohort N = %d (%s)\n', nSubj, strjoin(subjects, ', '));

%% Load GED trial spectra
ged = load(fullfile(paths.features, 'GCP_eeg_GED.mat'), ...
    'trials_powratio', 'scan_freqs', 'subjects');
scan_freqs = ged.scan_freqs(:)';
nFreq = numel(scan_freqs);
ged_idx = match_subjects(ged.subjects, subjects);

slope_power = nan(nCond, nSubj);
slope_freq = nan(nCond, nSubj);
slope_pupil = nan(nCond, nSubj);
slope_spec = nan(nCond, nSubj, nFreq);
peak_freq_ref = nan(nCond, nSubj);
n_spec = nan(nCond, nSubj);
n_pow = nan(nCond, nSubj);
n_freq = nan(nCond, nSubj);
n_pupil = nan(nCond, nSubj);
n_zero = nan(nCond, nSubj);
n_usable = nan(nCond, nSubj);
ms_sd = nan(nCond, nSubj);

pow_y = cell(nCond, 1);
pow_x = cell(nCond, 1);
pow_id = cell(nCond, 1);
freq_y = cell(nCond, 1);
freq_x = cell(nCond, 1);
freq_id = cell(nCond, 1);
pup_y = cell(nCond, 1);
pup_x = cell(nCond, 1);
pup_id = cell(nCond, 1);
scatter_pow_x = cell(nCond, nSubj);
scatter_pow_y = cell(nCond, nSubj);
scatter_freq_x = cell(nCond, nSubj);
scatter_freq_y = cell(nCond, nSubj);
scatter_pup_x = cell(nCond, nSubj);
scatter_pup_y = cell(nCond, nSubj);
for c = 1:nCond
    pow_y{c} = [];
    pow_x{c} = [];
    pow_id{c} = [];
    freq_y{c} = [];
    freq_x{c} = [];
    freq_id{c} = [];
    pup_y{c} = [];
    pup_x{c} = [];
    pup_id{c} = [];
end

%% Within participant slopes
for s = 1:nSubj
    gaze_file = fullfile(paths.features, subjects{s}, 'gaze', 'gaze_matrix_trial.mat');
    if ~isfile(gaze_file)
        error('GCP_ms_regression_gamma:MissingGaze', ...
            'Trial gaze file missing for %s: %s', subjects{s}, gaze_file);
    end
    gaze = load(gaze_file);
    gi = ged_idx(s);

    for c = 1:nCond
        pr = ged.trials_powratio{c, gi};
        if isempty(pr)
            fprintf('[MS REG] %s %d%%: empty GED spectra, skipped\n', ...
                subjects{s}, condPct(c));
            continue
        end
        pr = orient_powratio(pr, nFreq);
        nTrl = size(pr, 1);
        if ~isfield(gaze, condFields{c})
            error('GCP_ms_regression_gamma:MissingCond', ...
                '%s has no %s', subjects{s}, condFields{c});
        end
        [ms, pup, aligned] = align_gaze_to_spectra(gaze.(condFields{c}), nTrl);
        if ~aligned
            error('GCP_ms_regression_gamma:TrialMismatch', ...
                ['%s %d%%: gaze trials do not align with GED spectra ' ...
                '(nEEG = %d).'], subjects{s}, condPct(c), nTrl);
        end

        spec_ok = any(isfinite(pr), 2);
        usable = isfinite(ms) & spec_ok;
        n_usable(c, s) = nnz(usable);
        n_zero(c, s) = nnz(usable & ms == 0);
        if n_usable(c, s) >= 2
            ms_sd(c, s) = std(ms(usable), 0);
        end
        if n_usable(c, s) < min_trials || numel(unique(ms(usable))) < 2
            fprintf('[MS REG] %s %d%%: usable %d, 0 MS %.0f%%, skipped\n', ...
                subjects{s}, condPct(c), n_usable(c, s), ...
                100 * n_zero(c, s) / max(n_usable(c, s), 1));
            continue
        end

        ms_u = ms(usable);
        pr_u = pr(usable, :);
        n_spec(c, s) = numel(ms_u);
        slope_spec(c, s, :) = reshape(column_slopes(ms_u, pr_u, min_trials), 1, 1, nFreq);

        avg = compute_condition_average_powratio_ft(pr_u, scan_freqs);
        [f0, ~] = pick_tallest_peak(avg, scan_freqs, 0, 0);
        peak_freq_ref(c, s) = f0;
        band = abs(scan_freqs - f0) <= peak_halfwidth_hz;
        if isfinite(f0) && any(band)
            band_pow = mean(pr_u(:, band), 2, 'omitnan');
            [slope_power(c, s), y_used, x_c, x_raw] = scalar_slope(ms_u, band_pow, min_trials);
            n_pow(c, s) = numel(y_used);
            scatter_pow_x{c, s} = x_raw;
            scatter_pow_y{c, s} = y_used;
            [pow_y{c}, pow_x{c}, pow_id{c}] = append_pool( ...
                pow_y{c}, pow_x{c}, pow_id{c}, y_used, x_c, s);
        end

        [trl_peaks, ~, ~] = compute_trial_peak_metrics_from_powratio( ...
            pr, scan_freqs, true(size(scan_freqs)), peak_smooth_n, peak_halfwidth_hz);
        [slope_freq(c, s), y_used, x_c, x_raw] = scalar_slope(ms_u, trl_peaks(usable), min_trials);
        n_freq(c, s) = numel(y_used);
        scatter_freq_x{c, s} = x_raw;
        scatter_freq_y{c, s} = y_used;
        [freq_y{c}, freq_x{c}, freq_id{c}] = append_pool( ...
            freq_y{c}, freq_x{c}, freq_id{c}, y_used, x_c, s);

        [slope_pupil(c, s), y_used, x_c, x_raw] = scalar_slope(ms_u, pup(usable), min_trials);
        n_pupil(c, s) = numel(y_used);
        scatter_pup_x{c, s} = x_raw;
        scatter_pup_y{c, s} = y_used;
        [pup_y{c}, pup_x{c}, pup_id{c}] = append_pool( ...
            pup_y{c}, pup_x{c}, pup_id{c}, y_used, x_c, s);

        fprintf(['[MS REG] %s %d%%: n=%.0f, 0 MS %.0f%%, SD %.2f Hz, ' ...
            'power n %.0f, frequency n %.0f, pupil n %.0f\n'], ...
            subjects{s}, condPct(c), n_spec(c, s), ...
            100 * n_zero(c, s) / n_spec(c, s), ms_sd(c, s), ...
            n_pow(c, s), n_freq(c, s), n_pupil(c, s));
    end
end

%% Group tests
t_power = slope_group_test(slope_power);
t_freq = slope_group_test(slope_freq);
t_pupil = slope_group_test(slope_pupil);

lme_power = fit_contrast_lmes(pow_y, pow_x, pow_id);
lme_freq = fit_contrast_lmes(freq_y, freq_x, freq_id);
lme_pupil = fit_contrast_lmes(pup_y, pup_x, pup_id);

group_summary = [
    metric_table("Power", condPct, t_power, lme_power)
    metric_table("Frequency", condPct, t_freq, lme_freq)
    metric_table("Pupil", condPct, t_pupil, lme_pupil)
    ];
report_metric('Power', condPct, t_power, lme_power);
report_metric('Frequency', condPct, t_freq, lme_freq);
report_metric('Pupil', condPct, t_pupil, lme_pupil);
disp(group_summary);

save(fullfile(fig_dir, 'GCP_ms_regression.mat'), ...
    'subjects', 'scan_freqs', 'condPct', 'min_trials', ...
    'peak_halfwidth_hz', 'peak_smooth_n', ...
    'slope_power', 'slope_freq', 'slope_pupil', 'slope_spec', 'peak_freq_ref', ...
    'n_spec', 'n_pow', 'n_freq', 'n_pupil', 'n_zero', 'n_usable', 'ms_sd', ...
    't_power', 't_freq', 't_pupil', 'lme_power', 'lme_freq', 'lme_pupil', ...
    'group_summary', ...
    'scatter_pow_x', 'scatter_pow_y', 'scatter_freq_x', 'scatter_freq_y', ...
    'scatter_pup_x', 'scatter_pup_y');
fprintf('[MS REG] Saved %s\n', fullfile(fig_dir, 'GCP_ms_regression.mat'));

%% Trial scatters, one panel per participant
plot_participant_scatters(scatter_pow_x, scatter_pow_y, slope_power, ...
    subjects, colors, xtickLabs, 'Microsaccade rate [Hz]', 'Power [dB]', ...
    fullfile(fig_dir, 'GCP_msreg_scatter_Power.png'));
plot_participant_scatters(scatter_freq_x, scatter_freq_y, slope_freq, ...
    subjects, colors, xtickLabs, 'Microsaccade rate [Hz]', 'Frequency [Hz]', ...
    fullfile(fig_dir, 'GCP_msreg_scatter_Frequency.png'));
plot_participant_scatters(scatter_pup_x, scatter_pup_y, slope_pupil, ...
    subjects, colors, xtickLabs, 'Microsaccade rate [Hz]', 'Pupil Size [%]', ...
    fullfile(fig_dir, 'GCP_msreg_scatter_PupilSize.png'));

%% Slope spectra
plot_slope_spectra(slope_spec, scan_freqs, colors, condTitles, ...
    fullfile(fig_dir, 'GCP_msreg_slope_powspctrm.png'));

%% Microsaccade rate check
plot_ms_rate_check(n_zero, n_usable, n_spec, ms_sd, subjects, xtickLabs, colors, ...
    fullfile(fig_dir, 'GCP_msreg_ms_rate_check.png'));

fprintf('[MS REG] Done %s\n', datestr(now));

%% Local functions
function idx = match_subjects(allSubjects, keepSubjects)
allSubjects = cellstr(string(allSubjects(:)));
keepSubjects = cellstr(string(keepSubjects(:)));
idx = nan(1, numel(keepSubjects));
for s = 1:numel(keepSubjects)
    hit = find(strcmp(allSubjects, keepSubjects{s}), 1);
    if isempty(hit)
        error('GCP_ms_regression_gamma:SubjectMismatch', ...
            'Included subject %s is missing from GCP_eeg_GED.mat.', keepSubjects{s});
    end
    idx(s) = hit;
end
end

function pr = orient_powratio(pr, nFreq)
pr = double(pr);
if size(pr, 2) ~= nFreq && size(pr, 1) == nFreq
    pr = pr.';
end
if size(pr, 2) ~= nFreq
    error('GCP_ms_regression_gamma:FreqMismatch', ...
        'Trial spectra have %d frequencies, scan grid has %d.', size(pr, 2), nFreq);
end
end

function [ms, pup, aligned] = align_gaze_to_spectra(G, nTrl)
ms = nan(nTrl, 1);
pup = nan(nTrl, 1);
aligned = false;
if ~isstruct(G) || ~isfield(G, 'MSRate') || ~isfield(G, 'Trial')
    return
end
trial = double(G.Trial(:));
rate = double(G.MSRate(:));
if isfield(G, 'PupilSize_bl')
    pb = double(G.PupilSize_bl(:));
else
    pb = nan(size(rate));
end
if numel(rate) ~= numel(trial) || numel(pb) ~= numel(trial)
    return
end
if numel(trial) == nTrl && isequal(trial(:).', 1:nTrl)
    ms = rate;
    pup = pb;
    aligned = true;
    return
end
trial = round(trial);
if numel(unique(trial)) == numel(trial) && all(isfinite(trial)) && ...
        all(trial >= 1) && all(trial <= nTrl)
    ms(trial) = rate;
    pup(trial) = pb;
    aligned = true;
end
end

function [b, y_used, x_centered, x_raw] = scalar_slope(x, y, min_trials)
% Ordinary least squares slope. x is centered on the trials that enter.
b = NaN;
y_used = [];
x_centered = [];
x_raw = [];
x = double(x(:));
y = double(y(:));
if numel(x) ~= numel(y)
    error('GCP_ms_regression_gamma:LengthMismatch', ...
        'Predictor and outcome lengths differ.');
end
m = isfinite(x) & isfinite(y);
if nnz(m) < min_trials
    return
end
xx = x(m);
yy = y(m);
xx = xx - mean(xx);
den = sum(xx .^ 2);
if den <= 0
    return
end
b = sum(xx .* (yy - mean(yy))) / den;
y_used = yy;
x_centered = xx;
x_raw = x(m);
end

function b = column_slopes(x, Y, min_trials)
nFreq = size(Y, 2);
b = nan(1, nFreq);
for f = 1:nFreq
    b(f) = scalar_slope(x, Y(:, f), min_trials);
end
end

function [y_all, x_all, id_all] = append_pool(y_all, x_all, id_all, y_used, x_centered, subj_idx)
if isempty(y_used)
    return
end
y_all = [y_all; y_used(:)];
x_all = [x_all; x_centered(:)];
id_all = [id_all; repmat(subj_idx, numel(y_used), 1)];
end

function S = slope_group_test(M)
nCond = size(M, 1);
S = repmat(blank_test(), nCond, 1);
for c = 1:nCond
    v = M(c, :);
    v = v(isfinite(v));
    S(c).n = numel(v);
    if isempty(v)
        continue
    end
    S(c).mean = mean(v);
    if numel(v) < 2
        continue
    end
    S(c).sd = std(v, 0);
    [~, p, ci, stats] = ttest(v);
    S(c).tstat = stats.tstat;
    S(c).df = stats.df;
    S(c).p = p;
    S(c).ci = ci(:)';
end
end

function S = blank_test()
S = struct('n', 0, 'mean', NaN, 'sd', NaN, 'tstat', NaN, 'df', NaN, ...
    'p', NaN, 'ci', [NaN NaN]);
end

function L = fit_contrast_lmes(y_cell, x_cell, id_cell)
nCond = numel(y_cell);
L = repmat(blank_lme(), nCond, 1);
for c = 1:nCond
    L(c) = fit_within_lme(y_cell{c}, x_cell{c}, id_cell{c});
end
end

function row = fit_within_lme(y, x, id)
row = blank_lme();
y = y(:);
x = x(:);
id = id(:);
row.n = numel(y);
if isempty(id)
    row.n_subj = 0;
else
    row.n_subj = numel(unique(id));
end
if row.n < 10 || row.n_subj < 2 || numel(x) ~= numel(y) || numel(id) ~= numel(y)
    return
end
tbl = table(y, x, categorical(id), 'VariableNames', {'Y', 'MSRate', 'ID'});
try
    mdl = fitlme(tbl, 'Y ~ MSRate + (1|ID)', 'FitMethod', 'REML');
    coef = mdl.Coefficients;
    ix = find(string(coef.Name) == "MSRate", 1);
    if isempty(ix)
        return
    end
    row.estimate = coef.Estimate(ix);
    row.se = coef.SE(ix);
    row.df = coef.DF(ix);
    row.tstat = coef.tStat(ix);
    row.p = coef.pValue(ix);
    row.ci = [coef.Lower(ix), coef.Upper(ix)];
    row.converged = true;
catch ME
    fprintf('[MS REG] LME failed: %s\n', ME.message);
end
end

function row = blank_lme()
row = struct('n', 0, 'n_subj', 0, 'estimate', NaN, 'se', NaN, 'df', NaN, ...
    'tstat', NaN, 'p', NaN, 'ci', [NaN NaN], 'converged', false);
end

function T = metric_table(metric, condPct, tS, lS)
n = numel(condPct);
T = table();
T.Metric = repmat(string(metric), n, 1);
T.Contrast = condPct(:);
T.NSubj = arrayfun(@(s) s.n, tS);
T.MeanSlope = arrayfun(@(s) s.mean, tS);
T.SD = arrayfun(@(s) s.sd, tS);
T.t = arrayfun(@(s) s.tstat, tS);
T.df = arrayfun(@(s) s.df, tS);
T.p = arrayfun(@(s) s.p, tS);
T.CI95_low = arrayfun(@(s) s.ci(1), tS);
T.CI95_high = arrayfun(@(s) s.ci(2), tS);
T.LME_beta = arrayfun(@(s) s.estimate, lS);
T.LME_SE = arrayfun(@(s) s.se, lS);
T.LME_df = arrayfun(@(s) s.df, lS);
T.LME_p = arrayfun(@(s) s.p, lS);
T.LME_CI95_low = arrayfun(@(s) s.ci(1), lS);
T.LME_CI95_high = arrayfun(@(s) s.ci(2), lS);
T.LME_N = arrayfun(@(s) s.n, lS);
T.LME_converged = arrayfun(@(s) s.converged, lS);
end

function report_metric(name, condPct, tS, lS)
for c = 1:numel(condPct)
    fprintf(['[MS REG] %s %d%%: mean slope %.4g, t(%g)=%.2f, p=%.4g, ' ...
        'LME %.4g (SE %.4g), p=%.4g, trials %d, subjects %d\n'], ...
        name, condPct(c), tS(c).mean, tS(c).df, tS(c).tstat, tS(c).p, ...
        lS(c).estimate, lS(c).se, lS(c).p, lS(c).n, lS(c).n_subj);
end
end

function plot_participant_scatters(xCell, yCell, slopes, subjects, colors, condLabels, xlabStr, ylabStr, outPath)
% One panel per participant. Points are trials. Each line is the within
% contrast ordinary least squares fit. xCell and yCell are contrast x subject.
nCond = size(xCell, 1);
nSubj = numel(subjects);
nCol = min(4, max(nSubj, 1));
nRow = max(1, ceil(nSubj / nCol));
fontSize = 15;
dotSize = 16;
dotAlpha = 0.35;

xMax = 0;
for s = 1:nSubj
    for c = 1:nCond
        xx = xCell{c, s};
        if ~isempty(xx)
            xMax = max(xMax, max(xx));
        end
    end
end
if ~isfinite(xMax) || xMax <= 0
    xMax = 1;
end
xLim = [0, xMax * 1.04];

fig = new_hidden_figure();
closeFig = onCleanup(@() close_if_open(fig)); %#ok<NASGU>

for s = 1:nSubj
    subplot(nRow, nCol, s);
    hold on
    yAll = [];
    lineY = [];
    for c = 1:nCond
        xx = xCell{c, s};
        yy = yCell{c, s};
        if isempty(xx) || isempty(yy)
            continue
        end
        xx = xx(:);
        yy = yy(:);
        scatter(xx, yy, dotSize, colors(c, :), 'filled', ...
            'MarkerFaceAlpha', dotAlpha, 'MarkerEdgeColor', 'none');
        b = slopes(c, s);
        if isfinite(b) && numel(xx) >= 2
            xg = linspace(min(xx), max(xx), 40);
            yg = mean(yy) + b .* (xg - mean(xx));
            plot(xg, yg, '-', 'Color', colors(c, :), 'LineWidth', 2.5);
            lineY = [lineY, yg]; %#ok<AGROW>
        end
        yAll = [yAll; yy]; %#ok<AGROW>
    end
    yFinite = yAll(isfinite(yAll));
    if numel(yFinite) >= 5
        yLo = prctile(yFinite, 2);
        yHi = prctile(yFinite, 98);
        yShow = [yLo; yHi; lineY(:)];
    else
        yShow = [yFinite; lineY(:)];
    end
    yShow = yShow(isfinite(yShow));
    if ~isempty(yShow)
        yPad = 0.08 * range(yShow);
        if ~isfinite(yPad) || yPad == 0
            yPad = 1;
        end
        ylim([min(yShow) - yPad, max(yShow) + yPad]);
    end
    xlim(xLim);
    title(sprintf('Participant %s', subjects{s}), 'FontSize', fontSize, 'FontWeight', 'normal');
    xlabel(xlabStr, 'FontSize', fontSize * 0.9);
    ylabel(ylabStr, 'FontSize', fontSize * 0.9);
    set(gca, 'FontSize', fontSize * 0.8, 'Box', 'off');
    box off
    hold off
end

if nSubj < nRow * nCol
    subplot(nRow, nCol, nSubj + 1);
    hold on
    legH = gobjects(nCond, 1);
    for c = 1:nCond
        legH(c) = plot(nan, nan, 'o-', 'Color', colors(c, :), ...
            'MarkerFaceColor', colors(c, :), 'LineWidth', 2, 'MarkerSize', 7);
    end
    axis off
    legend(legH, condLabels, 'Location', 'west', 'Box', 'off', 'FontSize', fontSize);
    hold off
end

save_png(fig, outPath);
fprintf('[MS REG] Saved %s\n', outPath);
end

function plot_slope_spectra(spec, scan_freqs, colors, condTitles, outPath)
% spec is contrast x subject x frequency. Values are dB per Hz MS rate.
nCond = size(spec, 1);
nFreq = numel(scan_freqs);
fontSize = 28;
fig = new_hidden_figure();
closeFig = onCleanup(@() close_if_open(fig)); %#ok<NASGU>

ym = nan(nCond, nFreq);
ye = nan(nCond, nFreq);
for c = 1:nCond
    [ym(c, :), ye(c, :)] = grand_mean_sem_ft(squeeze(spec(c, :, :)), scan_freqs);
end
yShow = [ym(:) - ye(:); ym(:) + ye(:); ym(:); 0];
yShow = yShow(isfinite(yShow));
if isempty(yShow)
    yLim = [-1 1];
else
    yPad = 0.05 * range(yShow);
    if yPad == 0
        yPad = 0.5;
    end
    yLim = [min(yShow) - yPad, max(yShow) + yPad];
end

for c = 1:nCond
    subplot(2, 2, c);
    hold on
    fx = scan_freqs(:)';
    mu = ym(c, :);
    sem = ye(c, :);
    if any(isfinite(mu))
        eb = shadedErrorBar(fx, mu, sem, 'lineProps', {'-'}, 'transparent', true);
        set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', 3);
        set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.15);
        set(eb.edge(1), 'Color', 'none');
        set(eb.edge(2), 'Color', 'none');
    end
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    xlim([scan_freqs(1), scan_freqs(end)]);
    ylim(yLim);
    xticks(30:10:90);
    title(condTitles{c}, 'FontSize', fontSize, 'FontWeight', 'bold');
    xlabel('Frequency [Hz]', 'FontSize', fontSize * 0.85);
    ylabel('Power slope [dB per Hz MS]', 'FontSize', fontSize * 0.85);
    set(gca, 'FontSize', fontSize * 0.75, 'Box', 'off');
    box off
    hold off
end

save_png(fig, outPath);
fprintf('[MS REG] Saved %s\n', outPath);
end

function [mu, sem] = grand_mean_sem_ft(curves, scan_freqs)
% Participant mean and SEM via ft_freqgrandaverage.
% curves is subject x frequency. Rows without a finite spectrum are omitted.
nFreq = numel(scan_freqs);
mu = nan(1, nFreq);
sem = nan(1, nFreq);
if isvector(curves)
    curves = curves(:)';
end
pow_per_subj = {};
for i = 1:size(curves, 1)
    row = curves(i, :);
    if ~any(isfinite(row))
        continue
    end
    pow_per_subj{end + 1} = ged_powcurve_to_freq_ft(row, scan_freqs); %#ok<AGROW>
end
if isempty(pow_per_subj)
    return
end
cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter = 'powspctrm';
ga = ft_freqgrandaverage(cfg, pow_per_subj{:});
nf = numel(ga.freq);
pow = reshape(double(ga.powspctrm), size(ga.powspctrm, 1), nf);
n_fin = sum(isfinite(pow), 1);
mu = mean(pow, 1, 'omitnan');
sem = std(pow, 0, 1, 'omitnan') ./ sqrt(max(1, n_fin));
sem(n_fin < 2) = NaN;
end

function plot_ms_rate_check(n_zero, n_usable, n_spec, ms_sd, subjects, xtickLabs, colors, outPath)
% One panel per participant. Bars are the share of usable trials with no
% microsaccade. A contrast omitted from the spectral regression is marked.
nCond = size(n_usable, 1);
nSubj = numel(subjects);
nCol = min(4, max(nSubj, 1));
nRow = max(1, ceil(nSubj / nCol));
fontSize = 16;

fig = new_hidden_figure();
closeFig = onCleanup(@() close_if_open(fig)); %#ok<NASGU>

for s = 1:nSubj
    subplot(nRow, nCol, s);
    hold on
    pct0 = nan(nCond, 1);
    for c = 1:nCond
        nu = n_usable(c, s);
        if ~isfinite(nu) || nu < 1
            continue
        end
        pct0(c) = 100 * n_zero(c, s) / nu;
        bar(c, pct0(c), 0.7, 'FaceColor', colors(c, :), 'EdgeColor', 'none');
        note = sprintf('n=%d', nu);
        if isfinite(ms_sd(c, s))
            note = sprintf('%s\n%.1f Hz', note, ms_sd(c, s));
        end
        if ~isfinite(n_spec(c, s))
            note = sprintf('%s\nexcluded', note);
        end
        text(c, pct0(c) + 1.5, note, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', [0.15 0.15 0.15]);
    end
    xlim([0.4, nCond + 0.6]);
    ylim([0, 100]);
    xticks(1:nCond);
    xticklabels(xtickLabs);
    if s > (nRow - 1) * nCol
        xlabel('Contrast', 'FontSize', fontSize);
    end
    if mod(s - 1, nCol) == 0
        ylabel('Trials with 0 MS [%]', 'FontSize', fontSize);
    end
    title(sprintf('Participant %s', subjects{s}), 'FontSize', fontSize, 'FontWeight', 'normal');
    set(gca, 'FontSize', fontSize * 0.85, 'Box', 'off');
    box off
    hold off
end

save_png(fig, outPath);
fprintf('[MS REG] Saved %s\n', outPath);
end

function fig = new_hidden_figure()
fig = figure('Visible', 'off', 'Position', [0 0 1512 982], 'Color', 'w', ...
    'MenuBar', 'none', 'ToolBar', 'none', 'DockControls', 'off');
set(fig, 'DefaultAxesToolbarVisible', 'off');
end

function close_if_open(fig)
if isgraphics(fig)
    close(fig);
end
end

function save_png(fig, outPath)
set(fig, 'Visible', 'off');
axs = findall(fig, 'Type', 'axes');
for iAx = 1:numel(axs)
    if isprop(axs(iAx), 'Toolbar') && ~isempty(axs(iAx).Toolbar)
        axs(iAx).Toolbar.Visible = 'off';
    end
end
set(fig, 'PaperPositionMode', 'auto');
print(fig, outPath, '-dpng', '-r600');
end
