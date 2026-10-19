%% GCP microsaccade rate split of GED gamma
%
% Within each participant and each contrast, trials are split at the
% median stimulus window microsaccade rate (0 to 2 s, Hz). Trials equal
% to that median are dropped. A contrast is omitted when either side has
% fewer than min_trials_per_bin trials with a finite GED spectrum.
%
% Gamma metrics are the peak of the bin averaged spectrum, using the same
% helpers as the confirmatory condition averages. The GED spatial filter
% is not refit.
%
% Figures (figures/controls/ms_rate_median_split):
%   Paired low vs high boxes: power, frequency, pupil size
%   High minus low boxes: power, frequency
%   Sanity check spectra: one panel per contrast
%   Sanity check of the microsaccade split: one panel per participant
%
% Inputs:
%   paths.features/GCP_eeg_GED.mat (trials_powratio, scan_freqs)
%   paths.features/<subj>/gaze/gaze_matrix_trial.mat (MSRate, PupilSize_bl)

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar');
addpath('/Users/Arne/Documents/GitHub/GCP/2_feature_extraction/GED-helpers');
nSubj = numel(subjects);
nCond = 4;
condPct = [25 50 75 100];
condFields = {'subj_data_gaze_trial_c25', 'subj_data_gaze_trial_c50', ...
    'subj_data_gaze_trial_c75', 'subj_data_gaze_trial_c100'};
xtickLabs = {'25%', '50%', '75%', '100%'};
condTitles = {'25% Contrast', '50% Contrast', '75% Contrast', '100% Contrast'};
min_trials_per_bin = 8;
peak_halfwidth_hz = 5;

fig_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/controls/ms_rate_median_split';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

fprintf('[MS SPLIT] Included cohort N = %d (%s)\n', nSubj, strjoin(subjects, ', '));

%% Load GED trial spectra
ged = load(fullfile(paths.features, 'GCP_eeg_GED.mat'), ...
    'trials_powratio', 'scan_freqs', 'subjects');
scan_freqs = ged.scan_freqs(:)';
nFreq = numel(scan_freqs);
ged_idx = match_subjects(ged.subjects, subjects);

peak_freq = nan(nCond, nSubj, 2);
peak_power = nan(nCond, nSubj, 2);
pupil_bl = nan(nCond, nSubj, 2);
spec = nan(nCond, nSubj, nFreq, 2);
n_low = nan(nCond, nSubj);
n_high = nan(nCond, nSubj);
n_tie = nan(nCond, nSubj);
n_zero = nan(nCond, nSubj);
n_usable = nan(nCond, nSubj);
ms_median = nan(nCond, nSubj);

%% Median split and bin averaged spectra
for s = 1:nSubj
    gaze_file = fullfile(paths.features, subjects{s}, 'gaze', 'gaze_matrix_trial.mat');
    if ~isfile(gaze_file)
        error('GCP_ms_split_gamma:MissingGaze', ...
            'Trial gaze file missing for %s: %s', subjects{s}, gaze_file);
    end
    gaze = load(gaze_file);
    gi = ged_idx(s);

    for c = 1:nCond
        pr = ged.trials_powratio{c, gi};
        if isempty(pr)
            fprintf('[MS SPLIT] %s %d%%: empty GED spectra, skipped\n', ...
                subjects{s}, condPct(c));
            continue
        end
        pr = orient_powratio(pr, nFreq);
        nTrl = size(pr, 1);
        if ~isfield(gaze, condFields{c})
            error('GCP_ms_split_gamma:MissingCond', ...
                '%s has no %s', subjects{s}, condFields{c});
        end
        [ms, pup, aligned] = align_gaze_to_spectra(gaze.(condFields{c}), nTrl);
        if ~aligned
            error('GCP_ms_split_gamma:TrialMismatch', ...
                ['%s %d%%: gaze trials do not align with GED spectra ' ...
                '(nEEG = %d).'], subjects{s}, condPct(c), nTrl);
        end

        spec_ok = any(isfinite(pr), 2);
        usable = isfinite(ms) & spec_ok;
        n_usable(c, s) = nnz(usable);
        n_zero(c, s) = nnz(usable & ms == 0);
        if n_usable(c, s) < 1
            fprintf('[MS SPLIT] %s %d%%: no usable trials, skipped\n', ...
                subjects{s}, condPct(c));
            continue
        end

        med = median(ms(usable));
        ms_median(c, s) = med;
        low = usable & ms < med;
        high = usable & ms > med;
        n_low(c, s) = nnz(low);
        n_high(c, s) = nnz(high);
        n_tie(c, s) = nnz(usable & ms == med);
        pct_zero = 100 * n_zero(c, s) / n_usable(c, s);

        if n_low(c, s) < min_trials_per_bin || n_high(c, s) < min_trials_per_bin
            fprintf('[MS SPLIT] %s %d%%: 0 MS %.0f%%, low %d, high %d, ties %d, skipped\n', ...
                subjects{s}, condPct(c), pct_zero, n_low(c, s), n_high(c, s), n_tie(c, s));
            continue
        end

        bins = {low, high};
        for b = 1:2
            avg = compute_condition_average_powratio_ft(pr(bins{b}, :), scan_freqs);
            spec(c, s, :, b) = avg;
            [peak_freq(c, s, b), peak_power(c, s, b)] = pick_tallest_peak( ...
                avg, scan_freqs, 0, peak_halfwidth_hz);
            pupil_bl(c, s, b) = mean(pup(bins{b}), 'omitnan');
        end
        fprintf('[MS SPLIT] %s %d%%: 0 MS %.0f%%, low %d, high %d, ties dropped %d\n', ...
            subjects{s}, condPct(c), pct_zero, n_low(c, s), n_high(c, s), n_tie(c, s));
    end
end

save(fullfile(fig_dir, 'GCP_ms_split.mat'), ...
    'subjects', 'scan_freqs', 'condPct', 'min_trials_per_bin', ...
    'peak_freq', 'peak_power', 'pupil_bl', 'spec', ...
    'n_low', 'n_high', 'n_tie', 'n_zero', 'n_usable', 'ms_median');
fprintf('[MS SPLIT] Saved %s\n', fullfile(fig_dir, 'GCP_ms_split.mat'));

%% Paired boxplots
plot_paired_boxes(peak_power, colors, xtickLabs, 'Power [dB]', false, ...
    fullfile(fig_dir, 'GCP_split_boxplot_Power.png'));
plot_paired_boxes(peak_freq, colors, xtickLabs, 'Frequency [Hz]', false, ...
    fullfile(fig_dir, 'GCP_split_boxplot_Frequency.png'));
plot_paired_boxes(pupil_bl, colors, xtickLabs, 'Pupil Size [%]', true, ...
    fullfile(fig_dir, 'GCP_split_boxplot_PupilSize.png'));

%% High minus low boxplots
plot_diff_boxes(peak_power(:, :, 2) - peak_power(:, :, 1), colors, xtickLabs, ...
    'Power (High-Low) [dB]', ...
    fullfile(fig_dir, 'GCP_split_boxplot_Power_diff.png'));
plot_diff_boxes(peak_freq(:, :, 2) - peak_freq(:, :, 1), colors, xtickLabs, ...
    'Frequency (High-Low) [Hz]', ...
    fullfile(fig_dir, 'GCP_split_boxplot_Frequency_diff.png'));

%% Sanity check spectra
plot_split_spectra(spec, peak_freq, scan_freqs, colors, condTitles, ...
    fullfile(fig_dir, 'GCP_split_powspctrm.png'));

%% Sanity check of the microsaccade split
plot_ms_split_check(n_zero, n_low, n_high, n_tie, n_usable, ms_median, ...
    subjects, xtickLabs, min_trials_per_bin, ...
    fullfile(fig_dir, 'GCP_split_ms_rate_check.png'));

fprintf('[MS SPLIT] Done %s\n', datestr(now));

%% Local functions
function idx = match_subjects(allSubjects, keepSubjects)
allSubjects = cellstr(string(allSubjects(:)));
keepSubjects = cellstr(string(keepSubjects(:)));
idx = nan(1, numel(keepSubjects));
for s = 1:numel(keepSubjects)
    hit = find(strcmp(allSubjects, keepSubjects{s}), 1);
    if isempty(hit)
        error('GCP_ms_split_gamma:SubjectMismatch', ...
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
    error('GCP_ms_split_gamma:FreqMismatch', ...
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

function plot_paired_boxes(Y, colors, xtickLabs, ylabStr, drawZero, outPath)
% Y is contrast x subject x bin, bin 1 = low MS rate, bin 2 = high MS rate.
nCond = size(Y, 1);
nSubj = size(Y, 2);
fontSize = 50;
axisLabelSize = fontSize * 0.9;
yTickFontSize = axisLabelSize;
xLabelFontSize = axisLabelSize;
xTickFontSize = axisLabelSize * (1 - 0.33);
dotSize = 250;
dotAlpha = 0.85;
pairOffset = 0.18;
jitterAmp = 0.08;
boxWidth = 0.28;

rng(1, 'twister');
jit = (rand(nSubj, nCond) - 0.5) * jitterAmp;

fig = new_hidden_figure();
closeFig = onCleanup(@() close_if_open(fig)); %#ok<NASGU>
hold on

drawnColors = zeros(0, 3);
for c = 1:nCond
    for b = 1:2
        yC = Y(c, :, b);
        yC = yC(isfinite(yC));
        if isempty(yC)
            continue
        end
        pos = c + (2 * b - 3) * pairOffset;
        boxplot(yC(:), ones(numel(yC), 1), 'Positions', pos, 'Symbol', '', ...
            'Widths', boxWidth, 'Colors', 'k');
        hold on
        drawnColors(end + 1, :) = colors(c, :);
    end
end
styleCurrentBoxplot(drawnColors);

for s = 1:nSubj
    for c = 1:nCond
        yLow = Y(c, s, 1);
        yHigh = Y(c, s, 2);
        if ~(isfinite(yLow) && isfinite(yHigh))
            continue
        end
        xLow = c - pairOffset + jit(s, c);
        xHigh = c + pairOffset + jit(s, c);
        plot([xLow, xHigh], [yLow, yHigh], '-', ...
            'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    end
end

for c = 1:nCond
    for b = 1:2
        y = Y(c, :, b);
        y = y(:);
        x = c + (2 * b - 3) * pairOffset + jit(:, c);
        x = x(:);
        keep = isfinite(y) & isfinite(x);
        if ~any(keep)
            continue
        end
        scatter(x(keep), y(keep), dotSize, colors(c, :), 'filled', ...
            'MarkerFaceAlpha', dotAlpha);
    end
end

yValid = Y(isfinite(Y));
if drawZero
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    yValid = [yValid(:); 0];
end
if ~isempty(yValid)
    yPad = 0.05 * range(yValid);
    if yPad == 0
        yPad = 1;
    end
    ylim([min(yValid) - yPad, max(yValid) + yPad]);
end

xlim([0.45, nCond + 0.55]);
xticks(1:nCond);
xticklabels(xtickLabs);
hXlab = xlabel('Contrast', 'FontSize', xLabelFontSize);
hYlab = ylabel(ylabStr, 'Interpreter', 'none', 'FontSize', axisLabelSize);
ylNow = ylim;
yTxt = ylNow(1) + 0.02 * (ylNow(2) - ylNow(1));
for c = 1:nCond
    text(c - pairOffset, yTxt, 'Low MS', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 16, 'Color', [0.15 0.15 0.15]);
    text(c + pairOffset, yTxt, 'High MS', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 16, 'Color', [0.15 0.15 0.15]);
end
ax = gca;
ax.XAxis.FontSize = xTickFontSize;
ax.YAxis.FontSize = yTickFontSize;
hXlab.FontSize = xLabelFontSize;
hYlab.FontSize = axisLabelSize;
set(ax, 'Box', 'off');
box off
hold off

save_png(fig, outPath);
fprintf('[MS SPLIT] Saved %s\n', outPath);
end

function plot_diff_boxes(M, colors, xtickLabs, ylabStr, outPath)
% M is contrast x subject. Points are not linked across contrast.
nCond = size(M, 1);
nSubj = size(M, 2);
fontSize = 50;
axisLabelSize = fontSize * 0.9;
yTickFontSize = axisLabelSize;
xLabelFontSize = axisLabelSize;
xTickFontSize = axisLabelSize * (1 - 0.33);
dotSize = 400;
dotAlpha = 0.85;
jitterAmp = 0.4;
boxWidth = 0.55;

rng(1, 'twister');
xJit = nan(nCond, nSubj);
for s = 1:nSubj
    for c = 1:nCond
        if isfinite(M(c, s))
            xJit(c, s) = c + jitterAmp * (rand - 0.5);
        end
    end
end

fig = new_hidden_figure();
closeFig = onCleanup(@() close_if_open(fig)); %#ok<NASGU>
hold on

drawnColors = zeros(0, 3);
for c = 1:nCond
    yC = M(c, isfinite(M(c, :)));
    if isempty(yC)
        continue
    end
    boxplot(yC(:), ones(numel(yC), 1), 'Positions', c, 'Symbol', '', ...
        'Widths', boxWidth, 'Colors', 'k');
    hold on
    drawnColors(end + 1, :) = colors(c, :);
end
styleCurrentBoxplot(drawnColors);
yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

for c = 1:nCond
    keep = isfinite(M(c, :)) & isfinite(xJit(c, :));
    if ~any(keep)
        continue
    end
    scatter(xJit(c, keep), M(c, keep), dotSize, colors(c, :), 'filled', ...
        'MarkerFaceAlpha', dotAlpha);
end

yValid = M(isfinite(M));
yValid = [yValid(:); 0];
if ~isempty(yValid)
    yPad = 0.05 * range(yValid);
    if yPad == 0
        yPad = 1;
    end
    ylim([min(yValid) - yPad, max(yValid) + yPad]);
end
xlim([0.5, nCond + 0.5]);
xticks(1:nCond);
xticklabels(xtickLabs);
hXlab = xlabel('Contrast', 'FontSize', xLabelFontSize);
hYlab = ylabel(ylabStr, 'Interpreter', 'none', 'FontSize', axisLabelSize);
ax = gca;
ax.XAxis.FontSize = xTickFontSize;
ax.YAxis.FontSize = yTickFontSize;
hXlab.FontSize = xLabelFontSize;
hYlab.FontSize = axisLabelSize;
set(ax, 'Box', 'off');
box off
hold off

save_png(fig, outPath);
fprintf('[MS SPLIT] Saved %s\n', outPath);
end

function plot_split_spectra(spec, peak_freq, scan_freqs, colors, condTitles, outPath)
% spec is contrast x subject x frequency x bin (1 low, 2 high).
nCond = size(spec, 1);
nFreq = numel(scan_freqs);
fontSize = 28;
fig = new_hidden_figure();
closeFig = onCleanup(@() close_if_open(fig));

ym = nan(nCond, 2, nFreq);
ye = nan(nCond, 2, nFreq);
for c = 1:nCond
    for b = 1:2
        [ym(c, b, :), ye(c, b, :)] = grand_mean_sem_ft( ...
            squeeze(spec(c, :, :, b)), scan_freqs);
    end
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

lineStyles = {'--', '-'};
binNames = {'Low MS rate', 'High MS rate'};
for c = 1:nCond
    subplot(2, 2, c);
    hold on
    leg_p = gobjects(0);
    leg_lbl = {};
    for b = 1:2
        fx = scan_freqs(:)';
        mu = squeeze(ym(c, b, :))';
        sem = squeeze(ye(c, b, :))';
        if ~any(isfinite(mu))
            continue
        end
        eb = shadedErrorBar(fx, mu, sem, 'lineProps', {lineStyles{b}}, 'transparent', true);
        set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', 3);
        set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.15);
        set(eb.edge(1), 'Color', 'none');
        set(eb.edge(2), 'Color', 'none');
        pk = mean(peak_freq(c, :, b), 'omitnan');
        leg_p(end + 1) = eb.mainLine;
        leg_lbl{end + 1} = binNames{b};
    end
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xlim([scan_freqs(1), scan_freqs(end)]);
    ylim(yLim);
    xticks(30:10:90);
    title(condTitles{c}, 'FontSize', fontSize, 'FontWeight', 'bold');
    xlabel('Frequency [Hz]', 'FontSize', fontSize * 0.85);
    ylabel('Power [dB]', 'FontSize', fontSize * 0.85);
    set(gca, 'FontSize', fontSize * 0.75, 'Box', 'off');
    if ~isempty(leg_p)
        legend(leg_p, leg_lbl, 'Location', 'northeast', 'Box', 'off', ...
            'FontSize', fontSize * 0.7);
    end
    box off
    hold off
end

save_png(fig, outPath);
fprintf('[MS SPLIT] Saved %s\n', outPath);
end

function [mu, sem] = grand_mean_sem_ft(curves, scan_freqs)
% Participant mean and SEM via ft_freqgrandaverage.
% curves is subject x frequency. Rows without a finite spectrum are omitted.
nFreq = numel(scan_freqs);
mu = nan(1, nFreq);
sem = nan(1, nFreq);
pow_per_subj = {};
for i = 1:size(curves, 1)
    row = curves(i, :);
    if ~any(isfinite(row))
        continue
    end
    pow_per_subj{end + 1} = ged_powcurve_to_freq_ft(row, scan_freqs);
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

function styleCurrentBoxplot(boxColors)
hBoxes = findobj(gca, 'Tag', 'Box');
if isempty(hBoxes)
    return
end
nB = numel(hBoxes);
mx = zeros(nB, 1);
for ii = 1:nB
    xd = get(hBoxes(ii), 'XData');
    mx(ii) = mean(xd(:), 'omitnan');
end
[~, ord] = sort(mx);
hBoxes = hBoxes(ord);
for bi = 1:min(numel(hBoxes), size(boxColors, 1))
    patch(get(hBoxes(bi), 'XData'), get(hBoxes(bi), 'YData'), ...
        boxColors(bi, :), 'FaceAlpha', 0.25, ...
        'EdgeColor', boxColors(bi, :), 'LineWidth', 1.5);
end
end

function plot_ms_split_check(n_zero, n_low, n_high, n_tie, n_usable, ms_median, ...
    subjects, xtickLabs, min_trials_per_bin, outPath)
% One panel per participant. Each bar is one contrast.
% Segments are mutually exclusive shares of usable trials:
%   0 MS, low with rate above 0, median ties with rate above 0, high.
% 0 MS trials sit in the low bin when the median is above 0, and in the
% dropped tie when the median is 0 Hz.
nCond = size(n_usable, 1);
nSubj = numel(subjects);
nCol = min(4, max(nSubj, 1));
nRow = max(1, ceil(nSubj / nCol));
stackColors = [
    0.75 0.75 0.75
    0.55 0.75 0.90
    0.90 0.75 0.45
    0.20 0.20 0.20
    ];
stackNames = {'0 MS', 'Low, above 0', 'Tied at median, above 0', 'High'};
fontSize = 16;

fig = new_hidden_figure();
closeFig = onCleanup(@() close_if_open(fig)); %#ok<NASGU>
lastBh = gobjects(0);

for s = 1:nSubj
    subplot(nRow, nCol, s);
    hold on
    pct = nan(nCond, 4);
    excluded = false(nCond, 1);
    for c = 1:nCond
        nu = n_usable(c, s);
        if ~isfinite(nu) || nu < 1
            continue
        end
        z = n_zero(c, s);
        lo = n_low(c, s);
        hi = n_high(c, s);
        ti = n_tie(c, s);
        med = ms_median(c, s);
        zLow = 0;
        zTie = 0;
        zHigh = 0;
        if isfinite(med) && med > 0
            zLow = z;
        elseif isfinite(med) && med == 0
            zTie = z;
        elseif isfinite(med) && med < 0
            zHigh = z;
        end
        if ~isfinite(lo) || ~isfinite(hi) || ~isfinite(ti)
            pct(c, :) = [100 * z / nu, nan, nan, nan];
        else
            pct(c, :) = 100 * [z, lo - zLow, ti - zTie, hi - zHigh] / nu;
        end
        excluded(c) = ~isfinite(lo) || ~isfinite(hi) || ~isfinite(ti) || ...
            lo < min_trials_per_bin || hi < min_trials_per_bin;
    end

    bh = bar(1:nCond, pct, 'stacked', 'BarWidth', 0.72);
    for k = 1:numel(bh)
        bh(k).FaceColor = stackColors(k, :);
        bh(k).EdgeColor = 'none';
    end
    lastBh = bh;
    for c = 1:nCond
        if ~isfinite(pct(c, 1))
            continue
        end
        note = sprintf('%.0f%%\nn=%d', pct(c, 1), n_usable(c, s));
        if excluded(c)
            note = sprintf('%s\nexcluded', note);
        end
        text(c, 103, note, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', [0.15 0.15 0.15]);
    end
    xlim([0.4, nCond + 0.6]);
    ylim([0 128]);
    xticks(1:nCond);
    xticklabels(xtickLabs);
    if s > (nRow - 1) * nCol
        xlabel('Contrast', 'FontSize', fontSize);
    end
    if mod(s - 1, nCol) == 0
        ylabel('Usable trials [%]', 'FontSize', fontSize);
    end
    title(sprintf('Participant %s', subjects{s}), 'FontSize', fontSize, 'FontWeight', 'normal');
    set(gca, 'FontSize', fontSize * 0.85, 'Box', 'off');
    box off
    hold off
end

if nSubj < nRow * nCol
    subplot(nRow, nCol, nSubj + 1);
    hold on
    legH = gobjects(4, 1);
    for k = 1:4
        legH(k) = patch(nan, nan, stackColors(k, :), 'EdgeColor', 'none');
    end
    axis off
    legend(legH, stackNames, 'Location', 'west', 'Box', 'off', 'FontSize', fontSize);
    hold off
elseif ~isempty(lastBh)
    legend(lastBh, stackNames, 'Location', 'northeast', 'Box', 'off', ...
        'FontSize', fontSize * 0.8);
end

save_png(fig, outPath);
fprintf('[MS SPLIT] Saved %s\n', outPath);
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
