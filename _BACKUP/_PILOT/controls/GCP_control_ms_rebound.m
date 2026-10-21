%% Control: microsaccade rate from 350 to 900 ms
% Same raw rate series as
% 3_visualization/gaze/microsaccades/GCP_gaze_microsaccades_TC.m.
% No baseline correction. The boxplot window is fixed at [0.35, 0.90] s.
%
% Boxplots are the subject mean of the display-smoothed rate inside that
% window.
%
% Output: figures/controls/
%   GCP_control_ms_rebound_rate.png
%   GCP_control_ms_rebound_boxplot.png

%% Setup
startup
set(0, 'DefaultFigureVisible', 'off');
[subjects, paths, colors, ~] = setup('GCP');
if ispc
    addpath('C:\Users\Administrator\Documents\GitHub\GCP\_BACKUP\_PILOT');
else
    addpath('/Users/Arne/Documents/GitHub/GCP/_BACKUP/_PILOT');
end
[paths, subjects] = gcp_pilot_paths(paths);
try
    subjects = gcp_subject_inclusion(subjects, paths);
catch
    subjects = subjects(~ismember(str2double(string(subjects)), [602, 604, 608]));
    fprintf('[CTRL MS RB] Inclusion file was locked. Excluding 602, 604 and 608.\n');
end
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar')

datapath = paths.features;
figpath  = fullfile(paths.figures, 'controls');
if ~exist(figpath, 'dir')
    mkdir(figpath);
end
nSubj = length(subjects);

%% Display smoothing (same kernel as the confirmatory time course)
fsample = 500;
sigma_ms_disp = 20;
sigma_samp = max(1, round(sigma_ms_disp / (1000 / fsample)));
kHalf = 3 * sigma_samp;
x_kern = -kHalf:kHalf;
gKernel = exp(-x_kern.^2 / (2 * sigma_samp^2));
gKernel = gKernel ./ sum(gKernel);

t_win = [-0.5 2];
lineW = 4;
fontSize = 50;

condFields = {'msTS_c25', 'msTS_c50', 'msTS_c75', 'msTS_c100'};
condLabels = {' 25% Contrast', ' 50% Contrast', ' 75% Contrast', ' 100% Contrast'};
nConds = length(condFields);

%% Load subject time courses
fprintf('\n[CTRL MS RB] Loading microsaccade rate\n');
subjCurves = cell(nSubj, nConds);
t_vec = [];

for subj = 1:nSubj
    fprintf('[CTRL MS RB] Subject %d/%d (%s)\n', subj, nSubj, subjects{subj});
    spath = fullfile(datapath, subjects{subj}, 'gaze');
    dat = load(fullfile(spath, 'gaze_microsaccade_timeseries.mat'));

    for c = 1:nConds
        msTS = dat.(condFields{c});
        ch = find(strcmp(msTS.label, 'MSRate'), 1, 'first');
        thisTime = msTS.time(:)';
        thisRate = msTS.avg(ch, :);
        thisRateSm = conv(thisRate, gKernel, 'same');
        idxDisp = thisTime >= t_win(1) & thisTime <= t_win(2);
        thisTimeDisp = thisTime(idxDisp);
        thisRateDisp = thisRateSm(idxDisp);

        if isempty(t_vec)
            t_vec = thisTimeDisp;
        elseif numel(thisTimeDisp) ~= numel(t_vec) || any(abs(thisTimeDisp - t_vec) > 1e-6)
            error('GCP_control_ms_rebound:TimeMismatch', ...
                'Microsaccade time axes differ across subjects or conditions.');
        end
        subjCurves{subj, c} = thisRateDisp;
    end
end

%% Assemble subject x time x condition matrix
n_disp = numel(t_vec);
subjRates = nan(nSubj, n_disp, nConds);
for subj = 1:nSubj
    for c = 1:nConds
        subjRates(subj, :, c) = subjCurves{subj, c};
    end
end

grandMean = squeeze(nanmean(subjRates, 1));
nValid_ts = squeeze(sum(~isnan(subjRates), 1));
grandSEM  = squeeze(nanstd(subjRates, 0, 1)) ./ sqrt(max(nValid_ts, 1));
grandSEM(nValid_ts < 2) = NaN;

%% Fixed window, raw rate, no baseline
rebound_win = [0.35 0.90];
winIdx = t_vec >= rebound_win(1) & t_vec <= rebound_win(2);
if ~any(winIdx)
    error('GCP_control_ms_rebound:NoWindow', ...
        'The saved rate does not cover 0.35 to 0.90 s.');
end
fprintf('[CTRL MS RB] Raw rate, no baseline. Window [%.3f, %.3f] s (%d samples).\n', ...
    rebound_win(1), rebound_win(2), nnz(winIdx));

reboundRate = nan(nConds, nSubj);
for subj = 1:nSubj
    for c = 1:nConds
        reboundRate(c, subj) = nanmean(subjRates(subj, winIdx, c));
    end
end

%% Time course with the rebound window marked
close all
figure('Position', [0 0 1512 982], 'Color', 'w', 'Visible', 'off');
hold on

for c = 1:nConds
    eb = shadedErrorBar(t_vec, grandMean(:, c), grandSEM(:, c), ...
        'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', lineW);
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.2);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end

yl = ylim;
rbPatch = patch( ...
    [rebound_win(1) rebound_win(2) rebound_win(2) rebound_win(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.45, 'HandleVisibility', 'off');
uistack(rbPatch, 'bottom');
ylim(yl);

xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--', ...
    'HandleVisibility', 'off');
xlim(t_win);
xlabel('Time [s]', 'FontSize', fontSize * 0.8);
ylabel('Microsaccade Rate [Hz]', 'FontSize', fontSize * 0.8);
leg_p = gobjects(nConds, 1);
for c = 1:nConds
    leg_p(c) = patch(nan, nan, colors(c, :), 'FaceAlpha', 0.33, ...
        'EdgeColor', colors(c, :), 'LineWidth', 1.5);
end
set(gca, 'FontSize', fontSize * 0.8);
legend(leg_p, condLabels, 'Location', 'northeast', 'FontSize', fontSize * 0.65, 'Box', 'off');
box off
hold off
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(figpath, 'GCP_control_ms_rebound_rate.png'), '-dpng', '-r600');

%% Rebound window boxplots
xtickLabs = {'25%', '50%', '75%', '100%'};
axisLabelSize  = fontSize * 0.9;
yTickFontSize  = axisLabelSize;
xLabelFontSize = axisLabelSize;
xTickFontSize  = axisLabelSize * (1 - 0.33);
dotSize  = 400;
dotAlpha = 0.85;
jitter   = 0.4;
boxWidth = 0.55;
subjIDs  = str2double(string(subjects(:)));

y = nan(nConds * nSubj, 1);
condIdx = nan(nConds * nSubj, 1);
idVec = nan(nConds * nSubj, 1);
row = 0;
for s = 1:nSubj
    for c = 1:nConds
        row = row + 1;
        idVec(row) = subjIDs(s);
        condIdx(row) = c;
        y(row) = reboundRate(c, s);
    end
end

close all
figure('Position', [0 0 1512 982], 'Color', 'w', 'Visible', 'off');
hold on

xJit = nan(size(y));
for s = 1:nSubj
    idxSubj = idVec == subjIDs(s);
    for c = 1:nConds
        idxPt = idxSubj & condIdx == c & isfinite(y);
        if any(idxPt)
            xJit(idxPt) = c + jitter * (rand - 0.5);
        end
    end
end

for c = 1:nConds
    idxC = condIdx == c & isfinite(y);
    yC = y(idxC);
    if isempty(yC)
        continue
    end
    boxplot(yC, ones(numel(yC), 1), 'Positions', c, 'Symbol', '', ...
        'Widths', boxWidth, 'Colors', 'k');
end
styleCurrentBoxplot(colors(1:nConds, :));

for s = 1:nSubj
    idxSubj = idVec == subjIDs(s);
    xSubj = xJit(idxSubj);
    ySubj = y(idxSubj);
    condSubjIdx = condIdx(idxSubj);
    valid = isfinite(ySubj) & isfinite(xSubj);
    xSubj = xSubj(valid);
    ySubj = ySubj(valid);
    condSubjIdx = condSubjIdx(valid);
    if numel(ySubj) < 2
        continue
    end
    [~, sortIdx] = sort(condSubjIdx);
    plot(xSubj(sortIdx), ySubj(sortIdx), '-', ...
        'Color', [0.8 0.8 0.8], 'LineWidth', 1);
end

for c = 1:nConds
    idxC = condIdx == c & isfinite(y);
    scatter(xJit(idxC), y(idxC), dotSize, colors(c, :), 'filled', ...
        'MarkerFaceAlpha', dotAlpha);
end

xlim([0.5 nConds + 0.5]);
yValid = y(isfinite(y));
if ~isempty(yValid)
    yPad = 0.05 * range(yValid);
    ylim([min(yValid) - yPad, max(yValid) + yPad]);
end
xticks(1:nConds);
xticklabels(xtickLabs);

hXlab = xlabel('Contrast', 'FontSize', xLabelFontSize);
hYlab = ylabel('Microsaccade Rate [Hz]', 'Interpreter', 'none', ...
    'FontSize', axisLabelSize);
ax = gca;
ax.XAxis.FontSize = xTickFontSize;
ax.YAxis.FontSize = yTickFontSize;
hXlab.FontSize = xLabelFontSize;
hYlab.FontSize = axisLabelSize;
set(ax, 'Box', 'off');
box off
hold off

set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(figpath, 'GCP_control_ms_rebound_boxplot.png'), '-dpng', '-r600');

for c = 1:nConds
    fprintf('[CTRL MS RB] %s 350 to 900 ms mean %.3f Hz (SD %.3f)\n', ...
        strtrim(condLabels{c}), nanmean(reboundRate(c, :)), ...
        nanstd(reboundRate(c, :)));
end
fprintf('[CTRL MS RB] Saved figures in %s\n', figpath);

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
