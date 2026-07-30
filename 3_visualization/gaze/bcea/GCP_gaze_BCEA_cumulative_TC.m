%% GCP Gaze BCEA Time Course
% Computes cumulative 95% BCEA from gaze samples within each trial.
% Accumulation starts at the baseline onset (-1.5 s) so that the first
% displayed sample (t = -0.5) is a real BCEA estimate over the full
% baseline window [-1.5, -0.5] (0% change vs that same reference).
% Later samples continue the cumulative window [-1.5, t].
%
% Outlier trials are excluded within each subject x condition using a
% robust MAD rule on the trial-wise peak |% change| after baseline end,
% plus a hard absolute |%| cap.

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar')

figpath = fullfile(paths.figures, 'gaze', 'bcea');
if ~isfolder(figpath), mkdir(figpath); end

condVars = {'dataET_c25', 'dataET_c50', 'dataET_c75', 'dataET_c100'};
condLabels = {' 25% Contrast', ' 50% Contrast', ' 75% Contrast', ' 100% Contrast'};
condValues = [25 50 75 100];
nCond = numel(condVars);
nSubj = numel(subjects);

computeWindow = [-1.5 2];
baselineWindow = [-1.5 -0.5];
cumulStart = baselineWindow(1);   % start before display so t=-0.5 is real data
analysisWindow = [0 2];
displayWindow = [-0.5 2];
screenW = 800;
screenH = 600;
blinkWin = 25;
bceaK95 = 2.291;
minValidSamples = 30;
outlierMadThresh = 5;     % modified z = |x-median| / (1.4826*MAD)
outlierAbsPct = 500;      % hard cap on trial peak |% change|

fontSize = 50;
lineW = 4;

%% Compute participant time courses
fprintf('\n=== GCP BCEA time course (cumulative) ===\n');
fprintf('Cumulative window: [%.2f, t] (starts before display)\n', cumulStart);
fprintf('Display: %.2f to %.2f s; at t=-0.5 BCEA equals full baseline\n', ...
    displayWindow);
fprintf('Baseline reference: %.2f to %.2f s\n', baselineWindow);
fprintf('Minimum samples per estimate: %d\n', minValidSamples);
fprintf('Outlier trials: MAD thresh=%.1f, abs |%%| cap=%.0f\n', ...
    outlierMadThresh, outlierAbsPct);

tVec = [];
bceaRaw = [];
bceaDB = [];
baselineBCEA = nan(nSubj, nCond);
BCEA_bl_summary = nan(nSubj, nCond);
nTrialsKept = nan(nSubj, nCond);
nTrialsOutlier = nan(nSubj, nCond);

for subj = 1:nSubj
    gazePath = fullfile(paths.features, subjects{subj}, 'gaze');
    dataPath = fullfile(gazePath, 'dataET.mat');
    if ~isfile(dataPath)
        dataPath = fullfile(gazePath, 'dataET');
    end
    fprintf('Subject %d/%d: %s\n', subj, nSubj, subjects{subj});
    dat = load(dataPath, condVars{:});

    for c = 1:nCond
        [thisTime, thisBCEA, thisDB, thisBaseline, nKeep, nOut] = ...
            conditionBceaTimeCourse( ...
            dat.(condVars{c}), computeWindow, baselineWindow, cumulStart, ...
            screenW, screenH, blinkWin, bceaK95, minValidSamples, ...
            outlierMadThresh, outlierAbsPct);

        if isempty(tVec)
            tVec = thisTime;
            bceaRaw = nan(nSubj, nCond, numel(tVec));
            bceaDB = nan(nSubj, nCond, numel(tVec));
        elseif numel(thisTime) ~= numel(tVec) || ...
                max(abs(thisTime - tVec)) > 1e-9
            thisBCEA = interp1(thisTime, thisBCEA, tVec, 'linear', NaN);
            thisDB = interp1(thisTime, thisDB, tVec, 'linear', NaN);
        end

        bceaRaw(subj, c, :) = thisBCEA;
        bceaDB(subj, c, :) = thisDB;
        baselineBCEA(subj, c) = thisBaseline;
        nTrialsKept(subj, c) = nKeep;
        nTrialsOutlier(subj, c) = nOut;
        if nOut > 0
            fprintf('  Cond %d%%: excluded %d outlier trial(s), kept %d\n', ...
                condValues(c), nOut, nKeep);
        end

        analysisIdx = tVec >= analysisWindow(1) & tVec <= analysisWindow(2);
        endIdx = find(analysisIdx, 1, 'last');
        if ~isempty(endIdx)
            BCEA_bl_summary(subj, c) = thisDB(endIdx);
        end
    end
end

if isempty(tVec)
    error('No valid BCEA time course data were found.');
end

%% Save participant time courses and TC-derived summaries for boxplots
outData = fullfile(paths.features, 'GCP_gaze_BCEA_timeseries.mat');
save(outData, 'bceaRaw', 'bceaDB', 'baselineBCEA', 'BCEA_bl_summary', 'tVec', ...
    'subjects', 'condValues', 'baselineWindow', 'analysisWindow', 'cumulStart', ...
    'nTrialsKept', 'nTrialsOutlier', 'outlierMadThresh', 'outlierAbsPct');

outSum = fullfile(paths.features, 'GCP_gaze_BCEA_trace_summaries.mat');
save(outSum, 'BCEA_bl_summary', 'subjects', 'condValues', 'analysisWindow');

%% Grand average and SEM
grandMeanDB = squeeze(mean(bceaDB, 1, 'omitnan'));
nValid = squeeze(sum(isfinite(bceaDB), 1));
grandSemDB = squeeze(std(bceaDB, 0, 1, 'omitnan')) ./ sqrt(max(nValid, 1));
grandSemDB(nValid < 2) = NaN;

%% Figure
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on

displayIdx = tVec >= displayWindow(1) & tVec <= displayWindow(2);
x = tVec(displayIdx);

% Ensure display starts on a finite sample (should already be true at -0.5)
for c = 1:nCond
    mu = grandMeanDB(c, displayIdx);
    sem = grandSemDB(c, displayIdx);
    if ~isfinite(mu(1))
        warning('Condition %d: first display sample is non-finite.', c);
    end
    eb = shadedErrorBar(x, mu, sem, 'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', lineW);
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.2);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end

xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
    'LineStyle', '--', 'HandleVisibility', 'off');
yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
    'LineStyle', '--', 'HandleVisibility', 'off');
xlim(displayWindow);
xlabel('Time [s]', 'FontSize', fontSize * 0.8);
ylabel('Cumulative BCEA [%]', 'FontSize', fontSize * 0.8);

legendHandles = gobjects(nCond, 1);
for c = 1:nCond
    legendHandles(c) = patch(NaN, NaN, colors(c, :), 'FaceAlpha', 0.33, ...
        'EdgeColor', colors(c, :), 'LineWidth', 1.5);
end
set(gca, 'FontSize', fontSize * 0.8);
legend(legendHandles, condLabels, 'Location', 'best', ...
    'FontSize', fontSize * 0.65, 'Box', 'off');
box off
hold off

outFigure = fullfile(figpath, 'GCP_gaze_BCEA_TC_db.png');
exportgraphics(gcf, outFigure, 'Resolution', 600, 'BackgroundColor', 'white');

fprintf('\nSaved figure: %s\n', outFigure);
fprintf('Saved participant time courses: %s\n', outData);
fprintf('Saved TC summaries for boxplots: %s\n', outSum);
fprintf('Median BCEA_bl summary (cumulative endpoint) by condition: %s\n', ...
    mat2str(median(BCEA_bl_summary, 1, 'omitnan'), 4));
fprintf('Grand-mean at t=-0.5 by condition: %s\n', ...
    mat2str(grandMeanDB(:, nearestTimeIndex(tVec, -0.5))', 4));
fprintf('Valid participants at t=-0.5 by condition: %s\n', ...
    mat2str(squeeze(nValid(:, nearestTimeIndex(tVec, -0.5)))'));
fprintf('Baseline BCEA medians [px^2]: %s\n', ...
    mat2str(median(baselineBCEA, 1, 'omitnan'), 5));
fprintf('Outlier trials excluded (sum over subjects) by condition: %s\n', ...
    mat2str(sum(nTrialsOutlier, 1, 'omitnan')));
fprintf('=== GCP BCEA time course done ===\n\n');

%% Local functions
function [tVec, bceaMean, bceaPctMean, baselineMean, nKeep, nOut] = ...
    conditionBceaTimeCourse( ...
    dataET, computeWindow, baselineWindow, cumulStart, ...
    screenW, screenH, blinkWin, bceaK95, minValidSamples, ...
    outlierMadThresh, outlierAbsPct)

fsample = dataET.fsample;
tVec = computeWindow(1):(1 / fsample):computeWindow(2);
nTrials = numel(dataET.trial);
nTime = numel(tVec);
startIdx = nearestTimeIndex(tVec, cumulStart);
baselineEndIdx = nearestTimeIndex(tVec, baselineWindow(2));
baselineIdx = tVec >= baselineWindow(1) & tVec <= baselineWindow(2);
postBaseIdx = tVec > baselineWindow(2);

bceaTrials = nan(nTrials, nTime);
pctTrials = nan(nTrials, nTime);
baselineTrials = nan(nTrials, 1);

for trl = 1:nTrials
    raw = double(dataET.trial{trl});
    trialTime = dataET.time{trl};
    if size(raw, 1) < 2 || isempty(raw) || isempty(trialTime)
        continue
    end

    valid = raw(1, :) >= 0 & raw(1, :) <= screenW & ...
        raw(2, :) >= 0 & raw(2, :) <= screenH;
    gaze = nan(3, size(raw, 2));
    gaze(1, valid) = raw(1, valid);
    gaze(2, valid) = screenH - raw(2, valid);
    gaze = remove_blinks(gaze, blinkWin);

    x = interp1(trialTime, gaze(1, :), tVec, 'linear', NaN);
    y = interp1(trialTime, gaze(2, :), tVec, 'linear', NaN);
    good = isfinite(x) & isfinite(y);

    % Full baseline-window BCEA (same definition as GCP_gaze_fex)
    xBl = x(baselineIdx & good);
    yBl = y(baselineIdx & good);
    baselineVal = bceaFromXY(xBl, yBl, bceaK95, minValidSamples);
    baselineTrials(trl) = baselineVal;
    if ~(isfinite(baselineVal) && baselineVal > 0)
        continue
    end

    % Cumulative BCEA from cumulStart (-1.5 s): at t=-0.5 this equals baseline
    n = 0;
    sumX = 0;
    sumY = 0;
    sumX2 = 0;
    sumY2 = 0;
    sumXY = 0;
    for iTime = startIdx:nTime
        xi = x(iTime);
        yi = y(iTime);
        if isfinite(xi) && isfinite(yi)
            n = n + 1;
            sumX = sumX + xi;
            sumY = sumY + yi;
            sumX2 = sumX2 + xi * xi;
            sumY2 = sumY2 + yi * yi;
            sumXY = sumXY + xi * yi;
        end

        bceaVal = bceaFromMoments( ...
            n, sumX, sumY, sumX2, sumY2, sumXY, bceaK95, minValidSamples);
        if ~(isfinite(bceaVal) && bceaVal > 0)
            continue
        end
        bceaTrials(trl, iTime) = bceaVal;
        pctTrials(trl, iTime) = 100 * (bceaVal - baselineVal) / baselineVal;
    end

    % Pin the baseline-end sample to exact 0% (same window as the reference)
    if isfinite(pctTrials(trl, baselineEndIdx))
        pctTrials(trl, baselineEndIdx) = 0;
        bceaTrials(trl, baselineEndIdx) = baselineVal;
    end
end

% Exclude outlier trials: peak |% change| after baseline end
trialPeak = max(abs(pctTrials(:, postBaseIdx)), [], 2, 'omitnan');
keep = isfinite(trialPeak);
if any(keep)
    medPeak = median(trialPeak(keep), 'omitnan');
    madPeak = median(abs(trialPeak(keep) - medPeak), 'omitnan');
    if isfinite(madPeak) && madPeak > 0
        modZ = abs(trialPeak - medPeak) / (1.4826 * madPeak);
        keep = keep & (modZ <= outlierMadThresh);
    end
    keep = keep & (trialPeak <= outlierAbsPct);
end

nKeep = sum(keep);
nOut = sum(isfinite(trialPeak) & ~keep);
pctTrials(~keep, :) = NaN;
bceaTrials(~keep, :) = NaN;
baselineTrials(~keep) = NaN;

bceaMean = mean(bceaTrials, 1, 'omitnan');
bceaPctMean = mean(pctTrials, 1, 'omitnan');
baselineMean = mean(baselineTrials, 1, 'omitnan');
end

function val = bceaFromXY(x, y, bceaK95, minValidSamples)
val = NaN;
if numel(x) < minValidSamples
    return
end
sx = std(x(:), 'omitnan');
sy = std(y(:), 'omitnan');
if ~(isfinite(sx) && isfinite(sy)) || sx <= 0 || sy <= 0
    return
end
rho = corr(x(:), y(:));
if ~isfinite(rho)
    rho = 0;
end
val = 2 * bceaK95 * pi * sx * sy * sqrt(max(0, 1 - rho^2));
end

function val = bceaFromMoments(n, sumX, sumY, sumX2, sumY2, sumXY, bceaK95, minValidSamples)
val = NaN;
if n < minValidSamples
    return
end
vx = (sumX2 - sumX * sumX / n) / (n - 1);
vy = (sumY2 - sumY * sumY / n) / (n - 1);
cxy = (sumXY - sumX * sumY / n) / (n - 1);
if ~(vx > 0 && vy > 0 && isfinite(vx) && isfinite(vy) && isfinite(cxy))
    return
end
sx = sqrt(vx);
sy = sqrt(vy);
rho = cxy / (sx * sy);
if ~isfinite(rho)
    rho = 0;
end
rho = max(-1, min(1, rho));
val = 2 * bceaK95 * pi * sx * sy * sqrt(max(0, 1 - rho^2));
end

function idx = nearestTimeIndex(time, target)
[~, idx] = min(abs(time - target));
end
