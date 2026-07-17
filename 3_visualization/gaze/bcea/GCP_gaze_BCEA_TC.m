%% GCP Gaze BCEA Time Course
% Computes time resolved 95% BCEA from gaze samples across trials within
% a centered moving window. One time course is estimated per participant
% and condition, baseline normalized in dB, and plotted as mean plus SEM.

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
baselineWindow = [-1.5 -0.25];
displayWindow = [-0.5 2];
movingWindowSec = 0.100;
screenW = 800;
screenH = 600;
blinkWin = 25;
bceaK95 = 2.291;
minValidPairs = 30;
minValidTrials = 3;

fontSize = 50;
lineW = 4;

%% Compute participant time courses
fprintf('\n=== GCP BCEA time course ===\n');
fprintf('Moving window: %.0f ms\n', movingWindowSec * 1000);
fprintf('Baseline: %.2f to %.2f s\n', baselineWindow);
fprintf('Minimum data per estimate: %d gaze pairs across at least %d trials\n', ...
    minValidPairs, minValidTrials);

tVec = [];
bceaRaw = [];
bceaDB = [];
baselineBCEA = nan(nSubj, nCond);

for subj = 1:nSubj
    gazePath = fullfile(paths.features, subjects{subj}, 'gaze');
    dataPath = fullfile(gazePath, 'dataET.mat');
    if ~isfile(dataPath)
        dataPath = fullfile(gazePath, 'dataET');
    end
    fprintf('Subject %d/%d: %s\n', subj, nSubj, subjects{subj});
    dat = load(dataPath, condVars{:});

    for c = 1:nCond
        [thisTime, thisBCEA] = conditionBceaTimeCourse( ...
            dat.(condVars{c}), computeWindow, movingWindowSec, ...
            screenW, screenH, blinkWin, bceaK95, ...
            minValidPairs, minValidTrials);

        if isempty(tVec)
            tVec = thisTime;
            bceaRaw = nan(nSubj, nCond, numel(tVec));
            bceaDB = nan(nSubj, nCond, numel(tVec));
        elseif numel(thisTime) ~= numel(tVec) || ...
                max(abs(thisTime - tVec)) > 1e-9
            thisBCEA = interp1(thisTime, thisBCEA, tVec, 'linear', NaN);
        end

        baselineIdx = tVec >= baselineWindow(1) & tVec <= baselineWindow(2);
        thisBaseline = mean(thisBCEA(baselineIdx), 'omitnan');
        thisDB = 10 * log10(thisBCEA ./ thisBaseline);
        thisDB(~isfinite(thisDB) | thisBCEA <= 0 | thisBaseline <= 0) = NaN;

        bceaRaw(subj, c, :) = thisBCEA;
        bceaDB(subj, c, :) = thisDB;
        baselineBCEA(subj, c) = thisBaseline;
    end
end

if isempty(tVec)
    error('No valid BCEA time course data were found.');
end

%% Save participant time courses
outData = fullfile(paths.features, 'GCP_gaze_BCEA_timeseries.mat');
save(outData, 'bceaRaw', 'bceaDB', 'baselineBCEA', 'tVec', ...
    'subjects', 'condValues', 'movingWindowSec', 'baselineWindow');

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

for c = 1:nCond
    mu = grandMeanDB(c, displayIdx);
    sem = grandSemDB(c, displayIdx);
    eb = shadedErrorBar(x, mu, sem, 'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', lineW);
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.125);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end

xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
    'LineStyle', '--', 'HandleVisibility', 'off');
yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
    'LineStyle', '--', 'HandleVisibility', 'off');
xlim(displayWindow);
xlabel('Time [s]', 'FontSize', fontSize * 0.8);
ylabel('BCEA [dB]', 'FontSize', fontSize * 0.8);

legendHandles = gobjects(nCond, 1);
for c = 1:nCond
    legendHandles(c) = patch(NaN, NaN, colors(c, :), 'FaceAlpha', 0.25, ...
        'EdgeColor', colors(c, :), 'LineWidth', 1.5);
end
set(gca, 'FontSize', fontSize * 0.8);
legend(legendHandles, condLabels, 'Location', 'northeast', ...
    'FontSize', fontSize * 0.65, 'Box', 'off');
box off
hold off

outFigure = fullfile(figpath, 'GCP_gaze_BCEA_TC_db.png');
exportgraphics(gcf, outFigure, 'Resolution', 600, 'BackgroundColor', 'white');

fprintf('\nSaved figure: %s\n', outFigure);
fprintf('Saved participant time courses: %s\n', outData);
fprintf('Valid participants at time zero by condition: %s\n', ...
    mat2str(squeeze(nValid(:, nearestTimeIndex(tVec, 0)))'));
fprintf('Baseline BCEA medians [px^2]: %s\n', ...
    mat2str(median(baselineBCEA, 1, 'omitnan'), 5));
fprintf('=== GCP BCEA time course done ===\n\n');

%% Local functions
function [tVec, bcea] = conditionBceaTimeCourse(dataET, computeWindow, movingWindowSec, ...
    screenW, screenH, blinkWin, bceaK95, minValidPairs, minValidTrials)

fsample = dataET.fsample;
tVec = computeWindow(1):(1 / fsample):computeWindow(2);
nTrials = numel(dataET.trial);
nTime = numel(tVec);
xTrials = nan(nTrials, nTime);
yTrials = nan(nTrials, nTime);

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

    xTrials(trl, :) = interp1(trialTime, gaze(1, :), tVec, 'linear', NaN);
    yTrials(trl, :) = interp1(trialTime, gaze(2, :), tVec, 'linear', NaN);
end

halfWindow = max(1, round(movingWindowSec * fsample / 2));
bcea = nan(1, nTime);

for iTime = 1:nTime
    idx = max(1, iTime - halfWindow):min(nTime, iTime + halfWindow);
    xWindow = xTrials(:, idx);
    yWindow = yTrials(:, idx);
    validPairs = isfinite(xWindow) & isfinite(yWindow);
    validTrials = sum(any(validPairs, 2));

    x = xWindow(validPairs);
    y = yWindow(validPairs);
    if numel(x) < minValidPairs || validTrials < minValidTrials
        continue
    end

    covariance = cov([x(:) y(:)], 'omitrows');
    covarianceDet = det(covariance);
    if all(isfinite(covariance), 'all') && covarianceDet >= 0
        bcea(iTime) = 2 * bceaK95 * pi * sqrt(covarianceDet);
    end
end
end

function idx = nearestTimeIndex(time, target)
[~, idx] = min(abs(time - target));
end
