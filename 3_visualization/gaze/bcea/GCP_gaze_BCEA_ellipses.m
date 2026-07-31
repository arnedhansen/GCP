%% GCP Gaze BCEA Ellipses
% Computes one gaze dispersion estimate per participant and condition.
% The screen plot shows the group mean BCEA95 ellipse only
% (Mahalanobis radius sqrt(2*k) with k = 2.291, ~2.14 SD; not 2 SD or 3 SD),
% on full screen coordinates, with a textbox of BCEA95 areas.
% Same formula as GCP_gaze_fex.m; values are also saved for LMMs.

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP', 0);
subjects = gcp_subject_inclusion(subjects, paths);
figpath = fullfile(paths.figures, 'gaze', 'bcea');

%% Settings
condVars = {'dataET_c25', 'dataET_c50', 'dataET_c75', 'dataET_c100'};
condValues = [25 50 75 100];
condLabels = {' 25% Contrast', ' 50% Contrast', ' 75% Contrast', ' 100% Contrast'};
nCond = numel(condVars);
nSubj = numel(subjects);

stimWindow = [0 2];
baselineWindow = [-1.5 -0.5];
screenW = 800;
screenH = 600;
blinkWin = 25;
bceaK95 = 2.291;

fontSize = 40;
lineW = 2;

%% Participant and condition estimates
centroidStim = nan(nSubj, nCond, 2);
covStim = nan(2, 2, nSubj, nCond);
bceaStim = nan(nSubj, nCond);
bceaBase = nan(nSubj, nCond);
nValidStim = zeros(nSubj, nCond);
nValidBase = zeros(nSubj, nCond);

fprintf('\n[VIZ GAZE BCEA] Participant and condition estimates\n');
fprintf('[VIZ GAZE BCEA] Stimulus window: %.2f to %.2f s\n', stimWindow);
fprintf('[VIZ GAZE BCEA] Baseline window: %.2f to %.2f s\n', baselineWindow);
fprintf('[VIZ GAZE BCEA] BCEA95 formula: 2 x %.3f x pi x sqrt(det(covariance))\n', bceaK95);

for subj = 1:nSubj
    gazePath = fullfile(paths.features, subjects{subj}, 'gaze');
    dataPath = fullfile(gazePath, 'dataET.mat');
    if ~isfile(dataPath)
        dataPath = fullfile(gazePath, 'dataET');
    end
    clc; fprintf('[VIZ GAZE BCEA] Subject %d/%d: %s\n', subj, nSubj, subjects{subj});
    dat = load(dataPath, condVars{:});

    for c = 1:nCond
        [muStim, thisCovStim, nStim] = conditionMoments( ...
            dat.(condVars{c}), stimWindow, screenW, screenH, blinkWin);
        [~, thisCovBase, nBase] = conditionMoments( ...
            dat.(condVars{c}), baselineWindow, screenW, screenH, blinkWin);

        centroidStim(subj, c, :) = muStim;
        covStim(:, :, subj, c) = thisCovStim;
        nValidStim(subj, c) = nStim;
        nValidBase(subj, c) = nBase;

        if all(isfinite(thisCovStim), 'all') && det(thisCovStim) >= 0
            bceaStim(subj, c) = 2 * bceaK95 * pi * sqrt(det(thisCovStim));
        end
        if all(isfinite(thisCovBase), 'all') && det(thisCovBase) >= 0
            bceaBase(subj, c) = 2 * bceaK95 * pi * sqrt(det(thisCovBase));
        end
    end
end

BCEA_bl = 10 * log10(bceaStim ./ bceaBase);
BCEA_bl(~isfinite(BCEA_bl) | bceaStim <= 0 | bceaBase <= 0) = NaN;

%% Save long format data for LMMs
nRows = nSubj * nCond;
Subject = cell(nRows, 1);
Condition = nan(nRows, 1);
BCEA = nan(nRows, 1);
BaselineBCEA = nan(nRows, 1);
BCEA_bl_long = nan(nRows, 1);
CentroidX = nan(nRows, 1);
CentroidY = nan(nRows, 1);
SDX = nan(nRows, 1);
SDY = nan(nRows, 1);
RhoXY = nan(nRows, 1);
ValidSamples = nan(nRows, 1);
BaselineValidSamples = nan(nRows, 1);

row = 0;
for subj = 1:nSubj
    for c = 1:nCond
        row = row + 1;
        Subject{row} = subjects{subj};
        Condition(row) = condValues(c);
        BCEA(row) = bceaStim(subj, c);
        BaselineBCEA(row) = bceaBase(subj, c);
        BCEA_bl_long(row) = BCEA_bl(subj, c);
        CentroidX(row) = centroidStim(subj, c, 1);
        CentroidY(row) = centroidStim(subj, c, 2);
        thisCov = covStim(:, :, subj, c);
        SDX(row) = sqrt(thisCov(1, 1));
        SDY(row) = sqrt(thisCov(2, 2));
        RhoXY(row) = thisCov(1, 2) / sqrt(thisCov(1, 1) * thisCov(2, 2));
        ValidSamples(row) = nValidStim(subj, c);
        BaselineValidSamples(row) = nValidBase(subj, c);
    end
end

bceaTable = table(Subject, Condition, BCEA, BaselineBCEA, BCEA_bl_long, ...
    CentroidX, CentroidY, SDX, SDY, RhoXY, ValidSamples, BaselineValidSamples, ...
    'VariableNames', {'Subject', 'Condition', 'BCEA', 'BaselineBCEA', 'BCEA_bl', ...
    'CentroidX', 'CentroidY', 'SDX', 'SDY', 'RhoXY', ...
    'ValidSamples', 'BaselineValidSamples'});

csvPath = fullfile(paths.features, 'GCP_gaze_BCEA_subject_condition.csv');
matPath = fullfile(paths.features, 'GCP_gaze_BCEA_subject_condition.mat');
writetable(bceaTable, csvPath);
save(matPath, 'bceaTable', 'bceaStim', 'bceaBase', 'BCEA_bl', ...
    'centroidStim', 'covStim', 'subjects', 'condValues');

%% Group mean within participant BCEA95 ellipses on screen coordinates
% Area = 2*k*pi*sqrt(det(cov)) equals pi*c^2*sqrt(det(cov)) with c = sqrt(2*k).
bceaRadius = sqrt(2 * bceaK95);  % ~2.14 SD for k = 2.291
groupCentroid = nan(nCond, 2);
groupCov = nan(2, 2, nCond);
areaBCEA95 = nan(nCond, 1);

for c = 1:nCond
    validSubj = squeeze(all(isfinite(covStim(:, :, :, c)), [1 2]));
    if ~any(validSubj)
        error('No valid covariance estimates for the %d%% condition.', condValues(c));
    end
    groupCentroid(c, :) = squeeze(mean(centroidStim(validSubj, c, :), 1, 'omitnan'));
    groupCov(:, :, c) = mean(covStim(:, :, validSubj, c), 3, 'omitnan');
    areaBCEA95(c) = 2 * bceaK95 * pi * sqrt(max(det(groupCov(:, :, c)), 0));
end

close all
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on

theta = linspace(0, 2 * pi, 361);
unitCircle = [cos(theta); sin(theta)];
for c = 1:nCond
    [vectors, values] = eig(groupCov(:, :, c));
    axisTransform = vectors * sqrt(max(values, 0));
    ellipse95 = groupCentroid(c, :)' + bceaRadius * axisTransform * unitCircle;

    patch(ellipse95(1, :), ellipse95(2, :), colors(c, :), ...
        'FaceAlpha', 0.125, 'EdgeColor', colors(c, :), ...
        'LineStyle', '--', 'LineWidth', lineW, 'HandleVisibility', 'off');
    plot(groupCentroid(c, 1), groupCentroid(c, 2), 'o', ...
        'MarkerSize', 8, 'MarkerFaceColor', colors(c, :), ...
        'MarkerEdgeColor', 'w', 'HandleVisibility', 'off');
end

plot(400, 300, '+', 'MarkerSize', 18, 'LineWidth', 2, ...
    'Color', 'k', 'HandleVisibility', 'off');

conditionHandles = gobjects(nCond, 1);
for c = 1:nCond
    conditionHandles(c) = patch(nan, nan, colors(c, :), ...
        'FaceAlpha', 0.25, 'EdgeColor', colors(c, :), 'LineWidth', 1.5);
end

axis equal
pbaspect([4 3 1]);
xlim([300 500])
ylim([225 375])
%xlim([0 screenW]);
%ylim([0 screenH]);
xline(400, '--', 'LineWidth', 0.25)
yline(300, '--', 'LineWidth', 0.25)
axis manual
%xticks(0:200:800);
%yticks(0:150:600);
box off
set(gca, 'FontSize', fontSize);
xlabel('Screen Width [px]', 'FontSize', fontSize);
ylabel('Screen Height [px]', 'FontSize', fontSize);
legend(conditionHandles, condLabels, 'Location', 'northeast', ...
    'FontSize', fontSize * 0.55, 'Box', 'off');

areaLines = cell(nCond + 1, 1);
areaLines{1} = 'BCEA95 [px^2]';
for c = 1:nCond
    areaLines{c + 1} = sprintf('%d%%: %.0f', condValues(c), areaBCEA95(c));
end
% annotation('textbox', [0.5 0.5 0.7 0.7], 'String', areaLines, ...
%     'FontName', 'FixedWidth', 'FontSize', fontSize * 0.45, ...
%     'EdgeColor', 'k', 'LineWidth', 1.2, 'BackgroundColor', 'w', ...
%     'Margin', 10, 'VerticalAlignment', 'middle', 'FitBoxToText', 'on');

outFigure = fullfile(figpath, 'GCP_gaze_BCEA_ellipses.png');
exportgraphics(gcf, outFigure, 'Resolution', 600, 'BackgroundColor', 'white');

fprintf('\n[VIZ GAZE BCEA] Saved figure: %s\n', outFigure);
fprintf('[VIZ GAZE BCEA] Saved LMM table: %s\n', csvPath);
fprintf('[VIZ GAZE BCEA] Valid participant estimates by condition: %s\n', ...
    mat2str(sum(isfinite(bceaStim), 1)));
fprintf('[VIZ GAZE BCEA] BCEA95 Mahalanobis radius: %.3f SD (neither 2 nor 3)\n', bceaRadius);
fprintf('[VIZ GAZE BCEA] Group BCEA95 areas [px^2]: %s\n', mat2str(areaBCEA95', 5));
fprintf('[VIZ GAZE BCEA] Median BCEA95 by condition [px^2]: %s\n', ...
    mat2str(median(bceaStim, 1, 'omitnan'), 5));
fprintf('[VIZ GAZE BCEA] Median BCEA_bl by condition [dB]: %s\n', ...
    mat2str(median(BCEA_bl, 1, 'omitnan'), 4));
fprintf('[VIZ GAZE BCEA] Done.\n\n');

%% Local functions
function [centroid, covariance, nValid] = conditionMoments(dataET, timeWindow, screenW, screenH, blinkWin)
allX = cell(1, numel(dataET.trial));
allY = cell(1, numel(dataET.trial));

for trl = 1:numel(dataET.trial)
    idx = dataET.time{trl} >= timeWindow(1) & dataET.time{trl} <= timeWindow(2);
    raw = double(dataET.trial{trl}(:, idx));
    if size(raw, 1) < 2 || isempty(raw)
        continue
    end

    valid = raw(1, :) >= 0 & raw(1, :) <= screenW & ...
        raw(2, :) >= 0 & raw(2, :) <= screenH;
    gaze = nan(3, size(raw, 2));
    gaze(1, valid) = raw(1, valid);
    gaze(2, valid) = screenH - raw(2, valid);
    gaze = remove_blinks(gaze, blinkWin);

    keep = isfinite(gaze(1, :)) & isfinite(gaze(2, :));
    allX{trl} = gaze(1, keep);
    allY{trl} = gaze(2, keep);
end

x = [allX{:}];
y = [allY{:}];
nValid = numel(x);
if nValid < 3
    centroid = [NaN NaN];
    covariance = nan(2);
    return
end

centroid = [mean(x, 'omitnan') mean(y, 'omitnan')];
covariance = cov([x(:) y(:)], 'omitrows');
end
