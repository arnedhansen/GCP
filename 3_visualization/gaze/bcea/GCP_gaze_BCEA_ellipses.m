%% GCP Gaze BCEA Ellipses
% Computes one gaze dispersion estimate per participant and condition.
% Group mean BCEA95 ellipses are saved for full / early / late stimulus windows
% (consumed by GCP_assemble_manuscript_figures.m Figure 6).
% Figure 2: single-subject BCEA95 ellipses for the full window only.
% (Mahalanobis radius sqrt(k) with k = 5.991 for 95% coverage)

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

baselineWindow = [-1.5 -0.5];
screenW = 800;
screenH = 600;
blinkWin = 25;
bceaK95 = 5.991;  % k = -2*ln(1-P) for P = 0.95
bceaRadius = sqrt(bceaK95);  % Mahalanobis radius for 95% ellipse

fontSize = 40;
lineW = 2;

% Window tag, latency [s], and PNG name suffix (full keeps legacy unsuffixed name).
winDefs = { ...
    'full',  [0 2], ''; ...
    'early', [0 1], '_early'; ...
    'late',  [1 2], '_late'};
nWin = size(winDefs, 1);

%% Participant and condition estimates (all windows)
centroidStim = nan(nSubj, nCond, 2, nWin);
covStim = nan(2, 2, nSubj, nCond, nWin);
bceaStim = nan(nSubj, nCond, nWin);
bceaBase = nan(nSubj, nCond);
nValidStim = zeros(nSubj, nCond, nWin);
nValidBase = zeros(nSubj, nCond);

fprintf('\n[VIZ GAZE BCEA] Participant and condition estimates\n');
fprintf('[VIZ GAZE BCEA] Baseline window: %.2f to %.2f s\n', baselineWindow);
fprintf('[VIZ GAZE BCEA] BCEA95 formula: %.3f x pi x sqrt(det(covariance))\n', bceaK95);

for subj = 1:nSubj
    gazePath = fullfile(paths.features, subjects{subj}, 'gaze');
    dataPath = fullfile(gazePath, 'dataET.mat');
    if ~isfile(dataPath)
        dataPath = fullfile(gazePath, 'dataET');
    end
    clc; fprintf('[VIZ GAZE BCEA] Subject %d/%d: %s\n', subj, nSubj, subjects{subj});
    dat = load(dataPath, condVars{:});

    for c = 1:nCond
        [~, thisCovBase, nBase] = conditionMoments( ...
            dat.(condVars{c}), baselineWindow, screenW, screenH, blinkWin);
        nValidBase(subj, c) = nBase;
        if all(isfinite(thisCovBase), 'all') && det(thisCovBase) >= 0
            bceaBase(subj, c) = bceaK95 * pi * sqrt(det(thisCovBase));
        end

        for iWin = 1:nWin
            stimWindow = winDefs{iWin, 2};
            [muStim, thisCovStim, nStim] = conditionMoments( ...
                dat.(condVars{c}), stimWindow, screenW, screenH, blinkWin);
            centroidStim(subj, c, :, iWin) = muStim;
            covStim(:, :, subj, c, iWin) = thisCovStim;
            nValidStim(subj, c, iWin) = nStim;
            if all(isfinite(thisCovStim), 'all') && det(thisCovStim) >= 0
                bceaStim(subj, c, iWin) = bceaK95 * pi * sqrt(det(thisCovStim));
            end
        end
    end
end

% Full-window baselined table (legacy export for LMM / QC)
iFull = 1;
BCEA_bl = 100 * (bceaStim(:, :, iFull) - bceaBase) ./ bceaBase;
BCEA_bl(~isfinite(BCEA_bl) | bceaStim(:, :, iFull) <= 0 | bceaBase <= 0) = NaN;

%% Save BCEA data (full window)
nRows = nSubj * nCond;
Subject = cell(nRows, 1);
Condition = nan(nRows, 1);
BCEA = nan(nRows, 1);
BaselineBCEA = nan(nRows, 1);
BCEA_bl_long = nan(nRows, 1);
CentroidX = nan(nRows, 1);
CentroidY = nan(nRows, 1);
BCEA_Direction = nan(nRows, 1);
BCEA_Eccentricity = nan(nRows, 1);
SDX = nan(nRows, 1);
SDY = nan(nRows, 1);
RhoXY = nan(nRows, 1);
ValidSamples = nan(nRows, 1);
BaselineValidSamples = nan(nRows, 1);

fixXY = [400 300];
ppd = 50;

row = 0;
for subj = 1:nSubj
    for c = 1:nCond
        row = row + 1;
        Subject{row} = subjects{subj};
        Condition(row) = condValues(c);
        BCEA(row) = bceaStim(subj, c, iFull);
        BaselineBCEA(row) = bceaBase(subj, c);
        BCEA_bl_long(row) = BCEA_bl(subj, c);
        CentroidX(row) = centroidStim(subj, c, 1, iFull);
        CentroidY(row) = centroidStim(subj, c, 2, iFull);
        if all(isfinite([CentroidX(row) CentroidY(row)]))
            dx = CentroidX(row) - fixXY(1);
            dy = CentroidY(row) - fixXY(2);
            BCEA_Direction(row) = atan2(dy, dx) * (180 / pi);
            BCEA_Eccentricity(row) = hypot(dx, dy) / ppd;
        end
        thisCov = covStim(:, :, subj, c, iFull);
        SDX(row) = sqrt(thisCov(1, 1));
        SDY(row) = sqrt(thisCov(2, 2));
        RhoXY(row) = thisCov(1, 2) / sqrt(thisCov(1, 1) * thisCov(2, 2));
        ValidSamples(row) = nValidStim(subj, c, iFull);
        BaselineValidSamples(row) = nValidBase(subj, c);
    end
end

bceaTable = table(Subject, Condition, BCEA, BaselineBCEA, BCEA_bl_long, ...
    CentroidX, CentroidY, BCEA_Direction, BCEA_Eccentricity, SDX, SDY, RhoXY, ...
    ValidSamples, BaselineValidSamples, ...
    'VariableNames', {'Subject', 'Condition', 'BCEA', 'BaselineBCEA', 'BCEA_bl', ...
    'CentroidX', 'CentroidY', 'BCEA_Direction', 'BCEA_Eccentricity', 'SDX', 'SDY', 'RhoXY', ...
    'ValidSamples', 'BaselineValidSamples'});

csvPath = fullfile(paths.features, 'GCP_gaze_BCEA_subject_condition.csv');
matPath = fullfile(paths.features, 'GCP_gaze_BCEA_subject_condition.mat');
writetable(bceaTable, csvPath);
save(matPath, 'bceaTable', 'bceaStim', 'bceaBase', 'BCEA_bl', ...
    'centroidStim', 'covStim', 'subjects', 'condValues', 'winDefs');

%% Group mean BCEA95 ellipses per analysis window
for iWin = 1:nWin
    winTag = winDefs{iWin, 1};
    nameSuffix = winDefs{iWin, 3};
    stimWindow = winDefs{iWin, 2};
    fprintf('[VIZ GAZE BCEA] Group ellipses: %s (%.0f-%.0f ms)\n', ...
        winTag, stimWindow(1) * 1000, stimWindow(2) * 1000);

    groupCentroid = nan(nCond, 2);
    groupCov = nan(2, 2, nCond);
    areaBCEA95 = nan(nCond, 1);
    for c = 1:nCond
        validSubj = squeeze(all(isfinite(covStim(:, :, :, c, iWin)), [1 2]));
        if ~any(validSubj)
            error('No valid covariance estimates for the %d%% condition (%s).', ...
                condValues(c), winTag);
        end
        groupCentroid(c, :) = squeeze(mean(centroidStim(validSubj, c, :, iWin), 1, 'omitnan'));
        groupCov(:, :, c) = mean(covStim(:, :, validSubj, c, iWin), 3, 'omitnan');
        areaBCEA95(c) = bceaK95 * pi * sqrt(max(det(groupCov(:, :, c)), 0));
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
            'MarkerSize', 12, 'MarkerFaceColor', colors(c, :), ...
            'MarkerEdgeColor', 'w', 'HandleVisibility', 'off');
    end

    plot(400, 300, '+', 'MarkerSize', 18, 'LineWidth', 2, ...
        'Color', 'k', 'HandleVisibility', 'off');

    conditionHandles = gobjects(nCond, 1);
    for c = 1:nCond
        conditionHandles(c) = patch(nan, nan, colors(c, :), ...
            'FaceAlpha', 0.25, 'EdgeColor', colors(c, :), 'LineWidth', 1.5);
    end

    % Match manuscript sibling panel aspect (~1.83) while keeping true px geometry.
    xlimMain = [265 535];
    ylimMain = [225 375];
    axis equal
    xlim(xlimMain)
    ylim(ylimMain)
    xline(400, '--', 'LineWidth', 0.25)
    yline(300, '--', 'LineWidth', 0.25)
    axis manual
    box on
    set(gca, 'FontSize', fontSize);
    xlabel('Screen Width [px]', 'FontSize', fontSize);
    ylabel('Screen Height [px]', 'FontSize', fontSize);
    axMain = gca;
    legend(axMain, conditionHandles, condLabels, 'Location', 'northeast', ...
        'FontSize', fontSize * 0.75, 'Box', 'off');

    % Top-left screen overview inset inside the visible plot box
    drawnow;
    axMain.Units = 'normalized';
    pos = [0.15 0.2 0.77 0.82];
    dataAsp = diff(xlimMain) / diff(ylimMain);
    axAsp = pos(3) / pos(4);
    if axAsp > dataAsp
        pbH = pos(4);
        pbW = pos(4) * dataAsp;
        pbL = pos(1) + (pos(3) - pbW) / 2;
        pbB = pos(2);
    else
        pbW = pos(3);
        pbH = pos(3) / dataAsp;
        pbL = pos(1);
        pbB = pos(2) + (pos(4) - pbH) / 2;
    end
    insetW = pbW * 0.30;
    insetH = insetW * (screenH / screenW);
    if insetH > pbH * 0.42
        insetH = pbH * 0.42;
        insetW = insetH * (screenW / screenH);
    end
    marginX = pbW * 0.04;
    marginY = pbH * 0.04;
    axInset = axes( ...
        'Parent', gcf, ...
        'Units', 'normalized', ...
        'Position', [pbL + marginX, pbB + pbH - marginY - insetH, insetW, insetH], ...
        'Color', 'w', ...
        'Box', 'on');
    hold(axInset, 'on');
    zoomFace = [0.70 0.70 0.70];
    zoomEdge = [0.35 0.35 0.35];
    patch(axInset, ...
        [xlimMain(1) xlimMain(2) xlimMain(2) xlimMain(1)], ...
        [ylimMain(1) ylimMain(1) ylimMain(2) ylimMain(2)], ...
        zoomFace, 'FaceAlpha', 0.25, 'EdgeColor', zoomEdge, ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(axInset, 400, '--', 'LineWidth', 0.25, 'Color', [0.5 0.5 0.5], ...
        'HandleVisibility', 'off');
    yline(axInset, 300, '--', 'LineWidth', 0.25, 'Color', [0.5 0.5 0.5], ...
        'HandleVisibility', 'off');
    plot(axInset, 400, 300, '+', 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'Color', 'k', 'HandleVisibility', 'off');
    set(axInset, ...
        'XLim', [0 screenW], ...
        'YLim', [0 screenH], ...
        'DataAspectRatio', [1 1 1], ...
        'PlotBoxAspectRatio', [screenW screenH 1], ...
        'XTick', 0:200:screenW, ...
        'YTick', 0:150:screenH, ...
        'FontSize', fontSize * 0.28, ...
        'Color', 'w', ...
        'Box', 'on');
    xlabel(axInset, 'Screen Width [px]', 'FontSize', fontSize * 0.28);
    ylabel(axInset, 'Screen Height [px]', 'FontSize', fontSize * 0.28);
    uistack(axInset, 'top');

    drawnow; pause(0.05);
    outFigure = fullfile(figpath, sprintf('GCP_gaze_BCEA_ellipses%s.png', nameSuffix));
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, outFigure, '-dpng', '-r600');
    fprintf('[VIZ GAZE BCEA] Saved figure: %s\n', outFigure);
    fprintf('[VIZ GAZE BCEA] Group BCEA95 areas [%s, px^2]: %s\n', ...
        winTag, mat2str(areaBCEA95', 5));
end

%% Single-subject BCEA95 ellipses (full window; one subplot per participant)
nCols = 5;
nRows = 2;
nSubjTiles = 8;  % tiles 1-8 for subjects; tiles 9-10 combined = screen overview
fontSizeSub = max(14, round(fontSize * 0.45));
lineWSub = max(1, lineW * 0.75);
limPad = 0.15;  % +/- 15% of ellipse span per axis
theta = linspace(0, 2 * pi, 361);
unitCircle = [cos(theta); sin(theta)];

if nSubj > nSubjTiles
    error('Subject layout reserves tiles 1-%d; nSubj = %d exceeds that.', ...
        nSubjTiles, nSubj);
end

% Shared axis limits from all subject ellipses (+/- 15% of span)
allEllipseX = [];
allEllipseY = [];
ellipseBySubj = cell(nSubj, nCond);
centroidBySubj = nan(nSubj, nCond, 2);
for subj = 1:nSubj
    for c = 1:nCond
        thisCov = covStim(:, :, subj, c, iFull);
        thisMu = squeeze(centroidStim(subj, c, :, iFull))';
        centroidBySubj(subj, c, :) = thisMu;
        if ~(all(isfinite(thisCov), 'all') && all(isfinite(thisMu)))
            continue
        end
        [vectors, values] = eig(thisCov);
        axisTransform = vectors * sqrt(max(values, 0));
        ellipse95 = thisMu' + bceaRadius * axisTransform * unitCircle;
        ellipseBySubj{subj, c} = ellipse95;
        allEllipseX = [allEllipseX, ellipse95(1, :)]; %#ok<AGROW>
        allEllipseY = [allEllipseY, ellipse95(2, :)]; %#ok<AGROW>
    end
end
if isempty(allEllipseX)
    error('No valid subject ellipses to set shared axis limits.');
end
xMin = min(allEllipseX);
xMax = max(allEllipseX);
yMin = min(allEllipseY);
yMax = max(allEllipseY);
xPad = limPad * max(xMax - xMin, eps);
yPad = limPad * max(yMax - yMin, eps);
xlimShared = [xMin - xPad, xMax + xPad];
ylimShared = [yMin - yPad, yMax + yPad];

figure('Position', [0 0 1512 982], 'Color', 'w');
tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');

for subj = 1:nSubj
    nexttile(subj)
    hold on

    for c = 1:nCond
        ellipse95 = ellipseBySubj{subj, c};
        thisMu = squeeze(centroidBySubj(subj, c, :))';
        if isempty(ellipse95)
            continue
        end
        patch(ellipse95(1, :), ellipse95(2, :), colors(c, :), ...
            'FaceAlpha', 0.125, 'EdgeColor', colors(c, :), ...
            'LineStyle', '--', 'LineWidth', lineWSub, 'HandleVisibility', 'off');
        plot(thisMu(1), thisMu(2), 'o', ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(c, :), ...
            'MarkerEdgeColor', 'w', 'HandleVisibility', 'off');
    end
    plot(400, 300, '+', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', 'k', ...
        'HandleVisibility', 'off');
    axis equal
    xlim(xlimShared)
    ylim(ylimShared)
    box on
    set(gca, 'FontSize', fontSizeSub * 0.7);
    title(sprintf('Participant %s', subjects{subj}), 'FontSize', fontSizeSub);
    if subj > 5
        xlabel('Screen Width [px]', 'FontSize', fontSizeSub * 0.7);
    end
    if mod(subj - 1, 5) == 0
        ylabel('Screen Height [px]', 'FontSize', fontSizeSub * 0.7);
    end
end

% Combined screen overview on remaining tiles
axOverview = nexttile(nSubjTiles + 1, [1 2]);
hold(axOverview, 'on');
for c = 1:nCond
    validSubj = squeeze(all(isfinite(covStim(:, :, :, c, iFull)), [1 2]));
    groupCentroid = squeeze(mean(centroidStim(validSubj, c, :, iFull), 1, 'omitnan'));
    groupCov = mean(covStim(:, :, validSubj, c, iFull), 3, 'omitnan');
    [vectors, values] = eig(groupCov);
    axisTransform = vectors * sqrt(max(values, 0));
    ellipse95 = groupCentroid(:) + bceaRadius * axisTransform * unitCircle;
    patch(axOverview, ellipse95(1, :), ellipse95(2, :), colors(c, :), ...
        'FaceAlpha', 0.125, 'EdgeColor', colors(c, :), ...
        'LineStyle', '--', 'LineWidth', lineWSub, 'HandleVisibility', 'off');
    plot(axOverview, groupCentroid(1), groupCentroid(2), 'o', ...
        'MarkerSize', 8, 'MarkerFaceColor', colors(c, :), ...
        'MarkerEdgeColor', 'w', 'HandleVisibility', 'off');
end
plot(axOverview, 400, 300, '+', 'MarkerSize', 12, 'LineWidth', 1.5, ...
    'Color', 'k', 'HandleVisibility', 'off');
axis(axOverview, 'equal');
xlim(axOverview, [0 screenW]);
ylim(axOverview, [0 screenH]);
box(axOverview, 'on');
set(axOverview, 'FontSize', fontSizeSub * 0.7);
title(axOverview, 'Group overview', 'FontSize', fontSizeSub);
xlabel(axOverview, 'Screen Width [px]', 'FontSize', fontSizeSub * 0.7);
ylabel(axOverview, 'Screen Height [px]', 'FontSize', fontSizeSub * 0.7);

conditionHandles = gobjects(nCond, 1);
if nSubj >= nSubjTiles
    for c = 1:nCond
        conditionHandles(c) = patch(nan, nan, colors(c, :), ...
            'FaceAlpha', 0.25, 'EdgeColor', colors(c, :), 'LineWidth', 1.5);
    end
    lgd = legend(conditionHandles, condLabels, 'Location', 'northeast', ...
        'FontSize', fontSizeSub * 0.7, 'Box', 'off');
    lgd.ItemTokenSize = [12 12];
end

outFigureSubj = fullfile(figpath, 'GCP_gaze_BCEA_ellipses_subjects.png');
exportgraphics(gcf, outFigureSubj, 'Resolution', 600, 'BackgroundColor', 'white');

fprintf('[VIZ GAZE BCEA] Saved figure: %s\n', outFigureSubj);
fprintf('[VIZ GAZE BCEA] Saved LMM table: %s\n', csvPath);
fprintf('[VIZ GAZE BCEA] Valid participant estimates by condition (full): %s\n', ...
    mat2str(sum(isfinite(bceaStim(:, :, iFull)), 1)));
fprintf('[VIZ GAZE BCEA] BCEA95 Mahalanobis radius: %.3f (sqrt(k), k = %.3f)\n', bceaRadius, bceaK95);
fprintf('[VIZ GAZE BCEA] Median BCEA95 by condition [px^2]: %s\n', ...
    mat2str(median(bceaStim(:, :, iFull), 1, 'omitnan'), 5));
fprintf('[VIZ GAZE BCEA] Median BCEA_bl by condition [%%]: %s\n', ...
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
