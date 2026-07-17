%% GCP Gaze BCEA Ellipses
% Computes one gaze dispersion estimate per participant and condition.
% The screen plot shows group mean within participant covariance ellipses
% at exactly 2 SD and 3 SD. BCEA95 uses the same 95% formula as
% GCP_gaze_fex.m.

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);

figpath = fullfile(paths.figures, 'gaze', 'bcea');
if ~isfolder(figpath), mkdir(figpath); end

condVars = {'dataET_c25', 'dataET_c50', 'dataET_c75', 'dataET_c100'};
condValues = [25 50 75 100];
condLabels = {' 25% Contrast', ' 50% Contrast', ' 75% Contrast', ' 100% Contrast'};
nCond = numel(condVars);
nSubj = numel(subjects);

stimWindow = [0 2];
baselineWindow = [-1.5 -0.25];
screenW = 800;
screenH = 600;
blinkWin = 25;
bceaK95 = 2.291;

fontSize = 40;
lineW = 4;

%% Participant and condition estimates
centroidStim = nan(nSubj, nCond, 2);
covStim = nan(2, 2, nSubj, nCond);
bceaStim = nan(nSubj, nCond);
bceaBase = nan(nSubj, nCond);
nValidStim = zeros(nSubj, nCond);
nValidBase = zeros(nSubj, nCond);

fprintf('\n=== GCP BCEA participant and condition estimates ===\n');
fprintf('Stimulus window: %.2f to %.2f s\n', stimWindow);
fprintf('Baseline window: %.2f to %.2f s\n', baselineWindow);
fprintf('BCEA95 formula: 2 x %.3f x pi x sqrt(det(covariance))\n', bceaK95);

for subj = 1:nSubj
    gazePath = fullfile(paths.features, subjects{subj}, 'gaze');
    dataPath = fullfile(gazePath, 'dataET.mat');
    if ~isfile(dataPath)
        dataPath = fullfile(gazePath, 'dataET');
    end
    fprintf('Subject %d/%d: %s\n', subj, nSubj, subjects{subj});
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

dBBCEA = 10 * log10(bceaStim ./ bceaBase);
dBBCEA(~isfinite(dBBCEA) | bceaStim <= 0 | bceaBase <= 0) = NaN;

%% Save long format data for LMMs
nRows = nSubj * nCond;
Subject = cell(nRows, 1);
Condition = nan(nRows, 1);
BCEA = nan(nRows, 1);
BaselineBCEA = nan(nRows, 1);
dBBCEA_long = nan(nRows, 1);
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
        dBBCEA_long(row) = dBBCEA(subj, c);
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

bceaTable = table(Subject, Condition, BCEA, BaselineBCEA, dBBCEA_long, ...
    CentroidX, CentroidY, SDX, SDY, RhoXY, ValidSamples, BaselineValidSamples, ...
    'VariableNames', {'Subject', 'Condition', 'BCEA', 'BaselineBCEA', 'dBBCEA', ...
    'CentroidX', 'CentroidY', 'SDX', 'SDY', 'RhoXY', ...
    'ValidSamples', 'BaselineValidSamples'});

csvPath = fullfile(paths.features, 'GCP_gaze_BCEA_subject_condition.csv');
matPath = fullfile(paths.features, 'GCP_gaze_BCEA_subject_condition.mat');
writetable(bceaTable, csvPath);
save(matPath, 'bceaTable', 'bceaStim', 'bceaBase', 'dBBCEA', ...
    'centroidStim', 'covStim', 'subjects', 'condValues');

%% Group mean within participant ellipses on screen coordinates
groupCentroid = nan(nCond, 2);
groupCov = nan(2, 2, nCond);

for c = 1:nCond
    validSubj = squeeze(all(isfinite(covStim(:, :, :, c)), [1 2]));
    if ~any(validSubj)
        error('No valid covariance estimates for the %d%% condition.', condValues(c));
    end
    groupCentroid(c, :) = squeeze(mean(centroidStim(validSubj, c, :), 1, 'omitnan'));
    groupCov(:, :, c) = mean(covStim(:, :, validSubj, c), 3, 'omitnan');
end

close all
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on

theta = linspace(0, 2 * pi, 361);
unitCircle = [cos(theta); sin(theta)];
for c = 1:nCond
    [vectors, values] = eig(groupCov(:, :, c));
    axisTransform = vectors * sqrt(max(values, 0));

    ellipse2 = groupCentroid(c, :)' + 2 * axisTransform * unitCircle;
    ellipse3 = groupCentroid(c, :)' + 3 * axisTransform * unitCircle;

    patch(ellipse3(1, :), ellipse3(2, :), colors(c, :), ...
        'FaceAlpha', 0.035, 'EdgeColor', colors(c, :), ...
        'LineStyle', ':', 'LineWidth', lineW, 'HandleVisibility', 'off');
    patch(ellipse2(1, :), ellipse2(2, :), colors(c, :), ...
        'FaceAlpha', 0.125, 'EdgeColor', colors(c, :), ...
        'LineStyle', '-', 'LineWidth', lineW, 'HandleVisibility', 'off');
    plot(groupCentroid(c, 1), groupCentroid(c, 2), 'o', ...
        'MarkerSize', 10, 'MarkerFaceColor', colors(c, :), ...
        'MarkerEdgeColor', 'w', 'HandleVisibility', 'off');
end

plot(400, 300, '+', 'MarkerSize', 18, 'LineWidth', 2, ...
    'Color', 'k', 'HandleVisibility', 'off');

conditionHandles = gobjects(nCond, 1);
for c = 1:nCond
    conditionHandles(c) = patch(nan, nan, colors(c, :), ...
        'FaceAlpha', 0.25, 'EdgeColor', colors(c, :), 'LineWidth', 1.5);
end

xlim([0 screenW]);
ylim([0 screenH]);
xticks(0:200:800);
yticks(0:150:600);
axis equal
box off
set(gca, 'FontSize', fontSize);
xlabel('Screen Width [px]', 'FontSize', fontSize);
ylabel('Screen Height [px]', 'FontSize', fontSize);
legend(conditionHandles, condLabels, 'Location', 'northeast', ...
    'FontSize', fontSize * 0.65, 'Box', 'off');

outFigure = fullfile(figpath, 'GCP_gaze_BCEA_ellipses_2SD_3SD.png');
exportgraphics(gcf, outFigure, 'Resolution', 600, 'BackgroundColor', 'white');

fprintf('\nSaved figure: %s\n', outFigure);
fprintf('Saved LMM table: %s\n', csvPath);
fprintf('Valid participant estimates by condition: %s\n', ...
    mat2str(sum(isfinite(bceaStim), 1)));
fprintf('Median BCEA95 by condition [px^2]: %s\n', ...
    mat2str(median(bceaStim, 1, 'omitnan'), 5));
fprintf('Median dBBCEA by condition [dB]: %s\n', ...
    mat2str(median(dBBCEA, 1, 'omitnan'), 4));
fprintf('=== GCP BCEA done ===\n\n');

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
