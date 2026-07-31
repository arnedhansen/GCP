%% GCP BCEA Explainer Panels
% Real 100% contrast gaze data with BCEA95 ellipses
% (Mahalanobis radius sqrt(2*k) with k = 2.291; same as GCP_gaze_fex).
% Figure 1 (1 x 3): 1 trial / 10 trials / all valid trials for one subject.
% Figure 2 (2 x 5): one randomly chosen trial per included subject.
% Stimulus window [0 2] s.

%% Setup
startup
[subjects, paths, ~, ~] = setup('GCP', 0);
subjects = gcp_subject_inclusion(subjects, paths);

figpath = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/gaze/bcea/explainer';
if ~isfolder(figpath), mkdir(figpath); end

screenW = 800;
screenH = 600;
stimWindow = [0 2];
blinkWin = 25;
minValidSamples = 100;
nTrialsMid = 10;
condVar = 'dataET_c100';
condValue = 100;
bceaK95 = 2.291;  % same as GCP_gaze_fex / BCEA_ellipses
bceaRadius = sqrt(2 * bceaK95);  % Mahalanobis radius ~2.14 SD

fontSize = 40;
lineW = 2;
pointSize = 10;

colPoints = [0.16 0.39 0.71];
colEllipse = [0.84 0.19 0.15];

%% Pick a random participant with enough valid 100% trials
rng('shuffle')
subjOrder = randperm(numel(subjects));

trialPts = {};
pickedSubj = '';

for iSubj = 1:numel(subjOrder)
    subj = subjects{subjOrder(iSubj)};
    gazePath = fullfile(paths.features, subj, 'gaze');
    dataPath = fullfile(gazePath, 'dataET.mat');
    if ~isfile(dataPath)
        dataPath = fullfile(gazePath, 'dataET');
    end
    if ~isfile(dataPath) && ~exist(dataPath, 'file')
        continue
    end

    dat = load(dataPath, condVar);
    dataET = dat.(condVar);
    nTrials = numel(dataET.trial);
    trlOrder = randperm(nTrials);

    ptsList = {};
    for iTrl = 1:nTrials
        pts = trialGazePoints(dataET, trlOrder(iTrl), ...
            stimWindow, screenW, screenH, blinkWin);
        if size(pts, 1) >= minValidSamples
            ptsList{end + 1} = pts; %#ok<AGROW>
        end
    end

    if numel(ptsList) >= nTrialsMid
        trialPts = ptsList;
        pickedSubj = subj;
        break
    end
end

if isempty(trialPts)
    error(['No participant with at least %d valid 100%% trials ' ...
        '(>=%d samples in [%.1f %.1f] s) was found.'], ...
        nTrialsMid, minValidSamples, stimWindow(1), stimWindow(2));
end

nTrialsAll = numel(trialPts);
pts1 = trialPts{1};
pts10 = vertcat(trialPts{1:nTrialsMid});
ptsAll = vertcat(trialPts{:});

panelPts = {pts1, pts10, ptsAll};
panelNTrials = [1, nTrialsMid, nTrialsAll];

fprintf(['[VIZ GAZE BCEA] Subject %s, %d%% contrast: ' ...
    '1 / %d / %d trials (N = %d / %d / %d samples)\n'], ...
    pickedSubj, condValue, nTrialsMid, nTrialsAll, ...
    size(pts1, 1), size(pts10, 1), size(ptsAll, 1));

%% Figure: 1 x 3
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

theta = linspace(0, 2 * pi, 361);
unitCircle = [cos(theta); sin(theta)];

for iPanel = 1:3
    pts = panelPts{iPanel};
    nTrials = panelNTrials(iPanel);
    nSamples = size(pts, 1);

    thisMu = mean(pts, 1, 'omitnan');
    thisCov = cov(pts, 'omitrows');
    [vectors, values] = eig(thisCov);
    axisTransform = vectors * sqrt(max(values, 0));
    ellipse95 = thisMu' + bceaRadius * axisTransform * unitCircle;

    ax = nexttile;
    hold(ax, 'on')

    scatter(ax, pts(:, 1), pts(:, 2), pointSize, ...
        'MarkerFaceColor', colPoints, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.1);

    patch(ax, ellipse95(1, :), ellipse95(2, :), colEllipse, ...
        'FaceAlpha', 0.14, ...
        'EdgeColor', colEllipse, ...
        'LineStyle', '-', ...
        'LineWidth', lineW);

    axis(ax, 'equal')
    pbaspect(ax, [4 3 1]);
    xlim(ax, [0 screenW]);
    ylim(ax, [0 screenH]);
    axis(ax, 'manual')
    xticks(ax, 0:200:800);
    yticks(ax, 0:150:600);
    box(ax, 'on')
    ax.LineWidth = 1.2;
    ax.FontSize = fontSize * 0.36;

    xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize * 0.42);
    if iPanel == 1
        ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize * 0.42);
    end
    title(ax, {sprintf('N_{samples} = %d', nSamples); ...
        sprintf('N_{trials} = %d', nTrials)}, ...
        'FontSize', fontSize * 0.40, 'FontWeight', 'bold');
end

hp = scatter(nan, nan, pointSize, 'MarkerFaceColor', colPoints, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
h95 = plot(nan, nan, 'Color', colEllipse, 'LineWidth', lineW, 'LineStyle', '-');
lgd = legend([hp h95], {'Gaze points', 'BCEA95 ellipse'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off', ...
    'FontSize', fontSize * 0.34);
lgd.Layout.Tile = 'south';

outFigure = fullfile(figpath, sprintf( ...
    'GCP_gaze_BCEA_explainer_3panels_%s_c%d.png', pickedSubj, condValue));
drawnow; pause(0.05);
exportgraphics(gcf, outFigure, 'Resolution', 600, 'BackgroundColor', 'white');
fprintf('[VIZ GAZE BCEA] Saved figure: %s\n', outFigure);

%% Figure: 2 x 5, one random trial per subject
nSubj = numel(subjects);
if nSubj ~= 10
    warning('[VIZ GAZE BCEA] Expected 10 subjects, found %d. Using 2 x 5 layout.', nSubj);
end

subjTrialPts = cell(nSubj, 1);
subjTrialIdx = nan(nSubj, 1);
for iSubj = 1:nSubj
    subj = subjects{iSubj};
    gazePath = fullfile(paths.features, subj, 'gaze');
    dataPath = fullfile(gazePath, 'dataET.mat');
    if ~isfile(dataPath)
        dataPath = fullfile(gazePath, 'dataET');
    end
    if ~isfile(dataPath) && ~exist(dataPath, 'file')
        error('Missing ET data for subject %s.', subj);
    end

    dat = load(dataPath, condVar);
    dataET = dat.(condVar);
    nTrials = numel(dataET.trial);
    trlOrder = randperm(nTrials);

    found = false;
    for iTrl = 1:nTrials
        trl = trlOrder(iTrl);
        pts = trialGazePoints(dataET, trl, stimWindow, screenW, screenH, blinkWin);
        if size(pts, 1) >= minValidSamples
            subjTrialPts{iSubj} = pts;
            subjTrialIdx(iSubj) = trl;
            found = true;
            break
        end
    end
    if ~found
        error('No valid 100%% trial for subject %s.', subj);
    end
    fprintf('[VIZ GAZE BCEA] Subject %s: trial %d, N = %d\n', ...
        subj, subjTrialIdx(iSubj), size(subjTrialPts{iSubj}, 1));
end

figure('Position', [0 0 1512 982], 'Color', 'w');
tiledlayout(2, 5, 'Padding', 'compact', 'TileSpacing', 'compact');

for iSubj = 1:min(nSubj, 10)
    pts = subjTrialPts{iSubj};
    nSamples = size(pts, 1);

    thisMu = mean(pts, 1, 'omitnan');
    thisCov = cov(pts, 'omitrows');
    [vectors, values] = eig(thisCov);
    axisTransform = vectors * sqrt(max(values, 0));
    ellipse95 = thisMu' + bceaRadius * axisTransform * unitCircle;

    ax = nexttile;
    hold(ax, 'on')

    scatter(ax, pts(:, 1), pts(:, 2), pointSize, ...
        'MarkerFaceColor', colPoints, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.1);

    patch(ax, ellipse95(1, :), ellipse95(2, :), colEllipse, ...
        'FaceAlpha', 0.14, ...
        'EdgeColor', colEllipse, ...
        'LineStyle', '-', ...
        'LineWidth', lineW);

    axis(ax, 'equal')
    pbaspect(ax, [4 3 1]);
    xlim(ax, [0 screenW]);
    ylim(ax, [0 screenH]);
    axis(ax, 'manual')
    xticks(ax, 0:400:800);
    yticks(ax, 0:300:600);
    box(ax, 'on')
    ax.LineWidth = 1.0;
    ax.FontSize = fontSize * 0.28;

    if iSubj > 5
        xlabel(ax, 'Width [px]', 'FontSize', fontSize * 0.30);
    end
    if mod(iSubj - 1, 5) == 0
        ylabel(ax, 'Height [px]', 'FontSize', fontSize * 0.30);
    end
    title(ax, {subjects{iSubj}; sprintf('N_{samples} = %d', nSamples)}, ...
        'FontSize', fontSize * 0.30, 'FontWeight', 'bold');
end

hp = scatter(nan, nan, pointSize, 'MarkerFaceColor', colPoints, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
h95 = plot(nan, nan, 'Color', colEllipse, 'LineWidth', lineW, 'LineStyle', '-');
lgd = legend([hp h95], {'Gaze points', 'BCEA95 ellipse'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off', ...
    'FontSize', fontSize * 0.30);
lgd.Layout.Tile = 'south';

outFigureAll = fullfile(figpath, sprintf( ...
    'GCP_gaze_BCEA_explainer_10subj_1trial_c%d.png', condValue));
drawnow; pause(0.05);
exportgraphics(gcf, outFigureAll, 'Resolution', 600, 'BackgroundColor', 'white');
fprintf('[VIZ GAZE BCEA] Saved figure: %s\n', outFigureAll);

%% Local functions
function pts = trialGazePoints(dataET, trl, timeWindow, screenW, screenH, blinkWin)
pts = zeros(0, 2);
raw = double(dataET.trial{trl});
trialTime = dataET.time{trl};
if size(raw, 1) < 2 || isempty(raw) || isempty(trialTime)
    return
end

idx = trialTime >= timeWindow(1) & trialTime <= timeWindow(2);
raw = raw(:, idx);
if isempty(raw)
    return
end

valid = raw(1, :) >= 0 & raw(1, :) <= screenW & ...
    raw(2, :) >= 0 & raw(2, :) <= screenH;
gaze = nan(3, size(raw, 2));
gaze(1, valid) = raw(1, valid);
gaze(2, valid) = screenH - raw(2, valid);
gaze = remove_blinks(gaze, blinkWin);

keep = isfinite(gaze(1, :)) & isfinite(gaze(2, :));
pts = [gaze(1, keep)' gaze(2, keep)'];
end
