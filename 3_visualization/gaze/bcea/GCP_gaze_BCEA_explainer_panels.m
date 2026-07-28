%% GCP BCEA Explainer Panels
% Static 3 panel visualization for explaining convex hull versus BCEA.

%% Setup
startup
[~, paths, ~, ~] = setup('GCP', 0);

figpath = fullfile(paths.figures, 'gaze', 'bcea');
screenW = 800;
screenH = 600;
nThis = 1000;

fontSize = 40;
lineW = 2;
pointSize = 34;

colPoints = [0.16 0.39 0.71];
colHull = [0.98 0.72 0.12];
colEllipse = [0.84 0.19 0.15];
colCentroid = [0.65 0.08 0.09];

% Deterministic synthetic gaze points with central cluster plus sparse outliers.
rng(19)
nTotal = nThis;
nCore = round(0.82 * nTotal);
nOut = nTotal - nCore;

muCore = [400 300];
covCore = [95^2 0.35 * 95 * 72; 0.35 * 95 * 72 72^2];
ptsCore = mvnrnd(muCore, covCore, nCore);
ptsOut = [screenW * rand(nOut, 1), screenH * rand(nOut, 1)];
ptsAll = [ptsCore; ptsOut];
ptsAll(:, 1) = min(max(ptsAll(:, 1), 0), screenW);
ptsAll(:, 2) = min(max(ptsAll(:, 2), 0), screenH);
ptsAll = ptsAll(randperm(size(ptsAll, 1)), :);

%% Figure
close all
figure('Position', [0 0 1512 982*0.6], 'Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

theta = linspace(0, 2 * pi, 361);
unitCircle = [cos(theta); sin(theta)];

pts = ptsAll(1:nThis, :);
thisMu = mean(pts, 1, 'omitnan');
thisCov = cov(pts, 'omitrows');
[vectors, values] = eig(thisCov);
axisTransform = vectors * sqrt(max(values, 0));

ellipse2 = thisMu' + 2 * axisTransform * unitCircle;
ellipse3 = thisMu' + 3 * axisTransform * unitCircle;

hullIdx = convhull(pts(:, 1), pts(:, 2));
hullX = pts(hullIdx, 1);
hullY = pts(hullIdx, 2);

for iPanel = 1:3

    ax = nexttile;
    hold(ax, 'on')

    scatter(ax, pts(:, 1), pts(:, 2), pointSize, ...
        'MarkerFaceColor', colPoints, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.75);
    if iPanel >= 2
        plot(ax, hullX, hullY, ...
            'Color', colHull, ...
            'LineWidth', lineW, ...
            'LineStyle', '-');
    end

    if iPanel == 3
        patch(ax, ellipse3(1, :), ellipse3(2, :), colEllipse, ...
            'FaceAlpha', 0.00, ...
            'EdgeColor', colEllipse, ...
            'LineStyle', ':', ...
            'LineWidth', lineW);

        patch(ax, ellipse2(1, :), ellipse2(2, :), colEllipse, ...
            'FaceAlpha', 0.14, ...
            'EdgeColor', colEllipse, ...
            'LineStyle', '-', ...
            'LineWidth', lineW);

        plot(ax, thisMu(1), thisMu(2), 'o', ...
            'MarkerSize', 9, ...
            'MarkerFaceColor', colCentroid, ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 1.2);

        plot(ax, 400, 300, '+', ...
            'Color', [0.1 0.1 0.1], ...
            'LineWidth', 1.5, ...
            'MarkerSize', 13);
    end

    axis(ax, 'equal')
    pbaspect(ax, [4 3 1]);
    xlim(ax, [0 screenW]);
    ylim(ax, [0 screenH]);
    axis(ax, 'manual')
    xticks(ax, 0:200:800);
    yticks(ax, 0:150:600);
    box(ax, 'on')
    ax.LineWidth = 1.2;
    ax.FontSize = fontSize * 0.40;

    xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize * 0.48);
    ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize * 0.48);
    title(ax, sprintf('N = %d', nThis), 'FontSize', fontSize * 0.46, 'FontWeight', 'bold');
end

hp = scatter(nan, nan, pointSize, 'MarkerFaceColor', colPoints, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.75);
hh = plot(nan, nan, 'Color', colHull, 'LineWidth', lineW, 'LineStyle', '-');
h3 = plot(nan, nan, 'Color', colEllipse, 'LineWidth', lineW, 'LineStyle', ':');
h2 = plot(nan, nan, 'Color', colEllipse, 'LineWidth', lineW, 'LineStyle', '-');
lgd = legend([hp hh h3 h2], {'Gaze points', 'Convex hull', '3 SD ellipse', '2 SD ellipse'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off', ...
    'FontSize', fontSize * 0.34);
lgd.Layout.Tile = 'south';

outFigure = fullfile(figpath, 'GCP_gaze_BCEA_explainer_3panels_N1000.png');
exportgraphics(gcf, outFigure, 'Resolution', 600, 'BackgroundColor', 'white');

fprintf('Saved figure: %s\n', outFigure);
