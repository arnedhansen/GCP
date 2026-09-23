%% GCP Assemble Manuscript Figures

%% Setup
startup
[~, paths, ~, ~] = setup('GCP', 0);
outDir = fullfile(paths.figures, 'manuscript');
powerDir = fullfile(paths.figures, 'power_analysis');
boxplotDir = fullfile(paths.figures, 'stats', 'boxplots');
gazeDir = fullfile(paths.figures, 'gaze');
eegDir = fullfile(paths.figures, 'eeg');

exportDpi = 300;
titleFontSize = 20;
denseTitleFontSize = 16;
overallTitleFontSize = 24;
statsDir = fullfile(paths.figures, 'stats');
fprintf('[VIZ MANUSCRIPT] Assembling GCP manuscript figures...\n');

%% Figure 1: Contrast detection paradigm
fig1 = figureSpec('Figure1_paradigm', 1, 1, [0 0 1512 982]);
fig1.titleFontSize = titleFontSize;
fig1.overallTitle = 'Grating Task Paradigm';
fig1.overallTitleFontSize = overallTitleFontSize;
fig1.panels = {panelSpec(fullfile(paths.figures, 'paradigm', 'GCP_paradigm.png'), '', '')};

%% Figure 2: Power analyses
fig2 = figureSpec('Figure2_power_analysis', 2, 2, [0 0 1512 982]);
fig2.colGap = 0.035;
fig2.rowGap = 0.02;
fig2.forceFill = false;
fig2.equalPanelHeight = true;
fig2.titleFontSize = denseTitleFontSize;
fig2.overallTitle = 'Gamma Power and Frequency Power Analysis';
fig2.overallTitleFontSize = overallTitleFontSize;
fig2.panels = {
    panelSpec(fullfile(powerDir, 'GCP_power_analysis_power_heatmap.png'), ...
        'A', 'Gamma Power Statistical Power Heatmap', 1);
    panelSpec(fullfile(powerDir, 'GCP_power_analysis_power.png'), ...
        'B', 'Gamma Power Statistical Power Curve', 2);
    panelSpec(fullfile(powerDir, 'GCP_power_analysis_frequency_heatmap.png'), ...
        'C', 'Gamma Frequency Statistical Power Heatmap', 3);
    panelSpec(fullfile(powerDir, 'GCP_power_analysis_frequency.png'), ...
        'D', 'Gamma Frequency Statistical Power Curve', 4)
    };
for iPanel = 1:numel(fig2.panels)
    fig2.panels{iPanel}.trimWhite = true;
end

%% Figure 3: Gaze TCs and Boxplots (full window)
fig3 = figureSpec('Figure3_gaze', 4, 2, [0 0 1512 round(982 * 2)]);
fig3.titleFontSize = denseTitleFontSize;
fig3.overallTitle = 'Gaze Measures Time Courses and Boxplots';
fig3.overallTitleFontSize = overallTitleFontSize;
fig3.colGap = 0.018;
fig3.rowGap = 0.015;
fig3.packAdjacent = true;
fig3.fixedColumnLayout = true;
fig3.panels = {
    panelSpec(fullfile(gazeDir, 'pupil', 'GCP_gaze_pupil_size_TC.png'), ...
        'A', 'Pupil Size Time Course', 1);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_PupilSize_bl.png'), ...
        'B', 'Pupil Size Boxplots', 2);
    panelSpec(fullfile(gazeDir, 'microsaccades', 'GCP_gaze_microsaccades_rate.png'), ...
        'C', 'Microsaccade Rate Time Course', 3);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_MSRate_bl.png'), ...
        'D', 'Microsaccade Rate Boxplots', 4);
    panelSpec(fullfile(gazeDir, 'bcea', 'GCP_gaze_BCEA_ellipses.png'), ...
        'E', 'Gaze Dispersion Ellipses', 5);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_BCEA_bl.png'), ...
        'F', 'Gaze Dispersion Boxplots', 6);
    panelSpec(fullfile(gazeDir, 'velocity', 'GCP_gaze_velocity_Vel2D_TC.png'), ...
        'G', 'Eye Velocity Time Course', 7);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_Vel2D_bl.png'), ...
        'H', 'Eye Velocity Boxplots', 8)
    };
for iPanel = 1:numel(fig3.panels)
    fig3.panels{iPanel}.trimWhite = true;
end

%% Figure 4: F-tests on Gaze Measures
fig4 = figureSpec('Figure4_gaze_ftests', 1, 1, [0 0 1512 982]);
fig4.titleFontSize = titleFontSize;
fig4.overallTitle = 'F-tests on Gaze Measure Time Courses';
fig4.overallTitleFontSize = overallTitleFontSize;
fig4.panels = {
    panelSpec(fullfile(statsDir, 'f-tests', 'GCP_stats_gaze_ftests.png'), ...
        '', '')
    };

%% Figure 5: Power Spectrum
fig5 = figureSpec('Figure5_powspctrm', 1, 1, [0 0 1512 982]);
fig5.titleFontSize = titleFontSize;
fig5.overallTitle = 'Power Spectrum';
fig5.overallTitleFontSize = overallTitleFontSize;
fig5.panels = {
    panelSpec(fullfile(eegDir, 'powspctrm', 'GCP_eeg_GED_powspctrm_grand_average.png'), ...
        '', '')
    };

%% Figure 6: Gamma power and frequency boxplots (full window)
fig6 = figureSpec('Figure6_gamma_boxplots', 1, 2, [0 0 1512 982]);
fig6.titleFontSize = titleFontSize;
fig6.overallTitle = 'Gamma Power and Frequency Boxplots';
fig6.overallTitleFontSize = overallTitleFontSize;
fig6.panels = {
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_Power.png'), ...
        'A', 'Gamma Power', 1);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_Frequency.png'), ...
        'B', 'Gamma Peak Frequency', 2)
    };
for iPanel = 1:numel(fig6.panels)
    fig6.panels{iPanel}.trimWhite = true;
end

%% Figure 7: Gamma GED time frequency representation
fig7 = figureSpec('Figure7_tfr', 1, 1, [0 0 1512 982]);
fig7.titleFontSize = titleFontSize;
fig7.overallTitle = 'Time-Frequency Representations';
fig7.overallTitleFontSize = overallTitleFontSize;
fig7.panels = {
    panelSpec(fullfile(eegDir, 'tfr', 'GCP_eeg_tfr_GED.png'), ...
        '', '')
    };

%% Supplementary Figure S1: Combined GED components
figS1 = figureSpec('FigureS1_ged_components', 1, 1, [0 0 1512 982]);
figS1.titleFontSize = titleFontSize;
figS1.titleGap = 0.05;
figS1.panels = {
    panelSpec(fullfile(eegDir, 'ged', 'component_selection', ...
        'GCP_eeg_GED_components_allsubjects.png'), ...
        '', 'Combined GED Components')
    };

%% Supplementary Figure S2: Single participant power spectra
figS2 = figureSpec('FigureS2_powspctrm_subjects', 1, 1, [0 0 1512 982]);
figS2.titleFontSize = titleFontSize;
figS2.titleGap = 0.005;
figS2.panels = {
    panelSpec(fullfile(eegDir, 'powspctrm', ...
        'GCP_eeg_GED_powspctrm_overview_subjects.png'), ...
        '', 'Single Subject Powerspectra')
    };

%% Assemble all figures
figSpecs = {fig1, fig2, fig3, fig4, fig5, fig6, fig7, figS1, figS2};
for iFig = 1:numel(figSpecs)
    spec = figSpecs{iFig};
    outPng = fullfile(outDir, ['GCP_manuscript_' spec.name '.png']);
    assembleManuscriptFigure(spec, outPng, exportDpi);
    fprintf('[VIZ MANUSCRIPT] Saved %s\n', outPng);
end
disp(datestr(now))

%% Local functions
function spec = figureSpec(name, nrow, ncol, figSize)
spec = struct();
spec.name = name;
spec.nrow = nrow;
spec.ncol = ncol;
spec.figSize = figSize;
spec.panels = {};
end

function p = panelSpec(imagePath, letter, titleText, tile, trimWhite)
if nargin < 4 || isempty(tile)
    tile = 1;
end
p = struct( ...
    'file', imagePath, ...
    'letter', letter, ...
    'title', titleText, ...
    'tile', tile, ...
    'trimWhite', false, ...
    'imageScale', 1, ...
    'displayScale', 1);
if nargin >= 5 && ~isempty(trimWhite)
    p.trimWhite = trimWhite;
end
end

function assembleManuscriptFigure(spec, outPng, exportDpi)
missing = {};
for iPanel = 1:numel(spec.panels)
    if ~isfile(spec.panels{iPanel}.file)
        missing{end + 1} = spec.panels{iPanel}.file; %#ok<AGROW>
    end
end

if ~isempty(missing)
    fprintf('\n[VIZ MANUSCRIPT] Missing panel files for %s:\n', spec.name);
    for iMissing = 1:numel(missing)
        fprintf('  %s\n', missing{iMissing});
    end
    error('GCP_assemble_manuscript_figures:MissingPanels', ...
        'Cannot assemble %s until all panel PNG files exist.', spec.name);
end

fig = figure('Visible', 'off', 'Color', 'w');
fig.Units = 'pixels';
fig.Position = spec.figSize;
drawnow;
titleFontSize = 20;
if isfield(spec, 'titleFontSize')
    titleFontSize = spec.titleFontSize;
end
overallTitleFontSize = 24;
if isfield(spec, 'overallTitleFontSize')
    overallTitleFontSize = spec.overallTitleFontSize;
end
% Draw the overall title first, with the same annotation mechanism used for
% the panel titles, so it cannot be affected by later axes creation.
hasOverallTitle = isfield(spec, 'overallTitle') && ~isempty(spec.overallTitle);
overallTitleBand = 0;
overallTitleGap = 0.035;
% Hold the band and the gap below it at a constant pixel height, otherwise a
% taller canvas gets a proportionally larger gap under the title.
bandScale = 982 / spec.figSize(4);
if hasOverallTitle
    overallTitleBand = 0.055 * bandScale;
    overallTitleGap = 0.035 * bandScale;
    annotation(fig, 'textbox', ...
        [0, 1 - overallTitleBand, 1, overallTitleBand], ...
        'String', spec.overallTitle, ...
        'Interpreter', 'tex', ...
        'FontSize', overallTitleFontSize, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FitBoxToText', 'off', ...
        'Margin', 0, ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
end

if isfield(spec, 'colWidths') || spec.ncol == 2
    left = 0.025;
    right = 0.025;
    top = overallTitleGap + overallTitleBand;
    bottom = 0.025;
    colGap = 0.010;
    if isfield(spec, 'colGap')
        colGap = spec.colGap;
    end
    rowGap = 0.065;
    if isfield(spec, 'rowGap')
        rowGap = spec.rowGap;
    end
    usableWidth = 1 - left - right - (spec.ncol - 1) * colGap;
    usableHeight = 1 - top - bottom - (spec.nrow - 1) * rowGap;
    if isfield(spec, 'colWidths')
        columnWeights = spec.colWidths;
    else
        columnWeights = ones(1, spec.ncol);
    end
    colWidths = usableWidth .* columnWeights ./ sum(columnWeights);
    rowHeight = usableHeight / spec.nrow;
    forceFill = false;
    if isfield(spec, 'forceFill')
        forceFill = spec.forceFill;
    end
    equalPanelHeight = false;
    if isfield(spec, 'equalPanelHeight')
        equalPanelHeight = spec.equalPanelHeight;
    end
    packAdjacent = ~equalPanelHeight && ~forceFill;
    if isfield(spec, 'packAdjacent')
        packAdjacent = spec.packAdjacent;
    end

    panelImages = cell(size(spec.panels));
    panelAspects = zeros(1, numel(spec.panels));
    for iPanel = 1:numel(spec.panels)
        panelImages{iPanel} = loadPanelImage(spec.panels{iPanel});
        panelAspects(iPanel) = size(panelImages{iPanel}, 2) / ...
            size(panelImages{iPanel}, 1);
    end
    if equalPanelHeight && ~isfield(spec, 'colWidths')
        columnWeights = zeros(1, spec.ncol);
        colCounts = zeros(1, spec.ncol);
        for iPanel = 1:numel(spec.panels)
            col = mod(spec.panels{iPanel}.tile - 1, spec.ncol) + 1;
            columnWeights(col) = columnWeights(col) + panelAspects(iPanel);
            colCounts(col) = colCounts(col) + 1;
        end
        columnWeights = columnWeights ./ max(1, colCounts);
        colWidths = usableWidth .* columnWeights ./ sum(columnWeights);
    end

    % Use the requested canvas size for layout math, not a screen-clamped
    % on-screen window (that mismatch was stretching force-filled panels).
    figWidth = spec.figSize(3);
    figHeight = spec.figSize(4);
    fig.Position = [spec.figSize(1:2) figWidth figHeight];
    rowTargetHeight = zeros(1, spec.nrow);
    if equalPanelHeight
        for iRow = 1:spec.nrow
            maxHeights = rowHeight;
            for iPanel = 1:numel(spec.panels)
                p = spec.panels{iPanel};
                row = ceil(p.tile / spec.ncol);
                if row ~= iRow
                    continue
                end
                col = mod(p.tile - 1, spec.ncol) + 1;
                maxHeights = min(maxHeights, ...
                    (colWidths(col) * figWidth) / ...
                    (panelAspects(iPanel) * figHeight));
            end
            rowTargetHeight(iRow) = maxHeights;
        end
    end

    if packAdjacent
        % Place panels in pixel units so width/height match image aspect
        % exactly, independent of any on-screen window resizing.
        titleBand = 0.028;
        if isfield(spec, 'titleBand')
            titleBand = spec.titleBand;
        end
        fixedColumnLayout = false;
        if isfield(spec, 'fixedColumnLayout')
            fixedColumnLayout = spec.fixedColumnLayout;
        end
        packWidth = 1 - left - right;
        rowPanels = cell(1, spec.nrow);
        for iPanel = 1:numel(spec.panels)
            row = ceil(spec.panels{iPanel}.tile / spec.ncol);
            rowPanels{row}(end + 1) = iPanel; %#ok<AGROW>
        end
        columnX = [];
        columnWidths = [];
        slotHeightFixed = [];
        if fixedColumnLayout && spec.ncol == 2
            colAspect = zeros(1, 2);
            colCount = zeros(1, 2);
            for iPanel = 1:numel(spec.panels)
                col = mod(spec.panels{iPanel}.tile - 1, spec.ncol) + 1;
                colAspect(col) = colAspect(col) + panelAspects(iPanel);
                colCount(col) = colCount(col) + 1;
            end
            colAspect = colAspect ./ max(1, colCount);
            heightBudget = max(rowHeight - titleBand, eps);
            widthBudget = packWidth - colGap;
            heightFromWidth = widthBudget / ...
                (sum(colAspect) * figHeight / figWidth);
            slotHeightFixed = min(heightBudget, heightFromWidth);
            columnWidths = slotHeightFixed * colAspect * figHeight / figWidth;
            rowWidthFixed = sum(columnWidths) + colGap;
            columnX = [left, left] + (packWidth - rowWidthFixed) / 2;
            columnX(2) = columnX(1) + columnWidths(1) + colGap;
        end
        for iRow = 1:spec.nrow
            idxs = rowPanels{iRow};
            if isempty(idxs)
                continue
            end
            [~, order] = sort(cellfun(@(p) p.tile, spec.panels(idxs)));
            idxs = idxs(order);
            aspects = panelAspects(idxs);
            nInRow = numel(idxs);
            if fixedColumnLayout && spec.ncol == 2
                slotHeight = slotHeightFixed;
                slotWidths = zeros(1, nInRow);
                xSlots = zeros(1, nInRow);
                for iInRow = 1:nInRow
                    col = mod(spec.panels{idxs(iInRow)}.tile - 1, spec.ncol) + 1;
                    slotWidths(iInRow) = columnWidths(col);
                    xSlots(iInRow) = columnX(col);
                end
            else
                heightBudget = max(rowHeight - titleBand, eps);
                widthBudget = packWidth - (nInRow - 1) * colGap;
                heightFromWidth = widthBudget / ...
                    (sum(aspects) * figHeight / figWidth);
                slotHeight = min(heightBudget, heightFromWidth);
                slotWidths = slotHeight * aspects * figHeight / figWidth;
                rowWidth = sum(slotWidths) + (nInRow - 1) * colGap;
                xSlots = zeros(1, nInRow);
                xSlots(1) = left + (packWidth - rowWidth) / 2;
                for iInRow = 2:nInRow
                    xSlots(iInRow) = xSlots(iInRow - 1) + ...
                        slotWidths(iInRow - 1) + colGap;
                end
            end
            yTile = 1 - top - iRow * rowHeight - (iRow - 1) * rowGap;
            ySlot = yTile + rowHeight - titleBand - slotHeight;
            for iInRow = 1:nInRow
                iPanel = idxs(iInRow);
                p = spec.panels{iPanel};
                displayScale = 1;
                if isfield(p, 'displayScale') && ~isempty(p.displayScale)
                    displayScale = p.displayScale;
                end
                imgHeight = slotHeight * displayScale;
                imgWidth = slotWidths(iInRow) * displayScale;
                xImg = xSlots(iInRow) + (slotWidths(iInRow) - imgWidth) / 2;
                % Keep top edge aligned with the full-size neighbour (F).
                yImg = ySlot + slotHeight - imgHeight;

                ax = axes(fig, 'Units', 'pixels', 'Position', ...
                    panelPixels(figWidth, figHeight, ...
                    xImg, yImg, imgWidth, imgHeight));
                renderPanel(ax, p, titleFontSize, false, ...
                    panelImages{iPanel}, false);

                panelTitleStr = panelTitleString(p);
                if ~isempty(panelTitleStr)
                    annotation(fig, 'textbox', ...
                        [xSlots(iInRow), ySlot + slotHeight, ...
                        slotWidths(iInRow), titleBand], ...
                        'String', panelTitleStr, ...
                        'Interpreter', 'tex', ...
                        'FontSize', titleFontSize, ...
                        'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'EdgeColor', 'none', ...
                        'BackgroundColor', 'none');
                end
            end
        end
    else
        for iPanel = 1:numel(spec.panels)
            p = spec.panels{iPanel};
            row = ceil(p.tile / spec.ncol);
            col = mod(p.tile - 1, spec.ncol) + 1;
            xTile = left + sum(colWidths(1:col - 1)) + (col - 1) * colGap;
            yTile = 1 - top - row * rowHeight - (row - 1) * rowGap;
            if equalPanelHeight
                panelHeight = rowTargetHeight(row);
                panelWidth = panelHeight * panelAspects(iPanel) * ...
                    figHeight / figWidth;
                xPosition = xTile + (colWidths(col) - panelWidth) / 2;
                yPosition = yTile + rowHeight - panelHeight;
                ax = axes(fig, 'Position', ...
                    [xPosition yPosition panelWidth panelHeight]);
                renderPanel(ax, p, titleFontSize, true, panelImages{iPanel});
            else
                ax = axes(fig, 'Position', ...
                    [xTile yTile colWidths(col) rowHeight]);
                renderPanel(ax, p, titleFontSize, forceFill, ...
                    panelImages{iPanel});
            end
        end
    end
elseif spec.ncol == 1 && spec.nrow == 1
    % Place the single panel explicitly. Using tiledlayout here left a large
    % unpredictable gap between the overall title and the image.
    p = spec.panels{1};
    img = loadPanelImage(p);
    imgAspect = size(img, 2) / size(img, 1);
    figWidth = spec.figSize(3);
    figHeight = spec.figSize(4);
    left = 0.03;
    right = 0.03;
    bottom = 0.03;
    panelTitleStr = panelTitleString(p);
    subtitleBand = 0;
    if ~isempty(panelTitleStr)
        subtitleBand = 0.05;
    end
    titleGap = 0.012;
    if isfield(spec, 'titleGap')
        titleGap = spec.titleGap;
    end
    usableWidth = 1 - left - right;
    contentTop = 1 - overallTitleBand - subtitleBand - titleGap;
    usableHeight = contentTop - bottom;
    heightFromWidth = usableWidth * figWidth / (imgAspect * figHeight);
    imgHeight = min(usableHeight, heightFromWidth);
    imgWidth = imgHeight * imgAspect * figHeight / figWidth;
    xImg = left + (usableWidth - imgWidth) / 2;
    yImg = contentTop - imgHeight;
    ax = axes(fig, 'Position', [xImg yImg imgWidth imgHeight]);
    renderPanel(ax, p, titleFontSize, true, img, false);
    if ~isempty(panelTitleStr)
        annotation(fig, 'textbox', ...
            [xImg, yImg + imgHeight + titleGap, imgWidth, subtitleBand], ...
            'String', panelTitleStr, ...
            'Interpreter', 'tex', ...
            'FontSize', titleFontSize, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FitBoxToText', 'off', ...
            'Margin', 0, ...
            'EdgeColor', 'none', ...
            'BackgroundColor', 'none');
    end
else
    tl = tiledlayout(fig, spec.nrow, spec.ncol, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    if hasOverallTitle
        tl.OuterPosition = [0, 0, 1, 1 - overallTitleBand];
    end
    for iPanel = 1:numel(spec.panels)
        p = spec.panels{iPanel};
        ax = nexttile(tl, p.tile);
        renderPanel(ax, p, titleFontSize, false);
    end
end

drawnow; pause(0.05);
% print renders the whole canvas. exportgraphics crops to what it considers
% content and at 600 dpi that dropped the overall title band on multi-row
% layouts, so the surrounding white is trimmed here instead.
fig.PaperPositionMode = 'auto';
fig.InvertHardcopy = 'off';
print(fig, outPng, '-dpng', sprintf('-r%d', exportDpi));
close(fig);
imwrite(trimWhiteBorders(imread(outPng), round(exportDpi / 12)), outPng);
end

function str = panelTitleString(p)
if ~isempty(p.letter) && ~isempty(p.title)
    str = sprintf('\\bf{%s} | %s', p.letter, p.title);
elseif ~isempty(p.letter)
    str = sprintf('\\bf{%s}', p.letter);
else
    str = p.title;
end
end

function img = loadPanelImage(p)
img = imread(p.file);
if p.trimWhite
    img = trimWhiteBorders(img);
end
if p.imageScale ~= 1
    img = scaleImageOnCanvas(img, p.imageScale);
end
end

function posPx = panelPixels(figWidth, figHeight, xNorm, yNorm, wNorm, hNorm)
% Convert normalized figure coordinates to pixel axes Position, snapping
% width from height so rounding cannot introduce aspect distortion.
x = max(1, round(xNorm * figWidth) + 1);
y = max(1, round(yNorm * figHeight) + 1);
h = max(1, round(hNorm * figHeight));
imgAspect = (wNorm / max(hNorm, eps)) * (figWidth / figHeight);
w = max(1, round(h * imgAspect));
posPx = [x, y, w, h];
end

function renderPanel(ax, p, titleFontSize, forceFill, img, showTitle)
if nargin < 5 || isempty(img)
    img = loadPanelImage(p);
end
if nargin < 6
    showTitle = true;
end

image(ax, img);
axis(ax, 'image');
axis(ax, 'off');
set(ax, 'Color', 'w', 'DataAspectRatio', [1 1 1]);
if showTitle
    panelTitleStr = panelTitleString(p);
    if ~isempty(panelTitleStr)
        title(ax, panelTitleStr, ...
            'Interpreter', 'tex', 'FontSize', titleFontSize, ...
            'FontWeight', 'bold');
    end
end
end

function imgOut = scaleImageOnCanvas(imgIn, imageScale)
scaledImage = imresize(imgIn, imageScale);
if isinteger(imgIn)
    whiteValue = intmax(class(imgIn));
else
    whiteValue = 1;
end
imgOut = repmat(cast(whiteValue, 'like', imgIn), size(imgIn));
rowStart = floor((size(imgIn, 1) - size(scaledImage, 1)) / 2) + 1;
colStart = floor((size(imgIn, 2) - size(scaledImage, 2)) / 2) + 1;
rowEnd = rowStart + size(scaledImage, 1) - 1;
colEnd = colStart + size(scaledImage, 2) - 1;
imgOut(rowStart:rowEnd, colStart:colEnd, :) = scaledImage;
end

function imgOut = trimWhiteBorders(imgIn, padding)
if nargin < 2 || isempty(padding)
    padding = 6;
end
imgOut = imgIn;
if isempty(imgIn)
    return
end

if size(imgIn, 3) == 1
    gray = imgIn;
else
    gray = rgb2gray(imgIn(:, :, 1:3));
end

mask = gray < 248;
if ~any(mask(:))
    return
end

[rows, cols] = find(mask);
r1 = max(1, min(rows) - padding);
r2 = min(size(imgIn, 1), max(rows) + padding);
c1 = max(1, min(cols) - padding);
c2 = min(size(imgIn, 2), max(cols) + padding);
imgOut = imgIn(r1:r2, c1:c2, :);
end
