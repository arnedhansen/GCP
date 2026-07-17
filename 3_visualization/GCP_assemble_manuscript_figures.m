%% GCP Assemble Manuscript Figures
% Combines exported PNG panels into manuscript figures.

%% Setup
startup
[~, paths, ~, ~] = setup('GCP', 0);

exportDpi = 600;
outDir = fullfile(paths.figures, 'manuscript');
if ~isfolder(outDir)
    mkdir(outDir);
end

powerDir = fullfile(paths.figures, 'power_analysis');
boxplotDir = fullfile(paths.figures, 'stats', 'boxplots');
gazeDir = fullfile(paths.figures, 'gaze');
eegDir = fullfile(paths.figures, 'eeg');

fprintf('Assembling GCP manuscript figures...\n');
fprintf('Output directory: %s\n', outDir);

%% Figure 1: Contrast detection paradigm
fig1 = figureSpec('Figure1', 1, 1, [0 0 1512 982]);
fig1.titleFontSize = 16;
fig1.panels = {
    panelSpec(fullfile(paths.figures, 'paradigm', 'GCP_paradigm.png'), ...
        'A', 'Grating Task Paradigm', 1, [1 1])
    };

%% Figure 2: Power analyses
% Frequency is shown in the first row and power in the second row.
% Heatmaps are placed on the left and power curves on the right.
fig2 = figureSpec('Figure2', 2, 2, [0 0 1512 982]);
fig2.colWidths = [2 1];
fig2.forceFill = false;
fig2.titleFontSize = 16;
fig2.panels = {
    panelSpec(fullfile(powerDir, 'GCP_power_analysis_frequency_heatmap.png'), ...
        'A', 'Gamma Frequency Statistical Power Heatmap', 1, [1 1]);
    panelSpec(fullfile(powerDir, 'GCP_power_analysis_frequency.png'), ...
        'B', 'Gamma Frequency Statistical Power Curve', 2, [1 1]);
    panelSpec(fullfile(powerDir, 'GCP_power_analysis_power_heatmap.png'), ...
        'C', 'Gamma Power Statistical Power Heatmap', 3, [1 1]);
    panelSpec(fullfile(powerDir, 'GCP_power_analysis_power.png'), ...
        'D', 'Gamma Power Statistical Power Curve', 4, [1 1])
    };
for iPanel = 1:numel(fig2.panels)
    fig2.panels{iPanel}.trimWhite = true;
end

%% Figure 3: Pupil size
fig3 = figureSpec('Figure3', 1, 2, [0 0 1512 982]);
fig3.panels = {
    panelSpec(fullfile(gazeDir, 'pupil', 'GCP_gaze_pupil_size_TC_db.png'), ...
        'A', 'Pupil Size Time Course', 1, [1 1]);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_dBPupilSize.png'), ...
        'B', 'Pupil Size', 2, [1 1])
    };
fig3.panels{2}.imageScale = 0.90;
for iPanel = 1:numel(fig3.panels)
    fig3.panels{iPanel}.trimWhite = true;
end

%% Figure 4: Microsaccade rate
fig4 = figureSpec('Figure4', 1, 2, [0 0 1512 982]);
fig4.panels = {
    panelSpec(fullfile(gazeDir, 'microsaccades', 'GCP_gaze_microsaccades_rate_db.png'), ...
        'A', 'Microsaccade Rate Time Course', 1, [1 1]);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_dBMSRate.png'), ...
        'B', 'Microsaccade Rate', 2, [1 1])
    };
fig4.panels{2}.imageScale = 0.90;
for iPanel = 1:numel(fig4.panels)
    fig4.panels{iPanel}.trimWhite = true;
end

%% Figure 5: Eye velocity
fig5 = figureSpec('Figure5', 1, 2, [0 0 1512 982]);
fig5.panels = {
    panelSpec(fullfile(gazeDir, 'velocity', 'GCP_gaze_velocity_Vel2D_TC_db.png'), ...
        'A', 'Combined Eye Velocity Time Course', 1, [1 1]);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_dBVel2D.png'), ...
        'B', 'Combined Eye Velocity', 2, [1 1])
    };
fig5.panels{2}.imageScale = 0.90;
for iPanel = 1:numel(fig5.panels)
    fig5.panels{iPanel}.trimWhite = true;
end

%% Figure 6: Bivariate contour ellipse area
fig6 = figureSpec('Figure6', 1, 2, [0 0 1512 982]);
fig6.panels = {
    panelSpec(fullfile(gazeDir, 'bcea', 'GCP_gaze_BCEA_TC_db.png'), ...
        'A', 'Gaze Dispersion Time Course', 1, [1 1]);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_dBBCEA.png'), ...
        'B', 'Gaze Dispersion', 2, [1 1])
    };
fig6.panels{2}.imageScale = 0.90;
for iPanel = 1:numel(fig6.panels)
    fig6.panels{iPanel}.trimWhite = true;
end

%% Gaze composite
figGazeComp = figureSpec('figGazeComp', 4, 2, ...
    [0 0 1512 round(982 * 2)]);
figGazeComp.titleFontSize = 16;
figGazeComp.panels = {
    panelSpec(fullfile(gazeDir, 'pupil', 'GCP_gaze_pupil_size_TC_db.png'), ...
        'A', 'Pupil Size Time Course', 1, [1 1]);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_dBPupilSize.png'), ...
        'B', 'Pupil Size', 2, [1 1]);
    panelSpec(fullfile(gazeDir, 'microsaccades', 'GCP_gaze_microsaccades_rate_db.png'), ...
        'C', 'Microsaccade Rate Time Course', 3, [1 1]);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_dBMSRate.png'), ...
        'D', 'Microsaccade Rate', 4, [1 1]);
    panelSpec(fullfile(gazeDir, 'velocity', 'GCP_gaze_velocity_Vel2D_TC_db.png'), ...
        'E', 'Combined Eye Velocity Time Course', 5, [1 1]);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_dBVel2D.png'), ...
        'F', 'Combined Eye Velocity', 6, [1 1]);
    panelSpec(fullfile(gazeDir, 'bcea', 'GCP_gaze_BCEA_TC_db.png'), ...
        'G', 'Gaze Dispersion Time Course', 7, [1 1]);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_dBBCEA.png'), ...
        'H', 'Gaze Dispersion', 8, [1 1])
    };
for iPanel = 2:2:numel(figGazeComp.panels)
    figGazeComp.panels{iPanel}.imageScale = 0.90;
end
for iPanel = 1:numel(figGazeComp.panels)
    figGazeComp.panels{iPanel}.trimWhite = true;
end

%% Figure 7: Gamma power and frequency
fig7 = figureSpec('Figure7', 2, 2, [0 0 1512 982]);
fig7.panels = {
    panelSpec(fullfile(eegDir, 'ersd', 'GCP_eeg_ersd_GED_timecourse.png'), ...
        'A', 'Gamma ERSD Time Course', 1, [1 1]);
    panelSpec(fullfile(eegDir, 'topos', 'GCP_eeg_topos_ga.png'), ...
        'B', 'Gamma Power Topographies', 2, [1 1]);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_Power.png'), ...
        'C', 'Gamma Power', 3, [1 1]);
    panelSpec(fullfile(boxplotDir, 'GCP_stats_boxplot_Frequency.png'), ...
        'D', 'Gamma Peak Frequency', 4, [1 1])
    };
fig7.panels{2}.imageScale = 0.92;
fig7.panels{3}.imageScale = 0.88;
fig7.panels{4}.imageScale = 0.88;
for iPanel = 1:numel(fig7.panels)
    fig7.panels{iPanel}.trimWhite = true;
end

%% Figure 8: Gamma GED time frequency representation
fig8 = figureSpec('Figure8', 1, 1, [0 0 1512 982]);
fig8.panels = {
    panelSpec(fullfile(eegDir, 'tfr', 'GCP_eeg_tfr_GED.png'), ...
        'A', 'Time-Frequency Representations', 1, [1 1])
    };

%% Supplementary Figure S1: Single participant power spectra
figS1 = figureSpec('FigureS1', 1, 1, [0 0 1512 982]);
figS1.titleFontSize = 16;
figS1.panels = {
    panelSpec(fullfile(eegDir, 'powspctrm', ...
        'GCP_eeg_GED_powspctrm_overview_subjects.png'), ...
        'S1', 'Single Subject Powerspectra', 1, [1 1])
    };

%% Assemble all figures
figSpecs = {fig1, fig2, fig3, fig4, fig5, fig6, ...
    figGazeComp, fig7, fig8, figS1};
for iFig = 1:numel(figSpecs)
    spec = figSpecs{iFig};
    outPng = fullfile(outDir, ['GCP_manuscript_' spec.name '.png']);
    assembleManuscriptFigure(spec, outPng, exportDpi);
    fprintf('Saved %s\n', outPng);
end

fprintf('\nDone. Manuscript composites saved to:\n  %s\n', outDir);

%% Local functions
function spec = figureSpec(name, nrow, ncol, figSize)
spec = struct();
spec.name = name;
spec.nrow = nrow;
spec.ncol = ncol;
spec.figSize = figSize;
spec.panels = {};
end

function p = panelSpec(imagePath, letter, titleText, tile, span, trimWhite)
p = struct( ...
    'file', imagePath, ...
    'letter', letter, ...
    'title', titleText, ...
    'tile', tile, ...
    'span', span, ...
    'trimWhite', false, ...
    'imageScale', 1);
if nargin >= 6 && ~isempty(trimWhite)
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
    fprintf('\nMissing panel files for %s:\n', spec.name);
    for iMissing = 1:numel(missing)
        fprintf('  %s\n', missing{iMissing});
    end
    error('GCP_assemble_manuscript_figures:MissingPanels', ...
        'Cannot assemble %s until all panel PNG files exist.', spec.name);
end

fig = figure('Position', spec.figSize, 'Color', 'w');
titleFontSize = 20;
if isfield(spec, 'titleFontSize')
    titleFontSize = spec.titleFontSize;
end

if isfield(spec, 'colWidths') || spec.ncol == 2
    left = 0.025;
    right = 0.025;
    top = 0.035;
    bottom = 0.025;
    colGap = 0.010;
    rowGap = 0.065;
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

    for iPanel = 1:numel(spec.panels)
        p = spec.panels{iPanel};
        row = ceil(p.tile / spec.ncol);
        col = mod(p.tile - 1, spec.ncol) + 1;
        xPosition = left + sum(colWidths(1:col - 1)) + (col - 1) * colGap;
        yPosition = 1 - top - row * rowHeight - (row - 1) * rowGap;
        ax = axes(fig, 'Position', ...
            [xPosition yPosition colWidths(col) rowHeight]);
        renderPanel(ax, p, titleFontSize, forceFill);
    end
else
    tl = tiledlayout(fig, spec.nrow, spec.ncol, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    for iPanel = 1:numel(spec.panels)
        p = spec.panels{iPanel};
        ax = nexttile(tl, p.tile, p.span);
        renderPanel(ax, p, titleFontSize, false);
    end
end

drawnow;
exportgraphics(fig, outPng, 'Resolution', exportDpi, ...
    'BackgroundColor', 'white');
close(fig);
end

function renderPanel(ax, p, titleFontSize, forceFill)
img = imread(p.file);
if p.trimWhite
    img = trimWhiteBorders(img);
end
if p.imageScale ~= 1
    img = scaleImageOnCanvas(img, p.imageScale);
end

image(ax, img);
if forceFill
    axis(ax, 'normal');
else
    axis(ax, 'image');
end
axis(ax, 'off');
set(ax, 'Color', 'w');
title(ax, sprintf('\\bf{%s} | %s', p.letter, p.title), ...
    'Interpreter', 'tex', 'FontSize', titleFontSize, ...
    'FontWeight', 'bold');
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

function imgOut = trimWhiteBorders(imgIn)
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
padding = 6;
r1 = max(1, min(rows) - padding);
r2 = min(size(imgIn, 1), max(rows) + padding);
c1 = max(1, min(cols) - padding);
c2 = min(size(imgIn, 2), max(cols) + padding);
imgOut = imgIn(r1:r2, c1:c2, :);
end
