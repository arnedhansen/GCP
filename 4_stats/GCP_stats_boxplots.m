%% GCP Stats Overview
% One figure per variable: boxplots + subject lines + jittered dots
startup
[~, paths, colors, ~] = setup('GCP', 0);

% Load merged data
dat = load(fullfile(paths.features, 'GCP_merged_data.mat'));
tbl = struct2table(dat.GCP_merged_data);

if ismember('Include', tbl.Properties.VariableNames)
    nRowsBefore = height(tbl);
    tbl = tbl(tbl.Include, :);
    fprintf('GED cohort filter: %d rows kept (%d excluded).\n', ...
        height(tbl), nRowsBefore - height(tbl));
end

% Identify numeric variables only (excluding ID and Condition)
varNames    = tbl.Properties.VariableNames;
numericVars = {};
for i = 1:numel(varNames)
    v = tbl.(varNames{i});
    if isnumeric(v) && ~strcmp(varNames{i}, 'ID') && ~strcmp(varNames{i}, 'Condition') ...
            && ~strcmp(varNames{i}, 'Include')
        numericVars{end+1} = varNames{i};
    end
end

% Condition coding
condRaw = tbl.Condition;
[condLevels, ~, condIdx] = unique(condRaw, 'stable');  % 1..K
nCond   = numel(condLevels);
xtickLabs = strcat(num2str(condLevels*25), "% Contrast");

% Subject IDs
subjIDs = unique(tbl.ID);

% Outlier exclusion settings (applied per variable and per condition)
useOutlierExclusion = false;
madThreshold = 5; % robust z-threshold based on MAD

% Plot aesthetics
fontSize        = 50;
xTickFontSize   = fontSize *0.9;
yTickFontSize   = fontSize *0.9;
ylabelFontSize  = fontSize *0.9;
titleFontSize   = fontSize;
legendFontSize  = fontSize *0.65;
dotSize         = 400;
dotAlpha        = 0.85;
jitter          = 0.4;
boxWidth        = 0.55;
showXTickLabels = true;
showLegend      = false;
showTitle       = false;

%% Loop through variables
for iVar = 1:numel(numericVars)
    close all

    varName = numericVars{iVar};
    yRaw    = tbl.(varName);
    y       = yRaw;

    if useOutlierExclusion
        isOutlier = false(size(yRaw));
        for c = 1:nCond
            idxC = condIdx == c;
            yC = yRaw(idxC);

            medC = median(yC, 'omitnan');
            madC = median(abs(yC - medC), 'omitnan');

            % If MAD is undefined/zero, robust z-scores cannot be computed.
            if isnan(madC) || madC == 0
                continue
            end

            robustZ = abs(yC - medC) ./ (1.4826 * madC);
            outC = robustZ > madThreshold;
            outC(isnan(robustZ)) = false;

            idxAll = find(idxC);
            isOutlier(idxAll(outC)) = true;
        end

        y(isOutlier) = NaN;
    end

    % New figure for this variable
    figure('Position', [0 0 1512 982], 'Color', 'w');
    hold on

    % One jittered x position per subject per condition (shared by lines and dots)
    xJit = nan(size(y));
    for s = 1:numel(subjIDs)
        thisID  = subjIDs(s);
        idxSubj = tbl.ID == thisID;
        for c = 1:nCond
            idxPt = idxSubj & condIdx == c & ~isnan(y);
            if any(idxPt)
                xJit(idxPt) = c + jitter * (rand - 0.5);
            end
        end
    end

    % Boxplot per condition (separate calls avoid grouped-box misalignment)
    for c = 1:nCond
        idxC = condIdx == c & ~isnan(y);
        yC   = y(idxC);
        if isempty(yC)
            continue
        end
        boxplot(yC, ones(numel(yC), 1), 'Positions', c, 'Symbol', '', ...
            'Widths', boxWidth, 'Colors', 'k');
    end
    styleCurrentBoxplot(colors(1:nCond, :));

    % Zero line
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

    % Subject-wise lines across conditions (connect jittered dot positions)
    for s = 1:numel(subjIDs)
        thisID  = subjIDs(s);
        idxSubj = tbl.ID == thisID;

        xSubj = xJit(idxSubj);
        ySubj = y(idxSubj);
        condSubjIdx = condIdx(idxSubj);

        valid = ~isnan(ySubj) & ~isnan(xSubj);
        xSubj = xSubj(valid);
        ySubj = ySubj(valid);
        condSubjIdx = condSubjIdx(valid);

        if numel(ySubj) < 2
            continue
        end

        [condSubjIdx, sortIdx] = sort(condSubjIdx);
        xSubj = xSubj(sortIdx);
        ySubj = ySubj(sortIdx);

        plot(xSubj, ySubj, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    end

    % Jittered dots on top (foreground)
    for c = 1:nCond
        idxC = condIdx == c & ~isnan(y);
        scatter(xJit(idxC), y(idxC), dotSize, colors(c, :), 'filled', ...
            'MarkerFaceAlpha', dotAlpha);
    end

    % Axes formatting
    xlabel('');
    yMin = min(y, [], 'omitnan');
    yMax = max(y, [], 'omitnan');
    yPadFrac = 0.1;
    if isempty(yMin) || isnan(yMin) || isempty(yMax) || isnan(yMax)
        ylim([0 1]);
    else
        yRange = yMax - yMin;
        if yRange == 0
            yRange = max(abs(yMax), abs(yMin), 1) * 0.1;
        end
        yPad = yRange * yPadFrac;
        ylim([yMin - yPad, yMax + yPad]);
    end
    xlim([0.5 nCond + 0.5]);
    xticks(1:nCond);
    if showXTickLabels
        xticklabels(xtickLabs);
    else
        xticklabels({});
    end
    prettyName = strrep(varName, '_', ' ');
    ylabStr = prettyName;
    titleStr = prettyName;
    if strcmp(varName, 'dBMSRate')
        ylabStr = 'Microsaccade Rate [%]';
        titleStr = 'Microsaccade Rate';
    elseif strcmp(varName, 'dBPupilSize')
        ylabStr = 'Pupil Size [%]';
        titleStr = 'Pupil Size';
    elseif strcmp(varName, 'Power')
        ylabStr = 'Power [dB]';
        titleStr = 'Gamma Power';
        ylim([1 4.5])
    elseif strcmp(varName, 'dBGazeDeviation')
        ylabStr = 'Gaze Deviation [dB]';
        titleStr = 'Gaze Deviation';
    elseif strcmp(varName, 'dBVel2D')
        ylabStr = 'Eye Velocity [%]';
        titleStr = 'Combined Eye Velocity';
    elseif strcmp(varName, 'dBBCEA')
        ylabStr = 'BCEA [%]';
        titleStr = 'BCEA';
    elseif strcmp(varName, 'dBFixations')
        ylabStr = 'Fixations [dB]';
        titleStr = 'Fixations';
    elseif strcmp(varName, 'dBBlinks')
        ylabStr = 'Blinks [dB]';
        titleStr = 'Blinks';
    elseif strcmp(varName, 'dBSaccades')
        ylabStr = 'Saccades [dB]';
        titleStr = 'Saccades';
    elseif strcmp(varName, 'Frequency')
        ylabStr = 'Frequency [Hz]';
        titleStr = 'Gamma Frequency';
    end

    hYlab = ylabel(ylabStr, 'Interpreter', 'none', 'FontSize', ylabelFontSize);
    if showTitle
        hTitle = title(titleStr, 'FontSize', titleFontSize, 'FontWeight', 'bold');
    else
        hTitle = title('');
    end

    if showLegend
        hLeg = gobjects(nCond, 1);
        for c = 1:nCond
            hLeg(c) = patch(nan, nan, colors(c, :), 'FaceAlpha', 0.25, ...
                'EdgeColor', colors(c, :), 'LineWidth', 1.5);
        end
        legend(hLeg, xtickLabs, 'FontSize', legendFontSize, 'Location', 'northeast', 'Box', 'off');
    end

    ax = gca;
    ax.XAxis.FontSize = xTickFontSize;
    ax.YAxis.FontSize = yTickFontSize;
    hYlab.FontSize = ylabelFontSize;
    if showTitle
        hTitle.FontSize = titleFontSize;
    end
    set(ax, 'Box', 'off');
    box off;
    hold off;

    % Save
    drawnow;
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, fullfile(paths.figures, 'stats', 'boxplots', ['GCP_stats_boxplot_' varName '.png']), '-dpng', '-r600');
end

%%
function styleCurrentBoxplot(boxColors)
% Color box faces and edges to match condition scatter colors
hBoxes = findobj(gca, 'Tag', 'Box');
if isempty(hBoxes)
    return;
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
        boxColors(bi, :), 'FaceAlpha', 0.25, 'EdgeColor', boxColors(bi, :), 'LineWidth', 1.5);
end
end