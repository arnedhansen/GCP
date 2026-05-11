%% GCP Stats Overview
% One figure per variable: boxplots + subject lines + jittered dots
clear
[~, paths, colors, ~] = setup('GCP', 0);

% Load and convert
S = load(fullfile(paths.features, 'GCP_merged_data.mat'));
if isfield(S, 'GCP_merged_data')
    T = struct2table(S.GCP_merged_data);
elseif isfield(S, 'merged_data')
    T = struct2table(S.merged_data);
else
    error('GCP_merged_data.mat does not contain GCP_merged_data or merged_data.');
end

% Identify numeric variables only (excluding ID and Condition)
varNames    = T.Properties.VariableNames;
numericVars = {};

for i = 1:numel(varNames)
    v = T.(varNames{i});
    if isnumeric(v) && ~strcmp(varNames{i}, 'ID') && ~strcmp(varNames{i}, 'Condition')
        numericVars{end+1} = varNames{i};
    end
end

% Condition coding
condRaw = T.Condition;
[condLevels, ~, condIdx] = unique(condRaw, 'stable');  % 1..K
nCond   = numel(condLevels);

% Nice x-tick labels for contrast (change here if Condition is something else)
xtickLabs = strcat(num2str(condLevels*25), "% Contrast");

% Subject IDs
subjIDs = unique(T.ID);

% Outlier exclusion settings (applied per variable and per condition)
useOutlierExclusion = true;
madThreshold = 3.5; % robust z-threshold based on MAD

%% Loop through variables
for iVar = 35%%%%%%1:numel(numericVars)
    close all

    varName = numericVars{iVar};
    yRaw    = T.(varName);
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

    % First: subject-wise lines across conditions (light grey)
    for s = 1:numel(subjIDs)
        thisID  = subjIDs(s);
        idxSubj = T.ID == thisID;

        ySubj       = y(idxSubj);
        condSubjIdx = condIdx(idxSubj); % 1..nCond per row

        % sort within subject by condition order
        [condSubjIdx, sortIdx] = sort(condSubjIdx);
        ySubj = ySubj(sortIdx);

        % skip if fewer than 2 conditions for this subject
        if numel(ySubj) < 2
            continue
        end

        plot(condSubjIdx, ySubj, '-', ...
            'Color', [0.8 0.8 0.8], ...
            'LineWidth', 1);
    end

    % Boxplot per condition (using condIdx so positions are 1..nCond)
    boxplot(y, condIdx, 'Colors', 'k', 'Symbol', '');

    hold on
    % Overlay jittered dots + mean/range markers per condition
    for c = 1:nCond
        idxC = condIdx == c & ~isnan(y);
        yC   = y(idxC);

        xJit = c + (rand(sum(idxC),1)-0.5)*0.1;
        scatter(xJit, yC, 250, colors(c,:), 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth',0.5)
    end

    % Zero line
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

    % Axes formatting
    xlabel('');
    maxval = max(abs(y), [], 'omitnan');
    if isempty(maxval) || isnan(maxval)
        maxval = 1;
    end
    if contains(varName, 'Pct')
        if maxval == 0
            ylim([-1 1])
        else
            ylim([-maxval*1.1 maxval*1.1])
        end
    elseif strcmp(varName, 'Frequency')
        ylim([30 90])
    % elseif strcmp(varName, 'Power')
        %set(gca, 'YScale', 'log')
    else
        if maxval == 0
            ylim([0 1])
        else
            ylim([0 maxval*1.1])
        end
    end
    xlim([0.5 nCond + 0.5]);
    xticks(1:nCond);
    xticklabels(xtickLabs);
    prettyName = strrep(varName, '_', ' ');
    ylabel(prettyName, 'Interpreter', 'none');
    set(gca, 'FontSize', 20, 'Box', 'off');
    title(prettyName, 'FontSize', 30, 'FontWeight', 'bold');
    if strcmp(varName, 'PctMSRate')
        ylabel('Microsaccade Rate [%]')
        title('Microsaccade Rate', 'FontSize', 30, 'FontWeight', 'bold')
    elseif strcmp(varName, 'PctPupilSize')
        ylabel('Pupil Size [%]')
        title('Pupil Size', 'FontSize', 30, 'FontWeight', 'bold')
    elseif strcmp(varName, 'Power')
        ylabel('Power [dB]')
        title('Gamma Power', 'FontSize', 30, 'FontWeight', 'bold')
    elseif strcmp(varName, 'PctGazeDeviation')
        ylabel('Gaze Deviation [%]')
        title('Gaze Deviation', 'FontSize', 30, 'FontWeight', 'bold')
    elseif strcmp(varName, 'PctVel2D')
        ylabel('Eye Velocity [%]')
        title('Combined Eye Velocity', 'FontSize', 30, 'FontWeight', 'bold')
    elseif strcmp(varName, 'PctFixations')
        ylabel('Fixations [%]')
        title('Fixations', 'FontSize', 30, 'FontWeight', 'bold')
    elseif strcmp(varName, 'PctBlinks')
        ylabel('Blinks [%]')
        title('Blinks', 'FontSize', 30, 'FontWeight', 'bold')
    elseif strcmp(varName, 'PctSaccades')
        ylabel('Saccades [%]')
        title('Saccades', 'FontSize', 30, 'FontWeight', 'bold')
    end
    % Save
    drawnow;
    exportgraphics(gcf, fullfile(paths.figures, 'stats', 'boxplots', ['GCP_stats_boxplot_' varName '.png']), 'Resolution', 600);
end