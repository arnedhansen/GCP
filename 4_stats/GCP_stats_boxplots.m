%% GCP Stats Overview
% One figure per variable: boxplots + subject lines + jittered dots
clear
clc
close all
[~, ~, colors, ~] = setup('GCP', 0);

% Load and convert
load('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/merged_data.mat')
T = struct2table(merged_data);

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

%% Loop through variables
for iVar = 22%%%%%1:numel(numericVars)
    close all

    varName = numericVars{iVar};
    y       = T.(varName);

    % New figure for this variable
    figure;
    set(gcf, 'Position', [200 200 900 700]);
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
    boxplot(y, condIdx, ...
        'Positions', 1:nCond, ...
        'Symbol', '');  % no outlier symbols
    set(findobj(gca, 'Tag', 'Box'), 'LineWidth', 1.5);

    % Overlay jittered dots + mean/range markers per condition
    jitterWidth = 0.15;

    for c = 1:nCond
        idxC = condIdx == c;
        yC   = y(idxC);

        % Jittered x positions
        xJit = c + (rand(sum(idxC), 1) - 0.5) * 2 * jitterWidth;

        scatter(xJit, yC, 125, colors(c, :), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

        % Vertical range line (min–max)
        if ~isempty(yC)
            yMin = min(yC);
            yMax = max(yC);
            plot([c c], [yMin yMax], ':', ...
                'Color', [0.5 0.5 0.5], 'LineWidth', 1);

            % Mean marker
            plot(c, mean(yC), 'o', ...
                'MarkerFaceColor', colors(c, :), ...
                'MarkerEdgeColor', 'k', ...
                'MarkerSize', 7, ...
                'LineWidth', 1);
        end
    end

    % Zero line
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

    % Axes formatting
    xlim([0.5 nCond + 0.5]);
    xticks(1:nCond);
    xticklabels(xtickLabs);

    % Titles and labels (pretty-print variable name)
    prettyName = strrep(varName, '_', ' ');
    title(prettyName, 'FontSize', 20, 'FontWeight', 'bold');
    xlabel('Contrast');
    ylabel(prettyName, 'Interpreter', 'none');

    set(gca, 'FontSize', 12, 'Box', 'off');
    saveas(gcf, ['/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/stats/boxplots/GCP_stats_boxplot_', varName, '.png'])
end
