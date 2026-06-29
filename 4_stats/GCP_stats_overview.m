%% GCP Stats Overview
% Quick overview plots: boxplots per variable across conditions
[~, paths] = setup('GCP', 0);

% Load and convert
load(fullfile(paths.features, 'GCP_merged_data.mat'))
T = struct2table(merged_data);

if ismember('Include', T.Properties.VariableNames)
    T = T(T.Include, :);
end

% Identify numeric variables only (excluding ID and Condition)
var_names = T.Properties.VariableNames;
numeric_vars = {};

for i = 1:numel(var_names)
    v = T.(var_names{i});
    if isnumeric(v) && ~strcmp(var_names{i}, 'ID') && ~strcmp(var_names{i}, 'Condition') ...
            && ~strcmp(var_names{i}, 'Include')
        numeric_vars{end+1} = var_names{i};
    end
end

% Determine subplot grid
nVars = numel(numeric_vars);
nCols = ceil(sqrt(nVars));
nRows = ceil(nVars / nCols);

% Prepare figure
figure('Position', [0 0 1512 982], 'Color', 'w');

% Loop through variables
for i = 1:nVars
    subplot(nRows, nCols, i);
    
    % Get current variable
    var = numeric_vars{i};
    y = T.(var);
    cond = T.Condition;
    
    % Boxplot
    boxplot(y, cond, 'Symbol', '');
    title(var, 'Interpreter', 'none');
    xlabel('Condition');
    ylabel(var, 'Interpreter', 'none');
    
    grid on;
end

sgtitle('Merged Data Overview: Boxplots per Condition');
exportgraphics(gcf, fullfile(paths.figures, 'stats', 'overview', 'GCP_stats_overview.png'), 'Resolution', 600);