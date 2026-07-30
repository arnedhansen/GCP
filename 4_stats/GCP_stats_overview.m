%% GCP Stats Overview
% Quick overview plots: boxplots per variable across conditions
[~, paths] = setup('GCP', 0);

% Load and convert
load(fullfile(paths.features, 'GCP_merged_data.mat'))
T = struct2table(merged_data);

if ismember('Include', T.Properties.VariableNames)
    T = T(T.Include, :);
end

% Keep all numeric variables, but exclude non baselined gaze measures
var_names = T.Properties.VariableNames;
gaze_roots = {'MSRate', 'BCEA', 'Vel2D', 'VelV', 'PupilSize', ...
    'Blinks', 'Fixations', 'Saccades'};
numeric_vars = {};

for i = 1:numel(var_names)
    vn = var_names{i};
    v = T.(vn);
    if ~isnumeric(v)
        continue
    end
    if strcmp(vn, 'ID') || strcmp(vn, 'Condition') || strcmp(vn, 'Include')
        continue
    end

    vn_core = regexprep(vn, '^Gaze_', '');
    is_gaze_metric = false;
    for gi = 1:numel(gaze_roots)
        if startsWith(vn_core, gaze_roots{gi})
            is_gaze_metric = true;
            break
        end
    end

    if is_gaze_metric
        % Keep gaze only when baselined (full, early, or late windows)
        if contains(vn, '_bl')
            numeric_vars{end+1} = vn; %#ok<AGROW>
        end
    else
        numeric_vars{end+1} = vn; %#ok<AGROW>
    end
end

if isempty(numeric_vars)
    error('GCP_stats_overview:NoVariablesSelected', ...
        'No variables available after filtering.');
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

sgtitle('All Variables with Baselined Gaze Measures: Boxplots per Condition');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(paths.figures, 'stats', 'overview', 'GCP_stats_overview.png'), '-dpng', '-r600');