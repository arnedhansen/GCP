%% GCP Stats Overview
% Quick overview plots: boxplots per variable across conditions
[~, paths] = setup('GCP', 0);

% Load and convert
load(fullfile(paths.features, 'GCP_merged_data.mat'))
T = struct2table(merged_data);

if ismember('Include', T.Properties.VariableNames)
    T = T(T.Include == 1, :);
end

% Keep all numeric variables, but exclude non baselined gaze measures.
% Velocity: Vel2D only (drop axis components VelH/VelV and Baseline* scalars).
var_names = T.Properties.VariableNames;
gaze_roots = {'MSRate', 'BCEA', 'Vel2D', 'PupilSize', ...
    'Blinks', 'Fixations', 'Saccades'};
drop_roots = {'VelH', 'VelV', 'Baseline'};
numeric_vars = {};

for i = 1:numel(var_names)
    vn = var_names{i};
    v = T.(vn);
    if ~isnumeric(v)
        continue
    end
    if strcmp(vn, 'ID') || strcmp(vn, 'Condition') || strcmp(vn, 'Include') || strcmp(vn, 'Age')
        continue
    end

    vn_core = regexprep(vn, '^Gaze_', '');

    drop_var = false;
    for di = 1:numel(drop_roots)
        if startsWith(vn_core, drop_roots{di})
            drop_var = true;
            break
        end
    end
    if drop_var
        continue
    end

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

set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(paths.figures, 'stats', 'overview', 'GCP_stats_overview.png'), '-dpng', '-r600');