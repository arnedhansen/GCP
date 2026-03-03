%% GCP Trial-Level Mixed Models
%
% Dedicated trial-level analysis script using the merged trial table.
% This script intentionally does not depend on GCP_hypotheses_trials.m.
%
% Input:
%   - /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_merged_data_trials.mat
%
% Output:
%   - /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_trial_level_model_results.mat
%   - /Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_trial_level_model_results.csv

%% Setup
clear
clc
close all

features_root = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features';
merged_path = fullfile(features_root, 'GCP_merged_data_trials.mat');

if ~isfile(merged_path)
    error('Merged trial-level file not found: %s', merged_path);
end

S = load(merged_path);
if isfield(S, 'GCP_merged_table_trials')
    T = S.GCP_merged_table_trials;
elseif isfield(S, 'GCP_merged_data_trials')
    T = struct2table(S.GCP_merged_data_trials);
else
    error('No expected merged variable found in %s', merged_path);
end

required_keys = {'ID','Condition','Trial'};
if ~all(ismember(required_keys, T.Properties.VariableNames))
    error('Merged table is missing one or more key columns: ID, Condition, Trial');
end

T.ID = categorical(T.ID);
T.Condition = categorical(T.Condition, [1 2 3 4], {'25','50','75','100'}, 'Ordinal', true);

results = table('Size', [0, 6], ...
    'VariableTypes', {'string','double','double','double','double','double'}, ...
    'VariableNames', {'Model','N','NumSubjects','BetaCondition','PCondition','Converged'});

%% Behavioral trial-level analyses
results = [results; run_behavior_models(T)]; %#ok<AGROW>

%% Gaze trial-level analyses
results = [results; run_gaze_models(T)]; %#ok<AGROW>

%% GED trial-level analyses (if available)
results = [results; run_ged_models(T)]; %#ok<AGROW>

%% Save
GCP_trial_level_model_results = results; %#ok<NASGU>
save(fullfile(features_root, 'GCP_trial_level_model_results.mat'), ...
    'GCP_trial_level_model_results');
writetable(results, fullfile(features_root, 'GCP_trial_level_model_results.csv'));

disp('Trial-level model summary:')
disp(results)

%% Local functions
function out = run_behavior_models(T)
out = empty_result_table();

% Accuracy (binomial GLME)
if has_vars(T, {'Behavior_Accuracy'})
    Tb = T(~isnan(T.Behavior_Accuracy), :);
    if ~isempty(Tb)
        Tb.Behavior_Accuracy = double(Tb.Behavior_Accuracy ~= 0);
        model_name = "Behavior_Accuracy_binomial";
        out = [out; fit_glme_row(Tb, ...
            'Behavior_Accuracy ~ Condition + (1|ID)', ...
            'Binomial', 'Logit', model_name)]; %#ok<AGROW>
    end
end

% Reaction time (log Gaussian LME)
if has_vars(T, {'Behavior_ReactionTime'})
    Tb = T(~isnan(T.Behavior_ReactionTime) & T.Behavior_ReactionTime > 0, :);
    if ~isempty(Tb)
        Tb.LogRT = log(double(Tb.Behavior_ReactionTime));
        model_name = "Behavior_LogRT_gaussian";
        out = [out; fit_lme_row(Tb, ...
            'LogRT ~ Condition + (1|ID)', ...
            model_name)]; %#ok<AGROW>
    end
end
end

function out = run_gaze_models(T)
out = empty_result_table();
gaze_specs = { ...
    'Gaze_PctMSRate', 'Gaze_PctMSRate_lme'; ...
    'Gaze_PctVel2D', 'Gaze_PctVel2D_lme'; ...
    'Gaze_PctPupilSize', 'Gaze_PctPupilSize_lme'; ...
    'Gaze_PctBCEA', 'Gaze_PctBCEA_lme'};

for i = 1:size(gaze_specs, 1)
    varname = gaze_specs{i, 1};
    model_name = string(gaze_specs{i, 2});
    if ~has_vars(T, {varname})
        continue
    end
    Tg = T(~isnan(T.(varname)), :);
    if isempty(Tg)
        continue
    end
    formula = sprintf('%s ~ Condition + (1|ID)', varname);
    out = [out; fit_lme_row(Tg, formula, model_name)]; %#ok<AGROW>
end
end

function out = run_ged_models(T)
out = empty_result_table();
ged_specs = { ...
    'GED_PeakFrequency', 'GED_PeakFrequency_lme'; ...
    'GED_PeakAmplitude', 'GED_PeakAmplitude_lme'; ...
    'GED_BroadbandPower', 'GED_BroadbandPower_lme'};

for i = 1:size(ged_specs, 1)
    varname = ged_specs{i, 1};
    model_name = string(ged_specs{i, 2});
    if ~has_vars(T, {varname})
        continue
    end
    Te = T(~isnan(T.(varname)), :);
    if isempty(Te)
        continue
    end
    formula = sprintf('%s ~ Condition + (1|ID)', varname);
    out = [out; fit_lme_row(Te, formula, model_name)]; %#ok<AGROW>
end
end

function row = fit_lme_row(T, formula, model_name)
try
    mdl = fitlme(T, formula, 'FitMethod', 'REML');
    [beta, p] = condition_beta_and_p(mdl.Coefficients);
    n_subj = numel(unique(T.ID));
    row = make_row(model_name, height(T), n_subj, beta, p, 1);
catch
    n_subj = numel(unique(T.ID));
    row = make_row(model_name, height(T), n_subj, NaN, NaN, 0);
end
end

function row = fit_glme_row(T, formula, dist, link, model_name)
try
    mdl = fitglme(T, formula, 'Distribution', dist, 'Link', link, 'FitMethod', 'Laplace');
    [beta, p] = condition_beta_and_p(mdl.Coefficients);
    n_subj = numel(unique(T.ID));
    row = make_row(model_name, height(T), n_subj, beta, p, 1);
catch
    n_subj = numel(unique(T.ID));
    row = make_row(model_name, height(T), n_subj, NaN, NaN, 0);
end
end

function [beta, p] = condition_beta_and_p(coef_tbl)
coef_names = string(coef_tbl.Name);
is_cond = startsWith(coef_names, "Condition_");
if any(is_cond)
    beta = coef_tbl.Estimate(find(is_cond, 1, 'first'));
    p = coef_tbl.pValue(find(is_cond, 1, 'first'));
else
    beta = NaN;
    p = NaN;
end
end

function tf = has_vars(T, vars)
tf = all(ismember(vars, T.Properties.VariableNames));
end

function row = make_row(name, n, n_subj, beta, p, conv)
row = table(string(name), double(n), double(n_subj), double(beta), double(p), double(conv), ...
    'VariableNames', {'Model','N','NumSubjects','BetaCondition','PCondition','Converged'});
end

function T = empty_result_table()
T = table('Size', [0, 6], ...
    'VariableTypes', {'string','double','double','double','double','double'}, ...
    'VariableNames', {'Model','N','NumSubjects','BetaCondition','PCondition','Converged'});
end
