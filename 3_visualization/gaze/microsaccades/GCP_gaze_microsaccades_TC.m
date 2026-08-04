%% GCP Gaze Microsaccade Suppression Contrast Conditions
% Loads already processed microsaccade percentage-change time courses from
% 2_feature_extraction/GCP_gaze_fex.m (msTS_cXX_bl) and plots
% per-condition grand-average traces with SEM shading.
%
% Important:
%   - Trials with MS rate < 0.1 Hz in the baseline or stimulus window are excluded.
%   - Remaining trials use per-trial % change, then are averaged for the TC.
%   - Subject boxplot scalars are the mean of trial-level % values per condition.
%   - A light additional display smoothing is applied here.

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar')

datapath = paths.features;
figpath  = fullfile(paths.figures, 'gaze', 'microsaccades');
mkdir(figpath);
nSubj = length(subjects);

%% Display-only smoothing
fsample = 500;                    % Hz
sigma_ms_disp = 20;               % light display smoothing
sigma_samp = max(1, round(sigma_ms_disp / (1000 / fsample)));
kHalf = 3 * sigma_samp;
x_kern = -kHalf:kHalf;
gKernel = exp(-x_kern.^2 / (2 * sigma_samp^2));
gKernel = gKernel ./ sum(gKernel);

t_win = [-0.5 2];
lineW = 4;
fontSize = 50;

%% Condition definitions
condFields = {'msTS_c25_bl', 'msTS_c50_bl', 'msTS_c75_bl', 'msTS_c100_bl'};
condLabels = {' 25% Contrast', ' 50% Contrast', ' 75% Contrast', ' 100% Contrast'};
nConds     = length(condFields);

%% Process all conditions
fprintf('\n[VIZ GAZE MS] Processing contrast conditions\n');

subjCurves = cell(nSubj, nConds);
t_vec = [];

for subj = 1:nSubj
    clc; fprintf('[VIZ GAZE MS] Subject %d/%d (%s)\n', subj, nSubj, subjects{subj});
    spath = fullfile(datapath, subjects{subj}, 'gaze');
    dat = load(fullfile(spath, 'gaze_microsaccade_timeseries.mat'));

    for c = 1:nConds
        msTS = dat.(condFields{c});
        ch = find(strcmp(msTS.label, 'MSRate'), 1, 'first');
        thisTime = msTS.time(:)';
        thisRate = msTS.avg(ch, :);
        idxDisp = thisTime >= t_win(1) & thisTime <= t_win(2);
        thisTimeDisp = thisTime(idxDisp);
        thisRateDisp = conv(thisRate(idxDisp), gKernel, 'same');

        if isempty(t_vec)
            t_vec = thisTimeDisp;
        end
        subjCurves{subj, c} = thisRateDisp;
    end
end

%% Assemble subject x time x condition matrix
n_disp = numel(t_vec);
subjRates = nan(nSubj, n_disp, nConds);
for subj = 1:nSubj
    for c = 1:nConds
        subjRates(subj, :, c) = subjCurves{subj, c};
    end
end

%% Grand averages
grandMean = squeeze(nanmean(subjRates, 1));                                % n_disp x nConds
nValid_ts = squeeze(sum(~isnan(subjRates), 1));                            % n_disp x nConds
grandSEM  = squeeze(nanstd(subjRates, 0, 1)) ./ sqrt(max(nValid_ts, 1));  % n_disp x nConds
grandSEM(nValid_ts < 2) = NaN;

%% FIGURE % change MS rate time courses per condition
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on

% Per-condition lines with SEM shading
for c = 1:nConds
    mu  = grandMean(:, c);
    sem = grandSEM(:, c);

    eb = shadedErrorBar(t_vec, mu, sem, 'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', lineW);
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.2);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end

xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
xlim(t_win);
%ylim([-1.75 1.25])
xlabel('Time [s]', 'FontSize', fontSize*0.8);
ylabel('Microsaccade Rate [%]', 'FontSize', fontSize*0.8);
leg_p = gobjects(nConds, 1);
for c = 1:nConds
    leg_p(c) = patch(nan, nan, colors(c, :), 'FaceAlpha', 0.33, ...
        'EdgeColor', colors(c, :), 'LineWidth', 1.5);
end
set(gca, 'FontSize', fontSize*0.8);
legend(leg_p, condLabels, 'Location', 'northeast', 'FontSize', fontSize*0.65, 'Box', 'off');
box off
hold off
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(figpath, 'GCP_gaze_microsaccades_rate.png'), '-dpng', '-r600');