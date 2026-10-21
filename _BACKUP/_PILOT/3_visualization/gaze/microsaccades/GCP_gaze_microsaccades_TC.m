%% GCP Gaze Microsaccade Rate Contrast Conditions
% Loads the raw smoothed microsaccade rate in Hz from
% 2_feature_extraction/GCP_gaze_fex.m (msTS_cXX) and plots
% per-condition grand-average traces with SEM shading.
%
% Trials with no microsaccades remain in the average. A light additional
% display smoothing is applied here.

%% Setup
startup
set(0, 'DefaultFigureVisible', 'off');
[subjects, paths, colors, ~] = setup('GCP');
if ispc
    addpath('C:\Users\Administrator\Documents\GitHub\GCP\_BACKUP\_PILOT');
else
    addpath('/Users/Arne/Documents/GitHub/GCP/_BACKUP/_PILOT');
end
[paths, subjects] = gcp_pilot_paths(paths);
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
condFields = {'msTS_c25', 'msTS_c50', 'msTS_c75', 'msTS_c100'};
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
        % The saved series extends past the plotted window. Smooth that
        % longer series, then crop, so both ends of [-0.5, 2] use real samples.
        pad_s = kHalf / fsample;
        if thisTime(1) > t_win(1) - pad_s || thisTime(end) < t_win(2) + pad_s
            error('GCP_gaze_microsaccades_TC:ShortSeries', ...
                'Saved microsaccade series does not extend past the plotted window.');
        end
        thisRateSm = conv(thisRate, gKernel, 'same');
        idxDisp = thisTime >= t_win(1) & thisTime <= t_win(2);
        thisTimeDisp = thisTime(idxDisp);
        thisRateDisp = thisRateSm(idxDisp);

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
figure('Position', [0 0 1512 982], 'Color', 'w', 'Visible', 'off');
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
xlim(t_win);
%ylim([-1.75 1.25])
xlabel('Time [s]', 'FontSize', fontSize*0.8);
ylabel('Microsaccade Rate [Hz]', 'FontSize', fontSize*0.8);
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