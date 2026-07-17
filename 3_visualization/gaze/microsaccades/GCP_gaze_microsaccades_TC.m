%% GCP Gaze Microsaccade Suppression Contrast Conditions
% Loads already processed microsaccade percentage-change time courses from
% 2_feature_extraction/GCP_gaze_fex.m (msTS_cXX_bl_db) and plots
% per-condition grand-average traces with SEM shading.
%
% Important:
%   - Baseline correction has already been applied during feature extraction
%     as percentage change: 100*(stim/baseline - 1). Variable names retain *_bl_db.
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
condFields = {'msTS_c25_bl_db', 'msTS_c50_bl_db', 'msTS_c75_bl_db', 'msTS_c100_bl_db'};
condLabels = {' 25% Contrast', ' 50% Contrast', ' 75% Contrast', ' 100% Contrast'};
nConds     = length(condFields);

%% Process all conditions
fprintf('\n=== Processing GCP Contrast Conditions ===\n');

% Store per-subject, per-condition display traces
subjCurves = cell(nSubj, nConds);
t_vec = [];

for subj = 1:nSubj
    fprintf('  Subject %d/%d (%s)\n', subj, nSubj, subjects{subj});
    spath = fullfile(datapath, subjects{subj}, 'gaze');

    % Load already processed microsaccade time courses from feature extraction
    fpath = fullfile(spath, 'gaze_microsaccade_timeseries.mat');
    if ~isfile(fpath)
        warning('File not found for subject %s', subjects{subj});
        continue
    end
    dat = load(fpath);

    for c = 1:nConds
        cField = condFields{c};

        if ~isfield(dat, cField)
            warning('No %s for subject %s', cField, subjects{subj});
            continue
        end

        msTS = dat.(cField);
        if ~isfield(msTS, 'avg') || ~isfield(msTS, 'time') || isempty(msTS.avg) || isempty(msTS.time)
            warning('Invalid timelock structure in %s for subject %s', cField, subjects{subj});
            continue
        end

        % Select the microsaccade-rate channel explicitly by label when available.
        ch = 1;
        if isfield(msTS, 'label') && ~isempty(msTS.label)
            iCh = find(strcmp(msTS.label, 'MSRate'), 1, 'first');
            if ~isempty(iCh)
                ch = iCh;
            end
        end

        thisTime = msTS.time(:)';
        thisRate = msTS.avg(ch, :);
        idxDisp = thisTime >= t_win(1) & thisTime <= t_win(2);
        if ~any(idxDisp)
            continue
        end

        thisTimeDisp = thisTime(idxDisp);
        thisRateDisp = thisRate(idxDisp);
        thisRateDisp = conv(thisRateDisp, gKernel, 'same');

        if isempty(t_vec)
            t_vec = thisTimeDisp;
            subjCurves{subj, c} = thisRateDisp;
        else
            if numel(thisTimeDisp) == numel(t_vec) && max(abs(thisTimeDisp - t_vec)) < 1e-9
                subjCurves{subj, c} = thisRateDisp;
            else
                % Rare fallback if sampling grids differ across inputs.
                subjCurves{subj, c} = interp1(thisTimeDisp, thisRateDisp, t_vec, 'linear', NaN);
            end
        end
    end

    clear dat
end

%% Assemble subject x time x condition matrix
if isempty(t_vec)
    error('No valid microsaccade time-course data found.');
end
n_disp = numel(t_vec);
subjRates_db = nan(nSubj, n_disp, nConds);

for subj = 1:nSubj
    for c = 1:nConds
        if ~isempty(subjCurves{subj, c})
            subjRates_db(subj, :, c) = subjCurves{subj, c};
        end
    end
end

%% Grand averages
grandMean_db = squeeze(nanmean(subjRates_db, 1));                                % n_disp x nConds
nValid_ts    = squeeze(sum(~isnan(subjRates_db), 1));                            % n_disp x nConds
grandSEM_db  = squeeze(nanstd(subjRates_db, 0, 1)) ./ sqrt(max(nValid_ts, 1));  % n_disp x nConds
grandSEM_db(nValid_ts < 2) = NaN;

%% FIGURE dB-baselined MS rate time courses per condition
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on

% Per-condition lines with SEM shading
for c = 1:nConds
    mu  = grandMean_db(:, c);
    sem = grandSEM_db(:, c);

    eb = shadedErrorBar(t_vec, mu, sem, 'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', lineW);
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.125);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end

xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
xlim(t_win);
ylim([-1.75 1.25])
xlabel('Time [s]', 'FontSize', fontSize*0.8);
ylabel('Microsaccade Rate [%]', 'FontSize', fontSize*0.8);
leg_p_db = gobjects(nConds, 1);
for c = 1:nConds
    leg_p_db(c) = patch(nan, nan, colors(c, :), 'FaceAlpha', 0.25, ...
        'EdgeColor', colors(c, :), 'LineWidth', 1.5);
end
set(gca, 'FontSize', fontSize*0.8);
legend(leg_p_db, condLabels, 'Location', 'northeast', 'FontSize', fontSize*0.65, 'Box', 'off');
box off
hold off
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(figpath, 'GCP_gaze_microsaccades_rate_db.png'), '-dpng', '-r600');