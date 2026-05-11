%% GCP Gaze Microsaccade Suppression — Contrast Conditions
% Loads already processed microsaccade percentage-change time courses from
% 2_feature_extraction/GCP_gaze_fex.m (msTS_cXX_bl_pct) and plots
% per-condition grand-average traces with SEM shading.
%
% Important:
%   - Baseline correction has already been applied during feature extraction
%     (FieldTrip relativechange, then multiplied by 100).
%   - A light additional display smoothing is applied here.

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');

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

t_win = [-1 2];
fontSize = 20;

%% Condition definitions
condFields = {'msTS_c25_bl_pct', 'msTS_c50_bl_pct', 'msTS_c75_bl_pct', 'msTS_c100_bl_pct'};
condLabels = {'25% contrast', '50% contrast', '75% contrast', '100% contrast'};
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
subjRates_pct = nan(nSubj, n_disp, nConds);

for subj = 1:nSubj
    for c = 1:nConds
        if ~isempty(subjCurves{subj, c})
            subjRates_pct(subj, :, c) = subjCurves{subj, c};
        end
    end
end

%% Grand averages
grandMean_pct = squeeze(nanmean(subjRates_pct, 1));                                % n_disp x nConds
nValid_ts     = squeeze(sum(~isnan(subjRates_pct), 1));                            % n_disp x nConds
grandSEM_pct  = squeeze(nanstd(subjRates_pct, 0, 1)) ./ sqrt(max(nValid_ts, 1));  % n_disp x nConds
grandSEM_pct(nValid_ts < 2) = NaN;

%% FIGURE Percentage-Change MS Rate Time Courses per Condition
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on

% Per-condition lines with SEM shading
for c = 1:nConds
    mu  = grandMean_pct(:, c);
    sem = grandSEM_pct(:, c);

    % SEM ribbon
    fill([t_vec, fliplr(t_vec)], ...
         [(mu + sem)', fliplr((mu - sem)')], ...
         colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');

    plot(t_vec, mu, '-', 'Color', colors(c, :), 'LineWidth', 2.5);
end

xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlim(t_win);
xlabel('Time [s]');
ylabel('Microsaccade Rate [%]');
leg_p_pct = gobjects(nConds, 1);
for c = 1:nConds
    leg_p_pct(c) = patch(NaN, NaN, colors(c, :), 'EdgeColor', 'none');
end
legend(leg_p_pct, condLabels, 'Location', 'best', 'FontSize', fontSize - 4, 'Box', 'off');
set(gca, 'FontSize', fontSize);
hold off

exportgraphics(gcf, fullfile(figpath, 'GCP_gaze_microsaccades_rate_pct.png'), 'Resolution', 600);