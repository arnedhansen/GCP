%% GCP Stats Boxplots (Full / Early / Late)
%
% Subject-level boxplots across contrast for gamma and gaze metrics.
% Windows: full (0-2 s), early (0-1 s), late (1-2 s).
%
% Metrics (all windows):
%   Frequency, Power
%   MSRate_bl, BCEA_bl, Vel2D_bl, PupilSize_bl
%   Blinks_bl, Fixations_bl, Saccades_bl
%
% Data (precomputed; no window recomputation here):
%   GCP_eeg_GED.mat
%   GCP_gaze_window_summaries.mat  (from GCP_gaze_fex.m)
%
% Output: figures/stats/boxplots/
%   Full window: GCP_stats_boxplot_<Metric>.png (assembly-compatible)
%   Early/late:  GCP_stats_boxplot_<Metric>_<early|late>.png

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP', 0);
subjects = gcp_subject_inclusion(subjects, paths);
nSubj = numel(subjects);
fprintf('Included GED cohort: N = %d (%s)\n', nSubj, strjoin(subjects, ', '));

nCond = 4;
xtickLabs = {'25%', '50%', '75%', '100%'};
winNames = {'full', 'early', 'late'};

out_dir = fullfile(paths.figures, 'stats', 'boxplots');
if ~isfolder(out_dir)
    mkdir(out_dir);
end

%% Aesthetics
fontSize       = 50;
axisLabelSize  = fontSize * 0.9;
yTickFontSize  = axisLabelSize;
xLabelFontSize = axisLabelSize;
xTickFontSize  = axisLabelSize * (1 - 0.33);
dotSize        = 400;
dotAlpha       = 0.85;
jitter         = 0.4;
boxWidth       = 0.55;

%% Load gamma data
ged = load(fullfile(paths.features, 'GCP_eeg_GED.mat'), ...
    'trials_median', 'trials_median_early', 'trials_median_late', ...
    'trials_gamma_power', 'trials_gamma_power_early', 'trials_gamma_power_late', ...
    'subjects');
ged_idx = match_subjects(ged.subjects, subjects);

gamma = struct();
gamma.Frequency.full  = pick_matrix(ged.trials_median, ged_idx);
gamma.Frequency.early = pick_matrix(ged.trials_median_early, ged_idx);
gamma.Frequency.late  = pick_matrix(ged.trials_median_late, ged_idx);
gamma.Power.full      = pick_matrix(ged.trials_gamma_power, ged_idx);
gamma.Power.early     = pick_matrix(ged.trials_gamma_power_early, ged_idx);
gamma.Power.late      = pick_matrix(ged.trials_gamma_power_late, ged_idx);

%% Gaze
gazePath = fullfile(paths.features, 'GCP_gaze_window_summaries.mat');
if ~isfile(gazePath)
    error('GCP_stats_boxplots:MissingGazeSummaries', ...
        ['Missing %s. Run GCP_gaze_fex.m first so full/early/late ', ...
         'subject x condition scalars exist.'], gazePath);
end
gazeSum = load(gazePath);
gaze_idx = match_subjects(gazeSum.subjects, subjects);

gazeMetrics = {'MSRate_bl','BCEA_bl','Vel2D_bl','PupilSize_bl', ...
    'Blinks_bl','Fixations_bl','Saccades_bl'};
gaze = struct();
for mi = 1:numel(gazeMetrics)
    name = gazeMetrics{mi};
    if ~isfield(gazeSum, name)
        error('GCP_stats_boxplots:MissingMetric', ...
            'Metric %s missing in %s', name, gazePath);
    end
    for wi = 1:numel(winNames)
        wn = winNames{wi};
        if ~isfield(gazeSum.(name), wn)
            error('GCP_stats_boxplots:MissingWindow', ...
                'Metric %s window %s missing in %s', name, wn, gazePath);
        end
        gaze.(name).(wn) = pick_matrix(gazeSum.(name).(wn), gaze_idx);
    end
end

%% Plot specs: {varName, windowsStruct, yLabel, drawZero}
plotSpecs = {
    'Frequency',   gamma.Frequency,    'Frequency [Hz]',        false
    'Power',       gamma.Power,        'Power [dB]',            true
    'MSRate_bl',    gaze.MSRate_bl,      'Microsaccade Rate [%]', true
    'BCEA_bl',      gaze.BCEA_bl,        'BCEA [%]',              true
    'Vel2D_bl',     gaze.Vel2D_bl,       'Eye Velocity [%]',      true
    'PupilSize_bl', gaze.PupilSize_bl,   'Pupil Size [%]',        true
    'Blinks_bl',    gaze.Blinks_bl,      'Blinks [%]',            true
    'Fixations_bl', gaze.Fixations_bl,   'Fixations [%]',         true
    'Saccades_bl',  gaze.Saccades_bl,    'Saccades [%]',          true
    };

subjIDs = str2double(string(subjects(:)));

%% Loop
for iMetric = 1:size(plotSpecs, 1)
    varName = plotSpecs{iMetric, 1};
    winMats = plotSpecs{iMetric, 2};
    ylabStr = plotSpecs{iMetric, 3};
    drawZero = plotSpecs{iMetric, 4};

    for iWin = 1:numel(winNames)
        winName = winNames{iWin};
        if ~isfield(winMats, winName) || isempty(winMats.(winName))
            continue
        end
        M = winMats.(winName);
        if ~any(isfinite(M(:)))
            continue
        end

        y = nan(nCond * nSubj, 1);
        condIdx = nan(nCond * nSubj, 1);
        idVec = nan(nCond * nSubj, 1);
        row = 0;
        for s = 1:nSubj
            for c = 1:nCond
                row = row + 1;
                idVec(row) = subjIDs(s);
                condIdx(row) = c;
                if c <= size(M, 1) && s <= size(M, 2)
                    y(row) = M(c, s);
                end
            end
        end

        fprintf('%s %s: %d finite of %d cells\n', ...
            varName, winName, nnz(isfinite(y)), numel(y));

        close all
        figure('Position', [0 0 1512 982], 'Color', 'w');
        hold on

        xJit = nan(size(y));
        for s = 1:nSubj
            idxSubj = idVec == subjIDs(s);
            for c = 1:nCond
                idxPt = idxSubj & condIdx == c & isfinite(y);
                if any(idxPt)
                    xJit(idxPt) = c + jitter * (rand - 0.5);
                end
            end
        end

        for c = 1:nCond
            idxC = condIdx == c & isfinite(y);
            yC = y(idxC);
            if isempty(yC)
                continue
            end
            boxplot(yC, ones(numel(yC), 1), 'Positions', c, 'Symbol', '', ...
                'Widths', boxWidth, 'Colors', 'k');
        end
        styleCurrentBoxplot(colors(1:nCond, :));

        yl = ylim;
        if drawZero && yl(1) <= 0 && yl(2) >= 0
            yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        end

        for s = 1:nSubj
            idxSubj = idVec == subjIDs(s);
            xSubj = xJit(idxSubj);
            ySubj = y(idxSubj);
            condSubjIdx = condIdx(idxSubj);
            valid = isfinite(ySubj) & isfinite(xSubj);
            xSubj = xSubj(valid);
            ySubj = ySubj(valid);
            condSubjIdx = condSubjIdx(valid);
            if numel(ySubj) < 2
                continue
            end
            [~, sortIdx] = sort(condSubjIdx);
            plot(xSubj(sortIdx), ySubj(sortIdx), '-', ...
                'Color', [0.8 0.8 0.8], 'LineWidth', 1);
        end

        for c = 1:nCond
            idxC = condIdx == c & isfinite(y);
            scatter(xJit(idxC), y(idxC), dotSize, colors(c, :), 'filled', ...
                'MarkerFaceAlpha', dotAlpha);
        end

        xlim([0.5 nCond + 0.5]);
        yValid = y(isfinite(y));
        if ~isempty(yValid)
            yPad = 0.05 * range(yValid);
            ylim([min(yValid) - yPad, max(yValid) + yPad]);
        end
        xticks(1:nCond);
        xticklabels(xtickLabs);

        hXlab = xlabel('Contrast', 'FontSize', xLabelFontSize);
        hYlab = ylabel(ylabStr, 'Interpreter', 'none', 'FontSize', axisLabelSize);

        ax = gca;
        ax.XAxis.FontSize = xTickFontSize;
        ax.YAxis.FontSize = yTickFontSize;
        hXlab.FontSize = xLabelFontSize;
        hYlab.FontSize = axisLabelSize;
        set(ax, 'Box', 'off');
        box off;
        hold off;

        if strcmp(winName, 'full')
            outName = sprintf('GCP_stats_boxplot_%s.png', varName);
        else
            outName = sprintf('GCP_stats_boxplot_%s_%s.png', varName, winName);
        end
        outPath = fullfile(out_dir, outName);
        drawnow;
        set(gcf, 'PaperPositionMode', 'auto');
        print(gcf, outPath, '-dpng', '-r600');
        fprintf('Saved %s\n', outPath);
    end
end

fprintf('\nDone. Boxplots in:\n  %s\n', out_dir);

%% Local functions
function idx = match_subjects(allSubjects, keepSubjects)
idx = arrayfun(@(s) find(strcmp(allSubjects, keepSubjects{s}), 1), ...
    1:numel(keepSubjects));
if any(cellfun(@isempty, num2cell(idx)))
    error('GCP_stats_boxplots:SubjectMismatch', ...
        'One or more included subjects are missing from a data source.');
end
end

function M = pick_matrix(src, subj_idx)
M = nan(size(src, 1), numel(subj_idx));
for s = 1:numel(subj_idx)
    si = subj_idx(s);
    if si <= size(src, 2)
        M(:, s) = src(:, si);
    end
end
end

function styleCurrentBoxplot(boxColors)
hBoxes = findobj(gca, 'Tag', 'Box');
if isempty(hBoxes)
    return
end
nB = numel(hBoxes);
mx = zeros(nB, 1);
for ii = 1:nB
    xd = get(hBoxes(ii), 'XData');
    mx(ii) = mean(xd(:), 'omitnan');
end
[~, ord] = sort(mx);
hBoxes = hBoxes(ord);
for bi = 1:min(numel(hBoxes), size(boxColors, 1))
    patch(get(hBoxes(bi), 'XData'), get(hBoxes(bi), 'YData'), ...
        boxColors(bi, :), 'FaceAlpha', 0.25, ...
        'EdgeColor', boxColors(bi, :), 'LineWidth', 1.5);
end
end
