%% GCP Gaze Fixation Time Course — Contrast Conditions
% Builds fixation-rate time courses per contrast condition from
% EyeLink fixation events. Fixation rate is reconstructed from fixation
% onset events (L_fixation, R_fixation) as fixations/s, then baseline-
% normalised as percentage change and plotted with SEM.
%
% This script mirrors the microsaccade plotting workflow:
%   1. Fixation-rate time courses (raw fixations/s)
%   1b. Percentage-change baseline-corrected fixation rate
%   2. Raster plots of fixation onsets
%   3. Combined figure (rate + raster)
%   3b. Combined figure with percentage-change rate

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');

figpath  = fullfile(paths.figures, 'gaze', 'fixations');
mkdir(figpath);

nSubj = length(subjects);

%% Parameters
fsample = 500;          % Hz (updated from data when available)

% Computation window (wider than display to avoid edge artefacts)
t_comp     = [-2.0 2.5];                                % seconds
n_comp     = round(diff(t_comp) * fsample) + 1;         % samples in computation window
t_comp_vec = linspace(t_comp(1), t_comp(2), n_comp);    % computation time axis

% Storage window (includes baseline + display)
t_store = [-1.5 2];                                      % seconds
[~, store_start] = min(abs(t_comp_vec - t_store(1)));
[~, store_end]   = min(abs(t_comp_vec - t_store(2)));
store_idx   = store_start : store_end;
n_store     = length(store_idx);
t_store_vec = t_comp_vec(store_idx);

% Display time window (for plotting)
t_win = [-1 2];                                          % seconds
disp_idx = t_store_vec >= t_win(1) & t_store_vec <= t_win(2);
n_disp   = sum(disp_idx);
t_vec    = t_store_vec(disp_idx);

% Raster crop (display window in computation axis)
[~, raster_crop_start] = min(abs(t_comp_vec - t_win(1)));
[~, raster_crop_end]   = min(abs(t_comp_vec - t_win(2)));
raster_crop_idx = raster_crop_start : raster_crop_end;

% Baseline period for percentage-change normalisation
bl_win = [-1 -0.25];                                     % seconds
bl_idx = t_store_vec >= bl_win(1) & t_store_vec <= bl_win(2);

% Plotting
fontSize = 20;

%% Condition definitions
condCodes  = {'61', '62', '63', '64'};
condLabels = {'25% contrast', '50% contrast', '75% contrast', '100% contrast'};
nConds     = length(condCodes);

%% ====================================================================
%  Process all conditions
%  ====================================================================
fprintf('\n=== Processing GCP Fixation Rate Conditions ===\n');

% Preallocate: subjects x timepoints (storage window) x conditions
subjFix = nan(nSubj, n_store, nConds);
rasterAll = cell(nConds, 1);

for subj = 1:nSubj
    fprintf('  Subject %d/%d (%s)\n', subj, nSubj, subjects{subj});
    subjMergedPath = fullfile(paths.merged, subjects{subj});

    for c = 1:nConds
        condCode = condCodes{c};

        condOnset = [];

        for block = 1:4
            mergedFile = fullfile(subjMergedPath, ...
                sprintf('%s_EEG_ET_GCP_block%d_merged.mat', subjects{subj}, block));
            if ~isfile(mergedFile)
                continue
            end

            B = load(mergedFile, 'EEG');
            if ~isfield(B, 'EEG')
                error('Merged file does not contain EEG: %s', mergedFile);
            end

            EEG = B.EEG;
            if ~isfield(EEG, 'event') || isempty(EEG.event)
                error('No event array in merged EEG for %s block %d.', subjects{subj}, block);
            end

            EEG_ep = pop_epoch(EEG, {condCode}, t_comp);
            if EEG_ep.trials < 1
                continue
            end

            fsample = EEG_ep.srate;
            nTrialSamples = EEG_ep.pnts;
            nTrialsBlock = EEG_ep.trials;

            if ~all(isfield(EEG_ep.event, {'type', 'latency', 'epoch', 'duration'}))
                error(['Missing required event fields in %s block %d, cond %s. ' ...
                       'Expected: type, latency, epoch, duration.'], ...
                      subjects{subj}, block, condCode);
            end

            eventTypes = cell(1, numel(EEG_ep.event));
            for ev = 1:numel(EEG_ep.event)
                thisType = EEG_ep.event(ev).type;
                if ischar(thisType)
                    eventTypes{ev} = thisType;
                elseif isstring(thisType)
                    eventTypes{ev} = char(thisType);
                elseif iscell(thisType) && ~isempty(thisType)
                    if ischar(thisType{1})
                        eventTypes{ev} = thisType{1};
                    elseif isstring(thisType{1})
                        eventTypes{ev} = char(thisType{1});
                    else
                        error('Unsupported event.type cell content in %s block %d.', subjects{subj}, block);
                    end
                else
                    error('Unsupported event.type format in %s block %d.', subjects{subj}, block);
                end
            end

            isFixEvent = strcmp(eventTypes, 'L_fixation') | strcmp(eventTypes, 'R_fixation');

            for trl = 1:nTrialsBlock
                onsetVec = zeros(1, nTrialSamples);

                trlMask = [EEG_ep.event.epoch] == trl;
                trlFixEvents = find(isFixEvent & trlMask);

                for ev = 1:numel(trlFixEvents)
                    e = EEG_ep.event(trlFixEvents(ev));

                    latencyVal = double(e.latency);
                    durationVal = double(e.duration);
                    if ~isfinite(latencyVal) || ~isfinite(durationVal) || durationVal <= 0
                        error('Invalid fixation latency/duration in %s block %d, cond %s, trial %d.', ...
                            subjects{subj}, block, condCode, trl);
                    end

                    onsetTrial = mod(round(latencyVal) - 1, nTrialSamples) + 1;
                    durSamples = round(durationVal);
                    offsetTrial = min(nTrialSamples, onsetTrial + durSamples - 1); %#ok<NASGU>
                    onsetVec(onsetTrial) = 1;
                end

                if nTrialSamples >= n_comp
                    onsetVec = onsetVec(1:n_comp);
                else
                    onsetVec(end+1 : n_comp) = 0;
                end

                condOnset(end+1, :) = onsetVec;
            end
        end

        if size(condOnset, 1) >= 3
            % Mean onset probability per sample -> fixation rate (Hz),
            % then crop to storage window [-1.5, 2] for baseline/display.
            rateComp = mean(condOnset, 1) * fsample;
            subjFix(subj, :, c) = rateComp(store_idx);
            if size(condOnset, 2) < max(raster_crop_idx)
                error('Raster crop exceeds available samples in %s, cond %s.', subjects{subj}, condCode);
            end
            rasterAll{c} = [rasterAll{c}; condOnset(:, raster_crop_idx)];
        end
    end
end

%% Percentage-change baseline normalisation (per subject, per condition)
subjFix_pct = nan(size(subjFix));
for subj = 1:nSubj
    for c = 1:nConds
        fix_ts = subjFix(subj, :, c);
        if all(isnan(fix_ts)); continue; end
        bl_mean = nanmean(fix_ts(bl_idx));
        if bl_mean == 0 || isnan(bl_mean); continue; end
        subjFix_pct(subj, :, c) = (fix_ts - bl_mean) ./ bl_mean * 100;
    end
end

%% Grand averages — display portion only (for plotting)
% Raw fixation rate (fixations/s)
subjFix_disp = subjFix(:, disp_idx, :);
grandMean = squeeze(nanmean(subjFix_disp, 1));                        % n_disp x nConds
nValid    = squeeze(sum(~isnan(subjFix_disp(:, 1, :)), 1));           % nConds x 1
grandSEM  = squeeze(nanstd(subjFix_disp, 0, 1)) ./ sqrt(nValid');     % n_disp x nConds

% Percentage change
subjFix_pct_disp = subjFix_pct(:, disp_idx, :);
grandMean_pct = squeeze(nanmean(subjFix_pct_disp, 1));                          % n_disp x nConds
nValid_pct    = squeeze(sum(~isnan(subjFix_pct_disp(:, 1, :)), 1));             % nConds x 1
grandSEM_pct  = squeeze(nanstd(subjFix_pct_disp, 0, 1)) ./ sqrt(nValid_pct');   % n_disp x nConds

%% ================================================================
%  FIGURE 1: Fixation Rate Time Courses per Condition (raw fixations/s)
%  ================================================================
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on

for c = 1:nConds
    mu  = grandMean(:, c);
    sem = grandSEM(:, c);

    fill([t_vec, fliplr(t_vec)], ...
         [(mu + sem)', fliplr((mu - sem)')], ...
         colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');

    plot(t_vec, mu, '-', 'Color', colors(c, :), 'LineWidth', 2.5);
end

xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlim(t_win);
xlabel('Time [s]');
ylabel('Fixations/s');
title('Fixation Rate Time Course');
leg_p = gobjects(nConds, 1);
for c = 1:nConds
    leg_p(c) = patch(NaN, NaN, colors(c, :), 'EdgeColor', 'none');
end
legend(leg_p, condLabels, 'Location', 'northeast', 'FontSize', fontSize - 4, 'Box', 'off');
set(gca, 'FontSize', fontSize);
hold off

exportgraphics(gcf, fullfile(figpath, 'GCP_gaze_fixations_rate.png'), 'Resolution', 600);

%% ================================================================
%  FIGURE 1b: Baseline-Relative Fixation Rate (% change)
%  ================================================================
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on

for c = 1:nConds
    mu  = grandMean_pct(:, c);
    sem = grandSEM_pct(:, c);

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
ylabel('Fixations [%]');
title('Fixation Rate Time Course');
leg_p_pct = gobjects(nConds, 1);
for c = 1:nConds
    leg_p_pct(c) = patch(NaN, NaN, colors(c, :), 'EdgeColor', 'none');
end
legend(leg_p_pct, condLabels, 'Location', 'northeast', 'FontSize', fontSize - 4, 'Box', 'off');
set(gca, 'FontSize', fontSize);
hold off

exportgraphics(gcf, fullfile(figpath, 'GCP_gaze_fixations_rate_pct.png'), 'Resolution', 600);

%% ================================================================
%  FIGURE 2: Raster Plots of Fixation Onsets
%  ================================================================
figure('Position', [0 0 1512 982], 'Color', 'w');

for c = 1:nConds
    subplot(1, nConds, c);
    hold on

    raster = rasterAll{c};
    if isempty(raster)
        title(condLabels{c}, 'FontSize', fontSize);
        continue
    end

    nR = size(raster, 1);

    maxShow = 300;
    if nR > maxShow
        rIdx = sort(randperm(nR, maxShow));
        raster = raster(rIdx, :);
        nR = maxShow;
    end

    nSampPlot = min(size(raster, 2), n_disp);
    [trialIdx, sampleIdx] = find(raster(:, 1:nSampPlot));
    if ~isempty(trialIdx)
        plot(t_vec(sampleIdx), trialIdx, '.', ...
            'Color', colors(c, :), 'MarkerSize', 12);
    end

    xline(0, 'k--', 'LineWidth', 1.5);
    xlim(t_win);
    ylim([0 nR + 1]);
    xlabel('Time [s]');
    if c == 1; ylabel('Trial'); end
    title(condLabels{c}, 'FontSize', fontSize);
    set(gca, 'FontSize', fontSize - 4, 'YDir', 'reverse');
    hold off
end

sgtitle('Fixation Onset Raster', 'FontSize', fontSize + 2);

exportgraphics(gcf, fullfile(figpath, 'GCP_gaze_fixations_onset_raster.png'), 'Resolution', 600);

%% ================================================================
%  FIGURE 3: Combined (Fixation Rate + Raster)
%  ================================================================
figure('Position', [0 0 1512 982], 'Color', 'w');

% --- Top panel: Fixation rate ---
ax1 = subplot(2, 1, 1);
hold on
for c = 1:nConds
    mu  = grandMean(:, c);
    sem = grandSEM(:, c);
    fill([t_vec, fliplr(t_vec)], ...
         [(mu + sem)', fliplr((mu - sem)')], ...
         colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
    plot(t_vec, mu, '-', 'Color', colors(c, :), 'LineWidth', 2.5);
end
xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlim(t_win);
ylabel('Fixations/s');
title('Fixation Dynamics');
leg_p_comb = gobjects(nConds, 1);
for c = 1:nConds
    leg_p_comb(c) = patch(NaN, NaN, colors(c, :), 'EdgeColor', 'none');
end
legend(leg_p_comb, condLabels, 'Location', 'northeast', 'FontSize', fontSize - 4, 'Box', 'off');
set(gca, 'FontSize', fontSize, 'XTickLabel', []);
hold off

% --- Bottom panel: Raster grouped by condition ---
ax2 = subplot(2, 1, 2);
hold on
yOff      = 0;
yTicks    = [];
yTickLbls = {};

for c = 1:nConds
    raster = rasterAll{c};
    if isempty(raster); continue; end

    nR = min(size(raster, 1), 100);
    raster = raster(1:nR, :);

    nSampPlot = min(size(raster, 2), n_disp);
    [trialIdx, sampleIdx] = find(raster(:, 1:nSampPlot));
    if ~isempty(trialIdx)
        plot(t_vec(sampleIdx), trialIdx + yOff, '.', ...
            'Color', colors(c, :), 'MarkerSize', 12);
    end

    yTicks(end+1)    = yOff + nR / 2;
    yTickLbls{end+1} = condLabels{c};
    yOff = yOff + nR + 5;
end

xline(0, 'k--', 'LineWidth', 1.5);
xlim(t_win);
ylim([0 yOff]);
xlabel('Time [s]');
ylabel('Trials');
set(gca, 'FontSize', fontSize, 'YDir', 'reverse', ...
    'YTick', yTicks, 'YTickLabel', yTickLbls);
hold off

linkaxes([ax1, ax2], 'x');

exportgraphics(gcf, fullfile(figpath, 'GCP_gaze_fixations_combined.png'), 'Resolution', 600);

%% ================================================================
%  FIGURE 3b: Combined Baseline-Relative Fixation Rate (% change) + Raster
%  ================================================================
figure('Position', [0 0 1512 982], 'Color', 'w');

% --- Top panel: Percentage-change fixation rate ---
ax1p = subplot(2, 1, 1);
hold on
for c = 1:nConds
    mu  = grandMean_pct(:, c);
    sem = grandSEM_pct(:, c);
    fill([t_vec, fliplr(t_vec)], ...
         [(mu + sem)', fliplr((mu - sem)')], ...
         colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
    plot(t_vec, mu, '-', 'Color', colors(c, :), 'LineWidth', 2.5);
end
xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlim(t_win);
ylabel('Fixations [%]');
title('Fixation Dynamics');
leg_p_comb_pct = gobjects(nConds, 1);
for c = 1:nConds
    leg_p_comb_pct(c) = patch(NaN, NaN, colors(c, :), 'EdgeColor', 'none');
end
legend(leg_p_comb_pct, condLabels, 'Location', 'southeast', 'FontSize', fontSize - 4, 'Box', 'off');
set(gca, 'FontSize', fontSize, 'XTickLabel', []);
hold off

% --- Bottom panel: Raster grouped by condition ---
ax2p = subplot(2, 1, 2);
hold on
yOff      = 0;
yTicks    = [];
yTickLbls = {};

for c = 1:nConds
    raster = rasterAll{c};
    if isempty(raster); continue; end

    nR = min(size(raster, 1), 100);
    raster = raster(1:nR, :);

    nSampPlot = min(size(raster, 2), n_disp);
    [trialIdx, sampleIdx] = find(raster(:, 1:nSampPlot));
    if ~isempty(trialIdx)
        plot(t_vec(sampleIdx), trialIdx + yOff, '.', ...
            'Color', colors(c, :), 'MarkerSize', 12);
    end

    yTicks(end+1)    = yOff + nR / 2;
    yTickLbls{end+1} = condLabels{c};
    yOff = yOff + nR + 5;
end

xline(0, 'k--', 'LineWidth', 1.5);
xlim(t_win);
ylim([0 yOff]);
xlabel('Time [s]');
ylabel('Trials');
set(gca, 'FontSize', fontSize, 'YDir', 'reverse', ...
    'YTick', yTicks, 'YTickLabel', yTickLbls);
hold off

linkaxes([ax1p, ax2p], 'x');

exportgraphics(gcf, fullfile(figpath, 'GCP_gaze_fixations_combined_pct.png'), 'Resolution', 600);

fprintf('\n=== All fixation figures saved to %s ===\n', figpath);
