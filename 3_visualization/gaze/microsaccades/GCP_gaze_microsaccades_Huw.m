%% GCP Gaze Microsaccade Suppression — Contrast Conditions
% Adapted from AOC_gaze_microsaccades_Huw.m (Huw Swanborough, InfCry project).
%
% Detects microsaccades per trial using Engbert & Kliegl (2003), builds
% binary spike trains, convolves with a Gaussian kernel to produce smoothed
% microsaccade rate time courses, and plots per-condition rate curves with
% SEM shading and raster plots.
%
% Edge effects from the Gaussian smoothing kernel are handled by computing
% over a wider window (t_comp) and cropping to the display window (t_win).
%
% Key outputs:
%   1. Smoothed MS rate time courses per contrast condition (± SEM)
%   1b. Percentage-change baseline-corrected rate time courses
%   2. Raster plots of microsaccade onsets per condition
%   3. Combined figure (rate + raster)
%   3b. Combined figure with percentage-change rate

%% Setup
startup
[subjects, path, colors, ~] = setup('GCP');
close all
clc

datapath = path;
figpath  = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/gaze/microsaccades/';
mkdir(figpath);

nSubj = length(subjects);

%% Parameters
fsample   = 500;        % Hz
screenW   = 800;
screenH   = 600;
blink_win = 50;         % blink removal window (samples, ±100 ms)

% Engbert & Kliegl (2003) microsaccade detection
velthres = 6;           % velocity threshold multiplier
mindur   = 6;           % minimum duration (samples) = 12 ms at 500 Hz

% Gaussian smoothing kernel for continuous rate estimation
sigma_ms   = 50;                                     % kernel SD in ms
sigma_samp = round(sigma_ms / (1000 / fsample));     % SD in samples
kHalf      = 3 * sigma_samp;
x_kern     = -kHalf : kHalf;
gKernel    = exp(-x_kern.^2 / (2 * sigma_samp^2));
gKernel    = gKernel / sum(gKernel);                 % normalise

% Computation window (wider than display to avoid convolution edge effects)
% Kernel half-width = 3*sigma = 150 ms; padding per side must exceed this.
% Must also cover the baseline period [-1.5, -0.25] for pct-change normalisation.
t_comp     = [-2.0 2.5];                                % seconds
n_comp     = round(diff(t_comp) * fsample) + 1;         % samples in computation window
t_comp_vec = linspace(t_comp(1), t_comp(2), n_comp);    % computation time axis

% Storage window (includes baseline + display, edge-effect-free)
t_store = [-1.5 2];                                      % seconds
[~, store_start] = min(abs(t_comp_vec - t_store(1)));
[~, store_end]   = min(abs(t_comp_vec - t_store(2)));
store_idx  = store_start : store_end;
n_store    = length(store_idx);
t_store_vec = t_comp_vec(store_idx);

% Display time window (what is shown in figures)
t_win = [-0.5 2];                                        % seconds
disp_idx = t_store_vec >= t_win(1) & t_store_vec <= t_win(2);  % logical into t_store_vec
n_disp   = sum(disp_idx);
t_vec    = t_store_vec(disp_idx);                        % display time axis

% Raster crop: indices in t_comp_vec corresponding to display window
[~, raster_crop_start] = min(abs(t_comp_vec - t_win(1)));
[~, raster_crop_end]   = min(abs(t_comp_vec - t_win(2)));
raster_crop_idx = raster_crop_start : raster_crop_end;

% Baseline period for percentage-change normalisation (matches GCP_gaze_fex.m)
bl_win = [-1.5 -0.25];                                  % seconds
bl_idx = t_store_vec >= bl_win(1) & t_store_vec <= bl_win(2);  % logical into t_store_vec

% Plotting
fontSize = 20;

%% Condition definitions
condFields = {'dataET_c25', 'dataET_c50', 'dataET_c75', 'dataET_c100'};
condLabels = {'25% contrast', '50% contrast', '75% contrast', '100% contrast'};
nConds     = length(condFields);

%% ====================================================================
%  Process all conditions
%  ====================================================================
fprintf('\n=== Processing GCP Contrast Conditions ===\n');

% Preallocate: subjects x timepoints (storage window) x conditions
subjRates = nan(nSubj, n_store, nConds);
rasterAll = cell(nConds, 1);

for subj = 1:nSubj
    fprintf('  Subject %d/%d (%s)\n', subj, nSubj, subjects{subj});
    spath = fullfile(datapath, subjects{subj}, 'gaze');

    % Load ET data (contains dataET_c25, dataET_c50, dataET_c75, dataET_c100)
    fpath = fullfile(spath, 'dataET.mat');
    if ~isfile(fpath)
        warning('File not found for subject %s', subjects{subj});
        continue
    end
    S = load(fpath);

    for c = 1:nConds
        cField = condFields{c};

        if ~isfield(S, cField)
            warning('No %s for subject %s', cField, subjects{subj});
            continue
        end
        et = S.(cField);

        nTrials    = length(et.trial);
        condSpikes = [];

        for trl = 1:nTrials
            raw = et.trial{trl};
            t   = et.time{trl};

            % Keep x, y, pupil; invert Y to screen coordinates
            raw = raw(1:min(3, size(raw, 1)), :);
            raw(2, :) = screenH - raw(2, :);

            % Out-of-bounds → NaN
            oob = raw(1,:) < 0 | raw(1,:) > screenW | ...
                  raw(2,:) < 0 | raw(2,:) > screenH;
            raw(:, oob) = NaN;

            % Remove blinks
            raw = remove_blinks(raw, blink_win);

            % Extract COMPUTATION window (wider than display for edge-free smoothing)
            idx_win = t >= t_comp(1) & t <= t_comp(2);
            gx = raw(1, idx_win);
            gy = raw(2, idx_win);

            if sum(isfinite(gx) & isfinite(gy)) < 50
                continue
            end

            % Detect microsaccade onsets → binary spike train
            spikeVec = detect_ms_spike_train(gx, gy, fsample, velthres, mindur);

            % Pad/trim to fixed computation-window length
            if length(spikeVec) >= n_comp
                spikeVec = spikeVec(1:n_comp);
            else
                spikeVec(end+1 : n_comp) = 0;
            end

            condSpikes(end+1, :) = spikeVec;
        end

        % Per-subject: average spike trains → rate → smooth → crop
        if size(condSpikes, 1) >= 3
            % Mean spike probability × fsample = rate (Hz)
            rate = mean(condSpikes, 1) * fsample;

            % Smooth with Gaussian kernel (over full computation window)
            smoothed = conv(rate, gKernel, 'same');

            % Crop to storage window (baseline + display, edge-effect-free)
            subjRates(subj, :, c) = smoothed(store_idx);

            % Accumulate rasters (cropped to display window only)
            rasterAll{c} = [rasterAll{c}; condSpikes(:, raster_crop_idx)];
        end
    end

    clear S et
end

%% Percentage-change baseline normalisation (per subject, per condition)
% Baseline is computed over [-1.5, -0.25] within the storage window,
% then the display portion [-0.5, 2] is extracted for plotting.
subjRates_pct = nan(size(subjRates));
for subj = 1:nSubj
    for c = 1:nConds
        rate_ts = subjRates(subj, :, c);
        if all(isnan(rate_ts)); continue; end
        bl_mean = nanmean(rate_ts(bl_idx));
        if bl_mean == 0 || isnan(bl_mean); continue; end
        subjRates_pct(subj, :, c) = (rate_ts - bl_mean) ./ bl_mean * 100;
    end
end

%% Grand averages — display portion only (for plotting)
% Raw Hz
subjRates_disp = subjRates(:, disp_idx, :);
grandMean = squeeze(nanmean(subjRates_disp, 1));                        % n_disp x nConds
nValid    = squeeze(sum(~isnan(subjRates_disp(:, 1, :)), 1));           % nConds x 1
grandSEM  = squeeze(nanstd(subjRates_disp, 0, 1)) ./ sqrt(nValid');    % n_disp x nConds

% Percentage change
subjRates_pct_disp = subjRates_pct(:, disp_idx, :);
grandMean_pct = squeeze(nanmean(subjRates_pct_disp, 1));                          % n_disp x nConds
nValid_pct    = squeeze(sum(~isnan(subjRates_pct_disp(:, 1, :)), 1));             % nConds x 1
grandSEM_pct  = squeeze(nanstd(subjRates_pct_disp, 0, 1)) ./ sqrt(nValid_pct');  % n_disp x nConds

%% ================================================================
%  FIGURE 1: Smoothed MS Rate Time Courses per Condition (raw Hz)
%  ================================================================
close all
figure('Color', 'w', 'Position', [0 0 1512 982]);
hold on

% Per-condition lines with SEM shading
h = gobjects(nConds, 1);
for c = 1:nConds
    mu  = grandMean(:, c);
    sem = grandSEM(:, c);

    % SEM ribbon (excluded from legend)
    fill([t_vec, fliplr(t_vec)], ...
         [(mu + sem)', fliplr((mu - sem)')], ...
         colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');

    h(c) = plot(t_vec, mu, '-', 'Color', colors(c, :), 'LineWidth', 2.5);
end

xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlim(t_win);
xlabel('Time [s]');
ylabel('Microsaccade Rate [Hz]');
title('GCP — Microsaccade Rate Time Course');
legend(h, condLabels, 'Location', 'northeast', 'FontSize', fontSize - 4);
set(gca, 'FontSize', fontSize);
hold off

saveas(gcf, fullfile(figpath, 'GCP_gaze_microsaccades_Huw_rate.png'));

%% ================================================================
%  FIGURE 1b: Percentage-Change MS Rate Time Courses per Condition
%  ================================================================
figure('Color', 'w', 'Position', [0 0 1512 982]);
hold on

% Per-condition lines with SEM shading
h_pct = gobjects(nConds, 1);
for c = 1:nConds
    mu  = grandMean_pct(:, c);
    sem = grandSEM_pct(:, c);

    % SEM ribbon
    fill([t_vec, fliplr(t_vec)], ...
         [(mu + sem)', fliplr((mu - sem)')], ...
         colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');

    h_pct(c) = plot(t_vec, mu, '-', 'Color', colors(c, :), 'LineWidth', 2.5);
end

xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlim(t_win);
xlabel('Time [s]');
ylabel('Microsaccade Rate [% change]');
title('GCP — Microsaccade Rate Time Course (Percentage Change)');
legend(h_pct, condLabels, 'Location', 'northeast', 'FontSize', fontSize - 4);
set(gca, 'FontSize', fontSize);
hold off

saveas(gcf, fullfile(figpath, 'GCP_gaze_microsaccades_Huw_rate_pct.png'));

%% ================================================================
%  FIGURE 2: Raster Plots per Condition
%  ================================================================
figure('Color', 'w', 'Position', [0 0 1512 982]);

for c = 1:nConds
    subplot(1, nConds, c);
    hold on

    raster = rasterAll{c};
    if isempty(raster)
        title(condLabels{c}, 'FontSize', fontSize);
        continue
    end

    nR = size(raster, 1);

    % Subsample if too many trials for readability
    maxShow = 300;
    if nR > maxShow
        rIdx = sort(randperm(nR, maxShow));
        raster = raster(rIdx, :);
        nR = maxShow;
    end

    % Plot each microsaccade onset as a dot
    nSampPlot = min(size(raster, 2), n_disp);
    [trialIdx, sampleIdx] = find(raster(:, 1:nSampPlot));
    if ~isempty(trialIdx)
        plot(t_vec(sampleIdx), trialIdx, '.', ...
            'Color', colors(c, :), 'MarkerSize', 10);
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

sgtitle('GCP — Microsaccade Raster', 'FontSize', fontSize + 2);

saveas(gcf, fullfile(figpath, 'GCP_gaze_microsaccades_Huw_raster.png'));

%% ================================================================
%  FIGURE 2b: Collapsed Raster (all conditions pooled)
%  ================================================================
figure('Color', 'w', 'Position', [0 0 1512 982]);
hold on

% Pool all conditions
rasterPooled = [];
for c = 1:nConds
    if ~isempty(rasterAll{c})
        rasterPooled = [rasterPooled; rasterAll{c}];
    end
end

if ~isempty(rasterPooled)
    nR = size(rasterPooled, 1);

    % Subsample if too many trials
    maxShow = 500;
    if nR > maxShow
        rIdx = sort(randperm(nR, maxShow));
        rasterPooled = rasterPooled(rIdx, :);
        nR = maxShow;
    end

    nSampPlot = min(size(rasterPooled, 2), n_disp);
    [trialIdx, sampleIdx] = find(rasterPooled(:, 1:nSampPlot));
    if ~isempty(trialIdx)
        plot(t_vec(sampleIdx), trialIdx, '.', ...
            'Color', [0.3 0.3 0.3], 'MarkerSize', 10);
    end
end

xline(0, 'k--', 'LineWidth', 1.5);
xlim(t_win);
ylim([0 nR + 1]);
xlabel('Time [s]');
ylabel('Trial');
title('GCP — Microsaccade Raster (All Conditions)');
set(gca, 'FontSize', fontSize, 'YDir', 'reverse');
hold off

saveas(gcf, fullfile(figpath, 'GCP_gaze_microsaccades_Huw_raster_collapsed.png'));

%% ================================================================
%  FIGURE 3: Combined (Rate on top + Raster on bottom)
%  ================================================================
figure('Color', 'w', 'Position', [0 0 1512 982]);

% --- Top panel: Smoothed rate ---
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
ylabel('MS Rate [Hz]');
title('GCP — Microsaccade Dynamics');
legend(condLabels, 'Location', 'northeast', 'FontSize', fontSize - 4);
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
            'Color', colors(c, :), 'MarkerSize', 10);
    end

    yTicks(end+1)    = yOff + nR / 2;
    yTickLbls{end+1} = condLabels{c};
    yOff = yOff + nR + 5;  % gap between conditions
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

saveas(gcf, fullfile(figpath, 'GCP_gaze_microsaccades_Huw_combined.png'));

%% ================================================================
%  FIGURE 3b: Combined Percentage-Change (Rate on top + Raster on bottom)
%  ================================================================
figure('Color', 'w', 'Position', [0 0 1512 982]);

% --- Top panel: Percentage-change rate ---
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
ylabel('MS Rate [% change]');
title('GCP — Microsaccade Dynamics (Percentage Change)');
legend(condLabels, 'Location', 'southeast', 'FontSize', fontSize - 4);
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
            'Color', colors(c, :), 'MarkerSize', 10);
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

saveas(gcf, fullfile(figpath, 'GCP_gaze_microsaccades_Huw_combined_pct.png'));

fprintf('\n=== All figures saved to %s ===\n', figpath);

%% ====================================================================
%  LOCAL FUNCTION: Microsaccade spike train detection
%  ====================================================================
function spikeVec = detect_ms_spike_train(gx, gy, fsample, velthres, mindur)
% DETECT_MS_SPIKE_TRAIN  Engbert & Kliegl (2003) microsaccade detection.
%   Returns a binary vector with 1s at microsaccade onset samples.
%
%   Adapted from detect_microsaccades.m and ICET_PreProcessingv2.m.
%   Uses FieldTrip's ft_preproc_padding for edge handling.

    nSamp    = length(gx);
    spikeVec = zeros(1, nSamp);
    velData  = [gx; gy];

    % Interpolate NaNs for velocity computation
    for ch = 1:2
        nanIdx = isnan(velData(ch, :));
        if all(nanIdx)
            return
        end
        if any(nanIdx)
            goodIdx = find(~nanIdx);
            velData(ch, nanIdx) = interp1(goodIdx, velData(ch, goodIdx), ...
                find(nanIdx), 'nearest', 'extrap');
        end
    end

    % Velocity via convolution kernel (Engbert et al. 2003, Eqn. 1)
    kernel = [1 1 0 -1 -1] .* (fsample / 6);
    n   = size(kernel, 2);
    pad = ceil(n / 2);
    dat = ft_preproc_padding(velData, 'localmean', pad);
    vel = convn(dat, kernel, 'same');
    vel = ft_preproc_padding(vel, 'remove', pad);

    % Velocity thresholds (Eqn. 2)
    medianstd = sqrt(median(vel.^2, 2) - (median(vel, 2)).^2);
    if any(medianstd == 0 | ~isfinite(medianstd))
        return
    end
    radius = velthres * medianstd;

    % Ellipse test (Eqn. 3)
    test   = sum((vel ./ radius(:, ones(1, nSamp))).^2, 1);
    sacsmp = find(test > 1);
    if isempty(sacsmp); return; end

    % Find consecutive above-threshold segments
    j = find(diff(sacsmp) == 1);
    if isempty(j); return; end

    j1  = [j; j + 1];
    com = intersect(j, j + 1);
    cut = ~ismember(j1, com);
    if sum(cut) < 2; return; end

    sacidx = reshape(j1(cut), 2, []);

    % Mark onset sample for segments meeting minimum duration
    for k = 1:size(sacidx, 2)
        dur = sacidx(1, k):sacidx(2, k);
        if length(dur) >= mindur
            onset = sacsmp(dur(1));
            if onset >= 1 && onset <= nSamp
                spikeVec(onset) = 1;
            end
        end
    end
end
