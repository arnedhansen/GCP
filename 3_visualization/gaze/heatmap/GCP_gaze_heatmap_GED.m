%% GCP Gaze Pixel GED: pooled S_stim vs S_base, then condition contrast
%
% Per subject, trial-level gaze density maps are built at screen pixel
% resolution (800 x 600; one bin per pixel). Stimulus [0 2] and baseline
% [-1.5 -0.5] maps are duration-normalized (rate). All conditions are
% pooled to form second-moment matrices S_stim and S_base on
% support-masked pixels, reduced by uncentered SVD, then GED:
%   S_stim * w = lambda * S_base * w
%
% Second-moment (uncentered) matrices are used instead of demeaned
% covariance so mean gaze topography is retained. Demeaned covariance
% would emphasize trial-to-trial map variance and discard the mean
% spatial structure that answers where gaze sits.
%
% After the common stim-vs-base filter:
%   - Haufe activation pattern (S_stim * w) is mapped to pixel space
%   - Each trial map is projected to a scalar GED score
%   - Condition means of scores and of stim maps are contrasted
%   - Linear contrast maps (and |Haufe|-weighted versions) target where
%     gaze density changes with stimulus contrast
%
% Full-pixel GED without reduction is not tractable (~4.8e5 dimensions).
% Support masking plus uncentered SVD is required; with ~160 trials per
% condition the reduced-rank GED is well posed.
%
% Gaze preprocessing matches GCP_gaze_fex.m / GCP_gaze_heatmap_CBPT_outline.m:
% in-bounds mask, Y -> screen via (600 - y), remove_blinks (25 samples).
%
% Outputs (figures/gaze/heatmap/ged/ and data/features/):
%   subject Haufe patterns, linear contrast maps, condition scores,
%   group averages, and GCP_gaze_heatmap_GED.mat

%% Setup
startup
[subjects, paths, colors] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
pathFeat = paths.features;
figDir = fullfile(paths.figures, 'gaze', 'heatmap', 'ged');
if ~isfolder(figDir), mkdir(figDir); end

stimWindow = [0 2];
baselineWindow = [-1.5 -0.5];
stimDur = stimWindow(2) - stimWindow(1);
baseDur = baselineWindow(2) - baselineWindow(1);

screenW = 800;
screenH = 600;
x_edges = 0:screenW;   % 800 pixel bins
y_edges = 0:screenH;   % 600 pixel bins
nX = numel(x_edges) - 1;
nY = numel(y_edges) - 1;
nPix = nX * nY;
x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;

blink_win = 25;
smoothing_sigma = 8;          % Gaussian sigma in pixels (sparse sample maps)
support_frac = 0.05;          % keep pixels above this fraction of max mean rate
n_pc_max = 150;               % SVD / GED rank cap
ged_reg = 0.05;               % ridge on both second-moment matrices
min_trials_cond = 20;         % skip subject/condition if fewer valid trials
min_support_pix = 500;

condVars = {'dataET_c25', 'dataET_c50', 'dataET_c75', 'dataET_c100'};
condFileTags = {'c25', 'c50', 'c75', 'c100'};
condLabels = {'25%', '50%', '75%', '100%'};
nCond = numel(condVars);
% Same coding as GCP_gaze_heatmap_CBPT_outline.m (per 25 pp contrast step)
linearWeights = [-0.3, -0.1, 0.1, 0.3];

nSub = numel(subjects);
haufe_pix = nan(nSub, nY, nX);
filter_pix = nan(nSub, nY, nX);
linear_map = nan(nSub, nY, nX);
linear_map_haufe_weighted = nan(nSub, nY, nX);
cond_mean_map = nan(nSub, nCond, nY, nX);
score_stim_cond = nan(nSub, nCond);   % mean stim projection per condition
score_base_cond = nan(nSub, nCond);
score_stim_pool = nan(nSub, 1);
score_base_pool = nan(nSub, 1);
lambda1 = nan(nSub, 1);
n_pc_used = nan(nSub, 1);
n_support = nan(nSub, 1);
n_trials_stim = nan(nSub, nCond);
n_trials_base = nan(nSub, nCond);
valid_subj = false(nSub, 1);

fprintf('Gaze pixel GED: %d subjects, map %d x %d (%d pixels)\n', ...
    nSub, nX, nY, nPix);

%% Subject loop
for subj = 1:nSub
    datapath = fullfile(pathFeat, subjects{subj}, 'gaze');
    etFile = fullfile(datapath, 'dataET.mat');
    if ~isfile(etFile)
        fprintf('[skip] Subject %s: missing dataET\n', subjects{subj});
        continue
    end
    T = load(etFile, condVars{:});
    clc
    fprintf('[GED] Subject %s (%d/%d)\n', subjects{subj}, subj, nSub);

    % --- Pass 1: counts + support (no full map storage) ---
    supportAcc = zeros(nX, nY);
    nPerStim = zeros(1, nCond);
    nPerBase = zeros(1, nCond);
    okCond = true(1, nCond);

    for c = 1:nCond
        if ~isfield(T, condVars{c}) || isempty(T.(condVars{c}))
            okCond(c) = false;
            continue
        end
        [sumStim, nSt] = sumGazeRateMaps( ...
            T.(condVars{c}), stimWindow, stimDur, ...
            x_edges, y_edges, blink_win, screenH, smoothing_sigma);
        [sumBase, nBa] = sumGazeRateMaps( ...
            T.(condVars{c}), baselineWindow, baseDur, ...
            x_edges, y_edges, blink_win, screenH, smoothing_sigma);
        nPerStim(c) = nSt;
        nPerBase(c) = nBa;
        n_trials_stim(subj, c) = nSt;
        n_trials_base(subj, c) = nBa;
        if nSt < min_trials_cond || nBa < min_trials_cond
            okCond(c) = false;
            continue
        end
        supportAcc = supportAcc + sumStim + sumBase;
    end

    nStimAcc = sum(nPerStim);
    nBaseAcc = sum(nPerBase);
    if ~all(okCond) || nStimAcc < min_trials_cond || nBaseAcc < min_trials_cond
        fprintf('[skip] Subject %s: insufficient trials in one or more conditions\n', ...
            subjects{subj});
        continue
    end

    supportMean = supportAcc / max(nStimAcc + nBaseAcc, 1);
    thr = support_frac * max(supportMean(:));
    supportMask = supportMean > thr;   % [nX x nY]
    pixIdx = find(supportMask);
    nSup = numel(pixIdx);
    n_support(subj) = nSup;
    if nSup < min_support_pix
        fprintf('[skip] Subject %s: support too small (%d pixels)\n', ...
            subjects{subj}, nSup);
        continue
    end
    fprintf('  support: %d / %d pixels (%.1f%%), stim trials=%d, base trials=%d\n', ...
        nSup, nPix, 100 * nSup / nPix, nStimAcc, nBaseAcc);

    % --- Pass 2: support-pixel trial vectors only ---
    Xstim = zeros(nStimAcc, nSup, 'single');
    Xbase = zeros(nBaseAcc, nSup, 'single');
    condStimIdx = cell(nCond, 1);
    condBaseIdx = cell(nCond, 1);
    iS = 0;
    iB = 0;
    for c = 1:nCond
        [Xstim_c, nSt] = trialGazeRateVectors( ...
            T.(condVars{c}), stimWindow, stimDur, ...
            x_edges, y_edges, blink_win, screenH, smoothing_sigma, pixIdx);
        [Xbase_c, nBa] = trialGazeRateVectors( ...
            T.(condVars{c}), baselineWindow, baseDur, ...
            x_edges, y_edges, blink_win, screenH, smoothing_sigma, pixIdx);
        if nSt ~= nPerStim(c) || nBa ~= nPerBase(c)
            fprintf('[skip] Subject %s: trial count mismatch on pass 2\n', subjects{subj});
            okCond(c) = false;
            break
        end
        rowsS = iS + (1:nSt);
        rowsB = iB + (1:nBa);
        condStimIdx{c} = rowsS;
        condBaseIdx{c} = rowsB;
        Xstim(rowsS, :) = Xstim_c;
        Xbase(rowsB, :) = Xbase_c;
        iS = iS + nSt;
        iB = iB + nBa;
    end
    if ~all(okCond)
        continue
    end

    % --- Uncentered SVD dimensionality reduction ---
    Xall = [Xstim; Xbase];
    nObs = size(Xall, 1);
    nPC = min([n_pc_max, nObs - 2, nSup]);
    % Truncated SVD of uncentered data (principal axes of spatial energy)
    Xall_d = double(Xall);
    [~, ~, V] = svds(Xall_d, nPC);
    Zstim = double(Xstim) * V;
    Zbase = double(Xbase) * V;
    n_pc_used(subj) = nPC;

    % --- Second-moment GED in PC space ---
    Sstim = (Zstim' * Zstim) / size(Zstim, 1);
    Sbase = (Zbase' * Zbase) / size(Zbase, 1);
    trS = mean(diag(Sstim));
    trB = mean(diag(Sbase));
    Sstim_reg = (1 - ged_reg) * Sstim + ged_reg * trS * eye(nPC);
    Sbase_reg = (1 - ged_reg) * Sbase + ged_reg * trB * eye(nPC);

    [W, D] = eig(Sstim_reg, Sbase_reg);
    [evals, ord] = sort(real(diag(D)), 'descend');
    W = real(W(:, ord));
    w_pc = W(:, 1);
    lambda1(subj) = evals(1);

    % Orient: mean stim score > mean base score
    scoreStimAll = Zstim * w_pc;
    scoreBaseAll = Zbase * w_pc;
    if mean(scoreStimAll) < mean(scoreBaseAll)
        w_pc = -w_pc;
        scoreStimAll = -scoreStimAll;
        scoreBaseAll = -scoreBaseAll;
    end
    score_stim_pool(subj) = mean(scoreStimAll);
    score_base_pool(subj) = mean(scoreBaseAll);

    % Pixel-space filter and Haufe pattern (activation)
    w_sup = V * w_pc;
    a_pc = Sstim_reg * w_pc;
    a_sup = V * a_pc;

    wMap = zeros(nX, nY);
    aMap = zeros(nX, nY);
    wMap(pixIdx) = w_sup;
    aMap(pixIdx) = a_sup;
    % Store as [y x] for plotting (freq=y, time=x)
    filter_pix(subj, :, :) = wMap';
    haufe_pix(subj, :, :) = aMap';

    % Condition mean scores and mean stim maps
    linPix = zeros(nX, nY);
    for c = 1:nCond
        scS = scoreStimAll(condStimIdx{c});
        scB = scoreBaseAll(condBaseIdx{c});
        score_stim_cond(subj, c) = mean(scS);
        score_base_cond(subj, c) = mean(scB);

        meanVec = mean(double(Xstim(condStimIdx{c}, :)), 1);
        mMap = zeros(nX, nY);
        mMap(pixIdx) = meanVec;
        cond_mean_map(subj, c, :, :) = mMap';
        linPix = linPix + linearWeights(c) * mMap;
    end
    linear_map(subj, :, :) = linPix';
    % Emphasize pixels that carry the stim-vs-base GED pattern
    aAbs = abs(aMap);
    if max(aAbs(:)) > 0
        aAbs = aAbs / max(aAbs(:));
    end
    linear_map_haufe_weighted(subj, :, :) = linPix .* aAbs;

    valid_subj(subj) = true;
    fprintf('  lambda1=%.3f | nPC=%d | score stim=%.4g base=%.4g\n', ...
        lambda1(subj), nPC, score_stim_pool(subj), score_base_pool(subj));
end

subjKeep = find(valid_subj);
nKeep = numel(subjKeep);
if nKeep < 2
    error('Fewer than 2 subjects with valid gaze pixel GED.');
end
fprintf('\nValid subjects: %d / %d\n', nKeep, nSub);

%% Group averages (sign-align Haufe to first valid subject)
ref = squeeze(haufe_pix(subjKeep(1), :, :));
haufe_ga = zeros(nY, nX);
filt_ga = zeros(nY, nX);
lin_ga = zeros(nY, nX);
linw_ga = zeros(nY, nX);
cond_ga = zeros(nCond, nY, nX);
sign_flip = ones(nKeep, 1);

for ii = 1:nKeep
    s = subjKeep(ii);
    h = squeeze(haufe_pix(s, :, :));
    if corr(h(:), ref(:), 'rows', 'complete') < 0
        sign_flip(ii) = -1;
        h = -h;
        filter_pix(s, :, :) = -filter_pix(s, :, :);
        haufe_pix(s, :, :) = -haufe_pix(s, :, :);
        score_stim_cond(s, :) = -score_stim_cond(s, :);
        score_base_cond(s, :) = -score_base_cond(s, :);
        score_stim_pool(s) = -score_stim_pool(s);
        score_base_pool(s) = -score_base_pool(s);
    end
    haufe_ga = haufe_ga + h;
    filt_ga = filt_ga + squeeze(filter_pix(s, :, :));
    lin_ga = lin_ga + squeeze(linear_map(s, :, :));
    linw_ga = linw_ga + squeeze(linear_map_haufe_weighted(s, :, :));
    for c = 1:nCond
        cond_ga(c, :, :) = squeeze(cond_ga(c, :, :)) + ...
            squeeze(cond_mean_map(s, c, :, :));
    end
end
haufe_ga = haufe_ga / nKeep;
filt_ga = filt_ga / nKeep;
lin_ga = lin_ga / nKeep;
linw_ga = linw_ga / nKeep;
cond_ga = cond_ga / nKeep;

%% Condition score summaries
scores = score_stim_cond(subjKeep, :);
score_linear = scores * linearWeights(:);
[~, p_lin, ~, stats_lin] = ttest(score_linear);
fprintf('Linear contrast of GED stim scores: t(%d)=%.3f, p=%.4g, mean=%.4g\n', ...
    stats_lin.df, stats_lin.tstat, p_lin, mean(score_linear));

%% Figures
centerX = 400;
centerY = 300;
fontSize = 28;
cmapDiv = customcolormap_preset('red-white-blue');

% 1) Group Haufe pattern (stim-related topography)
plotPixelMap(haufe_ga, x_centers, y_centers, cmapDiv, fontSize, centerX, centerY, ...
    'GED Haufe activation (stim vs base)', 'Activation [a.u.]', ...
    fullfile(figDir, 'GCP_gaze_GED_haufe_GA.png'));

% 2) Group filter
plotPixelMap(filt_ga, x_centers, y_centers, cmapDiv, fontSize, centerX, centerY, ...
    'GED filter (stim vs base)', 'Filter weight [a.u.]', ...
    fullfile(figDir, 'GCP_gaze_GED_filter_GA.png'));

% 3) Condition mean stim maps
for c = 1:nCond
    plotPixelMap(squeeze(cond_ga(c, :, :)), x_centers, y_centers, ...
        parula(256), fontSize, centerX, centerY, ...
        sprintf('Mean stim gaze rate (%s)', condLabels{c}), ...
        'Rate [counts/s]', ...
        fullfile(figDir, sprintf('GCP_gaze_GED_stimmap_%s.png', condFileTags{c})));
end

% 4) Linear contrast of condition mean stim maps (primary "where with contrast")
plotPixelMap(lin_ga, x_centers, y_centers, cmapDiv, fontSize, centerX, centerY, ...
    'Linear contrast of stim gaze maps', ...
    'Density change per 25 pp contrast', ...
    fullfile(figDir, 'GCP_gaze_GED_linear_contrast.png'));

% 5) Same contrast weighted by |Haufe| (GED-relevant pixels)
plotPixelMap(linw_ga, x_centers, y_centers, cmapDiv, fontSize, centerX, centerY, ...
    'Linear contrast x |Haufe|', ...
    'Weighted density change', ...
    fullfile(figDir, 'GCP_gaze_GED_linear_contrast_haufe_weighted.png'));

% 6) Condition GED scores
figure('Position', [0 0 1512 982], 'Color', 'w');
mu = mean(scores, 1);
se = std(scores, 0, 1) / sqrt(nKeep);
hold on
if exist('colors', 'var') && isstruct(colors) && isfield(colors, 'cond')
    barCols = colors.cond;
else
    barCols = lines(nCond);
end
for c = 1:nCond
    bar(c, mu(c), 'FaceColor', barCols(min(c, size(barCols, 1)), :), ...
        'EdgeColor', 'none', 'BarWidth', 0.7);
end
errorbar(1:nCond, mu, se, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 12);
% Subject lines
for ii = 1:nKeep
    plot(1:nCond, scores(ii, :), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
end
plot(1:nCond, mu, 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
set(gca, 'XTick', 1:nCond, 'XTickLabel', condLabels, 'FontSize', fontSize);
ylabel('GED stim score [a.u.]', 'FontSize', fontSize);
xlabel('Stimulus contrast', 'FontSize', fontSize);
title(sprintf('Condition GED scores (linear t=%.2f, p=%.3g)', ...
    stats_lin.tstat, p_lin), 'FontSize', fontSize);
box off
exportgraphics(gcf, fullfile(figDir, 'GCP_gaze_GED_condition_scores.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');

% 7) Stim vs base pooled scores
figure('Position', [0 0 1512 982], 'Color', 'w');
sb = [score_stim_pool(subjKeep), score_base_pool(subjKeep)];
mu_sb = mean(sb, 1);
se_sb = std(sb, 0, 1) / sqrt(nKeep);
bar(1:2, mu_sb, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none'); hold on
errorbar(1:2, mu_sb, se_sb, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 12);
for ii = 1:nKeep
    plot(1:2, sb(ii, :), '-o', 'Color', [0.6 0.6 0.6], 'MarkerSize', 4);
end
set(gca, 'XTick', 1:2, 'XTickLabel', {'Stimulus', 'Baseline'}, 'FontSize', fontSize);
ylabel('GED score [a.u.]', 'FontSize', fontSize);
title('Pooled GED projection (filter check)', 'FontSize', fontSize);
box off
exportgraphics(gcf, fullfile(figDir, 'GCP_gaze_GED_stim_vs_base_scores.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');

%% Save
outMat = fullfile(pathFeat, 'GCP_gaze_heatmap_GED.mat');
save(outMat, ...
    'subjects', 'subjKeep', 'valid_subj', ...
    'haufe_pix', 'filter_pix', 'linear_map', 'linear_map_haufe_weighted', ...
    'cond_mean_map', 'score_stim_cond', 'score_base_cond', ...
    'score_stim_pool', 'score_base_pool', 'score_linear', ...
    'lambda1', 'n_pc_used', 'n_support', 'n_trials_stim', 'n_trials_base', ...
    'haufe_ga', 'filt_ga', 'lin_ga', 'linw_ga', 'cond_ga', ...
    'stimWindow', 'baselineWindow', 'linearWeights', ...
    'smoothing_sigma', 'support_frac', 'n_pc_max', 'ged_reg', ...
    'x_centers', 'y_centers', 'nX', 'nY', ...
    '-v7.3');
fprintf('Saved: %s\n', outMat);
fprintf('Figures: %s\n', figDir);
fprintf('DONE.\n');

%% Local functions
function [x, y] = cleanGazeTrial(trialMat, blink_win, screenH)
x = [];
y = [];
data = double(trialMat);
if size(data, 1) < 3 || isempty(data)
    return
end
valid_tr = data(1, :) >= 0 & data(1, :) <= 800 & ...
           data(2, :) >= 0 & data(2, :) <= screenH;
data = data(1:3, valid_tr);
if isempty(data)
    return
end
data(2, :) = screenH - data(2, :);
data = remove_blinks(data, blink_win);
x_positions = data(1, :);
y_positions = data(2, :);
fin = isfinite(x_positions) & isfinite(y_positions);
x_positions = x_positions(fin);
y_positions = y_positions(fin);
ok = ~(x_positions == 0 & y_positions == 0);
x = x_positions(ok);
y = y_positions(ok);
end

function [sumMap, nKeep] = sumGazeRateMaps(dataET, latencyWindow, winDur, ...
    x_edges, y_edges, blink_win, screenH, smoothing_sigma)
% Sum of trial rate maps [nX x nY] for support estimation.

cfg = [];
cfg.avgovertime = 'no';
cfg.keeptrials = 'yes';
cfg.latency = latencyWindow;
dataSel = ft_selectdata(cfg, dataET);

nX = numel(x_edges) - 1;
nY = numel(y_edges) - 1;
sumMap = zeros(nX, nY);
nKeep = 0;
for tr = 1:numel(dataSel.trial)
    [x_positions, y_positions] = cleanGazeTrial(dataSel.trial{tr}, blink_win, screenH);
    if numel(x_positions) < 5
        continue
    end
    binned = histcounts2(x_positions, y_positions, x_edges, y_edges);
    if smoothing_sigma > 0
        binned = imgaussfilt(binned, smoothing_sigma);
    end
    sumMap = sumMap + binned / winDur;
    nKeep = nKeep + 1;
end
end

function [X, nKeep] = trialGazeRateVectors(dataET, latencyWindow, winDur, ...
    x_edges, y_edges, blink_win, screenH, smoothing_sigma, pixIdx)
% X: [nTrials x nSupport] duration-normalized, smoothed gaze rate on support pixels.

cfg = [];
cfg.avgovertime = 'no';
cfg.keeptrials = 'yes';
cfg.latency = latencyWindow;
dataSel = ft_selectdata(cfg, dataET);

nSup = numel(pixIdx);
nTr = numel(dataSel.trial);
X = zeros(nTr, nSup, 'single');
keep = false(nTr, 1);

for tr = 1:nTr
    [x_positions, y_positions] = cleanGazeTrial(dataSel.trial{tr}, blink_win, screenH);
    if numel(x_positions) < 5
        continue
    end
    binned = histcounts2(x_positions, y_positions, x_edges, y_edges); % [nX x nY]
    if smoothing_sigma > 0
        binned = imgaussfilt(binned, smoothing_sigma);
    end
    rate = binned / winDur;
    X(tr, :) = single(rate(pixIdx)');
    keep(tr) = true;
end

X = X(keep, :);
nKeep = size(X, 1);
end

function plotPixelMap(mapYX, x_centers, y_centers, cmap, fontSize, centerX, centerY, ...
    titleStr, cbLabel, outPath)
% mapYX: [nY x nX]

freq = [];
freq.label = {'et'};
freq.dimord = 'chan_freq_time';
freq.time = x_centers;
freq.freq = y_centers;
freq.powspctrm = zeros(1, numel(y_centers), numel(x_centers));
freq.powspctrm(1, :, :) = mapYX;

vals = mapYX(isfinite(mapYX));
if isempty(vals)
    zlimVals = [-1 1];
elseif all(vals >= 0)
    zlimVals = [0, max(prctile(vals, 99.5), eps)];
else
    lim = max(prctile(abs(vals), 99.5), eps);
    zlimVals = [-lim lim];
end

figure('Position', [0 0 1512 982], 'Color', 'w');
cfg = [];
cfg.figure = 'gcf';
cfg.zlim = zlimVals;
cfg.colormap = cmap;
ft_singleplotTFR(cfg, freq);
title(titleStr, 'FontSize', fontSize);
xlim([0 800]);
ylim([0 600]);
yticks([0 150 300 450 600]);
set(gca, 'FontSize', fontSize);
xlabel('Screen Width [px]', 'FontSize', fontSize);
ylabel('Screen Height [px]', 'FontSize', fontSize);
cb = colorbar;
set(cb, 'FontSize', fontSize);
ylabel(cb, cbLabel, 'FontSize', fontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
exportgraphics(gcf, outPath, 'Resolution', 300, 'BackgroundColor', 'white');
end
