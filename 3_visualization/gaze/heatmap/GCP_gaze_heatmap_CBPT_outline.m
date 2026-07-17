%% GCP Gaze Heatmap: baselined conditions, linear contrast, and CBPT
% Builds duration-normalized gaze density heatmaps for stimulus [0 2] and
% baseline [-1.5 -0.25] (GCP EEG baseline). Per-bin percentage change is
% unstable for sparse fixation maps, so baselining uses a global denominator:
%   bl = 100 * (stim_rate - base_rate) / mean(base_rate > 0)
% i.e. density change as % of mean baseline density. Plots all contrasts,
% then a linear contrast and CBPT (100% vs 25%) restricted to bins with
% sufficient gaze support. The CBPT figure shows the grand-average
% difference with a significant-cluster outline.
% Gaze preprocessing matches GCP_gaze_fex.m: in-bounds mask on raw
% coordinates, Y -> screen coordinates via (600 - y), remove_blinks
% (25 samples).

%% Setup
startup
[subjects, paths, ~, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
path = paths.features;
figDir = fullfile(paths.figures, 'gaze', 'heatmap');
if ~isfolder(figDir), mkdir(figDir); end

stimWindow = [0 2];
baselineWindow = [-1.5 -0.25];
stimDur = stimWindow(2) - stimWindow(1);
baseDur = baselineWindow(2) - baselineWindow(1);

num_bins = 50;
smoothing_factor = 2;
x_edges = linspace(0, 800, num_bins);
y_edges = linspace(0, 600, num_bins);
x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;

blink_win = 25;
screenH = 600;

condVars = {'dataET_c25', 'dataET_c50', 'dataET_c75', 'dataET_c100'};
condFileTags = {'c25', 'c50', 'c75', 'c100'};
nCond = numel(condVars);

%% Subject-level baselined heatmaps
nSub = length(subjects);
dataBlAll = cell(nCond, nSub);
rateSupportAll = cell(nCond, nSub);

for subj = 1:nSub
    datapath = fullfile(path, subjects{subj}, 'gaze');
    T = load(fullfile(datapath, 'dataET'), condVars{:});
    clc
    disp(upper(['Loading ET data for subject ' num2str(subj) '/' num2str(nSub) '...']))

    for c = 1:nCond
        hmStim = buildGazeHeatmap(T.(condVars{c}), stimWindow, ...
            x_edges, y_edges, blink_win, screenH, smoothing_factor);
        hmBase = buildGazeHeatmap(T.(condVars{c}), baselineWindow, ...
            x_edges, y_edges, blink_win, screenH, smoothing_factor);

        rateStim = hmStim / stimDur;
        rateBase = hmBase / baseDur;
        rateSupportAll{c, subj} = rateStim + rateBase;

        posBase = rateBase(rateBase > 0);
        if isempty(posBase)
            blMap = zeros(size(rateStim));
        else
            meanBase = mean(posBase);
            blMap = 100 * (rateStim - rateBase) / meanBase;
        end

        freq = [];
        freq.label = {'et'};
        freq.dimord = 'chan_freq_time';
        freq.time = x_centers;
        freq.freq = y_centers;
        freq.powspctrm(1, :, :) = blMap.';
        dataBlAll{c, subj} = freq;
    end
end

%% Keep subjects with all four conditions
valid = true(1, nSub);
for c = 1:nCond
    valid = valid & ~cellfun(@isempty, dataBlAll(c, :));
end
dataBlAll = dataBlAll(:, valid);
rateSupportAll = rateSupportAll(:, valid);
nSub = sum(valid);
if nSub < 2
    error('Fewer than 2 subjects with valid baselined gaze heatmaps.');
end
fprintf('Sample size: %d subjects\n', nSub);

%% Gaze-support mask (exclude sparse periphery from CBPT)
% Keep bins where the group-mean rate (stim + baseline, all contrasts)
% exceeds 5% of the map maximum. Unsupported bins are set to NaN before
% ft_freqstatistics so empty periphery cannot inflate variance or clusters.
supportMean = zeros(size(rateSupportAll{1, 1}));
for c = 1:nCond
    for s = 1:nSub
        supportMean = supportMean + rateSupportAll{c, s};
    end
end
supportMean = supportMean / (nCond * nSub);
supportThr = 0.05 * max(supportMean(:));
gazeSupport = supportMean > supportThr;
gazeSupportFt = gazeSupport.';
nBinTot = numel(gazeSupport);
nBinKeep = nnz(gazeSupport);
fprintf('Gaze support: %d / %d bins (%.1f%%), threshold = %.5g\n', ...
    nBinKeep, nBinTot, 100 * nBinKeep / nBinTot, supportThr);
if nBinKeep < 10
    error('Gaze support mask too small (%d bins). Check heatmap / threshold.', nBinKeep);
end

%% Grand averages per contrast
datGA = cell(1, nCond);
for c = 1:nCond
    datGA{c} = ft_freqgrandaverage([], dataBlAll{c, :});
end

%% Plot baselined heatmaps for all contrasts
close all
overallFontSize = 50;
centerX = 400;
centerY = 300;
colMapBl = customcolormap_preset('red-white-blue');
powPool = [];
for c = 1:nCond
    powPool = [powPool; datGA{c}.powspctrm(:)]; %#ok<AGROW>
end
robustLim = prctile(abs(powPool(isfinite(powPool))), 99.5);
robustLim = max(robustLim, 1);

for c = 1:nCond
    plotBlDensityMap(datGA{c}, colMapBl, [-robustLim robustLim], overallFontSize, ...
        centerX, centerY, ...
        fullfile(figDir, sprintf('GCP_gaze_heatmap_bl_%s.png', condFileTags{c})));
end

%% Linear contrast across 25%, 50%, 75%, and 100%
% Regression slope per one condition step, corresponding to a 25 percentage
% point increase in stimulus contrast.
linearWeights = [-0.3 -0.1 0.1 0.3];
linearGA = datGA{1};
linearGA.powspctrm = zeros(size(datGA{1}.powspctrm));
for c = 1:nCond
    linearGA.powspctrm = linearGA.powspctrm + ...
        linearWeights(c) * datGA{c}.powspctrm;
end
linearLim = prctileFinite(abs(linearGA.powspctrm(:)), 99.5);
linearLim = max(linearLim, 1);
plotBlDensityMap(linearGA, colMapBl, [-linearLim linearLim], overallFontSize, ...
    centerX, centerY, ...
    fullfile(figDir, 'GCP_gaze_heatmap_bl_linear_contrast.png'));

fprintf('Linear contrast weights: [%.1f %.1f %.1f %.1f]\n', linearWeights);
fprintf('Linear map units: gaze density change per 25 percentage point contrast increase.\n');

%% Baselined difference: 100% vs 25%
diffGA = datGA{4};
diffGA.powspctrm = datGA{4}.powspctrm - datGA{1}.powspctrm;
diffLim = prctile(abs(diffGA.powspctrm(:)), 99.5);
diffLim = max(diffLim, 1);
plotBlDensityMap(diffGA, colMapBl, [-diffLim diffLim], overallFontSize, ...
    centerX, centerY, ...
    fullfile(figDir, 'GCP_gaze_heatmap_bl_diff_c100_c25.png'));

%% CBPT: 100% vs 25% on gaze-supported bins only
dataCbpt = cell(size(dataBlAll));
for c = [1 4]
    for s = 1:nSub
        tmp = dataBlAll{c, s};
        P = squeeze(tmp.powspctrm);
        P(~gazeSupportFt) = NaN;
        tmp.powspctrm(1, :, :) = P;
        dataCbpt{c, s} = tmp;
    end
end

cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 'all';
cfg.neighbours = [];
cfg.minnbchan = 0;

design = zeros(2, 2 * nSub);
for i = 1:nSub
    design(1, i) = i;
end
for i = 1:nSub
    design(1, nSub + i) = i;
end
design(2, 1:nSub) = 1;
design(2, nSub + 1:2 * nSub) = 2;

cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;

clc
fprintf('\n========== CBPT: 100%% vs 25%% (baselined, gaze-supported bins) ==========\n');
fprintf('Subjects: %d\n', nSub);
fprintf('Supported bins: %d / %d\n', nBinKeep, nBinTot);
fprintf('Statistic: %s | Correction: %s | Randomizations: all\n', ...
    cfg.statistic, cfg.correctm);
fprintf('clusteralpha = %.3f | alpha = %.3f | tail = %d\n', ...
    cfg.clusteralpha, cfg.alpha, cfg.tail);
disp('Computing ft_freqstatistics...')
statDIFF = ft_freqstatistics(cfg, dataCbpt{4, :}, dataCbpt{1, :});

%% Effect-size summary (supported bins); plot uses grand-average difference
statD = statDIFF;
statD.stat = statDIFF.stat ./ sqrt(nSub);
reportCbptSummary(statDIFF, statD, cfg.alpha);

%% Plot grand-average difference with significant-cluster outline
close all
diffPlot = diffGA;
Pdiff = squeeze(diffPlot.powspctrm);
diffLimCbpt = prctileFinite(abs(Pdiff(gazeSupportFt)), 99.5);
diffLimCbpt = max(diffLimCbpt, 1);
Pdiff(~gazeSupportFt) = 0;
diffPlot.powspctrm(1, :, :) = Pdiff;
if isfield(statDIFF, 'mask')
    diffPlot.mask = statDIFF.mask;
else
    diffPlot.mask = false(size(diffPlot.powspctrm));
end
diffColMap = customcolormap_preset('red-white-blue');
outCbpt = fullfile(figDir, 'GCP_gaze_heatmap_bl_CBPT_outline_diff.png');
plotDiffWithClusterOutline(diffPlot, diffColMap, [-diffLimCbpt diffLimCbpt], ...
    overallFontSize, centerX, centerY, outCbpt);
fprintf('Saved: %s\n', outCbpt);
fprintf('========== CBPT done ==========\n\n');

%%
function hm = buildGazeHeatmap(dataET, latencyWindow, x_edges, y_edges, blink_win, screenH, smoothing_factor)
cfg = [];
cfg.avgovertime = 'no';
cfg.keeptrials = 'yes';
cfg.latency = latencyWindow;
dataSel = ft_selectdata(cfg, dataET);

parts_x = cell(1, numel(dataSel.trial));
parts_y = cell(1, numel(dataSel.trial));
for tr = 1:numel(dataSel.trial)
    data = double(dataSel.trial{tr});
    if size(data, 1) < 3 || isempty(data)
        continue
    end
    valid_tr = data(1, :) >= 0 & data(1, :) <= 800 & ...
               data(2, :) >= 0 & data(2, :) <= screenH;
    data = data(1:3, valid_tr);
    if isempty(data)
        continue
    end
    data(2, :) = screenH - data(2, :);
    data = remove_blinks(data, blink_win);
    x_positions = data(1, :);
    y_positions = data(2, :);
    fin = isfinite(x_positions) & isfinite(y_positions);
    x_positions = x_positions(fin);
    y_positions = y_positions(fin);
    ok = ~(x_positions == 0 & y_positions == 0);
    parts_x{tr} = x_positions(ok);
    parts_y{tr} = y_positions(ok);
end
x_positions = [parts_x{:}];
y_positions = [parts_y{:}];
binned_data = histcounts2(x_positions, y_positions, x_edges, y_edges);
hm = imgaussfilt(binned_data, smoothing_factor);
end

function plotBlDensityMap(freqData, cmap, zlimVals, fontSize, centerX, centerY, outPath)
figure('Position', [0 0 1512 982], 'Color', 'w');
cfg = [];
cfg.figure = 'gcf';
cfg.zlim = zlimVals;
cfg.colormap = cmap;
ft_singleplotTFR(cfg, freqData);
title('');
xlim([0 800]);
ylim([0 600]);
yticks([0 150 300 450 600]);
set(gca, 'FontSize', fontSize);
xlabel('Screen Width [px]', 'FontSize', fontSize);
ylabel('Screen Height [px]', 'FontSize', fontSize);
cb = colorbar;
set(cb, 'FontSize', fontSize);
ylabel(cb, 'Gaze Density [%]', 'FontSize', fontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
fitHeatmapLayout(gca, cb, fontSize);
exportgraphics(gcf, outPath, 'Resolution', 300, 'BackgroundColor', 'white');
end

function lim = prctileFinite(x, p)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    lim = 1;
elseif numel(x) < 5
    lim = max(x);
else
    lim = prctile(x, p);
end
end

function reportCbptSummary(statT, statD, alphaThr)
t = statT.stat(:);
t = t(isfinite(t));
d = statD.stat(:);
d = d(isfinite(d));
fprintf('\n--- Observed effects (gaze-supported bins, unthresholded) ---\n');
fprintf('t:  min = %8.3f | max = %8.3f | mean |t| = %.3f\n', ...
    min(t), max(t), mean(abs(t)));
fprintf('d:  min = %8.3f | max = %8.3f | mean |d| = %.3f\n', ...
    min(d), max(d), mean(abs(d)));

if isfield(statT, 'mask') && ~isempty(statT.mask)
    nMask = nnz(statT.mask);
    nFinite = nnz(isfinite(statT.stat));
    fprintf('Significant mask voxels: %d / %d finite supported (%.2f%%)\n', ...
        nMask, nFinite, 100 * nMask / max(nFinite, 1));
else
    fprintf('No mask field on stat structure.\n');
    nMask = 0;
end

printClusterSide('Positive', getfield_or(statT, 'posclusters'), alphaThr);
printClusterSide('Negative', getfield_or(statT, 'negclusters'), alphaThr);

if nMask == 0
    fprintf('RESULT: no cluster survived correction (outline will be empty).\n');
    fprintf('Grand-average difference is still plotted on supported bins.\n');
else
    fprintf('RESULT: at least one significant cluster (outline drawn).\n');
end
fprintf('---------------------------------------\n');
end

function v = getfield_or(s, name)
if isfield(s, name)
    v = s.(name);
else
    v = [];
end
end

function printClusterSide(sideLabel, clusters, alphaThr)
if isempty(clusters)
    fprintf('%s clusters: none formed\n', sideLabel);
    return
end
nCl = numel(clusters);
fprintf('%s clusters formed: %d\n', sideLabel, nCl);
nShow = min(nCl, 10);
for i = 1:nShow
    cl = clusters(i);
    prob = NaN;
    cstat = NaN;
    if isfield(cl, 'prob'), prob = cl.prob; end
    if isfield(cl, 'clusterstat'), cstat = cl.clusterstat; end
    sigTag = '';
    if isfinite(prob) && prob < alphaThr
        sigTag = '  <-- significant';
    end
    fprintf('  #%d  clusterstat = %10.3f | prob = %.4f%s\n', ...
        i, cstat, prob, sigTag);
end
if nCl > nShow
    fprintf('  ... %d more not shown\n', nCl - nShow);
end
end

function plotDiffWithClusterOutline(freqData, cmap, zlimVals, fontSize, centerX, centerY, outPath)
figure('Position', [0 0 1512 982], 'Color', 'w');
cfg = [];
cfg.figure = 'gcf';
cfg.parameter = 'powspctrm';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.interactivecolor = [0 0 0];
cfg.zlim = zlimVals;
cfg.colormap = cmap;
ft_singleplotTFR(cfg, freqData);
title('');
xlim([0 800]);
ylim([0 600]);
yticks([0 150 300 450 600]);
set(gca, 'FontSize', fontSize);
xlabel('Screen Width [px]', 'FontSize', fontSize);
ylabel('Screen Height [px]', 'FontSize', fontSize);
cb = colorbar;
set(cb, 'FontSize', fontSize);
ylabel(cb, 'Gaze Density [%]', 'FontSize', fontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
fitHeatmapLayout(gca, cb, fontSize);
exportgraphics(gcf, outPath, 'Resolution', 300, 'BackgroundColor', 'white');
end

function fitHeatmapLayout(ax, cb, fontSize)
% Manual layout: TightInset/LooseInset do not reserve space for axis/colorbar labels.
s = fontSize / 50;
left = 0.15 + 0.06 * s;
bottom = 0.20 + 0.08 * s;
topPad = 0.06 + 0.03 * s;
cbW = 0.02;
rightPad = 0.06 + 0.10 * s;
cbX = 1 - rightPad - cbW;
axW = cbX - left - 0.015;
axH = 1 - bottom - topPad;
ax.Position = [left bottom axW axH];
cb.Position = [cbX bottom cbW axH];
end
