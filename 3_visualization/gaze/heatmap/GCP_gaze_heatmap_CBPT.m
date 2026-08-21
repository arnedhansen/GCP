%% GCP Gaze Heatmap: each contrast vs baseline
% Duration-normalized gaze density maps for stimulus [0 2] and baseline
% [-1.5 -0.5]. Baselining uses a global denominator:
%   bl = 100 * (stim_rate - base_rate) / mean(base_rate > 0)
% Gaze preprocessing matches GCP_gaze_fex.m.

%% Setup
startup
[subjects, paths, ~, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
path = paths.features;
figDir = fullfile(paths.figures, 'gaze', 'heatmap');
if ~isfolder(figDir)
    mkdir(figDir);
end

stimWindow = [0 2];
baselineWindow = [-1.5 -0.5];
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
condLabels = {'25% Contrast', '50% Contrast', '75% Contrast', '100% Contrast'};
nCond = numel(condVars);

%% Subject-level stim, baseline, and baselined heatmaps
nSub = length(subjects);
dataStimAll = cell(nCond, nSub);
dataBaseAll = cell(nCond, nSub);
dataBlAll = cell(nCond, nSub);

for subj = 1:nSub
    datapath = fullfile(path, subjects{subj}, 'gaze');
    T = load(fullfile(datapath, 'dataET'), condVars{:});
    clc; fprintf('[VIZ GAZE HEATMAP] Subject %d/%d (%s)\n', subj, nSub, subjects{subj})

    for c = 1:nCond
        hmStim = buildGazeHeatmap(T.(condVars{c}), stimWindow, ...
            x_edges, y_edges, blink_win, screenH, smoothing_factor);
        hmBase = buildGazeHeatmap(T.(condVars{c}), baselineWindow, ...
            x_edges, y_edges, blink_win, screenH, smoothing_factor);

        rateStim = hmStim / stimDur;
        rateBase = hmBase / baseDur;

        posBase = rateBase(rateBase > 0);
        if isempty(posBase)
            blMap = zeros(size(rateStim));
        else
            meanBase = mean(posBase);
            blMap = 100 * (rateStim - rateBase) / meanBase;
        end

        dataStimAll{c, subj} = mapToFreq(rateStim, x_centers, y_centers);
        dataBaseAll{c, subj} = mapToFreq(rateBase, x_centers, y_centers);
        dataBlAll{c, subj} = mapToFreq(blMap, x_centers, y_centers);
    end
end

%% Grand averages per contrast
datGA = cell(1, nCond);
for c = 1:nCond
    datGA{c} = ft_freqgrandaverage([], dataBlAll{c, :});
end

%% CBPT: each contrast vs baseline
%{
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
design(1, 1:nSub) = 1:nSub;
design(1, nSub + 1:2 * nSub) = 1:nSub;
design(2, 1:nSub) = 1;
design(2, nSub + 1:2 * nSub) = 2;
cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;

statAll = cell(1, nCond);
for c = 1:nCond
    dataStim = cell(1, nSub);
    dataBase = cell(1, nSub);
    for s = 1:nSub
        dataStim{s} = dataStimAll{c, s};
        dataBase{s} = dataBaseAll{c, s};
    end

    clc
    fprintf('\n========== CBPT: %s vs baseline ==========\n', condLabels{c});
    fprintf('Subjects: %d\n', nSub);
    fprintf('Statistic: %s | Correction: %s | Randomizations: all\n', ...
        cfg.statistic, cfg.correctm);
    fprintf('clusteralpha = %.3f | alpha = %.3f | tail = %d\n', ...
        cfg.clusteralpha, cfg.alpha, cfg.tail);
    disp('Computing ft_freqstatistics...')
    statAll{c} = ft_freqstatistics(cfg, dataStim{:}, dataBase{:});

    statD = statAll{c};
    statD.stat = statAll{c}.stat ./ sqrt(nSub);
    reportCbptSummary(statAll{c}, statD, cfg.alpha, condLabels{c});
end
%}

%% Plot 2x2 grand-average percentage change
close all
overallFontSize = 20;
centerX = 400;
centerY = 300;
colMapBl = customcolormap_preset('red-white-blue');

powPool = [];
for c = 1:nCond
    powPool = [powPool; datGA{c}.powspctrm(:)]; %#ok<AGROW>
end
powPool = abs(powPool(isfinite(powPool)));
robustLim = prctile(powPool, 99.5);
if ~isfinite(robustLim) || robustLim <= 0
    robustLim = 1;
end
zlimVals = [-robustLim robustLim];

figure('Position', [0 0 1512 982], 'Color', 'w');
for c = 1:nCond
    subplot(2, 2, c);
    plotDat = datGA{c};
    % if exist('statAll', 'var') && isfield(statAll{c}, 'mask') && ~isempty(statAll{c}.mask)
    %     plotDat.mask = statAll{c}.mask;
    % end
    plotHeatmapPanel(plotDat, colMapBl, zlimVals, ...
        overallFontSize, centerX, centerY, condLabels{c});
end

exportgraphics(gcf, fullfile(figDir, 'GCP_gaze_heatmap_CBPT.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');

%%
function freq = mapToFreq(mapXY, x_centers, y_centers)
freq = [];
freq.label = {'et'};
freq.dimord = 'chan_freq_time';
freq.time = x_centers;
freq.freq = y_centers;
freq.powspctrm = zeros(1, numel(y_centers), numel(x_centers));
freq.powspctrm(1, :, :) = mapXY.';
end

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

function plotHeatmapPanel(freqData, cmap, zlimVals, fontSize, centerX, centerY, titleStr)
cfg = [];
cfg.figure = 'gcf';
cfg.parameter = 'powspctrm';
if isfield(freqData, 'mask')
    cfg.maskparameter = 'mask';
    cfg.maskstyle = 'outline';
    cfg.interactivecolor = [0 0 0];
end
cfg.zlim = zlimVals;
cfg.colormap = cmap;
ft_singleplotTFR(cfg, freqData);
title(titleStr, 'FontSize', fontSize);
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
plot(centerX, centerY, '+', 'MarkerSize', 12, 'LineWidth', 2, 'Color', 'k');
end

function reportCbptSummary(statT, statD, alphaThr, condLabel)
t = statT.stat(:);
t = t(isfinite(t));
d = statD.stat(:);
d = d(isfinite(d));
fprintf('\n--- Observed effects (%s, unthresholded) ---\n', condLabel);
fprintf('t:  min = %8.3f | max = %8.3f | mean |t| = %.3f\n', ...
    min(t), max(t), mean(abs(t)));
fprintf('d:  min = %8.3f | max = %8.3f | mean |d| = %.3f\n', ...
    min(d), max(d), mean(abs(d)));

if isfield(statT, 'mask') && ~isempty(statT.mask)
    nMask = nnz(statT.mask);
    nFinite = nnz(isfinite(statT.stat));
    fprintf('Significant mask voxels: %d / %d finite (%.2f%%)\n', ...
        nMask, nFinite, 100 * nMask / max(nFinite, 1));
else
    fprintf('No mask field on stat structure.\n');
    nMask = 0;
end

printClusterSide('Positive', getfield_or(statT, 'posclusters'), alphaThr);
printClusterSide('Negative', getfield_or(statT, 'negclusters'), alphaThr);

if nMask == 0
    fprintf('RESULT: no cluster survived correction (outline will be empty).\n');
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
