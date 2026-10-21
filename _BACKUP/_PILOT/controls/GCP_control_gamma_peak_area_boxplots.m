%% Boxplots of gamma peak frequency and power from the 50% area frequency
% Same subject-level layout as 4_stats/GCP_stats_boxplots.m.
% Frequency is the point at which half the area under the peak has
% accumulated. Power is the mean of the condition spectrum within +/- 5 Hz
% of that frequency, matching pick_tallest_peak.

startup
[subjects, paths, colors, ~] = setup('GCP', 0);
if ispc
    addpath('C:\Users\Administrator\Documents\GitHub\GCP\_BACKUP\_PILOT');
else
    addpath('/Users/Arne/Documents/GitHub/GCP/_BACKUP/_PILOT');
end
[paths, subjects] = gcp_pilot_paths(paths);
try
    subjects = gcp_subject_inclusion(subjects, paths);
catch
    subjects = subjects(~ismember(str2double(string(subjects)), [602, 604, 608]));
end
nSubj = numel(subjects);
nCond = 4;
fprintf('[AREA BOX] N = %d (%s)\n', nSubj, strjoin(subjects, ', '));

ged = load(fullfile(paths.features, 'GCP_eeg_GED.mat'), ...
    'all_condition_powspctrm', 'scan_freqs', 'subjects');
ged_idx = nan(1, nSubj);
for i = 1:nSubj
    ged_idx(i) = find(strcmp(ged.subjects, subjects{i}), 1, 'first');
end
x = ged.scan_freqs(:)';

freq = nan(nCond, nSubj);
pow = nan(nCond, nSubj);
for s = 1:nSubj
    for c = 1:nCond
        y = ged.all_condition_powspctrm{c, ged_idx(s)};
        if isempty(y) || numel(y) ~= numel(x)
            continue
        end
        f_area = peak_half_area(x, y(:)');
        freq(c, s) = f_area;
        if isfinite(f_area)
            band = abs(x - f_area) <= 5 & isfinite(y(:)');
            pow(c, s) = mean(y(band), 'omitnan');
        end
    end
end

out_dir = fullfile(paths.figures, 'stats', 'boxplots');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
ctrl_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/_PILOT/figures/controls';
if ~exist(ctrl_dir, 'dir')
    mkdir(ctrl_dir);
end

specs = {
    'Frequency_area50', freq, 'Frequency [Hz]', false
    'Power_area50',     pow,  'Power [dB]',     true
    };
plot_area_boxplots(specs, subjects, colors, out_dir, ctrl_dir);

save(fullfile(ctrl_dir, 'GCP_control_gamma_peak_area_boxplots.mat'), ...
    'subjects', 'freq', 'pow', '-v7.3');
fprintf('[AREA BOX] Frequency medians: %s\n', mat2str(median(freq, 2, 'omitnan'), 4));
fprintf('[AREA BOX] Power medians: %s\n', mat2str(median(pow, 2, 'omitnan'), 4));

function plot_area_boxplots(specs, subjects, colors, out_dir, ctrl_dir)
nCond = 4;
nSubj = numel(subjects);
xtickLabs = {'25%', '50%', '75%', '100%'};
fontSize = 50;
axisLabelSize = fontSize * 0.9;
yTickFontSize = axisLabelSize;
xLabelFontSize = axisLabelSize;
xTickFontSize = axisLabelSize * (1 - 0.33);
dotSize = 400;
dotAlpha = 0.85;
jitter = 0.4;
boxWidth = 0.55;
subjIDs = str2double(string(subjects(:)));

for iMetric = 1:size(specs, 1)
    varName = specs{iMetric, 1};
    M = specs{iMetric, 2};
    ylabStr = specs{iMetric, 3};
    drawZero = specs{iMetric, 4};

    y = nan(nCond * nSubj, 1);
    condIdx = nan(nCond * nSubj, 1);
    idVec = nan(nCond * nSubj, 1);
    row = 0;
    for s = 1:nSubj
        for c = 1:nCond
            row = row + 1;
            idVec(row) = subjIDs(s);
            condIdx(row) = c;
            y(row) = M(c, s);
        end
    end

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
    style_boxes(colors(1:nCond, :));

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

    outName = sprintf('GCP_stats_boxplot_%s.png', varName);
    drawnow; pause(0.05);
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, fullfile(out_dir, outName), '-dpng', '-r600');
    print(gcf, fullfile(ctrl_dir, outName), '-dpng', '-r600');
    fprintf('[AREA BOX] Saved %s\n', fullfile(out_dir, outName));
end
end

function style_boxes(boxColors)
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

function f_area = peak_half_area(x, y)
f_area = NaN;
x = x(:)';
y = y(:)';
mask = x >= 30 & x <= 90 & isfinite(x) & isfinite(y);
if nnz(mask) < 9
    return
end
xx = x(mask);
yy = y(mask);
ys = movmean(yy, 5, 'omitnan');
[~, ipeak] = max(ys);
ilo = ipeak;
while ilo > 2 && ys(ilo) >= ys(ilo - 1)
    ilo = ilo - 1;
end
ihi = ipeak;
while ihi < numel(ys) - 1 && ys(ihi) >= ys(ihi + 1)
    ihi = ihi + 1;
end
if ihi - ilo < 4
    return
end
xseg = xx(ilo:ihi);
yseg = ys(ilo:ihi);
chord = linspace(yseg(1), yseg(end), numel(yseg));
above = max(yseg - chord, 0);
area = trapz(xseg, above);
if ~isfinite(area) || area <= 0
    return
end
c = cumtrapz(xseg, above);
k = find(c >= 0.5 * area, 1, 'first');
if isempty(k)
    return
end
if k == 1
    f_area = xseg(1);
else
    f_area = xseg(k - 1) + (xseg(k) - xseg(k - 1)) * ...
        ((0.5 * area - c(k - 1)) / max(c(k) - c(k - 1), eps));
end
end
