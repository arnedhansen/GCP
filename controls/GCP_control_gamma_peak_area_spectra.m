%% 50% area frequency, one spectrum per subject
% Each panel is the 100% contrast GED stimulus/baseline spectrum.
% The peak is the tallest point of a 5-point moving mean (5 Hz on the
% 1 Hz grid). Its bounds are the flanking minima of that smoothed curve.
% Area is the region above the straight line joining those minima.
% The marked frequency is where that area reaches 50%.

startup
[subjects, paths, ~, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);

out_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/controls';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
set(0, 'DefaultFigureVisible', 'off');

ged = load(fullfile(paths.features, 'GCP_eeg_GED.mat'), ...
    'all_condition_powspctrm', 'scan_freqs', 'subjects');
idx = nan(1, numel(subjects));
for i = 1:numel(subjects)
    idx(i) = find(strcmp(ged.subjects, subjects{i}), 1, 'first');
end
x = ged.scan_freqs(:)';
nSubj = numel(subjects);
cond = 4; % 100% contrast
if nSubj ~= 10
    warning('[AREA SPEC] Expected 10 subjects, found %d. Layout stays 2 x 5.', nSubj);
end

colRaw = [0.62 0.62 0.62];
colLine = [0.16 0.32 0.62];
colArea = [0.84 0.28 0.16];
f_area = nan(1, nSubj);
f_max = nan(1, nSubj);
peaks = cell(1, nSubj);

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
tl = tiledlayout(2, 5, 'Padding', 'compact', 'TileSpacing', 'compact');

for s = 1:min(nSubj, 10)
    ax = nexttile(tl);
    hold(ax, 'on');
    y = ged.all_condition_powspctrm{cond, idx(s)};
    pk = peak_half_area(x, y);
    peaks{s} = pk;
    f_area(s) = pk.f_area;
    f_max(s) = pk.f_max;

    if isempty(pk.xx)
        title(ax, subjects{s}, 'FontSize', 13);
        axis(ax, 'off');
        continue
    end

    if isfinite(pk.f_area)
        [xL, topL, botL, xR, topR, botR] = split_peak_area(pk);
        shade_between(ax, xL, topL, botL, colArea, 0.50);
        shade_between(ax, xR, topR, botR, colArea, 0.16);
        plot(ax, pk.xseg, pk.chord, '--', 'Color', [0.15 0.15 0.15], ...
            'LineWidth', 1.1, 'HandleVisibility', 'off');
        plot(ax, pk.xx([pk.ilo, pk.ihi]), pk.ys([pk.ilo, pk.ihi]), 'o', ...
            'Color', [0.15 0.15 0.15], 'MarkerSize', 5, ...
            'MarkerFaceColor', 'w', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    plot(ax, pk.xx, pk.yy, 'Color', colRaw, 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    plot(ax, pk.xx, pk.ys, 'Color', colLine, 'LineWidth', 1.8, ...
        'HandleVisibility', 'off');

    yBand = [pk.yy, pk.ys];
    yBand = yBand(isfinite(yBand));
    pad = 0.08 * max(range(yBand), 0.5);
    ylim(ax, [min(yBand) - pad, max(yBand) + pad]);
    xlim(ax, [30 90]);
    xticks(ax, 30:20:90);

    if isfinite(pk.f_area)
        xline(ax, pk.f_area, 'Color', colArea, 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
        yl = ylim(ax);
        if pk.f_area > 72
            ha = 'right';
            hzLabel = sprintf('%.1f Hz  ', pk.f_area);
        else
            ha = 'left';
            hzLabel = sprintf('  %.1f Hz', pk.f_area);
        end
        text(ax, pk.f_area, yl(2) - 0.02 * range(yl), hzLabel, ...
            'Color', colArea, 'FontSize', 10, 'FontWeight', 'bold', ...
            'HorizontalAlignment', ha, 'VerticalAlignment', 'top', ...
            'Clipping', 'on');
    end

    if s > 5
        xlabel(ax, 'Frequency [Hz]', 'FontSize', 12);
    end
    if mod(s - 1, 5) == 0
        ylabel(ax, 'Power [dB]', 'FontSize', 12);
    end
    title(ax, subjects{s}, 'FontSize', 14, 'FontWeight', 'bold');
    box(ax, 'off');
    set(ax, 'FontSize', 11, 'FontName', 'Arial', 'TickDir', 'out');
end

hRaw = plot(nan, nan, 'Color', colRaw, 'LineWidth', 1);
hLine = plot(nan, nan, 'Color', colLine, 'LineWidth', 1.8);
hBound = plot(nan, nan, 'o', 'Color', [0.15 0.15 0.15], 'MarkerSize', 5, ...
    'MarkerFaceColor', 'w', 'LineWidth', 1);
hFirst = patch(nan, nan, colArea, 'FaceAlpha', 0.50, 'EdgeColor', 'none');
hSecond = patch(nan, nan, colArea, 'FaceAlpha', 0.16, 'EdgeColor', 'none');
hFreq = plot(nan, nan, 'Color', colArea, 'LineWidth', 1.5);
lgd = legend([hRaw, hLine, hBound, hFirst, hSecond, hFreq], ...
    {'Spectrum', 'Smoothed spectrum', 'Flanking minima', ...
    'First 50% of peak area', 'Remaining peak area', '50% area frequency'}, ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11, 'FontName', 'Arial');
lgd.Layout.Tile = 'south';
title(tl, '100% contrast, 50% area frequency of the gamma peak', ...
    'FontSize', 16, 'FontWeight', 'normal', 'FontName', 'Arial');

print(fig, fullfile(out_dir, 'GCP_control_gamma_peak_area_spectra.png'), '-dpng', '-r300');
close(fig);
save(fullfile(out_dir, 'GCP_control_gamma_peak_area_spectra.mat'), ...
    'subjects', 'cond', 'f_max', 'f_area', 'x', '-v7.3');
fprintf('[AREA SPEC] 100%% contrast. Median 50%% area frequency %.1f Hz. Median tallest sample %.1f Hz.\n', ...
    median(f_area, 'omitnan'), median(f_max, 'omitnan'));
fprintf('[AREA SPEC] Saved %s\n', fullfile(out_dir, 'GCP_control_gamma_peak_area_spectra.png'));

function pk = peak_half_area(x, y)
pk = struct('f_max', NaN, 'f_area', NaN, ...
    'xx', [], 'yy', [], 'ys', [], ...
    'ilo', NaN, 'ihi', NaN, ...
    'xseg', [], 'yseg', [], 'chord', [], 'above', []);
if isempty(y)
    return
end
x = x(:)';
y = y(:)';
if numel(y) ~= numel(x)
    return
end
mask = x >= 30 & x <= 90 & isfinite(x) & isfinite(y);
if nnz(mask) < 9
    return
end
xx = x(mask);
yy = y(mask);
ys = movmean(yy, 5, 'omitnan');
[~, ipeak] = max(ys);
[~, iraw] = max(yy);
pk.f_max = xx(iraw);
pk.xx = xx;
pk.yy = yy;
pk.ys = ys;

ilo = ipeak;
while ilo > 2 && ys(ilo) >= ys(ilo - 1)
    ilo = ilo - 1;
end
ihi = ipeak;
while ihi < numel(ys) - 1 && ys(ihi) >= ys(ihi + 1)
    ihi = ihi + 1;
end
pk.ilo = ilo;
pk.ihi = ihi;
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
cArea = cumtrapz(xseg, above);
k = find(cArea >= 0.5 * area, 1, 'first');
if isempty(k)
    return
end
if k == 1
    f50 = xseg(1);
else
    f50 = xseg(k - 1) + (xseg(k) - xseg(k - 1)) * ...
        ((0.5 * area - cArea(k - 1)) / max(cArea(k) - cArea(k - 1), eps));
end
pk.f_area = f50;
pk.xseg = xseg;
pk.yseg = yseg;
pk.chord = chord;
pk.above = above;
end

function [xL, topL, botL, xR, topR, botR] = split_peak_area(pk)
xseg = pk.xseg(:)';
top = pk.chord(:)' + pk.above(:)';
bot = pk.chord(:)';
f50 = pk.f_area;
[xL, topL, botL] = clip_to_frequency(xseg, top, bot, f50, 'left');
[xR, topR, botR] = clip_to_frequency(xseg, top, bot, f50, 'right');
end

function [xs, ts, bs] = clip_to_frequency(xseg, top, bot, f50, side)
if strcmp(side, 'left')
    keep = xseg < f50;
    xs = xseg(keep);
    ts = top(keep);
    bs = bot(keep);
    if isempty(xs) || xs(end) < f50
        xs = [xs, f50];
        ts = [ts, interp1(xseg, top, f50, 'linear')];
        bs = [bs, interp1(xseg, bot, f50, 'linear')];
    end
else
    keep = xseg > f50;
    xs = xseg(keep);
    ts = top(keep);
    bs = bot(keep);
    if isempty(xs) || xs(1) > f50
        xs = [f50, xs];
        ts = [interp1(xseg, top, f50, 'linear'), ts];
        bs = [interp1(xseg, bot, f50, 'linear'), bs];
    end
end
end

function shade_between(ax, xs, top, bot, color, alpha)
if numel(xs) < 2 || any(~isfinite([xs, top, bot]))
    return
end
fill(ax, [xs, fliplr(xs)], [top, fliplr(bot)], color, ...
    'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
