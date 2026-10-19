%% 50% area frequency of the gamma peak
% Analogous to 50% area latency: the frequency at which half the area of the
% peak has accumulated. Bounds are the flanking minima of a lightly smoothed
% spectrum. Area is above the straight line joining those minima, so a narrow
% spike adds little area and does not set the frequency by itself.

startup
[subjects, paths, colors, ~] = setup('GCP');
try
    subjects = gcp_subject_inclusion(subjects, paths);
catch
    subjects = subjects(~ismember(str2double(string(subjects)), [602, 604, 608]));
end

out_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/controls';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
set(0, 'DefaultFigureVisible', 'off');

ged = load(fullfile(paths.features, 'GCP_eeg_GED.mat'), ...
    'all_condition_powspctrm', 'all_condition_peak_freq', 'scan_freqs', 'subjects');
idx = nan(1, numel(subjects));
for i = 1:numel(subjects)
    idx(i) = find(strcmp(ged.subjects, subjects{i}), 1, 'first');
end
x = ged.scan_freqs(:)';
nCond = 4;
nSubj = numel(subjects);
condLabels = {'25%', '50%', '75%', '100%'};
f_max = nan(nCond, nSubj);
f_area = nan(nCond, nSubj);
spec = nan(nCond, nSubj, numel(x));

for s = 1:nSubj
    for c = 1:nCond
        y = ged.all_condition_powspctrm{c, idx(s)};
        if isempty(y) || numel(y) ~= numel(x)
            continue
        end
        y = y(:)';
        spec(c, s, :) = y;
        [f_max(c, s), f_area(c, s)] = peak_half_area(x, y);
    end
end

d = f_area - f_max;
near50 = abs(f_max - 50) <= 2;
fprintf('[AREA] Tallest peak median %.1f Hz\n', median(f_max(:), 'omitnan'));
fprintf('[AREA] 50%% area frequency median %.1f Hz\n', median(f_area(:), 'omitnan'));
fprintf('[AREA] Median |area minus tallest| %.2f Hz\n', median(abs(d(:)), 'omitnan'));
fprintf('[AREA] When tallest peak is within 2 Hz of 50: n=%d, median area frequency %.1f Hz, median shift %.2f Hz\n', ...
    sum(near50(:)), median(f_area(near50), 'omitnan'), median(d(near50), 'omitnan'));

saved = ged.all_condition_peak_freq(:, idx);
fprintf('[AREA] Saved pick_tallest median %.1f Hz. Median |area minus saved| %.2f Hz\n', ...
    median(saved(:), 'omitnan'), median(abs(f_area(:) - saved(:)), 'omitnan'));

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
examples = {'607', '603', '610', '606'};
for e = 1:numel(examples)
    s = find(strcmp(subjects, examples{e}), 1);
    subplot(2, 2, e); hold on;
    if isempty(s)
        title(sprintf('%s not included', examples{e}));
        continue
    end
    for c = 1:nCond
        plot(x, squeeze(spec(c, s, :)), 'Color', colors(c, :), 'LineWidth', 1.6);
        if isfinite(f_area(c, s))
            xline(f_area(c, s), 'Color', colors(c, :), 'LineStyle', '--', 'LineWidth', 1);
        end
    end
    xline(50, ':', 'Color', [0.5 0.5 0.5]);
    xlim([30 90]);
    xlabel('Frequency [Hz]');
    ylabel('Stim / base [dB]');
    title(sprintf('%s, dashed = 50%% area frequency', subjects{s}));
    box off;
    set(gca, 'FontSize', 12);
    if e == 1
        legend(condLabels, 'Box', 'off', 'Location', 'northeast');
    end
end
print(fig, fullfile(out_dir, 'GCP_control_gamma_peak_area.png'), '-dpng', '-r300');
close(fig);

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
for c = 1:nCond
    scatter(f_max(c, :), f_area(c, :), 80, colors(c, :), 'filled', 'MarkerFaceAlpha', 0.85);
end
plot([30 90], [30 90], 'k--');
xline(50, ':', 'Color', [0.6 0.6 0.6]);
yline(50, ':', 'Color', [0.6 0.6 0.6]);
xlim([30 90]); ylim([30 90]);
xlabel('Tallest peak [Hz]');
ylabel('50% area frequency [Hz]');
title('Peak frequency from half the area under the peak');
legend(condLabels, 'Box', 'off', 'Location', 'southeast');
axis square; box off; set(gca, 'FontSize', 14);
print(fig, fullfile(out_dir, 'GCP_control_gamma_peak_area_scatter.png'), '-dpng', '-r300');
close(fig);

save(fullfile(out_dir, 'GCP_control_gamma_peak_area.mat'), ...
    'subjects', 'f_max', 'f_area', 'spec', 'x', 'condLabels', '-v7.3');
fprintf('[AREA] Saved figures to %s\n', out_dir);

function [f_max, f_area] = peak_half_area(x, y)
f_max = NaN;
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
[~, iraw] = max(yy);
f_max = xx(iraw);

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
