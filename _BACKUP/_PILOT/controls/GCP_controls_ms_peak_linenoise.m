%% Controls: microsaccade line vs box, gamma peak definitions, 40 to 80 Hz GED, line noise
% Results go to figures/controls. Does not overwrite confirmatory features.

startup
[subjects, paths, colors, headmodel] = setup('GCP');
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
    fprintf('[CTRL] Inclusion file was locked. Excluding 602, 604 and 608.\n');
end
addpath('/Users/Arne/Documents/GitHub/GCP/_BACKUP/_PILOT/2_feature_extraction/GED-helpers');
global ft_default
ft_default.showcallinfo = 'no';

out_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/_PILOT/figures/controls';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
set(0, 'DefaultFigureVisible', 'off');

nSubj = numel(subjects);
condFieldsHz = {'msTS_c25', 'msTS_c50', 'msTS_c75', 'msTS_c100'};
condFieldsBl = {'msTS_c25_bl', 'msTS_c50_bl', 'msTS_c75_bl', 'msTS_c100_bl'};
evFields = {'ms_events_c25', 'ms_events_c50', 'ms_events_c75', 'ms_events_c100'};
condLabels = {'25%', '50%', '75%', '100%'};
nCond = 4;

%% Microsaccades: boxplot scalar vs line plot
fprintf('\n[MS] Comparing boxplot scalars with line plot time courses\n');
gazeSum = load(fullfile(paths.features, 'GCP_gaze_window_summaries.mat'));
gaze_idx = match_subjects_local(gazeSum.subjects, subjects);
box = gazeSum.MSRate_bl.full(:, gaze_idx); % condition x subject, percent

fs_et = 500;
sigma_ms_disp = 20;
sigma_samp = max(1, round(sigma_ms_disp / (1000 / fs_et)));
kHalf = 3 * sigma_samp;
x_kern = -kHalf:kHalf;
gKernel = exp(-x_kern.^2 / (2 * sigma_samp^2));
gKernel = gKernel ./ sum(gKernel);

line_stim = nan(nCond, nSubj);
line_pre = nan(nCond, nSubj);
line_base = nan(nCond, nSubj);
hz_base = nan(nCond, nSubj);
hz_pre = nan(nCond, nSubj);
hz_stim = nan(nCond, nSubj);
pct_from_hz = nan(nCond, nSubj);
evt_all = nan(nCond, nSubj);
evt_kept = nan(nCond, nSubj);
n_trials = nan(nCond, nSubj);
n_kept = nan(nCond, nSubj);
t_ref = [];
curves_bl = nan(nSubj, 0, nCond);

for s = 1:nSubj
    spath = fullfile(paths.features, subjects{s}, 'gaze');
    dat = load(fullfile(spath, 'gaze_microsaccade_timeseries.mat'));
    ev = load(fullfile(spath, 'gaze_microsaccade_events.mat'));
    for c = 1:nCond
        msB = dat.(condFieldsBl{c});
        msH = dat.(condFieldsHz{c});
        t = msB.time(:)';
        yb = msB.avg(1, :);
        yh = msH.avg(1, :);
        yb_s = conv(yb, gKernel, 'same');
        if isempty(t_ref)
            t_ref = t;
            curves_bl = nan(nSubj, numel(t_ref), nCond);
        end
        curves_bl(s, :, c) = interp1(t, yb_s, t_ref, 'linear', NaN);
        line_stim(c, s) = mean(yb_s(t >= 0 & t <= 2), 'omitnan');
        line_pre(c, s) = mean(yb_s(t >= -0.5 & t < 0), 'omitnan');
        line_base(c, s) = mean(yb_s(t >= -1.5 & t <= -0.5), 'omitnan');
        hz_base(c, s) = mean(yh(t >= -1.5 & t <= -0.5), 'omitnan');
        hz_pre(c, s) = mean(yh(t >= -0.5 & t < 0), 'omitnan');
        hz_stim(c, s) = mean(yh(t >= 0 & t <= 2), 'omitnan');
        pct_from_hz(c, s) = 100 * (hz_stim(c, s) - hz_base(c, s)) ./ hz_base(c, s);

        onsets = ev.(evFields{c}).Onset;
        pct_trials = nan(1, numel(onsets));
        keep = false(1, numel(onsets));
        for tr = 1:numel(onsets)
            ton = onsets{tr};
            if isempty(ton)
                ton = [];
            end
            ton = ton(:)';
            rb = sum(ton >= -1.5 & ton <= -0.5) / 1.0;
            rs = sum(ton >= 0 & ton <= 2) / 2.0;
            if rb > 0
                pct_trials(tr) = 100 * (rs - rb) / rb;
            end
            keep(tr) = rb >= 0.1 && rs >= 0.1;
        end
        n_trials(c, s) = numel(onsets);
        n_kept(c, s) = sum(keep);
        evt_all(c, s) = mean(pct_trials, 'omitnan');
        evt_kept(c, s) = mean(pct_trials(keep), 'omitnan');
    end
end

fprintf('[MS] Boxplot scalar: median %.1f%%, fraction negative %.2f\n', ...
    median(box(:), 'omitnan'), mean(box(:) < 0, 'omitnan'));
fprintf('[MS] Line mean 0 to 2 s: median %.1f%%, fraction positive %.2f\n', ...
    median(line_stim(:), 'omitnan'), mean(line_stim(:) > 0, 'omitnan'));
fprintf('[MS] Line mean -0.5 to 0 s: median %.1f%%, fraction positive %.2f\n', ...
    median(line_pre(:), 'omitnan'), mean(line_pre(:) > 0, 'omitnan'));
fprintf('[MS] Line mean -1.5 to -0.5 s (should be near 0 if the denominator matched): median %.1f%%\n', ...
    median(line_base(:), 'omitnan'));
fprintf('[MS] Percent from subject mean Hz, same windows as the boxplot: median %.1f%%, fraction negative %.2f\n', ...
    median(pct_from_hz(:), 'omitnan'), mean(pct_from_hz(:) < 0, 'omitnan'));
fprintf('[MS] Full trace event counts, all trials: median %.1f%%, fraction negative %.2f\n', ...
    median(evt_all(:), 'omitnan'), mean(evt_all(:) < 0, 'omitnan'));
fprintf('[MS] Full trace event counts, trials with both windows >= 0.1 Hz: median %.1f%%, fraction negative %.2f\n', ...
    median(evt_kept(:), 'omitnan'), mean(evt_kept(:) < 0, 'omitnan'));
fprintf('[MS] Median fraction of trials kept by the 0.1 Hz rule: %.2f\n', ...
    median(n_kept ./ n_trials, 'all', 'omitnan'));
fprintf('[MS] Mean Hz baseline / pre / stim: %.3f / %.3f / %.3f\n', ...
    mean(hz_base(:), 'omitnan'), mean(hz_pre(:), 'omitnan'), mean(hz_stim(:), 'omitnan'));

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
subplot(2, 2, 1); hold on;
for c = 1:nCond
    scatter(box(c, :), line_stim(c, :), 80, colors(c, :), 'filled', ...
        'MarkerFaceAlpha', 0.85);
end
lims = [min([box(:); line_stim(:)], [], 'omitnan'), max([box(:); line_stim(:)], [], 'omitnan')];
plot(lims, lims, 'k--');
xline(0, 'Color', [0.6 0.6 0.6]); yline(0, 'Color', [0.6 0.6 0.6]);
xlabel('Boxplot scalar [%]');
ylabel('Mean of line, 0 to 2 s [%]');
title('Same subjects, different summaries');
axis square; box off; set(gca, 'FontSize', 14);

subplot(2, 2, 2); hold on;
grand = squeeze(mean(curves_bl, 1, 'omitnan'));
n_ok = squeeze(sum(isfinite(curves_bl), 1));
sem = squeeze(std(curves_bl, 0, 1, 'omitnan')) ./ sqrt(max(n_ok, 1));
for c = 1:nCond
    plot(t_ref, grand(:, c), 'Color', colors(c, :), 'LineWidth', 2.5);
    plot(t_ref, grand(:, c) + sem(:, c), 'Color', colors(c, :), 'LineWidth', 0.5);
    plot(t_ref, grand(:, c) - sem(:, c), 'Color', colors(c, :), 'LineWidth', 0.5);
end
xline(0, '--', 'Color', [0.5 0.5 0.5]);
yline(0, '--', 'Color', [0.5 0.5 0.5]);
xline(-0.5, ':', 'Color', [0.3 0.3 0.3]);
xlim([-1.5 2]);
xlabel('Time [s]');
ylabel('Microsaccade rate [%]');
title('Line plot. Baseline used for % is -1.5 to -0.5 s');
box off; set(gca, 'FontSize', 14);
legend(condLabels, 'Location', 'northeast', 'Box', 'off');

subplot(2, 2, 3); hold on;
vals = {box(:), evt_all(:), evt_kept(:), line_stim(:), pct_from_hz(:)};
labs = {'Box all trials', 'Events all', 'Events >= 0.1 Hz', 'Line 0 to 2 s', 'Hz windows'};
for i = 1:numel(vals)
    v = vals{i};
    v = v(isfinite(v));
    boxplot(v, ones(numel(v), 1) * i, 'Positions', i, 'Symbol', '', 'Widths', 0.55, 'Colors', 'k');
end
yline(0, '--', 'Color', [0.5 0.5 0.5]);
xlim([0.4 5.6]);
set(gca, 'XTick', 1:5, 'XTickLabel', labs, 'XTickLabelRotation', 20, 'FontSize', 12);
ylabel('Percent change');
title('What moves the sign');
box off;

subplot(2, 2, 4); hold on;
win_mu = [mean(hz_base, 2, 'omitnan'), mean(hz_pre, 2, 'omitnan'), mean(hz_stim, 2, 'omitnan')];
b = bar(win_mu');
for c = 1:nCond
    b(c).FaceColor = colors(c, :);
end
set(gca, 'XTick', 1:3, 'XTickLabel', {'-1.5 to -0.5 s', '-0.5 to 0 s', '0 to 2 s'}, 'FontSize', 14);
ylabel('Rate [Hz]');
title('Unbaselined line, subject mean');
legend(condLabels, 'Location', 'northeast', 'Box', 'off');
box off;
print(fig, fullfile(out_dir, 'GCP_control_ms_line_vs_box.png'), '-dpng', '-r300');
close(fig);

%% Gamma peak definitions on the saved 30 to 90 Hz GED spectra
fprintf('\n[PEAK] Refitting saved GED condition spectra\n');
ged = load(fullfile(paths.features, 'GCP_eeg_GED.mat'), ...
    'all_condition_powspctrm', 'all_condition_peak_freq', 'all_condition_peak_power', ...
    'scan_freqs', 'subjects');
ged_idx = match_subjects_local(ged.subjects, subjects);
scan_freqs = ged.scan_freqs(:)';
nF = numel(scan_freqs);
spec = nan(nCond, nSubj, nF);
peak_saved = ged.all_condition_peak_freq(:, ged_idx);
for s = 1:nSubj
    for c = 1:nCond
        y = ged.all_condition_powspctrm{c, ged_idx(s)};
        if ~isempty(y) && numel(y) == nF
            spec(c, s, :) = y(:)';
        end
    end
end

methods = {'Tallest', 'FOOOF', 'Laplace', 'HalfArea'};
bands = {[30 90], [40 80]};
band_names = {'30 to 90', '40 to 80'};
peak_hz = nan(numel(bands), numel(methods), nCond, nSubj);
fit_r2 = nan(nCond, nSubj);
for s = 1:nSubj
    for c = 1:nCond
        y = squeeze(spec(c, s, :))';
        if all(~isfinite(y))
            continue
        end
        for b = 1:numel(bands)
            [hz, r2] = peak_definitions(y, scan_freqs, bands{b});
            peak_hz(b, :, c, s) = hz;
            if b == 1
                fit_r2(c, s) = r2;
            end
        end
    end
    fprintf('[PEAK] %s done\n', subjects{s});
end

for b = 1:numel(bands)
    fprintf('[PEAK] Search band %s\n', band_names{b});
    for m = 1:numel(methods)
        v = squeeze(peak_hz(b, m, :, :));
        near50 = mean(abs(v(:) - 50) <= 2, 'omitnan');
        fprintf('  %s: median %.1f Hz, fraction within 2 Hz of 50 Hz %.2f, n finite %d\n', ...
            methods{m}, median(v(:), 'omitnan'), near50, sum(isfinite(v(:))));
    end
end
fprintf('[PEAK] Saved pick_tallest median %.1f Hz. Median |FOOOF minus tallest| %.2f Hz. Median |Laplace minus tallest| %.2f Hz. Median |half area minus tallest| %.2f Hz.\n', ...
    median(peak_saved(:), 'omitnan'), ...
    median(abs(squeeze(peak_hz(1, 2, :, :)) - peak_saved), 'all', 'omitnan'), ...
    median(abs(squeeze(peak_hz(1, 3, :, :)) - peak_saved), 'all', 'omitnan'), ...
    median(abs(squeeze(peak_hz(1, 4, :, :)) - peak_saved), 'all', 'omitnan'));
fprintf('[PEAK] Median FOOOF r2 on the stim/base ratio %.2f\n', median(fit_r2(:), 'omitnan'));

% 50 percent area latency of the sustained GED gamma time course (dB)
lat50 = nan(nCond, nSubj);
tfr_file = fullfile(paths.features, 'GCP_eeg_GED_TFR.mat');
if isfile(tfr_file)
    tfr = load(tfr_file, 'tfr_cond_avg', 'subjects');
    tfr_idx = match_subjects_local(tfr.subjects, subjects);
    for s = 1:nSubj
        for c = 1:nCond
            lat50(c, s) = half_area_latency(tfr.tfr_cond_avg{c, tfr_idx(s)});
        end
    end
    fprintf('[PEAK] 50%% area latency of GED gamma (0 to 2 s, positive dB only): median %.3f s\n', ...
        median(lat50(:), 'omitnan'));
end

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
subplot(2, 2, 1); hold on;
ga = squeeze(mean(spec, 2, 'omitnan'));
for c = 1:nCond
    plot(scan_freqs, ga(c, :), 'Color', colors(c, :), 'LineWidth', 2.2);
end
xline(50, '--', 'Color', [0.4 0.4 0.4]);
xlabel('Frequency [Hz]'); ylabel('Stim / base [dB]');
title('Saved GED spectra');
legend(condLabels, 'Box', 'off', 'Location', 'northeast');
box off; set(gca, 'FontSize', 14); xlim([30 90]);

subplot(2, 2, 2); hold on;
for m = 1:numel(methods)
    v = squeeze(peak_hz(1, m, :, :));
    for c = 1:nCond
        jitter = (rand(1, nSubj) - 0.5) * 0.15;
        scatter(m + (c - 2.5) * 0.12 + jitter, v(c, :), 36, colors(c, :), 'filled', ...
            'MarkerFaceAlpha', 0.8);
    end
    plot(m, median(v(:), 'omitnan'), 'k_', 'MarkerSize', 28, 'LineWidth', 2);
end
yline(50, '--', 'Color', [0.5 0.5 0.5]);
xlim([0.5 4.5]);
set(gca, 'XTick', 1:4, 'XTickLabel', methods, 'FontSize', 14);
ylabel('Peak frequency [Hz]');
title('Definitions on 30 to 90 Hz');
box off;

subplot(2, 2, 3); hold on;
for m = 1:numel(methods)
    v = squeeze(peak_hz(2, m, :, :));
    for c = 1:nCond
        jitter = (rand(1, nSubj) - 0.5) * 0.15;
        scatter(m + (c - 2.5) * 0.12 + jitter, v(c, :), 36, colors(c, :), 'filled', ...
            'MarkerFaceAlpha', 0.8);
    end
    plot(m, median(v(:), 'omitnan'), 'k_', 'MarkerSize', 28, 'LineWidth', 2);
end
yline(50, '--', 'Color', [0.5 0.5 0.5]);
xlim([0.5 4.5]); ylim([35 85]);
set(gca, 'XTick', 1:4, 'XTickLabel', methods, 'FontSize', 14);
ylabel('Peak frequency [Hz]');
title('Same spectra, search restricted to 40 to 80 Hz');
box off;

subplot(2, 2, 4); hold on;
if any(isfinite(lat50(:)))
    for c = 1:nCond
        scatter(c * ones(1, nSubj), lat50(c, :), 70, colors(c, :), 'filled', ...
            'MarkerFaceAlpha', 0.85);
    end
    xlim([0.5 4.5]);
    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 14);
    ylabel('Latency [s]');
    title('50% area latency is a time, not a frequency');
    box off;
else
    axis off;
    text(0.1, 0.5, 'TFR file missing', 'FontSize', 16);
end
print(fig, fullfile(out_dir, 'GCP_control_gamma_peak_methods.png'), '-dpng', '-r300');
close(fig);

%% Line noise and a real 40 to 80 Hz GED refit
fprintf('\n[EEG] Line noise and 40 to 80 Hz GED. One subject at a time.\n');
lambda = 0.05;
baseline_window = [-1.5, -0.5];
stim_window = [0, 2.0];
scan_score = 30:1:90;
bands_ged = {[30 90], [40 80]};
ged_names = {'30 to 90', '40 to 80'};

line_excess = nan(nSubj, 2, 2); % subj, region (occ, front), window (stim, base)
line_excess_100 = nan(nSubj, 2);
spec_occ_stim = [];
spec_occ_base = [];
freq_axis = [];
ged_out = cell(nSubj, 2);

for s = 1:nSubj
    t_subj = tic;
    eeg_path = fullfile(paths.features, subjects{s}, 'eeg', 'dataEEG.mat');
    fprintf('[EEG] Loading %s\n', subjects{s});
    try
        eeg = load(eeg_path, 'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
        dat_per_cond = {eeg.dataEEG_c25, eeg.dataEEG_c50, eeg.dataEEG_c75, eeg.dataEEG_c100};
        clear eeg
        fs = dat_per_cond{1}.fsample;
        labels = dat_per_cond{1}.label;
        occ = labels(occ_mask_of(labels));
        front = labels(front_mask_of(labels));
        try
            [f_axis, p_stim_occ, p_base_occ, ex_stim_occ, ex_base_occ] = ...
                region_line_spectra(dat_per_cond, occ, stim_window, baseline_window);
            [~, ~, ~, ex_stim_front, ex_base_front] = ...
                region_line_spectra(dat_per_cond, front, stim_window, baseline_window);
            if isempty(freq_axis)
                freq_axis = f_axis;
                spec_occ_stim = nan(nSubj, numel(freq_axis));
                spec_occ_base = nan(nSubj, numel(freq_axis));
            end
            spec_occ_stim(s, :) = p_stim_occ;
            spec_occ_base(s, :) = p_base_occ;
            line_excess(s, 1, :) = [ex_stim_occ, ex_base_occ];
            line_excess(s, 2, :) = [ex_stim_front, ex_base_front];
            line_excess_100(s, 1) = harmonic_excess_db(f_axis, p_stim_occ, 100);
            line_excess_100(s, 2) = harmonic_excess_db(f_axis, p_base_occ, 100);
        catch ME
            fprintf('[EEG] %s line spectra failed: %s\n', subjects{s}, ME.message);
            ex_stim_occ = NaN;
            ex_base_occ = NaN;
        end

        for b = 1:2
            ged_out{s, b} = summarize_ged_band( ...
                dat_per_cond, bands_ged{b}, baseline_window, stim_window, ...
                lambda, headmodel, scan_score);
        end
        clear dat_per_cond
        fprintf('[EEG] %s done in %.1f min. 50 Hz excess stim occ %.2f dB, base %.2f dB. Eligible %d vs %d.\n', ...
            subjects{s}, toc(t_subj) / 60, ex_stim_occ, ex_base_occ, ...
            ged_out{s, 1}.n_eligible, ged_out{s, 2}.n_eligible);
    catch ME
        fprintf('[EEG] %s failed: %s\n', subjects{s}, ME.message);
    end
end

fprintf('[LINE] Occipital 50 Hz excess, stim median %.2f dB, baseline median %.2f dB\n', ...
    median(line_excess(:, 1, 1), 'omitnan'), median(line_excess(:, 1, 2), 'omitnan'));
fprintf('[LINE] Frontal 50 Hz excess, stim median %.2f dB, baseline median %.2f dB\n', ...
    median(line_excess(:, 2, 1), 'omitnan'), median(line_excess(:, 2, 2), 'omitnan'));
fprintf('[LINE] Occipital 100 Hz excess, stim median %.2f dB, baseline median %.2f dB\n', ...
    median(line_excess_100(:, 1), 'omitnan'), median(line_excess_100(:, 2), 'omitnan'));

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
subplot(2, 2, 1); hold on;
if ~isempty(freq_axis)
    mu_s = mean(10 * log10(spec_occ_stim), 1, 'omitnan');
    mu_b = mean(10 * log10(spec_occ_base), 1, 'omitnan');
    plot(freq_axis, mu_b, 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
    plot(freq_axis, mu_s, 'Color', [0.1 0.35 0.7], 'LineWidth', 2);
    xline(50, '--', 'Color', [0.8 0.2 0.2]);
    xline(100, '--', 'Color', [0.8 0.2 0.2]);
    xlim([1 120]);
    xlabel('Frequency [Hz]'); ylabel('Occipital power [dB]');
    legend({'Baseline', 'Stimulus'}, 'Box', 'off');
    title('Hanning spectra, no notch in preprocessing');
    box off; set(gca, 'FontSize', 14);
end

subplot(2, 2, 2); hold on;
for s = 1:nSubj
    plot([1 2], squeeze(line_excess(s, 1, :)), '-o', 'Color', [0.7 0.7 0.7]);
end
plot([1 2], squeeze(median(line_excess(:, 1, :), 1, 'omitnan')), '-o', ...
    'Color', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'k');
xlim([0.7 2.3]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Stimulus', 'Baseline'}, 'FontSize', 14);
ylabel('50 Hz excess [dB]');
title('Occipital 50 Hz above 45 to 48 and 52 to 55 Hz');
box off;

subplot(2, 2, 3); hold on;
for s = 1:nSubj
    plot([1 2], squeeze(line_excess(s, 2, :)), '-o', 'Color', [0.7 0.7 0.7]);
end
plot([1 2], squeeze(median(line_excess(:, 2, :), 1, 'omitnan')), '-o', ...
    'Color', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'k');
xlim([0.7 2.3]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Stimulus', 'Baseline'}, 'FontSize', 14);
ylabel('50 Hz excess [dB]');
title('Frontal 50 Hz excess');
box off;

subplot(2, 2, 4); hold on;
scatter(line_excess(:, 1, 1), line_excess_100(:, 1), 80, [0.1 0.35 0.7], 'filled');
xline(0, '--', 'Color', [0.5 0.5 0.5]); yline(0, '--', 'Color', [0.5 0.5 0.5]);
xlabel('50 Hz excess, stimulus [dB]');
ylabel('100 Hz excess, stimulus [dB]');
title('Harmonics, occipital stimulus');
box off; set(gca, 'FontSize', 14);
print(fig, fullfile(out_dir, 'GCP_control_line_noise.png'), '-dpng', '-r300');
close(fig);

% GED band comparison figure
hf = nan(nSubj, 2);
ex50 = nan(nSubj, 2);
prom = nan(nSubj, 2);
lam = nan(nSubj, 2);
nelig = nan(nSubj, 2);
pk = nan(nSubj, 2);
spec_ged = nan(nSubj, 2, numel(scan_score));
for s = 1:nSubj
    for b = 1:2
        g = ged_out{s, b};
        if isempty(g), continue, end
        hf(s, b) = g.hf_slope;
        ex50(s, b) = g.excess50_db;
        prom(s, b) = g.peak_prom_db;
        lam(s, b) = g.lambda_top;
        nelig(s, b) = g.n_eligible;
        pk(s, b) = g.peak_hz;
        if ~isempty(g.spectrum) && numel(g.spectrum) == numel(scan_score)
            spec_ged(s, b, :) = g.spectrum(:)';
        end
    end
end
fprintf('[GED] Median top lambda 30 to 90: %.3f; 40 to 80: %.3f\n', ...
    median(lam(:, 1), 'omitnan'), median(lam(:, 2), 'omitnan'));
fprintf('[GED] Median eligible count 30 to 90: %.1f; 40 to 80: %.1f\n', ...
    median(nelig(:, 1), 'omitnan'), median(nelig(:, 2), 'omitnan'));
fprintf('[GED] Median HF slope 70 to 90, 30 to 90 filter: %.3f; 40 to 80 filter: %.3f\n', ...
    median(hf(:, 1), 'omitnan'), median(hf(:, 2), 'omitnan'));
fprintf('[GED] Median 50 Hz excess on ratio spectrum, 30 to 90: %.2f dB; 40 to 80: %.2f dB\n', ...
    median(ex50(:, 1), 'omitnan'), median(ex50(:, 2), 'omitnan'));
fprintf('[GED] Median peak prominence, 30 to 90: %.2f dB; 40 to 80: %.2f dB\n', ...
    median(prom(:, 1), 'omitnan'), median(prom(:, 2), 'omitnan'));
fprintf('[GED] Median peak Hz, 30 to 90: %.1f; 40 to 80: %.1f\n', ...
    median(pk(:, 1), 'omitnan'), median(pk(:, 2), 'omitnan'));

fig = figure('Position', [0 0 1512 982], 'Color', 'w');
subplot(2, 2, 1); hold on;
cols_b = [0.2 0.2 0.2; 0.85 0.35 0.1];
for b = 1:2
    mu = squeeze(mean(spec_ged(:, b, :), 1, 'omitnan'));
    plot(scan_score, mu, 'Color', cols_b(b, :), 'LineWidth', 2.4);
end
xline(50, '--', 'Color', [0.6 0.6 0.6]);
xlim([30 90]);
xlabel('Frequency [Hz]'); ylabel('Stim / base [dB]');
legend(ged_names, 'Box', 'off');
title('Top eligible component, or top component if none passed');
box off; set(gca, 'FontSize', 13);

subplot(2, 2, 2); hold on;
paired_dots(hf, cols_b);
set(gca, 'XTick', [1 2], 'XTickLabel', ged_names, 'FontSize', 13);
ylabel('HF slope, 70 to 90 Hz');
title('Positive slope is the EMG gate failure');
yline(0, '--', 'Color', [0.5 0.5 0.5]);
box off;

subplot(2, 2, 3); hold on;
paired_dots(ex50, cols_b);
set(gca, 'XTick', [1 2], 'XTickLabel', ged_names, 'FontSize', 13);
ylabel('50 Hz excess on ratio [dB]');
title('Line noise left in the GED ratio');
yline(0, '--', 'Color', [0.5 0.5 0.5]);
box off;

subplot(2, 2, 4); hold on;
paired_dots(prom, cols_b);
set(gca, 'XTick', [1 2], 'XTickLabel', ged_names, 'FontSize', 13);
ylabel('Peak minus median [dB]');
title('Peak prominence on 30 to 90 Hz');
box off;
print(fig, fullfile(out_dir, 'GCP_control_ged_40to80.png'), '-dpng', '-r300');
close(fig);

save(fullfile(out_dir, 'GCP_control_ms_peak_linenoise.mat'), ...
    'subjects', 'box', 'line_stim', 'line_pre', 'line_base', ...
    'hz_base', 'hz_pre', 'hz_stim', 'pct_from_hz', 'evt_all', 'evt_kept', ...
    'n_trials', 'n_kept', 't_ref', 'curves_bl', ...
    'peak_hz', 'peak_saved', 'fit_r2', 'lat50', 'methods', 'scan_freqs', ...
    'line_excess', 'line_excess_100', 'freq_axis', 'spec_occ_stim', 'spec_occ_base', ...
    'hf', 'ex50', 'prom', 'lam', 'nelig', 'pk', 'spec_ged', 'scan_score', '-v7.3');
fprintf('\n[DONE] Figures written to %s\n', out_dir);

%% Local functions
function idx = match_subjects_local(file_subjects, subjects)
idx = nan(1, numel(subjects));
for i = 1:numel(subjects)
    hit = find(strcmp(file_subjects, subjects{i}), 1, 'first');
    if ~isempty(hit)
        idx(i) = hit;
    end
end
if any(~isfinite(idx))
    error('Subject mismatch between inclusion list and saved features.');
end
end

function [hz, r2] = peak_definitions(y, x, band)
hz = nan(1, 4);
r2 = NaN;
mask = x >= band(1) & x <= band(2) & isfinite(y) & isfinite(x);
if nnz(mask) < 7
    return
end
xx = x(mask);
yy = y(mask);
[hz(1), ~] = pick_tallest_peak(yy, xx, 0, 0);
lin = 10 .^ (yy / 10);
[hz(2), r2] = fooof_center(xx, lin, band);
hz(3) = laplace_center(xx, yy);
hz(4) = half_area_frequency(xx, yy);
end

function [mu, r2] = fooof_center(x, lin_power, band)
mu = NaN;
r2 = NaN;
settings = struct();
settings.peak_width_limits = [2, 20];
settings.max_n_peaks = 3;
settings.min_peak_height = 0.05;
settings.peak_threshold = 1.5;
settings.aperiodic_mode = 'fixed';
settings.verbose = false;
try
    fr = fooof(x, lin_power, band, settings, true);
    r2 = fr.r_squared;
    if ~isempty(fr.peak_params)
        pp = double(fr.peak_params);
        if size(pp, 2) >= 2
            [~, ix] = max(pp(:, 2));
            mu = pp(ix, 1);
        elseif numel(pp) >= 1
            mu = pp(1);
        end
    end
catch
    mu = NaN;
end
if ~isfinite(mu)
    mu = gaussian_center(x, 10 * log10(lin_power));
    r2 = NaN;
end
end

function mu = gaussian_center(x, y_db)
mu = NaN;
y = y_db(:);
x = x(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
if numel(x) < 7
    return
end
[~, i0] = max(y);
p0 = [max(y) - median(y), x(i0), 8, median(y), 0];
lb = [0, min(x), 2, min(y) - 5, -2];
ub = [max(20, range(y) * 2), max(x), 25, max(y) + 5, 2];
obj = @(p) mean((y - (p(1) * exp(-0.5 * ((x - p(2)) / p(3)).^2) + p(4) + p(5) * log10(x))).^2, 'omitnan');
opts = optimset('Display', 'off', 'MaxIter', 800, 'MaxFunEvals', 2000);
try
    p = fminsearch(@(pp) bounded_obj(pp, lb, ub, obj), p0, opts);
    p = min(max(p, lb), ub);
    mu = p(2);
catch
    mu = x(i0);
end
end

function mu = laplace_center(x, y)
mu = NaN;
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
if numel(x) < 7
    return
end
[x, ix] = sort(x);
y = y(ix);
y = movmean(y, 5, 'omitnan');
floor_y = prctile(y, 20);
y_pos = max(y - floor_y, 0);
scale = max(y_pos);
if ~isfinite(scale) || scale <= eps
    return
end
y_norm = y_pos / scale;
[~, i0] = max(y_pos);
p0 = [1, x(i0), 6];
lb = [0.2, min(x), 1];
ub = [1.5, max(x), 20];
obj = @(p) mean((y_norm(:)' - laplace_template_model(x(:)', p)).^2, 'all', 'omitnan');
opts = optimset('Display', 'off', 'MaxIter', 800, 'MaxFunEvals', 2000);
try
    p = fminsearch(@(pp) bounded_obj(pp, lb, ub, obj), p0, opts);
    p = min(max(p, lb), ub);
    mu = p(2);
catch
    mu = x(i0);
end
end

function f50 = half_area_frequency(x, y)
f50 = NaN;
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
if numel(x) < 5
    return
end
floor_y = prctile(y, 20);
y_pos = max(y - floor_y, 0);
area = trapz(x, y_pos);
if ~isfinite(area) || area <= 0
    return
end
c = cumtrapz(x, y_pos);
k = find(c >= 0.5 * area, 1, 'first');
if isempty(k)
    return
end
if k == 1
    f50 = x(1);
else
    c0 = c(k - 1);
    f50 = x(k - 1) + (x(k) - x(k - 1)) * ((0.5 * area - c0) / max(c(k) - c0, eps));
end
end

function lat = half_area_latency(freq)
lat = NaN;
if ~isstruct(freq) || ~isfield(freq, 'powspctrm') || ~isfield(freq, 'time') || ~isfield(freq, 'freq')
    return
end
p = double(freq.powspctrm);
t = freq.time(:)';
f = freq.freq(:)';
if ndims(p) == 3
    p = squeeze(mean(p, 1, 'omitnan'));
end
if isvector(p)
    p = p(:)';
end
if size(p, 1) == numel(f) && size(p, 2) == numel(t)
    % freq x time already
elseif size(p, 2) == numel(f) && size(p, 1) == numel(t)
    p = p';
elseif size(p, 1) == numel(t) && size(p, 2) == numel(f)
    p = p';
else
    if size(p, 2) ~= numel(t)
        return
    end
end
fmask = f >= 40 & f <= 80;
tmask = t >= 0 & t <= 2;
if ~any(fmask) || nnz(tmask) < 5
    return
end
yt = mean(p(fmask, tmask), 1, 'omitnan');
tt = t(tmask);
yp = max(yt, 0);
area = trapz(tt, yp);
if ~isfinite(area) || area <= 0
    return
end
c = cumtrapz(tt, yp);
k = find(c >= 0.5 * area, 1, 'first');
if isempty(k)
    return
end
if k == 1
    lat = tt(1);
else
    lat = tt(k - 1) + (tt(k) - tt(k - 1)) * ((0.5 * area - c(k - 1)) / max(c(k) - c(k - 1), eps));
end
end

function m = occ_mask_of(labels)
m = cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO)', 'once')), labels);
end

function m = front_mask_of(labels)
m = cellfun(@(l) ~isempty(regexp(l, '^(Fp|AF|F)', 'once')), labels);
end

function [f_axis, p_stim, p_base, ex_stim, ex_base] = region_line_spectra(dat_per_cond, chans, stim_window, baseline_window)
f_axis = [];
p_stim = [];
p_base = [];
ex_stim = NaN;
ex_base = NaN;
if isempty(chans)
    return
end
cfg = [];
cfg.channel = chans;
parts = cell(1, numel(dat_per_cond));
for c = 1:numel(dat_per_cond)
    if isempty(dat_per_cond{c})
        continue
    end
    parts{c} = ft_selectdata(cfg, dat_per_cond{c});
end
parts = parts(~cellfun('isempty', parts));
dat = ft_appenddata([], parts{:});
fs = dat.fsample;
foi = 1:1:min(140, floor(fs / 2) - 1);
p_stim = hanning_mean_spectrum(dat, stim_window, foi);
p_base = hanning_mean_spectrum(dat, baseline_window, foi);
f_axis = foi;
ex_stim = harmonic_excess_db(foi, p_stim, 50);
ex_base = harmonic_excess_db(foi, p_base, 50);
end

function pow = hanning_mean_spectrum(dat, latency, foi)
cfg = [];
cfg.latency = latency;
d = ft_selectdata(cfg, dat);
cfgf = [];
cfgf.method = 'mtmfft';
cfgf.taper = 'hanning';
cfgf.output = 'pow';
cfgf.foi = foi;
cfgf.keeptrials = 'yes';
cfgf.pad = 4;
cfgf.feedback = 'none';
freq = ft_freqanalysis(cfgf, d);
pow = squeeze(mean(mean(double(freq.powspctrm), 1, 'omitnan'), 2, 'omitnan'));
pow = pow(:)';
end

function ex = harmonic_excess_db(freq, pow, f0)
ex = NaN;
if isempty(freq) || isempty(pow)
    return
end
freq = freq(:)';
pow = pow(:)';
i0 = find(abs(freq - f0) < 0.51, 1, 'first');
neigh = (freq >= f0 - 5 & freq <= f0 - 2) | (freq >= f0 + 2 & freq <= f0 + 5);
if isempty(i0) || nnz(neigh) < 4 || ~isfinite(pow(i0)) || pow(i0) <= 0
    return
end
bg = mean(pow(neigh), 'omitnan');
if ~isfinite(bg) || bg <= 0
    return
end
ex = 10 * log10(pow(i0) / bg);
end

function out = summarize_ged_band(dat_per_cond, band, baseline_window, stim_window, lambda, headmodel, scan_freqs)
out = empty_ged();
ref = [];
for c = 1:numel(dat_per_cond)
    if ~isempty(dat_per_cond{c})
        ref = dat_per_cond{c};
        break
    end
end
if isempty(ref)
    return
end
labels = ref.label;
nChans = numel(labels);
fs = ref.fsample;
cfg_filt = [];
cfg_filt.bpfilter = 'yes';
cfg_filt.bpfreq = band;
cfg_filt.bpfilttype = 'fir';
cfg_filt.bpfiltord = round(3 * fs / band(1));
cfg_filt.feedback = 'none';
covStim = zeros(nChans);
covBase = zeros(nChans);
nTot = 0;
for c = 1:numel(dat_per_cond)
    dat = dat_per_cond{c};
    if isempty(dat) || ~isfield(dat, 'trial') || isempty(dat.trial)
        continue
    end
    dat_f = ft_preprocessing(cfg_filt, dat);
    nTrl = numel(dat_f.trial);
    if nTrl < 1
        continue
    end
    covBase = covBase + compute_pooled_covariance_window(dat_f, baseline_window) * nTrl;
    covStim = covStim + compute_pooled_covariance_window(dat_f, stim_window) * nTrl;
    nTot = nTot + nTrl;
    clear dat_f
end
if nTot < 1
    return
end
covStim = covStim / nTot;
covBase = covBase / nTot;
covStim_reg = (1 - lambda) * covStim + lambda * mean(diag(covStim)) * eye(nChans);
covBase_reg = (1 - lambda) * covBase + lambda * mean(diag(covBase)) * eye(nChans);
[W, D] = eig(covStim_reg, covBase_reg);
[evals, ord] = sort(real(diag(D)), 'descend');
W = W(:, ord);
nSearch = min(10, size(W, 2));
[post_idx, front_idx, temp_idx] = region_indices(labels);
sim_template = occipital_template(labels, headmodel, post_idx, front_idx);
filters = nan(nChans, nSearch);
topos = nan(nChans, nSearch);
post_front = nan(nSearch, 1);
post_temp = nan(nSearch, 1);
topo_frac = nan(nSearch, 1);
for ci = 1:nSearch
    w = W(:, ci);
    topo = covStim_reg * w;
    r = corr(topo, sim_template, 'rows', 'complete');
    if ~isnan(r) && r < 0
        w = -w;
        topo = -topo;
    end
    filters(:, ci) = w;
    topos(:, ci) = topo;
    post_s = mean(abs(topo(post_idx)));
    front_s = mean(abs(topo(front_idx)));
    temp_s = mean(abs(topo(temp_idx)));
    post_front(ci) = post_s / max(front_s, eps);
    post_temp(ci) = post_s / max(temp_s, eps);
    topo_abs = abs(topo);
    topo_frac(ci) = max(topo_abs) / max(sum(topo_abs), eps);
end
proxies = estimate_components_artifact_proxies( ...
    filters, dat_per_cond, stim_window, baseline_window, ...
    fs, scan_freqs, 3, sprintf('[GED %g to %g]', band(1), band(2)));
spectra = nan(nSearch, numel(scan_freqs));
hf_slope = zeros(nSearch, 1);
for ci = 1:nSearch
    spectra(ci, :) = proxies(ci).mean_pr_spectrum(:)';
    if isfinite(proxies(ci).hf_slope)
        hf_slope(ci) = proxies(ci).hf_slope;
    end
end
[pf, ~] = compute_powspctrm_form_laplacian_score_from_spectra(spectra, scan_freqs, band);
elig = evals(1:nSearch) >= 1.1 & pf >= 0.80 & post_front > 1 & post_temp > 1 & ...
    topo_frac <= 0.25 & hf_slope <= 0 & isfinite(pf);
out.n_eligible = sum(elig);
out.lambda_top = evals(1);
if any(elig)
    idx = find(elig);
    wts = evals(idx);
    wts = wts(:)' / sum(wts);
    Wsel = filters(:, idx);
    w_comb = Wsel * wts(:);
    spec = wts * spectra(idx, :);
    out.used = 'eligible';
else
    w_comb = filters(:, 1);
    spec = spectra(1, :);
    out.used = 'top';
end
out.spectrum = spec;
out.hf_slope = hf_slope_of(scan_freqs, spec);
out.excess50_db = band_excess_db(scan_freqs, spec, 50);
[out.peak_hz, out.peak_prom_db] = prominence_of(scan_freqs, spec);
out.lambda_used = evals(1);
end

function out = empty_ged()
out = struct('n_eligible', NaN, 'lambda_top', NaN, 'lambda_used', NaN, ...
    'hf_slope', NaN, 'excess50_db', NaN, 'peak_hz', NaN, 'peak_prom_db', NaN, ...
    'spectrum', [], 'used', '');
end

function [post_idx, front_idx, temp_idx] = region_indices(labels)
post_idx = find(cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P)', 'once')), labels));
front_idx = find(cellfun(@(l) ~isempty(regexp(l, '^(Fp|AF|F)', 'once')), labels));
temp_idx = find(cellfun(@(l) ~isempty(regexp(l, '^(T|TP|FT)', 'once')), labels));
end

function tmpl = occipital_template(labels, headmodel, post_idx, front_idx)
nChans = numel(labels);
tmpl = zeros(nChans, 1);
use_lay = isstruct(headmodel) && isfield(headmodel, 'layANThead') && ...
    isfield(headmodel.layANThead, 'label') && isfield(headmodel.layANThead, 'pos');
if use_lay
    lay_labels = headmodel.layANThead.label;
    lay_pos = headmodel.layANThead.pos;
    pos = nan(nChans, 2);
    for ch = 1:nChans
        li = find(strcmp(lay_labels, labels{ch}), 1, 'first');
        if ~isempty(li)
            pos(ch, :) = lay_pos(li, :);
        end
    end
    has = ~any(isnan(pos), 2);
    occ_idx = find(cellfun(@(l) ~isempty(regexp(l, '^(O|I|PO|PPO|P10|P9)', 'once')), labels));
    if any(has(occ_idx)) && any(has(front_idx))
        occ_ctr = mean(pos(occ_idx(has(occ_idx)), :), 1);
        front_ctr = mean(pos(front_idx(has(front_idx)), :), 1);
        d_occ = sum((pos - occ_ctr).^2, 2);
        d_front = sum((pos - front_ctr).^2, 2);
        tmpl = exp(-d_occ / (2 * 0.12^2)) - 0.75 * exp(-d_front / (2 * 0.25^2));
        tmpl(~has) = 0;
    end
end
if all(tmpl == 0)
    tmpl(post_idx) = 1;
    tmpl(front_idx) = -0.75;
end
if std(tmpl) > 0
    tmpl = (tmpl - mean(tmpl)) / std(tmpl);
end
end

function s = hf_slope_of(freq, spec)
s = NaN;
mask = freq >= 70 & freq <= 90;
if nnz(mask) < 3
    return
end
x = log(freq(mask));
y = log(abs(spec(mask)) + eps);
ok = isfinite(x) & isfinite(y);
if nnz(ok) < 3
    return
end
p = polyfit(x(ok), y(ok), 1);
s = p(1);
end

function ex = band_excess_db(freq, spec_db, f0)
ex = NaN;
i0 = find(abs(freq - f0) < 0.51, 1, 'first');
neigh = (abs(freq - f0) >= 4 & abs(freq - f0) <= 8);
if isempty(i0) || nnz(neigh) < 4 || ~isfinite(spec_db(i0))
    return
end
ex = spec_db(i0) - mean(spec_db(neigh), 'omitnan');
end

function [pk, prom] = prominence_of(freq, spec)
pk = NaN;
prom = NaN;
mask = freq >= 30 & freq <= 90 & isfinite(spec);
if nnz(mask) < 5
    return
end
yy = spec(mask);
xx = freq(mask);
[pk, ~] = pick_tallest_peak(yy, xx, 0, 0);
if ~isfinite(pk)
    [~, ix] = max(yy);
    pk = xx(ix);
end
prom = max(yy) - median(yy, 'omitnan');
end

function paired_dots(M, cols)
hold on;
for s = 1:size(M, 1)
    plot(1:2, M(s, :), '-', 'Color', [0.75 0.75 0.75]);
end
for b = 1:2
    scatter(b * ones(size(M, 1), 1), M(:, b), 70, cols(b, :), 'filled', ...
        'MarkerFaceAlpha', 0.9);
end
xlim([0.6 2.4]);
end
