%% GCP Hypothesis Testing — Full Analysis & Visualization
%
% Tests all hypotheses by loading pre-computed variables and generating
% results figures. Each figure addresses one hypothesis.
%
% Oculomotor:
%   H1: Microsaccade rate decreases with increasing contrast
%   H2: Eye velocity increases post-stimulus, suppressed by contrast
%   H3: Pupil constriction amplitude increases monotonically with contrast
%   H4: Gaze dispersion increases with stimulus contrast
%
% Gamma oscillations:
%   H5: Gamma peak frequency is higher at higher contrast
%   H6: Gamma peak amplitude peaks at 75% contrast (inverted U)
%
% Gamma as cortical index of visual gain:
%   H7: Gamma frequency relates to oculomotor dynamics;
%       Early (0–600 ms): MS suppression + reduced velocity with contrast
%       Late  (1–2 s):    greater gaze dispersion with contrast
%
% Dependencies:
%   - Per-subject gaze time series  (gaze_fex output)
%   - gaze_matrix.mat               (group-level scalar gaze)
%   - eeg_matrix.mat                 (group-level EEG peaks)
%   - GCP_eeg_GED_gamma_metrics.mat  (GED gamma metrics)
%   - FieldTrip on path

clear; close all; clc

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');
nSubj = length(subjects);

condLabels    = {'25%', '50%', '75%', '100%'};
contrast_vals = [25, 50, 75, 100];
fontSize      = 16;

if ispc
    data_dir = 'W:\Students\Arne\GCP\data\features';
    fig_dir  = 'W:\Students\Arne\GCP\figures\hypotheses';
else
    data_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features';
    fig_dir  = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/hypotheses';
end
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% ====================================================================
%  LOAD DATA
%  ====================================================================

% --- Scalar gaze measures (group-level) ---
load(fullfile(data_dir, 'gaze_matrix.mat'));   % → gaze_data
load(fullfile(data_dir, 'eeg_matrix.mat'));     % → eeg_data

% GED gamma metrics (load into struct to avoid overwriting 'subjects')
ged = load(fullfile(data_dir, 'GCP_eeg_GED_gamma_metrics.mat'));

% Organise scalar measures into [4 x nSubj] matrices
msRate_subj    = nan(4, nSubj);
vel2D_subj     = nan(4, nSubj);
pupilSize_subj = nan(4, nSubj);
gazeDev_subj   = nan(4, nSubj);
gazeStdX_subj  = nan(4, nSubj);
gazeStdY_subj  = nan(4, nSubj);
pctPupil_subj  = nan(4, nSubj);
pctGazeDev_subj = nan(4, nSubj);

for row = 1:numel(gaze_data)
    sid  = gaze_data(row).ID;
    cond = gaze_data(row).Condition;
    si   = find(strcmp(subjects, num2str(sid)));
    if isempty(si), continue; end

    msRate_subj(cond, si)    = gaze_data(row).MSRate;
    vel2D_subj(cond, si)     = gaze_data(row).Vel2D;
    pupilSize_subj(cond, si) = gaze_data(row).PupilSize;
    gazeDev_subj(cond, si)   = gaze_data(row).GazeDeviation;
    gazeStdX_subj(cond, si)  = gaze_data(row).GazeStdX;
    gazeStdY_subj(cond, si)  = gaze_data(row).GazeStdY;

    if isfield(gaze_data, 'PctPupilSize')
        pctPupil_subj(cond, si) = gaze_data(row).PctPupilSize;
    end
    if isfield(gaze_data, 'PctGazeDeviation')
        pctGazeDev_subj(cond, si) = gaze_data(row).PctGazeDeviation;
    end
end

% EEG peak frequency (FOOOF pipeline)
peakFreq_eeg = nan(4, nSubj);
for row = 1:numel(eeg_data)
    sid  = eeg_data(row).ID;
    cond = eeg_data(row).Condition;
    si   = find(strcmp(subjects, num2str(sid)));
    if isempty(si), continue; end
    peakFreq_eeg(cond, si) = eeg_data(row).Frequency;
end

% GED-based metrics [4 x nSubj]
ged_peakAmp  = ged.metric_peak_amp;
ged_bbPower  = ged.metric_bb_power;
ged_centroid = ged.metric_centroid;
ged_auc      = ged.metric_auc;

ged_peakFreq = nan(4, nSubj);
for cond = 1:4
    for s = 1:nSubj
        pf = ged.trl_peak_freq{cond, s};
        if ~isempty(pf)
            ged_peakFreq(cond, s) = nanmean(pf);
        end
    end
end

% --- Per-subject time series ---
fprintf('Loading per-subject time series...\n');

t_ms  = [];   % time axis for MS rate
t_vel = [];   % time axis for velocity
t_pup = [];   % time axis for pupil

ms_mat  = [];  % [nSubj x nTime] per condition
vel_mat = [];
pup_mat = [];

for subj = 1:nSubj
    gazepath = fullfile(path, subjects{subj}, 'gaze');

    % Microsaccade rate (raw Hz)
    msfile = fullfile(gazepath, 'gaze_microsaccade_timeseries.mat');
    if exist(msfile, 'file')
        S = load(msfile, 'msTS_c25', 'msTS_c50', 'msTS_c75', 'msTS_c100');
        ms_structs = {S.msTS_c25, S.msTS_c50, S.msTS_c75, S.msTS_c100};
        if isempty(t_ms)
            t_ms   = ms_structs{1}.time;
            ms_mat = nan(nSubj, length(t_ms), 4);
        end
        for cond = 1:4
            ch = find(strcmp(ms_structs{cond}.label, 'MSRate'));
            if ~isempty(ch)
                ms_mat(subj, :, cond) = ms_structs{cond}.avg(ch, :);
            end
        end
    end

    % Velocity (raw px/s)
    velfile = fullfile(gazepath, 'gaze_velocity_timeseries.mat');
    if exist(velfile, 'file')
        S = load(velfile, 'velTS_c25', 'velTS_c50', 'velTS_c75', 'velTS_c100');
        vel_structs = {S.velTS_c25, S.velTS_c50, S.velTS_c75, S.velTS_c100};
        if isempty(t_vel)
            t_vel   = vel_structs{1}.time;
            vel_mat = nan(nSubj, length(t_vel), 4);
        end
        ch = find(strcmp(vel_structs{1}.label, 'Vel2D'));
        for cond = 1:4
            vel_mat(subj, :, cond) = vel_structs{cond}.avg(ch, :);
        end
    end

    % Pupil (baselined, subtractive)
    pupfile = fullfile(gazepath, 'gaze_pupil_timeseries.mat');
    if exist(pupfile, 'file')
        S = load(pupfile, 'pupTS_c25_bl', 'pupTS_c50_bl', 'pupTS_c75_bl', 'pupTS_c100_bl');
        pup_structs = {S.pupTS_c25_bl, S.pupTS_c50_bl, S.pupTS_c75_bl, S.pupTS_c100_bl};
        if isempty(t_pup)
            t_pup   = pup_structs{1}.time;
            pup_mat = nan(nSubj, length(t_pup), 4);
        end
        for cond = 1:4
            pup_mat(subj, :, cond) = pup_structs{cond}.avg(1, :);
        end
    end

    fprintf('  Subject %s (%d/%d)\n', subjects{subj}, subj, nSubj);
end

% --- H7: Extract early/late window means from time series ---
early_win = [0, 0.6];
late_win  = [1, 2];

ms_early  = nan(4, nSubj);
ms_late   = nan(4, nSubj);
vel_early = nan(4, nSubj);
vel_late  = nan(4, nSubj);

if ~isempty(t_ms)
    e_idx_ms = t_ms >= early_win(1) & t_ms <= early_win(2);
    l_idx_ms = t_ms >= late_win(1)  & t_ms <= late_win(2);
    for cond = 1:4
        ms_early(cond, :) = nanmean(ms_mat(:, e_idx_ms, cond), 2)';
        ms_late(cond, :)  = nanmean(ms_mat(:, l_idx_ms, cond), 2)';
    end
end

if ~isempty(t_vel)
    e_idx_vel = t_vel >= early_win(1) & t_vel <= early_win(2);
    l_idx_vel = t_vel >= late_win(1)  & t_vel <= late_win(2);
    for cond = 1:4
        vel_early(cond, :) = nanmean(vel_mat(:, e_idx_vel, cond), 2)';
        vel_late(cond, :)  = nanmean(vel_mat(:, l_idx_vel, cond), 2)';
    end
end

fprintf('\nData loading complete. %d subjects.\n\n', nSubj);

%% ====================================================================
%  FIGURE 1: H1 — Microsaccade Rate Decreases with Contrast
%  ====================================================================
close all

fig1 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H1] Microsaccade Rate Decreases with Increasing Contrast', ...
    'FontSize', 18, 'FontWeight', 'bold');

% --- Panel a: Time course ---
subplot(1, 2, 1); hold on;
t_plot = [-0.5 2];
if ~isempty(t_ms)
    t_idx = t_ms >= t_plot(1) & t_ms <= t_plot(2);
    h = gobjects(4, 1);
    for c = 1:4
        dat = ms_mat(:, t_idx, c);
        n_valid = sum(~isnan(dat(:, 1)));
        mu  = nanmean(dat, 1);
        sem = nanstd(dat, [], 1) / sqrt(n_valid);

        fill([t_ms(t_idx), fliplr(t_ms(t_idx))], ...
            [mu + sem, fliplr(mu - sem)], ...
            colors(c,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        h(c) = plot(t_ms(t_idx), mu, '-', 'Color', colors(c,:), 'LineWidth', 2.5);
    end
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    legend(h, condLabels, 'Location', 'northeast', 'FontSize', fontSize - 4);
end
xlabel('Time [s]'); ylabel('Microsaccade Rate [Hz]');
title('Time Course'); set(gca, 'FontSize', fontSize - 2); xlim(t_plot);

% --- Panel b: Contrast response function ---
subplot(1, 2, 2); hold on;
plot_crf(msRate_subj, contrast_vals, colors, nSubj);
xlabel('Contrast [%]'); ylabel('Microsaccade Rate [Hz]');
title('Contrast Response'); set(gca, 'FontSize', fontSize - 2);

r_trend = report_trend('H1 MS Rate', msRate_subj, contrast_vals, nSubj);

saveas(fig1, fullfile(fig_dir, 'GCP_H1_microsaccade_rate.png'));

%% ====================================================================
%  FIGURE 2: H2 — Eye Velocity Post-Stimulus Dynamics
%  ====================================================================

fig2 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H2] Eye Velocity: Post-Stimulus Increase, Contrast-Dependent Suppression', ...
    'FontSize', 18, 'FontWeight', 'bold');

% --- Panel a: Time course ---
subplot(1, 2, 1); hold on;
if ~isempty(t_vel)
    t_idx = t_vel >= t_plot(1) & t_vel <= t_plot(2);
    h = gobjects(4, 1);
    for c = 1:4
        dat = vel_mat(:, t_idx, c);
        n_valid = sum(~isnan(dat(:, 1)));
        mu  = nanmean(dat, 1);
        sem = nanstd(dat, [], 1) / sqrt(n_valid);

        fill([t_vel(t_idx), fliplr(t_vel(t_idx))], ...
            [mu + sem, fliplr(mu - sem)], ...
            colors(c,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        h(c) = plot(t_vel(t_idx), mu, '-', 'Color', colors(c,:), 'LineWidth', 2.5);
    end
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    legend(h, condLabels, 'Location', 'northeast', 'FontSize', fontSize - 4);
end
xlabel('Time [s]'); ylabel('Combined Eye Velocity [px/s]');
title('Time Course'); set(gca, 'FontSize', fontSize - 2); xlim(t_plot);

% --- Panel b: CRF of mean velocity ---
subplot(1, 2, 2); hold on;
plot_crf(vel2D_subj, contrast_vals, colors, nSubj);
xlabel('Contrast [%]'); ylabel('Combined Eye Velocity [px/s]');
title('Contrast Response'); set(gca, 'FontSize', fontSize - 2);

report_trend('H2 Eye Velocity', vel2D_subj, contrast_vals, nSubj);

saveas(fig2, fullfile(fig_dir, 'GCP_H2_eye_velocity.png'));

%% ====================================================================
%  FIGURE 3: H3 — Pupil Constriction Scales with Contrast
%  ====================================================================

fig3 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H3] Pupil Constriction Amplitude Increases with Contrast', ...
    'FontSize', 18, 'FontWeight', 'bold');

% --- Panel a: Baselined time course ---
subplot(1, 2, 1); hold on;
if ~isempty(t_pup)
    t_idx = t_pup >= t_plot(1) & t_pup <= t_plot(2);
    h = gobjects(4, 1);
    for c = 1:4
        dat = pup_mat(:, t_idx, c);
        n_valid = sum(~isnan(dat(:, 1)));
        mu  = nanmean(dat, 1);
        sem = nanstd(dat, [], 1) / sqrt(n_valid);

        fill([t_pup(t_idx), fliplr(t_pup(t_idx))], ...
            [mu + sem, fliplr(mu - sem)], ...
            colors(c,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        h(c) = plot(t_pup(t_idx), mu, '-', 'Color', colors(c,:), 'LineWidth', 2.5);
    end
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    legend(h, condLabels, 'Location', 'southeast', 'FontSize', fontSize - 4);
end
xlabel('Time [s]'); ylabel('Pupil Size [baselined, a.u.]');
title('Time Course'); set(gca, 'FontSize', fontSize - 2); xlim(t_plot);

% --- Panel b: Constriction amplitude (min of baselined pupil) ---
subplot(1, 2, 2); hold on;
constriction = nan(4, nSubj);
if ~isempty(t_pup)
    cwin = t_pup >= 0.3 & t_pup <= 1.5;
    for cond = 1:4
        for s = 1:nSubj
            ts = pup_mat(s, cwin, cond);
            if any(~isnan(ts))
                constriction(cond, s) = min(ts);
            end
        end
    end
elseif ~all(isnan(pctPupil_subj(:)))
    constriction = pctPupil_subj;
end

plot_crf(constriction, contrast_vals, colors, nSubj);
xlabel('Contrast [%]'); ylabel('Max Constriction [a.u.]');
title('Constriction Amplitude'); set(gca, 'FontSize', fontSize - 2);

report_trend('H3 Pupil Constriction', constriction, contrast_vals, nSubj);

saveas(fig3, fullfile(fig_dir, 'GCP_H3_pupil_constriction.png'));

%% ====================================================================
%  FIGURE 4: H4 — Gaze Dispersion Increases with Contrast
%  ====================================================================

fig4 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H4] Gaze Dispersion Increases with Stimulus Contrast', ...
    'FontSize', 18, 'FontWeight', 'bold');

% --- Panel a: Gaze deviation CRF ---
subplot(1, 2, 1); hold on;
plot_crf(gazeDev_subj, contrast_vals, colors, nSubj);
xlabel('Contrast [%]'); ylabel('Gaze Deviation [px]');
title('Gaze Deviation'); set(gca, 'FontSize', fontSize - 2);

% --- Panel b: Gaze Std (X + Y combined) ---
subplot(1, 2, 2); hold on;
gazeStd_subj = sqrt(gazeStdX_subj.^2 + gazeStdY_subj.^2);
plot_crf(gazeStd_subj, contrast_vals, colors, nSubj);
xlabel('Contrast [%]'); ylabel('Gaze Std [px]');
title('Gaze Spatial Spread'); set(gca, 'FontSize', fontSize - 2);

report_trend('H4 Gaze Deviation', gazeDev_subj, contrast_vals, nSubj);
report_trend('H4 Gaze Std', gazeStd_subj, contrast_vals, nSubj);

saveas(fig4, fullfile(fig_dir, 'GCP_H4_gaze_dispersion.png'));

%% ====================================================================
%  FIGURE 5: H5 + H6 — Gamma Peak Frequency & Peak Amplitude
%  ====================================================================

fig5 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H5] Gamma Peak Frequency & [H6] Peak Amplitude vs Contrast', ...
    'FontSize', 18, 'FontWeight', 'bold');

% --- Panel a: H5 — Peak frequency increases with contrast ---
subplot(1, 2, 1); hold on;
plot_crf(ged_peakFreq, contrast_vals, colors, nSubj);
xlabel('Contrast [%]'); ylabel('GED Peak Frequency [Hz]');
title('[H5] Peak Frequency (linear increase)');
set(gca, 'FontSize', fontSize - 2);

% --- Panel b: H6 — Peak amplitude peaks at 75% (inverted U) ---
subplot(1, 2, 2); hold on;
plot_crf(ged_peakAmp, contrast_vals, colors, nSubj);
xlabel('Contrast [%]'); ylabel('GED Peak Amplitude [\Delta PR]');
title('[H6] Peak Amplitude (max at 75%)');
set(gca, 'FontSize', fontSize - 2);

report_trend('H5 Gamma Peak Freq', ged_peakFreq, contrast_vals, nSubj);
report_trend('H6 Gamma Peak Amp', ged_peakAmp, contrast_vals, nSubj);

% H6: also test quadratic (inverted U) — check if 75% > mean(25%,50%,100%)
amp75  = ged_peakAmp(3, :);
ampOth = nanmean(ged_peakAmp([1 2 4], :), 1);
[~, p_quad] = ttest(amp75 - ampOth);
fprintf('  H6 inverted-U test: 75%% vs mean(25%%,50%%,100%%): p = %.4f\n\n', p_quad);

saveas(fig5, fullfile(fig_dir, 'GCP_H5H6_gamma_freq_amp.png'));

%% ====================================================================
%  FIGURE 6: H7 — Gamma Frequency as Index of Visual Gain
%  ====================================================================

fig6 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H7] Gamma Frequency Relates to Oculomotor Dynamics', ...
    'FontSize', 18, 'FontWeight', 'bold');

% --- Panel (1,1): Gamma freq vs MS rate ---
subplot(2, 3, 1); hold on;
plot_scatter_by_cond(ged_peakFreq, msRate_subj, colors, condLabels, nSubj);
xlabel('GED Peak Frequency [Hz]'); ylabel('MS Rate [Hz]');
title('\gamma Freq vs MS Rate'); set(gca, 'FontSize', fontSize - 3);

% --- Panel (1,2): Gamma freq vs Velocity ---
subplot(2, 3, 2); hold on;
plot_scatter_by_cond(ged_peakFreq, vel2D_subj, colors, condLabels, nSubj);
xlabel('GED Peak Frequency [Hz]'); ylabel('Eye Velocity [px/s]');
title('\gamma Freq vs Eye Velocity'); set(gca, 'FontSize', fontSize - 3);

% --- Panel (1,3): Gamma freq vs Gaze deviation ---
subplot(2, 3, 3); hold on;
plot_scatter_by_cond(ged_peakFreq, gazeDev_subj, colors, condLabels, nSubj);
xlabel('GED Peak Frequency [Hz]'); ylabel('Gaze Deviation [px]');
title('\gamma Freq vs Gaze Dispersion'); set(gca, 'FontSize', fontSize - 3);

% --- Panel (2,1): Early window MS rate ---
subplot(2, 3, 4); hold on;
plot_early_late_bars(ms_early, ms_late, colors, condLabels, nSubj);
ylabel('MS Rate [Hz]'); title('MS Rate: Early vs Late');
set(gca, 'FontSize', fontSize - 3);

% --- Panel (2,2): Early window velocity ---
subplot(2, 3, 5); hold on;
plot_early_late_bars(vel_early, vel_late, colors, condLabels, nSubj);
ylabel('Eye Velocity [px/s]'); title('Velocity: Early vs Late');
set(gca, 'FontSize', fontSize - 3);

% --- Panel (2,3): CRF summary (gamma freq + MS rate normalised overlay) ---
subplot(2, 3, 6); hold on;
mu_freq = nanmean(ged_peakFreq, 2);
mu_ms   = nanmean(msRate_subj, 2);
sem_freq = nanstd(ged_peakFreq, [], 2) / sqrt(nSubj);
sem_ms   = nanstd(msRate_subj, [], 2) / sqrt(nSubj);

yyaxis left
errorbar(contrast_vals, mu_freq, sem_freq, 'b-o', 'LineWidth', 2.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'b', 'CapSize', 8);
ylabel('\gamma Peak Freq [Hz]');
set(gca, 'YColor', 'b');

yyaxis right
errorbar(contrast_vals, mu_ms, sem_ms, 'r-s', 'LineWidth', 2.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r', 'CapSize', 8);
ylabel('MS Rate [Hz]');
set(gca, 'YColor', 'r');

xlabel('Contrast [%]');
title('Dual CRF: \gamma Freq & MS Rate');
set(gca, 'XTick', contrast_vals, 'FontSize', fontSize - 3);
xlim([15 110]); grid on; box on;
legend({'\gamma Peak Freq', 'MS Rate'}, 'Location', 'best', 'FontSize', 10);

% Report correlations for H7
fprintf('H7 Correlations (across all subject-condition pairs):\n');
gf = ged_peakFreq(:); ms = msRate_subj(:); vl = vel2D_subj(:); gd = gazeDev_subj(:);
valid = ~isnan(gf) & ~isnan(ms);
if sum(valid) > 5
    [r, p] = corr(gf(valid), ms(valid));
    fprintf('  Gamma freq vs MS rate:     r = %.3f, p = %.4f\n', r, p);
end
valid = ~isnan(gf) & ~isnan(vl);
if sum(valid) > 5
    [r, p] = corr(gf(valid), vl(valid));
    fprintf('  Gamma freq vs Velocity:    r = %.3f, p = %.4f\n', r, p);
end
valid = ~isnan(gf) & ~isnan(gd);
if sum(valid) > 5
    [r, p] = corr(gf(valid), gd(valid));
    fprintf('  Gamma freq vs Gaze Dev:    r = %.3f, p = %.4f\n', r, p);
end
fprintf('\n');

saveas(fig6, fullfile(fig_dir, 'GCP_H7_gamma_oculomotor.png'));

%% ====================================================================
%  FIGURE 7: Grand Summary Dashboard
%  ====================================================================

fig7 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Hypothesis Summary — All Measures', ...
    'FontSize', 16, 'FontWeight', 'bold');

summary_data = {msRate_subj, vel2D_subj, constriction, ...
                gazeDev_subj, ged_peakFreq, ged_peakAmp};
summary_names = {'[H1] MS Rate', '[H2] Eye Velocity', '[H3] Constriction', ...
                 '[H4] Gaze Dev', '[H5] \gamma Peak Freq', '[H6] \gamma Peak Amp'};
summary_expect = {'decrease', 'decrease', 'decrease', ...
                  'increase', 'increase', 'inverted-U'};

for mi = 1:6
    subplot(2, 3, mi); hold on;
    dat = summary_data{mi};
    mu  = nanmean(dat, 2);
    sem = nanstd(dat, [], 2) / sqrt(nSubj);

    for c = 1:4
        bar(c, mu(c), 0.6, 'FaceColor', colors(c,:), 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    end
    errorbar(1:4, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);

    for s = 1:nSubj
        plot(1:4, dat(:, s), '-', 'Color', [0.6 0.6 0.6 0.3], 'LineWidth', 0.6);
        for c = 1:4
            if ~isnan(dat(c, s))
                scatter(c + (rand-0.5)*0.2, dat(c, s), 25, [0.3 0.3 0.3], ...
                    'filled', 'MarkerFaceAlpha', 0.5);
            end
        end
    end

    set(gca, 'XTick', 1:4, 'XTickLabel', {'25', '50', '75', '100'}, 'FontSize', 11);
    xlabel('Contrast [%]');
    title(sprintf('%s (%s)', summary_names{mi}, summary_expect{mi}), ...
        'FontSize', 11, 'FontWeight', 'bold');
    box on; grid on;
end

saveas(fig7, fullfile(fig_dir, 'GCP_hypotheses_summary.png'));

%% Done
fprintf('\n=== GCP Hypothesis Testing Complete ===\n');
fprintf('Figures saved to:\n  %s\n\n', fig_dir);
fprintf('Figures generated:\n');
fprintf('  1. H1_microsaccade_rate.png\n');
fprintf('  2. H2_eye_velocity.png\n');
fprintf('  3. H3_pupil_constriction.png\n');
fprintf('  4. H4_gaze_dispersion.png\n');
fprintf('  5. H5H6_gamma_freq_amp.png\n');
fprintf('  6. H7_gamma_oculomotor.png\n');
fprintf('  7. hypotheses_summary.png\n');

%% ====================================================================
%  HELPER FUNCTIONS
%  ====================================================================

function plot_crf(dat, contrast_vals, colors, nSubj)
    for s = 1:nSubj
        plot(contrast_vals, dat(:, s), '-o', 'Color', [0.7 0.7 0.7 0.5], ...
            'MarkerSize', 4, 'MarkerFaceColor', [0.7 0.7 0.7], 'LineWidth', 0.8);
    end
    mu  = nanmean(dat, 2);
    sem = nanstd(dat, [], 2) / sqrt(nSubj);
    errorbar(contrast_vals, mu, sem, 'k-o', 'LineWidth', 2.5, ...
        'MarkerSize', 10, 'MarkerFaceColor', 'k', 'CapSize', 8);
    for c = 1:4
        scatter(contrast_vals(c), mu(c), 120, colors(c,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
    set(gca, 'XTick', contrast_vals);
    xlim([15 110]); grid on; box on;
end

function r_mean = report_trend(label, dat, contrast_vals, nSubj)
    mu  = nanmean(dat, 2);
    sem = nanstd(dat, [], 2) / sqrt(nSubj);
    fprintf('%s:\n', label);
    for c = 1:4
        fprintf('  %3d%%: %.3f +/- %.3f\n', contrast_vals(c), mu(c), sem(c));
    end

    r_per_subj = nan(nSubj, 1);
    for s = 1:nSubj
        vals = dat(:, s);
        if sum(~isnan(vals)) >= 3
            r_per_subj(s) = corr(contrast_vals(:), vals, 'rows', 'complete');
        end
    end
    r_mean = nanmean(r_per_subj);
    [~, p] = ttest(r_per_subj);
    fprintf('  Linear trend: mean r = %.3f, t-test p = %.4f\n\n', r_mean, p);
end

function plot_scatter_by_cond(x_data, y_data, colors, condLabels, nSubj)
    all_x = []; all_y = [];

    % Subject lines first (behind scatter points)
    for s = 1:nSubj
        xv = x_data(:, s);
        yv = y_data(:, s);
        if sum(~isnan(xv) & ~isnan(yv)) >= 2
            plot(xv, yv, '-', 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.6, ...
                'HandleVisibility', 'off');
        end
    end

    % Scatter by condition
    h = gobjects(4, 1);
    for c = 1:4
        xv = x_data(c, :);
        yv = y_data(c, :);
        valid = ~isnan(xv) & ~isnan(yv);
        h(c) = scatter(xv(valid), yv(valid), 80, colors(c,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        all_x = [all_x, xv(valid)];
        all_y = [all_y, yv(valid)];
    end

    % Regression line
    valid = ~isnan(all_x) & ~isnan(all_y);
    if sum(valid) > 5
        [r, pval] = corr(all_x(valid)', all_y(valid)');
        p = polyfit(all_x(valid), all_y(valid), 1);
        xl = [min(all_x(valid)), max(all_x(valid))];
        xfit = linspace(xl(1), xl(2), 100);
        plot(xfit, polyval(p, xfit), 'k--', 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
        yl = ylim;
        text(xl(1) + 0.05*diff(xl), yl(2) - 0.08*diff(yl), ...
            sprintf('r = %.2f, p = %.3f', r, pval), 'FontSize', 10);
    end

    legend(h, condLabels, 'Location', 'best', 'FontSize', 9);
    grid on; box on;
end

function plot_early_late_bars(early_dat, late_dat, colors, condLabels, nSubj)
    bar_w = 0.7;
    h_leg = gobjects(2, 1);

    for c = 1:4
        idx_e = (c-1)*3 + 1;
        idx_l = (c-1)*3 + 2;

        mu_e  = nanmean(early_dat(c, :));
        sem_e = nanstd(early_dat(c, :)) / sqrt(nSubj);
        mu_l  = nanmean(late_dat(c, :));
        sem_l = nanstd(late_dat(c, :)) / sqrt(nSubj);

        hb_e = bar(idx_e, mu_e, bar_w, 'FaceColor', colors(c,:), ...
            'FaceAlpha', 0.4, 'EdgeColor', 'k', 'HandleVisibility', 'off');
        hb_l = bar(idx_l, mu_l, bar_w, 'FaceColor', colors(c,:), ...
            'FaceAlpha', 1.0, 'EdgeColor', 'k', 'HandleVisibility', 'off');
        errorbar(idx_e, mu_e, sem_e, 'k', 'LineWidth', 1.5, 'CapSize', 6, ...
            'HandleVisibility', 'off');
        errorbar(idx_l, mu_l, sem_l, 'k', 'LineWidth', 1.5, 'CapSize', 6, ...
            'HandleVisibility', 'off');

        if c == 1
            h_leg(1) = hb_e;
            h_leg(2) = hb_l;
            set(h_leg(1), 'HandleVisibility', 'on');
            set(h_leg(2), 'HandleVisibility', 'on');
        end
    end

    set(gca, 'XTick', [1.5, 4.5, 7.5, 10.5], 'XTickLabel', condLabels);
    legend(h_leg, {'Early (0-0.6s)', 'Late (1-2s)'}, 'Location', 'best', 'FontSize', 9);
    grid on; box on;
end
