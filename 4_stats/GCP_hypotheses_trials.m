%% GCP Hypothesis Testing — Trial-Level Analysis & Visualization
%
% Same hypotheses as GCP_hypotheses.m but using trial-level data instead
% of subject-level condition means. Visualises distributions with
% raincloud plots (violin + box + jittered scatter).
%
% Oculomotor (from gaze_matrix_trial.mat):
%   H1: Microsaccade rate decreases with increasing contrast
%   H2: Eye velocity increases post-stimulus, suppressed by contrast
%   H3: Pupil constriction amplitude increases monotonically with contrast
%   H4: BCEA increases with stimulus contrast
%
% Gamma oscillations (from GED trial-level data):
%   H5: Gamma peak frequency is higher at higher contrast
%   H6: Gamma peak amplitude peaks at 75% contrast (inverted U)
%
% Gamma as cortical index of visual gain:
%   H7: Gamma frequency relates to oculomotor dynamics
%       (uses subject-level means since gaze and EEG trials cannot
%        be reliably matched across preprocessing pipelines)
%
% Dependencies:
%   - Per-subject gaze_matrix_trial.mat  (gaze_fex output)
%   - GCP_eeg_GED_gamma_metrics.mat      (GED trial-level cell arrays)
%   - FieldTrip on path

clear; close all; clc

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP', 0);

nSubj = length(subjects);
fprintf('Subjects (N=%d): %s\n', nSubj, strjoin(subjects, ', '));
fprintf('No manual exclusions — outlier trials are rejected data-driven (median +/- 3 MAD).\n');

condLabels    = {'25%', '50%', '75%', '100%'};
contrast_vals = [25, 50, 75, 100];
fontSize      = 16;

data_dir = paths.features;
fig_dir = fullfile(paths.figures, 'hypotheses');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% ====================================================================
%  LOAD TRIAL-LEVEL GAZE DATA
%  ====================================================================
fprintf('Loading trial-level gaze data...\n');

gaze_Subject   = [];
gaze_Condition = [];
gaze_Trial     = [];
gaze_PctMSRate = [];
gaze_PctVel2D  = [];
gaze_PctPupil  = [];
gaze_PctBCEA   = [];
gaze_MSRate    = [];
gaze_Vel2D     = [];
gaze_PupilSize = [];
gaze_BCEA      = [];

for subj = 1:nSubj
    gazepath = fullfile(paths.features, subjects{subj}, 'gaze');
    trialfile = fullfile(gazepath, 'gaze_matrix_trial.mat');
    if ~exist(trialfile, 'file')
        fprintf('  %s: gaze_matrix_trial.mat not found, skipping.\n', subjects{subj});
        continue;
    end

    S = load(trialfile);
    sid = str2double(subjects{subj});

    cond_names = {'subj_data_gaze_trial_c25', 'subj_data_gaze_trial_c50', ...
                  'subj_data_gaze_trial_c75', 'subj_data_gaze_trial_c100'};

    for ci = 1:4
        if ~isfield(S, cond_names{ci}), continue; end
        d = S.(cond_names{ci});
        nTrl = numel(d.Trial);

        gaze_Subject   = [gaze_Subject;   repmat(sid, nTrl, 1)];
        gaze_Condition = [gaze_Condition;  repmat(ci, nTrl, 1)];
        gaze_Trial     = [gaze_Trial;     d.Trial(:)];

        gaze_MSRate    = [gaze_MSRate;    d.MSRate(:)];
        gaze_Vel2D     = [gaze_Vel2D;     d.Vel2D(:)];
        gaze_PupilSize = [gaze_PupilSize; d.PupilSize(:)];
        gaze_BCEA      = [gaze_BCEA;      d.BCEA(:)];

        if isfield(d, 'PctMSRate')
            gaze_PctMSRate = [gaze_PctMSRate; d.PctMSRate(:)];
        else
            gaze_PctMSRate = [gaze_PctMSRate; nan(nTrl, 1)];
        end
        if isfield(d, 'PctVel2D')
            gaze_PctVel2D = [gaze_PctVel2D; d.PctVel2D(:)];
        else
            gaze_PctVel2D = [gaze_PctVel2D; nan(nTrl, 1)];
        end
        if isfield(d, 'PctPupilSize')
            gaze_PctPupil = [gaze_PctPupil; d.PctPupilSize(:)];
        else
            gaze_PctPupil = [gaze_PctPupil; nan(nTrl, 1)];
        end
        if isfield(d, 'PctBCEA')
            gaze_PctBCEA = [gaze_PctBCEA; d.PctBCEA(:)];
        else
            gaze_PctBCEA = [gaze_PctBCEA; nan(nTrl, 1)];
        end
    end
    fprintf('  Subject %s (%d/%d)\n', subjects{subj}, subj, nSubj);
end

fprintf('  Gaze: %d total trials loaded.\n', numel(gaze_Subject));

%% ====================================================================
%  LOAD TRIAL-LEVEL GED DATA
%  ====================================================================
fprintf('Loading trial-level GED gamma data...\n');

ged = load(fullfile(data_dir, 'GCP_eeg_GED_gamma_metrics.mat'));

ged_Subject   = [];
ged_Condition = [];
ged_Trial     = [];
ged_PeakFreq  = [];
ged_PeakAmp   = [];
ged_BBPower   = [];
ged_Centroid  = [];
ged_AUC       = [];

ged_subj_order = ged.subjects;

for subj = 1:nSubj
    gi = find(strcmp(ged_subj_order, subjects{subj}));
    if isempty(gi), continue; end
    sid = str2double(subjects{subj});

    for cond = 1:4
        pf = ged.trl_peak_freq{cond, gi};
        if isempty(pf), continue; end
        nTrl = numel(pf);

        ged_Subject   = [ged_Subject;   repmat(sid, nTrl, 1)];
        ged_Condition = [ged_Condition;  repmat(cond, nTrl, 1)];
        ged_Trial     = [ged_Trial;     (1:nTrl)'];

        ged_PeakFreq  = [ged_PeakFreq;  pf(:)];
        ged_PeakAmp   = [ged_PeakAmp;   ged.trl_peak_amp{cond, gi}(:)];
        ged_BBPower   = [ged_BBPower;   ged.trl_bb_power{cond, gi}(:)];
        ged_Centroid  = [ged_Centroid;  ged.trl_centroid{cond, gi}(:)];
        ged_AUC       = [ged_AUC;       ged.trl_auc{cond, gi}(:)];
    end
end

fprintf('  GED: %d total trials loaded.\n', numel(ged_Subject));

%% ====================================================================
%  OUTLIER REJECTION (median ± 3 MAD per condition)
%  ====================================================================
fprintf('Removing outlier trials (median +/- 3 MAD per condition)...\n');

% Gaze metrics
[gaze_PctMSRate, n1] = reject_outliers_by_cond(gaze_PctMSRate, gaze_Condition);
[gaze_PctVel2D,  n2] = reject_outliers_by_cond(gaze_PctVel2D,  gaze_Condition);
[gaze_PctPupil,  n3] = reject_outliers_by_cond(gaze_PctPupil,  gaze_Condition);
[gaze_PctBCEA,   n4] = reject_outliers_by_cond(gaze_PctBCEA,   gaze_Condition);
[gaze_MSRate,    n5] = reject_outliers_by_cond(gaze_MSRate,    gaze_Condition);
[gaze_Vel2D,     n6] = reject_outliers_by_cond(gaze_Vel2D,     gaze_Condition);
[gaze_PupilSize, n7] = reject_outliers_by_cond(gaze_PupilSize, gaze_Condition);
[gaze_BCEA,      n8] = reject_outliers_by_cond(gaze_BCEA,      gaze_Condition);
fprintf('  Gaze outliers removed: PctMS=%d, PctVel=%d, PctPup=%d, PctBCEA=%d, MS=%d, Vel=%d, Pup=%d, BCEA=%d\n', ...
    n1, n2, n3, n4, n5, n6, n7, n8);

% GED metrics
[ged_PeakFreq, n1] = reject_outliers_by_cond(ged_PeakFreq, ged_Condition);
[ged_PeakAmp,  n2] = reject_outliers_by_cond(ged_PeakAmp,  ged_Condition);
[ged_BBPower,  n3] = reject_outliers_by_cond(ged_BBPower,  ged_Condition);
[ged_Centroid, n4] = reject_outliers_by_cond(ged_Centroid, ged_Condition);
[ged_AUC,      n5] = reject_outliers_by_cond(ged_AUC,      ged_Condition);
fprintf('  GED outliers removed: PeakFreq=%d, PeakAmp=%d, BBPow=%d, Centroid=%d, AUC=%d\n', ...
    n1, n2, n3, n4, n5);

%% ====================================================================
%  COMPUTE SUBJECT-LEVEL MEANS (for H7 cross-modal analysis)
%  ====================================================================

subj_pctMSRate = nan(4, nSubj);
subj_pctVel2D  = nan(4, nSubj);
subj_pctBCEA   = nan(4, nSubj);
subj_gedFreq   = nan(4, nSubj);

for s = 1:nSubj
    sid = str2double(subjects{s});
    for c = 1:4
        idx_g = gaze_Subject == sid & gaze_Condition == c;
        if any(idx_g)
            subj_pctMSRate(c, s) = nanmean(gaze_PctMSRate(idx_g));
            subj_pctVel2D(c, s)  = nanmean(gaze_PctVel2D(idx_g));
            subj_pctBCEA(c, s)   = nanmean(gaze_PctBCEA(idx_g));
        end

        idx_e = ged_Subject == sid & ged_Condition == c;
        if any(idx_e)
            subj_gedFreq(c, s) = nanmean(ged_PeakFreq(idx_e));
        end
    end
end

% No second-pass outlier rejection on subject means — trial-level MAD is sufficient

%% ====================================================================
%  FIGURE 1: H1 — Microsaccade Rate (trial-level)
%  ====================================================================
close all
fprintf('\nGenerating trial-level figures...\n');

fig1 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H1] Microsaccade Rate Decreases with Contrast — Trial-Level', ...
    'FontSize', 18, 'FontWeight', 'bold');

subplot(1, 2, 1);
plot_raincloud_trial(gaze_PctMSRate, gaze_Condition, colors, condLabels, ...
    'MS Rate [% change]');

subplot(1, 2, 2);
plot_raincloud_trial(gaze_MSRate, gaze_Condition, colors, condLabels, ...
    'MS Rate [Hz]');

report_trend_trial('H1 MS Rate (% change)', gaze_PctMSRate, gaze_Condition, ...
    gaze_Subject, contrast_vals);

saveas(fig1, fullfile(fig_dir, 'GCP_H1_microsaccade_rate_trials.png'));

%% ====================================================================
%  FIGURE 2: H2 — Eye Velocity (trial-level)
%  ====================================================================

fig2 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H2] Eye Velocity: Contrast-Dependent Suppression — Trial-Level', ...
    'FontSize', 18, 'FontWeight', 'bold');

subplot(1, 2, 1);
plot_raincloud_trial(gaze_PctVel2D, gaze_Condition, colors, condLabels, ...
    'Eye Velocity [% change]');

subplot(1, 2, 2);
plot_raincloud_trial(gaze_Vel2D, gaze_Condition, colors, condLabels, ...
    'Eye Velocity [px/s]');

report_trend_trial('H2 Velocity (% change)', gaze_PctVel2D, gaze_Condition, ...
    gaze_Subject, contrast_vals);

saveas(fig2, fullfile(fig_dir, 'GCP_H2_eye_velocity_trials.png'));

%% ====================================================================
%  FIGURE 3: H3 — Pupil Size (trial-level)
%  ====================================================================

fig3 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H3] Pupil Constriction Scales with Contrast — Trial-Level', ...
    'FontSize', 18, 'FontWeight', 'bold');

subplot(1, 2, 1);
plot_raincloud_trial(gaze_PctPupil, gaze_Condition, colors, condLabels, ...
    'Pupil Size [% change]');

subplot(1, 2, 2);
plot_raincloud_trial(gaze_PupilSize, gaze_Condition, colors, condLabels, ...
    'Pupil Size [a.u.]');

report_trend_trial('H3 Pupil (% change)', gaze_PctPupil, gaze_Condition, ...
    gaze_Subject, contrast_vals);

saveas(fig3, fullfile(fig_dir, 'GCP_H3_pupil_trials.png'));

%% ====================================================================
%  FIGURE 4: H4 — BCEA (trial-level)
%  ====================================================================

fig4 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H4] BCEA Increases with Stimulus Contrast — Trial-Level', ...
    'FontSize', 18, 'FontWeight', 'bold');

subplot(1, 2, 1);
plot_raincloud_trial(gaze_PctBCEA, gaze_Condition, colors, condLabels, ...
    'BCEA [% change]');

subplot(1, 2, 2);
plot_raincloud_trial(gaze_BCEA, gaze_Condition, colors, condLabels, ...
    'BCEA [px^2]');

report_trend_trial('H4 BCEA (% change)', gaze_PctBCEA, gaze_Condition, ...
    gaze_Subject, contrast_vals);

saveas(fig4, fullfile(fig_dir, 'GCP_H4_bcea_trials.png'));

%% ====================================================================
%  FIGURE 5: H5 + H6 — Gamma Peak Frequency & Power (trial-level)
%  ====================================================================

fig5 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H5] Gamma Peak Frequency & [H6] Peak Power — Trial-Level', ...
    'FontSize', 18, 'FontWeight', 'bold');

subplot(2, 3, 1);
plot_raincloud_trial(ged_PeakFreq, ged_Condition, colors, condLabels, ...
    'GED Peak Frequency [Hz]');
title('[H5] GED Peak Frequency', 'FontSize', 13, 'FontWeight', 'bold');

subplot(2, 3, 2);
plot_raincloud_trial(ged_Centroid, ged_Condition, colors, condLabels, ...
    'Spectral Centroid [Hz]');
title('[H5] Spectral Centroid', 'FontSize', 13, 'FontWeight', 'bold');

subplot(2, 3, 3);
plot_raincloud_trial(ged_AUC, ged_Condition, colors, condLabels, ...
    'AUC [a.u.]');
title('GED AUC', 'FontSize', 13, 'FontWeight', 'bold');

subplot(2, 3, 4);
plot_raincloud_trial(ged_PeakAmp, ged_Condition, colors, condLabels, ...
    'Peak Amplitude [\Delta PR]');
title('[H6] GED Peak Amplitude', 'FontSize', 13, 'FontWeight', 'bold');

subplot(2, 3, 5);
plot_raincloud_trial(ged_BBPower, ged_Condition, colors, condLabels, ...
    'Broadband Power [ratio]');
title('[H6] GED Broadband Power', 'FontSize', 13, 'FontWeight', 'bold');

% Panel 6: empty — inverted-U test summary
subplot(2, 3, 6); hold on;
axis off;
txt = {'\bf{Trend tests (trial-level):}'};
for pp = 1:3
    switch pp
        case 1, pdat = ged_PeakAmp;  plbl = 'GED Peak Amp';  pcond = ged_Condition;
        case 2, pdat = ged_BBPower;   plbl = 'GED BB Power';  pcond = ged_Condition;
        case 3, pdat = ged_Centroid;  plbl = 'Centroid';       pcond = ged_Condition;
    end
    val75  = pdat(pcond == 3);
    valOth = [pdat(pcond == 1); pdat(pcond == 2); pdat(pcond == 4)];
    val75  = val75(~isnan(val75));
    valOth = valOth(~isnan(valOth));
    [~, p_rs] = ranksum(val75, valOth);
    txt{end+1} = sprintf('%s: 75%% vs others p = %.4f', plbl, p_rs);
end
text(0.1, 0.5, txt, 'FontSize', 12, 'VerticalAlignment', 'middle');

report_trend_trial('H5 GED Peak Freq', ged_PeakFreq, ged_Condition, ...
    ged_Subject, contrast_vals);
report_trend_trial('H6 GED Peak Amp', ged_PeakAmp, ged_Condition, ...
    ged_Subject, contrast_vals);
report_trend_trial('H6 GED BB Power', ged_BBPower, ged_Condition, ...
    ged_Subject, contrast_vals);

saveas(fig5, fullfile(fig_dir, 'GCP_H5H6_gamma_trials.png'));

%% ====================================================================
%  FIGURE 6: H7 — Gamma–Oculomotor Relationship (subject-level means)
%  ====================================================================

fig6 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H7] Gamma Frequency vs Oculomotor Dynamics (subject means from trials)', ...
    'FontSize', 18, 'FontWeight', 'bold');

subplot(2, 3, 1); hold on;
plot_scatter_by_cond(subj_gedFreq, subj_pctMSRate, colors, condLabels, nSubj);
xlabel('GED Peak Frequency [Hz]'); ylabel('MS Rate [% change]');
title('\gamma Freq vs MS Rate'); set(gca, 'FontSize', fontSize - 3);

subplot(2, 3, 2); hold on;
plot_scatter_by_cond(subj_gedFreq, subj_pctVel2D, colors, condLabels, nSubj);
xlabel('GED Peak Frequency [Hz]'); ylabel('Eye Velocity [% change]');
title('\gamma Freq vs Velocity'); set(gca, 'FontSize', fontSize - 3);

subplot(2, 3, 3); hold on;
plot_scatter_by_cond(subj_gedFreq, subj_pctBCEA, colors, condLabels, nSubj);
xlabel('GED Peak Frequency [Hz]'); ylabel('BCEA [% change]');
title('\gamma Freq vs BCEA'); set(gca, 'FontSize', fontSize - 3);

% Dual CRF overlay (from trial-level data, consistent with rainclouds)
subplot(2, 3, 4); hold on;
mu_freq  = nan(4,1); sem_freq = nan(4,1);
mu_ms    = nan(4,1); sem_ms   = nan(4,1);
for c = 1:4
    vf = ged_PeakFreq(ged_Condition == c);  vf = vf(~isnan(vf));
    mu_freq(c)  = mean(vf);  sem_freq(c) = std(vf)/sqrt(numel(vf));
    vm = gaze_PctMSRate(gaze_Condition == c); vm = vm(~isnan(vm));
    mu_ms(c)    = mean(vm);  sem_ms(c)   = std(vm)/sqrt(numel(vm));
end

yyaxis left
errorbar(contrast_vals, mu_freq, sem_freq, 'b-o', 'LineWidth', 2.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'b', 'CapSize', 8);
ylabel('\gamma Peak Freq [Hz]');
set(gca, 'YColor', 'b');

yyaxis right
errorbar(contrast_vals, mu_ms, sem_ms, 'r-s', 'LineWidth', 2.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r', 'CapSize', 8);
ylabel('MS Rate [% change]');
set(gca, 'YColor', 'r');

xlabel('Contrast [%]');
title('Dual CRF: \gamma Freq & MS Rate');
set(gca, 'XTick', contrast_vals, 'FontSize', fontSize - 3);
xlim([15 110]);  box on;
legend({'\gamma Peak Freq', 'MS Rate'}, 'Location', 'best', 'FontSize', 10);

% Trial counts per condition
subplot(2, 3, 5); hold on;
for c = 1:4
    n_gaze = sum(gaze_Condition == c & ~isnan(gaze_PctMSRate));
    n_ged  = sum(ged_Condition == c & ~isnan(ged_PeakFreq));
    bar(c - 0.15, n_gaze, 0.3, 'FaceColor', colors(c,:), 'FaceAlpha', 0.5, ...
        'EdgeColor', 'k');
    bar(c + 0.15, n_ged, 0.3, 'FaceColor', colors(c,:), 'FaceAlpha', 1.0, ...
        'EdgeColor', 'k');
end
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 12);
ylabel('N trials'); title('Trial Counts'); legend({'Gaze', 'GED'}, 'FontSize', 10);
 box on;

% Correlation summary text
subplot(2, 3, 6); hold on; axis off;
gf = subj_gedFreq(:); ms = subj_pctMSRate(:);
vl = subj_pctVel2D(:); bc = subj_pctBCEA(:);
txt = {'\bf{H7 Correlations (subject-condition pairs):}'};

valid = ~isnan(gf) & ~isnan(ms);
if sum(valid) > 5
    [r, p] = corr(gf(valid), ms(valid));
    txt{end+1} = sprintf('\\gamma freq vs MS rate:  r = %.3f, p = %.4f', r, p);
end
valid = ~isnan(gf) & ~isnan(vl);
if sum(valid) > 5
    [r, p] = corr(gf(valid), vl(valid));
    txt{end+1} = sprintf('\\gamma freq vs Velocity: r = %.3f, p = %.4f', r, p);
end
valid = ~isnan(gf) & ~isnan(bc);
if sum(valid) > 5
    [r, p] = corr(gf(valid), bc(valid));
    txt{end+1} = sprintf('\\gamma freq vs BCEA:     r = %.3f, p = %.4f', r, p);
end
text(0.05, 0.5, txt, 'FontSize', 12, 'VerticalAlignment', 'middle');

saveas(fig6, fullfile(fig_dir, 'GCP_H7_gamma_oculomotor_trials.png'));

%% ====================================================================
%  FIGURE 7: Summary Dashboard (trial-level distributions)
%  ====================================================================

fig7 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Hypothesis Summary — Trial-Level Distributions', ...
    'FontSize', 16, 'FontWeight', 'bold');

summary_data   = {gaze_PctMSRate, gaze_PctVel2D, gaze_PctPupil, ...
                  gaze_PctBCEA, ged_PeakFreq, ged_PeakAmp};
summary_cond   = {gaze_Condition, gaze_Condition, gaze_Condition, ...
                  gaze_Condition, ged_Condition, ged_Condition};
summary_names  = {'[H1] MS Rate [%\Delta]', '[H2] Velocity [%\Delta]', ...
                  '[H3] Pupil [%\Delta]', '[H4] BCEA [%\Delta]', ...
                  '[H5] \gamma Peak Freq', '[H6] \gamma Peak Amp'};
summary_expect = {'decrease', 'decrease', 'decrease', ...
                  'increase', 'increase', 'inverted-U'};

for mi = 1:6
    subplot(2, 3, mi); hold on;
    dat  = summary_data{mi};
    cnd  = summary_cond{mi};

    for c = 1:4
        vals = dat(cnd == c);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end

        mu  = mean(vals);
        sem = std(vals) / sqrt(numel(vals));

        xJit = c + (rand(size(vals)) - 0.5) * 0.3;
        scatter(xJit, vals, 8, colors(c,:), 'filled', ...
            'MarkerFaceAlpha', 0.15, 'HandleVisibility', 'off');

        bar(c, mu, 0.6, 'FaceColor', colors(c,:), 'EdgeColor', 'k', ...
            'FaceAlpha', 0.6);
        errorbar(c, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, ...
            'CapSize', 8);
    end

    set(gca, 'XTick', 1:4, 'XTickLabel', {'25', '50', '75', '100'}, 'FontSize', 11);
    xlabel('Contrast [%]');
    title(sprintf('%s (%s)', summary_names{mi}, summary_expect{mi}), ...
        'FontSize', 11, 'FontWeight', 'bold');
    box on; 
end

saveas(fig7, fullfile(fig_dir, 'GCP_hypotheses_summary_trials.png'));

%% Done
fprintf('\n=== GCP Trial-Level Hypothesis Testing Complete ===\n');
fprintf('Figures saved to:\n  %s\n\n', fig_dir);
fprintf('Trial counts:\n');
fprintf('  Gaze trials:  %d\n', numel(gaze_Subject));
fprintf('  GED trials:   %d\n', numel(ged_Subject));
fprintf('  Subjects:     %d\n', nSubj);

%% ====================================================================
%  HELPER FUNCTIONS
%  ====================================================================

function plot_raincloud_trial(data, cond_vec, colors, condLabels, y_label)
    hold on;
    for c = 1:4
        vals = data(cond_vec == c);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end

        % Kernel density (violin on left side)
        if numel(vals) >= 3
            [f_dens, xi] = ksdensity(vals);
            f_dens = f_dens / max(f_dens) * 0.3;
            patch(c - f_dens - 0.05, xi, colors(c,:), ...
                'FaceAlpha', 0.3, 'EdgeColor', colors(c,:), 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end

        % Jittered scatter (right side)
        xJit = c + 0.15 + (rand(size(vals)) - 0.5) * 0.15;
        scatter(xJit, vals, 15, colors(c,:), 'filled', ...
            'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % Box: median, IQR, whiskers
        q25 = prctile(vals, 25);
        q50 = prctile(vals, 50);
        q75 = prctile(vals, 75);
        iqr_val = q75 - q25;
        whi_lo  = max(min(vals), q25 - 1.5 * iqr_val);
        whi_hi  = min(max(vals), q75 + 1.5 * iqr_val);

        bw = 0.12;
        fill([c-bw c+bw c+bw c-bw], [q25 q25 q75 q75], colors(c,:), ...
            'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
        plot([c-bw c+bw], [q50 q50], 'k-', 'LineWidth', 2, ...
            'HandleVisibility', 'off');
        plot([c c], [whi_lo q25], 'k-', 'LineWidth', 1, ...
            'HandleVisibility', 'off');
        plot([c c], [q75 whi_hi], 'k-', 'LineWidth', 1, ...
            'HandleVisibility', 'off');

        % Mean marker
        scatter(c, mean(vals), 60, 'k', 'diamond', 'filled', ...
            'HandleVisibility', 'off');
    end

    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 13);
    ylabel(y_label);
     box on;
end

function report_trend_trial(label, data, cond_vec, subj_vec, contrast_vals)
    fprintf('%s (trial-level):\n', label);
    for c = 1:4
        vals = data(cond_vec == c);
        vals = vals(~isnan(vals));
        fprintf('  %3d%%: mean = %.3f, median = %.3f, N = %d\n', ...
            contrast_vals(c), mean(vals), median(vals), numel(vals));
    end

    % Per-subject linear trend
    uSubj = unique(subj_vec);
    r_per_subj = nan(numel(uSubj), 1);
    for si = 1:numel(uSubj)
        subj_means = nan(4, 1);
        for c = 1:4
            idx = subj_vec == uSubj(si) & cond_vec == c;
            vals = data(idx);
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                subj_means(c) = mean(vals);
            end
        end
        if sum(~isnan(subj_means)) >= 3
            r_per_subj(si) = corr(contrast_vals(:), subj_means, 'rows', 'complete');
        end
    end
    r_mean = nanmean(r_per_subj);
    [~, p] = ttest(r_per_subj);
    fprintf('  Linear trend (per-subj means): mean r = %.3f, t-test p = %.4f\n\n', r_mean, p);
end

function plot_scatter_by_cond(x_data, y_data, colors, condLabels, nSubj)
    all_x = []; all_y = [];

    for s = 1:nSubj
        xv = x_data(:, s);
        yv = y_data(:, s);
        if sum(~isnan(xv) & ~isnan(yv)) >= 2
            plot(xv, yv, '-', 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.6, ...
                'HandleVisibility', 'off');
        end
    end

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
     box on;
end

function [data_out, n_removed] = reject_outliers_by_cond(data, cond_vec)
    data_out  = data;
    n_removed = 0;
    for c = 1:4
        idx = cond_vec == c;
        vals = data(idx);
        valid = ~isnan(vals);
        if sum(valid) < 5, continue; end
        med_val = median(vals(valid));
        mad_val = mad(vals(valid), 1);
        if mad_val == 0, continue; end
        outlier = valid & (abs(vals - med_val) > 3 * mad_val);
        n_removed = n_removed + sum(outlier);
        vals(outlier) = NaN;
        data_out(idx) = vals;
    end
end

function [data_out, n_removed] = reject_outliers_by_cond_mat(data_mat)
    data_out  = data_mat;
    n_removed = 0;
    for c = 1:size(data_mat, 1)
        vals  = data_mat(c, :);
        valid = ~isnan(vals);
        if sum(valid) < 3, continue; end
        med_val = median(vals(valid));
        mad_val = mad(vals(valid), 1);
        if mad_val == 0, continue; end
        outlier = valid & (abs(vals - med_val) > 3 * mad_val);
        n_removed = n_removed + sum(outlier);
        data_out(c, outlier) = NaN;
    end
end
