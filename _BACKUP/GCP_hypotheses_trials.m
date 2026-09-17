%% GCP Hypothesis Testing: Trial-Level Analysis & Visualization
%
% Trial-level tests of the current gamma / oculomotor hypotheses.
% Prefers GCP_merged_data_trials.csv so gamma and gaze share ID/Condition/Trial.
%
% Oculomotor (baselined % change; Vel2D_bl = eye velocity):
%   H1: Microsaccade rate (MSRate_bl) decreases with increasing contrast
%   H2: Eye velocity (Vel2D_bl) increases with contrast
%   H3: BCEA (BCEA_bl) increases with stimulus contrast
%
% Gamma oscillations (GED per-trial peaks):
%   H4: Peak gamma frequency increases with contrast
%   H5: Peak gamma power increases with contrast
%
% Gamma as cortical index of visual gain:
%   H6: Higher gamma frequency/power co-occurs with lower MS rate,
%       higher BCEA, and higher eye velocity (matched trials)
%
% Dependencies:
%   - GCP_merged_data_trials.csv  (preferred; from GCP_master_matrix_trials.m)
%   - Fallback: gaze_matrix_trial.mat + GCP_eeg_GED.mat
%   - FieldTrip on path (via startup)

clear; close all; clc

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP', 0);
subjects = gcp_subject_inclusion(subjects, paths);

nSubj = length(subjects);
fprintf('[STATS HYP TRIALS] Subjects (N=%d, GED cohort): %s\n', nSubj, strjoin(subjects, ', '));
fprintf('[STATS HYP TRIALS] Outlier trials rejected per condition (IQR, 1.5x; matches rainclouds).\n');

condLabels    = {'25%', '50%', '75%', '100%'};
contrast_vals = [25, 50, 75, 100];
fontSize      = 16;

data_dir = paths.features;
fig_dir = fullfile(paths.figures, 'hypotheses');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% LOAD TRIAL DATA
merged_csv = fullfile(data_dir, 'GCP_merged_data_trials.csv');
merged_mat = fullfile(data_dir, 'GCP_merged_data_trials.mat');

if isfile(merged_csv)
    fprintf('[STATS HYP TRIALS] Loading merged trial table:\n  %s\n', merged_csv);
    T = readtable(merged_csv);
elseif isfile(merged_mat)
    fprintf('[STATS HYP TRIALS] Loading merged trial MAT:\n  %s\n', merged_mat);
    S = load(merged_mat);
    if isfield(S, 'GCP_merged_table_trials')
        T = S.GCP_merged_table_trials;
    elseif isfield(S, 'GCP_merged_data_trials')
        T = struct2table(S.GCP_merged_data_trials);
    else
        error('GCP_hypotheses_trials:BadMergeMat', ...
            'No expected merged variable in %s', merged_mat);
    end
else
    fprintf('[STATS HYP TRIALS] Merged table missing; building trial vectors from gaze + GED files.\n');
    T = build_trial_table_from_sources(subjects, paths, data_dir);
end

T = restrict_to_subjects(T, subjects);
T = standardize_condition_codes(T);

ms_var   = pick_var(T, {'Gaze_MSRate_bl', 'MSRate_bl', 'Gaze_dBMSRate', 'dBMSRate'});
vel_var  = pick_var(T, {'Gaze_Vel2D_bl', 'Vel2D_bl', 'Gaze_dBVel2D', 'dBVel2D'});
bcea_var = pick_var(T, {'Gaze_BCEA_bl', 'BCEA_bl', 'Gaze_dBBCEA', 'dBBCEA'});
freq_var = pick_var(T, {'GED_GammaFrequency', 'GammaFrequency', 'GED_PeakFrequency'});
pow_var  = pick_var(T, {'GED_GammaPower', 'GammaPower', 'GED_PeakAmplitude'});
id_var   = pick_var(T, {'ID'});
cond_var = pick_var(T, {'Condition'});

need = {id_var, cond_var, ms_var, vel_var, bcea_var, freq_var, pow_var};
if any(cellfun(@isempty, need))
    error('GCP_hypotheses_trials:MissingCols', ...
        ['Merged/source table missing required columns. Need ID, Condition, ', ...
         'MSRate_bl, Vel2D_bl, BCEA_bl, GammaFrequency, GammaPower.']);
end

fprintf('Using columns: MS=%s, Vel=%s, BCEA=%s, Freq=%s, Power=%s\n', ...
    ms_var, vel_var, bcea_var, freq_var, pow_var);

gaze_Subject   = double(T.(id_var));
gaze_Condition = double(T.(cond_var));
gaze_MSRate_bl = double(T.(ms_var));
gaze_Vel2D_bl  = double(T.(vel_var));
gaze_BCEA_bl   = double(T.(bcea_var));
ged_Subject    = gaze_Subject;
ged_Condition  = gaze_Condition;
ged_PeakFreq   = double(T.(freq_var));
ged_PeakPower  = double(T.(pow_var));

fprintf('  Loaded %d trial rows.\n', height(T));

%% OUTLIER REJECTION (IQR per condition; same spirit as GCP_stats_rainclouds.py)
fprintf('Removing outlier trials (IQR 1.5x per condition)...\n');
[gaze_MSRate_bl, n1] = reject_outliers_iqr_by_cond(gaze_MSRate_bl, gaze_Condition);
[gaze_Vel2D_bl,  n2] = reject_outliers_iqr_by_cond(gaze_Vel2D_bl,  gaze_Condition);
[gaze_BCEA_bl,   n3] = reject_outliers_iqr_by_cond(gaze_BCEA_bl,   gaze_Condition);
[ged_PeakFreq,   n4] = reject_outliers_iqr_by_cond(ged_PeakFreq,   ged_Condition);
[ged_PeakPower,  n5] = reject_outliers_iqr_by_cond(ged_PeakPower,  ged_Condition);
fprintf('  Outliers removed: MS=%d, Vel=%d, BCEA=%d, Freq=%d, Power=%d\n', ...
    n1, n2, n3, n4, n5);

%% SUBJECT x CONDITION MEANS (for H6 subject-level panels / dual CRF)
subj_MSRate_bl = nan(4, nSubj);
subj_Vel2D_bl  = nan(4, nSubj);
subj_BCEA_bl   = nan(4, nSubj);
subj_gedFreq   = nan(4, nSubj);
subj_gedPower  = nan(4, nSubj);

for s = 1:nSubj
    sid = str2double(subjects{s});
    for c = 1:4
        idx_g = gaze_Subject == sid & gaze_Condition == c;
        if any(idx_g)
            subj_MSRate_bl(c, s) = mean(gaze_MSRate_bl(idx_g), 'omitnan');
            subj_Vel2D_bl(c, s)  = mean(gaze_Vel2D_bl(idx_g), 'omitnan');
            subj_BCEA_bl(c, s)   = mean(gaze_BCEA_bl(idx_g), 'omitnan');
        end
        idx_e = ged_Subject == sid & ged_Condition == c;
        if any(idx_e)
            subj_gedFreq(c, s)  = mean(ged_PeakFreq(idx_e), 'omitnan');
            subj_gedPower(c, s) = mean(ged_PeakPower(idx_e), 'omitnan');
        end
    end
end

%% TRIAL-LEVEL LMEs (contrast coded as continuous percent)
fprintf('\n=== Trial-level LMEs (Y ~ Contrast + (1|ID)) ===\n');
lme_rows = empty_lme_table();
lme_rows = [lme_rows; fit_contrast_lme('H1_MSRate_bl', gaze_MSRate_bl, gaze_Condition, gaze_Subject, contrast_vals)]; %#ok<AGROW>
lme_rows = [lme_rows; fit_contrast_lme('H2_Vel2D_bl',  gaze_Vel2D_bl,  gaze_Condition, gaze_Subject, contrast_vals)]; %#ok<AGROW>
lme_rows = [lme_rows; fit_contrast_lme('H3_BCEA_bl',   gaze_BCEA_bl,   gaze_Condition, gaze_Subject, contrast_vals)]; %#ok<AGROW>
lme_rows = [lme_rows; fit_contrast_lme('H4_GammaFreq', ged_PeakFreq,   ged_Condition,  ged_Subject,  contrast_vals)]; %#ok<AGROW>
lme_rows = [lme_rows; fit_contrast_lme('H5_GammaPower',ged_PeakPower,  ged_Condition,  ged_Subject,  contrast_vals)]; %#ok<AGROW>
disp(lme_rows);

% H6: gamma ~ oculomotor + contrast on matched trials
fprintf('\n=== H6 trial-level LMEs (Gamma ~ Gaze + Contrast + (1|ID)) ===\n');
h6_rows = empty_lme_table();
h6_specs = { ...
    'Freq_vs_MSRate',   ged_PeakFreq,  gaze_MSRate_bl; ...
    'Freq_vs_BCEA',     ged_PeakFreq,  gaze_BCEA_bl; ...
    'Freq_vs_Vel2D',    ged_PeakFreq,  gaze_Vel2D_bl; ...
    'Power_vs_MSRate',  ged_PeakPower, gaze_MSRate_bl; ...
    'Power_vs_BCEA',    ged_PeakPower, gaze_BCEA_bl; ...
    'Power_vs_Vel2D',   ged_PeakPower, gaze_Vel2D_bl};
for i = 1:size(h6_specs, 1)
    h6_rows = [h6_rows; fit_gamma_gaze_lme( ...
        h6_specs{i, 1}, h6_specs{i, 2}, h6_specs{i, 3}, ...
        gaze_Condition, gaze_Subject, contrast_vals)]; %#ok<AGROW>
end
disp(h6_rows);

writetable([lme_rows; h6_rows], fullfile(data_dir, 'GCP_hypotheses_trials_lme.csv'));

%% FIGURE 1: H1 Microsaccade Rate
close all
fprintf('\nGenerating trial-level figures...\n');

fig1 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H1] Microsaccade Rate Decreases with Contrast (Trial-Level)', ...
    'FontSize', 18, 'FontWeight', 'bold');
plot_raincloud_trial(gaze_MSRate_bl, gaze_Condition, colors, condLabels, ...
    'Microsaccade Rate [%]');
report_trend_trial('H1 MS Rate (%)', gaze_MSRate_bl, gaze_Condition, ...
    gaze_Subject, contrast_vals);
exportgraphics(fig1, fullfile(fig_dir, 'GCP_H1_microsaccade_rate_trials.png'), 'Resolution', 600);

%% FIGURE 2: H2 Eye velocity
fig2 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H2] Eye Velocity Increases with Contrast (Trial-Level)', ...
    'FontSize', 18, 'FontWeight', 'bold');
plot_raincloud_trial(gaze_Vel2D_bl, gaze_Condition, colors, condLabels, ...
    'Eye Velocity [%]');
report_trend_trial('H2 Velocity (%)', gaze_Vel2D_bl, gaze_Condition, ...
    gaze_Subject, contrast_vals);
exportgraphics(fig2, fullfile(fig_dir, 'GCP_H2_eye_velocity_trials.png'), 'Resolution', 600);

%% FIGURE 3: H3 BCEA
fig3 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H3] BCEA Increases with Stimulus Contrast (Trial-Level)', ...
    'FontSize', 18, 'FontWeight', 'bold');
plot_raincloud_trial(gaze_BCEA_bl, gaze_Condition, colors, condLabels, ...
    'BCEA [%]');
report_trend_trial('H3 BCEA (%)', gaze_BCEA_bl, gaze_Condition, ...
    gaze_Subject, contrast_vals);
exportgraphics(fig3, fullfile(fig_dir, 'GCP_H3_bcea_trials.png'), 'Resolution', 600);

%% FIGURE 4: H4 + H5 Gamma frequency & power
fig4 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H4] Gamma Peak Frequency & [H5] Peak Power (Trial-Level)', ...
    'FontSize', 18, 'FontWeight', 'bold');

subplot(1, 2, 1);
plot_raincloud_trial(ged_PeakFreq, ged_Condition, colors, condLabels, ...
    'GED Peak Frequency [Hz]');
title('[H4] Peak Frequency', 'FontSize', 14, 'FontWeight', 'bold');

subplot(1, 2, 2);
plot_raincloud_trial(ged_PeakPower, ged_Condition, colors, condLabels, ...
    'GED Peak Power [dB]');
title('[H5] Peak Power', 'FontSize', 14, 'FontWeight', 'bold');

report_trend_trial('H4 GED Peak Freq', ged_PeakFreq, ged_Condition, ...
    ged_Subject, contrast_vals);
report_trend_trial('H5 GED Peak Power', ged_PeakPower, ged_Condition, ...
    ged_Subject, contrast_vals);

exportgraphics(fig4, fullfile(fig_dir, 'GCP_H4H5_gamma_trials.png'), 'Resolution', 600);

%% FIGURE 5: H6 Gamma vs oculomotor (matched trial means + trial-level scatter)
fig5 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('[H6] Gamma vs Oculomotor Dynamics (subject means from matched trials)', ...
    'FontSize', 18, 'FontWeight', 'bold');

subplot(2, 3, 1); hold on;
plot_scatter_by_cond(subj_gedFreq, subj_MSRate_bl, colors, condLabels, nSubj);
xlabel('GED Peak Frequency [Hz]'); ylabel('MS Rate [%]');
title('\gamma Freq vs MS Rate'); set(gca, 'FontSize', fontSize - 3);

subplot(2, 3, 2); hold on;
plot_scatter_by_cond(subj_gedFreq, subj_Vel2D_bl, colors, condLabels, nSubj);
xlabel('GED Peak Frequency [Hz]'); ylabel('Eye Velocity [%]');
title('\gamma Freq vs Velocity'); set(gca, 'FontSize', fontSize - 3);

subplot(2, 3, 3); hold on;
plot_scatter_by_cond(subj_gedFreq, subj_BCEA_bl, colors, condLabels, nSubj);
xlabel('GED Peak Frequency [Hz]'); ylabel('BCEA [%]');
title('\gamma Freq vs BCEA'); set(gca, 'FontSize', fontSize - 3);

subplot(2, 3, 4); hold on;
plot_scatter_by_cond(subj_gedPower, subj_MSRate_bl, colors, condLabels, nSubj);
xlabel('GED Peak Power [dB]'); ylabel('MS Rate [%]');
title('\gamma Power vs MS Rate'); set(gca, 'FontSize', fontSize - 3);

subplot(2, 3, 5); hold on;
plot_scatter_by_cond(subj_gedPower, subj_Vel2D_bl, colors, condLabels, nSubj);
xlabel('GED Peak Power [dB]'); ylabel('Eye Velocity [%]');
title('\gamma Power vs Velocity'); set(gca, 'FontSize', fontSize - 3);

subplot(2, 3, 6); hold on;
plot_scatter_by_cond(subj_gedPower, subj_BCEA_bl, colors, condLabels, nSubj);
xlabel('GED Peak Power [dB]'); ylabel('BCEA [%]');
title('\gamma Power vs BCEA'); set(gca, 'FontSize', fontSize - 3);

exportgraphics(fig5, fullfile(fig_dir, 'GCP_H6_gamma_oculomotor_trials.png'), 'Resolution', 600);

%% FIGURE 6: Summary dashboard
fig6 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Hypothesis Summary (Trial-Level Distributions)', ...
    'FontSize', 16, 'FontWeight', 'bold');

summary_data   = {gaze_MSRate_bl, gaze_Vel2D_bl, gaze_BCEA_bl, ...
                  ged_PeakFreq, ged_PeakPower};
summary_cond   = {gaze_Condition, gaze_Condition, gaze_Condition, ...
                  ged_Condition, ged_Condition};
summary_names  = {'[H1] MS Rate [%]', '[H2] Fix. Velocity [%]', ...
                  '[H3] BCEA [%]', '[H4] \gamma Peak Freq', '[H5] \gamma Peak Power'};
summary_expect = {'decrease', 'increase', 'increase', 'increase', 'increase'};

for mi = 1:5
    subplot(2, 3, mi); hold on;
    dat = summary_data{mi};
    cnd = summary_cond{mi};

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

% Expected direction checklist from LMEs
subplot(2, 3, 6); axis off;
txt = {'\bf{Contrast LME direction check:}', ...
    sprintf('H1 MSRate:  beta=%.3g, p=%.3g (expect <0)', lme_rows.BetaContrast(1), lme_rows.PContrast(1)), ...
    sprintf('H2 Vel2D:   beta=%.3g, p=%.3g (expect >0)', lme_rows.BetaContrast(2), lme_rows.PContrast(2)), ...
    sprintf('H3 BCEA:    beta=%.3g, p=%.3g (expect >0)', lme_rows.BetaContrast(3), lme_rows.PContrast(3)), ...
    sprintf('H4 Freq:    beta=%.3g, p=%.3g (expect >0)', lme_rows.BetaContrast(4), lme_rows.PContrast(4)), ...
    sprintf('H5 Power:   beta=%.3g, p=%.3g (expect >0)', lme_rows.BetaContrast(5), lme_rows.PContrast(5))};
text(0.05, 0.55, txt, 'FontSize', 12, 'VerticalAlignment', 'middle', ...
    'Interpreter', 'none');

exportgraphics(fig6, fullfile(fig_dir, 'GCP_hypotheses_summary_trials.png'), 'Resolution', 600);

%% Done
fprintf('\n[STATS HYP TRIALS] Trial-Level Hypothesis Testing Complete\n');
fprintf('[STATS HYP TRIALS] Figures saved to:\n  %s\n', fig_dir);
fprintf('[STATS HYP TRIALS] LME table saved to:\n  %s\n', fullfile(data_dir, 'GCP_hypotheses_trials_lme.csv'));
fprintf('[STATS HYP TRIALS] Trial rows: %d | Subjects: %d\n', height(T), nSubj);

%% HELPER FUNCTIONS

function T = restrict_to_subjects(T, subjects)
id_var = pick_var(T, {'ID'});
if isempty(id_var)
    return
end
ids = T.(id_var);
if iscell(ids) || isstring(ids) || ischar(ids)
    id_str = string(ids);
else
    id_str = string(ids);
end
keep = ismember(id_str, string(subjects(:)));
T = T(keep, :);
end

function T = standardize_condition_codes(T)
cond_var = pick_var(T, {'Condition'});
if isempty(cond_var)
    return
end
c = T.(cond_var);
if iscell(c) || isstring(c) || ischar(c) || iscategorical(c)
    lab = string(c);
    out = nan(height(T), 1);
    out(lab == "25%" | lab == "25") = 1;
    out(lab == "50%" | lab == "50") = 2;
    out(lab == "75%" | lab == "75") = 3;
    out(lab == "100%" | lab == "100") = 4;
    T.(cond_var) = out;
else
    T.(cond_var) = double(c);
end
end

function name = pick_var(T, candidates)
name = '';
for i = 1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        name = candidates{i};
        return
    end
end
end

function T = build_trial_table_from_sources(subjects, paths, data_dir)
% Fallback when merged table is absent.
nSubj = numel(subjects);
ID = [];
Condition = [];
Trial = [];
MSRate_bl = [];
Vel2D_bl = [];
BCEA_bl = [];
GammaFrequency = [];
GammaPower = [];

for subj = 1:nSubj
    gazepath = fullfile(paths.features, subjects{subj}, 'gaze');
    trialfile = fullfile(gazepath, 'gaze_matrix_trial.mat');
    if ~exist(trialfile, 'file')
        fprintf('  %s: gaze_matrix_trial.mat not found, skipping gaze.\n', subjects{subj});
        continue
    end
    dat = load(trialfile);
    sid = str2double(subjects{subj});
    cond_names = {'subj_data_gaze_trial_c25', 'subj_data_gaze_trial_c50', ...
                  'subj_data_gaze_trial_c75', 'subj_data_gaze_trial_c100'};
    for ci = 1:4
        if ~isfield(dat, cond_names{ci}), continue; end
        d = dat.(cond_names{ci});
        nTrl = numel(d.Trial);
        ID = [ID; repmat(sid, nTrl, 1)]; %#ok<AGROW>
        Condition = [Condition; repmat(ci, nTrl, 1)]; %#ok<AGROW>
        Trial = [Trial; d.Trial(:)]; %#ok<AGROW>
        MSRate_bl = [MSRate_bl; pick_field_vec(d, {'MSRate_bl', 'dBMSRate'}, nTrl)]; %#ok<AGROW>
        Vel2D_bl  = [Vel2D_bl;  pick_field_vec(d, {'Vel2D_bl', 'dBVel2D'}, nTrl)]; %#ok<AGROW>
        BCEA_bl   = [BCEA_bl;   pick_field_vec(d, {'BCEA_bl', 'dBBCEA'}, nTrl)]; %#ok<AGROW>
        GammaFrequency = [GammaFrequency; nan(nTrl, 1)]; %#ok<AGROW>
        GammaPower = [GammaPower; nan(nTrl, 1)]; %#ok<AGROW>
    end
end

ged_path = fullfile(data_dir, 'GCP_eeg_GED.mat');
if isfile(ged_path)
    ged = load(ged_path, ...
        'trials_peaks', 'trials_powratio_fullscan', ...
        'trials_outlier_mask_power_full', 'scan_freqs', 'subjects');
    scan_freqs = ged.scan_freqs(:);
    for subj = 1:nSubj
        gi = find(strcmp(ged.subjects, subjects{subj}));
        if isempty(gi), continue; end
        sid = str2double(subjects{subj});
        for cond = 1:4
            pf = ged.trials_peaks{cond, gi};
            if isempty(pf), continue; end
            pf = pf(:);
            nTrl = numel(pf);
            pr = ged.trials_powratio_fullscan{cond, gi};
            pp = reconstruct_trial_peak_power(pf, pr, scan_freqs, 5);
            if isfield(ged, 'trials_outlier_mask_power_full')
                mask = ged.trials_outlier_mask_power_full{cond, gi};
                if ~isempty(mask) && numel(mask) == nTrl
                    pp(logical(mask(:))) = NaN;
                end
            end
            for t = 1:nTrl
                idx = ID == sid & Condition == cond & Trial == t;
                if any(idx)
                    GammaFrequency(idx) = pf(t);
                    GammaPower(idx) = pp(t);
                else
                    ID = [ID; sid]; %#ok<AGROW>
                    Condition = [Condition; cond]; %#ok<AGROW>
                    Trial = [Trial; t]; %#ok<AGROW>
                    MSRate_bl = [MSRate_bl; NaN]; %#ok<AGROW>
                    Vel2D_bl = [Vel2D_bl; NaN]; %#ok<AGROW>
                    BCEA_bl = [BCEA_bl; NaN]; %#ok<AGROW>
                    GammaFrequency = [GammaFrequency; pf(t)]; %#ok<AGROW>
                    GammaPower = [GammaPower; pp(t)]; %#ok<AGROW>
                end
            end
        end
    end
end

T = table(ID, Condition, Trial, MSRate_bl, Vel2D_bl, BCEA_bl, GammaFrequency, GammaPower);
end

function v = pick_field_vec(d, names, nTrl)
v = nan(nTrl, 1);
for i = 1:numel(names)
    if isfield(d, names{i})
        x = d.(names{i});
        v = x(:);
        if numel(v) ~= nTrl
            v = nan(nTrl, 1);
        end
        return
    end
end
end

function plot_raincloud_trial(data, cond_vec, colors, condLabels, y_label)
    hold on;
    for c = 1:4
        vals = data(cond_vec == c);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end

        if numel(vals) >= 3
            [f_dens, xi] = ksdensity(vals);
            f_dens = f_dens / max(f_dens) * 0.3;
            patch(c - f_dens - 0.05, xi, colors(c,:), ...
                'FaceAlpha', 0.3, 'EdgeColor', colors(c,:), 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end

        xJit = c + 0.15 + (rand(size(vals)) - 0.5) * 0.15;
        scatter(xJit, vals, 15, colors(c,:), 'filled', ...
            'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');

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
    valid_r = ~isnan(r_per_subj);
    r_mean = mean(r_per_subj(valid_r));
    if sum(valid_r) >= 2
        [~, p] = ttest(r_per_subj(valid_r));
    else
        p = NaN;
    end
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

function pp = reconstruct_trial_peak_power(pf, pr, scan_freqs, halfwidth_hz)
nTrl = numel(pf);
pp = nan(nTrl, 1);
scan_freqs = scan_freqs(:);
if isempty(pr) || size(pr, 1) ~= nTrl
    return
end
for t = 1:nTrl
    if ~isfinite(pf(t))
        continue
    end
    band = abs(scan_freqs - pf(t)) <= halfwidth_hz;
    if ~any(band)
        continue
    end
    pp(t) = mean(pr(t, band), 'omitnan');
end
end

function [data_out, n_removed] = reject_outliers_iqr_by_cond(data, cond_vec)
    data_out  = data;
    n_removed = 0;
    for c = 1:4
        idx = cond_vec == c;
        vals = data(idx);
        valid = ~isnan(vals);
        if sum(valid) < 5, continue; end
        q1 = prctile(vals(valid), 25);
        q3 = prctile(vals(valid), 75);
        iqr_val = q3 - q1;
        if iqr_val == 0, continue; end
        lo = q1 - 1.5 * iqr_val;
        hi = q3 + 1.5 * iqr_val;
        outlier = valid & (vals < lo | vals > hi);
        n_removed = n_removed + sum(outlier);
        vals(outlier) = NaN;
        data_out(idx) = vals;
    end
end

function out = empty_lme_table()
out = table('Size', [0, 7], ...
    'VariableTypes', {'string','double','double','double','double','double','double'}, ...
    'VariableNames', {'Model','N','NumSubjects','BetaContrast','PContrast','BetaGaze','PGaze'});
end

function row = fit_contrast_lme(model_name, y, cond_vec, subj_vec, contrast_vals)
ok = isfinite(y) & isfinite(cond_vec) & isfinite(subj_vec) & cond_vec >= 1 & cond_vec <= 4;
y = y(ok); cond_vec = cond_vec(ok); subj_vec = subj_vec(ok);
if numel(y) < 10
    row = make_lme_row(model_name, numel(y), numel(unique(subj_vec)), NaN, NaN, NaN, NaN);
    return
end
tbl = table();
tbl.Y = y(:);
tbl.Contrast = contrast_vals(cond_vec(:))';
tbl.ID = categorical(subj_vec(:));
try
    mdl = fitlme(tbl, 'Y ~ Contrast + (1|ID)', 'FitMethod', 'REML');
    [b, p] = coef_beta_p(mdl, 'Contrast');
    row = make_lme_row(model_name, height(tbl), numel(unique(subj_vec)), b, p, NaN, NaN);
catch
    row = make_lme_row(model_name, height(tbl), numel(unique(subj_vec)), NaN, NaN, NaN, NaN);
end
end

function row = fit_gamma_gaze_lme(model_name, y_gamma, x_gaze, cond_vec, subj_vec, contrast_vals)
ok = isfinite(y_gamma) & isfinite(x_gaze) & isfinite(cond_vec) & isfinite(subj_vec) ...
    & cond_vec >= 1 & cond_vec <= 4;
y_gamma = y_gamma(ok); x_gaze = x_gaze(ok);
cond_vec = cond_vec(ok); subj_vec = subj_vec(ok);
if numel(y_gamma) < 10
    row = make_lme_row(model_name, numel(y_gamma), numel(unique(subj_vec)), NaN, NaN, NaN, NaN);
    return
end
% z-score gaze within sample so beta is per SD
mu = mean(x_gaze); sd = std(x_gaze);
if sd > 0
    x_gaze = (x_gaze - mu) / sd;
end
tbl = table();
tbl.Y = y_gamma(:);
tbl.Gaze = x_gaze(:);
tbl.Contrast = contrast_vals(cond_vec(:))';
tbl.ID = categorical(subj_vec(:));
try
    mdl = fitlme(tbl, 'Y ~ Gaze + Contrast + (1|ID)', 'FitMethod', 'REML');
    [bC, pC] = coef_beta_p(mdl, 'Contrast');
    [bG, pG] = coef_beta_p(mdl, 'Gaze');
    row = make_lme_row(model_name, height(tbl), numel(unique(subj_vec)), bC, pC, bG, pG);
catch
    row = make_lme_row(model_name, height(tbl), numel(unique(subj_vec)), NaN, NaN, NaN, NaN);
end
end

function [beta, p] = coef_beta_p(mdl, name)
cn = string(mdl.Coefficients.Name);
ix = find(cn == string(name), 1);
if isempty(ix)
    beta = NaN; p = NaN;
    return
end
beta = mdl.Coefficients.Estimate(ix);
p = mdl.Coefficients.pValue(ix);
end

function row = make_lme_row(model_name, n, n_subj, bC, pC, bG, pG)
row = table(string(model_name), n, n_subj, bC, pC, bG, pG, ...
    'VariableNames', {'Model','N','NumSubjects','BetaContrast','PContrast','BetaGaze','PGaze'});
end
