%% GCP GED Split Gamma (Low/High) with Single-Peak Handling
%
% Purpose:
%   Build a participant-specific low/high gamma split from the trough between
%   the two strongest peaks in the participant-averaged detrended GED spectrum.
%   Then apply that split to all trials while keeping single-peak trials.
%
% Input:
%   GCP_eeg_GED_trials.mat (must contain all_trial_powratio, scan_freqs, subjects)
%
% Output (all contained in splitGammaLowHigh folders):
%   - Subject split table (.csv + .mat)
%   - Trial-level long table with low/high metrics and peak-presence flags
%   - Subject x condition summary table
%   - QC / summary figures

clear; close all; clc

%% Setup
startup
[subjects_setup, ~, colors, ~] = setup('GCP');

condLabels_default = {'25%', '50%', '75%', '100%'};
contrast_vals = [25, 50, 75, 100];

% Parameters for split QC
params.poly_order = 2;
params.smooth_win = 5;
params.min_peak_distance = 5;       % Hz
params.min_sep_hz = 8;              % minimum distance between top two peaks
params.min_prom_rel = 0.10;         % each selected peak >= 10% of max positive peak
params.min_trough_drop_rel = 0.05;  % trough must be >=5% below smaller selected peak
params.fallback_split_hz = 50;

% IO paths (new isolated folders)
if ispc
    data_root = 'W:\Students\Arne\GCP\data\features';
    fig_root  = 'W:\Students\Arne\GCP\figures\eeg\ged';
else
    data_root = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features';
    fig_root  = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged';
end

data_dir = fullfile(data_root, 'splitGammaLowHigh');
fig_dir  = fullfile(fig_root,  'splitGammaLowHigh');
if ~exist(data_dir, 'dir'), mkdir(data_dir); end
if ~exist(fig_dir,  'dir'), mkdir(fig_dir); end

%% Load trial-level GED source data
in_path = fullfile(data_root, 'GCP_eeg_GED_trials.mat');
fprintf('Loading GED trial data:\n  %s\n', in_path);
S = load(in_path);

required_vars = {'all_trial_powratio', 'scan_freqs'};
for vi = 1:numel(required_vars)
    if ~isfield(S, required_vars{vi})
        error('Missing required variable "%s" in %s.', required_vars{vi}, in_path);
    end
end

all_trial_powratio = S.all_trial_powratio;
scan_freqs = S.scan_freqs(:)';

if isfield(S, 'subjects')
    subjects = S.subjects;
else
    subjects = subjects_setup;
end
nSubj = numel(subjects);

if isfield(S, 'condLabels')
    condLabels = S.condLabels;
else
    condLabels = condLabels_default;
end

nCond = 4;

%% Preallocate subject split outputs
split_hz = nan(nSubj, 1);
split_method = strings(nSubj, 1);
split_qc_pass = false(nSubj, 1);
split_pk1_hz = nan(nSubj, 1);
split_pk2_hz = nan(nSubj, 1);
split_pk1_amp = nan(nSubj, 1);
split_pk2_amp = nan(nSubj, 1);
split_trough_amp = nan(nSubj, 1);
split_sep_hz = nan(nSubj, 1);

% Keep spectra for plotting
subject_mean_dt = nan(nSubj, numel(scan_freqs));

%% Trial-level long-format rows
row_Subject = [];
row_Condition = [];
row_Contrast = [];
row_Trial = [];
row_SplitHz = [];
row_SplitMethod = strings(0,1);
row_SplitQCPass = [];
row_nPeaks = [];
row_PeakClass = strings(0,1);    % none/single_low/single_high/dual
row_HasLowPeak = [];
row_HasHighPeak = [];
row_DominantPeakFreq = [];
row_DominantPeakAmp = [];
row_DominantBand = strings(0,1); % none/low/high
row_LowPeakFreq = [];
row_LowPeakAmp = [];
row_HighPeakFreq = [];
row_HighPeakAmp = [];
row_LowPower = [];
row_HighPower = [];
row_BroadbandPower = [];

fprintf('Computing participant-specific split + trial metrics...\n');

for subj = 1:nSubj
    % ---------- Build participant mean detrended spectrum ----------
    all_dt = [];
    for cond = 1:nCond
        pr_mat = all_trial_powratio{cond, subj};
        if isempty(pr_mat), continue; end
        nTrl = size(pr_mat, 1);
        dt = nan(nTrl, numel(scan_freqs));
        for trl = 1:nTrl
            pr = pr_mat(trl, :);
            if all(isnan(pr)), continue; end
            p = polyfit(scan_freqs, pr, params.poly_order);
            dt(trl, :) = pr - polyval(p, scan_freqs);
        end
        all_dt = [all_dt; dt]; %#ok<AGROW>
    end

    if isempty(all_dt) || all(all(isnan(all_dt)))
        split_hz(subj) = params.fallback_split_hz;
        split_method(subj) = "fallback50_noData";
        split_qc_pass(subj) = false;
    else
        mu_dt = nanmean(all_dt, 1);
        mu_dt_s = movmean(mu_dt, params.smooth_win);
        subject_mean_dt(subj, :) = mu_dt_s;

        [pks, locs] = findpeaks(mu_dt_s, scan_freqs, ...
            'MinPeakDistance', params.min_peak_distance);
        pos = pks > 0;
        pks_pos = pks(pos);
        locs_pos = locs(pos);

        do_fallback = true;
        if numel(pks_pos) >= 2
            [p_sorted, idx_sorted] = sort(pks_pos, 'descend');
            pk1_amp = p_sorted(1);
            pk2_amp = p_sorted(2);
            pk1_hz = locs_pos(idx_sorted(1));
            pk2_hz = locs_pos(idx_sorted(2));

            lo_pk = min(pk1_hz, pk2_hz);
            hi_pk = max(pk1_hz, pk2_hz);
            lo_idx = find(scan_freqs == lo_pk, 1, 'first');
            hi_idx = find(scan_freqs == hi_pk, 1, 'first');

            sep_hz = hi_pk - lo_pk;
            max_pos_amp = max(pks_pos);
            amp_ok = (pk1_amp >= params.min_prom_rel * max_pos_amp) && ...
                     (pk2_amp >= params.min_prom_rel * max_pos_amp);
            sep_ok = sep_hz >= params.min_sep_hz;

            if ~isempty(lo_idx) && ~isempty(hi_idx) && hi_idx > lo_idx
                seg = mu_dt_s(lo_idx:hi_idx);
                [trough_amp, rel_min_idx] = min(seg);
                trough_idx = lo_idx + rel_min_idx - 1;
                trough_hz = scan_freqs(trough_idx);

                min_pk_amp = min(pk1_amp, pk2_amp);
                trough_drop = min_pk_amp - trough_amp;
                trough_ok = trough_drop >= params.min_trough_drop_rel * max_pos_amp;

                if amp_ok && sep_ok && trough_ok
                    split_hz(subj) = trough_hz;
                    split_method(subj) = "trough";
                    split_qc_pass(subj) = true;
                    do_fallback = false;
                end

                split_pk1_hz(subj) = pk1_hz;
                split_pk2_hz(subj) = pk2_hz;
                split_pk1_amp(subj) = pk1_amp;
                split_pk2_amp(subj) = pk2_amp;
                split_trough_amp(subj) = trough_amp;
                split_sep_hz(subj) = sep_hz;
            end
        end

        if do_fallback
            split_hz(subj) = params.fallback_split_hz;
            if numel(pks_pos) < 2
                split_method(subj) = "fallback50_lt2Peaks";
            else
                split_method(subj) = "fallback50_qcFail";
            end
            split_qc_pass(subj) = false;
        end
    end

    % ---------- Per-subject QC figure ----------
    fig_sub = figure('Position', [0 0 1512 982], 'Color', 'w');
    tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile; hold on;
    if ~all(isnan(subject_mean_dt(subj,:)))
        plot(scan_freqs, subject_mean_dt(subj,:), 'k-', 'LineWidth', 2);
        yline(0, 'k:', 'LineWidth', 1);
        xline(split_hz(subj), '--', 'LineWidth', 2, 'Color', [0.85 0.2 0.2]);
        if ~isnan(split_pk1_hz(subj)), xline(split_pk1_hz(subj), ':', 'Color', [0.2 0.2 0.7]); end
        if ~isnan(split_pk2_hz(subj)), xline(split_pk2_hz(subj), ':', 'Color', [0.2 0.2 0.7]); end
    end
    xlim([scan_freqs(1), scan_freqs(end)]);
    xlabel('Frequency [Hz]'); ylabel('\Delta Power ratio');
    title('Participant mean detrended spectrum');
    grid on; box on;

    nexttile; axis off;
    txt = {
        sprintf('Subject: %s', subjects{subj})
        sprintf('split_hz: %.1f', split_hz(subj))
        sprintf('method: %s', split_method(subj))
        sprintf('qc_pass: %d', split_qc_pass(subj))
        sprintf('pk1_hz/pk2_hz: %.1f / %.1f', split_pk1_hz(subj), split_pk2_hz(subj))
        sprintf('pk1_amp/pk2_amp: %.3f / %.3f', split_pk1_amp(subj), split_pk2_amp(subj))
        sprintf('sep_hz: %.1f', split_sep_hz(subj))
    };
    text(0.05, 0.9, txt, 'FontSize', 12, 'VerticalAlignment', 'top');

    nexttile([1 2]); hold on;
    edges = scan_freqs(1):2:scan_freqs(end);
    centers = edges(1:end-1) + diff(edges)/2;
    hmat = zeros(nCond, numel(centers));
    for cond = 1:nCond
        pr_mat = all_trial_powratio{cond, subj};
        if isempty(pr_mat), continue; end
        cond_peaks = [];
        for trl = 1:size(pr_mat, 1)
            pr = pr_mat(trl, :);
            if all(isnan(pr)), continue; end
            p = polyfit(scan_freqs, pr, params.poly_order);
            dt = movmean(pr - polyval(p, scan_freqs), params.smooth_win);
            [pks_t, locs_t] = findpeaks(dt, scan_freqs, 'MinPeakDistance', params.min_peak_distance);
            cond_peaks = [cond_peaks, locs_t(pks_t > 0)]; %#ok<AGROW>
        end
        hmat(cond, :) = histcounts(cond_peaks, edges);
    end
    bh = bar(centers, hmat', 'stacked', 'EdgeColor', 'none', 'BarWidth', 1);
    for cond = 1:nCond
        bh(cond).FaceColor = colors(cond, :);
    end
    xline(split_hz(subj), '--', 'Color', [0 0 0], 'LineWidth', 2);
    xlim([scan_freqs(1), scan_freqs(end)]);
    xlabel('Peak frequency [Hz]'); ylabel('Count');
    title('All detected positive peaks (stacked by condition)');
    legend(condLabels, 'Location', 'best');
    grid on; box on;

    sgtitle(sprintf('SplitGammaLowHigh QC — Subject %s', subjects{subj}), ...
        'FontSize', 16, 'FontWeight', 'bold');
    saveas(fig_sub, fullfile(fig_dir, sprintf('splitGammaLowHigh_subj_%s.png', subjects{subj})));
    close(fig_sub);

    % ---------- Trial-level metrics using participant split ----------
    lo_mask_subj = scan_freqs < split_hz(subj);
    hi_mask_subj = scan_freqs >= split_hz(subj);

    for cond = 1:nCond
        pr_mat = all_trial_powratio{cond, subj};
        if isempty(pr_mat), continue; end
        nTrl = size(pr_mat, 1);

        for trl = 1:nTrl
            pr = pr_mat(trl, :);

            sid_num = str2double(subjects{subj});
            if isnan(sid_num), sid_num = subj; end

            row_Subject(end+1,1) = sid_num; %#ok<SAGROW>
            row_Condition(end+1,1) = cond; %#ok<SAGROW>
            row_Contrast(end+1,1) = contrast_vals(cond); %#ok<SAGROW>
            row_Trial(end+1,1) = trl; %#ok<SAGROW>
            row_SplitHz(end+1,1) = split_hz(subj); %#ok<SAGROW>
            row_SplitMethod(end+1,1) = split_method(subj); %#ok<SAGROW>
            row_SplitQCPass(end+1,1) = split_qc_pass(subj); %#ok<SAGROW>

            if all(isnan(pr))
                row_nPeaks(end+1,1) = 0; %#ok<SAGROW>
                row_PeakClass(end+1,1) = "none"; %#ok<SAGROW>
                row_HasLowPeak(end+1,1) = false; %#ok<SAGROW>
                row_HasHighPeak(end+1,1) = false; %#ok<SAGROW>
                row_DominantPeakFreq(end+1,1) = nan; %#ok<SAGROW>
                row_DominantPeakAmp(end+1,1) = nan; %#ok<SAGROW>
                row_DominantBand(end+1,1) = "none"; %#ok<SAGROW>
                row_LowPeakFreq(end+1,1) = nan; %#ok<SAGROW>
                row_LowPeakAmp(end+1,1) = nan; %#ok<SAGROW>
                row_HighPeakFreq(end+1,1) = nan; %#ok<SAGROW>
                row_HighPeakAmp(end+1,1) = nan; %#ok<SAGROW>
                row_LowPower(end+1,1) = nan; %#ok<SAGROW>
                row_HighPower(end+1,1) = nan; %#ok<SAGROW>
                row_BroadbandPower(end+1,1) = nan; %#ok<SAGROW>
                continue;
            end

            % Power metrics are always computed
            low_pow = nanmean(pr(lo_mask_subj));
            high_pow = nanmean(pr(hi_mask_subj));
            bb_pow = nanmean(pr);

            % Detrended peak detection
            p = polyfit(scan_freqs, pr, params.poly_order);
            dt = movmean(pr - polyval(p, scan_freqs), params.smooth_win);
            [pks_t, locs_t] = findpeaks(dt, scan_freqs, 'MinPeakDistance', params.min_peak_distance);
            pos = pks_t > 0;
            pks_pos = pks_t(pos);
            locs_pos = locs_t(pos);

            nPeaks = numel(pks_pos);

            has_lo = any(locs_pos < split_hz(subj));
            has_hi = any(locs_pos >= split_hz(subj));

            low_pf = nan; low_pa = nan;
            high_pf = nan; high_pa = nan;
            dom_pf = nan; dom_pa = nan; dom_band = "none";

            if has_lo
                [low_pa, bi] = max(pks_pos(locs_pos < split_hz(subj)));
                tmp = locs_pos(locs_pos < split_hz(subj));
                low_pf = tmp(bi);
            end
            if has_hi
                [high_pa, bi] = max(pks_pos(locs_pos >= split_hz(subj)));
                tmp = locs_pos(locs_pos >= split_hz(subj));
                high_pf = tmp(bi);
            end

            if nPeaks > 0
                [dom_pa, bi] = max(pks_pos);
                dom_pf = locs_pos(bi);
                if dom_pf < split_hz(subj)
                    dom_band = "low";
                else
                    dom_band = "high";
                end
            end

            if nPeaks == 0
                peak_class = "none";
            elseif has_lo && has_hi
                peak_class = "dual";
            elseif has_lo
                peak_class = "single_low";
            else
                peak_class = "single_high";
            end

            row_nPeaks(end+1,1) = nPeaks; %#ok<SAGROW>
            row_PeakClass(end+1,1) = peak_class; %#ok<SAGROW>
            row_HasLowPeak(end+1,1) = has_lo; %#ok<SAGROW>
            row_HasHighPeak(end+1,1) = has_hi; %#ok<SAGROW>
            row_DominantPeakFreq(end+1,1) = dom_pf; %#ok<SAGROW>
            row_DominantPeakAmp(end+1,1) = dom_pa; %#ok<SAGROW>
            row_DominantBand(end+1,1) = dom_band; %#ok<SAGROW>
            row_LowPeakFreq(end+1,1) = low_pf; %#ok<SAGROW>
            row_LowPeakAmp(end+1,1) = low_pa; %#ok<SAGROW>
            row_HighPeakFreq(end+1,1) = high_pf; %#ok<SAGROW>
            row_HighPeakAmp(end+1,1) = high_pa; %#ok<SAGROW>
            row_LowPower(end+1,1) = low_pow; %#ok<SAGROW>
            row_HighPower(end+1,1) = high_pow; %#ok<SAGROW>
            row_BroadbandPower(end+1,1) = bb_pow; %#ok<SAGROW>
        end
    end

    clc
    fprintf('Subject %s (%d/%d): split %.1f Hz (%s)\n', ...
        subjects{subj}, subj, nSubj, split_hz(subj), split_method(subj));
end

%% Build and save subject split table
SubjectNumeric = nan(nSubj, 1);
for s = 1:nSubj
    v = str2double(subjects{s});
    if isnan(v), v = s; end
    SubjectNumeric(s) = v;
end

T_split = table( ...
    SubjectNumeric, split_hz, split_method, split_qc_pass, ...
    split_pk1_hz, split_pk2_hz, split_pk1_amp, split_pk2_amp, ...
    split_trough_amp, split_sep_hz, ...
    'VariableNames', {'Subject','SplitHz','Method','QCPass', ...
    'Peak1Hz','Peak2Hz','Peak1Amp','Peak2Amp','TroughAmp','PeakSeparationHz'});

split_csv = fullfile(data_dir, 'GCP_eeg_GED_splitGammaLowHigh_subjectSplits.csv');
split_mat = fullfile(data_dir, 'GCP_eeg_GED_splitGammaLowHigh_subjectSplits.mat');
writetable(T_split, split_csv);
save(split_mat, 'T_split', 'params', 'scan_freqs', 'subjects', 'condLabels');

%% Build and save trial-level long table
T_trials = table( ...
    row_Subject, row_Condition, row_Contrast, row_Trial, ...
    row_SplitHz, row_SplitMethod, row_SplitQCPass, ...
    row_nPeaks, row_PeakClass, row_HasLowPeak, row_HasHighPeak, ...
    row_DominantPeakFreq, row_DominantPeakAmp, row_DominantBand, ...
    row_LowPeakFreq, row_LowPeakAmp, row_HighPeakFreq, row_HighPeakAmp, ...
    row_LowPower, row_HighPower, row_BroadbandPower, ...
    'VariableNames', {'Subject','Condition','Contrast','Trial', ...
    'SplitHz','SplitMethod','SplitQCPass', ...
    'NumPeaks','PeakClass','HasLowPeak','HasHighPeak', ...
    'DominantPeakFreq','DominantPeakAmp','DominantBand', ...
    'LowPeakFreq','LowPeakAmp','HighPeakFreq','HighPeakAmp', ...
    'LowGammaPower','HighGammaPower','BroadbandPower'});

trials_csv = fullfile(data_dir, 'GCP_eeg_GED_splitGammaLowHigh_trials.csv');
trials_mat = fullfile(data_dir, 'GCP_eeg_GED_splitGammaLowHigh_trials.mat');
writetable(T_trials, trials_csv);
save(trials_mat, 'T_trials', 'params', 'scan_freqs', 'subjects', 'condLabels');

%% Build subject x condition summary
nRows = nSubj * nCond;
sum_Subject = nan(nRows, 1);
sum_Condition = nan(nRows, 1);
sum_Contrast = nan(nRows, 1);
sum_SplitHz = nan(nRows, 1);
sum_QCPass = false(nRows, 1);
sum_NTrials = zeros(nRows, 1);
sum_DualRate = nan(nRows, 1);
sum_SingleLowRate = nan(nRows, 1);
sum_SingleHighRate = nan(nRows, 1);
sum_NoPeakRate = nan(nRows, 1);
sum_HasLowRate = nan(nRows, 1);
sum_HasHighRate = nan(nRows, 1);
sum_MeanLowPeakFreq = nan(nRows, 1);
sum_MeanHighPeakFreq = nan(nRows, 1);
sum_MeanLowPower = nan(nRows, 1);
sum_MeanHighPower = nan(nRows, 1);
sum_MeanBroadbandPower = nan(nRows, 1);

ri = 0;
for subj = 1:nSubj
    sid = SubjectNumeric(subj);
    for cond = 1:nCond
        ri = ri + 1;
        idx = (T_trials.Subject == sid) & (T_trials.Condition == cond);
        tt = T_trials(idx, :);
        nT = height(tt);
        if nT == 0, continue; end

        sum_Subject(ri) = sid;
        sum_Condition(ri) = cond;
        sum_Contrast(ri) = contrast_vals(cond);
        sum_SplitHz(ri) = split_hz(subj);
        sum_QCPass(ri) = split_qc_pass(subj);
        sum_NTrials(ri) = nT;

        sum_DualRate(ri) = mean(tt.PeakClass == "dual");
        sum_SingleLowRate(ri) = mean(tt.PeakClass == "single_low");
        sum_SingleHighRate(ri) = mean(tt.PeakClass == "single_high");
        sum_NoPeakRate(ri) = mean(tt.PeakClass == "none");
        sum_HasLowRate(ri) = mean(tt.HasLowPeak);
        sum_HasHighRate(ri) = mean(tt.HasHighPeak);

        sum_MeanLowPeakFreq(ri) = nanmean(tt.LowPeakFreq);
        sum_MeanHighPeakFreq(ri) = nanmean(tt.HighPeakFreq);
        sum_MeanLowPower(ri) = nanmean(tt.LowGammaPower);
        sum_MeanHighPower(ri) = nanmean(tt.HighGammaPower);
        sum_MeanBroadbandPower(ri) = nanmean(tt.BroadbandPower);
    end
end

T_summary = table( ...
    sum_Subject, sum_Condition, sum_Contrast, sum_SplitHz, sum_QCPass, sum_NTrials, ...
    sum_DualRate, sum_SingleLowRate, sum_SingleHighRate, sum_NoPeakRate, ...
    sum_HasLowRate, sum_HasHighRate, ...
    sum_MeanLowPeakFreq, sum_MeanHighPeakFreq, ...
    sum_MeanLowPower, sum_MeanHighPower, sum_MeanBroadbandPower, ...
    'VariableNames', {'Subject','Condition','Contrast','SplitHz','SplitQCPass','NTrials', ...
    'DualPeakRate','SingleLowRate','SingleHighRate','NoPeakRate', ...
    'HasLowPeakRate','HasHighPeakRate', ...
    'MeanLowPeakFreq','MeanHighPeakFreq', ...
    'MeanLowGammaPower','MeanHighGammaPower','MeanBroadbandPower'});

summary_csv = fullfile(data_dir, 'GCP_eeg_GED_splitGammaLowHigh_summary.csv');
summary_mat = fullfile(data_dir, 'GCP_eeg_GED_splitGammaLowHigh_summary.mat');
writetable(T_summary, summary_csv);
save(summary_mat, 'T_summary');

%% Grand figure 1: split boundaries + fallback rate
fig1 = figure('Position', [0 0 1512 982], 'Color', 'w');
tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
histogram(split_hz, scan_freqs(1):2:scan_freqs(end), ...
    'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'w');
xline(params.fallback_split_hz, '--k', 'LineWidth', 2);
xlabel('Participant split [Hz]'); ylabel('Count');
title('Distribution of participant split boundaries');
grid on; box on;

nexttile; hold on;
is_fallback = ~split_qc_pass;
bar(1, mean(~is_fallback) * 100, 0.5, 'FaceColor', [0.2 0.6 0.2]);
bar(2, mean(is_fallback) * 100, 0.5, 'FaceColor', [0.8 0.3 0.3]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Trough QC pass', 'Fallback 50 Hz'});
ylabel('Participants [%]');
title('Split method usage');
ylim([0 100]); grid on; box on;

nexttile([1 2]); hold on;
for s = 1:nSubj
    if all(isnan(subject_mean_dt(s,:))), continue; end
    c = [0.6 0.6 0.6];
    if split_qc_pass(s), c = [0.1 0.5 0.1]; end
    plot(scan_freqs, subject_mean_dt(s,:), '-', 'Color', [c 0.45], 'LineWidth', 1.5);
    xline(split_hz(s), ':', 'Color', [c 0.35]);
end
yline(0, 'k-', 'LineWidth', 1);
xlim([scan_freqs(1), scan_freqs(end)]);
xlabel('Frequency [Hz]'); ylabel('\Delta Power ratio');
title('Participant mean detrended spectra with individual split lines');
grid on; box on;

sgtitle('SplitGammaLowHigh — Participant Split QC', 'FontSize', 18, 'FontWeight', 'bold');
saveas(fig1, fullfile(fig_dir, 'splitGammaLowHigh_splitQC_overview.png'));
close(fig1);

%% Grand figure 2: peak class composition by condition
fig2 = figure('Position', [0 0 1512 982], 'Color', 'w');
cls = ["none","single_low","single_high","dual"];
cls_labels = {'No peak', 'Single low', 'Single high', 'Dual'};
cls_colors = [0.75 0.75 0.75; 0.2 0.4 0.8; 0.8 0.3 0.3; 0.2 0.65 0.4];

pct = nan(nCond, numel(cls));
for cond = 1:nCond
    idx = T_trials.Condition == cond;
    for ci = 1:numel(cls)
        pct(cond, ci) = mean(T_trials.PeakClass(idx) == cls(ci)) * 100;
    end
end

b = bar(1:nCond, pct, 'stacked', 'EdgeColor', 'none');
for ci = 1:numel(cls)
    b(ci).FaceColor = cls_colors(ci, :);
end
set(gca, 'XTick', 1:nCond, 'XTickLabel', condLabels);
ylabel('Trials [%]');
xlabel('Contrast condition');
ylim([0 100]);
legend(cls_labels, 'Location', 'eastoutside');
title('Trial composition: none vs single-peak vs dual-peak');
grid on; box on;
saveas(fig2, fullfile(fig_dir, 'splitGammaLowHigh_peakClass_byCondition.png'));
close(fig2);

%% Grand figure 3: low/high power CRF and low/high peak detection rates
fig3 = figure('Position', [0 0 1512 982], 'Color', 'w');
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
mu_lo = nan(nCond,1); sem_lo = nan(nCond,1);
mu_hi = nan(nCond,1); sem_hi = nan(nCond,1);
for cond = 1:nCond
    v_lo = T_summary.MeanLowGammaPower(T_summary.Condition == cond);
    v_hi = T_summary.MeanHighGammaPower(T_summary.Condition == cond);
    mu_lo(cond) = nanmean(v_lo); sem_lo(cond) = nanstd(v_lo) / sqrt(sum(~isnan(v_lo)));
    mu_hi(cond) = nanmean(v_hi); sem_hi(cond) = nanstd(v_hi) / sqrt(sum(~isnan(v_hi)));
end
errorbar(contrast_vals, mu_lo, sem_lo, '-o', 'Color', [0.2 0.4 0.85], ...
    'LineWidth', 2.5, 'MarkerFaceColor', [0.2 0.4 0.85]);
errorbar(contrast_vals, mu_hi, sem_hi, '-s', 'Color', [0.85 0.25 0.25], ...
    'LineWidth', 2.5, 'MarkerFaceColor', [0.85 0.25 0.25]);
xlabel('Contrast [%]'); ylabel('Power ratio');
title('Low/high gamma power (all trials)');
legend({'Low gamma', 'High gamma'}, 'Location', 'best');
grid on; box on;

nexttile; hold on;
mu_has_lo = nan(nCond,1); sem_has_lo = nan(nCond,1);
mu_has_hi = nan(nCond,1); sem_has_hi = nan(nCond,1);
for cond = 1:nCond
    v_lo = T_summary.HasLowPeakRate(T_summary.Condition == cond);
    v_hi = T_summary.HasHighPeakRate(T_summary.Condition == cond);
    mu_has_lo(cond) = nanmean(v_lo) * 100; sem_has_lo(cond) = nanstd(v_lo) / sqrt(sum(~isnan(v_lo))) * 100;
    mu_has_hi(cond) = nanmean(v_hi) * 100; sem_has_hi(cond) = nanstd(v_hi) / sqrt(sum(~isnan(v_hi))) * 100;
end
errorbar(contrast_vals, mu_has_lo, sem_has_lo, '-o', 'Color', [0.2 0.4 0.85], ...
    'LineWidth', 2.5, 'MarkerFaceColor', [0.2 0.4 0.85]);
errorbar(contrast_vals, mu_has_hi, sem_has_hi, '-s', 'Color', [0.85 0.25 0.25], ...
    'LineWidth', 2.5, 'MarkerFaceColor', [0.85 0.25 0.25]);
xlabel('Contrast [%]'); ylabel('Trials with peak [%]');
title('Peak presence rate (single-peak trials retained)');
legend({'Has low peak', 'Has high peak'}, 'Location', 'best');
ylim([0 100]);
grid on; box on;

sgtitle('SplitGammaLowHigh — Main Outcomes', 'FontSize', 18, 'FontWeight', 'bold');
saveas(fig3, fullfile(fig_dir, 'splitGammaLowHigh_mainOutcomes.png'));
close(fig3);

%% Grand figure 4: boxplots of low/high peak frequencies by condition
fig4 = figure('Position', [0 0 1512 982], 'Color', 'w');
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Build grouped vectors for boxplots (trial-level pooled)
g_low = [];
y_low = [];
g_high = [];
y_high = [];
for cond = 1:nCond
    v_lo = T_trials.LowPeakFreq(T_trials.Condition == cond);
    v_lo = v_lo(~isnan(v_lo));
    g_low = [g_low; repmat(cond, numel(v_lo), 1)]; %#ok<AGROW>
    y_low = [y_low; v_lo]; %#ok<AGROW>

    v_hi = T_trials.HighPeakFreq(T_trials.Condition == cond);
    v_hi = v_hi(~isnan(v_hi));
    g_high = [g_high; repmat(cond, numel(v_hi), 1)]; %#ok<AGROW>
    y_high = [y_high; v_hi]; %#ok<AGROW>
end

nexttile; hold on;
if ~isempty(y_low)
    boxplot(y_low, g_low, 'Labels', condLabels, 'Symbol', '');
    for cond = 1:nCond
        vals = y_low(g_low == cond);
        if isempty(vals), continue; end
        xj = cond + (rand(size(vals)) - 0.5) * 0.18;
        scatter(xj, vals, 10, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.20);
    end
end
xlabel('Contrast condition');
ylabel('Low gamma peak frequency [Hz]');
title('Low gamma peak frequencies across conditions');
ylim([scan_freqs(1), max(params.fallback_split_hz + 5, 55)]);
grid on; box on;

nexttile; hold on;
if ~isempty(y_high)
    boxplot(y_high, g_high, 'Labels', condLabels, 'Symbol', '');
    for cond = 1:nCond
        vals = y_high(g_high == cond);
        if isempty(vals), continue; end
        xj = cond + (rand(size(vals)) - 0.5) * 0.18;
        scatter(xj, vals, 10, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.20);
    end
end
xlabel('Contrast condition');
ylabel('High gamma peak frequency [Hz]');
title('High gamma peak frequencies across conditions');
ylim([max(params.fallback_split_hz - 5, 45), scan_freqs(end)]);
grid on; box on;

sgtitle('SplitGammaLowHigh — Peak Frequency Boxplots', 'FontSize', 18, 'FontWeight', 'bold');
saveas(fig4, fullfile(fig_dir, 'splitGammaLowHigh_peakFrequency_boxplots.png'));
close(fig4);

%% Save one consolidated MAT bundle
all_out = fullfile(data_dir, 'GCP_eeg_GED_splitGammaLowHigh_all.mat');
save(all_out, ...
    'T_split', 'T_trials', 'T_summary', 'params', ...
    'scan_freqs', 'subjects', 'condLabels', 'contrast_vals', ...
    'subject_mean_dt');

fprintf('\nDone: SplitGammaLowHigh outputs created.\n');
fprintf('Data folder:\n  %s\n', data_dir);
fprintf('Figures folder:\n  %s\n', fig_dir);
fprintf('Key files:\n');
fprintf('  %s\n', split_csv);
fprintf('  %s\n', trials_csv);
fprintf('  %s\n', summary_csv);
