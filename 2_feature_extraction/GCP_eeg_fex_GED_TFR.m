%% GCP GED-based Time-Frequency Representation (TFR)
%
% This script reconstructs subject-specific GED spatial filters from the
% original EEG epochs, applies the selected GED component weights saved by
% the GED pipeline, and computes trial-resolved TFRs for each condition.
%
% Output:
%   data/features/GCP_eeg_GED_TFR.mat
%
% Notes:
%   - GED component indices/weights are loaded from GCP_eeg_GED.mat.
%   - Baseline correction follows an ERS/ERD-style subtraction in dB:
%       ERS/ERD(t) = power_dB(t) - mean(power_dB(baseline)).

%% Setup
startup
[subjects, paths, ~] = setup('GCP');

%% Parameters
baseline_window = [-1.5, -0.25];
stim_window = [0, 2.0];
gamma_range = [30, 90];
lambda = 0.05;

% TFR settings
tfr_foi = 30:1:90;
tfr_toi = -1.75:0.05:2.00;
tfr_win_sec = 0.50;
tfr_tapsmofrq = 5;

condNames = {'c25', 'c50', 'c75', 'c100'};
condCodes = [61, 62, 63, 64];
condLabels = {'25%', '50%', '75%', '100%'};

%% Paths
ged_path = fullfile(paths.features, 'GCP_eeg_GED.mat');
if ~isfile(ged_path)
    error('Missing GED features file: %s', ged_path);
end

Sged = load(ged_path, ...
    'all_component_selection_stats_full', ...
    'all_component_selection_stats', ...
    'subjects');
if isfield(Sged, 'all_component_selection_stats_full') && ~isempty(Sged.all_component_selection_stats_full)
    comp_stats_all = Sged.all_component_selection_stats_full;
elseif isfield(Sged, 'all_component_selection_stats') && ~isempty(Sged.all_component_selection_stats)
    comp_stats_all = Sged.all_component_selection_stats;
else
    error('Component selection stats are missing in %s', ged_path);
end

out_path = fullfile(paths.features, 'GCP_eeg_GED_TFR.mat');

%% Storage
nSubj = numel(subjects);
nCond = numel(condNames);

tfr_cond_trials = cell(nCond, nSubj);
tfr_cond_avg = cell(nCond, nSubj);
ged_filter_meta = cell(1, nSubj);

%% Subject loop
for subj = 1:nSubj
    fprintf('GED-TFR Subject %s (%d/%d)\n', subjects{subj}, subj, nSubj);

    % Load subject EEG
    subj_eeg_path = fullfile(paths.features, subjects{subj}, 'eeg', 'dataEEG.mat');
    if ~isfile(subj_eeg_path)
        warning('Missing EEG file for %s: %s', subjects{subj}, subj_eeg_path);
        continue;
    end
    E = load(subj_eeg_path, 'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    dat_by_cond = {E.dataEEG_c25, E.dataEEG_c50, E.dataEEG_c75, E.dataEEG_c100};

    fs = dat_by_cond{1}.fsample;
    nChan = numel(dat_by_cond{1}.label);

    % Recover selected GED component indices and weights (full window)
    stat_full = comp_stats_all{subj};
    if isempty(stat_full) || ~isfield(stat_full, 'selected_idx') || isempty(stat_full.selected_idx)
        warning('No selected GED components for %s. Skipping subject.', subjects{subj});
        continue;
    end
    sel_idx = stat_full.selected_idx(:)';
    sel_w = stat_full.selected_weights(:)';
    if isempty(sel_w) || any(~isfinite(sel_w))
        sel_w = ones(size(sel_idx));
    end
    sel_w = sel_w / max(sum(sel_w), eps);

    % Reconstruct full-window GED filters (deterministic, same settings)
    [W_sorted, evals_sorted] = reconstruct_ged_filters_fullwindow( ...
        dat_by_cond, condCodes, baseline_window, stim_window, gamma_range, lambda, fs, nChan);
    if isempty(W_sorted)
        warning('Failed GED reconstruction for %s. Skipping subject.', subjects{subj});
        continue;
    end
    sel_idx = sel_idx(sel_idx >= 1 & sel_idx <= size(W_sorted, 2));
    if isempty(sel_idx)
        warning('Selected GED indices are out of bounds for %s. Skipping subject.', subjects{subj});
        continue;
    end
    sel_w = sel_w(1:numel(sel_idx));
    sel_w = sel_w / max(sum(sel_w), eps);
    W_sel = W_sorted(:, sel_idx);

    ged_filter_meta{subj} = struct( ...
        'subject', subjects{subj}, ...
        'selected_idx', sel_idx, ...
        'selected_weights', sel_w, ...
        'selected_eigenvalues', evals_sorted(sel_idx));

    % Condition loop
    for c = 1:nCond
        dat = dat_by_cond{c};
        if isempty(dat) || ~isfield(dat, 'trial') || isempty(dat.trial)
            continue;
        end

        % Keep only task trials matching condition code
        trl_idx = find(dat.trialinfo == condCodes(c));
        if isempty(trl_idx)
            continue;
        end
        cfg_sel = [];
        cfg_sel.trials = trl_idx;
        dat = ft_selectdata(cfg_sel, dat);

        % Project each trial through selected GED components and combine
        ged_dat = dat;
        ged_dat.label = {'GED'};
        for tr = 1:numel(dat.trial)
            x = double(dat.trial{tr}); % chan x time
            y_comp = W_sel' * x;       % comp x time
            y = sel_w(:)' * y_comp;    % 1 x time
            ged_dat.trial{tr} = y;
        end

        % Trial-resolved TFR
        cfg_tfr = [];
        cfg_tfr.method = 'mtmconvol';
        cfg_tfr.output = 'pow';
        cfg_tfr.taper = 'dpss';
        cfg_tfr.foi = tfr_foi;
        cfg_tfr.toi = tfr_toi;
        cfg_tfr.t_ftimwin = tfr_win_sec * ones(size(tfr_foi));
        cfg_tfr.tapsmofrq = tfr_tapsmofrq * ones(size(tfr_foi));
        cfg_tfr.keeptrials = 'yes';
        tfr_trials = ft_freqanalysis(cfg_tfr, ged_dat);
        tfr_trials = baseline_correct_db_subtract(tfr_trials, baseline_window);

        % Subject-condition average
        cfg_avg = [];
        cfg_avg.keeptrials = 'no';
        tfr_avg = ft_freqdescriptives(cfg_avg, tfr_trials);

        tfr_cond_trials{c, subj} = tfr_trials;
        tfr_cond_avg{c, subj} = tfr_avg;
    end
end

%% Save
save(out_path, ...
    'tfr_cond_trials', 'tfr_cond_avg', ...
    'ged_filter_meta', ...
    'subjects', 'condNames', 'condLabels', ...
    'baseline_window', 'stim_window', ...
    'tfr_foi', 'tfr_toi', 'tfr_win_sec', 'tfr_tapsmofrq', ...
    '-v7.3');

fprintf('Saved GED-TFR data to: %s\n', out_path);

%% Local functions
function [W_sorted, evals_sorted] = reconstruct_ged_filters_fullwindow(dat_by_cond, condCodes, baseline_window, stim_window, gamma_range, lambda, fs, nChan)
cov_stim = zeros(nChan);
cov_base = zeros(nChan);
n_total = 0;

for c = 1:numel(dat_by_cond)
    dat = dat_by_cond{c};
    if isempty(dat) || ~isfield(dat, 'trial') || isempty(dat.trial)
        continue;
    end
    trl_idx = find(dat.trialinfo == condCodes(c));
    if isempty(trl_idx)
        continue;
    end

    cfg_sel = [];
    cfg_sel.trials = trl_idx;
    dat = ft_selectdata(cfg_sel, dat);

    cfg_filt = [];
    cfg_filt.bpfilter = 'yes';
    cfg_filt.bpfreq = gamma_range;
    cfg_filt.bpfilttype = 'fir';
    cfg_filt.bpfiltord = round(3 * fs / gamma_range(1));
    dat_gamma = ft_preprocessing(cfg_filt, dat);

    cfg_t = [];
    cfg_t.latency = baseline_window;
    dat_base = ft_selectdata(cfg_t, dat_gamma);
    cfg_t.latency = stim_window;
    dat_stim = ft_selectdata(cfg_t, dat_gamma);

    nTr = numel(dat_stim.trial);
    for tr = 1:nTr
        xb = double(dat_base.trial{tr});
        xs = double(dat_stim.trial{tr});
        xb = xb - mean(xb, 2);
        xs = xs - mean(xs, 2);
        cov_base = cov_base + (xb * xb') / size(xb, 2);
        cov_stim = cov_stim + (xs * xs') / size(xs, 2);
    end
    n_total = n_total + nTr;
end

if n_total < 1
    W_sorted = [];
    evals_sorted = [];
    return;
end

cov_base = cov_base / n_total;
cov_stim = cov_stim / n_total;
cov_stim_reg = (1 - lambda) * cov_stim + lambda * mean(diag(cov_stim)) * eye(nChan);
cov_base_reg = (1 - lambda) * cov_base + lambda * mean(diag(cov_base)) * eye(nChan);

[W, D] = eig(cov_stim_reg, cov_base_reg);
[evals_sorted, idx] = sort(real(diag(D)), 'descend');
W_sorted = W(:, idx);
end

function tfr_out = baseline_correct_db_subtract(tfr_in, baseline_window)
tfr_out = tfr_in;
pow = tfr_in.powspctrm;
pow = max(pow, eps);
pow_db = 10 * log10(pow);

t = tfr_in.time;
bmask = t >= baseline_window(1) & t <= baseline_window(2);
if ~any(bmask)
    warning('Baseline window not found in TFR time axis. Returning dB-converted power only.');
    tfr_out.powspctrm = pow_db;
    return;
end

% Supports keeptrials='yes' (rpt x chan x freq x time) and no-trial dims.
if ndims(pow_db) == 4
    base = mean(pow_db(:, :, :, bmask), 4, 'omitnan');
    tfr_out.powspctrm = pow_db - base;
elseif ndims(pow_db) == 3
    base = mean(pow_db(:, :, bmask), 3, 'omitnan');
    tfr_out.powspctrm = pow_db - base;
else
    tfr_out.powspctrm = pow_db;
end
end
