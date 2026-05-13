%% GCP Test GED microsaccade covariance regression
% Trial weighted GED that prioritizes gamma covariance patterns during trials
% with higher microsaccade rate.

clear; clc; close all;
startup;
[subjects, paths, colors] = setup('GCP', 0);

out_data_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/tests';
out_fig_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/tests';
if ~exist(out_data_dir, 'dir'), mkdir(out_data_dir); end
if ~exist(out_fig_dir, 'dir'), mkdir(out_fig_dir); end

analysis_freq_range = [30 90];
stim_window = [0 2.0];
lambda = 0.01;
cond_codes = [61 62 63 64];

merged_csv = fullfile(paths.features, 'GCP_merged_data_trials.csv');
if ~isfile(merged_csv)
    error('Merged trial table not found: %s', merged_csv);
end
T = readtable(merged_csv);

id_var = pick_first_var(T, {'ID'});
cond_var = pick_first_var(T, {'Condition'});
trial_var = pick_first_var(T, {'Trial'});
ms_var = pick_first_var(T, {'Gaze_MSRate', 'Gaze_PctMSRate', 'MSRate', 'PctMSRate'});
if isempty(ms_var)
    error('No microsaccade metric column found in merged table.');
end

n_subj = numel(subjects);
subject_id = nan(n_subj, 1);
coupling_r = nan(n_subj, 1);
coupling_p = nan(n_subj, 1);
n_trials_used = zeros(n_subj, 1);
cond_component_power = nan(4, n_subj);
cond_ms_metric = nan(4, n_subj);

for si = 1:n_subj
    subj_str = subjects{si};
    sid = str2double(subj_str);
    subject_id(si) = sid;

    eeg_path = fullfile(paths.features, subj_str, 'eeg', 'dataEEG.mat');
    if ~isfile(eeg_path)
        continue;
    end
    S = load(eeg_path, 'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    data_structs = {S.dataEEG_c25, S.dataEEG_c50, S.dataEEG_c75, S.dataEEG_c100};

    n_ch = numel(data_structs{1}.label);
    cov_ref = zeros(n_ch);
    cov_sig = zeros(n_ch);
    trial_store = struct('x', {}, 'y', {}, 'cond', {});

    for ci = 1:4
        dat = data_structs{ci};
        trl_idx = find(dat.trialinfo == cond_codes(ci));
        if isempty(trl_idx), continue; end

        cfg = [];
        cfg.trials = trl_idx;
        dat = ft_selectdata(cfg, dat);

        cfgf = [];
        cfgf.bpfilter = 'yes';
        cfgf.bpfreq = analysis_freq_range;
        cfgf.bpfilttype = 'fir';
        cfgf.bpfiltord = round(3 * dat.fsample / analysis_freq_range(1));
        dat_g = ft_preprocessing(cfgf, dat);

        Ts = T(T.(id_var) == sid & T.(cond_var) == ci, :);
        y_by_trial = nan(numel(dat_g.trial), 1);
        if ~isempty(Ts)
            idx_ok = Ts.(trial_var) >= 1 & Ts.(trial_var) <= numel(dat_g.trial);
            if any(idx_ok)
                y_by_trial(Ts.(trial_var)(idx_ok)) = Ts.(ms_var)(idx_ok);
            end
        end

        y_valid = isfinite(y_by_trial);
        if sum(y_valid) < 5
            continue;
        end
        y = y_by_trial(y_valid);
        y = y - mean(y, 'omitnan');
        y = y ./ max(std(y, 0, 'omitnan'), eps);
        y01 = y - min(y);
        y01 = y01 ./ max(max(y01), eps);

        kept = find(y_valid);
        for ki = 1:numel(kept)
            trl = kept(ki);
            x = double(dat_g.trial{trl});
            t = dat_g.time{trl};
            idx_t = t >= stim_window(1) & t <= stim_window(2);
            x = x(:, idx_t);
            if isempty(x) || size(x, 2) < 5
                continue;
            end
            x = x - mean(x, 2);
            C = (x * x') / max(size(x, 2) - 1, 1);
            cov_ref = cov_ref + C;
            cov_sig = cov_sig + y01(ki) * C;

            trial_store(end+1).x = x; %#ok<SAGROW>
            trial_store(end).y = y(ki);
            trial_store(end).cond = ci;
        end
    end

    if numel(trial_store) < 20
        continue;
    end

    cov_ref = cov_ref / numel(trial_store);
    cov_sig = cov_sig / numel(trial_store);
    cov_ref_reg = (1 - lambda) * cov_ref + lambda * mean(diag(cov_ref)) * eye(n_ch);
    cov_sig_reg = (1 - lambda) * cov_sig + lambda * mean(diag(cov_sig)) * eye(n_ch);

    [W, D] = eig(cov_sig_reg, cov_ref_reg);
    [~, ix] = max(real(diag(D)));
    w = real(W(:, ix));
    w = w / max(norm(w), eps);

    comp_pow = nan(numel(trial_store), 1);
    y_all = nan(numel(trial_store), 1);
    cond_all = nan(numel(trial_store), 1);
    for ti = 1:numel(trial_store)
        z = w' * trial_store(ti).x;
        comp_pow(ti) = mean(z.^2, 'omitnan');
        y_all(ti) = trial_store(ti).y;
        cond_all(ti) = trial_store(ti).cond;
    end

    good = isfinite(comp_pow) & isfinite(y_all);
    if sum(good) < 10
        continue;
    end
    [r, p] = corr(comp_pow(good), y_all(good), 'type', 'Spearman');
    coupling_r(si) = r;
    coupling_p(si) = p;
    n_trials_used(si) = sum(good);

    for ci = 1:4
        m = good & cond_all == ci;
        if any(m)
            cond_component_power(ci, si) = mean(comp_pow(m), 'omitnan');
            cond_ms_metric(ci, si) = mean(y_all(m), 'omitnan');
        end
    end
end

R = table(subject_id, coupling_r, coupling_p, n_trials_used, ...
    'VariableNames', {'subject', 'coupling_r_spearman', 'coupling_p', 'n_trials_used'});
writetable(R, fullfile(out_data_dir, 'GCP_test_GED_microsaccade_spoc_subject_summary.csv'));
save(fullfile(out_data_dir, 'GCP_test_GED_microsaccade_spoc.mat'), ...
    'R', 'cond_component_power', 'cond_ms_metric', 'subjects', ...
    'analysis_freq_range', 'stim_window', 'lambda');

fig = figure('Position', [0 0 1512 982]);
set(fig, 'Color', 'w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
vals = coupling_r(isfinite(coupling_r));
if isempty(vals), vals = NaN; end
histogram(vals, 12, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none');
xlabel('Subject level coupling (Spearman r)');
ylabel('Count');
title('GED microsaccade coupling');

nexttile;
means = mean(cond_component_power, 2, 'omitnan');
sems = std(cond_component_power, 0, 2, 'omitnan') ./ sqrt(sum(isfinite(cond_component_power), 2));
hold on;
for ci = 1:4
    bar(ci, means(ci), 'FaceColor', colors(ci, :), 'EdgeColor', 'none');
end
errorbar(1:4, means, sems, 'k.', 'LineWidth', 1.2);
xlim([0.5 4.5]);
set(gca, 'XTick', 1:4, 'XTickLabel', {'25', '50', '75', '100'});
xlabel('Contrast');
ylabel('Component power');
title('Group mean component power by contrast');
hold off;

exportgraphics(fig, fullfile(out_fig_dir, 'GCP_test_GED_microsaccade_spoc.png'), 'Resolution', 300);

function name = pick_first_var(T, candidates)
name = '';
for i = 1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        name = candidates{i};
        return;
    end
end
end
