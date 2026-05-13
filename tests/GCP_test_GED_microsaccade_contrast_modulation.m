%% GCP Test GED contrast modulation with microsaccades
% GED optimized for high versus low contrast, then tested for
% microsaccade by contrast interaction on component gamma power.

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
beta_ms = nan(n_subj, 1);
beta_contrast = nan(n_subj, 1);
beta_interaction = nan(n_subj, 1);
p_interaction = nan(n_subj, 1);
n_trials_model = zeros(n_subj, 1);

slopes_by_cond = nan(4, n_subj);

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
    cov_low = zeros(n_ch);
    cov_high = zeros(n_ch);
    n_low = 0;
    n_high = 0;
    trial_rows = struct('pow', {}, 'ms', {}, 'cond_idx', {});

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

        for trl = 1:numel(dat_g.trial)
            x = double(dat_g.trial{trl});
            t = dat_g.time{trl};
            idx_t = t >= stim_window(1) & t <= stim_window(2);
            x = x(:, idx_t);
            if isempty(x) || size(x, 2) < 5
                continue;
            end
            x = x - mean(x, 2);
            C = (x * x') / max(size(x, 2) - 1, 1);

            if ci == 1
                cov_low = cov_low + C;
                n_low = n_low + 1;
            elseif ci == 4
                cov_high = cov_high + C;
                n_high = n_high + 1;
            end

            trial_rows(end+1).pow = NaN; %#ok<SAGROW>
            trial_rows(end).ms = y_by_trial(trl);
            trial_rows(end).cond_idx = ci;
            trial_rows(end).x = x;
        end
    end

    if n_low < 10 || n_high < 10 || isempty(trial_rows)
        continue;
    end

    cov_low = cov_low / n_low;
    cov_high = cov_high / n_high;
    cov_ref = (1 - lambda) * cov_low + lambda * mean(diag(cov_low)) * eye(n_ch);
    cov_sig = (1 - lambda) * cov_high + lambda * mean(diag(cov_high)) * eye(n_ch);

    [W, D] = eig(cov_sig, cov_ref);
    [~, ix] = max(real(diag(D)));
    w = real(W(:, ix));
    w = w / max(norm(w), eps);

    n_rows = numel(trial_rows);
    comp_pow = nan(n_rows, 1);
    ms_all = nan(n_rows, 1);
    cond_all = nan(n_rows, 1);
    for ri = 1:n_rows
        z = w' * trial_rows(ri).x;
        comp_pow(ri) = mean(z.^2, 'omitnan');
        ms_all(ri) = trial_rows(ri).ms;
        cond_all(ri) = trial_rows(ri).cond_idx;
    end

    good = isfinite(comp_pow) & isfinite(ms_all) & isfinite(cond_all);
    if sum(good) < 20
        continue;
    end

    x_ms = zscore(ms_all(good));
    x_con = zscore(cond_all(good));
    X = [ones(sum(good), 1), x_ms, x_con, x_ms .* x_con];
    y = comp_pow(good);
    [b, ~, ~, ~, stats] = regress(y, X);

    beta_ms(si) = b(2);
    beta_contrast(si) = b(3);
    beta_interaction(si) = b(4);
    p_interaction(si) = stats(3);
    n_trials_model(si) = sum(good);

    for ci = 1:4
        m = good & cond_all == ci;
        if sum(m) >= 5
            xm = zscore(ms_all(m));
            ym = comp_pow(m);
            if all(isfinite(xm)) && all(isfinite(ym))
                b_ci = regress(ym, [ones(sum(m), 1), xm]);
                slopes_by_cond(ci, si) = b_ci(2);
            end
        end
    end
end

R = table(subject_id, beta_ms, beta_contrast, beta_interaction, p_interaction, n_trials_model, ...
    'VariableNames', {'subject', 'beta_ms', 'beta_contrast', 'beta_ms_x_contrast', 'p_model', 'n_trials'});
writetable(R, fullfile(out_data_dir, 'GCP_test_GED_microsaccade_contrast_modulation_subject_summary.csv'));
save(fullfile(out_data_dir, 'GCP_test_GED_microsaccade_contrast_modulation.mat'), ...
    'R', 'slopes_by_cond', 'subjects', 'analysis_freq_range', 'stim_window', 'lambda');

fig = figure('Position', [0 0 1512 982]);
set(fig, 'Color', 'w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
vals = beta_interaction(isfinite(beta_interaction));
if isempty(vals), vals = NaN; end
histogram(vals, 12, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none');
xlabel('Subject interaction beta (MS x contrast)');
ylabel('Count');
title('Contrast modulation of MS coupling');

nexttile;
m = mean(slopes_by_cond, 2, 'omitnan');
s = std(slopes_by_cond, 0, 2, 'omitnan') ./ sqrt(sum(isfinite(slopes_by_cond), 2));
hold on;
for ci = 1:4
    bar(ci, m(ci), 'FaceColor', colors(ci, :), 'EdgeColor', 'none');
end
errorbar(1:4, m, s, 'k.', 'LineWidth', 1.2);
xlim([0.5 4.5]);
set(gca, 'XTick', 1:4, 'XTickLabel', {'25', '50', '75', '100'});
xlabel('Contrast');
ylabel('Within contrast slope');
title('MS to power slope by contrast');
hold off;

exportgraphics(fig, fullfile(out_fig_dir, 'GCP_test_GED_microsaccade_contrast_modulation.png'), 'Resolution', 300);

function name = pick_first_var(T, candidates)
name = '';
for i = 1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        name = candidates{i};
        return;
    end
end
end
