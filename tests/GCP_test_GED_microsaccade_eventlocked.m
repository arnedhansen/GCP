%% GCP Test event locked GED for microsaccades by contrast
% GED is built from peri microsaccade segments versus matched control segments.
% The resulting component is evaluated with event locked gamma power curves.

clear; clc; close all;
startup;
[subjects, paths, colors] = setup('GCP', 0);
rng(123, 'twister');

out_data_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/tests';
out_fig_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/tests';
if ~exist(out_data_dir, 'dir'), mkdir(out_data_dir); end
if ~exist(out_fig_dir, 'dir'), mkdir(out_fig_dir); end

analysis_freq_range = [30 90];
stim_window = [0 2.0];
event_window = [-0.15 0.15];
guard_s = 0.20;
lambda = 0.01;
cond_codes = [61 62 63 64];

n_subj = numel(subjects);
subject_id = nan(n_subj, 1);
n_events_ms = zeros(n_subj, 1);
n_events_ctl = zeros(n_subj, 1);
peak_event_power = nan(4, n_subj);
time_vec_common = [];
subject_curves = [];

for si = 1:n_subj
    subj_str = subjects{si};
    sid = str2double(subj_str);
    subject_id(si) = sid;

    eeg_path = fullfile(paths.features, subj_str, 'eeg', 'dataEEG.mat');
    gaze_path = fullfile(paths.features, subj_str, 'gaze', 'dataET.mat');
    if ~isfile(eeg_path) || ~isfile(gaze_path)
        continue;
    end

    E = load(eeg_path, 'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    G = load(gaze_path, 'dataET_c25', 'dataET_c50', 'dataET_c75', 'dataET_c100');
    data_eeg = {E.dataEEG_c25, E.dataEEG_c50, E.dataEEG_c75, E.dataEEG_c100};
    data_gaze = {G.dataET_c25, G.dataET_c50, G.dataET_c75, G.dataET_c100};

    n_ch = numel(data_eeg{1}.label);
    cov_ms = zeros(n_ch);
    cov_ctl = zeros(n_ch);
    n_seg_ms = 0;
    n_seg_ctl = 0;
    event_signals_by_cond = cell(1, 4);

    for ci = 1:4
        dat = data_eeg{ci};
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
        dat_gam = ft_preprocessing(cfgf, dat);

        gaze = data_gaze{ci};
        n_trials = min(numel(dat_gam.trial), numel(gaze.trial));
        if n_trials < 1
            continue;
        end

        for trl = 1:n_trials
            x_eeg = double(dat_gam.trial{trl});
            t_eeg = dat_gam.time{trl};
            x_eeg = x_eeg - mean(x_eeg, 2);

            g = double(gaze.trial{trl}(1:2, :));
            t_g = gaze.time{trl};
            ms_onsets = detect_ms_onsets(g, t_g, gaze.fsample, stim_window);
            if isempty(ms_onsets)
                continue;
            end

            valid_onsets = [];
            for oi = 1:numel(ms_onsets)
                c = ms_onsets(oi);
                if c + event_window(1) < min(t_eeg) || c + event_window(2) > max(t_eeg)
                    continue;
                end
                valid_onsets(end+1) = c; %#ok<SAGROW>
            end
            if isempty(valid_onsets)
                continue;
            end

            ctl_centers = sample_control_centers(t_eeg, valid_onsets, stim_window, event_window, guard_s, numel(valid_onsets));
            if isempty(ctl_centers)
                continue;
            end

            for oi = 1:numel(valid_onsets)
                seg = extract_segment(x_eeg, t_eeg, valid_onsets(oi), event_window);
                if isempty(seg), continue; end
                C = (seg * seg') / max(size(seg, 2) - 1, 1);
                cov_ms = cov_ms + C;
                n_seg_ms = n_seg_ms + 1;
            end

            n_pair = min(numel(valid_onsets), numel(ctl_centers));
            for oi = 1:n_pair
                seg = extract_segment(x_eeg, t_eeg, ctl_centers(oi), event_window);
                if isempty(seg), continue; end
                C = (seg * seg') / max(size(seg, 2) - 1, 1);
                cov_ctl = cov_ctl + C;
                n_seg_ctl = n_seg_ctl + 1;
            end
        end
    end

    if n_seg_ms < 20 || n_seg_ctl < 20
        continue;
    end

    cov_ms = cov_ms / n_seg_ms;
    cov_ctl = cov_ctl / n_seg_ctl;
    cov_ref = (1 - lambda) * cov_ctl + lambda * mean(diag(cov_ctl)) * eye(n_ch);
    cov_sig = (1 - lambda) * cov_ms + lambda * mean(diag(cov_ms)) * eye(n_ch);

    [W, D] = eig(cov_sig, cov_ref);
    [~, ix] = max(real(diag(D)));
    w = real(W(:, ix));
    w = w / max(norm(w), eps);

    for ci = 1:4
        dat = data_eeg{ci};
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
        dat_gam = ft_preprocessing(cfgf, dat);
        gaze = data_gaze{ci};
        n_trials = min(numel(dat_gam.trial), numel(gaze.trial));

        event_mat = [];
        for trl = 1:n_trials
            x_eeg = double(dat_gam.trial{trl});
            t_eeg = dat_gam.time{trl};
            z = w' * x_eeg;
            g = double(gaze.trial{trl}(1:2, :));
            t_g = gaze.time{trl};
            ms_onsets = detect_ms_onsets(g, t_g, gaze.fsample, stim_window);
            for oi = 1:numel(ms_onsets)
                seg = extract_segment(z, t_eeg, ms_onsets(oi), event_window);
                if isempty(seg), continue; end
                event_mat = [event_mat; seg.^2]; %#ok<AGROW>
            end
        end

        if ~isempty(event_mat)
            event_signals_by_cond{ci} = event_mat;
            avg_curve = mean(event_mat, 1, 'omitnan');
            [mx, ~] = max(avg_curve);
            peak_event_power(ci, si) = mx;
            if isempty(time_vec_common)
                time_vec_common = linspace(event_window(1), event_window(2), size(event_mat, 2));
            end
        end
    end

    if isempty(subject_curves)
        if ~isempty(time_vec_common)
            subject_curves = nan(n_subj, 4, numel(time_vec_common));
        end
    end
    if ~isempty(subject_curves)
        for ci = 1:4
            if ~isempty(event_signals_by_cond{ci})
                subject_curves(si, ci, :) = mean(event_signals_by_cond{ci}, 1, 'omitnan');
            end
        end
    end

    n_events_ms(si) = n_seg_ms;
    n_events_ctl(si) = n_seg_ctl;
end

R = table(subject_id, n_events_ms, n_events_ctl, ...
    peak_event_power(1, :)', peak_event_power(2, :)', peak_event_power(3, :)', peak_event_power(4, :)', ...
    'VariableNames', {'subject', 'n_event_segments', 'n_control_segments', ...
    'peak_pow_c25', 'peak_pow_c50', 'peak_pow_c75', 'peak_pow_c100'});
writetable(R, fullfile(out_data_dir, 'GCP_test_GED_microsaccade_eventlocked_subject_summary.csv'));
save(fullfile(out_data_dir, 'GCP_test_GED_microsaccade_eventlocked.mat'), ...
    'R', 'subject_curves', 'time_vec_common', 'subjects', ...
    'analysis_freq_range', 'stim_window', 'event_window', 'lambda');

fig = figure('Position', [0 0 1512 982]);
set(fig, 'Color', 'w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
if ~isempty(subject_curves) && ~isempty(time_vec_common)
    hold on;
    for ci = 1:4
        Y = squeeze(subject_curves(:, ci, :));
        mu = mean(Y, 1, 'omitnan');
        se = std(Y, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(Y), 1));
        fill([time_vec_common, fliplr(time_vec_common)], ...
             [mu + se, fliplr(mu - se)], colors(ci, :), ...
             'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(time_vec_common, mu, 'Color', colors(ci, :), 'LineWidth', 2.0);
    end
    xline(0, 'k:', 'LineWidth', 1.1);
    xlabel('Time around microsaccade (s)');
    ylabel('Component gamma power');
    title('Event locked GED power');
    legend({'25', '50', '75', '100'}, 'Location', 'best', 'Box', 'off');
    hold off;
else
    text(0.1, 0.5, 'No valid event locked curves', 'FontSize', 14);
    axis off;
end

nexttile;
m = mean(peak_event_power, 2, 'omitnan');
s = std(peak_event_power, 0, 2, 'omitnan') ./ sqrt(sum(isfinite(peak_event_power), 2));
hold on;
for ci = 1:4
    bar(ci, m(ci), 'FaceColor', colors(ci, :), 'EdgeColor', 'none');
end
errorbar(1:4, m, s, 'k.', 'LineWidth', 1.2);
xlim([0.5 4.5]);
set(gca, 'XTick', 1:4, 'XTickLabel', {'25', '50', '75', '100'});
xlabel('Contrast');
ylabel('Peak event locked power');
title('Contrast effect on event locked power');
hold off;

exportgraphics(fig, fullfile(out_fig_dir, 'GCP_test_GED_microsaccade_eventlocked.png'), 'Resolution', 300);

function ms_onsets = detect_ms_onsets(gaze_xy, t, fs, stim_window)
x = gaze_xy(1, :);
y = gaze_xy(2, :);
valid = x >= 0 & x <= 800 & y >= 0 & y <= 600 & isfinite(x) & isfinite(y) & isfinite(t);
x = x(valid);
y = y(valid);
t = t(valid);
if numel(x) < 10
    ms_onsets = [];
    return;
end

vx = [0, diff(x)] * fs;
vy = [0, diff(y)] * fs;
sx = median(abs(vx - median(vx))) / 0.6745;
sy = median(abs(vy - median(vy))) / 0.6745;
sx = max(sx, eps);
sy = max(sy, eps);
vnorm = sqrt((vx / sx).^2 + (vy / sy).^2);

thr = 6;
is_ms = vnorm > thr;
is_ms = is_ms & t >= stim_window(1) & t <= stim_window(2);

d = diff([0, is_ms, 0]);
starts = find(d == 1);
stops = find(d == -1) - 1;
min_dur = max(3, round(0.006 * fs));
ms_onsets = [];
for i = 1:numel(starts)
    if stops(i) - starts(i) + 1 >= min_dur
        ms_onsets(end+1) = t(starts(i)); %#ok<AGROW>
    end
end
end

function centers = sample_control_centers(t, onsets, stim_window, event_window, guard_s, n_need)
if isempty(t)
    centers = [];
    return;
end
cand = t(t >= stim_window(1) & t <= stim_window(2));
cand = cand(cand + event_window(1) >= min(t) & cand + event_window(2) <= max(t));
if isempty(cand)
    centers = [];
    return;
end
keep = true(size(cand));
for i = 1:numel(cand)
    if any(abs(cand(i) - onsets) <= guard_s)
        keep(i) = false;
    end
end
cand = cand(keep);
if isempty(cand)
    centers = [];
    return;
end
n_take = min(n_need, numel(cand));
ix = randperm(numel(cand), n_take);
centers = cand(ix);
end

function seg = extract_segment(sig, t, center_t, win)
if isvector(sig)
    idx = t >= center_t + win(1) & t <= center_t + win(2);
    seg = sig(idx);
    if numel(seg) < 10
        seg = [];
    end
else
    idx = t >= center_t + win(1) & t <= center_t + win(2);
    seg = sig(:, idx);
    if size(seg, 2) < 10
        seg = [];
    end
end
end
