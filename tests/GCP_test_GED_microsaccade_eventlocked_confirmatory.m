%% GCP event-locked confirmatory GED around Engbert microsaccades
%
% Uses the FIXED occipital GED component from GCP_eeg_GED.mat (not a new
% MS-versus-control GED). Locks component gamma power to Engbert MS onsets
% from gaze_microsaccade_events.mat (fallback: Engbert redetect on dataET).
%
% This is a diagnostic for neural vs artifactual accounts: brief broadband
% MS-locked bursts favor spike artifact; sustained narrowband stimulus-locked
% gamma that survives MS-free analysis favors neural gamma.

clear; clc; close all;
startup;
[subjects, paths, colors] = setup('GCP', 0);

out_data_dir = fullfile(paths.features, 'tests');
out_fig_dir = fullfile(paths.figures, 'tests');
if ~exist(out_data_dir, 'dir'), mkdir(out_data_dir); end
if ~exist(out_fig_dir, 'dir'), mkdir(out_fig_dir); end

ged_path = fullfile(paths.features, 'GCP_eeg_GED.mat');
if ~isfile(ged_path)
    error('Missing %s. Run GCP_eeg_fex_GED.m first.', ged_path);
end
GED = load(ged_path, 'all_combined_filter_full', 'subjects');
if isfield(GED, 'subjects') && ~isempty(GED.subjects)
    subjects = GED.subjects;
end

event_window = [-0.15 0.15];
stim_window = [0 2.0];
condNames = {'c25', 'c50', 'c75', 'c100'};
n_subj = numel(subjects);

subject_id = nan(n_subj, 1);
n_events = zeros(n_subj, 1);
peak_event_power = nan(4, n_subj);
time_vec_common = linspace(event_window(1), event_window(2), round(diff(event_window) * 500) + 1);
subject_curves = nan(n_subj, 4, numel(time_vec_common));

for si = 1:n_subj
    subj = subjects{si};
    sid = str2double(subj);
    subject_id(si) = sid;

    w = [];
    if numel(GED.all_combined_filter_full) >= si
        w = GED.all_combined_filter_full{si};
    end
    if isempty(w) || ~all(isfinite(w(:)))
        continue
    end
    w = w(:);

    eeg_path = fullfile(paths.features, subj, 'eeg', 'dataEEG.mat');
    if ~isfile(eeg_path), continue; end
    E = load(eeg_path, 'dataEEG_c25', 'dataEEG_c50', 'dataEEG_c75', 'dataEEG_c100');
    data_eeg = {E.dataEEG_c25, E.dataEEG_c50, E.dataEEG_c75, E.dataEEG_c100};

    ms_events = load_ms_events_local(paths.features, subj, condNames);

    for ci = 1:4
        dat = data_eeg{ci};
        if isempty(dat), continue; end
        if size(dat.trial{1}, 1) ~= numel(w), continue; end

        onsets_cond = {};
        if ~isempty(ms_events) && numel(ms_events) >= ci && isfield(ms_events{ci}, 'Onset')
            onsets_cond = ms_events{ci}.Onset;
        end

        event_mat = [];
        n_trials = numel(dat.trial);
        for trl = 1:n_trials
            x = double(dat.trial{trl});
            t = dat.time{trl};
            z = w' * x;
            onsets = [];
            if iscell(onsets_cond) && trl <= numel(onsets_cond)
                onsets = onsets_cond{trl}(:);
            end
            onsets = onsets(isfinite(onsets) & onsets >= stim_window(1) & onsets <= stim_window(2));
            for oi = 1:numel(onsets)
                seg = extract_segment(z, t, onsets(oi), event_window);
                if isempty(seg), continue; end
                % interpolate onto common time grid via length match
                if numel(seg) ~= numel(time_vec_common)
                    seg = interp1(linspace(event_window(1), event_window(2), numel(seg)), ...
                        seg, time_vec_common, 'linear', NaN);
                end
                event_mat = [event_mat; seg.^2]; %#ok<AGROW>
            end
        end

        if ~isempty(event_mat)
            n_events(si) = n_events(si) + size(event_mat, 1);
            avg_curve = mean(event_mat, 1, 'omitnan');
            subject_curves(si, ci, :) = avg_curve;
            peak_event_power(ci, si) = max(avg_curve);
        end
    end
end

R = table(subject_id, n_events, ...
    peak_event_power(1, :)', peak_event_power(2, :)', ...
    peak_event_power(3, :)', peak_event_power(4, :)', ...
    'VariableNames', {'subject', 'n_events', ...
    'peak_pow_c25', 'peak_pow_c50', 'peak_pow_c75', 'peak_pow_c100'});
writetable(R, fullfile(out_data_dir, 'GCP_test_GED_microsaccade_eventlocked_confirmatory.csv'));
save(fullfile(out_data_dir, 'GCP_test_GED_microsaccade_eventlocked_confirmatory.mat'), ...
    'R', 'subject_curves', 'time_vec_common', 'subjects', 'event_window', 'stim_window');

fig = figure('Position', [0 0 1512 982]);
set(fig, 'Color', 'w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold on;
for ci = 1:4
    Y = squeeze(subject_curves(:, ci, :));
    mu = mean(Y, 1, 'omitnan');
    se = std(Y, 0, 1, 'omitnan') ./ sqrt(max(sum(isfinite(Y), 1), 1));
    fill([time_vec_common, fliplr(time_vec_common)], ...
         [mu + se, fliplr(mu - se)], colors(ci, :), ...
         'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(time_vec_common, mu, 'Color', colors(ci, :), 'LineWidth', 2.0);
end
xline(0, 'k:', 'LineWidth', 1.1);
xlabel('Time around microsaccade (s)');
ylabel('Confirmatory GED gamma power');
title('Event locked confirmatory GED');
legend({'25', '50', '75', '100'}, 'Location', 'best', 'Box', 'off');
hold off;

nexttile;
m = mean(peak_event_power, 2, 'omitnan');
s = std(peak_event_power, 0, 2, 'omitnan') ./ sqrt(max(sum(isfinite(peak_event_power), 2), 1));
hold on;
for ci = 1:4
    bar(ci, m(ci), 'FaceColor', colors(ci, :), 'EdgeColor', 'none');
end
errorbar(1:4, m, s, 'k.', 'LineWidth', 1.2);
xlim([0.5 4.5]);
set(gca, 'XTick', 1:4, 'XTickLabel', {'25', '50', '75', '100'});
xlabel('Contrast');
ylabel('Peak event locked power');
title('Contrast effect on MS-locked power');
hold off;

exportgraphics(fig, fullfile(out_fig_dir, 'GCP_test_GED_microsaccade_eventlocked_confirmatory.png'), 'Resolution', 300);
fprintf('Saved confirmatory event-locked outputs to %s\n', out_fig_dir);

function ms_events = load_ms_events_local(features_root, subj, condNames)
ms_events = cell(1, 4);
ev_path = fullfile(features_root, subj, 'gaze', 'gaze_microsaccade_events.mat');
if isfile(ev_path)
    S = load(ev_path);
    for ci = 1:4
        fn = sprintf('ms_events_%s', condNames{ci});
        if isfield(S, fn)
            ms_events{ci} = S.(fn);
        end
    end
    return
end
et_path = fullfile(features_root, subj, 'gaze', 'dataET.mat');
if ~isfile(et_path)
    return
end
G = load(et_path);
data_et = {G.dataET_c25, G.dataET_c50, G.dataET_c75, G.dataET_c100};
for ci = 1:4
    gaze = data_et{ci};
    nTrl = numel(gaze.trial);
    onsets = cell(1, nTrl);
    offsets = cell(1, nTrl);
    for trl = 1:nTrl
        [on_s, off_s] = detect_engbert_ms_times_local(gaze.trial{trl}, gaze.time{trl}, gaze.fsample);
        onsets{trl} = on_s;
        offsets{trl} = off_s;
    end
    ms_events{ci} = struct();
    ms_events{ci}.Onset = onsets;
    ms_events{ci}.Offset = offsets;
end
end

function [onset_s, offset_s] = detect_engbert_ms_times_local(raw, tVec, fsample)
onset_s = [];
offset_s = [];
if isempty(raw) || isempty(tVec), return; end
full_idx = tVec >= -1.5 & tVec <= 2.0;
t_full = tVec(full_idx);
x = raw(1, full_idx);
y = 600 - raw(2, full_idx);
valid = x >= 0 & x <= 800 & y >= 0 & y <= 600 & isfinite(x) & isfinite(y);
if sum(valid) < round(0.5 * fsample), return; end
x_clean = x(valid);
y_clean = y(valid);
[~, ms_det] = detect_microsaccades(fsample, [x_clean; y_clean], numel(x_clean));
if isempty(ms_det.Onset), return; end
idx_full_valid = find(valid);
nEv = min(numel(ms_det.Onset), numel(ms_det.Offset));
on_idx = ms_det.Onset(1:nEv);
off_idx = ms_det.Offset(1:nEv);
keep = on_idx >= 1 & on_idx <= numel(idx_full_valid) & ...
       off_idx >= 1 & off_idx <= numel(idx_full_valid);
onset_s = t_full(idx_full_valid(on_idx(keep)))';
offset_s = t_full(idx_full_valid(off_idx(keep)))';
end

function seg = extract_segment(sig, t, center_t, win)
idx = t >= center_t + win(1) & t <= center_t + win(2);
seg = sig(idx);
if numel(seg) < 10
    seg = [];
end
end
