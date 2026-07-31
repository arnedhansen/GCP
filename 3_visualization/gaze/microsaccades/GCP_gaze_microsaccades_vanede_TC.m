%% GCP Gaze Microsaccade Suppression Using the van Ede Method
% Detects gaze shifts from preprocessed trial-level gaze positions, creates
% baseline-normalised microsaccade-rate time courses, and plots condition
% time courses and participant-level condition summaries.
%
% Detection follows Liu, Nobre, and van Ede (2022):
% https://github.com/BaiweiLiu/VWMMS_EEG_alpha-covertAttention
%
% The released method is one-dimensional because the original task was
% horizontally lateralised. Here, two-dimensional speed is used because
% the GCP contrast stimulus is central and total microsaccade rate is the
% outcome of interest.

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
addpath('/Users/Arne/Documents/GitHub/functions')
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar')

datapath = paths.features;
figpath = fullfile(paths.figures, 'gaze', 'microsaccades', 'vanede');
mkdir(figpath);

condition_names = {'c25', 'c50', 'c75', 'c100'};
condition_labels = {'25% Contrast', '50% Contrast', '75% Contrast', '100% Contrast'};
n_conditions = numel(condition_names);
n_subjects = numel(subjects);

variant_names = { ...
    'bl_exclude_zero', ...
    'bl_include_zero', ...
    'raw_exclude_zero', ...
    'raw_include_zero'};
variant_ylabels = { ...
    'Microsaccade rate change [%]', ...
    'Microsaccade rate change [%]', ...
    'Microsaccade rate [Hz]', ...
    'Microsaccade rate [Hz]'};
variant_summary_ylabels = { ...
    'Mean microsaccade rate change, 0 to 2 s [%]', ...
    'Mean microsaccade rate change, 0 to 2 s [%]', ...
    'Mean microsaccade rate, 0 to 2 s [Hz]', ...
    'Mean microsaccade rate, 0 to 2 s [Hz]'};
variant_titles = { ...
    'Baselined, zero-stim excluded', ...
    'Baselined, zero-stim included', ...
    'Non-baselined, zero-stim excluded', ...
    'Non-baselined, zero-stim included'};
n_variants = numel(variant_names);

%% Analysis parameters
baseline_period = [-1.5 -0.5];
analysis_period = [0 2];
display_period = [-0.5 2];
blink_window_samples = 25;

detector_cfg = [];
detector_cfg.threshold = 5;
detector_cfg.smooth_ms = 7;
detector_cfg.before_ms = [-50 0];
detector_cfg.after_ms = [50 100];
detector_cfg.min_interval_ms = 100;

rate_window_ms = 100;
display_smooth_ms = 20;
minimum_valid_fraction = 0.60;
maximum_blink_loss_fraction = 0.40;
minimum_valid_duration_s = 0.50;
maximum_rate_hz = 8;
minimum_baseline_rate_hz = 0.1;

subject_curves = cell(n_subjects, n_conditions, n_variants);
subject_summary = nan(n_subjects, n_conditions, n_variants);
common_display_time = [];

%% Detect gaze shifts and construct subject time courses
fprintf('\n[VIZ GAZE MS VANEDE] Gaze-shift detection\n');

for subject_index = 1:n_subjects
    subject_id = subjects{subject_index};
    clc; fprintf('[VIZ GAZE MS VANEDE] Subject %d/%d (%s)\n', subject_index, n_subjects, subject_id);

    input_file = fullfile(datapath, subject_id, 'gaze', 'dataET.mat');
    if ~isfile(input_file)
        warning('GCP:MissingEyeData', 'File not found: %s', input_file);
        continue
    end
    eye_data = load(input_file);

    output_file = fullfile(datapath, subject_id, 'gaze', ...
        'gaze_microsaccade_timeseries_vanede.mat');
    output_data = struct();

    for condition_index = 1:n_conditions
        condition_name = condition_names{condition_index};
        data_field = ['dataET_' condition_name];
        if ~isfield(eye_data, data_field)
            warning('GCP:MissingCondition', 'No %s in %s', data_field, input_file);
            continue
        end

        data_condition = eye_data.(data_field);
        fsample = data_condition.fsample;
        common_time = baseline_period(1):(1 / fsample):analysis_period(2);
        n_trials = numel(data_condition.trial);
        trial_variants = cell(1, n_variants);
        for variant_index = 1:n_variants
            trial_variants{variant_index} = nan(n_trials, numel(common_time));
        end
        event_details = cell(n_trials, 1);

        rate_window_samples = max(1, round(rate_window_ms * fsample / 1000));
        rate_kernel = ones(1, rate_window_samples) / rate_window_samples;
        minimum_valid_samples = round(minimum_valid_duration_s * fsample);

        for trial_index = 1:n_trials
            trial_time = data_condition.time{trial_index};
            trial_data = data_condition.trial{trial_index};
            selected = trial_time >= baseline_period(1) & ...
                trial_time <= analysis_period(2);
            if ~any(selected) || size(trial_data, 1) < 2
                continue
            end

            time_selected = trial_time(selected);
            gaze = trial_data(1:2, selected);
            gaze(2, :) = 600 - gaze(2, :);

            in_screen = gaze(1, :) >= 0 & gaze(1, :) <= 800 & ...
                gaze(2, :) >= 0 & gaze(2, :) <= 600;
            gaze(:, ~in_screen) = NaN;
            valid_before_blinks = all(isfinite(gaze), 1);
            gaze = remove_blinks(gaze, blink_window_samples);
            valid_after_blinks = all(isfinite(gaze), 1);

            valid_fraction = sum(valid_after_blinks) / numel(valid_after_blinks);
            blink_loss_fraction = max(0, ...
                (sum(valid_before_blinks) - sum(valid_after_blinks)) / ...
                max(sum(valid_before_blinks), 1));

            if sum(valid_after_blinks) < minimum_valid_samples || ...
                    valid_fraction < minimum_valid_fraction || ...
                    blink_loss_fraction > maximum_blink_loss_fraction
                continue
            end

            [gaze_shift, details] = detect_gaze_shifts_van_ede( ...
                fsample, gaze, detector_cfg);

            event_impulses = gaze_shift ~= 0;
            local_valid_mass = conv(double(valid_after_blinks), ...
                rate_kernel, 'same');
            rate = conv(double(event_impulses), rate_kernel, 'same');
            rate = rate ./ local_valid_mass * fsample;
            rate(local_valid_mass == 0 | ~valid_after_blinks) = NaN;

            valid_duration = sum(valid_after_blinks) / fsample;
            trial_event_rate = sum(event_impulses) / valid_duration;
            if ~isfinite(trial_event_rate) || trial_event_rate > maximum_rate_hz
                continue
            end

            baseline_index = time_selected >= baseline_period(1) & ...
                time_selected <= baseline_period(2);
            analysis_index = time_selected >= analysis_period(1) & ...
                time_selected <= analysis_period(2);
            baseline_rate = mean(rate(baseline_index), 'omitnan');
            stimulus_rate = mean(rate(analysis_index), 'omitnan');
            stimulus_event_count = sum(event_impulses(analysis_index));

            if ~isfinite(stimulus_rate)
                continue
            end

            rate_common = interp1(time_selected, rate, ...
                common_time, 'linear', NaN);
            stimulus_is_nonzero = stimulus_event_count > 0;

            % Raw variants do not require a nonzero baseline.
            trial_variants{4}(trial_index, :) = rate_common;
            if stimulus_is_nonzero
                trial_variants{3}(trial_index, :) = rate_common;
            end

            % Percentage normalisation requires a finite nonzero baseline.
            if isfinite(baseline_rate) && ...
                    baseline_rate >= minimum_baseline_rate_hz
                rate_pct = 100 * (rate - baseline_rate) / baseline_rate;
                rate_pct_common = interp1(time_selected, rate_pct, ...
                    common_time, 'linear', NaN);
                trial_variants{2}(trial_index, :) = rate_pct_common;
                if stimulus_is_nonzero
                    trial_variants{1}(trial_index, :) = rate_pct_common;
                end
            end

            event_details{trial_index} = rmfield(details, 'Speed');
        end

        output_data.(['eventsVE_' condition_name]) = event_details;

        display_index = common_time >= display_period(1) & ...
            common_time <= display_period(2);
        summary_index = common_time >= analysis_period(1) & ...
            common_time <= analysis_period(2);
        display_samples = max(1, round(display_smooth_ms * fsample / 1000));
        this_display_time = common_time(display_index);
        if isempty(common_display_time)
            common_display_time = this_display_time;
        end

        for variant_index = 1:n_variants
            condition_curve = mean(trial_variants{variant_index}, 1, 'omitnan');
            output_structure = struct( ...
                'label', {{'MSRate'}}, ...
                'fsample', fsample, ...
                'time', common_time, ...
                'avg', condition_curve);
            output_field = ['msVE_' condition_name '_' ...
                variant_names{variant_index}];
            output_data.(output_field) = output_structure;

            display_curve = condition_curve(display_index);
            if display_samples > 1
                display_curve = smoothdata(display_curve, 2, ...
                    'gaussian', display_samples);
            end
            if numel(this_display_time) ~= numel(common_display_time) || ...
                    max(abs(this_display_time - common_display_time)) > 1e-9
                display_curve = interp1(this_display_time, display_curve, ...
                    common_display_time, 'linear', NaN);
            end

            subject_curves{subject_index, condition_index, variant_index} = ...
                display_curve;
            subject_summary(subject_index, condition_index, variant_index) = ...
                mean(condition_curve(summary_index), 'omitnan');
        end
    end

    output_data.detector_cfg = detector_cfg;
    output_data.rate_window_ms = rate_window_ms;
    output_data.minimum_baseline_rate_hz = minimum_baseline_rate_hz;
    save(output_file, '-struct', 'output_data');
end

if isempty(common_display_time)
    error('GCP:NoVanEdeData', 'No valid van Ede gaze-shift data were found.');
end

%% Assemble participant time courses
n_display_samples = numel(common_display_time);
subject_rate_matrix = nan(n_subjects, n_display_samples, ...
    n_conditions, n_variants);
for subject_index = 1:n_subjects
    for condition_index = 1:n_conditions
        for variant_index = 1:n_variants
            if ~isempty(subject_curves{subject_index, condition_index, variant_index})
                subject_rate_matrix(subject_index, :, condition_index, variant_index) = ...
                    subject_curves{subject_index, condition_index, variant_index};
            end
        end
    end
end

%% Combined figure: top row time courses, bottom row boxplots
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
tiled = tiledlayout(2, n_variants, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
axis_font = 16;

for variant_index = 1:n_variants
    variant_matrix = subject_rate_matrix(:, :, :, variant_index);
    grand_mean = squeeze(mean(variant_matrix, 1, 'omitnan'));
    n_valid = squeeze(sum(isfinite(variant_matrix), 1));
    grand_sem = squeeze(std(variant_matrix, 0, 1, 'omitnan')) ./ ...
        sqrt(max(n_valid, 1));
    grand_sem(n_valid < 2) = NaN;

    % Top row: grand-average time courses
    nexttile(variant_index);
    hold on
    for condition_index = 1:n_conditions
        error_bar = shadedErrorBar(common_display_time, ...
            grand_mean(:, condition_index), grand_sem(:, condition_index), ...
            'lineProps', {'-'}, 'transparent', true);
        set(error_bar.mainLine, 'Color', colors(condition_index, :), ...
            'LineWidth', 2);
        set(error_bar.patch, 'FaceColor', colors(condition_index, :), ...
            'FaceAlpha', 0.2);
        set(error_bar.edge, 'Color', 'none');
    end

    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
        'LineStyle', '--', 'HandleVisibility', 'off');
    yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
        'LineStyle', '--', 'HandleVisibility', 'off');
    xlim(display_period);
    xlabel('Time [s]');
    ylabel(variant_ylabels{variant_index});
    title(variant_titles{variant_index}, 'FontSize', axis_font);
    set(gca, 'FontSize', axis_font, 'Box', 'off');

    if variant_index == n_variants
        legend_handles = gobjects(n_conditions, 1);
        for condition_index = 1:n_conditions
            legend_handles(condition_index) = patch(nan, nan, ...
                colors(condition_index, :), 'FaceAlpha', 0.33, ...
                'EdgeColor', colors(condition_index, :), 'LineWidth', 1.5);
        end
        legend(legend_handles, condition_labels, 'Location', 'northeast', ...
            'FontSize', axis_font * 0.7, 'Box', 'off');
    end
    hold off

    % Bottom row: participant-level boxplots with raw points
    nexttile(n_variants + variant_index);
    hold on
    variant_summary = subject_summary(:, :, variant_index);
    for subject_index = 1:n_subjects
        valid_conditions = isfinite(variant_summary(subject_index, :));
        plot(find(valid_conditions), ...
            variant_summary(subject_index, valid_conditions), ...
            'Color', [0.8 0.8 0.8], 'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end

    boxplot(variant_summary, 'Labels', condition_labels, ...
        'Colors', 'k', 'Symbol', '', 'Widths', 0.45);
    for condition_index = 1:n_conditions
        values = variant_summary(:, condition_index);
        values = values(isfinite(values));
        if isempty(values)
            continue
        end
        jitter = linspace(-0.10, 0.10, numel(values))';
        scatter(condition_index + jitter, values, 30, ...
            colors(condition_index, :), 'filled', ...
            'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'none');
    end

    yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
        'LineStyle', '--', 'HandleVisibility', 'off');
    xlim([0.5 n_conditions + 0.5]);
    ylabel(variant_summary_ylabels{variant_index});
    set(gca, 'FontSize', axis_font, 'Box', 'off');
    hold off
end

title(tiled, 'van Ede microsaccade detection', ...
    'FontSize', axis_font * 1.3, 'FontWeight', 'bold');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(figpath, ...
    'GCP_gaze_microsaccades_vanede_overview.png'), '-dpng', '-r600');

%% Save participant summaries
summary_file = fullfile(datapath, ...
    'GCP_gaze_microsaccade_vanede_summary.mat');
save(summary_file, 'subjects', 'condition_names', 'condition_labels', ...
    'variant_names', 'subject_summary', 'subject_rate_matrix', ...
    'common_display_time', 'detector_cfg', 'rate_window_ms', ...
    'minimum_baseline_rate_hz');

fprintf('[VIZ GAZE MS VANEDE] Saved all four van Ede analysis variants.\n');
