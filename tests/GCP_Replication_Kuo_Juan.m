%% GCP Replication of Kuo & Juan (2026) Figures
% This script recreates all figures from:
% Kuo, C.-Y., & Juan, C.-H. (2026). Adaptive modulation of microsaccades 
% and saccade dynamics by global luminance. Frontiers in Systems Neuroscience, 19:1735778.
%
% Adapted for GCP project with contrast conditions: 25%, 50%, 75%, 100%
%
% Figures recreated:
% - Figure 2: Microsaccade rate, peak velocity, amplitude across time, main sequence
% - Figure 3: Averaged metrics during pre-target and target epochs
% - Figure 4: Effect of contrast and microsaccade occurrence on saccadic performance
% - Figure 5: Effect of target contrast on microsaccades during target epoch
% - Figure 6: Effect of target contrast and microsaccade occurrence on saccadic performance

%% Setup
clear; clc; close all;
startup;
[subjects, path, colors, ~] = setup('GCP');

% Data paths
data_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/';
gaze_path = fullfile(data_path, 'gaze_raw.mat');

% Time windows (relative to target onset at t=0)
pre_target_epoch = [-0.8, -0.2];  % -800 to -200 ms before target
target_epoch = [0.05, 0.14];      % 50 to 140 ms after target
peri_target_window = [-0.4, 0.1];  % -400 to 100 ms around target (for suppression analysis)
full_window = [-1.5, 2.0];         % Full analysis window

% Microsaccade detection parameters (Engbert & Kliegl, 2003)
ms_velocity_threshold = 6;  % median SDs
ms_amplitude_min = 0.1;     % degrees
ms_amplitude_max = 2.0;     % degrees
fsample = 500;              % Hz

% Saccade detection parameters
saccade_velocity_threshold = 30;  % deg/s
saccade_amplitude_min = 3;        % degrees

% Contrast conditions
contrast_levels = [25, 50, 75, 100];
n_contrasts = length(contrast_levels);

% Figure save directory
fig_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/tests/rep_kuo_juan/';
if ~exist(fig_save_dir, 'dir')
    mkdir(fig_save_dir);
end

%% Load raw gaze data (optional - for reference, but we'll use full dataET for analysis)
fprintf('Loading data...\n');
% Note: gaze_raw.mat contains windowed data [0.3, 2.0] s, but we need full epoch data
% So we'll extract directly from dataET instead

% Load preprocessed ET data for timing information
fprintf('Loading preprocessed ET data for timing...\n');
all_dataET = cell(length(subjects), 4);
for subj = 1:length(subjects)
    et_path = fullfile(data_path, subjects{subj}, 'gaze', 'dataET.mat');
    if exist(et_path, 'file')
        load(et_path);
        all_dataET{subj, 1} = dataET_c25;
        all_dataET{subj, 2} = dataET_c50;
        all_dataET{subj, 3} = dataET_c75;
        all_dataET{subj, 4} = dataET_c100;
    end
end

% Load behavioral data for saccade reaction times and accuracy
fprintf('Loading behavioral data...\n');
behav_path = fullfile(data_path, 'behavioral_matrix.mat');
if exist(behav_path, 'file')
    load(behav_path, 'behav_data');
else
    warning('Behavioral data not found. Some analyses will be skipped.');
    behav_data = [];
end

%% Extract microsaccade and saccade data
fprintf('Extracting microsaccade and saccade data...\n');

% Storage structures
all_ms_data = struct('subject', [], 'condition', [], 'trial', [], ...
    'pre_target_rate', [], 'pre_target_peak_vel', [], 'pre_target_amplitude', [], ...
    'target_rate', [], 'target_peak_vel', [], 'target_amplitude', [], ...
    'has_peri_target_ms', [], 'saccade_rt', [], 'saccade_peak_vel', [], ...
    'saccade_amplitude', [], 'saccade_endpoint_dev', [], 'saccade_accuracy', []);

ms_idx = 1;

for subj = 1:length(subjects)
    fprintf('  Processing subject %d/%d...\n', subj, length(subjects));
    
    for cond = 1:4
        % Get dataET for this condition
        if isempty(all_dataET{subj, cond})
            continue;
        end
        
        dataET = all_dataET{subj, cond};
        n_trials = length(dataET.trial);
        
        for trl = 1:n_trials
            % Extract full gaze data from dataET (not windowed)
            raw = dataET.trial{trl};
            tVec = dataET.time{trl};
            
            % Extract gaze coordinates
            x_raw = raw(1, :);
            y_raw = raw(2, :);
            
            % Flip Y axis to screen coordinates
            y_raw = 600 - y_raw;
            
            % Apply validity filter (remove points outside screen)
            valid = x_raw >= 0 & x_raw <= 800 & y_raw >= 0 & y_raw <= 600;
            x = x_raw(valid);
            y = y_raw(valid);
            tVec = tVec(valid);
            
            if length(x) < 10 || length(y) < 10
                continue;
            end
            
            % Ensure time vector matches data length
            if length(tVec) ~= length(x)
                tVec = (0:length(x)-1) / fsample + tVec(1);
            end
            
            gaze_data = [x; y];
            
            % Detect microsaccades in pre-target epoch
            [ms_rate_pre, ms_details_pre] = detect_microsaccades(fsample, gaze_data, tVec, pre_target_epoch, ...
                ms_velocity_threshold, ms_amplitude_min, ms_amplitude_max);
            
            % Calculate average peak velocity and amplitude for pre-target microsaccades
            if ~isempty(ms_details_pre)
                pre_peak_vel = mean([ms_details_pre.peak_velocity]);
                pre_amplitude = mean([ms_details_pre.amplitude]);
            else
                pre_peak_vel = NaN;
                pre_amplitude = NaN;
            end
            
            % Detect microsaccades in target epoch
            [ms_rate_target, ms_details_target] = detect_microsaccades(fsample, gaze_data, tVec, target_epoch, ...
                ms_velocity_threshold, ms_amplitude_min, ms_amplitude_max);
            
            % Calculate average peak velocity and amplitude for target microsaccades
            if ~isempty(ms_details_target)
                target_peak_vel = mean([ms_details_target.peak_velocity]);
                target_amplitude = mean([ms_details_target.amplitude]);
            else
                target_peak_vel = NaN;
                target_amplitude = NaN;
            end
            
            % Check for peri-target microsaccades (-400 to 100 ms)
            [~, ms_details_peri] = detect_microsaccades(fsample, gaze_data, tVec, peri_target_window, ...
                ms_velocity_threshold, ms_amplitude_min, ms_amplitude_max);
            has_peri_ms = ~isempty(ms_details_peri);
            
            % Extract saccade metrics (first saccade after target onset)
            % Calculate velocity
            if length(x) > 1
                px_to_deg = 0.05;
                vx = diff(x) * fsample * px_to_deg;
                vy = diff(y) * fsample * px_to_deg;
                v = sqrt(vx.^2 + vy.^2);
                
                % Find first saccade after target onset (velocity > 30 deg/s, amplitude > 3 deg)
                target_onset_idx = find(tVec >= 0, 1);
                if ~isempty(target_onset_idx) && target_onset_idx < length(v)
                    saccade_idx = find(v(target_onset_idx:end) > saccade_velocity_threshold, 1);
                    if ~isempty(saccade_idx)
                        saccade_start = target_onset_idx + saccade_idx - 1;
                        % Find saccade end
                        saccade_end = find(v(saccade_start:end) < saccade_velocity_threshold, 1);
                        if isempty(saccade_end)
                            saccade_end = length(v);
                        else
                            saccade_end = saccade_start + saccade_end - 1;
                        end
                        
                        % Calculate saccade metrics
                        saccade_rt = tVec(saccade_start);
                        saccade_peak_vel = max(v(saccade_start:saccade_end));
                        
                        dx_sacc = x(saccade_end+1) - x(saccade_start);
                        dy_sacc = y(saccade_end+1) - y(saccade_start);
                        saccade_amplitude = sqrt(dx_sacc^2 + dy_sacc^2) * px_to_deg;
                        
                        % Endpoint deviation (distance from target)
                        % Target is at 5 deg eccentricity (assume rightward for now)
                        target_x = 400 + 5/px_to_deg;  % Approximate pixel position
                        target_y = 300;
                        endpoint_x = x(saccade_end+1);
                        endpoint_y = y(saccade_end+1);
                        endpoint_dev = sqrt((endpoint_x - target_x)^2 + (endpoint_y - target_y)^2) * px_to_deg;
                        
                        % Accuracy (direction correct - simplified)
                        saccade_direction = atan2(dy_sacc, dx_sacc);
                        target_direction = 0;  % Assume rightward
                        saccade_accuracy = abs(saccade_direction - target_direction) < pi/4;  % Within 45 degrees
                    else
                        saccade_rt = NaN;
                        saccade_peak_vel = NaN;
                        saccade_amplitude = NaN;
                        endpoint_dev = NaN;
                        saccade_accuracy = NaN;
                    end
                else
                    saccade_rt = NaN;
                    saccade_peak_vel = NaN;
                    saccade_amplitude = NaN;
                    endpoint_dev = NaN;
                    saccade_accuracy = NaN;
                end
            else
                saccade_rt = NaN;
                saccade_peak_vel = NaN;
                saccade_amplitude = NaN;
                endpoint_dev = NaN;
                saccade_accuracy = NaN;
            end
            
            % Store data
            all_ms_data(ms_idx).subject = str2double(subjects{subj});
            all_ms_data(ms_idx).condition = cond;
            all_ms_data(ms_idx).trial = trl;
            all_ms_data(ms_idx).pre_target_rate = ms_rate_pre;
            all_ms_data(ms_idx).pre_target_peak_vel = pre_peak_vel;
            all_ms_data(ms_idx).pre_target_amplitude = pre_amplitude;
            all_ms_data(ms_idx).target_rate = ms_rate_target;
            all_ms_data(ms_idx).target_peak_vel = target_peak_vel;
            all_ms_data(ms_idx).target_amplitude = target_amplitude;
            all_ms_data(ms_idx).has_peri_target_ms = has_peri_ms;
            all_ms_data(ms_idx).saccade_rt = saccade_rt;
            all_ms_data(ms_idx).saccade_peak_vel = saccade_peak_vel;
            all_ms_data(ms_idx).saccade_amplitude = saccade_amplitude;
            all_ms_data(ms_idx).saccade_endpoint_dev = endpoint_dev;
            all_ms_data(ms_idx).saccade_accuracy = saccade_accuracy;
            
            ms_idx = ms_idx + 1;
        end
    end
end

fprintf('Data extraction complete.\n');

%% FIGURE 2: Microsaccade metrics across time and main sequence
fprintf('Creating Figure 2...\n');
fprintf('Computing time-resolved microsaccade metrics from raw data...\n');

% Time vector for time-resolved analysis
time_vec = -1.0:0.01:0.5;  % 10ms resolution
moving_window = 0.1;  % 100ms moving window as per paper
px_to_deg = 0.05;  % Pixel to degree conversion

% Storage for time-resolved data
ms_rate_time = nan(n_contrasts, length(time_vec));
ms_peak_vel_time = nan(n_contrasts, length(time_vec));
ms_amplitude_time = nan(n_contrasts, length(time_vec));

% Collect all microsaccades with timestamps for main sequence
all_ms_for_mainseq = cell(n_contrasts, 1);

for cond = 1:n_contrasts
    fprintf('  Processing contrast %d%%...\n', contrast_levels(cond));
    
    % Collect all microsaccades from all trials for this condition
    all_ms_times = [];
    all_ms_peak_vels = [];
    all_ms_amplitudes = [];
    
    % Storage for accumulating time-resolved data across trials
    ms_count_accum = zeros(length(time_vec), 1);  % Total count of microsaccades in window
    ms_peak_vel_sum = zeros(length(time_vec), 1);  % Sum of peak velocities
    ms_amplitude_sum = zeros(length(time_vec), 1);  % Sum of amplitudes
    trial_count_accum = zeros(length(time_vec), 1);  % Number of trials contributing to each time point
    
    % Process each subject
    for subj = 1:length(subjects)
        % Get dataET for this condition
        if isempty(all_dataET{subj, cond})
            continue;
        end
        
        dataET = all_dataET{subj, cond};
        n_trials = length(dataET.trial);
        
        % Process each trial
        for trl = 1:n_trials
            % Extract FULL gaze data from dataET (not windowed)
            raw = dataET.trial{trl};
            tVec = dataET.time{trl};
            
            % Extract gaze coordinates
            x_raw = raw(1, :);
            y_raw = raw(2, :);
            
            % Flip Y axis to screen coordinates
            y_raw = 600 - y_raw;
            
            % Apply validity filter (remove points outside screen)
            valid = x_raw >= 0 & x_raw <= 800 & y_raw >= 0 & y_raw <= 600;
            x = x_raw(valid);
            y = y_raw(valid);
            tVec = tVec(valid);
            
            if length(x) < 10 || length(y) < 10
                continue;
            end
            
            % Ensure time vector matches data length
            if length(tVec) ~= length(x)
                tVec = (0:length(x)-1) / fsample + tVec(1);
            end
            
            % Detect all microsaccades in this trial
            gaze_data = [x; y];
            [~, ms_details] = detect_microsaccades(fsample, gaze_data, tVec, [-inf, inf], ...
                ms_velocity_threshold, ms_amplitude_min, ms_amplitude_max);
            
            % Store microsaccade details for main sequence
            if ~isempty(ms_details)
                for ms = 1:length(ms_details)
                    all_ms_times(end+1) = ms_details(ms).onset;
                    all_ms_peak_vels(end+1) = ms_details(ms).peak_velocity;
                    all_ms_amplitudes(end+1) = ms_details(ms).amplitude;
                end
            end
            
            % Compute time-resolved metrics using moving window for this trial
            for t_idx = 1:length(time_vec)
                t_center = time_vec(t_idx);
                window_start = t_center - moving_window/2;
                window_end = t_center + moving_window/2;
                
                % Find microsaccades in this window
                ms_in_window = [];
                if ~isempty(ms_details)
                    for ms = 1:length(ms_details)
                        ms_time = ms_details(ms).onset;
                        if ms_time >= window_start && ms_time <= window_end
                            ms_in_window(end+1) = ms;
                        end
                    end
                end
                
                % Always increment trial count (even if no microsaccades)
                trial_count_accum(t_idx) = trial_count_accum(t_idx) + 1;
                
                if ~isempty(ms_in_window)
                    % Accumulate count of microsaccades
                    ms_count_accum(t_idx) = ms_count_accum(t_idx) + length(ms_in_window);
                    
                    % Accumulate peak velocity (sum, not mean)
                    peak_vels_window = [ms_details(ms_in_window).peak_velocity];
                    ms_peak_vel_sum(t_idx) = ms_peak_vel_sum(t_idx) + sum(peak_vels_window);
                    
                    % Accumulate amplitude (sum, not mean)
                    amps_window = [ms_details(ms_in_window).amplitude];
                    ms_amplitude_sum(t_idx) = ms_amplitude_sum(t_idx) + sum(amps_window);
                end
            end
        end
    end
    
    % Average across all trials
    for t_idx = 1:length(time_vec)
        if trial_count_accum(t_idx) > 0
            % Rate = total microsaccades / (window_duration * number_of_trials)
            ms_rate_time(cond, t_idx) = ms_count_accum(t_idx) / (moving_window * trial_count_accum(t_idx));
            
            % Average peak velocity and amplitude (divide sum by count of microsaccades)
            if ms_count_accum(t_idx) > 0
                ms_peak_vel_time(cond, t_idx) = ms_peak_vel_sum(t_idx) / ms_count_accum(t_idx);
                ms_amplitude_time(cond, t_idx) = ms_amplitude_sum(t_idx) / ms_count_accum(t_idx);
            end
        end
    end
    
    % Store for main sequence plot
    all_ms_for_mainseq{cond} = struct('amplitudes', all_ms_amplitudes, 'peak_vels', all_ms_peak_vels);
end

% Smooth the time-resolved data with moving average (100ms window)
smooth_win = round(moving_window / 0.01);  % Convert to samples
for cond = 1:n_contrasts
    ms_rate_time(cond, :) = movmean(ms_rate_time(cond, :), smooth_win, 'omitnan');
    ms_peak_vel_time(cond, :) = movmean(ms_peak_vel_time(cond, :), smooth_win, 'omitnan');
    ms_amplitude_time(cond, :) = movmean(ms_amplitude_time(cond, :), smooth_win, 'omitnan');
end

% Main sequence (peak velocity vs amplitude)
figure('Position', [0, 0, 1512, 982], 'Color', 'w');
sgtitle('Figure 2: Effect of Contrast on Microsaccadic Metrics', 'FontSize', 16, 'FontWeight', 'bold');

% Panel A: Microsaccade rate across time
subplot(2, 2, 1);
hold on;
for cond = 1:n_contrasts
    plot(time_vec, ms_rate_time(cond, :), 'LineWidth', 2, 'Color', colors(cond, :));
end
xline(0, 'k--', 'LineWidth', 1);
xlabel('Time relative to target (s)', 'FontSize', 12);
ylabel('Microsaccade rate (Hz)', 'FontSize', 12);
title('A: Microsaccade rate across time', 'FontSize', 12);
legend({'25%', '50%', '75%', '100%'}, 'Location', 'best');
grid on;
xlim([-1.0, 0.5]);

% Panel B: Peak velocity across time
subplot(2, 2, 2);
hold on;
for cond = 1:n_contrasts
    plot(time_vec, ms_peak_vel_time(cond, :), 'LineWidth', 2, 'Color', colors(cond, :));
end
xline(0, 'k--', 'LineWidth', 1);
xlabel('Time relative to target (s)', 'FontSize', 12);
ylabel('Peak velocity (deg/s)', 'FontSize', 12);
title('B: Peak velocity across time', 'FontSize', 12);
legend({'25%', '50%', '75%', '100%'}, 'Location', 'best');
grid on;
xlim([-1.0, 0.5]);

% Panel C: Amplitude across time
subplot(2, 2, 3);
hold on;
for cond = 1:n_contrasts
    plot(time_vec, ms_amplitude_time(cond, :), 'LineWidth', 2, 'Color', colors(cond, :));
end
xline(0, 'k--', 'LineWidth', 1);
xlabel('Time relative to target (s)', 'FontSize', 12);
ylabel('Amplitude (deg)', 'FontSize', 12);
title('C: Amplitude across time', 'FontSize', 12);
legend({'25%', '50%', '75%', '100%'}, 'Location', 'best');
grid on;
xlim([-1.0, 0.5]);

% Panel D: Main sequence (peak velocity vs amplitude)
subplot(2, 2, 4);
hold on;
for cond = 1:n_contrasts
    if ~isempty(all_ms_for_mainseq{cond}.amplitudes)
        amplitudes = all_ms_for_mainseq{cond}.amplitudes;
        peak_vels = all_ms_for_mainseq{cond}.peak_vels;
        valid_idx = ~isnan(amplitudes) & ~isnan(peak_vels) & amplitudes > 0 & peak_vels > 0;
        if sum(valid_idx) > 0
            h = scatter(amplitudes(valid_idx), peak_vels(valid_idx), 50, colors(cond, :), 'filled');
            if isprop(h, 'MarkerFaceAlpha')
                h.MarkerFaceAlpha = 0.6;
            end
        end
    end
end
xlabel('Amplitude (deg)', 'FontSize', 12);
ylabel('Peak velocity (deg/s)', 'FontSize', 12);
title('D: Main sequence', 'FontSize', 12);
legend({'25%', '50%', '75%', '100%'}, 'Location', 'best');
grid on;

saveas(gcf, fullfile(fig_save_dir, 'Figure2_Microsaccade_Metrics.png'));

%% FIGURE 3: Averaged metrics during pre-target and target epochs
fprintf('Creating Figure 3...\n');

figure('Position', [0, 0, 1400, 800], 'Color', 'w');
sgtitle('Figure 3: Effect of Contrast on Microsaccadic Metrics', 'FontSize', 16, 'FontWeight', 'bold');

% Pre-target epoch averages
pre_target_rate = nan(n_contrasts, 1);
pre_target_peak_vel = nan(n_contrasts, 1);
pre_target_amplitude = nan(n_contrasts, 1);

% Target epoch averages
target_rate = nan(n_contrasts, 1);
target_peak_vel = nan(n_contrasts, 1);
target_amplitude = nan(n_contrasts, 1);

for cond = 1:n_contrasts
    cond_data = all_ms_data([all_ms_data.condition] == cond);
    
    pre_target_rate(cond) = mean([cond_data.pre_target_rate], 'omitnan');
    pre_target_peak_vel(cond) = mean([cond_data.pre_target_peak_vel], 'omitnan');
    pre_target_amplitude(cond) = mean([cond_data.pre_target_amplitude], 'omitnan');
    
    target_rate(cond) = mean([cond_data.target_rate], 'omitnan');
    target_peak_vel(cond) = mean([cond_data.target_peak_vel], 'omitnan');
    target_amplitude(cond) = mean([cond_data.target_amplitude], 'omitnan');
end

% Pre-target epoch
subplot(2, 3, 1);
bar(pre_target_rate, 'FaceColor', 'flat', 'CData', colors(1:4, :));
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Rate (Hz)', 'FontSize', 12);
title('Pre-target: Rate', 'FontSize', 12);
grid on;

subplot(2, 3, 2);
bar(pre_target_peak_vel, 'FaceColor', 'flat', 'CData', colors(1:4, :));
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Peak velocity (deg/s)', 'FontSize', 12);
title('Pre-target: Peak velocity', 'FontSize', 12);
grid on;

subplot(2, 3, 3);
bar(pre_target_amplitude, 'FaceColor', 'flat', 'CData', colors(1:4, :));
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Amplitude (deg)', 'FontSize', 12);
title('Pre-target: Amplitude', 'FontSize', 12);
grid on;

% Target epoch
subplot(2, 3, 4);
bar(target_rate, 'FaceColor', 'flat', 'CData', colors(1:4, :));
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Rate (Hz)', 'FontSize', 12);
title('Target: Rate', 'FontSize', 12);
grid on;

subplot(2, 3, 5);
bar(target_peak_vel, 'FaceColor', 'flat', 'CData', colors(1:4, :));
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Peak velocity (deg/s)', 'FontSize', 12);
title('Target: Peak velocity', 'FontSize', 12);
grid on;

subplot(2, 3, 6);
bar(target_amplitude, 'FaceColor', 'flat', 'CData', colors(1:4, :));
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Amplitude (deg)', 'FontSize', 12);
title('Target: Amplitude', 'FontSize', 12);
grid on;

saveas(gcf, fullfile(fig_save_dir, 'Figure3_Epoch_Averages.png'));

%% FIGURE 4: Effect of contrast and microsaccade occurrence on saccadic performance
fprintf('Creating Figure 4...\n');

figure('Position', [0, 0, 1512, 982], 'Color', 'w');
sgtitle('Figure 4: Effect of Contrast and Microsaccade Occurrence on Saccadic Performance', 'FontSize', 16, 'FontWeight', 'bold');

% Separate trials with and without peri-target microsaccades
for cond = 1:n_contrasts
    cond_data = all_ms_data([all_ms_data.condition] == cond);
    
    % Trials without microsaccades
    no_ms_trials = cond_data([cond_data.has_peri_target_ms] == 0);
    % Trials with microsaccades
    ms_trials = cond_data([cond_data.has_peri_target_ms] == 1);
    
    % Calculate means
    accuracy_no_ms(cond) = mean([no_ms_trials.saccade_accuracy], 'omitnan') * 100;
    accuracy_ms(cond) = mean([ms_trials.saccade_accuracy], 'omitnan') * 100;
    
    rt_no_ms(cond) = mean([no_ms_trials.saccade_rt], 'omitnan') * 1000;  % Convert to ms
    rt_ms(cond) = mean([ms_trials.saccade_rt], 'omitnan') * 1000;
    
    peak_vel_no_ms(cond) = mean([no_ms_trials.saccade_peak_vel], 'omitnan');
    peak_vel_ms(cond) = mean([ms_trials.saccade_peak_vel], 'omitnan');
    
    amplitude_no_ms(cond) = mean([no_ms_trials.saccade_amplitude], 'omitnan');
    amplitude_ms(cond) = mean([ms_trials.saccade_amplitude], 'omitnan');
    
    endpoint_dev_no_ms(cond) = mean([no_ms_trials.saccade_endpoint_dev], 'omitnan');
    endpoint_dev_ms(cond) = mean([ms_trials.saccade_endpoint_dev], 'omitnan');
end

% Panel A: Accuracy
subplot(2, 3, 1);
x = 1:n_contrasts;
bar(x - 0.2, accuracy_no_ms, 0.4, 'FaceColor', [0.7, 0.7, 0.7], 'DisplayName', 'No MS');
hold on;
bar(x + 0.2, accuracy_ms, 0.4, 'FaceColor', [0.3, 0.3, 0.3], 'DisplayName', 'With MS');
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Accuracy (%)', 'FontSize', 12);
title('A: Saccadic directional accuracy', 'FontSize', 12);
legend('Location', 'best');
grid on;

% Panel B: Reaction time
subplot(2, 3, 2);
bar(x - 0.2, rt_no_ms, 0.4, 'FaceColor', [0.7, 0.7, 0.7], 'DisplayName', 'No MS');
hold on;
bar(x + 0.2, rt_ms, 0.4, 'FaceColor', [0.3, 0.3, 0.3], 'DisplayName', 'With MS');
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('SRT (ms)', 'FontSize', 12);
title('B: Saccade reaction time', 'FontSize', 12);
legend('Location', 'best');
grid on;

% Panel C: Peak velocity
subplot(2, 3, 3);
bar(x - 0.2, peak_vel_no_ms, 0.4, 'FaceColor', [0.7, 0.7, 0.7], 'DisplayName', 'No MS');
hold on;
bar(x + 0.2, peak_vel_ms, 0.4, 'FaceColor', [0.3, 0.3, 0.3], 'DisplayName', 'With MS');
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Peak velocity (deg/s)', 'FontSize', 12);
title('C: Peak velocity', 'FontSize', 12);
legend('Location', 'best');
grid on;

% Panel D: Amplitude
subplot(2, 3, 4);
bar(x - 0.2, amplitude_no_ms, 0.4, 'FaceColor', [0.7, 0.7, 0.7], 'DisplayName', 'No MS');
hold on;
bar(x + 0.2, amplitude_ms, 0.4, 'FaceColor', [0.3, 0.3, 0.3], 'DisplayName', 'With MS');
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Amplitude (deg)', 'FontSize', 12);
title('D: Amplitude', 'FontSize', 12);
legend('Location', 'best');
grid on;

% Panel E: Endpoint deviation
subplot(2, 3, 5);
bar(x - 0.2, endpoint_dev_no_ms, 0.4, 'FaceColor', [0.7, 0.7, 0.7], 'DisplayName', 'No MS');
hold on;
bar(x + 0.2, endpoint_dev_ms, 0.4, 'FaceColor', [0.3, 0.3, 0.3], 'DisplayName', 'With MS');
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Endpoint deviation (deg)', 'FontSize', 12);
title('E: Endpoint deviation', 'FontSize', 12);
legend('Location', 'best');
grid on;

saveas(gcf, fullfile(fig_save_dir, 'Figure4_Contrast_MS_Suppression.png'));

%% FIGURE 5: Effect of target contrast on microsaccades during target epoch
fprintf('Creating Figure 5...\n');

figure('Position', [0, 0, 1400, 800], 'Color', 'w');
sgtitle('Figure 5: Effect of Target Contrast on Microsaccades During Target Epoch', 'FontSize', 16, 'FontWeight', 'bold');

% Panel A: Microsaccade rate across time (target epoch highlighted)
subplot(2, 3, 1);
hold on;
for cond = 1:n_contrasts
    plot(time_vec, ms_rate_time(cond, :), 'LineWidth', 2, 'Color', colors(cond, :));
end
xline(0, 'k--', 'LineWidth', 1);
xline(target_epoch(1), 'r--', 'LineWidth', 1);
xline(target_epoch(2), 'r--', 'LineWidth', 1);
max_rate = max(ms_rate_time(:));
if ~isnan(max_rate) && max_rate > 0
    patch([target_epoch(1), target_epoch(2), target_epoch(2), target_epoch(1)], ...
          [0, 0, max_rate, max_rate], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
xlabel('Time relative to target (s)', 'FontSize', 12);
ylabel('Microsaccade rate (Hz)', 'FontSize', 12);
title('A: Microsaccade rate', 'FontSize', 12);
legend({'25%', '50%', '75%', '100%'}, 'Location', 'best');
grid on;
xlim([-0.2, 0.3]);

% Panel B: Averaged rate during target epoch
subplot(2, 3, 2);
bar(target_rate, 'FaceColor', 'flat', 'CData', colors(1:4, :));
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Rate (Hz)', 'FontSize', 12);
title('B: Averaged rate (target epoch)', 'FontSize', 12);
grid on;

% Panel C: Peak velocity across time
subplot(2, 3, 3);
hold on;
for cond = 1:n_contrasts
    plot(time_vec, ms_peak_vel_time(cond, :), 'LineWidth', 2, 'Color', colors(cond, :));
end
xline(0, 'k--', 'LineWidth', 1);
xline(target_epoch(1), 'r--', 'LineWidth', 1);
xline(target_epoch(2), 'r--', 'LineWidth', 1);
max_vel = max(ms_peak_vel_time(:));
if ~isnan(max_vel) && max_vel > 0
    patch([target_epoch(1), target_epoch(2), target_epoch(2), target_epoch(1)], ...
          [0, 0, max_vel, max_vel], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
xlabel('Time relative to target (s)', 'FontSize', 12);
ylabel('Peak velocity (deg/s)', 'FontSize', 12);
title('C: Peak velocity', 'FontSize', 12);
legend({'25%', '50%', '75%', '100%'}, 'Location', 'best');
grid on;
xlim([-0.2, 0.3]);

% Panel D: Averaged peak velocity during target epoch
subplot(2, 3, 4);
bar(target_peak_vel, 'FaceColor', 'flat', 'CData', colors(1:4, :));
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Peak velocity (deg/s)', 'FontSize', 12);
title('D: Averaged peak velocity (target epoch)', 'FontSize', 12);
grid on;

% Panel E: Amplitude across time
subplot(2, 3, 5);
hold on;
for cond = 1:n_contrasts
    plot(time_vec, ms_amplitude_time(cond, :), 'LineWidth', 2, 'Color', colors(cond, :));
end
xline(0, 'k--', 'LineWidth', 1);
xline(target_epoch(1), 'r--', 'LineWidth', 1);
xline(target_epoch(2), 'r--', 'LineWidth', 1);
max_amp = max(ms_amplitude_time(:));
if ~isnan(max_amp) && max_amp > 0
    patch([target_epoch(1), target_epoch(2), target_epoch(2), target_epoch(1)], ...
          [0, 0, max_amp, max_amp], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
xlabel('Time relative to target (s)', 'FontSize', 12);
ylabel('Amplitude (deg)', 'FontSize', 12);
title('E: Amplitude', 'FontSize', 12);
legend({'25%', '50%', '75%', '100%'}, 'Location', 'best');
grid on;
xlim([-0.2, 0.3]);

% Panel F: Averaged amplitude during target epoch
subplot(2, 3, 6);
bar(target_amplitude, 'FaceColor', 'flat', 'CData', colors(1:4, :));
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Amplitude (deg)', 'FontSize', 12);
title('F: Averaged amplitude (target epoch)', 'FontSize', 12);
grid on;

saveas(gcf, fullfile(fig_save_dir, 'Figure5_Contrast_Target_Epoch.png'));

%% FIGURE 6: Effect of target contrast and microsaccade occurrence on saccadic performance
fprintf('Creating Figure 6...\n');

figure('Position', [0, 0, 1512, 982], 'Color', 'w');
sgtitle('Figure 6: Effect of Target Contrast and Microsaccade Occurrence on Saccadic Performance', 'FontSize', 16, 'FontWeight', 'bold');

% Define x positions for bars
x = 1:n_contrasts;

% Panel A: Accuracy
subplot(2, 3, 1);
bar(x - 0.2, accuracy_no_ms, 0.4, 'FaceColor', [0.7, 0.7, 0.7], 'DisplayName', 'No MS');
hold on;
bar(x + 0.2, accuracy_ms, 0.4, 'FaceColor', [0.3, 0.3, 0.3], 'DisplayName', 'With MS');
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Accuracy (%)', 'FontSize', 12);
title('A: Saccadic directional accuracy', 'FontSize', 12);
legend('Location', 'best');
grid on;

% Panel B: Reaction time
subplot(2, 3, 2);
bar(x - 0.2, rt_no_ms, 0.4, 'FaceColor', [0.7, 0.7, 0.7], 'DisplayName', 'No MS');
hold on;
bar(x + 0.2, rt_ms, 0.4, 'FaceColor', [0.3, 0.3, 0.3], 'DisplayName', 'With MS');
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('SRT (ms)', 'FontSize', 12);
title('B: Saccade reaction time', 'FontSize', 12);
legend('Location', 'best');
grid on;

% Panel C: Peak velocity
subplot(2, 3, 3);
bar(x - 0.2, peak_vel_no_ms, 0.4, 'FaceColor', [0.7, 0.7, 0.7], 'DisplayName', 'No MS');
hold on;
bar(x + 0.2, peak_vel_ms, 0.4, 'FaceColor', [0.3, 0.3, 0.3], 'DisplayName', 'With MS');
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Peak velocity (deg/s)', 'FontSize', 12);
title('C: Peak velocity', 'FontSize', 12);
legend('Location', 'best');
grid on;

% Panel D: Amplitude
subplot(2, 3, 4);
bar(x - 0.2, amplitude_no_ms, 0.4, 'FaceColor', [0.7, 0.7, 0.7], 'DisplayName', 'No MS');
hold on;
bar(x + 0.2, amplitude_ms, 0.4, 'FaceColor', [0.3, 0.3, 0.3], 'DisplayName', 'With MS');
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Amplitude (deg)', 'FontSize', 12);
title('D: Amplitude', 'FontSize', 12);
legend('Location', 'best');
grid on;

% Panel E: Endpoint deviation
subplot(2, 3, 5);
bar(x - 0.2, endpoint_dev_no_ms, 0.4, 'FaceColor', [0.7, 0.7, 0.7], 'DisplayName', 'No MS');
hold on;
bar(x + 0.2, endpoint_dev_ms, 0.4, 'FaceColor', [0.3, 0.3, 0.3], 'DisplayName', 'With MS');
set(gca, 'XTickLabel', {'25%', '50%', '75%', '100%'});
ylabel('Endpoint deviation (deg)', 'FontSize', 12);
title('E: Endpoint deviation', 'FontSize', 12);
legend('Location', 'best');
grid on;

saveas(gcf, fullfile(fig_save_dir, 'Figure6_Contrast_MS_Saccadic_Performance.png'));

fprintf('All figures created successfully!\n');
fprintf('Figures saved in: %s\n', fig_save_dir);

%% Local Functions

function [ms_rate, ms_details] = detect_microsaccades(fsample, gaze_data, tVec, time_window, velocity_threshold, amplitude_min, amplitude_max)
    % Detect microsaccades using Engbert & Kliegl (2003) method
    % gaze_data: 2 x N matrix (x, y positions)
    % tVec: time vector
    % time_window: [t_start, t_end] in seconds
    % velocity_threshold: number of median SDs for threshold
    % amplitude_min: minimum amplitude in degrees
    % amplitude_max: maximum amplitude in degrees
    
    if nargin < 4
        time_window = [-inf, inf];
    end
    if nargin < 5
        velocity_threshold = 6;
    end
    if nargin < 6
        amplitude_min = 0.1;
    end
    if nargin < 7
        amplitude_max = 2.0;
    end
    
    % Ensure tVec and gaze_data have matching lengths
    if length(tVec) ~= size(gaze_data, 2)
        % If mismatch, use the shorter length
        min_len = min(length(tVec), size(gaze_data, 2));
        tVec = tVec(1:min_len);
        gaze_data = gaze_data(:, 1:min_len);
    end
    
    % Extract data within time window
    time_idx = tVec >= time_window(1) & tVec <= time_window(2);
    if sum(time_idx) < 10
        ms_rate = 0;
        ms_details = [];
        return;
    end
    
    x = gaze_data(1, time_idx);
    y = gaze_data(2, time_idx);
    t = tVec(time_idx);
    
    % Convert pixels to degrees (assuming screen distance ~70cm, screen size ~43x24 deg)
    % Approximate conversion: 1 pixel ≈ 0.05 degrees (adjust based on your setup)
    px_to_deg = 0.05;
    
    % Calculate velocity
    if length(x) < 2
        ms_rate = 0;
        ms_details = [];
        return;
    end
    
    vx = diff(x) * fsample * px_to_deg;  % deg/s
    vy = diff(y) * fsample * px_to_deg;  % deg/s
    v = sqrt(vx.^2 + vy.^2);
    
    % Velocity threshold (6 median SDs)
    med_v = median(v);
    mad_v = mad(v, 1);  % median absolute deviation
    threshold = med_v + velocity_threshold * mad_v;
    
    % Find microsaccades
    ms_idx = v > threshold;
    
    % Group consecutive samples into microsaccades
    ms_onsets = [];
    ms_offsets = [];
    in_ms = false;
    
    for i = 1:length(ms_idx)
        if ms_idx(i) && ~in_ms
            ms_onsets(end+1) = i;
            in_ms = true;
        elseif ~ms_idx(i) && in_ms
            ms_offsets(end+1) = i-1;
            in_ms = false;
        end
    end
    
    if in_ms
        ms_offsets(end+1) = length(ms_idx);
    end
    
    % Calculate microsaccade metrics
    ms_details = [];
    ms_count = 0;
    
    for i = 1:length(ms_onsets)
        onset_idx = ms_onsets(i);
        offset_idx = ms_offsets(i);
        
        % Amplitude
        dx = x(offset_idx+1) - x(onset_idx);
        dy = y(offset_idx+1) - y(onset_idx);
        amplitude = sqrt(dx^2 + dy^2) * px_to_deg;
        
        % Check amplitude criteria
        if amplitude >= amplitude_min && amplitude <= amplitude_max
            ms_count = ms_count + 1;
            ms_details(ms_count).onset = t(onset_idx);
            ms_details(ms_count).offset = t(offset_idx+1);
            ms_details(ms_count).amplitude = amplitude;
            ms_details(ms_count).peak_velocity = max(v(onset_idx:offset_idx));
            ms_details(ms_count).duration = (offset_idx - onset_idx + 1) / fsample;
        end
    end
    
    % Calculate rate (per second)
    window_duration = time_window(2) - time_window(1);
    if window_duration > 0
        ms_rate = ms_count / window_duration;
    else
        ms_rate = 0;
    end
end
