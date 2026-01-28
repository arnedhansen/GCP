%% GCP Generalized Eigendecomposition (GED) for Gamma Band Analysis
%
% This script uses GED to identify gamma oscillations by comparing
% baseline (pre-stimulus) vs stimulus (post-stimulus) periods.
% GED finds the spatial filter that maximizes the signal-to-noise ratio
% for gamma activity, helping to identify peak gamma frequency and power
% more reliably than single-subject power spectra.
%
% Implementation follows Cohen (2022) NeuroImage tutorial:
% "A tutorial on generalized eigendecomposition for denoising, contrast 
%  enhancement, and dimension reduction in multichannel electrophysiology"
%
% Input: Preprocessed and FOOOFed data (dataEEG)
% Output: GED components, power spectra, peak frequencies/powers, and visualizations
%
% For each condition separately:
%   1. Filter data to gamma band (30-90 Hz)
%   2. Compute baseline and stimulus covariance matrices (per trial, then average)
%   3. Perform GED (eig(covStim, covBase)) with shrinkage regularization
%   4. Extract top component(s) and handle sign uncertainty
%   5. Compute power spectrum of GED component time series
%   6. Identify peak gamma frequency and power
%   7. Visualize results

clear; close all; clc

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');

% Limit to first 10 pilot participants
nPilot = min(10, length(subjects));
subjects_pilot = subjects(1:nPilot);

% Time windows
baseline_window = [-1.5, -0.25];  % Pre-stimulus baseline
stimulus_window = [0, 2.0];        % Post-stimulus period

% Gamma frequency band
gamma_freq = [30, 90];  % Hz

% GED parameters
nComponents = 5;  % Extract top 5 components (but focus on top 1)

% Condition names and trial codes
condNames = {'c25', 'c50', 'c75', 'c100'};
trialCodes = [61, 62, 63, 64];
condLabels = {'25% Contrast', '50% Contrast', '75% Contrast', '100% Contrast'};

% Figure save directory
fig_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/tests/ged';
if ~exist(fig_save_dir, 'dir')
    mkdir(fig_save_dir);
end

% Load head model for topoplots
if ispc
    load('W:\Students\Arne\MA\headmodel\ant128lay.mat');
else
    load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/ant128lay.mat');
end

%% Preallocate storage for grand average
all_topos = cell(4, nPilot);  % [condition, subject]
all_powerspec_stim = cell(4, nPilot);  % Raw stimulus power
all_powerspec_base = cell(4, nPilot);  % Baseline power
all_powerspec_diff = cell(4, nPilot);  % Difference (stimulus - baseline)
all_powerspec_rel = cell(4, nPilot);   % Relative change ((stimulus - baseline) / baseline)
all_freqs = cell(4, nPilot);
all_peak_freqs = nan(4, nPilot);
all_peak_powers = nan(4, nPilot);
all_eigenvals = cell(4, nPilot);
all_component_ts_stim = cell(4, nPilot);
all_component_ts_base = cell(4, nPilot);
all_times_stim = cell(4, nPilot);
all_times_base = cell(4, nPilot);
chanlocs_all = {};  % Store channel labels

%% Process each subject
for subj = 1:nPilot
    
    fprintf('\n========================================\n');
    fprintf('Processing Subject %d/%d: %s\n', subj, nPilot, subjects_pilot{subj});
    fprintf('========================================\n');
    
    datapath = strcat(path, subjects_pilot{subj}, filesep, 'eeg');
    cd(datapath)
    
    % Load preprocessed data
    load dataEEG
    
    % Get sampling rate
    fsample = dataEEG_c25.fsample;
    
    % Identify trial indices for each condition
    ind61 = find(dataEEG_c25.trialinfo == 61);
    ind62 = find(dataEEG_c50.trialinfo == 62);
    ind63 = find(dataEEG_c75.trialinfo == 63);
    ind64 = find(dataEEG_c100.trialinfo == 64);
    
    trialIndices = {ind61, ind62, ind63, ind64};
    dataStructs = {dataEEG_c25, dataEEG_c50, dataEEG_c75, dataEEG_c100};
    
    %% Process each condition
    for cond = 1:4
        
        fprintf('  Condition %d/4: %s\n', cond, condLabels{cond});
        
        dat = dataStructs{cond};
        trlIdx = trialIndices{cond};
        
        if isempty(trlIdx)
            warning('No trials found for condition %d, subject %s', cond, subjects_pilot{subj});
            continue;
        end
        
        % Select trials
        cfg = [];
        cfg.trials = trlIdx;
        dat = ft_selectdata(cfg, dat);
        
        % Filter to gamma band
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = gamma_freq;
        cfg.bpfilttype = 'fir';
        cfg.bpfiltord = 4*round(fsample/gamma_freq(1));  % Filter order
        dat_gamma = ft_preprocessing(cfg, dat);
        
        % Extract baseline and stimulus periods
        cfg_base = [];
        cfg_base.latency = baseline_window;
        dat_base = ft_selectdata(cfg_base, dat_gamma);
        
        cfg_stim = [];
        cfg_stim.latency = stimulus_window;
        dat_stim = ft_selectdata(cfg_stim, dat_gamma);
        
        % Compute covariance matrices per trial, then average (as recommended in Cohen 2022)
        % This approach is more robust than concatenating all trials
        nTrials_base = length(dat_base.trial);
        nTrials_stim = length(dat_stim.trial);
        
        % Initialize covariance matrices
        nChans = length(dat_base.label);
        covBase_trials = zeros(nChans, nChans, nTrials_base);
        covStim_trials = zeros(nChans, nChans, nTrials_stim);
        
        % Compute covariance for each baseline trial
        for trl = 1:nTrials_base
            trial_data = double(dat_base.trial{trl});
            % Mean-center each trial separately (critical: mean-center per trial)
            trial_data = bsxfun(@minus, trial_data, mean(trial_data, 2));
            nTimePnts = size(trial_data, 2);
            covBase_trials(:, :, trl) = (trial_data * trial_data') / nTimePnts;
        end
        
        % Compute covariance for each stimulus trial
        for trl = 1:nTrials_stim
            trial_data = double(dat_stim.trial{trl});
            % Mean-center each trial separately (critical: mean-center per trial)
            trial_data = bsxfun(@minus, trial_data, mean(trial_data, 2));
            nTimePnts = size(trial_data, 2);
            covStim_trials(:, :, trl) = (trial_data * trial_data') / nTimePnts;
        end
        
        % Average covariance matrices across trials
        covBase = mean(covBase_trials, 3);
        covStim = mean(covStim_trials, 3);
        
        % Apply shrinkage regularization (Cohen 2022, Section 3.4)
        % R_reg = (1-λ)*R + λ*mean(diag(R))*I
        % This preserves the trace of the matrix
        lambda = 0.01;  % Regularization parameter
        meanDiagBase = mean(diag(covBase));
        meanDiagStim = mean(diag(covStim));
        covBase = (1-lambda) * covBase + lambda * meanDiagBase * eye(size(covBase));
        covStim = (1-lambda) * covStim + lambda * meanDiagStim * eye(size(covStim));
        
        % Perform GED: eig(covStim, covBase)
        % This finds components that maximize signal (stimulus) relative to noise (baseline)
        [eigenvecs, eigenvals] = eig(covStim, covBase);
        eigenvals = diag(eigenvals);
        
        % Sort by eigenvalue (descending)
        [eigenvals_sorted, sortIdx] = sort(eigenvals, 'descend');
        eigenvecs_sorted = eigenvecs(:, sortIdx);
        
        % Store eigenvalues
        all_eigenvals{cond, subj} = eigenvals_sorted(1:nComponents);
        
        % Extract top component
        topComp = eigenvecs_sorted(:, 1);
        
        % Handle sign uncertainty (Cohen 2022, Section 3.5)
        % Flip eigenvector so that the largest absolute value in the component map is positive
        topo_temp = covStim * topComp;
        [~, maxIdx] = max(abs(topo_temp));
        if topo_temp(maxIdx) < 0
            topComp = -topComp;
        end
        
        % Compute topography (forward model)
        % Component map: covStim * eigenvector (Cohen 2022, Section 3.6)
        topo = covStim * topComp;
        % Ensure topo is a column vector matching number of channels
        if size(topo, 1) ~= length(dat.label)
            topo = topo(:);  % Force column vector
            if length(topo) ~= length(dat.label)
                warning('Topo vector length (%d) does not match number of channels (%d)', ...
                    length(topo), length(dat.label));
            end
        end
        all_topos{cond, subj} = topo;
        
        % Store channel labels for topoplot (once per subject)
        if cond == 1
            chanlocs_all = dat.label;  % Store channel labels
        end
        
        % Project data onto top component to get component time series
        % Component time series: w' * data (Cohen 2022, Section 3.6)
        % Note: The spatial filter can be applied to any data, not just the data
        % used to construct the covariance matrices
        % Baseline period
        base_comp_ts = [];
        for trl = 1:length(dat_base.trial)
            base_comp_ts = [base_comp_ts, topComp' * double(dat_base.trial{trl})];
        end
        
        % Stimulus period
        stim_comp_ts = [];
        for trl = 1:length(dat_stim.trial)
            stim_comp_ts = [stim_comp_ts, topComp' * double(dat_stim.trial{trl})];
        end
        
        % Get time vectors
        base_times = [];
        for trl = 1:length(dat_base.time)
            base_times = [base_times, dat_base.time{trl}];
        end
        
        stim_times = [];
        for trl = 1:length(dat_stim.time)
            stim_times = [stim_times, dat_stim.time{trl}];
        end
        
        all_component_ts_base{cond, subj} = base_comp_ts;
        all_component_ts_stim{cond, subj} = stim_comp_ts;
        all_times_base{cond, subj} = base_times;
        all_times_stim{cond, subj} = stim_times;
        
        % Compute power spectra using FFT
        % Baseline power spectrum
        [pow_base, freq_base] = compute_power_spectrum(base_comp_ts, fsample);
        
        % Stimulus power spectrum
        [pow_stim, freq_stim] = compute_power_spectrum(stim_comp_ts, fsample);
        
        % Find gamma band indices
        gamma_idx = freq_stim >= gamma_freq(1) & freq_stim <= gamma_freq(2);
        freq_gamma = freq_stim(gamma_idx);
        pow_base_gamma = pow_base(gamma_idx);
        pow_stim_gamma = pow_stim(gamma_idx);
        
        % Baseline-corrected power (difference)
        pow_diff_gamma = pow_stim_gamma - pow_base_gamma;
        
        % Relative change: (stimulus - baseline) / baseline
        % Avoid division by zero
        pow_base_gamma_safe = pow_base_gamma;
        pow_base_gamma_safe(pow_base_gamma_safe == 0) = eps;  % Replace zeros with small value
        pow_rel_gamma = pow_diff_gamma ./ pow_base_gamma_safe;
        
        % Store power spectra
        all_powerspec_stim{cond, subj} = pow_stim_gamma;
        all_powerspec_base{cond, subj} = pow_base_gamma;
        all_powerspec_diff{cond, subj} = pow_diff_gamma;
        all_powerspec_rel{cond, subj} = pow_rel_gamma;
        all_freqs{cond, subj} = freq_gamma;
        
        % Find peak gamma frequency and power
        % Use baseline-corrected power for peak detection
        [peaks, locs] = findpeaks(pow_diff_gamma, freq_gamma, ...
            'MinPeakHeight', max(pow_diff_gamma) * 0.3, ...  % At least 30% of max
            'MinPeakDistance', 5);  % At least 5 Hz apart
        
        if ~isempty(peaks)
            [peak_power, peak_idx] = max(peaks);
            peak_freq = locs(peak_idx);
            
            % Compute mean power in ±5 Hz range around peak
            freq_range = peak_freq + [-5, 5];
            range_idx = freq_gamma >= freq_range(1) & freq_gamma <= freq_range(2);
            peak_power_mean = mean(pow_diff_gamma(range_idx));
            
            all_peak_freqs(cond, subj) = peak_freq;
            all_peak_powers(cond, subj) = peak_power_mean;
        else
            % If no peak found, use maximum
            [peak_power_mean, max_idx] = max(pow_diff_gamma);
            peak_freq = freq_gamma(max_idx);
            all_peak_freqs(cond, subj) = peak_freq;
            all_peak_powers(cond, subj) = peak_power_mean;
        end
        
    end  % End condition loop
    
    %% Create individual subject figure
    create_subject_figure(subj, subjects_pilot{subj}, condLabels, ...
        all_topos(:, subj), all_powerspec_stim(:, subj), ...
        all_powerspec_base(:, subj), all_powerspec_diff(:, subj), ...
        all_powerspec_rel(:, subj), all_freqs(:, subj), ...
        all_peak_freqs(:, subj), all_peak_powers(:, subj), ...
        all_eigenvals(:, subj), all_component_ts_stim(:, subj), ...
        all_component_ts_base(:, subj), all_times_stim(:, subj), ...
        all_times_base(:, subj), chanlocs_all, headmodel.layANThead, colors, ...
        fig_save_dir);
    
end  % End subject loop

%% Create grand average figure
fprintf('\n========================================\n');
fprintf('Creating Grand Average Figure\n');
fprintf('========================================\n');

% Get channel labels from first subject
chanlocs_ga = [];
for subj = 1:nPilot
    if ~isempty(all_topos{1, subj})
        % Load data to get channel labels
        datapath_temp = strcat(path, subjects_pilot{subj}, filesep, 'eeg');
        cd(datapath_temp);
        load dataEEG;
        chanlocs_ga = dataEEG_c25.label;
        break;
    end
end

create_grand_average_figure(condLabels, all_topos, all_powerspec_stim, ...
    all_powerspec_base, all_powerspec_diff, all_powerspec_rel, all_freqs, ...
    all_peak_freqs, all_peak_powers, all_eigenvals, all_component_ts_stim, ...
    all_component_ts_base, all_times_stim, all_times_base, layANThead, colors, ...
    fig_save_dir, nPilot, chanlocs_ga);

fprintf('\n========================================\n');
fprintf('GED Analysis Complete!\n');
fprintf('========================================\n');

%% Helper function: Compute power spectrum
function [power, freq] = compute_power_spectrum(signal, fsample)
    % Compute power spectrum using FFT
    N = length(signal);
    Y = fft(signal);
    P = abs(Y/N).^2;
    P = P(1:N/2+1);
    P(2:end-1) = 2*P(2:end-1);  % Single-sided spectrum
    
    freq = fsample * (0:(N/2)) / N;
    power = P;
end

%% Helper function: Create individual subject figure
function create_subject_figure(subj, subjID, condLabels, topos, powerspec_stim, ...
    powerspec_base, powerspec_diff, powerspec_rel, freqs, peak_freqs, peak_powers, ...
    eigenvals, comp_ts_stim, comp_ts_base, times_stim, times_base, chanlocs, ...
    layANThead, colors, fig_save_dir)
    
    figure('Position', [0, 0, 2000, 1400], 'Color', 'w');
    sgtitle(sprintf('GED Gamma Analysis: Subject %s', subjID), ...
        'FontSize', 24, 'FontWeight', 'bold');
    
    nCond = 4;
    
    for cond = 1:nCond
        
        % Row 1: Topography
        subplot(4, nCond, cond);
        if ~isempty(topos{cond})
            % Create FieldTrip structure for topoplot
            topo_data = [];
            topo_data.label = chanlocs;
            topo_data.avg = topos{cond};
            topo_data.dimord = 'chan';
            
            cfg = [];
            cfg.layout = layANThead;
            cfg.comment = 'no';
            cfg.marker = 'off';
            cfg.style = 'straight';
            cfg.gridscale = 300;
            cfg.zlim = 'maxabs';
            cfg.colormap = 'parula';
            
            try
                ft_topoplotER(cfg, topo_data);
            catch
                % Fallback: use imagesc if topoplot fails
                imagesc(topo_data.avg);
                colorbar;
                title('Topography (Top Component)', 'FontSize', 12);
            end
            title(sprintf('%s\nTopography (Top Component)', condLabels{cond}), ...
                'FontSize', 14, 'FontWeight', 'bold');
        end
        
        % Row 2: Power spectra (Raw stimulus, Difference, Relative change)
        subplot(4, nCond, cond + nCond);
        if ~isempty(powerspec_stim{cond}) && ~isempty(freqs{cond})
            % Left y-axis: Raw power values
            yyaxis left;
            hold on;
            % Raw stimulus power spectrum
            plot(freqs{cond}, powerspec_stim{cond}, 'Color', colors(cond, :), ...
                'LineWidth', 2, 'DisplayName', 'Stimulus (Raw)');
            % Baseline power
            plot(freqs{cond}, powerspec_base{cond}, 'Color', colors(cond, :), ...
                'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Baseline');
            % Difference (stimulus - baseline)
            plot(freqs{cond}, powerspec_diff{cond}, 'Color', colors(cond, :), ...
                'LineStyle', ':', 'LineWidth', 2.5, 'DisplayName', 'Difference');
            
            ylabel('Power', 'FontSize', 12);
            
            % Right y-axis: Relative change
            if ~isempty(powerspec_rel{cond})
                yyaxis right;
                plot(freqs{cond}, powerspec_rel{cond} * 100, 'Color', [0.5 0.5 0.5], ...
                    'LineStyle', '-.', 'LineWidth', 2, 'DisplayName', 'Relative Change (%)');
                ylabel('Relative Change (%)', 'FontSize', 12);
                yyaxis left;  % Switch back to left for peak marker
            end
            
            if ~isnan(peak_freqs(cond))
                xline(peak_freqs(cond), '--', 'LineWidth', 2, 'Color', colors(cond, :));
                text(peak_freqs(cond), max(powerspec_diff{cond}) * 0.9, ...
                    sprintf('%.1f Hz', peak_freqs(cond)), ...
                    'FontSize', 10, 'Color', colors(cond, :), ...
                    'HorizontalAlignment', 'center');
            end
            
            xlabel('Frequency (Hz)', 'FontSize', 12);
            title(sprintf('Power Spectrum\nPeak: %.1f Hz, Power: %.4f', ...
                peak_freqs(cond), peak_powers(cond)), 'FontSize', 12);
            legend('Location', 'best', 'FontSize', 9);
            grid on;
            xlim([30, 90]);
        end
        
        % Row 3: Component time series (stimulus period)
        subplot(4, nCond, cond + 2*nCond);
        if ~isempty(comp_ts_stim{cond}) && ~isempty(times_stim{cond})
            plot(times_stim{cond}, comp_ts_stim{cond}, 'Color', colors(cond, :), ...
                'LineWidth', 1);
            xlabel('Time (s)', 'FontSize', 12);
            ylabel('Amplitude', 'FontSize', 12);
            title('Component Time Series (Stimulus)', 'FontSize', 12);
            xlim([0, 2]);
            grid on;
        end
        
        % Row 4: Eigenvalues
        subplot(4, nCond, cond + 3*nCond);
        if ~isempty(eigenvals{cond})
            bar(1:length(eigenvals{cond}), eigenvals{cond}, ...
                'FaceColor', colors(cond, :), 'EdgeColor', 'k');
            xlabel('Component #', 'FontSize', 12);
            ylabel('Eigenvalue', 'FontSize', 12);
            title(sprintf('Top %d Components', length(eigenvals{cond})), ...
                'FontSize', 12);
            grid on;
        end
        
    end
    
    % Save figure
    saveas(gcf, fullfile(fig_save_dir, sprintf('GED_subj%s.png', subjID)), 'png');
    fprintf('  Saved figure for subject %s\n', subjID);
end

%% Helper function: Create grand average figure
function create_grand_average_figure(condLabels, all_topos, all_powerspec_stim, ...
    all_powerspec_base, all_powerspec_diff, all_powerspec_rel, all_freqs, ...
    all_peak_freqs, all_peak_powers, all_eigenvals, all_component_ts_stim, ...
    all_component_ts_base, all_times_stim, all_times_base, layANThead, colors, ...
    fig_save_dir, nSubj, chanlocs_ga)
    
    nCond = 4;
    
    % Compute grand averages
    % Topographies
    ga_topos = cell(nCond, 1);
    for cond = 1:nCond
        valid_topos = [];
        for subj = 1:nSubj
            if ~isempty(all_topos{cond, subj})
                if isempty(valid_topos)
                    valid_topos = all_topos{cond, subj};
                else
                    valid_topos = valid_topos + all_topos{cond, subj};
                end
            end
        end
        if ~isempty(valid_topos)
            ga_topos{cond} = valid_topos / nSubj;
        end
    end
    
    % Power spectra (interpolate to common frequency grid)
    freq_common = 30:0.5:90;  % Common frequency grid
    ga_powerspec_stim = cell(nCond, 1);
    ga_powerspec_base = cell(nCond, 1);
    ga_powerspec_diff = cell(nCond, 1);
    ga_powerspec_rel = cell(nCond, 1);
    
    for cond = 1:nCond
        all_pow_stim = [];
        all_pow_base = [];
        all_pow_diff = [];
        all_pow_rel = [];
        
        for subj = 1:nSubj
            if ~isempty(all_powerspec_stim{cond, subj}) && ~isempty(all_freqs{cond, subj})
                % Interpolate to common grid
                pow_stim_interp = interp1(all_freqs{cond, subj}, ...
                    all_powerspec_stim{cond, subj}, freq_common, 'linear', 'extrap');
                pow_base_interp = interp1(all_freqs{cond, subj}, ...
                    all_powerspec_base{cond, subj}, freq_common, 'linear', 'extrap');
                pow_diff_interp = interp1(all_freqs{cond, subj}, ...
                    all_powerspec_diff{cond, subj}, freq_common, 'linear', 'extrap');
                
                % Relative change
                if ~isempty(all_powerspec_rel{cond, subj})
                    pow_rel_interp = interp1(all_freqs{cond, subj}, ...
                        all_powerspec_rel{cond, subj}, freq_common, 'linear', 'extrap');
                    all_pow_rel = [all_pow_rel; pow_rel_interp];
                end
                
                all_pow_stim = [all_pow_stim; pow_stim_interp];
                all_pow_base = [all_pow_base; pow_base_interp];
                all_pow_diff = [all_pow_diff; pow_diff_interp];
            end
        end
        
        if ~isempty(all_pow_stim)
            ga_powerspec_stim{cond} = mean(all_pow_stim, 1);
            ga_powerspec_base{cond} = mean(all_pow_base, 1);
            ga_powerspec_diff{cond} = mean(all_pow_diff, 1);
            if ~isempty(all_pow_rel)
                ga_powerspec_rel{cond} = mean(all_pow_rel, 1);
            end
        end
    end
    
    % Mean peak frequencies and powers
    ga_peak_freqs = nanmean(all_peak_freqs, 2);
    ga_peak_powers = nanmean(all_peak_powers, 2);
    
    % Get channel locations from first available subject
    chanlocs = [];
    for subj = 1:nSubj
        if ~isempty(all_topos{1, subj})
            % Try to get from a data structure - we'll use a workaround
            break;
        end
    end
    
    % Create figure
    figure('Position', [0, 0, 2000, 1400], 'Color', 'w');
    sgtitle(sprintf('GED Gamma Analysis: Grand Average (N=%d)', nSubj), ...
        'FontSize', 24, 'FontWeight', 'bold');
    
    for cond = 1:nCond
        
        % Row 1: Topography
        subplot(4, nCond, cond);
        if ~isempty(ga_topos{cond}) && ~isempty(chanlocs_ga)
            % Create FieldTrip structure for topoplot
            topo_data = [];
            topo_data.label = chanlocs_ga;
            topo_data.avg = ga_topos{cond};
            topo_data.dimord = 'chan';
            
            cfg = [];
            cfg.layout = layANThead;
            cfg.comment = 'no';
            cfg.marker = 'off';
            cfg.style = 'straight';
            cfg.gridscale = 300;
            cfg.zlim = 'maxabs';
            cfg.colormap = 'parula';
            
            try
                ft_topoplotER(cfg, topo_data);
            catch ME
                % Fallback: use imagesc if topoplot fails
                warning('Topoplot failed: %s. Using imagesc instead.', ME.message);
                imagesc(ga_topos{cond});
                colorbar;
            end
            title(sprintf('%s\nTopography (Top Component)', condLabels{cond}), ...
                'FontSize', 14, 'FontWeight', 'bold');
        end
        
        % Row 2: Power spectra (Raw stimulus, Difference, Relative change)
        subplot(4, nCond, cond + nCond);
        if ~isempty(ga_powerspec_stim{cond})
            % Left y-axis: Raw power values
            yyaxis left;
            hold on;
            % Raw stimulus power spectrum
            plot(freq_common, ga_powerspec_stim{cond}, 'Color', colors(cond, :), ...
                'LineWidth', 2, 'DisplayName', 'Stimulus (Raw)');
            % Baseline power
            plot(freq_common, ga_powerspec_base{cond}, 'Color', colors(cond, :), ...
                'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Baseline');
            % Difference (stimulus - baseline)
            plot(freq_common, ga_powerspec_diff{cond}, 'Color', colors(cond, :), ...
                'LineStyle', ':', 'LineWidth', 2.5, 'DisplayName', 'Difference');
            
            ylabel('Power', 'FontSize', 12);
            
            % Right y-axis: Relative change
            if ~isempty(ga_powerspec_rel{cond})
                yyaxis right;
                plot(freq_common, ga_powerspec_rel{cond} * 100, 'Color', [0.5 0.5 0.5], ...
                    'LineStyle', '-.', 'LineWidth', 2, 'DisplayName', 'Relative Change (%)');
                ylabel('Relative Change (%)', 'FontSize', 12);
                yyaxis left;  % Switch back to left for peak marker
            end
            
            if ~isnan(ga_peak_freqs(cond))
                xline(ga_peak_freqs(cond), '--', 'LineWidth', 2, 'Color', colors(cond, :));
                text(ga_peak_freqs(cond), max(ga_powerspec_diff{cond}) * 0.9, ...
                    sprintf('%.1f Hz', ga_peak_freqs(cond)), ...
                    'FontSize', 10, 'Color', colors(cond, :), ...
                    'HorizontalAlignment', 'center');
            end
            
            xlabel('Frequency (Hz)', 'FontSize', 12);
            title(sprintf('Power Spectrum\nPeak: %.1f Hz, Power: %.4f', ...
                ga_peak_freqs(cond), ga_peak_powers(cond)), 'FontSize', 12);
            legend('Location', 'best', 'FontSize', 9);
            grid on;
            xlim([30, 90]);
        end
        
        % Row 3: Mean peak frequencies across subjects
        subplot(4, nCond, cond + 2*nCond);
        peak_freqs_cond = all_peak_freqs(cond, :);
        peak_freqs_cond = peak_freqs_cond(~isnan(peak_freqs_cond));
        if ~isempty(peak_freqs_cond)
            histogram(peak_freqs_cond, 10, 'FaceColor', colors(cond, :), ...
                'EdgeColor', 'k');
            xlabel('Peak Frequency (Hz)', 'FontSize', 12);
            ylabel('Number of Subjects', 'FontSize', 12);
            title(sprintf('Peak Frequency Distribution\nMean: %.1f Hz', ...
                mean(peak_freqs_cond)), 'FontSize', 12);
            grid on;
            xlim([30, 90]);
        end
        
        % Row 4: Mean peak powers across subjects
        subplot(4, nCond, cond + 3*nCond);
        peak_powers_cond = all_peak_powers(cond, :);
        peak_powers_cond = peak_powers_cond(~isnan(peak_powers_cond));
        if ~isempty(peak_powers_cond)
            histogram(peak_powers_cond, 10, 'FaceColor', colors(cond, :), ...
                'EdgeColor', 'k');
            xlabel('Peak Power', 'FontSize', 12);
            ylabel('Number of Subjects', 'FontSize', 12);
            title(sprintf('Peak Power Distribution\nMean: %.4f', ...
                mean(peak_powers_cond)), 'FontSize', 12);
            grid on;
        end
        
    end
    
    % Save figure
    saveas(gcf, fullfile(fig_save_dir, 'GED_grand_average.png'), 'png');
    fprintf('Saved grand average figure\n');
end
