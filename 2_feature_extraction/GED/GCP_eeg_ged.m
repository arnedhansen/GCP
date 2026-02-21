%% GCP Narrowband Scanning GED (Common Spatial Filter)
%
% Uses a COMMON spatial filter across all contrast conditions to guarantee
% the same gamma source is compared. The filter is derived from pooled
% covariance across all 4 conditions.
%
% Narrowband Power Scanning (adapted from Cohen, 2021):
%   - Pool all conditions -> broadband GED -> common spatial filter
%   - For each condition, at each candidate frequency (30-90 Hz):
%     bandpass filter, project through common filter, compute stimulus vs.
%     baseline power ratio
%   - Detrend the power-ratio spectrum (2nd-order polynomial)
%
% Two peak detection approaches:
%   1. Single peak — tallest peak in the full detrended spectrum
%   2. Dual peak   — assuming low + high gamma sub-bands:
%                     low gamma peak  (30-55 Hz)
%                     high gamma peak (50-85 Hz)
%                     overlap at 50-55 Hz so mid-range peaks count for both

clear; close all; clc

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');

nSubj = length(subjects);

% Time windows
baseline_window = [-1.5, -0.25];
stimulus_window = [0.3, 2.0];

% Gamma frequency range
gamma_range = [30, 90];

% Narrowband scanning parameters
scan_freqs = 30:1:90;
scan_width = 4; % +-2 Hz bandwidth at each step

% GED parameters
lambda = 0.01; % shrinkage regularization

% Detrending parameters for power-ratio spectrum
poly_order = 2;

% Condition info
condNames  = {'c25', 'c50', 'c75', 'c100'};
trialCodes = [61, 62, 63, 64];
condLabels = {'25%', '50%', '75%', '100%'};

% Figure save directory
if ispc
    fig_save_dir = 'W:\Students\Arne\GCP\figures\eeg\ged';
else
    fig_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged';
end
if ~exist(fig_save_dir, 'dir'), mkdir(fig_save_dir); end

%% Preallocate storage
% Method 1: Narrowband power scanning (through common filter)
all_powratio      = cell(4, nSubj);  % raw power-ratio spectrum
all_powratio_dt   = cell(4, nSubj);  % detrended power-ratio spectrum
all_scan_peakfreq = nan(4, nSubj);   % peak freq from detrended spectrum
all_scan_peakpow  = nan(4, nSubj);   % peak value (detrended)

% Dual-peak model (low + high gamma)
all_peak_low     = nan(4, nSubj);  % low gamma peak  (30-55 Hz)
all_peak_high    = nan(4, nSubj);  % high gamma peak (50-85 Hz)

% Common filter info (one per subject)
all_topos          = cell(1, nSubj);
all_eigenvalues    = nan(1, nSubj);

% GED component power spectra
all_powspec_stim  = cell(4, nSubj);
all_powspec_base  = cell(4, nSubj);
all_powspec_diff  = cell(4, nSubj);
all_powspec_freqs = cell(4, nSubj);

% Channel labels
chanlocs_all = {};

%% Process each subject
for subj = 1:nSubj

    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load dataEEG

    fsample = dataEEG_c25.fsample;

    trialIndices = { ...
        find(dataEEG_c25.trialinfo  == 61), ...
        find(dataEEG_c50.trialinfo  == 62), ...
        find(dataEEG_c75.trialinfo  == 63), ...
        find(dataEEG_c100.trialinfo == 64)};
    dataStructs = {dataEEG_c25, dataEEG_c50, dataEEG_c75, dataEEG_c100};

    if subj == 1
        chanlocs_all = dataEEG_c25.label;
    end

    nChans = length(dataEEG_c25.label);

    %% ================================================================
    %  PHASE 1: Build POOLED covariance across all conditions -> one GED
    %  ================================================================
    clc
    fprintf('Subject %s (%d/%d) — Phase 1: Pooled GED\n', ...
        subjects{subj}, subj, nSubj);

    covStim_pooled = zeros(nChans);
    covBase_pooled = zeros(nChans);
    nTrials_total  = 0;

    dat_per_cond = cell(1, 4);

    for cond = 1:4
        dat    = dataStructs{cond};
        trlIdx = trialIndices{cond};
        if isempty(trlIdx), continue; end

        cfg = [];
        cfg.trials = trlIdx;
        dat = ft_selectdata(cfg, dat);
        dat_per_cond{cond} = dat;

        % Broadband gamma filter for covariance
        cfg_filt = [];
        cfg_filt.bpfilter   = 'yes';
        cfg_filt.bpfreq     = gamma_range;
        cfg_filt.bpfilttype = 'fir';
        cfg_filt.bpfiltord  = round(3 * fsample / gamma_range(1));
        dat_gamma = ft_preprocessing(cfg_filt, dat);

        cfg_t = [];
        cfg_t.latency = baseline_window;
        dat_base = ft_selectdata(cfg_t, dat_gamma);

        cfg_t.latency = stimulus_window;
        dat_stim = ft_selectdata(cfg_t, dat_gamma);

        nTrl = length(dat_stim.trial);
        for trl = 1:nTrl
            d = double(dat_stim.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            covStim_pooled = covStim_pooled + (d * d') / size(d, 2);

            d = double(dat_base.trial{trl});
            d = bsxfun(@minus, d, mean(d, 2));
            covBase_pooled = covBase_pooled + (d * d') / size(d, 2);
        end
        nTrials_total = nTrials_total + nTrl;
    end

    covStim_pooled = covStim_pooled / nTrials_total;
    covBase_pooled = covBase_pooled / nTrials_total;

    % Shrinkage regularization
    covStim_pooled = (1-lambda)*covStim_pooled + lambda*mean(diag(covStim_pooled))*eye(nChans);
    covBase_pooled = (1-lambda)*covBase_pooled + lambda*mean(diag(covBase_pooled))*eye(nChans);

    % GED on pooled covariance
    [W, D] = eig(covStim_pooled, covBase_pooled);
    [evals_sorted, sortIdx] = sort(real(diag(D)), 'descend');
    W = W(:, sortIdx);

    topComp = W(:, 1);
    all_eigenvalues(subj) = evals_sorted(1);

    % Sign correction
    topo_temp = covStim_pooled * topComp;
    [~, mxI] = max(abs(topo_temp));
    if topo_temp(mxI) < 0, topComp = -topComp; end

    all_topos{subj} = covStim_pooled * topComp;

    %% ================================================================
    %  PHASE 2: Per condition — narrowband scanning + power spectrum
    %  ================================================================
    for cond = 1:4

        dat = dat_per_cond{cond};
        if isempty(dat), continue; end

        %% ============================================================
        %  METHOD 1: Narrowband Power Scanning through common filter
        %  ============================================================
        powratio = nan(1, length(scan_freqs));

        for fi = 1:length(scan_freqs)
            clc
            fprintf('Subject    %s (%d/%d)\nCondition  %d/4\nFrequency  %d/%d\n', ...
                subjects{subj}, subj, nSubj, cond, fi, length(scan_freqs));
            cf = scan_freqs(fi);
            bpfreq = [max(cf - scan_width/2, 1), cf + scan_width/2];

            cfg_filt = [];
            cfg_filt.bpfilter   = 'yes';
            cfg_filt.bpfreq     = bpfreq;
            cfg_filt.bpfilttype = 'fir';
            cfg_filt.bpfiltord  = round(3 * fsample / bpfreq(1));
            dat_nb = ft_preprocessing(cfg_filt, dat);

            cfg_t = [];
            cfg_t.latency = baseline_window;
            dat_base_nb = ft_selectdata(cfg_t, dat_nb);

            cfg_t.latency = stimulus_window;
            dat_stim_nb = ft_selectdata(cfg_t, dat_nb);

            % Project through common filter and compute power
            pow_stim = 0;
            for trl = 1:length(dat_stim_nb.trial)
                comp_ts = topComp' * double(dat_stim_nb.trial{trl});
                pow_stim = pow_stim + mean(comp_ts.^2);
            end
            pow_stim = pow_stim / length(dat_stim_nb.trial);

            pow_base = 0;
            for trl = 1:length(dat_base_nb.trial)
                comp_ts = topComp' * double(dat_base_nb.trial{trl});
                pow_base = pow_base + mean(comp_ts.^2);
            end
            pow_base = pow_base / length(dat_base_nb.trial);

            powratio(fi) = pow_stim / pow_base;
        end

        all_powratio{cond, subj} = powratio;

        % Detrend power-ratio spectrum
        p = polyfit(scan_freqs, powratio, poly_order);
        powratio_trend = polyval(p, scan_freqs);
        powratio_dt = powratio - powratio_trend;

        all_powratio_dt{cond, subj} = powratio_dt;

        % Find peak in detrended spectrum
        powratio_dt_smooth = movmean(powratio_dt, 5);
        [pks, locs] = findpeaks(powratio_dt_smooth, scan_freqs, ...
            'MinPeakProminence', max(powratio_dt_smooth) * 0.15, ...
            'MinPeakDistance', 5);

        if ~isempty(pks)
            [~, best_pk] = max(pks);
            scan_peak_freq = locs(best_pk);
        else
            [~, mi] = max(powratio_dt_smooth);
            scan_peak_freq = scan_freqs(mi);
        end
        scan_peak_pow = powratio_dt(scan_freqs == scan_peak_freq);

        all_scan_peakfreq(cond, subj) = scan_peak_freq;
        all_scan_peakpow(cond, subj)  = scan_peak_pow;

        %% --- Dual-peak model: low gamma + high gamma ---
        % Find all positive peaks in the full detrended spectrum
        [pks_all, locs_all] = findpeaks(powratio_dt_smooth, scan_freqs, ...
            'MinPeakDistance', 5);
        pos_mask = pks_all > 0;
        pks_pos  = pks_all(pos_mask);
        locs_pos = locs_all(pos_mask);
        nPosPeaks = length(pks_pos);

        if nPosPeaks >= 2
            % 2+ peaks: find tallest in each sub-range (overlapping middle)
            in_lo = locs_pos >= 30 & locs_pos <= 55;
            in_hi = locs_pos >= 50 & locs_pos <= 85;
            if any(in_lo)
                [~, bi] = max(pks_pos(in_lo));
                tmp = locs_pos(in_lo);
                all_peak_low(cond, subj) = tmp(bi);
            end
            if any(in_hi)
                [~, bi] = max(pks_pos(in_hi));
                tmp = locs_pos(in_hi);
                all_peak_high(cond, subj) = tmp(bi);
            end
        elseif nPosPeaks == 1
            % Single peak: assign based on middle margin (45-65 Hz)
            the_peak = locs_pos(1);
            if the_peak >= 45 && the_peak <= 65
                all_peak_low(cond, subj)  = the_peak;
                all_peak_high(cond, subj) = the_peak;
            elseif the_peak < 45
                all_peak_low(cond, subj) = the_peak;
            else
                all_peak_high(cond, subj) = the_peak;
            end
        end

        %% ============================================================
        %  GED Component Power Spectrum (through common filter)
        %  ============================================================
        cfg_t_pow = [];
        cfg_t_pow.latency = stimulus_window;
        dat_stim_raw = ft_selectdata(cfg_t_pow, dat);

        cfg_t_pow.latency = baseline_window;
        dat_base_raw = ft_selectdata(cfg_t_pow, dat);

        stim_comp = [];
        for trl = 1:length(dat_stim_raw.trial)
            stim_comp = [stim_comp, topComp' * double(dat_stim_raw.trial{trl})];
        end

        base_comp = [];
        for trl = 1:length(dat_base_raw.trial)
            base_comp = [base_comp, topComp' * double(dat_base_raw.trial{trl})];
        end

        N_stim = length(stim_comp);
        Y_stim = fft(stim_comp);
        P_stim = abs(Y_stim / N_stim).^2;
        P_stim = P_stim(1:floor(N_stim/2)+1);
        P_stim(2:end-1) = 2 * P_stim(2:end-1);
        freq_stim = fsample * (0:floor(N_stim/2)) / N_stim;

        N_base = length(base_comp);
        Y_base = fft(base_comp);
        P_base = abs(Y_base / N_base).^2;
        P_base = P_base(1:floor(N_base/2)+1);
        P_base(2:end-1) = 2 * P_base(2:end-1);
        freq_base = fsample * (0:floor(N_base/2)) / N_base;

        freq_common_pow = 30:0.5:90;
        pow_stim_interp = interp1(freq_stim, P_stim, freq_common_pow, 'linear', 0);
        pow_base_interp = interp1(freq_base, P_base, freq_common_pow, 'linear', 0);
        pow_diff_interp = pow_stim_interp - pow_base_interp;

        all_powspec_stim{cond, subj}  = pow_stim_interp;
        all_powspec_base{cond, subj}  = pow_base_interp;
        all_powspec_diff{cond, subj}  = pow_diff_interp;
        all_powspec_freqs{cond, subj} = freq_common_pow;

    end % condition loop

    %% ================================================================
    %  PER-SUBJECT FIGURE (3 rows × 4 columns)
    %  ================================================================
    close all
    fig = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('GED Gamma (Common Filter): Subject %s', subjects{subj}), ...
        'FontSize', 20, 'FontWeight', 'bold');

    % --- Row 1: Raw power-ratio spectra ---
    for cond = 1:4
        subplot(4, 4, cond); hold on;
        if ~isempty(all_powratio{cond, subj})
            plot(scan_freqs, all_powratio{cond, subj}, '-', ...
                'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            plot(scan_freqs, movmean(all_powratio{cond, subj}, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2);
            p = polyfit(scan_freqs, all_powratio{cond, subj}, poly_order);
            plot(scan_freqs, polyval(p, scan_freqs), 'k--', 'LineWidth', 1.5);
        end
        xlabel('Freq [Hz]'); ylabel('Power Ratio');
        title(sprintf('%s Raw', condLabels{cond}), 'FontSize', 12);
        set(gca, 'FontSize', 11); xlim([30 90]); grid on; box on;
    end

    % --- Row 2: Detrended with single peak ---
    for cond = 1:4
        subplot(4, 4, 4 + cond); hold on;
        if ~isempty(all_powratio_dt{cond, subj})
            plot(scan_freqs, all_powratio_dt{cond, subj}, '-', ...
                'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            plot(scan_freqs, movmean(all_powratio_dt{cond, subj}, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2.5);
            yline(0, 'k-', 'LineWidth', 0.5);
            xline(all_scan_peakfreq(cond, subj), '--', 'LineWidth', 2, ...
                'Color', colors(cond,:));
            text(all_scan_peakfreq(cond, subj) + 1, ...
                max(all_powratio_dt{cond, subj}) * 0.9, ...
                sprintf('%.0f Hz', all_scan_peakfreq(cond, subj)), ...
                'FontSize', 12, 'Color', colors(cond,:), 'FontWeight', 'bold');
        end
        xlabel('Freq [Hz]'); ylabel('\Delta PR');
        title(sprintf('%s Single Peak', condLabels{cond}), 'FontSize', 12);
        set(gca, 'FontSize', 11); xlim([30 90]); grid on; box on;
    end

    % --- Row 3: Detrended with low + high gamma peaks ---
    for cond = 1:4
        subplot(4, 4, 8 + cond); hold on;
        if ~isempty(all_powratio_dt{cond, subj})
            plot(scan_freqs, all_powratio_dt{cond, subj}, '-', ...
                'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            plot(scan_freqs, movmean(all_powratio_dt{cond, subj}, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2.5);
            yline(0, 'k-', 'LineWidth', 0.5);
            xline(50, 'k:', 'LineWidth', 0.8, 'Alpha', 0.4);
            xline(55, 'k:', 'LineWidth', 0.8, 'Alpha', 0.4);
            pf_lo = all_peak_low(cond, subj);
            pf_hi = all_peak_high(cond, subj);
            if ~isnan(pf_lo)
                xline(pf_lo, '--', 'LineWidth', 2, 'Color', [0 0 0.7]);
                text(pf_lo + 1, max(all_powratio_dt{cond, subj}) * 0.85, ...
                    sprintf('L:%.0f', pf_lo), 'FontSize', 10, ...
                    'Color', [0 0 0.7], 'FontWeight', 'bold');
            end
            if ~isnan(pf_hi)
                xline(pf_hi, '--', 'LineWidth', 2, 'Color', [0.7 0 0]);
                text(pf_hi + 1, max(all_powratio_dt{cond, subj}) * 0.7, ...
                    sprintf('H:%.0f', pf_hi), 'FontSize', 10, ...
                    'Color', [0.7 0 0], 'FontWeight', 'bold');
            end
        end
        xlabel('Freq [Hz]'); ylabel('\Delta PR');
        title(sprintf('%s Dual Peak', condLabels{cond}), 'FontSize', 12);
        set(gca, 'FontSize', 11); xlim([30 90]); grid on; box on;
    end

    % --- Row 4: Topoplot + all conditions overlay ---
    cfg_topo = [];
    cfg_topo.layout    = headmodel.layANThead;
    cfg_topo.comment   = 'no';
    cfg_topo.marker    = 'off';
    cfg_topo.style     = 'straight';
    cfg_topo.gridscale = 300;
    cfg_topo.zlim      = 'maxabs';
    cfg_topo.colormap  = '*RdBu';
    cfg_topo.figure    = 'gcf';

    subplot(4, 4, 13);
    if ~isempty(all_topos{subj})
        topo_data = [];
        topo_data.label  = chanlocs_all;
        topo_data.avg    = all_topos{subj};
        topo_data.dimord = 'chan';
        try
            ft_topoplotER(cfg_topo, topo_data);
            cb = colorbar; cb.FontSize = 9;
        catch
            imagesc(topo_data.avg); colorbar;
        end
        title(sprintf('Common Filter (\\lambda=%.1f)', all_eigenvalues(subj)), 'FontSize', 12);
    end

    subplot(4, 4, [14 15 16]); hold on;
    for cond = 1:4
        if ~isempty(all_powratio_dt{cond, subj})
            plot(scan_freqs, movmean(all_powratio_dt{cond, subj}, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2.5);
            xline(all_scan_peakfreq(cond, subj), '--', 'Color', colors(cond,:), 'LineWidth', 1.5);
        end
    end
    yline(0, 'k-', 'LineWidth', 0.5);
    xlabel('Freq [Hz]'); ylabel('\Delta PR (detrended)');
    title('All Conditions (Same Spatial Filter)', 'FontSize', 14);
    legend(condLabels, 'FontSize', 11, 'Location', 'best');
    set(gca, 'FontSize', 12); xlim([30 90]); grid on; box on;

    saveas(fig, fullfile(fig_save_dir, sprintf('GED_scan_subj%s.png', subjects{subj})));

end % subject loop

%% ====================================================================
%  GRAND AVERAGE FIGURE 1: Power-ratio spectra + peak scatter
%  ====================================================================
fprintf('\nCreating grand average figures...\n');
close all

fig_ga1 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Narrowband Scanning (Common Filter): Grand Average (N=%d)', nSubj), ...
    'FontSize', 20, 'FontWeight', 'bold');

% --- Top left: Raw power-ratio spectra ---
subplot(2, 2, 1); hold on;
for cond = 1:4
    pr_mat = nan(nSubj, length(scan_freqs));
    for s = 1:nSubj
        if ~isempty(all_powratio{cond, s})
            pr_mat(s,:) = all_powratio{cond, s};
        end
    end
    mu  = nanmean(pr_mat, 1);
    sem = nanstd(pr_mat, [], 1) / sqrt(sum(~isnan(pr_mat(:,1))));
    mu_s  = movmean(mu, 5);
    sem_s = movmean(sem, 5);
    faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
    patch([scan_freqs, fliplr(scan_freqs)], ...
          [mu_s - sem_s, fliplr(mu_s + sem_s)], ...
          colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(scan_freqs, mu_s, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);
end
xlabel('Frequency [Hz]'); ylabel('Power Ratio (stim/base)');
title('Raw Power-Ratio Spectra (mean +/- SEM)', 'FontSize', 14);
legend(condLabels, 'Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 13); xlim([30 90]); grid on; box on;

% --- Top right: Detrended power-ratio spectra ---
subplot(2, 2, 2); hold on;
for cond = 1:4
    prdt_mat = nan(nSubj, length(scan_freqs));
    for s = 1:nSubj
        if ~isempty(all_powratio_dt{cond, s})
            prdt_mat(s,:) = all_powratio_dt{cond, s};
        end
    end
    mu  = nanmean(prdt_mat, 1);
    sem = nanstd(prdt_mat, [], 1) / sqrt(sum(~isnan(prdt_mat(:,1))));
    mu_s  = movmean(mu, 5);
    sem_s = movmean(sem, 5);
    faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
    patch([scan_freqs, fliplr(scan_freqs)], ...
          [mu_s - sem_s, fliplr(mu_s + sem_s)], ...
          colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(scan_freqs, mu_s, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);
end
yline(0, 'k-', 'LineWidth', 0.5);
xlabel('Frequency [Hz]'); ylabel('\Delta PR (detrended)');
title('Detrended Power-Ratio Spectra (mean +/- SEM)', 'FontSize', 14);
legend(condLabels, 'Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 13); xlim([30 90]); grid on; box on;

% --- Bottom: Peak frequency scatter (Original method) ---
subplot(2, 2, [3 4]); hold on;
for cond = 1:4
    pf = all_scan_peakfreq(cond, :);
    pf = pf(~isnan(pf));
    scatter(ones(size(pf)) * cond + 0.15*(rand(size(pf))-0.5), pf, ...
        60, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.6);
    mu_pf  = nanmean(all_scan_peakfreq(cond,:));
    sem_pf = nanstd(all_scan_peakfreq(cond,:)) / sqrt(sum(~isnan(all_scan_peakfreq(cond,:))));
    errorbar(cond, mu_pf, sem_pf, 'k', 'LineWidth', 2, 'CapSize', 10);
    plot(cond, mu_pf, 'kd', 'MarkerSize', 12, ...
        'MarkerFaceColor', colors(cond,:), 'LineWidth', 1.5);
end
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 13);
ylabel('Peak Gamma Frequency [Hz]');
title('Individual Peak Frequencies (Original)', 'FontSize', 14);
ylim([30 90]); grid on; box on;

saveas(fig_ga1, fullfile(fig_save_dir, 'GED_scan_grand_average.png'));

%% ====================================================================
%  GRAND AVERAGE FIGURE 2: All subjects subplot
%  ====================================================================
nRows = ceil(nSubj / 5);
fig_all = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Detrended Power-Ratio Spectra: All Subjects (N=%d)', nSubj), ...
    'FontSize', 18, 'FontWeight', 'bold');

for s = 1:nSubj
    subplot(nRows, 5, s); hold on;
    for cond = 1:4
        if ~isempty(all_powratio_dt{cond, s})
            plot(scan_freqs, movmean(all_powratio_dt{cond, s}, 5), '-', ...
                'Color', colors(cond,:), 'LineWidth', 2);
            xline(all_scan_peakfreq(cond, s), '--', 'Color', colors(cond,:), ...
                'LineWidth', 1.2);
        end
    end
    yline(0, 'k-', 'LineWidth', 0.5);
    xlabel('Freq [Hz]'); ylabel('\Delta PR');
    title(sprintf('Subj %s', subjects{s}), 'FontSize', 11);
    set(gca, 'FontSize', 10); xlim([30 90]); grid on; box on;
    if s == 1
        legend(condLabels, 'FontSize', 8, 'Location', 'best');
    end
end
saveas(fig_all, fullfile(fig_save_dir, 'GED_scan_all_subjects.png'));

%% ====================================================================
%  BOXPLOT FIGURE 1: Single peak
%  ====================================================================
fig_box1 = figure('Position', [0 0 1200 982], 'Color', 'w');
hold on;

for s = 1:nSubj
    pf = all_scan_peakfreq(:, s);
    if sum(~isnan(pf)) >= 2
        plot(1:4, pf, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    end
end

y_box = all_scan_peakfreq(:);
g_box = repelem((1:4)', nSubj, 1);
valid_box = ~isnan(y_box);
boxplot(y_box(valid_box), g_box(valid_box), 'Colors', 'k', 'Symbol', '');

hold on;
for c = 1:4
    pf = all_scan_peakfreq(c, :);
    pf = pf(~isnan(pf));
    xJit = c + (rand(size(pf)) - 0.5) * 0.1;
    scatter(xJit, pf, 250, colors(c,:), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

xlim([0.5 4.5]); ylim([30 90]);
set(gca, 'XTick', 1:4, 'XTickLabel', {'25%', '50%', '75%', '100%'}, ...
    'FontSize', 20, 'Box', 'off');
ylabel('Peak Gamma Frequency [Hz]');
title('Peak Frequency: Single Peak', 'FontSize', 30, 'FontWeight', 'bold');

saveas(fig_box1, fullfile(fig_save_dir, 'GED_boxplot_peakfreq_SinglePeak.png'));

%% ====================================================================
%  BOXPLOT FIGURE 2: Low + High gamma side by side
%  ====================================================================
fig_box2 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Dual-Peak Gamma Frequency', 'FontSize', 30, 'FontWeight', 'bold');

dual_data   = {all_peak_low, all_peak_high};
dual_titles = {'Low Gamma', 'High Gamma'};
dual_ylims  = {[30 65], [40 90]};

for di = 1:2
    subplot(1, 2, di); hold on;
    peak_data = dual_data{di};

    for s = 1:nSubj
        pf = peak_data(:, s);
        if sum(~isnan(pf)) >= 2
            plot(1:4, pf, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
        end
    end

    y_box = peak_data(:);
    g_box = repelem((1:4)', nSubj, 1);
    valid_box = ~isnan(y_box);
    if any(valid_box)
        boxplot(y_box(valid_box), g_box(valid_box), 'Colors', 'k', 'Symbol', '');
    end

    hold on;
    for c = 1:4
        pf = peak_data(c, :);
        pf = pf(~isnan(pf));
        xJit = c + (rand(size(pf)) - 0.5) * 0.1;
        scatter(xJit, pf, 250, colors(c,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end

    xlim([0.5 4.5]); ylim(dual_ylims{di});
    set(gca, 'XTick', 1:4, 'XTickLabel', {'25%', '50%', '75%', '100%'}, ...
        'FontSize', 18, 'Box', 'off');
    ylabel('Peak Gamma Frequency [Hz]');
    title(dual_titles{di}, 'FontSize', 24, 'FontWeight', 'bold');
end

saveas(fig_box2, fullfile(fig_save_dir, 'GED_boxplot_peakfreq_DualGamma.png'));

%% ====================================================================
%  POWER SPECTRUM FIGURE 1: All subjects subplot
%  ====================================================================
freq_common_pow = 30:0.5:90;
nRows_pow = ceil(nSubj / 5);
fig_pow_all = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('GED Component Power Spectra: All Subjects (N=%d)', nSubj), ...
    'FontSize', 18, 'FontWeight', 'bold');

for s = 1:nSubj
    subplot(nRows_pow, 5, s); hold on;
    for cond = 1:4
        if ~isempty(all_powspec_diff{cond, s})
            plot(freq_common_pow, all_powspec_diff{cond, s}, '-', ...
                'Color', colors(cond,:), 'LineWidth', 2);
            xline(all_scan_peakfreq(cond, s), '--', 'Color', colors(cond,:), ...
                'LineWidth', 1.2);
        end
    end
    yline(0, 'k--', 'LineWidth', 0.5);
    xlabel('Freq [Hz]'); ylabel('Power');
    title(sprintf('Subj %s', subjects{s}), 'FontSize', 11);
    set(gca, 'FontSize', 10); xlim([30 90]); grid on; box on;
    if s == 1
        legend(condLabels, 'FontSize', 8, 'Location', 'best');
    end
end
saveas(fig_pow_all, fullfile(fig_save_dir, 'GED_powspec_all_subjects.png'));

%% ====================================================================
%  POWER SPECTRUM FIGURE 2: Grand average
%  ====================================================================
fig_pow_ga = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('GED Component Power Spectrum: Grand Average (N=%d)', nSubj), ...
    'FontSize', 20, 'FontWeight', 'bold');
hold on;

hl_ga = gobjects(1, 4);
for cond = 1:4
    pow_mat = nan(nSubj, length(freq_common_pow));
    for s = 1:nSubj
        if ~isempty(all_powspec_diff{cond, s})
            pow_mat(s,:) = all_powspec_diff{cond, s};
        end
    end
    mu  = nanmean(pow_mat, 1);
    sem = nanstd(pow_mat, [], 1) / sqrt(sum(~isnan(pow_mat(:,1))));

    faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
    patch([freq_common_pow, fliplr(freq_common_pow)], ...
          [mu - sem, fliplr(mu + sem)], ...
          colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hl_ga(cond) = plot(freq_common_pow, mu, '-', 'Color', colors(cond,:), 'LineWidth', 3);

    ga_pf = nanmean(all_scan_peakfreq(cond,:));
    xline(ga_pf, '--', 'Color', colors(cond,:), 'LineWidth', 1.5);
end
yline(0, 'k--', 'LineWidth', 1);
xlabel('Frequency [Hz]'); ylabel('Power (stim - base)');
title('Baseline-Corrected Power Spectrum (mean +/- SEM)', 'FontSize', 18);
legend(hl_ga, condLabels, 'Location', 'northeast', 'FontSize', 15);
set(gca, 'FontSize', 16); xlim([30 90]); grid on; box on;

saveas(fig_pow_ga, fullfile(fig_save_dir, 'GED_powspec_grand_average.png'));

%% Save results
if ispc
    save_path = 'W:\Students\Arne\GCP\data\features\ged_gamma_peaks.mat';
else
    save_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/ged_gamma_peaks.mat';
end
save(save_path, ...
    'all_powratio', 'all_powratio_dt', 'all_topos', 'all_eigenvalues', ...
    'all_scan_peakfreq', 'all_scan_peakpow', ...
    'all_peak_low', 'all_peak_high', ...
    'all_powspec_stim', 'all_powspec_base', 'all_powspec_diff', 'all_powspec_freqs', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames');

clc
fprintf('Done.\n');
