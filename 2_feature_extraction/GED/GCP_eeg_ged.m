%% GCP Narrowband Scanning GED + Frequency Sliding (Common Spatial Filter)
%
% Uses a COMMON spatial filter across all contrast conditions to guarantee
% the same gamma source is compared. The filter is derived from pooled
% covariance across all 4 conditions.
%
% METHOD 1: Narrowband Power Scanning (adapted from Cohen, 2021)
%   - Pool all conditions -> broadband GED -> common spatial filter
%   - For each condition, at each candidate frequency (30-90 Hz):
%     bandpass filter, project through common filter, compute stimulus vs.
%     baseline power ratio
%   - Detrend the power-ratio spectrum (2nd-order polynomial)
%   - Peak of detrended spectrum = individual gamma peak frequency
%
% METHOD 2: Frequency Sliding (Cohen, 2014, J Neurosci)
%   - Apply common spatial filter to broadband gamma data
%   - Hilbert transform -> instantaneous frequency distribution
%   - Mode of distribution = individual gamma peak frequency

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

% Method 2: Frequency sliding (through common filter)
all_fs_hist_freqs  = cell(4, nSubj);
all_fs_hist_counts = cell(4, nSubj);
all_fs_peakfreq    = nan(4, nSubj);   % mode of unweighted IF
all_fs_median      = nan(4, nSubj);
all_fs_wt_mean     = nan(4, nSubj);   % amplitude-weighted mean IF
all_fs_wt_peak     = nan(4, nSubj);   % mode of amplitude-weighted IF histogram

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
    %  PHASE 2: Per condition — scanning + frequency sliding + power spectrum
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

        %% ============================================================
        %  METHOD 2: Frequency Sliding through common filter
        %  ============================================================
        % Narrower bandpass (35-75 Hz) to reduce bias toward band center
        fs_range = [40, 80];
        cfg_broad = [];
        cfg_broad.bpfilter   = 'yes';
        cfg_broad.bpfreq     = fs_range;
        cfg_broad.bpfilttype = 'fir';
        cfg_broad.bpfiltord  = round(3 * fsample / fs_range(1));
        dat_broad = ft_preprocessing(cfg_broad, dat);

        cfg_t = [];
        cfg_t.latency = stimulus_window;
        dat_stim_broad = ft_selectdata(cfg_t, dat_broad);

        all_inst_freq = [];
        all_inst_amp  = [];
        for trl = 1:length(dat_stim_broad.trial)
            comp_ts = topComp' * double(dat_stim_broad.trial{trl});
            analytic = hilbert(comp_ts);
            phase_angles = angle(analytic);
            inst_amp  = abs(analytic);
            dphi = diff(unwrap(phase_angles));
            inst_freq = dphi * fsample / (2 * pi);
            inst_freq = medfilt1(inst_freq, 10);
            inst_amp  = inst_amp(1:end-1);
            valid = inst_freq >= fs_range(1) & inst_freq <= fs_range(2);
            all_inst_freq = [all_inst_freq, inst_freq(valid)];
            all_inst_amp  = [all_inst_amp,  inst_amp(valid)];
        end

        bin_edges   = fs_range(1):1:fs_range(2);
        bin_centers = bin_edges(1:end-1) + 0.5;
        counts = histcounts(all_inst_freq, bin_edges);
        counts = counts / sum(counts);

        % Amplitude-weighted histogram
        counts_wt = zeros(size(bin_centers));
        for bi = 1:length(bin_centers)
            in_bin = all_inst_freq >= bin_edges(bi) & all_inst_freq < bin_edges(bi+1);
            counts_wt(bi) = sum(all_inst_amp(in_bin));
        end
        if sum(counts_wt) > 0
            counts_wt = counts_wt / sum(counts_wt);
        end

        all_fs_hist_freqs{cond, subj}  = bin_centers;
        all_fs_hist_counts{cond, subj} = counts_wt;

        [~, mode_idx] = max(counts);
        all_fs_peakfreq(cond, subj) = bin_centers(mode_idx);
        all_fs_median(cond, subj)   = median(all_inst_freq);

        % Amplitude-weighted mean and peak
        all_fs_wt_mean(cond, subj) = sum(all_inst_freq .* all_inst_amp) / sum(all_inst_amp);
        [~, wt_mode_idx] = max(counts_wt);
        all_fs_wt_peak(cond, subj) = bin_centers(wt_mode_idx);

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
    %  PER-SUBJECT FIGURE
    %  ================================================================
    close all
    fig = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle(sprintf('GED Gamma (Common Filter): Subject %s', subjects{subj}), ...
        'FontSize', 20, 'FontWeight', 'bold');

    % --- Row 1: Raw + detrended power-ratio spectra (Method 1) ---
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

    % --- Row 2: Detrended power-ratio spectra with peak ---
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
        title(sprintf('%s Detrended', condLabels{cond}), 'FontSize', 12);
        set(gca, 'FontSize', 11); xlim([30 90]); grid on; box on;
    end

    % --- Row 3: Amplitude-weighted frequency sliding distributions (Method 2) ---
    for cond = 1:4
        subplot(4, 4, 8 + cond); hold on;
        if ~isempty(all_fs_hist_counts{cond, subj})
            bar(all_fs_hist_freqs{cond, subj}, all_fs_hist_counts{cond, subj}, 1, ...
                'FaceColor', colors(cond,:), 'FaceAlpha', 0.6, 'EdgeColor', 'none');
            xline(all_fs_wt_peak(cond, subj), '--', 'LineWidth', 2, ...
                'Color', colors(cond,:));
            text(all_fs_wt_peak(cond, subj) + 1, ...
                max(all_fs_hist_counts{cond, subj}) * 0.9, ...
                sprintf('%.0f Hz', all_fs_wt_peak(cond, subj)), ...
                'FontSize', 12, 'Color', colors(cond,:), 'FontWeight', 'bold');
        end
        xlabel('Freq [Hz]'); ylabel('P(f) [amp-wt]');
        title(sprintf('%s Freq Sliding', condLabels{cond}), 'FontSize', 12);
        set(gca, 'FontSize', 11); xlim([30 90]); grid on; box on;
    end

    % --- Row 4: Common topography + overlay plot ---
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

% --- Bottom left: Peak frequency scatter (Scanning) ---
subplot(2, 2, 3); hold on;
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
title('Individual Peak Frequencies (Power Scanning)', 'FontSize', 14);
ylim([30 90]); grid on; box on;

% --- Bottom right: Peak frequency scatter (Freq Sliding, amp-weighted) ---
subplot(2, 2, 4); hold on;
for cond = 1:4
    pf = all_fs_wt_peak(cond, :);
    pf = pf(~isnan(pf));
    scatter(ones(size(pf)) * cond + 0.15*(rand(size(pf))-0.5), pf, ...
        60, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.6);
    mu_pf  = nanmean(all_fs_wt_peak(cond,:));
    sem_pf = nanstd(all_fs_wt_peak(cond,:)) / sqrt(sum(~isnan(all_fs_wt_peak(cond,:))));
    errorbar(cond, mu_pf, sem_pf, 'k', 'LineWidth', 2, 'CapSize', 10);
    plot(cond, mu_pf, 'kd', 'MarkerSize', 12, ...
        'MarkerFaceColor', colors(cond,:), 'LineWidth', 1.5);
end
set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 13);
ylabel('Peak Gamma Frequency [Hz]');
title('Individual Peak Frequencies (Freq Sliding, Amp-Wt)', 'FontSize', 14);
ylim([30 90]); grid on; box on;

saveas(fig_ga1, fullfile(fig_save_dir, 'GED_scan_grand_average.png'));

%% ====================================================================
%  GRAND AVERAGE FIGURE 2: Method comparison
%  ====================================================================
fig_ga2 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle(sprintf('Method Comparison: Power Scanning vs Frequency Sliding (N=%d)', nSubj), ...
    'FontSize', 20, 'FontWeight', 'bold');

% --- Left: Correlation between methods ---
subplot(1, 2, 1); hold on;
all_scan_flat = all_scan_peakfreq(:);
all_fs_flat   = all_fs_wt_peak(:);
valid = ~isnan(all_scan_flat) & ~isnan(all_fs_flat);
for cond = 1:4
    scatter(all_scan_peakfreq(cond,:), all_fs_wt_peak(cond,:), ...
        80, colors(cond,:), 'filled', 'MarkerFaceAlpha', 0.7);
end
plot([30 90], [30 90], 'k--', 'LineWidth', 1.5);
if sum(valid) > 2
    [r, p_val] = corr(all_scan_flat(valid), all_fs_flat(valid), 'type', 'Spearman');
    text(32, 85, sprintf('r_s = %.2f, p = %.3f', r, p_val), 'FontSize', 14);
end
xlabel('Power Scanning Peak [Hz]'); ylabel('Freq Sliding Peak [Hz] (Amp-Wt)');
title('Agreement Between Methods', 'FontSize', 15);
set(gca, 'FontSize', 14); xlim([30 90]); ylim([30 90]);
axis square; grid on; box on;
legend(condLabels, 'Location', 'southeast', 'FontSize', 12);

% --- Right: Frequency sliding grand average distributions (amp-weighted) ---
subplot(1, 2, 2); hold on;
bin_centers_common = 40.5 : 1 : 79.5;
for cond = 1:4
    counts_mat = nan(nSubj, length(bin_centers_common));
    for s = 1:nSubj
        if ~isempty(all_fs_hist_counts{cond, s})
            counts_mat(s,:) = interp1(all_fs_hist_freqs{cond, s}, ...
                all_fs_hist_counts{cond, s}, bin_centers_common, 'linear', 0);
        end
    end
    mu = nanmean(counts_mat, 1);
    sem = nanstd(counts_mat, [], 1) / sqrt(sum(~isnan(counts_mat(:,1))));
    mu_s = movmean(mu, 3);
    sem_s = movmean(sem, 3);
    faceC = 0.8*colors(cond,:) + 0.2*[1 1 1];
    patch([bin_centers_common, fliplr(bin_centers_common)], ...
          [mu_s - sem_s, fliplr(mu_s + sem_s)], ...
          colors(cond,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot(bin_centers_common, mu_s, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);
end
xlabel('Frequency [Hz]'); ylabel('P(f)');
title('Grand Average Freq Sliding Distributions', 'FontSize', 15);
legend(condLabels, 'Location', 'northeast', 'FontSize', 12);
set(gca, 'FontSize', 14); xlim([30 90]); grid on; box on;

saveas(fig_ga2, fullfile(fig_save_dir, 'GED_method_comparison.png'));

%% ====================================================================
%  GRAND AVERAGE FIGURE 3: All subjects subplot
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
%  BOXPLOT FIGURE 1: Scanning peak frequencies across conditions
%  ====================================================================
fig_box1 = figure('Position', [0 0 1200 982], 'Color', 'w');
hold on;

for s = 1:nSubj
    pf = all_scan_peakfreq(:, s);
    if sum(~isnan(pf)) >= 2
        plot(1:4, pf, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    end
end

y_scan = all_scan_peakfreq(:);
g_scan = repelem((1:4)', nSubj, 1);
valid_scan = ~isnan(y_scan);
boxplot(y_scan(valid_scan), g_scan(valid_scan), 'Colors', 'k', 'Symbol', '');

hold on;
for c = 1:4
    pf = all_scan_peakfreq(c, :);
    pf = pf(~isnan(pf));
    xJit = c + (rand(size(pf)) - 0.5) * 0.1;
    scatter(xJit, pf, 250, colors(c,:), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

xlim([0.5 4.5]); ylim([30 90]);
set(gca, 'XTick', 1:4, 'XTickLabel', {'25%', '50%', '75%', '100%'}, 'FontSize', 20, 'Box', 'off');
ylabel('Peak Gamma Frequency [Hz]');
title('Power Scanning Peak Frequency', 'FontSize', 30, 'FontWeight', 'bold');

saveas(fig_box1, fullfile(fig_save_dir, 'GED_boxplot_peakfreq_GED.png'));

%% ====================================================================
%  BOXPLOT FIGURE 2: Frequency sliding peak frequencies across conditions
%  ====================================================================
fig_box2 = figure('Position', [0 0 1200 982], 'Color', 'w');
hold on;

for s = 1:nSubj
    pf = all_fs_wt_peak(:, s);
    if sum(~isnan(pf)) >= 2
        plot(1:4, pf, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    end
end

y_fs = all_fs_wt_peak(:);
g_fs = repelem((1:4)', nSubj, 1);
valid_fs = ~isnan(y_fs);
boxplot(y_fs(valid_fs), g_fs(valid_fs), 'Colors', 'k', 'Symbol', '');

hold on;
for c = 1:4
    pf = all_fs_wt_peak(c, :);
    pf = pf(~isnan(pf));
    xJit = c + (rand(size(pf)) - 0.5) * 0.1;
    scatter(xJit, pf, 250, colors(c,:), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

xlim([0.5 4.5]); ylim([30 90]);
set(gca, 'XTick', 1:4, 'XTickLabel', {'25%', '50%', '75%', '100%'}, 'FontSize', 20, 'Box', 'off');
ylabel('Peak Gamma Frequency [Hz]');
title('Frequency Sliding Peak Frequency (Amp-Weighted)', 'FontSize', 30, 'FontWeight', 'bold');

saveas(fig_box2, fullfile(fig_save_dir, 'GED_boxplot_peakfreq_FreqSliding.png'));

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
            xline(all_fs_wt_peak(cond, s), '--', 'Color', colors(cond,:), ...
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

    ga_pf = nanmean(all_fs_wt_peak(cond,:));
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
    'all_fs_hist_freqs', 'all_fs_hist_counts', ...
    'all_fs_peakfreq', 'all_fs_median', 'all_fs_wt_mean', 'all_fs_wt_peak', ...
    'all_powspec_stim', 'all_powspec_base', 'all_powspec_diff', 'all_powspec_freqs', ...
    'scan_freqs', 'subjects', 'condLabels', 'condNames');

clc
fprintf('Done.\n');
