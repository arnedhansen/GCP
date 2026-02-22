%% GCP GED Gamma Metrics — Power, Shape, and Coherence
%
% Extracts gamma metrics beyond peak frequency from GED-optimised data.
% Peak frequency is condition-invariant in human scalp EEG; this script
% tests whether contrast modulates gamma POWER, spectral SHAPE, and
% inter-electrode COHERENCE (Roberts et al. 2013).
%
% Part 1-2:  Load saved peak explorer data → power, spectral shape, and
%            count metrics (fast, no raw-data processing)
% Part 3:    Load raw EEG → recompute GED → inter-electrode coherence
%            and wPLI between occipital and parietal ROIs
% Part 4:    Visualisation (7 figures)
%
% Dependencies:
%   - GCP_eeg_GED_peak_explorer.mat  (from GCP_eeg_ged_peak_explorer.m)
%   - Raw EEG data (dataEEG)         (only if compute_coherence = true)
%   - FieldTrip on path

clear; close all; clc

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');

nSubj = length(subjects);

baseline_window = [-1.5, -0.25];
stimulus_window = [0.3, 2.0];
gamma_range     = [30, 90];
scan_freqs      = 30:1:90;
nFreqs          = length(scan_freqs);
poly_order      = 2;
lambda          = 0.01;

condNames  = {'c25', 'c50', 'c75', 'c100'};
condLabels = {'25%', '50%', '75%', '100%'};
contrast_vals = [25, 50, 75, 100];

compute_coherence = true;

if ispc
    fig_save_dir  = 'W:\Students\Arne\GCP\figures\eeg\ged\gamma_metrics';
    data_save_dir = 'W:\Students\Arne\GCP\data\features';
else
    fig_save_dir  = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/eeg/ged/gamma_metrics';
    data_save_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features';
end
if ~exist(fig_save_dir, 'dir'), mkdir(fig_save_dir); end

%% Load saved peak explorer data
explorer_path = fullfile(data_save_dir, 'GCP_eeg_GED_peak_explorer.mat');
fprintf('Loading peak explorer data from:\n  %s\n', explorer_path);
load(explorer_path);
fprintf('Loaded: %d subjects, freqs %d-%d Hz\n', nSubj, scan_freqs(1), scan_freqs(end));

%% ====================================================================
%  PART 1-2: Power, Spectral Shape, and Count Metrics (from saved data)
%  ====================================================================

% --- Preallocate subject-level matrices [4 x nSubj] ---
metric_peak_amp   = nan(4, nSubj);
metric_bb_power   = nan(4, nSubj);
metric_lo_power   = nan(4, nSubj);
metric_hi_power   = nan(4, nSubj);
metric_auc        = nan(4, nSubj);
metric_fwhm       = nan(4, nSubj);
metric_prominence = nan(4, nSubj);
metric_centroid   = nan(4, nSubj);
metric_det_rate   = nan(4, nSubj);
metric_peak_count = nan(4, nSubj);
metric_lohi_ratio = nan(4, nSubj);

% --- Trial-level storage ---
trl_peak_amp   = cell(4, nSubj);
trl_bb_power   = cell(4, nSubj);
trl_lo_power   = cell(4, nSubj);
trl_hi_power   = cell(4, nSubj);
trl_auc        = cell(4, nSubj);
trl_fwhm       = cell(4, nSubj);
trl_prominence = cell(4, nSubj);
trl_centroid   = cell(4, nSubj);
trl_peak_freq  = cell(4, nSubj);

lo_mask = scan_freqs <= 49;
hi_mask = scan_freqs >= 50;

for subj = 1:nSubj
    for cond = 1:4
        pr_mat = all_trial_powratio{cond, subj};
        if isempty(pr_mat), continue; end

        nTrl = size(pr_mat, 1);

        t_peak_amp  = nan(nTrl, 1);
        t_peak_freq = nan(nTrl, 1);
        t_bb_pow    = nan(nTrl, 1);
        t_lo_pow    = nan(nTrl, 1);
        t_hi_pow    = nan(nTrl, 1);
        t_auc       = nan(nTrl, 1);
        t_fwhm      = nan(nTrl, 1);
        t_prom      = nan(nTrl, 1);
        t_centroid  = nan(nTrl, 1);

        for trl = 1:nTrl
            pr = pr_mat(trl, :);
            if all(isnan(pr)), continue; end

            % --- Power metrics (raw power ratio) ---
            t_bb_pow(trl) = nanmean(pr);
            t_lo_pow(trl) = nanmean(pr(lo_mask));
            t_hi_pow(trl) = nanmean(pr(hi_mask));

            % --- Detrend + smooth ---
            p = polyfit(scan_freqs, pr, poly_order);
            pr_dt   = pr - polyval(p, scan_freqs);
            pr_dt_s = movmean(pr_dt, 5);

            % --- Peak detection ---
            [pks, locs, ~, proms] = findpeaks(pr_dt_s, scan_freqs, ...
                'MinPeakDistance', 5);
            pos = pks > 0;

            if any(pos)
                pks_p   = pks(pos);
                locs_p  = locs(pos);
                proms_p = proms(pos);
                [max_pk, bi] = max(pks_p);
                peak_freq = locs_p(bi);

                t_peak_amp(trl)  = max_pk;
                t_peak_freq(trl) = peak_freq;
                t_prom(trl)      = proms_p(bi);

                % FWHM via linear interpolation
                [~, peak_fi] = min(abs(scan_freqs - peak_freq));
                half_h = max_pk / 2;

                left_f = NaN;
                for ki = peak_fi:-1:2
                    if pr_dt_s(ki) >= half_h && pr_dt_s(ki-1) < half_h
                        frac   = (half_h - pr_dt_s(ki-1)) / (pr_dt_s(ki) - pr_dt_s(ki-1));
                        left_f = scan_freqs(ki-1) + frac;
                        break;
                    end
                end
                right_f = NaN;
                for ki = peak_fi:nFreqs-1
                    if pr_dt_s(ki) >= half_h && pr_dt_s(ki+1) < half_h
                        frac    = (pr_dt_s(ki) - half_h) / (pr_dt_s(ki) - pr_dt_s(ki+1));
                        right_f = scan_freqs(ki) + frac;
                        break;
                    end
                end
                if ~isnan(left_f) && ~isnan(right_f)
                    t_fwhm(trl) = right_f - left_f;
                end
            end

            % --- AUC above zero ---
            pos_spec   = max(pr_dt_s, 0);
            t_auc(trl) = trapz(scan_freqs, pos_spec);

            % --- Spectral centroid ---
            if any(pos_spec > 0)
                t_centroid(trl) = sum(scan_freqs .* pos_spec) / sum(pos_spec);
            end
        end

        % Store trial-level
        trl_peak_amp{cond, subj}   = t_peak_amp;
        trl_peak_freq{cond, subj}  = t_peak_freq;
        trl_bb_power{cond, subj}   = t_bb_pow;
        trl_lo_power{cond, subj}   = t_lo_pow;
        trl_hi_power{cond, subj}   = t_hi_pow;
        trl_auc{cond, subj}        = t_auc;
        trl_fwhm{cond, subj}       = t_fwhm;
        trl_prominence{cond, subj} = t_prom;
        trl_centroid{cond, subj}   = t_centroid;

        % Aggregate to subject level (mean across trials)
        metric_peak_amp(cond, subj)   = nanmean(t_peak_amp);
        metric_bb_power(cond, subj)   = nanmean(t_bb_pow);
        metric_lo_power(cond, subj)   = nanmean(t_lo_pow);
        metric_hi_power(cond, subj)   = nanmean(t_hi_pow);
        metric_auc(cond, subj)        = nanmean(t_auc);
        metric_fwhm(cond, subj)       = nanmean(t_fwhm);
        metric_prominence(cond, subj) = nanmean(t_prom);
        metric_centroid(cond, subj)   = nanmean(t_centroid);

        % Count/rate metrics
        pc = all_peak_counts{cond, subj};
        if ~isempty(pc)
            metric_det_rate(cond, subj)   = mean(pc > 0);
            metric_peak_count(cond, subj) = mean(pc);
        end

        % Low/high gamma ratio (fraction of peaks in low gamma)
        tpf = all_peak_freqs{cond, subj};
        if ~isempty(tpf)
            all_pf = [tpf{:}];
            if ~isempty(all_pf)
                metric_lohi_ratio(cond, subj) = sum(all_pf <= 50) / length(all_pf);
            end
        end
    end
    clc
    fprintf('Part 1-2: Subject %s (%d/%d) done\n', subjects{subj}, subj, nSubj);
end

fprintf('\nPart 1-2 complete. %d power/shape/count metrics extracted.\n', 11);

%% ====================================================================
%  PART 3: Coherence (raw EEG → GED → ft_connectivityanalysis)
%  ====================================================================

% Preallocate coherence storage
coh_occ_par   = nan(4, nSubj);   % mean gamma coherence: occ-par
wpli_occ_par  = nan(4, nSubj);   % mean gamma wPLI: occ-par
coh_interhemi = nan(4, nSubj);   % mean gamma coherence: left-right occ
wpli_interhemi = nan(4, nSubj);

coh_spectra_occ_par    = [];     % will be set after first subject (freq axis needed)
coh_spectra_interhemi  = [];
wpli_spectra_occ_par   = [];
wpli_spectra_interhemi = [];
coh_freq_axis          = [];

all_topComps = cell(1, nSubj);   % save GED filters for reuse

if compute_coherence

    % Define parietal ROI (must NOT overlap with occipital GED channels)
    par_candidates = {'Pz', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'};

    for subj = 1:nSubj
        tic
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        load dataEEG

        fsample = dataEEG_c25.fsample;
        labels  = dataEEG_c25.label;
        nChans  = length(labels);

        % --- Identify channel ROIs ---
        occ_mask = cellfun(@(l) ~isempty(regexp(l, '[OI]', 'once')), labels);
        occ_idx  = find(occ_mask);
        nOcc     = length(occ_idx);

        par_idx = find(ismember(labels, par_candidates));
        nPar    = length(par_idx);

        occ_labels = labels(occ_idx);
        left_occ  = [];
        right_occ = [];
        for ch = 1:length(occ_labels)
            nums = regexp(occ_labels{ch}, '\d+', 'match');
            if ~isempty(nums) && mod(str2double(nums{end}), 2) == 1
                left_occ = [left_occ, occ_idx(ch)];
            elseif ~isempty(nums) && mod(str2double(nums{end}), 2) == 0
                right_occ = [right_occ, occ_idx(ch)];
            end
        end

        % --- Recompute GED (Phase 1 from peak explorer) ---
        clc
        fprintf('Part 3 Coherence: Subject %s (%d/%d) — GED\n', ...
            subjects{subj}, subj, nSubj);

        trialIndices = { ...
            find(dataEEG_c25.trialinfo  == 61), ...
            find(dataEEG_c50.trialinfo  == 62), ...
            find(dataEEG_c75.trialinfo  == 63), ...
            find(dataEEG_c100.trialinfo == 64)};
        dataStructs = {dataEEG_c25, dataEEG_c50, dataEEG_c75, dataEEG_c100};

        covStim_full = zeros(nChans);
        covBase_full = zeros(nChans);
        nTrials_total = 0;
        dat_per_cond = cell(1, 4);

        for cond = 1:4
            dat    = dataStructs{cond};
            trlIdx = trialIndices{cond};
            if isempty(trlIdx), continue; end

            cfg = [];
            cfg.trials = trlIdx;
            dat = ft_selectdata(cfg, dat);
            dat_per_cond{cond} = dat;

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
                covStim_full = covStim_full + (d * d') / size(d, 2);

                d = double(dat_base.trial{trl});
                d = bsxfun(@minus, d, mean(d, 2));
                covBase_full = covBase_full + (d * d') / size(d, 2);
            end
            nTrials_total = nTrials_total + nTrl;
        end

        covStim_full = covStim_full / nTrials_total;
        covBase_full = covBase_full / nTrials_total;

        covStim_occ = covStim_full(occ_idx, occ_idx);
        covBase_occ = covBase_full(occ_idx, occ_idx);
        covStim_occ = (1-lambda)*covStim_occ + lambda*mean(diag(covStim_occ))*eye(nOcc);
        covBase_occ = (1-lambda)*covBase_occ + lambda*mean(diag(covBase_occ))*eye(nOcc);

        [W_occ, D_occ] = eig(covStim_occ, covBase_occ);
        [~, sortIdx] = sort(real(diag(D_occ)), 'descend');
        W_occ = W_occ(:, sortIdx);
        w_occ = W_occ(:, 1);

        topComp = zeros(nChans, 1);
        topComp(occ_idx) = w_occ;

        covStim_full_reg = (1-lambda)*covStim_full + lambda*mean(diag(covStim_full))*eye(nChans);
        topo_temp = covStim_full_reg * topComp;
        [~, mxI] = max(abs(topo_temp));
        if topo_temp(mxI) < 0, topComp = -topComp; end

        all_topComps{subj} = topComp;

        % --- Coherence per condition ---
        for cond = 1:4
            dat = dat_per_cond{cond};
            if isempty(dat), continue; end

            clc
            fprintf('Part 3 Coherence: Subject %s (%d/%d) — Cond %d/4\n', ...
                subjects{subj}, subj, nSubj, cond);

            cfg_t = [];
            cfg_t.latency = stimulus_window;
            dat_stim = ft_selectdata(cfg_t, dat);

            % Fourier analysis (multitaper)
            cfg_freq = [];
            cfg_freq.method     = 'mtmfft';
            cfg_freq.output     = 'fourier';
            cfg_freq.keeptrials = 'yes';
            cfg_freq.foilim     = [30 90];
            cfg_freq.tapsmofrq  = 4;
            cfg_freq.taper      = 'dpss';
            cfg_freq.pad        = 'nextpow2';
            freq = ft_freqanalysis(cfg_freq, dat_stim);

            % Inject GED component as virtual channel
            nRpt  = size(freq.fourierspctrm, 1);
            nFreq = size(freq.fourierspctrm, 3);
            ged_four = sum(freq.fourierspctrm .* reshape(topComp, 1, [], 1), 2);
            freq.fourierspctrm = cat(2, freq.fourierspctrm, ged_four);
            freq.label{end+1}  = 'GED_gamma';

            % Set up frequency axis on first pass
            if isempty(coh_freq_axis)
                coh_freq_axis = freq.freq;
                nCohFreq = length(coh_freq_axis);
                coh_spectra_occ_par    = nan(4, nSubj, nCohFreq);
                coh_spectra_interhemi  = nan(4, nSubj, nCohFreq);
                wpli_spectra_occ_par   = nan(4, nSubj, nCohFreq);
                wpli_spectra_interhemi = nan(4, nSubj, nCohFreq);
            end

            ged_chan_idx   = find(strcmp(freq.label, 'GED_gamma'));
            par_chan_idx   = find(ismember(freq.label, par_candidates));
            left_chan_idx  = find(ismember(freq.label, labels(left_occ)));
            right_chan_idx = find(ismember(freq.label, labels(right_occ)));

            % --- Coherence ---
            cfg_conn = [];
            cfg_conn.method  = 'coh';
            cfg_conn.complex = 'abs';
            conn_coh = ft_connectivityanalysis(cfg_conn, freq);

            % GED-parietal coherence spectrum
            if ~isempty(par_chan_idx)
                coh_gp = squeeze(mean(conn_coh.cohspctrm(ged_chan_idx, par_chan_idx, :), 2));
                coh_spectra_occ_par(cond, subj, :) = coh_gp;
                coh_occ_par(cond, subj) = mean(coh_gp);
            end

            % Inter-hemispheric coherence spectrum
            if ~isempty(left_chan_idx) && ~isempty(right_chan_idx)
                coh_lr = conn_coh.cohspctrm(left_chan_idx, right_chan_idx, :);
                coh_lr_mean = squeeze(mean(mean(coh_lr, 1), 2));
                coh_spectra_interhemi(cond, subj, :) = coh_lr_mean;
                coh_interhemi(cond, subj) = mean(coh_lr_mean);
            end

            % --- wPLI (debiased, robust to volume conduction) ---
            cfg_conn.method = 'wpli_debiased';
            conn_wpli = ft_connectivityanalysis(cfg_conn, freq);

            if ~isempty(par_chan_idx)
                wpli_gp = squeeze(mean(conn_wpli.wpli_debiasedspctrm(ged_chan_idx, par_chan_idx, :), 2));
                wpli_spectra_occ_par(cond, subj, :) = wpli_gp;
                wpli_occ_par(cond, subj) = mean(wpli_gp);
            end

            if ~isempty(left_chan_idx) && ~isempty(right_chan_idx)
                wpli_lr = conn_wpli.wpli_debiasedspctrm(left_chan_idx, right_chan_idx, :);
                wpli_lr_mean = squeeze(mean(mean(wpli_lr, 1), 2));
                wpli_spectra_interhemi(cond, subj, :) = wpli_lr_mean;
                wpli_interhemi(cond, subj) = mean(wpli_lr_mean);
            end
        end

        fprintf('  Coherence done (%.1f s)\n', toc);
    end
end

fprintf('\nAll metrics computed.\n');

%% ====================================================================
%  VISUALISATION
%  ====================================================================
close all

%% --- Figure 1: Power Metrics (2x2 rainclouds) ---
fig1 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Power Metrics by Contrast Level', 'FontSize', 18, 'FontWeight', 'bold');

metric_data   = {metric_peak_amp, metric_bb_power, metric_lo_power, metric_hi_power};
metric_titles = {'Peak Amplitude (\Delta PR)', 'Broadband Power (ratio)', ...
                 'Low-Gamma Power (30-49 Hz)', 'High-Gamma Power (50-90 Hz)'};

for mi = 1:4
    subplot(2, 2, mi);
    plot_raincloud(metric_data{mi}, colors, condLabels, metric_titles{mi});
end

saveas(fig1, fullfile(fig_save_dir, 'GCP_ged_metrics_power.png'));

%% --- Figure 2: Spectral Shape Metrics (2x2 rainclouds) ---
fig2 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Spectral Shape Metrics by Contrast Level', 'FontSize', 18, 'FontWeight', 'bold');

shape_data   = {metric_auc, metric_fwhm, metric_prominence, metric_centroid};
shape_titles = {'AUC Above Zero', 'FWHM [Hz]', 'Peak Prominence', 'Spectral Centroid [Hz]'};

for mi = 1:4
    subplot(2, 2, mi);
    plot_raincloud(shape_data{mi}, colors, condLabels, shape_titles{mi});
end

saveas(fig2, fullfile(fig_save_dir, 'GCP_ged_metrics_shape.png'));

%% --- Figure 3: Count Metrics (1x3 bar+scatter) ---
fig3 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Peak Count / Detection Metrics by Contrast Level', ...
    'FontSize', 18, 'FontWeight', 'bold');

count_data   = {metric_det_rate * 100, metric_peak_count, metric_lohi_ratio * 100};
count_titles = {'Detection Rate [%]', 'Mean Peak Count', 'Low-Gamma Peak Fraction [%]'};
count_ylims  = {[0 105], [], [0 105]};

for mi = 1:3
    subplot(1, 3, mi); hold on;
    dat = count_data{mi};

    mu  = nanmean(dat, 2);
    sem = nanstd(dat, [], 2) / sqrt(nSubj);

    b = bar(1:4, mu, 0.6);
    b.FaceColor = 'flat';
    for c = 1:4, b.CData(c,:) = colors(c,:); end
    errorbar(1:4, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 10);

    for s = 1:nSubj
        plot(1:4, dat(:, s), '-', 'Color', [0.6 0.6 0.6 0.4], 'LineWidth', 0.8);
        for c = 1:4
            if ~isnan(dat(c, s))
                scatter(c + (rand-0.5)*0.25, dat(c, s), 40, colors(c,:), ...
                    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
                    'MarkerFaceAlpha', 0.7);
            end
        end
    end

    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 14);
    ylabel(count_titles{mi}, 'FontSize', 13);
    title(count_titles{mi}, 'FontSize', 14, 'FontWeight', 'bold');
    if ~isempty(count_ylims{mi}), ylim(count_ylims{mi}); end
    grid on; box on;
end

saveas(fig3, fullfile(fig_save_dir, 'GCP_ged_metrics_counts.png'));

%% --- Figure 4: Coherence Spectra — Roberts et al. style ---
if compute_coherence && ~isempty(coh_freq_axis)

    fig4 = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle('Gamma Coherence Spectra by Contrast Level', ...
        'FontSize', 18, 'FontWeight', 'bold');

    spec_data   = {coh_spectra_occ_par, coh_spectra_interhemi, ...
                   wpli_spectra_occ_par, wpli_spectra_interhemi};
    spec_titles = {'Coherence: GED comp \leftrightarrow Parietal', ...
                   'Coherence: Left Occ \leftrightarrow Right Occ', ...
                   'wPLI: GED comp \leftrightarrow Parietal', ...
                   'wPLI: Left Occ \leftrightarrow Right Occ'};
    spec_ylabels = {'Coherence', 'Coherence', 'wPLI (debiased)', 'wPLI (debiased)'};

    for si = 1:4
        subplot(2, 2, si); hold on;
        for cond = 1:4
            spec = squeeze(spec_data{si}(cond, :, :));  % [nSubj x nFreq]
            mu  = nanmean(spec, 1);
            sem = nanstd(spec, [], 1) / sqrt(nSubj);

            fill([coh_freq_axis, fliplr(coh_freq_axis)], ...
                [mu - sem, fliplr(mu + sem)], ...
                colors(cond,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            plot(coh_freq_axis, mu, '-', 'Color', colors(cond,:), 'LineWidth', 2.5);
        end
        xlabel('Frequency [Hz]', 'FontSize', 12);
        ylabel(spec_ylabels{si}, 'FontSize', 12);
        title(spec_titles{si}, 'FontSize', 13, 'FontWeight', 'bold');
        set(gca, 'FontSize', 11); xlim([30 90]); grid on; box on;
        if si == 1
            legend(condLabels, 'FontSize', 10, 'Location', 'best');
        end
    end

    saveas(fig4, fullfile(fig_save_dir, 'GCP_ged_metrics_coherence_spectra.png'));

    %% --- Figure 5: Coherence Rainclouds ---
    fig5 = figure('Position', [0 0 1512 982], 'Color', 'w');
    sgtitle('Mean Gamma-Band Coherence by Contrast Level', ...
        'FontSize', 18, 'FontWeight', 'bold');

    coh_rc_data   = {coh_occ_par, coh_interhemi, wpli_occ_par, wpli_interhemi};
    coh_rc_titles = {'Coherence: GED\leftrightarrowParietal', ...
                     'Coherence: L-Occ\leftrightarrowR-Occ', ...
                     'wPLI: GED\leftrightarrowParietal', ...
                     'wPLI: L-Occ\leftrightarrowR-Occ'};

    for ci = 1:4
        subplot(2, 2, ci);
        plot_raincloud(coh_rc_data{ci}, colors, condLabels, coh_rc_titles{ci});
    end

    saveas(fig5, fullfile(fig_save_dir, 'GCP_ged_metrics_coherence_rainclouds.png'));
end

%% --- Figure 6: Contrast Response Functions ---
fig6 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Contrast Response Functions', 'FontSize', 18, 'FontWeight', 'bold');

if compute_coherence && ~isempty(coh_freq_axis)
    crf_data   = {metric_bb_power, metric_peak_amp, metric_auc, ...
                  coh_occ_par, wpli_occ_par, metric_centroid};
    crf_titles = {'Broadband Power', 'Peak Amplitude', 'AUC Above Zero', ...
                  'Coherence (GED-Par)', 'wPLI (GED-Par)', 'Spectral Centroid [Hz]'};
    nCRF = 6;
    nCols = 3; nRows = 2;
else
    crf_data   = {metric_bb_power, metric_peak_amp, metric_auc, ...
                  metric_lo_power, metric_hi_power, metric_centroid};
    crf_titles = {'Broadband Power', 'Peak Amplitude', 'AUC Above Zero', ...
                  'Low-Gamma Power', 'High-Gamma Power', 'Spectral Centroid [Hz]'};
    nCRF = 6;
    nCols = 3; nRows = 2;
end

for ci = 1:nCRF
    subplot(nRows, nCols, ci); hold on;
    dat = crf_data{ci};

    for s = 1:nSubj
        plot(contrast_vals, dat(:, s), '-o', 'Color', [0.7 0.7 0.7 0.5], ...
            'MarkerSize', 4, 'MarkerFaceColor', [0.7 0.7 0.7], 'LineWidth', 0.8);
    end

    mu  = nanmean(dat, 2);
    sem = nanstd(dat, [], 2) / sqrt(nSubj);
    errorbar(contrast_vals, mu, sem, 'k-o', 'LineWidth', 2.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', 'k', 'CapSize', 8);

    xlabel('Contrast [%]', 'FontSize', 12);
    ylabel(crf_titles{ci}, 'FontSize', 12);
    title(crf_titles{ci}, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'XTick', contrast_vals, 'FontSize', 11);
    xlim([15 110]); grid on; box on;
end

saveas(fig6, fullfile(fig_save_dir, 'GCP_ged_metrics_CRF.png'));

%% --- Figure 7: Summary Dashboard ---
fig7 = figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Gamma Metrics Summary — All Measures', ...
    'FontSize', 16, 'FontWeight', 'bold');

if compute_coherence && ~isempty(coh_freq_axis)
    all_metrics = {metric_peak_amp, metric_bb_power, metric_lo_power, ...
                   metric_hi_power, metric_auc, metric_fwhm, ...
                   metric_prominence, metric_centroid, ...
                   metric_det_rate*100, metric_peak_count, metric_lohi_ratio*100, ...
                   coh_occ_par, wpli_occ_par};
    all_names = {'Peak Amp', 'BB Power', 'Lo-\gamma Power', 'Hi-\gamma Power', ...
                 'AUC', 'FWHM [Hz]', 'Prominence', 'Centroid [Hz]', ...
                 'Det Rate [%]', 'Peak Count', 'Lo-\gamma Frac [%]', ...
                 'Coh (GED-Par)', 'wPLI (GED-Par)'};
    nMetrics = 13;
    nCols_d = 5; nRows_d = 3;
else
    all_metrics = {metric_peak_amp, metric_bb_power, metric_lo_power, ...
                   metric_hi_power, metric_auc, metric_fwhm, ...
                   metric_prominence, metric_centroid, ...
                   metric_det_rate*100, metric_peak_count, metric_lohi_ratio*100};
    all_names = {'Peak Amp', 'BB Power', 'Lo-\gamma Power', 'Hi-\gamma Power', ...
                 'AUC', 'FWHM [Hz]', 'Prominence', 'Centroid [Hz]', ...
                 'Det Rate [%]', 'Peak Count', 'Lo-\gamma Frac [%]'};
    nMetrics = 11;
    nCols_d = 4; nRows_d = 3;
end

for mi = 1:nMetrics
    subplot(nRows_d, nCols_d, mi); hold on;
    dat = all_metrics{mi};

    mu  = nanmean(dat, 2);
    sem = nanstd(dat, [], 2) / sqrt(nSubj);

    for c = 1:4
        bar(c, mu(c), 0.6, 'FaceColor', colors(c,:), 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    end
    errorbar(1:4, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 6);

    for s = 1:nSubj
        for c = 1:4
            if ~isnan(dat(c, s))
                scatter(c + (rand-0.5)*0.2, dat(c, s), 20, [0.3 0.3 0.3], ...
                    'filled', 'MarkerFaceAlpha', 0.5);
            end
        end
    end

    set(gca, 'XTick', 1:4, 'XTickLabel', {'25', '50', '75', '100'}, 'FontSize', 9);
    title(all_names{mi}, 'FontSize', 10, 'FontWeight', 'bold');
    box on;
end

saveas(fig7, fullfile(fig_save_dir, 'GCP_ged_metrics_summary.png'));

%% ====================================================================
%  TRIAL-LEVEL LONG-FORMAT TABLE (for R / lme4)
%  ====================================================================
fprintf('Building trial-level long-format table...\n');

trl_metric_names = {'PeakFrequency', 'PeakAmplitude', 'BroadbandPower', ...
                    'LowGammaPower', 'HighGammaPower', 'AUC', 'FWHM', ...
                    'Prominence', 'SpectralCentroid'};
trl_metric_cells = {trl_peak_freq, trl_peak_amp, trl_bb_power, trl_lo_power, ...
                    trl_hi_power, trl_auc, trl_fwhm, trl_prominence, trl_centroid};

row_Subject   = [];
row_Condition = [];
row_Contrast  = [];
row_Trial     = [];
row_metrics   = nan(0, length(trl_metric_names));

for subj = 1:nSubj
    for cond = 1:4
        nTrl = length(trl_peak_amp{cond, subj});
        if nTrl == 0, continue; end

        row_Subject   = [row_Subject;   repmat(str2double(subjects{subj}), nTrl, 1)];
        row_Condition = [row_Condition;  repmat(cond, nTrl, 1)];
        row_Contrast  = [row_Contrast;   repmat(contrast_vals(cond), nTrl, 1)];
        row_Trial     = [row_Trial;      (1:nTrl)'];

        block = nan(nTrl, length(trl_metric_names));
        for mi = 1:length(trl_metric_names)
            block(:, mi) = trl_metric_cells{mi}{cond, subj};
        end
        row_metrics = [row_metrics; block];
    end
end

% Also add trial-level peak count from peak explorer data
row_PeakCount = [];
for subj = 1:nSubj
    for cond = 1:4
        pc = all_peak_counts{cond, subj};
        if isempty(pc), continue; end
        row_PeakCount = [row_PeakCount; pc(:)];
    end
end

% Also add low/high ratio per trial (fraction of peaks <= 50 Hz)
row_LoHiRatio = [];
for subj = 1:nSubj
    for cond = 1:4
        tpf = all_peak_freqs{cond, subj};
        if isempty(tpf), continue; end
        nTrl = length(tpf);
        ratio_trl = nan(nTrl, 1);
        for trl = 1:nTrl
            pf = tpf{trl};
            if ~isempty(pf)
                ratio_trl(trl) = sum(pf <= 50) / length(pf);
            end
        end
        row_LoHiRatio = [row_LoHiRatio; ratio_trl];
    end
end

% Build table
T = table(row_Subject, row_Condition, row_Contrast, row_Trial, ...
    'VariableNames', {'Subject', 'Condition', 'Contrast', 'Trial'});

for mi = 1:length(trl_metric_names)
    T.(trl_metric_names{mi}) = row_metrics(:, mi);
end
T.PeakCount  = row_PeakCount;
T.LoHiRatio  = row_LoHiRatio;

fprintf('  Table: %d rows x %d columns\n', height(T), width(T));
fprintf('  Rows per condition: ');
for cond = 1:4
    fprintf('%s=%d  ', condLabels{cond}, sum(T.Condition == cond));
end
fprintf('\n');

% Save as CSV
if ispc
    csv_path = fullfile('W:\Students\Arne\GCP\data\features', 'GCP_eeg_GED_gamma_metrics_trials.csv');
else
    csv_path = fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features', 'GCP_eeg_GED_gamma_metrics_trials.csv');
end
writetable(T, csv_path);
fprintf('  CSV saved to:\n    %s\n', csv_path);

% Also save as .mat
if ispc
    tbl_mat_path = fullfile('W:\Students\Arne\GCP\data\features', 'GCP_eeg_GED_gamma_metrics_trials.mat');
else
    tbl_mat_path = fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features', 'GCP_eeg_GED_gamma_metrics_trials.mat');
end
save(tbl_mat_path, 'T');
fprintf('  MAT saved to:\n    %s\n', tbl_mat_path);

%% ====================================================================
%  SAVE (subject-level + coherence)
%  ====================================================================
save_path = fullfile(data_save_dir, 'GCP_eeg_GED_gamma_metrics.mat');
fprintf('Saving subject-level data to:\n  %s\n', save_path);

save(save_path, ...
    'metric_peak_amp', 'metric_bb_power', 'metric_lo_power', 'metric_hi_power', ...
    'metric_auc', 'metric_fwhm', 'metric_prominence', 'metric_centroid', ...
    'metric_det_rate', 'metric_peak_count', 'metric_lohi_ratio', ...
    'trl_peak_freq', 'trl_peak_amp', 'trl_bb_power', 'trl_lo_power', ...
    'trl_hi_power', 'trl_auc', 'trl_fwhm', 'trl_prominence', 'trl_centroid', ...
    'coh_occ_par', 'wpli_occ_par', 'coh_interhemi', 'wpli_interhemi', ...
    'coh_spectra_occ_par', 'coh_spectra_interhemi', ...
    'wpli_spectra_occ_par', 'wpli_spectra_interhemi', ...
    'coh_freq_axis', 'all_topComps', ...
    'scan_freqs', 'subjects', 'condLabels', 'contrast_vals');

clc
fprintf('GCP GED Gamma Metrics — Done.\n');
fprintf('Figures saved to:\n  %s\n', fig_save_dir);
fprintf('Subject-level data saved to:\n  %s\n', save_path);
fprintf('Trial-level CSV saved to:\n  %s\n', csv_path);
fprintf('Trial-level MAT saved to:\n  %s\n', tbl_mat_path);

%% ====================================================================
%  HELPER FUNCTION: Raincloud subplot
%  ====================================================================
function plot_raincloud(data_mat, colors, condLabels, title_str)
    hold on;
    nSubj = size(data_mat, 2);

    % Subject lines
    for s = 1:nSubj
        pf = data_mat(:, s);
        if sum(~isnan(pf)) >= 2
            plot(1:4, pf, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.8);
        end
    end

    % Violin (density) on left side
    for c = 1:4
        vals = data_mat(c, :);
        vals = vals(~isnan(vals));
        if length(vals) >= 3
            [f_dens, xi] = ksdensity(vals);
            f_dens = f_dens / max(f_dens) * 0.3;
            patch(c - f_dens - 0.05, xi, colors(c,:), ...
                'FaceAlpha', 0.3, 'EdgeColor', colors(c,:), 'LineWidth', 1);
        end
    end

    % Boxplot
    y_box = data_mat(:);
    g_box = repelem((1:4)', nSubj, 1);
    valid = ~isnan(y_box);
    if any(valid)
        boxplot(y_box(valid), g_box(valid), 'Colors', 'k', ...
            'Symbol', '', 'Widths', 0.15);
    end

    % Scatter (right side)
    hold on;
    for c = 1:4
        vals = data_mat(c, :);
        vals_valid = vals(~isnan(vals));
        xJit = c + 0.15 + (rand(size(vals_valid)) - 0.5) * 0.12;
        scatter(xJit, vals_valid, 180, colors(c,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end

    set(gca, 'XTick', 1:4, 'XTickLabel', condLabels, 'FontSize', 13, 'Box', 'off');
    title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
end
