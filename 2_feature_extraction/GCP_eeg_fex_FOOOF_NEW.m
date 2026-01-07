%% GCP EEG FOOOF
%
% Outputs per subject:
%   TFR (mtmconvol) raw + baseline (for visualisation)
%   TFR_FOOOF (sliding windows): powspctrm = model fit (aperiodic + peaks) in FOOOF/log space
%   plus: power_spectrum + aperiodic params (offset/exponent/error/r^2) per window

clear; close all; clc

%% Setup
startup
[subjects, path, ~, ~] = setup('GCP');

%% DataQueue for parfor progress
% function printProgress(s)
% fprintf('Subj %d | Cond %d | Time %d/%d finished\n', ...
%     s.subj, s.cond, s.time, s.nTimePnts);
% end
% D = parallel.pool.DataQueue;
% afterEach(D, @printProgress);

%% Loop subjects
for subj = 1:length(subjects)

    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    close all
    clc

    disp(['Processing GCP TFR + sliding-window FOOOF for subject ', num2str(subjects{subj})])

    load dataEEG
    if ispc
        load('W:\Students\Arne\MA\headmodel\ant128lay.mat');
    else
        load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/ant128lay.mat');
    end

    %% Identify trial indices (per dataset)
    ind61 = find(dataEEG_c25.trialinfo  == 61);
    ind62 = find(dataEEG_c50.trialinfo  == 62);
    ind63 = find(dataEEG_c75.trialinfo  == 63);
    ind64 = find(dataEEG_c100.trialinfo == 64);

    %% TFR (mtmconvol) for visualisation / later use
    cfg             = [];
    cfg.output      = 'pow';
    cfg.method      = 'mtmconvol';
    cfg.taper       = 'dpss';
    cfg.foi         = 1:1:100;
    cfg.tapsmofrq   = 5;
    cfg.t_ftimwin   = ones(length(cfg.foi),1) .* 0.5;
    cfg.toi         = -1.75:0.05:2;
    cfg.keeptrials  = 'no';
    cfg.pad         = 'nextpow2';

    cfg.trials = ind61;  tfr_c25  = ft_freqanalysis(cfg, dataEEG_c25);
    cfg.trials = ind62;  tfr_c50  = ft_freqanalysis(cfg, dataEEG_c50);
    cfg.trials = ind63;  tfr_c75  = ft_freqanalysis(cfg, dataEEG_c75);
    cfg.trials = ind64;  tfr_c100 = ft_freqanalysis(cfg, dataEEG_c100);

    %% Baseline (raw TFR)
    cfg              = [];
    cfg.baseline     = [-1.5 -.25];
    cfg.baselinetype = 'db';
    tfr_c25_bl       = ft_freqbaseline(cfg, tfr_c25);
    tfr_c50_bl       = ft_freqbaseline(cfg, tfr_c50);
    tfr_c75_bl       = ft_freqbaseline(cfg, tfr_c75);
    tfr_c100_bl      = ft_freqbaseline(cfg, tfr_c100);

    disp(upper('Raw TFR + baseline done...'))

    %% Sliding-window FOOOF on trial-averaged spectra (mtmfft)
    startWin_FOOOF = [-1.75 -1.25];   % 500 ms window
    steps_FOOOF    = 0.05;            % 50 ms step
    toi_FOOOF      = -1.75:0.05:2;    % desired centre times (matching your TFR toi)
    winLen         = 0.5;
    halfWin        = winLen/2;
    step           = 0.05;

    tfr_fooof = cell(1,4);
    for cond = 1:4

        if cond == 1
            dat      = dataEEG_c25;
            trlIdx   = ind61;
            condName = 'c25';
        elseif cond == 2
            dat      = dataEEG_c50;
            trlIdx   = ind62;
            condName = 'c50';
        elseif cond == 3
            dat      = dataEEG_c75;
            trlIdx   = ind63;
            condName = 'c75';
        else
            dat      = dataEEG_c100;
            trlIdx   = ind64;
            condName = 'c100';
        end

        tmin           = dat.time{1}(1);
        tmax           = dat.time{1}(end);
        toi_centres    = (tmin + halfWin) : step : (tmax - halfWin);
        nTimePnts      = numel(toi_centres);

        disp(' ')
        disp(['Running sliding-window FOOOF for condition ', condName])

        % FOOOF config (trial-averaged spectrum per window)
        cfg_fooof            = [];
        cfg_fooof.method     = 'mtmfft';
        cfg_fooof.taper      = 'hanning';
        cfg_fooof.foilim     = [1 100];
        cfg_fooof.pad        = 5;
        cfg_fooof.output     = 'fooof';
        cfg_fooof.keeptrials = 'no';

        % One test window for sizing / master freq grid
        cfg_sel0         = [];
        cfg_sel0.latency = startWin_FOOOF;
        cfg_sel0.trials  = trlIdx;
        dat_win0         = ft_selectdata(cfg_sel0, dat);

        fooof_test       = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat_win0);

        nChan            = numel(fooof_test.label);
        freq_master      = fooof_test.freq(:);
        nFreq            = numel(freq_master);

        fooof_powspctrm  = nan(nChan, nFreq, nTimePnts);   % model fit (FOOOF/log space)
        fooof_powspec    = nan(nChan, nFreq, nTimePnts);   % input power spectrum (FOOOF/log space)
        fooof_aperiodic  = nan(nChan, 4,    nTimePnts);    % [offset exponent error r^2]

        for timePnt = 1:nTimePnts

            centre = toi_centres(timePnt);
            cfg_sel = [];
            cfg_sel.latency = [centre - halfWin, centre + halfWin];
            cfg_sel.trials  = trlIdx;
            dat_win = ft_selectdata(cfg_sel, dat);

            % Tracker
            clc
            disp(['Running FOOOF for Subject ', num2str(subjects{subj})])
            disp(['Subject:    ', num2str(subj), ' / ', num2str(length(subjects))])
            disp(['Condition:  ', num2str(cond), ' / 4'])
            disp(['Time Point: ', num2str(timePnt), ' / ', num2str(nTimePnts)])

            % fooof_out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat_win);
            out = evalc('fooof_out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat_win);');

            if iscell(fooof_out.fooofparams)
                repdata = fooof_out.fooofparams{1};
            else
                repdata = fooof_out.fooofparams;
            end

            freq_now = fooof_out.freq(:);

            % Map freq_now onto the master grid
            [tf, loc] = ismembertol(freq_now, freq_master, 1e-10);

            local_model  = nan(nChan, nFreq);
            local_offset = nan(nChan, 1);
            local_expo   = nan(nChan, 1);
            local_err    = nan(nChan, 1);
            local_rsq    = nan(nChan, 1);

            local_peaks = nan(nChan, nFreq);   % peaks-only (aperiodic removed), on master grid
            local_ps    = nan(nChan, nFreq);

            for ch = 1:nChan

                % input spectrum (FOOOF/log space)
                ps_tmp = repdata(ch).power_spectrum(:);

                % build aperiodic fit on fooof_out.freq
                ap = repdata(ch).aperiodic_params(:);
                local_err(ch) = repdata(ch).error;
                local_rsq(ch) = repdata(ch).r_squared;

                if numel(ap) == 2
                    offset = ap(1);
                    expo   = ap(2);
                    ap_fit = offset - expo .* log10(freq_now);
                else
                    offset = ap(1);
                    knee   = ap(2);
                    expo   = ap(3);
                    ap_fit = offset - log10(knee + freq_now.^expo);
                end

                % enforce equal length between ps_tmp and ap_fit (THIS FIXES YOUR CRASH)
                nMin = min(numel(ps_tmp), numel(ap_fit));
                ps_use   = ps_tmp(1:nMin);
                ap_use   = ap_fit(1:nMin);
                freq_use = freq_now(1:nMin);

                % remap indices using the truncated grid (not freq_now)
                [tf_use, loc_use] = ismembertol(freq_use, freq_master, 1e-10);

                % store input spectrum on master grid
                local_ps(ch, loc_use(tf_use)) = ps_use(tf_use).';

                % aperiodic-corrected input spectrum on the same grid
                ps_corr = ps_use - ap_use;
                local_model(ch, loc_use(tf_use)) = ps_corr(tf_use).';

                local_offset(ch) = offset;
                local_expo(ch)   = expo;

                % peaks (Gaussian sum) on freq_use (NOT freq_now)
                pk = repdata(ch).peak_params;
                gauss_sum = zeros(numel(freq_use), 1);
                if ~isempty(pk)
                    for p = 1:size(pk,1)
                        cf  = pk(p,1);
                        amp = pk(p,2);
                        bw  = pk(p,3);
                        gauss_sum = gauss_sum + amp .* exp(-(freq_use - cf).^2 ./ (2*bw.^2));
                    end
                end

                local_peaks(ch, loc_use(tf_use)) = gauss_sum(tf_use).';
            end

            % store outputs
            fooof_powspctrm(:, :, timePnt) = local_model; % original powerspectrum - aperiodic fit
            fooof_powspec(:,   :, timePnt) = local_ps;
            fooof_aperiodic(:, 1, timePnt) = local_offset;
            fooof_aperiodic(:, 2, timePnt) = local_expo;
            fooof_aperiodic(:, 3, timePnt) = local_err;
            fooof_aperiodic(:, 4, timePnt) = local_rsq;

            s           = struct();
            s.subj      = subj;
            s.cond      = cond;
            s.time      = timePnt;
            s.nTimePnts = nTimePnts;
            %send(D, s);
        end

        % Build FieldTrip-like freq struct
        tfr_ff                = [];
        tfr_ff.label          = fooof_test.label;
        tfr_ff.freq           = fooof_test.freq;
        tfr_ff.time           = toi_centres(:).';
        tfr_ff.powspctrm      = fooof_powspctrm;
        tfr_ff.power_spectrum = fooof_powspec;
        tfr_ff.fooofparams    = fooof_aperiodic;
        tfr_ff.dimord         = 'chan_freq_time';

        tfr_fooof{cond}       = tfr_ff;
    end

    tfr_c25_fooof  = tfr_fooof{1};
    tfr_c50_fooof  = tfr_fooof{2};
    tfr_c75_fooof  = tfr_fooof{3};
    tfr_c100_fooof = tfr_fooof{4};

    disp(upper('FOOOF done on trial-averaged spectra (sliding windows)...'))

    %% Sanity Check
    %
    % 3-panel figure (1 row):
    %   left   = original input spectrum (FOOOF/log space)
    %   middle = input spectrum minus aperiodic fit (FOOOF/log space)
    %   right  = peaks-only model (aperiodic-removed fit; Gaussian sum)
    %
    % Notes:
    %   - Uses one condition only (set cond = 1..4).
    %   - Maps fooof_sc.freq onto your master freq grid for clean plotting.

    latWin  = [0.25 0.75];
    tCentre = mean(latWin);

    cond        = 2;  % 1=c25, 2=c50, 3=c75, 4=c100
    cond_titles = {'25% Contrast','50% Contrast','75% Contrast','100% Contrast'};

    dat_all = {dataEEG_c25, dataEEG_c50, dataEEG_c75, dataEEG_c100};
    trl_all = {ind61,      ind62,      ind63,      ind64};

    tfr_all     = {tfr_c25_fooof, tfr_c50_fooof, tfr_c75_fooof, tfr_c100_fooof};
    freq_master = tfr_all{cond}.freq(:);
    nFreqM      = numel(freq_master);

    dat    = dat_all{cond};
    trlIdx = trl_all{cond};

    % select the sanity-check window
    cfg_sel         = [];
    cfg_sel.latency = latWin;
    cfg_sel.trials  = trlIdx;
    dat_win_sc      = ft_selectdata(cfg_sel, dat);

    % run FOOOF
    fooof_sc = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat_win_sc);

    if iscell(fooof_sc.fooofparams)
        repdata_sc = fooof_sc.fooofparams{1};
    else
        repdata_sc = fooof_sc.fooofparams;
    end

    freq_sc = fooof_sc.freq(:);
    nChan   = numel(repdata_sc);

    % map sanity-check outputs onto master grid
    [tf, loc] = ismembertol(freq_sc, freq_master, 1e-10);

    ps_all    = nan(nFreqM, nChan);   % input spectrum (FOOOF/log space)
    ap_all    = nan(nFreqM, nChan);   % aperiodic fit (FOOOF/log space)
    gauss_all = nan(nFreqM, nChan);   % peaks-only fit (FOOOF/log space)

    for ch = 1:nChan

        ps_tmp = repdata_sc(ch).power_spectrum(:);

        ap = repdata_sc(ch).aperiodic_params(:);
        if numel(ap) == 2
            offset = ap(1);
            expo   = ap(2);
            ap_fit = offset - expo .* log10(freq_sc);
        else
            offset = ap(1);
            knee   = ap(2);
            expo   = ap(3);
            ap_fit = offset - log10(knee + freq_sc.^expo);
        end

        pk = repdata_sc(ch).peak_params;
        gauss_tmp = zeros(numel(freq_sc), 1);
        if ~isempty(pk)
            for p = 1:size(pk,1)
                cf  = pk(p,1);
                amp = pk(p,2);
                bw  = pk(p,3);
                gauss_tmp = gauss_tmp + amp .* exp(-(freq_sc - cf).^2 ./ (2*bw.^2));
            end
        end

        ps_all(loc(tf), ch)    = ps_tmp(tf);
        ap_all(loc(tf), ch)    = ap_fit(tf);
        gauss_all(loc(tf), ch) = gauss_tmp(tf);
    end

    % average across channels
    ps_in     = mean(ps_all,    2, 'omitnan');
    ap_fit    = mean(ap_all,    2, 'omitnan');
    peaks_fit = mean(gauss_all, 2, 'omitnan');

    ps_corr = ps_in - ap_fit;   % original minus aperiodic

    % plot: 1x3 layout
    figure('Position', [0 0 1500 420], 'Color', 'w');

    subplot(1,3,1); hold on
    plot(freq_master, ps_in, 'LineWidth', 2)
    title(sprintf('%s | Original Powspctrm', cond_titles{cond}))
    xlabel('Frequency (Hz)')
    ylabel('Power (FOOOF/log space)')
    set(gca, 'FontSize', 13)
    xlim([min(freq_master) max(freq_master)])

    subplot(1,3,2); hold on
    plot(freq_master, ps_corr, 'LineWidth', 2)
    title('Original Powspctrm - Aperiodic')
    xlabel('Frequency (Hz)')
    ylabel('Power (FOOOF/log space)')
    set(gca, 'FontSize', 13)
    xlim([min(freq_master) max(freq_master)])

    subplot(1,3,3); hold on
    plot(freq_master, peaks_fit, 'LineWidth', 2)
    title('Aperiodic-removed fit (peaks-only)')
    xlabel('Frequency (Hz)')
    ylabel('Power (FOOOF/log space)')
    set(gca, 'FontSize', 13)
    xlim([min(freq_master) max(freq_master)])

    sgtitle(sprintf('FOOOF sanity check: Subject %s | Window [%.2f %.2f] s | t = %.2f s', ...
        subjects{subj}, latWin(1), latWin(2), tCentre), 'FontSize', 18)

    if ispc
        savePathControls = 'W:\Students\Arne\GCP\data\controls\FOOOF\';
    else
        savePathControls = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/controls/FOOOF/';
    end
    if ~exist(savePathControls, 'dir')
        mkdir(savePathControls)
    end

    saveName = sprintf('GCP_controls_FOOOF_powspctrm_stern_subj%s.png', subjects{subj});
    saveas(gcf, fullfile(savePathControls, saveName));

    %% Baseline (FOOOFed TFR)
    cfg              = [];
    cfg.baseline     = [-1.5 -.25];
    cfg.baselinetype = 'absolute';
    tfr_c25_fooof_bl  = ft_freqbaseline(cfg, tfr_c25_fooof);
    tfr_c50_fooof_bl  = ft_freqbaseline(cfg, tfr_c50_fooof);
    tfr_c75_fooof_bl  = ft_freqbaseline(cfg, tfr_c75_fooof);
    tfr_c100_fooof_bl = ft_freqbaseline(cfg, tfr_c100_fooof);

    disp(upper('FOOOF baseline done...'))

    %% Save
    cd(datapath)

    save data_tfr_fooof_sw ...
        tfr_c25 tfr_c50 tfr_c75 tfr_c100 ...
        tfr_c25_bl tfr_c50_bl tfr_c75_bl tfr_c100_bl ...
        tfr_c25_fooof tfr_c50_fooof tfr_c75_fooof tfr_c100_fooof ...
        tfr_c25_fooof_bl tfr_c50_fooof_bl tfr_c75_fooof_bl tfr_c100_fooof_bl

    clc
    fprintf('Subject GCP %s (%.3d/%.3d) DONE (sliding-window FOOOF) \n', ...
        num2str(subjects{subj}), subj, length(subjects))
end
