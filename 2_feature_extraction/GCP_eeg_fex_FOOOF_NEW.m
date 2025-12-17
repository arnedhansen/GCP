%% GCP EEG Feature Extraction
%
% Outputs per subject:
%   TFR (mtmconvol) raw + baseline (for visualisation)
%   TFR_FOOOF (sliding windows): powspctrm = log10(FOOOFed spectrum) per window centre
%   plus: power_spectrum + aperiodic params (offset/slope/error/r^2) per window

%% Setup
startup
[subjects, path, ~, ~] = setup('GCP');

%% DataQueue for parfor progress
D = parallel.pool.DataQueue;
afterEach(D, @printProgress);
function printProgress(s)
fprintf('Subj %d | Cond %d | Time %d/%d finished\n', ...
    s.subj, s.cond, s.time, s.nTimePnts);
end

%% Loop subjects
for subj = 1 : length(subjects)

    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    close all
    clc

    disp(['Processing GCP TFR + sliding-window FOOOF for subject ', num2str(subjects{subj})])

    load dataEEG
    load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify trial indices (per dataset)
    % Note: adjust if your trialinfo is not a vector or codes differ per file.
    ind61 = find(dataEEG_c25.trialinfo  == 61);
    ind62 = find(dataEEG_c50.trialinfo  == 62);
    ind63 = find(dataEEG_c75.trialinfo  == 63);
    ind64 = find(dataEEG_c100.trialinfo == 64);

    %% TFR (mtmconvol) for visualisation / later use
    cfg             = [];
    cfg.output      = 'pow';
    cfg.method      = 'mtmconvol';
    cfg.taper       = 'dpss';
    cfg.foi         = 30:5:120;
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
    nTimePnts      = numel(toi_FOOOF);

    toi_centres    = startWin_FOOOF(1):steps_FOOOF:(startWin_FOOOF(1) + steps_FOOOF*(nTimePnts-1));
    toi_centres    = toi_centres + 0.25;  % shift by half window length

    % Conditions container
    tfr_fooof = cell(1,4);

    for cond = 1:4

        if cond == 1
            dat    = dataEEG_c25;
            trlIdx = ind61;
            condName = 'c25';
        elseif cond == 2
            dat    = dataEEG_c50;
            trlIdx = ind62;
            condName = 'c50';
        elseif cond == 3
            dat    = dataEEG_c75;
            trlIdx = ind63;
            condName = 'c75';
        elseif cond == 4
            dat    = dataEEG_c100;
            trlIdx = ind64;
            condName = 'c100';
        end

        disp(' ')
        disp(['Running sliding-window FOOOF for condition ', condName])

        % FOOOF config (trial-averaged spectrum per window)
        cfg_fooof            = [];
        cfg_fooof.method     = 'mtmfft';
        cfg_fooof.taper      = 'hanning';
        cfg_fooof.foilim     = [30 120];
        cfg_fooof.pad        = 5;
        cfg_fooof.output     = 'fooof';
        cfg_fooof.keeptrials = 'no';

        % One test window for sizing
        cfg_sel0         = [];
        cfg_sel0.latency = startWin_FOOOF;
        cfg_sel0.trials  = trlIdx;
        dat_win0         = ft_selectdata(cfg_sel0, dat);

        fooof_test       = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat_win0);

        nChan            = numel(fooof_test.label);
        nFreq            = numel(fooof_test.freq);

        fooof_powspctrm  = nan(nChan, nFreq, nTimePnts);
        fooof_powspec    = nan(nChan, nFreq, nTimePnts);
        fooof_aperiodic  = nan(nChan, 4,    nTimePnts);   % [offset slope error r^2]

        % parfor over timepoints
        parfor timePnt = 1 : nTimePnts

            cfg_sel         = [];
            cfg_sel.latency = startWin_FOOOF + steps_FOOOF * (timePnt-1);
            cfg_sel.trials  = trlIdx;
            dat_win         = ft_selectdata(cfg_sel, dat);

            fooof_out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat_win);

            freq_master = fooof_test.freq(:);
            freq_now    = fooof_out.freq(:);

            % Index mapping: where does freq_now land on freq_master?
            [tf, loc] = ismembertol(freq_now, freq_master, 1e-10);

            % Preallocate local arrays on the MASTER grid
            local_pow = nan(nChan, numel(freq_master));
            local_ps  = nan(nChan, numel(freq_master));

            if iscell(fooof_out.fooofparams)
                repdata = fooof_out.fooofparams{1};
            else
                repdata = fooof_out.fooofparams;
            end

            tmpaper      = {repdata.aperiodic_params};
            tmperror     = {repdata.error};
            tmpr_sq      = {repdata.r_squared};
            tmp_pwr_spec = {repdata.power_spectrum};

            local_aper = nan(nChan, 2);
            local_err  = nan(nChan, 1);
            local_rsq  = nan(nChan, 1);

            for ch = 1:nChan

                ap               = tmpaper{ch};
                local_aper(ch,:) = ap(:).';

                local_err(ch,1)  = tmperror{ch};
                local_rsq(ch,1)  = tmpr_sq{ch};

                ps_tmp = tmp_pwr_spec{ch}(:);

                % Write only matching bins into the master grid
                if numel(ps_tmp) == numel(freq_now)
                    local_ps(ch, loc(tf))  = ps_tmp(tf).';
                else
                    % defensive: if wrapper returns mismatch internally
                    nMin = min(numel(ps_tmp), numel(freq_now));
                    local_ps(ch, loc(tf(1:nMin))) = ps_tmp(1:nMin).';
                end
            end

            % Same alignment for powspctrm (FOOOF-space) if its freq grid matches freq_now
            pow_tmp = fooof_out.powspctrm;
            if size(pow_tmp,2) == numel(freq_now)
                local_pow(:, loc(tf)) = log10(pow_tmp(:, tf));
            else
                % last-resort: just log10 the available columns and align what you can
                nMin = min(size(pow_tmp,2), numel(freq_now));
                local_pow(:, loc(tf(1:nMin))) = log10(pow_tmp(:, 1:nMin));
            end

            fooof_powspctrm(:, :, timePnt) = local_pow;
            fooof_powspec(:,   :, timePnt) = local_ps;

            local_aper_all       = nan(nChan, 4);
            local_aper_all(:,1)  = local_aper(:,1);
            local_aper_all(:,2)  = local_aper(:,2);
            local_aper_all(:,3)  = local_err;
            local_aper_all(:,4)  = local_rsq;

            fooof_aperiodic(:, :, timePnt) = local_aper_all;

            s           = struct();
            s.subj      = subj;
            s.cond      = cond;
            s.time      = timePnt;
            s.nTimePnts = nTimePnts;
            send(D, s);
        end

        % Build FieldTrip-like freq struct
        tfr_ff                    = [];
        tfr_ff.label              = fooof_test.label;
        tfr_ff.freq               = fooof_test.freq;
        tfr_ff.time               = toi_centres(:).';
        tfr_ff.powspctrm          = fooof_powspctrm;
        tfr_ff.power_spectrum     = fooof_powspec;
        tfr_ff.fooofparams        = fooof_aperiodic;
        tfr_ff.dimord             = 'chan_freq_time';

        tfr_fooof{cond}           = tfr_ff;
    end

    tfr_c25_fooof  = tfr_fooof{1};
    tfr_c50_fooof  = tfr_fooof{2};
    tfr_c75_fooof  = tfr_fooof{3};
    tfr_c100_fooof = tfr_fooof{4};

    disp(upper('FOOOF done on trial-averaged spectra (sliding windows)...'))

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
