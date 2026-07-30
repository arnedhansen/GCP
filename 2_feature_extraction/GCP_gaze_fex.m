%% GCP Gaze Feature Extraction
%
% Extracted features:
%   BCEA (Bivariate Contour Ellipse Area, k=2.291, 95%)
%   Pupil size (Time Series Raw, Baselined % change)
%   Microsaccades (TC ground truth: boxplot scalar = mean of % TC)
%   Eye velocity (Raw, Baselined % change)
%   EyeLink blinks, fixations, saccades (rates from preprocessing gaze_metrics; baselined % change)
%
% Subject x condition scalars are saved for full [0 2], early [0 1], and
% late [1 2] windows: features/GCP_gaze_window_summaries.mat

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('GCP');

% across-subjects raw struct (fields filled from subj_data_gaze)
gaze_data = struct([]);

% time‐windows
baseline_period   = [-1.5, -0.5];
analysis_period   = [0 2];          % full
analysis_early    = [0 1];
analysis_late     = [1 2];
analysis_periodTS = [-1 2];
pupil_store_window = [-1.5, 2.5];
vel_store_window   = [-1.5, 2.5];  % extend past 2 s so t=2 is not a kernel edge
vel_smooth_s       = 0.050;        % light smoothing of eye speed TC
ms_min_rate_hz     = 0.1;          % exclude trials below this in baseline or stimulus
win_size          = 25;    % blink‐removal window (samples)
fsample           = 500;   % eye‐tracker sampling rate

% prepare raw gaze storage across all subjects
gaze_x_c25   = {};  gaze_y_c25   = {};
gaze_x_c50   = {};  gaze_y_c50   = {};
gaze_x_c75   = {};  gaze_y_c75   = {};
gaze_x_c100  = {};  gaze_y_c100  = {};

%% Loop over subjects
clc
for subj = 1:numel(subjects)

    % Load preprocessed ET data
    clc
    fprintf('Loading Subject %d/%d...\n', subj, numel(subjects));
    datapath = fullfile(paths.features, subjects{subj}, 'gaze');
    cd(datapath)
    load(fullfile(datapath,'dataET'));

    %% Loop over conditions
    for conds = {'c25','c50','c75','c100'}
        cond = conds{1};

        % Pick dataET
        switch cond
            case 'c25', dataET = dataET_c25;
            case 'c50', dataET = dataET_c50;
            case 'c75', dataET = dataET_c75;
            case 'c100',dataET = dataET_c100;
        end

        % initialise per‐trial arrays
        subject_id         = [];
        trial_num          = [];
        condition          = [];

        bcea               = [];  baselineBcea      = [];
        bcea_early         = [];  bcea_late         = [];
        pupilSize          = [];  baselinePupilSize = [];
        microsaccadeRate   = [];  baselineMSRate    = [];
        microsaccadeRate_early = [];  microsaccadeRate_late = [];

        % Prepare velocity time series container for this condition
        velocityData     = dataET;    % copy meta-info
        velocityData.label = {'VelH','VelV','Vel2D'};

        % trial-level % baseline (scalar baseline per trial; field names keep *_db)
        velocityData_db = dataET;
        velocityData_db.label = {'VelH_bl','VelV_bl','Vel2D_bl'};

        fsample          = dataET.fsample;
        vel_kernel       = [1 1 0 -1 -1] * (fsample / 6); % Engbert velocity kernel
        vel_kernel_pad   = ceil(numel(vel_kernel) / 2);
        nTrials          = numel(dataET.trial);

        % For linking baseline and trial-wise means
        velHorz        = nan(1, nTrials);
        velVert        = nan(1, nTrials);
        vel2D          = nan(1, nTrials);
        baselineVelH   = nan(1, nTrials);
        baselineVelV   = nan(1, nTrials);
        baselineVel2D  = nan(1, nTrials);

        % FieldTrip-ready velocity container (full baseline + analysis window)
        velFT = [];
        velFT.label     = {'VelH','VelV','Vel2D'};
        velFT.fsample   = fsample;
        velFT.trial     = cell(1, nTrials);
        velFT.time      = cell(1, nTrials);
        velFT.trialinfo = dataET.trialinfo;

        % Continuous rectified velocity for OCC (same as velFT; no MS zeroing)
        velOCCFT = velFT;

        % FieldTrip-ready pupil container (full baseline + analysis window)
        pupFT = [];
        pupFT.label     = {'Pupil'};
        pupFT.fsample   = fsample;
        pupFT.trial     = cell(1, nTrials);
        pupFT.time      = cell(1, nTrials);
        pupFT.trialinfo = dataET.trialinfo;

        % FieldTrip-ready microsaccade rate container (full baseline + analysis window)
        msFT = [];
        msFT.label     = {'MSRate'};   % Hz
        msFT.fsample   = fsample;
        msFT.trial     = cell(1, nTrials);
        msFT.time      = cell(1, nTrials);
        msFT.trialinfo = dataET.trialinfo;

        % Kernel for microsaccade rate (event density -> Hz)
        ms_sigma_s   = 0.05;                             % 50 ms
        ms_sigma_smp = max(1, round(ms_sigma_s*fsample));
        ker_half     = 4 * ms_sigma_smp;                 % +/- 4 sigma
        ker_x        = -ker_half:ker_half;
        ms_kernel    = exp(-0.5 * (ker_x/ms_sigma_smp).^2);
        ms_kernel    = ms_kernel ./ sum(ms_kernel);      % area = 1

        % Trial-level microsaccade QC (applied before msFT is filled)
        qc_min_valid_frac      = 0.60;                  % min valid fraction in [-1.5, 2]
        qc_max_blink_loss_frac = 0.40;                  % max blink-loss fraction after remove_blinks
        qc_min_valid_samples   = round(0.50 * fsample); % at least 500 ms valid data
        qc_max_ms_rate_hz      = 8;                     % reject implausibly high trial-level MS rates

        % Engbert MS event times (seconds on trial time axis)
        ms_onsets_per_trial  = cell(1, nTrials);
        ms_offsets_per_trial = cell(1, nTrials);

        %% Trial loop
        for trl = 1:numel(dataET.trialinfo)
            raw    = dataET.trial{trl};
            tVec   = dataET.time{trl};

            % Eye velocity on the original regular time axis (no MS sample removal).
            % Store past analysis end so samples near t=2 are not convolution edges.
            full_idx = tVec >= vel_store_window(1) & tVec <= vel_store_window(2);
            t_full   = tVec(full_idx);
            full_dat = raw(1:3, full_idx);

            % flip y axis to screen coordinates
            full_dat(2,:) = 600 - full_dat(2,:);

            % Mark invalid samples and blink contaminated intervals as missing.
            valid_full_position = full_dat(1,:) >= 0 & full_dat(1,:) <= 800 & ...
                                  full_dat(2,:) >= 0 & full_dat(2,:) <= 600;
            full_dat(:, ~valid_full_position) = NaN;
            full_dat = remove_blinks(full_dat, win_size);
            valid_full_position = isfinite(full_dat(1,:)) & isfinite(full_dat(2,:));

            if numel(t_full) > 1 && any(valid_full_position)
                xy_padded = ft_preproc_padding(full_dat(1:2,:), 'localmean', vel_kernel_pad);
                vel_signed = convn(xy_padded, vel_kernel, 'same');
                vel_signed = ft_preproc_padding(vel_signed, 'remove', vel_kernel_pad);
                vel_signed(:, ~valid_full_position) = NaN;

                vx_full    = abs(vel_signed(1,:));
                vy_full    = abs(vel_signed(2,:));
                speed_full = sqrt(sum(vel_signed.^2, 1));

                % Light temporal smoothing after differentiation (not before).
                vel_smooth_samp = max(1, round(vel_smooth_s * fsample));
                if vel_smooth_samp > 1
                    vx_full = movmean(vx_full, vel_smooth_samp, 'omitnan');
                    vy_full = movmean(vy_full, vel_smooth_samp, 'omitnan');
                    speed_full = movmean(speed_full, vel_smooth_samp, 'omitnan');
                end

                velFT.trial{trl} = [vx_full; vy_full; speed_full];
                velFT.time{trl}  = t_full;

                % OCC uses the same continuous rectified velocity.
                velOCCFT.trial{trl} = velFT.trial{trl};
                velOCCFT.time{trl} = t_full;
            else
                velFT.trial{trl} = nan(3, numel(t_full));
                velFT.time{trl}  = t_full;
                velOCCFT.trial{trl} = nan(3, numel(t_full));
                velOCCFT.time{trl}  = t_full;
            end

            % BASELINE WINDOW (position-based metrics)
            bl_idx = tVec >= baseline_period(1) & tVec <= baseline_period(2);
            bl_dat = raw(:,bl_idx);
            valid  = bl_dat(1,:)>=0 & bl_dat(1,:)<=800 & bl_dat(2,:)>=0 & bl_dat(2,:)<=600;
            bl_dat = bl_dat(1:3, valid);
            bl_dat(2,:) = 600 - bl_dat(2,:);
            bl_dat = remove_blinks(bl_dat, win_size);

            bl_x = bl_dat(1,:);  bl_y = bl_dat(2,:);
            baseline_std_x    = nanstd(bl_x);
            baseline_std_y    = nanstd(bl_y);
            if numel(bl_x) > 2
                rho_bl = corr(bl_x(:), bl_y(:));
            else
                rho_bl = 0;
            end
            baseline_bcea_val = 2 * 2.291 * pi * baseline_std_x * baseline_std_y * sqrt(1 - rho_bl^2);
            baseline_pupil    = mean(bl_dat(3,:),'omitnan')/1000;
            [baseline_msrate, ~] = detect_microsaccades(fsample, [bl_x; bl_y], numel(bl_x));

            % Baseline eye velocity
            if ~isempty(velFT.time{trl})
                vel_bl_idx = velFT.time{trl} >= baseline_period(1) & ...
                             velFT.time{trl} <= baseline_period(2);
                baselineVelH(trl)  = mean(velFT.trial{trl}(1,vel_bl_idx), 'omitnan');
                baselineVelV(trl)  = mean(velFT.trial{trl}(2,vel_bl_idx), 'omitnan');
                baselineVel2D(trl) = mean(velFT.trial{trl}(3,vel_bl_idx), 'omitnan');
            else
                baselineVelH(trl)  = NaN;
                baselineVelV(trl)  = NaN;
                baselineVel2D(trl) = NaN;
            end

            % ANALYSIS WINDOW [0.3, 2.0] (cleaned position, as before)
            an_idx = tVec >= analysis_period(1) & tVec <= analysis_period(2);
            an_dat = raw(:,an_idx);
            valid  = an_dat(1,:)>=0 & an_dat(1,:)<=800 & an_dat(2,:)>=0 & an_dat(2,:)<=600;
            an_dat = an_dat(1:3, valid);
            an_dat(2,:) = 600 - an_dat(2,:);
            an_dat = remove_blinks(an_dat, win_size);

            % extract the cleaned gaze trace
            x = an_dat(1,:);
            y = an_dat(2,:);

            % Store raw gaze traces for this subject/condition/trial
            switch cond
                case 'c25'
                    gaze_x_c25{subj,trl}  = x;
                    gaze_y_c25{subj,trl}  = y;
                case 'c50'
                    gaze_x_c50{subj,trl}  = x;
                    gaze_y_c50{subj,trl}  = y;
                case 'c75'
                    gaze_x_c75{subj,trl}  = x;
                    gaze_y_c75{subj,trl}  = y;
                case 'c100'
                    gaze_x_c100{subj,trl} = x;
                    gaze_y_c100{subj,trl} = y;
            end

            % BCEA, microsaccades
            std_x       = nanstd(x);
            std_y       = nanstd(y);
            if numel(x) > 2
                rho_an = corr(x(:), y(:));
            else
                rho_an = 0;
            end
            bcea_val    = 2 * 2.291 * pi * std_x * std_y * sqrt(1 - rho_an^2);
            bcea_early_val = bcea_in_window(raw, tVec, analysis_early, win_size);
            bcea_late_val  = bcea_in_window(raw, tVec, analysis_late, win_size);
            pupil       = mean(an_dat(3,:),'omitnan')/1000;
            [msrate, ~] = detect_microsaccades(fsample, [x; y], numel(x));
            msrate_early = ms_rate_in_window(raw, tVec, analysis_early, win_size, fsample);
            msrate_late  = ms_rate_in_window(raw, tVec, analysis_late, win_size, fsample);

            % Analysis window eye velocity from the full regular trace
            if ~isempty(velFT.time{trl})
                vel_an_idx = velFT.time{trl} >= analysis_period(1) & ...
                             velFT.time{trl} <= analysis_period(2);
                velocityData.trial{trl} = velFT.trial{trl}(:, vel_an_idx);
                velocityData.time{trl}  = velFT.time{trl}(vel_an_idx);

                % Trial-wise mean velocities (for scalar features)
                velHorz(trl) = mean(velocityData.trial{trl}(1,:), 'omitnan');
                velVert(trl) = mean(velocityData.trial{trl}(2,:), 'omitnan');
                vel2D(trl)   = mean(velocityData.trial{trl}(3,:), 'omitnan');

                % Baseline-normalised velocity time series (% change, scalar baseline)
                db_vx = compute_pct_baseline( ...
                    velocityData.trial{trl}(1,:), baselineVelH(trl));
                db_vy = compute_pct_baseline( ...
                    velocityData.trial{trl}(2,:), baselineVelV(trl));
                db_speed = compute_pct_baseline( ...
                    velocityData.trial{trl}(3,:), baselineVel2D(trl));

                velocityData_db.trial{trl} = zeros(3, numel(db_speed));
                velocityData_db.trial{trl}(1,:) = db_vx;
                velocityData_db.trial{trl}(2,:) = db_vy;
                velocityData_db.trial{trl}(3,:) = db_speed;
                velocityData_db.time{trl}       = velocityData.time{trl};

            else
                % Too few points → fill with NaNs
                velocityData.trial{trl}      = nan(3,0);
                velocityData.time{trl}       = [];
                velocityData_db.trial{trl}  = nan(3,0);
                velocityData_db.time{trl}   = [];
            end

            % Ensure FieldTrip format for pupil data
            % Wider pupil storage window for edge-safe smoothing later
            pup_idx = tVec >= pupil_store_window(1) & tVec <= pupil_store_window(2);
            t_full = tVec(pup_idx);
            p_full = raw(3, pup_idx);    % 1 x N already
            p_full = p_full ./ 1000;      % keep units consistent with scalar pupil features
            pupFT.trial{trl} = p_full;
            pupFT.time{trl}  = t_full;

            % Microsaccade rate time series in full window [-1.5, 2]
            t_full_ms = tVec(full_idx);

            x_full = raw(1, full_idx);
            y_full = raw(2, full_idx);

            % Apply same screen-coordinate flip as elsewhere
            y_full = 600 - y_full;

            % Validity mask (match your position pipeline)
            valid_full = x_full>=0 & x_full<=800 & y_full>=0 & y_full<=600;

            x_val = x_full(valid_full);
            y_val = y_full(valid_full);

            % Blink removal (same function you already rely on)
            full_dat = [x_val; y_val; nan(1, numel(x_val))];
            full_dat = remove_blinks(full_dat, win_size);
            x_val = full_dat(1,:);
            y_val = full_dat(2,:);

            % Trial-level QC before constructing msFT traces
            n_total_full = numel(t_full_ms);
            n_valid_before = sum(isfinite(x_full(valid_full)) & isfinite(y_full(valid_full)));
            valid_clean = isfinite(x_val) & isfinite(y_val);
            n_valid_after = sum(valid_clean);

            valid_frac_full = n_valid_after / max(n_total_full, 1);
            blink_loss_frac = max(0, (n_valid_before - n_valid_after) / max(n_valid_before, 1));

            if n_valid_after < qc_min_valid_samples || ...
               valid_frac_full < qc_min_valid_frac || ...
               blink_loss_frac > qc_max_blink_loss_frac
                msFT.trial{trl} = nan(1, n_total_full);
                msFT.time{trl}  = t_full_ms;
                ms_onsets_per_trial{trl} = [];
                ms_offsets_per_trial{trl} = [];
            else
                x_clean = x_val(valid_clean);
                y_clean = y_val(valid_clean);

                % Detect microsaccades in cleaned valid samples only
                [~, ms_det] = detect_microsaccades(fsample, [x_clean; y_clean], numel(x_clean));

                % Map Engbert sample indices (in x_clean) back to trial time (s)
                idx_full_valid = find(valid_full);
                idx_clean_in_full = idx_full_valid(valid_clean);
                if ~isempty(ms_det.Onset) && ~isempty(ms_det.Offset)
                    nEv = min(numel(ms_det.Onset), numel(ms_det.Offset));
                    on_idx = ms_det.Onset(1:nEv);
                    off_idx = ms_det.Offset(1:nEv);
                    keep_ev = on_idx >= 1 & on_idx <= numel(idx_clean_in_full) & ...
                              off_idx >= 1 & off_idx <= numel(idx_clean_in_full);
                    on_idx = on_idx(keep_ev);
                    off_idx = off_idx(keep_ev);
                    ms_onsets_per_trial{trl} = t_full_ms(idx_clean_in_full(on_idx))';
                    ms_offsets_per_trial{trl} = t_full_ms(idx_clean_in_full(off_idx))';
                else
                    ms_onsets_per_trial{trl} = [];
                    ms_offsets_per_trial{trl} = [];
                end

                % Impulse train on valid samples
                ms_imp = zeros(1, numel(x_clean));
                onsets = [];
                if ~isempty(ms_det.Onset)
                    onsets = ms_det.Onset(:)';
                    onsets = onsets(onsets >= 1 & onsets <= numel(ms_imp));
                    ms_imp(onsets) = 1;
                end

                ms_rate_trial = numel(onsets) / (numel(x_clean) / fsample);
                if ~isfinite(ms_rate_trial) || ms_rate_trial > qc_max_ms_rate_hz
                    msFT.trial{trl} = nan(1, n_total_full);
                    msFT.time{trl}  = t_full_ms;
                    % Keep detected event times even if TC rate QC fails
                else
                    % Smooth impulses -> event density per sample.
                    % Edge-effect correction: renormalise by local kernel mass so
                    % boundary samples are not artificially attenuated.
                    ms_kernel_mass = conv(ones(size(ms_imp)), ms_kernel, 'same');
                    ms_kernel_mass(ms_kernel_mass < eps) = NaN;
                    ms_rate_clean = conv(ms_imp, ms_kernel, 'same') ./ ms_kernel_mass * fsample;

                    % Expand back to full in-bounds vector (NaN at blink/invalid samples)
                    ms_rate_val = nan(1, numel(x_val));
                    ms_rate_val(valid_clean) = ms_rate_clean;

                    % Map back to full continuous time axis (NaN for invalid samples)
                    ms_rate_full = nan(1, numel(t_full_ms));
                    ms_rate_full(valid_full) = ms_rate_val;

                    % Store in FieldTrip struct with the original regular time axis
                    % (TC is for plotting only; boxplot scalars use direct counts)
                    msFT.trial{trl} = ms_rate_full;
                    msFT.time{trl}  = t_full_ms;

                    % Drop sparse trials from the TC average only
                    bl_ms_idx = t_full_ms >= baseline_period(1) & t_full_ms <= baseline_period(2);
                    an_ms_idx = t_full_ms >= analysis_period(1) & t_full_ms <= analysis_period(2);
                    bl_tc = mean(ms_rate_full(bl_ms_idx), 'omitnan');
                    an_tc = mean(ms_rate_full(an_ms_idx), 'omitnan');
                    if ~isfinite(bl_tc) || ~isfinite(an_tc) || ...
                            bl_tc < ms_min_rate_hz || an_tc < ms_min_rate_hz
                        msFT.trial{trl} = nan(1, n_total_full);
                        msFT.time{trl}  = t_full_ms;
                    end
                end
            end

            % append to trial‐wise arrays
            subject_id(end+1)       = str2double(subjects{subj});
            trial_num(end+1)        = trl;
            condition(end+1)        = dataET.trialinfo(trl) - 60;
            bcea(end+1)             = bcea_val;
            bcea_early(end+1)       = bcea_early_val;
            bcea_late(end+1)        = bcea_late_val;
            baselineBcea(end+1)     = baseline_bcea_val;
            pupilSize(end+1)        = pupil;
            baselinePupilSize(end+1)= baseline_pupil;
            microsaccadeRate(end+1) = msrate;
            microsaccadeRate_early(end+1) = msrate_early;
            microsaccadeRate_late(end+1)  = msrate_late;
            baselineMSRate(end+1)   = baseline_msrate;
        end

        % FieldTrip timelock and baseline on velFT
        % Force identical time axes: empty/unequal trials break ft_selectdata.
        t_vel_common = vel_store_window(1):(1 / fsample):vel_store_window(2);
        velFT = force_common_timeaxis(velFT, t_vel_common);
        velOCCFT = force_common_timeaxis(velOCCFT, t_vel_common);

        % Timelocked average without baseline
        cfg = [];
        cfg.latency    = analysis_periodTS;   % [-1 2]
        cfg.keeptrials = 'no';
        velTS_noBL = ft_timelockanalysis(cfg, velFT);

        % Average Hz first, then % change (avoids noisy trial-wise ratios)
        vel_bl_idx = velTS_noBL.time >= baseline_period(1) & ...
                     velTS_noBL.time <= baseline_period(2);
        velTS_BL_db = velTS_noBL;
        for ch = 1:size(velTS_noBL.avg, 1)
            bl_mu = mean(velTS_noBL.avg(ch, vel_bl_idx), 'omitnan');
            velTS_BL_db.avg(ch, :) = compute_pct_baseline(velTS_noBL.avg(ch, :), bl_mu);
        end

        % Also keep trial-level % traces for OCC / diagnostics
        velFT_bl_db = velFT;
        for trl = 1:nTrials
            if all(~isfinite(velFT.trial{trl}(:)))
                continue
            end
            velFT_bl_db.trial{trl}(1,:) = compute_pct_baseline(velFT.trial{trl}(1,:), baselineVelH(trl));
            velFT_bl_db.trial{trl}(2,:) = compute_pct_baseline(velFT.trial{trl}(2,:), baselineVelV(trl));
            velFT_bl_db.trial{trl}(3,:) = compute_pct_baseline(velFT.trial{trl}(3,:), baselineVel2D(trl));
        end
        velFT_bl_db = force_common_timeaxis(velFT_bl_db, t_vel_common);

        % Timelocked average without baseline
        t_pup_common = pupil_store_window(1):(1 / fsample):pupil_store_window(2);
        pupFT = force_common_timeaxis(pupFT, t_pup_common);

        cfg = [];
        cfg.latency    = pupil_store_window;   % edge-safe pupil window
        cfg.keeptrials = 'no';
        pupTS_noBL = ft_timelockanalysis(cfg, pupFT);

        % % change: scalar trial baseline on pupil traces
        pupFT_bl_db = pupFT;
        for trl = 1:nTrials
            if all(~isfinite(pupFT.trial{trl}(:)))
                continue
            end
            pupFT_bl_db.trial{trl} = compute_pct_baseline(pupFT.trial{trl}, baselinePupilSize(trl));
        end
        pupFT_bl_db = force_common_timeaxis(pupFT_bl_db, t_pup_common);

        cfg = [];
        cfg.latency    = pupil_store_window;
        cfg.keeptrials = 'no';
        pupTS_BL_db = ft_timelockanalysis(cfg, pupFT_bl_db);

        % Timelocked average without baseline
        t_ms_common = baseline_period(1):(1 / fsample):analysis_period(2);
        msFT = force_common_timeaxis(msFT, t_ms_common);

        cfg = [];
        cfg.latency    = analysis_periodTS;
        cfg.keeptrials = 'no';
        msTS_noBL = ft_timelockanalysis(cfg, msFT);

        % Per-trial % TC for plotting (sparse trials already NaN in msFT)
        msFT_bl_db = msFT;
        for trl = 1:nTrials
            if all(~isfinite(msFT.trial{trl}(:))) || ...
                    ~isfinite(baselineMSRate(trl)) || baselineMSRate(trl) <= 0
                continue
            end
            msFT_bl_db.trial{trl} = compute_pct_baseline(msFT.trial{trl}, baselineMSRate(trl));
        end
        msFT_bl_db = force_common_timeaxis(msFT_bl_db, t_ms_common);

        cfg = [];
        cfg.latency    = analysis_periodTS;
        cfg.keeptrials = 'no';
        msTS_BL_db = ft_timelockanalysis(cfg, msFT_bl_db);

        % Trial-level baselined scalars (% change). MS uses direct count rates.
        BCEA_bl          = compute_pct_baseline(bcea, baselineBcea);
        BCEA_bl_early    = compute_pct_baseline(bcea_early, baselineBcea);
        BCEA_bl_late     = compute_pct_baseline(bcea_late, baselineBcea);
        PupilSize_bl     = compute_pct_baseline(pupilSize, baselinePupilSize);
        VelH_bl          = compute_pct_baseline(velHorz, baselineVelH);
        VelV_bl          = compute_pct_baseline(velVert, baselineVelV);
        Vel2D_bl         = compute_pct_baseline(vel2D, baselineVel2D);
        MSRate_bl        = compute_pct_baseline(microsaccadeRate, baselineMSRate);
        MSRate_bl_early  = compute_pct_baseline(microsaccadeRate_early, baselineMSRate);
        MSRate_bl_late   = compute_pct_baseline(microsaccadeRate_late, baselineMSRate);

        % Subject boxplot scalars: TC window means for vel/pupil; trial means for MS/BCEA
        [vel_tc_full_H, vel_tc_early_H, vel_tc_late_H] = tc_window_means( ...
            velTS_BL_db.avg(1, :), velTS_BL_db.time, analysis_period, analysis_early, analysis_late);
        [vel_tc_full_V, vel_tc_early_V, vel_tc_late_V] = tc_window_means( ...
            velTS_BL_db.avg(2, :), velTS_BL_db.time, analysis_period, analysis_early, analysis_late);
        [vel_tc_full_2D, vel_tc_early_2D, vel_tc_late_2D] = tc_window_means( ...
            velTS_BL_db.avg(3, :), velTS_BL_db.time, analysis_period, analysis_early, analysis_late);
        [pup_tc_full, pup_tc_early, pup_tc_late] = tc_window_means( ...
            pupTS_BL_db.avg(1, :), pupTS_BL_db.time, analysis_period, analysis_early, analysis_late);

        %% SUBJECT‐BY‐CONDITION AVERAGES
        switch cond
            case 'c25'
                c25_bcea      = mean(bcea,'omitnan');
                c25_bl_bcea   = mean(baselineBcea,'omitnan');
                c25_pups      = mean(pupilSize,'omitnan');
                c25_bl_pups   = mean(baselinePupilSize,'omitnan');
                c25_msrate    = mean(microsaccadeRate,'omitnan');
                c25_bl_msrate = mean(baselineMSRate,'omitnan');

                % Velocity: raw and baseline means
                c25_velHorz      = mean(velHorz,'omitnan');
                c25_bl_velHorz   = mean(baselineVelH,'omitnan');
                c25_velVert      = mean(velVert,'omitnan');
                c25_bl_velVert   = mean(baselineVelV,'omitnan');
                c25_vel2D        = mean(vel2D,'omitnan');
                c25_bl_vel2D     = mean(baselineVel2D,'omitnan');

                % Condition scalars: trial means for MS/BCEA; TC means for vel/pupil
                c25_bcea_bl         = mean(BCEA_bl, 'omitnan');
                c25_bcea_bl_early   = mean(BCEA_bl_early, 'omitnan');
                c25_bcea_bl_late    = mean(BCEA_bl_late, 'omitnan');
                c25_pups_bl         = pup_tc_full;
                c25_pups_bl_early   = pup_tc_early;
                c25_pups_bl_late    = pup_tc_late;
                c25_msrate_bl       = mean(MSRate_bl, 'omitnan');
                c25_msrate_bl_early = mean(MSRate_bl_early, 'omitnan');
                c25_msrate_bl_late  = mean(MSRate_bl_late, 'omitnan');
                c25_velHorz_bl       = vel_tc_full_H;
                c25_velHorz_bl_early = vel_tc_early_H;
                c25_velHorz_bl_late  = vel_tc_late_H;
                c25_velVert_bl       = vel_tc_full_V;
                c25_velVert_bl_early = vel_tc_early_V;
                c25_velVert_bl_late  = vel_tc_late_V;
                c25_vel2D_bl         = vel_tc_full_2D;
                c25_vel2D_bl_early   = vel_tc_early_2D;
                c25_vel2D_bl_late    = vel_tc_late_2D;

                subj_data_gaze_trial_c25 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'BCEA',bcea, 'BCEA_early',bcea_early, 'BCEA_late',bcea_late, ...
                    'BaselineBCEA',baselineBcea, 'BCEA_bl', BCEA_bl, 'BCEA_bl_early', BCEA_bl_early, 'BCEA_bl_late', BCEA_bl_late, ...
                    'PupilSize',pupilSize, 'BaselinePupilSize',baselinePupilSize, 'PupilSize_bl', PupilSize_bl, ...
                    'MSRate',microsaccadeRate, 'MSRate_early',microsaccadeRate_early, 'MSRate_late',microsaccadeRate_late, ...
                    'BaselineMSRate',baselineMSRate, 'MSRate_bl', MSRate_bl, 'MSRate_bl_early', MSRate_bl_early, 'MSRate_bl_late', MSRate_bl_late, ...
                    'VelH',velHorz, 'BaselineVelH',baselineVelH, 'VelH_bl', VelH_bl, ...
                    'VelV',velVert, 'BaselineVelV',baselineVelV, 'VelV_bl', VelV_bl, ...
                    'Vel2D',vel2D, 'BaselineVel2D',baselineVel2D, 'Vel2D_bl', Vel2D_bl );

                % Store FieldTrip velocity for this condition
                velTS_c25        = velTS_noBL;
                velTS_c25_bl_db = velTS_BL_db;

                % Store trial-level velocity structs under explicit names
                velTS_trials_c25      = velocityData;
                velTS_db_trials_c25  = velocityData_db;
                velOCC_trials_c25     = velOCCFT;

                % Store pupil size
                pupTS_c25        = pupTS_noBL;
                pupTS_c25_bl_db = pupTS_BL_db;

                % Store MS
                msTS_c25        = msTS_noBL;
                msTS_c25_bl_db = msTS_BL_db;
                ms_events_c25 = struct();
                ms_events_c25.Onset = ms_onsets_per_trial;
                ms_events_c25.Offset = ms_offsets_per_trial;

            case 'c50'
                c50_bcea      = mean(bcea,'omitnan');
                c50_bl_bcea   = mean(baselineBcea,'omitnan');
                c50_pups      = mean(pupilSize,'omitnan');
                c50_bl_pups   = mean(baselinePupilSize,'omitnan');
                c50_msrate    = mean(microsaccadeRate,'omitnan');
                c50_bl_msrate = mean(baselineMSRate,'omitnan');

                % Velocity: raw and baseline means
                c50_velHorz      = mean(velHorz,'omitnan');
                c50_bl_velHorz   = mean(baselineVelH,'omitnan');
                c50_velVert      = mean(velVert,'omitnan');
                c50_bl_velVert   = mean(baselineVelV,'omitnan');
                c50_vel2D        = mean(vel2D,'omitnan');
                c50_bl_vel2D     = mean(baselineVel2D,'omitnan');

                % Condition scalars: trial means for MS/BCEA; TC means for vel/pupil
                c50_bcea_bl         = mean(BCEA_bl, 'omitnan');
                c50_bcea_bl_early   = mean(BCEA_bl_early, 'omitnan');
                c50_bcea_bl_late    = mean(BCEA_bl_late, 'omitnan');
                c50_pups_bl         = pup_tc_full;
                c50_pups_bl_early   = pup_tc_early;
                c50_pups_bl_late    = pup_tc_late;
                c50_msrate_bl       = mean(MSRate_bl, 'omitnan');
                c50_msrate_bl_early = mean(MSRate_bl_early, 'omitnan');
                c50_msrate_bl_late  = mean(MSRate_bl_late, 'omitnan');
                c50_velHorz_bl       = vel_tc_full_H;
                c50_velHorz_bl_early = vel_tc_early_H;
                c50_velHorz_bl_late  = vel_tc_late_H;
                c50_velVert_bl       = vel_tc_full_V;
                c50_velVert_bl_early = vel_tc_early_V;
                c50_velVert_bl_late  = vel_tc_late_V;
                c50_vel2D_bl         = vel_tc_full_2D;
                c50_vel2D_bl_early   = vel_tc_early_2D;
                c50_vel2D_bl_late    = vel_tc_late_2D;

                subj_data_gaze_trial_c50 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'BCEA',bcea, 'BCEA_early',bcea_early, 'BCEA_late',bcea_late, ...
                    'BaselineBCEA',baselineBcea, 'BCEA_bl', BCEA_bl, 'BCEA_bl_early', BCEA_bl_early, 'BCEA_bl_late', BCEA_bl_late, ...
                    'PupilSize',pupilSize, 'BaselinePupilSize',baselinePupilSize, 'PupilSize_bl', PupilSize_bl, ...
                    'MSRate',microsaccadeRate, 'MSRate_early',microsaccadeRate_early, 'MSRate_late',microsaccadeRate_late, ...
                    'BaselineMSRate',baselineMSRate, 'MSRate_bl', MSRate_bl, 'MSRate_bl_early', MSRate_bl_early, 'MSRate_bl_late', MSRate_bl_late, ...
                    'VelH',velHorz, 'BaselineVelH',baselineVelH, 'VelH_bl', VelH_bl, ...
                    'VelV',velVert, 'BaselineVelV',baselineVelV, 'VelV_bl', VelV_bl, ...
                    'Vel2D',vel2D, 'BaselineVel2D',baselineVel2D, 'Vel2D_bl', Vel2D_bl );

                velTS_c50        = velTS_noBL;
                velTS_c50_bl_db = velTS_BL_db;

                velTS_trials_c50      = velocityData;
                velTS_db_trials_c50  = velocityData_db;
                velOCC_trials_c50     = velOCCFT;

                pupTS_c50        = pupTS_noBL;
                pupTS_c50_bl_db = pupTS_BL_db;

                msTS_c50        = msTS_noBL;
                msTS_c50_bl_db = msTS_BL_db;
                ms_events_c50 = struct();
                ms_events_c50.Onset = ms_onsets_per_trial;
                ms_events_c50.Offset = ms_offsets_per_trial;

            case 'c75'
                c75_bcea      = mean(bcea,'omitnan');
                c75_bl_bcea   = mean(baselineBcea,'omitnan');
                c75_pups      = mean(pupilSize,'omitnan');
                c75_bl_pups   = mean(baselinePupilSize,'omitnan');
                c75_msrate    = mean(microsaccadeRate,'omitnan');
                c75_bl_msrate = mean(baselineMSRate,'omitnan');

                % Velocity: raw and baseline means
                c75_velHorz      = mean(velHorz,'omitnan');
                c75_bl_velHorz   = mean(baselineVelH,'omitnan');
                c75_velVert      = mean(velVert,'omitnan');
                c75_bl_velVert   = mean(baselineVelV,'omitnan');
                c75_vel2D        = mean(vel2D,'omitnan');
                c75_bl_vel2D     = mean(baselineVel2D,'omitnan');

                % Condition scalars: trial means for MS/BCEA; TC means for vel/pupil
                c75_bcea_bl         = mean(BCEA_bl, 'omitnan');
                c75_bcea_bl_early   = mean(BCEA_bl_early, 'omitnan');
                c75_bcea_bl_late    = mean(BCEA_bl_late, 'omitnan');
                c75_pups_bl         = pup_tc_full;
                c75_pups_bl_early   = pup_tc_early;
                c75_pups_bl_late    = pup_tc_late;
                c75_msrate_bl       = mean(MSRate_bl, 'omitnan');
                c75_msrate_bl_early = mean(MSRate_bl_early, 'omitnan');
                c75_msrate_bl_late  = mean(MSRate_bl_late, 'omitnan');
                c75_velHorz_bl       = vel_tc_full_H;
                c75_velHorz_bl_early = vel_tc_early_H;
                c75_velHorz_bl_late  = vel_tc_late_H;
                c75_velVert_bl       = vel_tc_full_V;
                c75_velVert_bl_early = vel_tc_early_V;
                c75_velVert_bl_late  = vel_tc_late_V;
                c75_vel2D_bl         = vel_tc_full_2D;
                c75_vel2D_bl_early   = vel_tc_early_2D;
                c75_vel2D_bl_late    = vel_tc_late_2D;

                subj_data_gaze_trial_c75 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'BCEA',bcea, 'BCEA_early',bcea_early, 'BCEA_late',bcea_late, ...
                    'BaselineBCEA',baselineBcea, 'BCEA_bl', BCEA_bl, 'BCEA_bl_early', BCEA_bl_early, 'BCEA_bl_late', BCEA_bl_late, ...
                    'PupilSize',pupilSize, 'BaselinePupilSize',baselinePupilSize, 'PupilSize_bl', PupilSize_bl, ...
                    'MSRate',microsaccadeRate, 'MSRate_early',microsaccadeRate_early, 'MSRate_late',microsaccadeRate_late, ...
                    'BaselineMSRate',baselineMSRate, 'MSRate_bl', MSRate_bl, 'MSRate_bl_early', MSRate_bl_early, 'MSRate_bl_late', MSRate_bl_late, ...
                    'VelH',velHorz, 'BaselineVelH',baselineVelH, 'VelH_bl', VelH_bl, ...
                    'VelV',velVert, 'BaselineVelV',baselineVelV, 'VelV_bl', VelV_bl, ...
                    'Vel2D',vel2D, 'BaselineVel2D',baselineVel2D, 'Vel2D_bl', Vel2D_bl );

                velTS_c75        = velTS_noBL;
                velTS_c75_bl_db = velTS_BL_db;

                velTS_trials_c75      = velocityData;
                velTS_db_trials_c75  = velocityData_db;
                velOCC_trials_c75     = velOCCFT;

                pupTS_c75        = pupTS_noBL;
                pupTS_c75_bl_db = pupTS_BL_db;

                msTS_c75        = msTS_noBL;
                msTS_c75_bl_db = msTS_BL_db;
                ms_events_c75 = struct();
                ms_events_c75.Onset = ms_onsets_per_trial;
                ms_events_c75.Offset = ms_offsets_per_trial;

            case 'c100'
                c100_bcea      = mean(bcea,'omitnan');
                c100_bl_bcea   = mean(baselineBcea,'omitnan');
                c100_pups      = mean(pupilSize,'omitnan');
                c100_bl_pups   = mean(baselinePupilSize,'omitnan');
                c100_msrate    = mean(microsaccadeRate,'omitnan');
                c100_bl_msrate = mean(baselineMSRate,'omitnan');

                % Velocity: raw and baseline means
                c100_velHorz      = mean(velHorz,'omitnan');
                c100_bl_velHorz   = mean(baselineVelH,'omitnan');
                c100_velVert      = mean(velVert,'omitnan');
                c100_bl_velVert   = mean(baselineVelV,'omitnan');
                c100_vel2D        = mean(vel2D,'omitnan');
                c100_bl_vel2D     = mean(baselineVel2D,'omitnan');

                % Condition scalars: trial means for MS/BCEA; TC means for vel/pupil
                c100_bcea_bl         = mean(BCEA_bl, 'omitnan');
                c100_bcea_bl_early   = mean(BCEA_bl_early, 'omitnan');
                c100_bcea_bl_late    = mean(BCEA_bl_late, 'omitnan');
                c100_pups_bl         = pup_tc_full;
                c100_pups_bl_early   = pup_tc_early;
                c100_pups_bl_late    = pup_tc_late;
                c100_msrate_bl       = mean(MSRate_bl, 'omitnan');
                c100_msrate_bl_early = mean(MSRate_bl_early, 'omitnan');
                c100_msrate_bl_late  = mean(MSRate_bl_late, 'omitnan');
                c100_velHorz_bl       = vel_tc_full_H;
                c100_velHorz_bl_early = vel_tc_early_H;
                c100_velHorz_bl_late  = vel_tc_late_H;
                c100_velVert_bl       = vel_tc_full_V;
                c100_velVert_bl_early = vel_tc_early_V;
                c100_velVert_bl_late  = vel_tc_late_V;
                c100_vel2D_bl         = vel_tc_full_2D;
                c100_vel2D_bl_early   = vel_tc_early_2D;
                c100_vel2D_bl_late    = vel_tc_late_2D;

                subj_data_gaze_trial_c100 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'BCEA',bcea, 'BCEA_early',bcea_early, 'BCEA_late',bcea_late, ...
                    'BaselineBCEA',baselineBcea, 'BCEA_bl', BCEA_bl, 'BCEA_bl_early', BCEA_bl_early, 'BCEA_bl_late', BCEA_bl_late, ...
                    'PupilSize',pupilSize, 'BaselinePupilSize',baselinePupilSize, 'PupilSize_bl', PupilSize_bl, ...
                    'MSRate',microsaccadeRate, 'MSRate_early',microsaccadeRate_early, 'MSRate_late',microsaccadeRate_late, ...
                    'BaselineMSRate',baselineMSRate, 'MSRate_bl', MSRate_bl, 'MSRate_bl_early', MSRate_bl_early, 'MSRate_bl_late', MSRate_bl_late, ...
                    'VelH',velHorz, 'BaselineVelH',baselineVelH, 'VelH_bl', VelH_bl, ...
                    'VelV',velVert, 'BaselineVelV',baselineVelV, 'VelV_bl', VelV_bl, ...
                    'Vel2D',vel2D, 'BaselineVel2D',baselineVel2D, 'Vel2D_bl', Vel2D_bl );

                velTS_c100        = velTS_noBL;
                velTS_c100_bl_db = velTS_BL_db;

                velTS_trials_c100      = velocityData;
                velTS_db_trials_c100  = velocityData_db;
                velOCC_trials_c100     = velOCCFT;

                pupTS_c100        = pupTS_noBL;
                pupTS_c100_bl_db = pupTS_BL_db;

                msTS_c100        = msTS_noBL;
                msTS_c100_bl_db = msTS_BL_db;
                ms_events_c100 = struct();
                ms_events_c100.Onset = ms_onsets_per_trial;
                ms_events_c100.Offset = ms_offsets_per_trial;
        end
    end

    %% CREATE SUBJECT‐LEVEL STRUCTS ACROSS CONDITIONS
    savepath = fullfile(paths.features, subjects{subj}, 'gaze');
    if ~exist(savepath,'dir'); mkdir(savepath); end

    % EyeLink event rates from preprocessing (full / early / late + baseline % change)
    load(fullfile(savepath, 'gaze_metrics'));
    c25_blinks_bl = c25_pct_blinks; c50_blinks_bl = c50_pct_blinks;
    c75_blinks_bl = c75_pct_blinks; c100_blinks_bl = c100_pct_blinks;
    c25_blinks_bl_early = c25_pct_blinks_early; c50_blinks_bl_early = c50_pct_blinks_early;
    c75_blinks_bl_early = c75_pct_blinks_early; c100_blinks_bl_early = c100_pct_blinks_early;
    c25_blinks_bl_late = c25_pct_blinks_late; c50_blinks_bl_late = c50_pct_blinks_late;
    c75_blinks_bl_late = c75_pct_blinks_late; c100_blinks_bl_late = c100_pct_blinks_late;

    c25_fixations_bl = c25_pct_fixations; c50_fixations_bl = c50_pct_fixations;
    c75_fixations_bl = c75_pct_fixations; c100_fixations_bl = c100_pct_fixations;
    c25_fixations_bl_early = c25_pct_fixations_early; c50_fixations_bl_early = c50_pct_fixations_early;
    c75_fixations_bl_early = c75_pct_fixations_early; c100_fixations_bl_early = c100_pct_fixations_early;
    c25_fixations_bl_late = c25_pct_fixations_late; c50_fixations_bl_late = c50_pct_fixations_late;
    c75_fixations_bl_late = c75_pct_fixations_late; c100_fixations_bl_late = c100_pct_fixations_late;

    c25_saccades_bl = c25_pct_saccades; c50_saccades_bl = c50_pct_saccades;
    c75_saccades_bl = c75_pct_saccades; c100_saccades_bl = c100_pct_saccades;
    c25_saccades_bl_early = c25_pct_saccades_early; c50_saccades_bl_early = c50_pct_saccades_early;
    c75_saccades_bl_early = c75_pct_saccades_early; c100_saccades_bl_early = c100_pct_saccades_early;
    c25_saccades_bl_late = c25_pct_saccades_late; c50_saccades_bl_late = c50_pct_saccades_late;
    c75_saccades_bl_late = c75_pct_saccades_late; c100_saccades_bl_late = c100_pct_saccades_late;

    subj_data_gaze = struct( ...
        'ID',        num2cell(subject_id(1:4))', ...
        'Condition', num2cell([1;2;3;4]), ...
        'BCEA',          num2cell([c25_bcea;   c50_bcea;   c75_bcea;   c100_bcea]), ...
        'PupilSize',     num2cell([c25_pups;   c50_pups;   c75_pups;   c100_pups]), ...
        'MSRate',        num2cell([c25_msrate; c50_msrate; c75_msrate; c100_msrate]), ...
        'VelH',          num2cell([c25_velHorz;   c50_velHorz;   c75_velHorz;   c100_velHorz]), ...
        'VelV',          num2cell([c25_velVert;   c50_velVert;   c75_velVert;   c100_velVert]), ...
        'Vel2D',         num2cell([c25_vel2D;  c50_vel2D;  c75_vel2D;  c100_vel2D]), ...
        'Blinks',        num2cell([c25_blinks; c50_blinks; c75_blinks; c100_blinks]), ...
        'Fixations',     num2cell([c25_fixations;c50_fixations;c75_fixations;c100_fixations]), ...
        'Saccades',      num2cell([c25_saccades; c50_saccades; c75_saccades; c100_saccades]), ...
        'BaselineBCEA',          num2cell([c25_bl_bcea;   c50_bl_bcea;   c75_bl_bcea;   c100_bl_bcea]), ...
        'BaselinePupilSize',     num2cell([c25_bl_pups;   c50_bl_pups;   c75_bl_pups;   c100_bl_pups]), ...
        'BaselineMSRate',        num2cell([c25_bl_msrate; c50_bl_msrate; c75_bl_msrate; c100_bl_msrate]), ...
        'BaselineVelH',          num2cell([c25_bl_velHorz;   c50_bl_velHorz;   c75_bl_velHorz;   c100_bl_velHorz]), ...
        'BaselineVelV',          num2cell([c25_bl_velVert;   c50_bl_velVert;   c75_bl_velVert;   c100_bl_velVert]), ...
        'BaselineVel2D',         num2cell([c25_bl_vel2D;  c50_bl_vel2D;  c75_bl_vel2D;  c100_bl_vel2D]), ...
        'BaselineBlinks',        num2cell([c25_bl_blinks; c50_bl_blinks; c75_bl_blinks; c100_bl_blinks]), ...
        'BaselineFixations',     num2cell([c25_bl_fixations; c50_bl_fixations; c75_bl_fixations; c100_bl_fixations]), ...
        'BaselineSaccades',      num2cell([c25_bl_saccades; c50_bl_saccades; c75_bl_saccades; c100_bl_saccades]), ...
        'BCEA_bl',          num2cell([c25_bcea_bl;   c50_bcea_bl;   c75_bcea_bl;   c100_bcea_bl]), ...
        'BCEA_bl_early',    num2cell([c25_bcea_bl_early; c50_bcea_bl_early; c75_bcea_bl_early; c100_bcea_bl_early]), ...
        'BCEA_bl_late',     num2cell([c25_bcea_bl_late;  c50_bcea_bl_late;  c75_bcea_bl_late;  c100_bcea_bl_late]), ...
        'PupilSize_bl',     num2cell([c25_pups_bl;   c50_pups_bl;   c75_pups_bl;   c100_pups_bl]), ...
        'PupilSize_bl_early', num2cell([c25_pups_bl_early; c50_pups_bl_early; c75_pups_bl_early; c100_pups_bl_early]), ...
        'PupilSize_bl_late',  num2cell([c25_pups_bl_late;  c50_pups_bl_late;  c75_pups_bl_late;  c100_pups_bl_late]), ...
        'MSRate_bl',        num2cell([c25_msrate_bl; c50_msrate_bl; c75_msrate_bl; c100_msrate_bl]), ...
        'MSRate_bl_early',  num2cell([c25_msrate_bl_early; c50_msrate_bl_early; c75_msrate_bl_early; c100_msrate_bl_early]), ...
        'MSRate_bl_late',   num2cell([c25_msrate_bl_late;  c50_msrate_bl_late;  c75_msrate_bl_late;  c100_msrate_bl_late]), ...
        'VelH_bl',          num2cell([c25_velHorz_bl;   c50_velHorz_bl;   c75_velHorz_bl;   c100_velHorz_bl]), ...
        'VelH_bl_early',    num2cell([c25_velHorz_bl_early; c50_velHorz_bl_early; c75_velHorz_bl_early; c100_velHorz_bl_early]), ...
        'VelH_bl_late',     num2cell([c25_velHorz_bl_late;  c50_velHorz_bl_late;  c75_velHorz_bl_late;  c100_velHorz_bl_late]), ...
        'VelV_bl',          num2cell([c25_velVert_bl;   c50_velVert_bl;   c75_velVert_bl;   c100_velVert_bl]), ...
        'VelV_bl_early',    num2cell([c25_velVert_bl_early; c50_velVert_bl_early; c75_velVert_bl_early; c100_velVert_bl_early]), ...
        'VelV_bl_late',     num2cell([c25_velVert_bl_late;  c50_velVert_bl_late;  c75_velVert_bl_late;  c100_velVert_bl_late]), ...
        'Vel2D_bl',         num2cell([c25_vel2D_bl;  c50_vel2D_bl;  c75_vel2D_bl;  c100_vel2D_bl]), ...
        'Vel2D_bl_early',   num2cell([c25_vel2D_bl_early; c50_vel2D_bl_early; c75_vel2D_bl_early; c100_vel2D_bl_early]), ...
        'Vel2D_bl_late',    num2cell([c25_vel2D_bl_late;  c50_vel2D_bl_late;  c75_vel2D_bl_late;  c100_vel2D_bl_late]), ...
        'Blinks_bl',        num2cell([c25_blinks_bl; c50_blinks_bl; c75_blinks_bl; c100_blinks_bl]), ...
        'Blinks_bl_early',  num2cell([c25_blinks_bl_early; c50_blinks_bl_early; c75_blinks_bl_early; c100_blinks_bl_early]), ...
        'Blinks_bl_late',   num2cell([c25_blinks_bl_late;  c50_blinks_bl_late;  c75_blinks_bl_late;  c100_blinks_bl_late]), ...
        'Fixations_bl',     num2cell([c25_fixations_bl; c50_fixations_bl; c75_fixations_bl; c100_fixations_bl]), ...
        'Fixations_bl_early', num2cell([c25_fixations_bl_early; c50_fixations_bl_early; c75_fixations_bl_early; c100_fixations_bl_early]), ...
        'Fixations_bl_late',  num2cell([c25_fixations_bl_late;  c50_fixations_bl_late;  c75_fixations_bl_late;  c100_fixations_bl_late]), ...
        'Saccades_bl',      num2cell([c25_saccades_bl; c50_saccades_bl; c75_saccades_bl; c100_saccades_bl]), ...
        'Saccades_bl_early', num2cell([c25_saccades_bl_early; c50_saccades_bl_early; c75_saccades_bl_early; c100_saccades_bl_early]), ...
        'Saccades_bl_late',  num2cell([c25_saccades_bl_late;  c50_saccades_bl_late;  c75_saccades_bl_late;  c100_saccades_bl_late]) );

    subj_data_gaze_baseline = struct( ...
        'ID',        num2cell(subject_id(1:4))', ...
        'Condition', num2cell([1;2;3;4]), ...
        'BaselineBCEA',          num2cell([c25_bl_bcea;   c50_bl_bcea;   c75_bl_bcea;   c100_bl_bcea]), ...
        'BaselinePupilSize',     num2cell([c25_bl_pups;   c50_bl_pups;   c75_bl_pups;   c100_bl_pups]), ...
        'BaselineMSRate',        num2cell([c25_bl_msrate; c50_bl_msrate; c75_bl_msrate; c100_bl_msrate]), ...
        'BaselineVelH',          num2cell([c25_bl_velHorz;   c50_bl_velHorz;   c75_bl_velHorz;   c100_bl_velHorz]), ...
        'BaselineVelV',          num2cell([c25_bl_velVert;   c50_bl_velVert;   c75_bl_velVert;   c100_bl_velVert]), ...
        'BaselineVel2D',         num2cell([c25_bl_vel2D;  c50_bl_vel2D;  c75_bl_vel2D;  c100_bl_vel2D]), ...
        'BaselineBlinks',        num2cell([c25_bl_blinks; c50_bl_blinks; c75_bl_blinks; c100_bl_blinks]), ...
        'BaselineFixations',     num2cell([c25_bl_fixations; c50_bl_fixations; c75_bl_fixations; c100_bl_fixations]), ...
        'BaselineSaccades',      num2cell([c25_bl_saccades; c50_bl_saccades; c75_bl_saccades; c100_bl_saccades]) );

    subj_data_gaze_bl = struct( ...
        'ID',        num2cell(subject_id(1:4))', ...
        'Condition', num2cell([1;2;3;4]), ...
        'BCEA_bl',          num2cell([c25_bcea_bl;   c50_bcea_bl;   c75_bcea_bl;   c100_bcea_bl]), ...
        'BCEA_bl_early',    num2cell([c25_bcea_bl_early; c50_bcea_bl_early; c75_bcea_bl_early; c100_bcea_bl_early]), ...
        'BCEA_bl_late',     num2cell([c25_bcea_bl_late;  c50_bcea_bl_late;  c75_bcea_bl_late;  c100_bcea_bl_late]), ...
        'PupilSize_bl',     num2cell([c25_pups_bl;   c50_pups_bl;   c75_pups_bl;   c100_pups_bl]), ...
        'PupilSize_bl_early', num2cell([c25_pups_bl_early; c50_pups_bl_early; c75_pups_bl_early; c100_pups_bl_early]), ...
        'PupilSize_bl_late',  num2cell([c25_pups_bl_late;  c50_pups_bl_late;  c75_pups_bl_late;  c100_pups_bl_late]), ...
        'MSRate_bl',        num2cell([c25_msrate_bl; c50_msrate_bl; c75_msrate_bl; c100_msrate_bl]), ...
        'MSRate_bl_early',  num2cell([c25_msrate_bl_early; c50_msrate_bl_early; c75_msrate_bl_early; c100_msrate_bl_early]), ...
        'MSRate_bl_late',   num2cell([c25_msrate_bl_late;  c50_msrate_bl_late;  c75_msrate_bl_late;  c100_msrate_bl_late]), ...
        'VelH_bl',          num2cell([c25_velHorz_bl;   c50_velHorz_bl;   c75_velHorz_bl;   c100_velHorz_bl]), ...
        'VelH_bl_early',    num2cell([c25_velHorz_bl_early; c50_velHorz_bl_early; c75_velHorz_bl_early; c100_velHorz_bl_early]), ...
        'VelH_bl_late',     num2cell([c25_velHorz_bl_late;  c50_velHorz_bl_late;  c75_velHorz_bl_late;  c100_velHorz_bl_late]), ...
        'VelV_bl',          num2cell([c25_velVert_bl;   c50_velVert_bl;   c75_velVert_bl;   c100_velVert_bl]), ...
        'VelV_bl_early',    num2cell([c25_velVert_bl_early; c50_velVert_bl_early; c75_velVert_bl_early; c100_velVert_bl_early]), ...
        'VelV_bl_late',     num2cell([c25_velVert_bl_late;  c50_velVert_bl_late;  c75_velVert_bl_late;  c100_velVert_bl_late]), ...
        'Vel2D_bl',         num2cell([c25_vel2D_bl;  c50_vel2D_bl;  c75_vel2D_bl;  c100_vel2D_bl]), ...
        'Vel2D_bl_early',   num2cell([c25_vel2D_bl_early; c50_vel2D_bl_early; c75_vel2D_bl_early; c100_vel2D_bl_early]), ...
        'Vel2D_bl_late',    num2cell([c25_vel2D_bl_late;  c50_vel2D_bl_late;  c75_vel2D_bl_late;  c100_vel2D_bl_late]), ...
        'Blinks_bl',        num2cell([c25_blinks_bl; c50_blinks_bl; c75_blinks_bl; c100_blinks_bl]), ...
        'Blinks_bl_early',  num2cell([c25_blinks_bl_early; c50_blinks_bl_early; c75_blinks_bl_early; c100_blinks_bl_early]), ...
        'Blinks_bl_late',   num2cell([c25_blinks_bl_late;  c50_blinks_bl_late;  c75_blinks_bl_late;  c100_blinks_bl_late]), ...
        'Fixations_bl',     num2cell([c25_fixations_bl; c50_fixations_bl; c75_fixations_bl; c100_fixations_bl]), ...
        'Fixations_bl_early', num2cell([c25_fixations_bl_early; c50_fixations_bl_early; c75_fixations_bl_early; c100_fixations_bl_early]), ...
        'Fixations_bl_late',  num2cell([c25_fixations_bl_late;  c50_fixations_bl_late;  c75_fixations_bl_late;  c100_fixations_bl_late]), ...
        'Saccades_bl',      num2cell([c25_saccades_bl; c50_saccades_bl; c75_saccades_bl; c100_saccades_bl]), ...
        'Saccades_bl_early', num2cell([c25_saccades_bl_early; c50_saccades_bl_early; c75_saccades_bl_early; c100_saccades_bl_early]), ...
        'Saccades_bl_late',  num2cell([c25_saccades_bl_late;  c50_saccades_bl_late;  c75_saccades_bl_late;  c100_saccades_bl_late]) );


    %% Save
    save(fullfile(savepath,'gaze_matrix_trial'),  ...
        'subj_data_gaze_trial_c25',  'subj_data_gaze_trial_c50',  ...
        'subj_data_gaze_trial_c75',  'subj_data_gaze_trial_c100');
    save(fullfile(savepath,'gaze_matrix_subj'),   'subj_data_gaze');
    save(fullfile(savepath,'gaze_matrix_baseline'),'subj_data_gaze_baseline');
    save(fullfile(savepath,'gaze_matrix_bl'),'subj_data_gaze_bl');

    % velocity time series and FieldTrip timelocked data
    save(fullfile(savepath, 'gaze_velocity_timeseries'), ...
        'velTS_c25','velTS_c50','velTS_c75','velTS_c100', ...
        'velTS_c25_bl_db','velTS_c50_bl_db','velTS_c75_bl_db','velTS_c100_bl_db', ...
        'velTS_trials_c25','velTS_trials_c50','velTS_trials_c75','velTS_trials_c100', ...
        'velTS_db_trials_c25','velTS_db_trials_c50','velTS_db_trials_c75','velTS_db_trials_c100', ...
        'velOCC_trials_c25','velOCC_trials_c50','velOCC_trials_c75','velOCC_trials_c100');

    % pupil size time series
    save(fullfile(savepath, 'gaze_pupil_timeseries'), ...
        'pupTS_c25','pupTS_c50','pupTS_c75','pupTS_c100', ...
        'pupTS_c25_bl_db','pupTS_c50_bl_db','pupTS_c75_bl_db','pupTS_c100_bl_db');

    % ms time series
    save(fullfile(savepath, 'gaze_microsaccade_timeseries'), ...
        'msTS_c25','msTS_c50','msTS_c75','msTS_c100', ...
        'msTS_c25_bl_db','msTS_c50_bl_db','msTS_c75_bl_db','msTS_c100_bl_db');

    % Engbert microsaccade onset/offset times (s) per trial
    save(fullfile(savepath, 'gaze_microsaccade_events'), ...
        'ms_events_c25', 'ms_events_c50', 'ms_events_c75', 'ms_events_c100');

    % Append to across-subjects raw struct
    if isempty(gaze_data)
        gaze_data = subj_data_gaze;
    else
        gaze_data = [gaze_data; subj_data_gaze]; %#ok<AGROW>
    end
end

%% Save files
save(fullfile(paths.features, 'GCP_gaze_raw.mat'), 'gaze_x_c25','gaze_y_c25','gaze_x_c50','gaze_y_c50','gaze_x_c75','gaze_y_c75','gaze_x_c100','gaze_y_c100');
save(fullfile(paths.features, 'GCP_gaze_matrix.mat'), 'gaze_data');
save_gaze_window_summaries(paths.features, subjects, gaze_data);
clc;
fprintf('[GCP] Gaze Fex done! %d/%d Subjects\n', subj, numel(subjects));


function [fullV, earlyV, lateV] = tc_window_means(avg, t, winFull, winEarly, winLate)
t = t(:)';
avg = avg(:)';
fullV  = mean(avg(t >= winFull(1)  & t <= winFull(2)),  'omitnan');
earlyV = mean(avg(t >= winEarly(1) & t <= winEarly(2)), 'omitnan');
lateV  = mean(avg(t >= winLate(1)  & t <= winLate(2)),  'omitnan');
end


function rate = ms_rate_in_window(raw, tVec, tw, win_size, fsample)
idx = tVec >= tw(1) & tVec <= tw(2);
dat = raw(1:3, idx);
valid = dat(1,:) >= 0 & dat(1,:) <= 800 & dat(2,:) >= 0 & dat(2,:) <= 600;
dat = dat(1:3, valid);
if isempty(dat)
    rate = NaN;
    return
end
dat(2,:) = 600 - dat(2,:);
dat = remove_blinks(dat, win_size);
x = dat(1,:); y = dat(2,:);
valid = isfinite(x) & isfinite(y);
x = x(valid); y = y(valid);
if numel(x) < 3
    rate = NaN;
    return
end
[rate, ~] = detect_microsaccades(fsample, [x; y], numel(x));
end

function val = bcea_in_window(raw, tVec, tw, win_size)
idx = tVec >= tw(1) & tVec <= tw(2);
dat = raw(1:3, idx);
valid = dat(1,:) >= 0 & dat(1,:) <= 800 & dat(2,:) >= 0 & dat(2,:) <= 600;
dat = dat(1:3, valid);
if isempty(dat)
    val = NaN;
    return
end
dat(2,:) = 600 - dat(2,:);
dat = remove_blinks(dat, win_size);
x = dat(1,:); y = dat(2,:);
if numel(x) < 3
    val = NaN;
    return
end
sx = nanstd(x); sy = nanstd(y);
rho = corr(x(:), y(:));
val = 2 * 2.291 * pi * sx * sy * sqrt(1 - rho^2);
end

function save_gaze_window_summaries(featuresRoot, subjects, gaze_data)
% Build cond x subject matrices for boxplots (no recomputation downstream).
metricBases = {'MSRate_bl','Vel2D_bl','PupilSize_bl','BCEA_bl', ...
        'Blinks_bl','Fixations_bl','Saccades_bl'};
winSuffix = {'', '_early', '_late'};
winName = {'full', 'early', 'late'};
nSubj = numel(subjects);
nCond = 4;
out = struct();
out.subjects = subjects;
out.windows = winName;
out.analysisWindows = struct('full', [0 2], 'early', [0 1], 'late', [1 2]);

G = struct2table(gaze_data);
for mi = 1:numel(metricBases)
    base = metricBases{mi};
    S = struct();
    for wi = 1:numel(winSuffix)
        field = [base winSuffix{wi}];
        M = nan(nCond, nSubj);
        if ~ismember(field, G.Properties.VariableNames)
            S.(winName{wi}) = M;
            continue
        end
        for s = 1:nSubj
            sid = str2double(subjects{s});
            for c = 1:nCond
                idx = G.ID == sid & G.Condition == c;
                if any(idx)
                    M(c, s) = G.(field)(find(idx, 1));
                end
            end
        end
        S.(winName{wi}) = M;
    end
    out.(base) = S;
end

outPath = fullfile(featuresRoot, 'GCP_gaze_window_summaries.mat');
save(outPath, '-struct', 'out');
fprintf('Saved %s\n', outPath);
end

function db = compute_db_baseline(stim, baseline)
% Power-style dB ratio: 10*log10(stim/baseline). Non-positive ratios -> NaN.
ratio = stim ./ baseline;
db = 10 * log10(ratio);
db(~isfinite(stim) | ~isfinite(baseline) | ~isfinite(ratio) | ratio <= 0) = NaN;
end

function pct = compute_pct_baseline(stim, baseline)
% Percentage change: 100*(stim-baseline)/baseline. Non-positive baselines -> NaN.
pct = 100 * (stim - baseline) ./ baseline;
pct(~isfinite(stim) | ~isfinite(baseline) | ~isfinite(pct) | baseline <= 0) = NaN;
end

function data = force_common_timeaxis(data, t_common)
% Put every trial on the same regular time axis so ft_timelockanalysis
% does not fail in ft_selectdata when trials are empty or unequal.
t_common = t_common(:)';
nTrials = numel(data.trial);
nChan = numel(data.label);
for trl = 1:nTrials
    t = data.time{trl};
    x = data.trial{trl};
    if isempty(t) || isempty(x) || size(x, 1) ~= nChan || size(x, 2) ~= numel(t)
        data.trial{trl} = nan(nChan, numel(t_common));
        data.time{trl} = t_common;
        continue
    end
    t = t(:)';
    if isequal(t, t_common)
        data.time{trl} = t_common;
        continue
    end
    newx = nan(nChan, numel(t_common));
    for ch = 1:nChan
        valid = isfinite(t) & isfinite(x(ch, :));
        if nnz(valid) >= 2
            newx(ch, :) = interp1(t(valid), x(ch, valid), t_common, 'linear', NaN);
        elseif nnz(valid) == 1
            [~, j] = min(abs(t_common - t(valid)));
            newx(ch, j) = x(ch, valid);
        end
    end
    data.trial{trl} = newx;
    data.time{trl} = t_common;
end
end
