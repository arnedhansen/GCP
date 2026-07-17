%% GCP Gaze Feature Extraction
%
% Extracted features:
%   Gaze deviation (Euclidean distances)
%   Gaze standard deviation
%   BCEA (Bivariate Contour Ellipse Area, k=2.291, 95%)
%   Pupil size (Time Series Raw, Baselined % change)
%   Microsaccades (TC ground truth: boxplot scalar = mean of % TC over [0, 2] s)
%   Saccade cleaned fixational eye velocity (Raw, Baselined % change)
%
% Gaze metrics labelled by eye-tracker (saccades, blinks and
% fixations) are extracted already in GCP_preprocessing.m
%
% Note: fields named dB* store percentage change for pupil, MS, velocity,
% and BCEA (downstream name compatibility).

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('GCP');

% preallocate across‐all‐subjects matrix of raw data
gaze_data = struct('ID',{},'Condition',{},'GazeDeviation',{},...
    'GazeStdX',{},'GazeStdY',{},'BCEA',{},'PupilSize',{},...
    'MSRate',{}, ...
    'VelH',{},'VelV',{},'Vel2D',{}, ...
    'Blinks',{},'Fixations',{},'Saccades',{}, ...
    'BaselineGazeDeviation',{},'BaselineGazeStdX',{},'BaselineGazeStdY',{},...
    'BaselineBCEA',{},'BaselinePupilSize',{},'BaselineMSRate',{},...
    'BaselineVelH',{},'BaselineVelV',{},'BaselineVel2D',{}, ...
    'BaselineBlinks',{},'BaselineFixations',{},'BaselineSaccades',{}, ...
    'dBGazeDeviation',{},'dBGazeStdX',{},'dBGazeStdY',{},...
    'dBBCEA',{},'dBPupilSize',{},'dBMSRate',{},...
    'dBVelH',{},'dBVelV',{},'dBVel2D',{}, ...
    'dBBlinks',{},'dBFixations',{},'dBSaccades',{});

% time‐windows
baseline_period   = [-1.5, -0.25];
analysis_period   = [0 2];
analysis_periodTS = [-1 2];
pupil_store_window = [-1.5, 2.5];
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

        gazeDev            = [];  baselineGazeDev   = [];
        gazeSDx            = [];  baselineGazeSDx   = [];
        gazeSDy            = [];  baselineGazeSDy   = [];
        bcea               = [];  baselineBcea      = [];
        pupilSize          = [];  baselinePupilSize = [];
        microsaccadeRate   = [];  baselineMSRate    = [];

        % Prepare velocity time series container for this condition
        velocityData     = dataET;    % copy meta-info
        velocityData.label = {'VelH','VelV','Vel2D'};

        % trial-level % baseline (scalar baseline per trial; field names keep *_db)
        velocityData_db = dataET;
        velocityData_db.label = {'dBVelH','dBVelV','dBVel2D'};

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

        % Continuous rectified velocity for OCC: saccadic samples are zeroed
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

        %% Trial loop
        for trl = 1:numel(dataET.trialinfo)
            raw    = dataET.trial{trl};
            tVec   = dataET.time{trl};

            % Saccade cleaned fixational velocity on the original regular time axis.
            % No temporal smoothing is applied before differentiation because a
            % 50 to 100 ms moving average would suppress activity in the gamma range.
            full_idx = tVec >= baseline_period(1) & tVec <= analysis_period(2);
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

                % Detect rapid saccadic events on the same position trace.
                [~, saccade_details_vel] = detect_microsaccades( ...
                    fsample, full_dat(1:2,:), numel(t_full));
                saccade_mask = false(1, numel(t_full));
                for sac = 1:numel(saccade_details_vel.Onset)
                    sac_start = max(1, saccade_details_vel.Onset(sac));
                    sac_end   = min(numel(t_full), saccade_details_vel.Offset(sac));
                    saccade_mask(sac_start:sac_end) = true;
                end

                vel_signed(:, ~valid_full_position) = NaN;

                % Exclude saccadic samples from mean fixational velocity so the
                % estimate is not mechanically lowered when more saccades occur.
                vel_fixational = vel_signed;
                vel_fixational(:, saccade_mask) = NaN;
                vx_full    = abs(vel_fixational(1,:));
                vy_full    = abs(vel_fixational(2,:));
                speed_full = sqrt(sum(vel_fixational.^2, 1));

                velFT.trial{trl} = [vx_full; vy_full; speed_full];
                velFT.time{trl}  = t_full;

                % OCC requires a continuous rectified signal. Following the
                % proposal, detected saccadic samples are therefore set to zero.
                vel_occ = vel_signed;
                vel_occ(:, saccade_mask) = 0;
                velOCCFT.trial{trl} = [abs(vel_occ(1,:)); ...
                    abs(vel_occ(2,:)); sqrt(sum(vel_occ.^2, 1))];
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
            dx_bl = bl_x - 400 - nanmean(bl_x-400);
            dy_bl = bl_y - 300 - nanmean(bl_y-300);
            baseline_eucdev   = nanmean( sqrt(dx_bl.^2 + dy_bl.^2) );
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

            % Baseline fixational velocity after removal of saccadic samples
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

            % Compute GazeDev & Microsaccades
            dx = x - 400 - nanmean(x-400);
            dy = y - 300 - nanmean(y-300);

            mean_eucdev = nanmean( sqrt(dx.^2 + dy.^2) );
            std_x       = nanstd(x);
            std_y       = nanstd(y);
            if numel(x) > 2
                rho_an = corr(x(:), y(:));
            else
                rho_an = 0;
            end
            bcea_val    = 2 * 2.291 * pi * std_x * std_y * sqrt(1 - rho_an^2);
            pupil       = mean(an_dat(3,:),'omitnan')/1000;
            [msrate, ~] = detect_microsaccades(fsample, [x; y], numel(x));

            % Analysis window fixational velocity from the full regular trace
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
                msrate = NaN;
                baseline_msrate = NaN;
            else
                x_clean = x_val(valid_clean);
                y_clean = y_val(valid_clean);

                % Detect microsaccades in cleaned valid samples only
                [~, ms_det] = detect_microsaccades(fsample, [x_clean; y_clean], numel(x_clean));

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
                    msrate = NaN;
                    baseline_msrate = NaN;
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
                    msFT.trial{trl} = ms_rate_full;
                    msFT.time{trl}  = t_full_ms;

                    % Raw rates from the same TC used for plotting / boxplots
                    bl_ms_idx = t_full_ms >= baseline_period(1) & t_full_ms <= baseline_period(2);
                    an_ms_idx = t_full_ms >= analysis_period(1) & t_full_ms <= analysis_period(2);
                    baseline_msrate = mean(ms_rate_full(bl_ms_idx), 'omitnan');
                    msrate = mean(ms_rate_full(an_ms_idx), 'omitnan');
                end
            end

            % append to trial‐wise arrays
            subject_id(end+1)       = str2double(subjects{subj});
            trial_num(end+1)        = trl;
            condition(end+1)        = dataET.trialinfo(trl) - 60;
            gazeDev(end+1)          = mean_eucdev;
            baselineGazeDev(end+1)  = baseline_eucdev;
            gazeSDx(end+1)          = std_x;
            baselineGazeSDx(end+1)  = baseline_std_x;
            gazeSDy(end+1)          = std_y;
            baselineGazeSDy(end+1)  = baseline_std_y;
            bcea(end+1)             = bcea_val;
            baselineBcea(end+1)     = baseline_bcea_val;
            pupilSize(end+1)        = pupil;
            baselinePupilSize(end+1)= baseline_pupil;
            microsaccadeRate(end+1) = msrate;
            baselineMSRate(end+1)   = baseline_msrate;
        end

        % FieldTrip timelock and baseline on velFT
        % Force identical time axes: empty/unequal trials break ft_selectdata.
        t_vel_common = baseline_period(1):(1 / fsample):analysis_period(2);
        velFT = force_common_timeaxis(velFT, t_vel_common);
        velOCCFT = force_common_timeaxis(velOCCFT, t_vel_common);

        % Timelocked average without baseline
        cfg = [];
        cfg.latency    = analysis_periodTS;   % [-1 2]
        cfg.keeptrials = 'no';
        velTS_noBL = ft_timelockanalysis(cfg, velFT);

        % % change: scalar trial baseline 100*(x(t)/mean_base - 1)
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

        cfg = [];
        cfg.latency    = analysis_periodTS;
        cfg.keeptrials = 'no';
        velTS_BL_db = ft_timelockanalysis(cfg, velFT_bl_db);

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
        msFT = force_common_timeaxis(msFT, t_vel_common);

        cfg = [];
        cfg.latency    = analysis_periodTS;
        cfg.keeptrials = 'no';
        msTS_noBL = ft_timelockanalysis(cfg, msFT);

        % % change on MS TC; trial scalar = mean of that TC over [0, 2] s
        msFT_bl_db = msFT;
        dBMSRate = nan(1, nTrials);
        for trl = 1:nTrials
            if all(~isfinite(msFT.trial{trl}(:)))
                continue
            end
            msFT_bl_db.trial{trl} = compute_pct_baseline(msFT.trial{trl}, baselineMSRate(trl));
            t_ms = msFT_bl_db.time{trl};
            an_ms = t_ms >= analysis_period(1) & t_ms <= analysis_period(2);
            dBMSRate(trl) = mean(msFT_bl_db.trial{trl}(an_ms), 'omitnan');
        end
        msFT_bl_db = force_common_timeaxis(msFT_bl_db, t_vel_common);

        cfg = [];
        cfg.latency    = analysis_periodTS;
        cfg.keeptrials = 'no';
        msTS_BL_db = ft_timelockanalysis(cfg, msFT_bl_db);

        % Subject MS boxplot value: mean of the saved/ploted TC over [0, 2] s
        ms_an_idx = msTS_BL_db.time >= analysis_period(1) & msTS_BL_db.time <= analysis_period(2);
        ms_tc_summary = mean(msTS_BL_db.avg(:, ms_an_idx), 'omitnan');

        % Trial-level baselines for scalar gaze metrics
        dBGazeDeviation = compute_db_baseline(gazeDev, baselineGazeDev);
        dBGazeStdX      = compute_db_baseline(gazeSDx, baselineGazeSDx);
        dBGazeStdY      = compute_db_baseline(gazeSDy, baselineGazeSDy);
        dBBCEA          = compute_pct_baseline(bcea, baselineBcea);
        dBPupilSize     = compute_pct_baseline(pupilSize, baselinePupilSize);
        dBVelH          = compute_pct_baseline(velHorz, baselineVelH);
        dBVelV          = compute_pct_baseline(velVert, baselineVelV);
        dBVel2D         = compute_pct_baseline(vel2D, baselineVel2D);

        %% SUBJECT‐BY‐CONDITION AVERAGES
        switch cond
            case 'c25'
                c25_gdev      = mean(gazeDev,'omitnan');
                c25_bl_gdev   = mean(baselineGazeDev,'omitnan');
                c25_gSDx      = mean(gazeSDx,'omitnan');
                c25_bl_gSDx   = mean(baselineGazeSDx,'omitnan');
                c25_gSDy      = mean(gazeSDy,'omitnan');
                c25_bl_gSDy   = mean(baselineGazeSDy,'omitnan');
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

                % Condition mean of trial-level dB metrics
                c25_db_gdev    = mean(dBGazeDeviation, 'omitnan');
                c25_db_gSDx    = mean(dBGazeStdX, 'omitnan');
                c25_db_gSDy    = mean(dBGazeStdY, 'omitnan');
                c25_db_bcea    = mean(dBBCEA, 'omitnan');
                c25_db_pups    = mean(dBPupilSize, 'omitnan');
                c25_db_msrate  = ms_tc_summary;
                c25_db_velHorz = mean(dBVelH, 'omitnan');
                c25_db_velVert = mean(dBVelV, 'omitnan');
                c25_db_vel2D   = mean(dBVel2D, 'omitnan');

                subj_data_gaze_trial_c25 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'dBGazeDeviation',dBGazeDeviation, ...
                    'GazeStdX',gazeSDx,        'BaselineGazeStdX',baselineGazeSDx,        'dBGazeStdX',      dBGazeStdX, ...
                    'GazeStdY',gazeSDy,        'BaselineGazeStdY',baselineGazeSDy,        'dBGazeStdY',      dBGazeStdY, ...
                    'BCEA',bcea,               'BaselineBCEA',baselineBcea,               'dBBCEA',          dBBCEA, ...
                    'PupilSize',pupilSize,     'BaselinePupilSize',baselinePupilSize,     'dBPupilSize',     dBPupilSize, ...
                    'MSRate',microsaccadeRate, 'BaselineMSRate',baselineMSRate,           'dBMSRate',        dBMSRate, ...
                    'VelH',velHorz,            'BaselineVelH',baselineVelH,               'dBVelH',          dBVelH, ...
                    'VelV',velVert,            'BaselineVelV',baselineVelV,               'dBVelV',          dBVelV, ...
                    'Vel2D',vel2D,             'BaselineVel2D',baselineVel2D,             'dBVel2D',         dBVel2D );

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

            case 'c50'
                c50_gdev      = mean(gazeDev,'omitnan');
                c50_bl_gdev   = mean(baselineGazeDev,'omitnan');
                c50_gSDx      = mean(gazeSDx,'omitnan');
                c50_bl_gSDx   = mean(baselineGazeSDx,'omitnan');
                c50_gSDy      = mean(gazeSDy,'omitnan');
                c50_bl_gSDy   = mean(baselineGazeSDy,'omitnan');
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

                % Condition mean of trial-level dB metrics
                c50_db_gdev    = mean(dBGazeDeviation, 'omitnan');
                c50_db_gSDx    = mean(dBGazeStdX, 'omitnan');
                c50_db_gSDy    = mean(dBGazeStdY, 'omitnan');
                c50_db_bcea    = mean(dBBCEA, 'omitnan');
                c50_db_pups    = mean(dBPupilSize, 'omitnan');
                c50_db_msrate  = ms_tc_summary;
                c50_db_velHorz = mean(dBVelH, 'omitnan');
                c50_db_velVert = mean(dBVelV, 'omitnan');
                c50_db_vel2D   = mean(dBVel2D, 'omitnan');

                subj_data_gaze_trial_c50 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'dBGazeDeviation',dBGazeDeviation, ...
                    'GazeStdX',gazeSDx,        'BaselineGazeStdX',baselineGazeSDx,        'dBGazeStdX',      dBGazeStdX, ...
                    'GazeStdY',gazeSDy,        'BaselineGazeStdY',baselineGazeSDy,        'dBGazeStdY',      dBGazeStdY, ...
                    'BCEA',bcea,               'BaselineBCEA',baselineBcea,               'dBBCEA',          dBBCEA, ...
                    'PupilSize',pupilSize,     'BaselinePupilSize',baselinePupilSize,     'dBPupilSize',     dBPupilSize, ...
                    'MSRate',microsaccadeRate, 'BaselineMSRate',baselineMSRate,           'dBMSRate',        dBMSRate, ...
                    'VelH',velHorz,            'BaselineVelH',baselineVelH,               'dBVelH',          dBVelH, ...
                    'VelV',velVert,            'BaselineVelV',baselineVelV,               'dBVelV',          dBVelV, ...
                    'Vel2D',vel2D,             'BaselineVel2D',baselineVel2D,             'dBVel2D',         dBVel2D );

                velTS_c50        = velTS_noBL;
                velTS_c50_bl_db = velTS_BL_db;

                velTS_trials_c50      = velocityData;
                velTS_db_trials_c50  = velocityData_db;
                velOCC_trials_c50     = velOCCFT;

                pupTS_c50        = pupTS_noBL;
                pupTS_c50_bl_db = pupTS_BL_db;

                msTS_c50        = msTS_noBL;
                msTS_c50_bl_db = msTS_BL_db;

            case 'c75'
                c75_gdev      = mean(gazeDev,'omitnan');
                c75_bl_gdev   = mean(baselineGazeDev,'omitnan');
                c75_gSDx      = mean(gazeSDx,'omitnan');
                c75_bl_gSDx   = mean(baselineGazeSDx,'omitnan');
                c75_gSDy      = mean(gazeSDy,'omitnan');
                c75_bl_gSDy   = mean(baselineGazeSDy,'omitnan');
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

                % Condition mean of trial-level dB metrics
                c75_db_gdev    = mean(dBGazeDeviation, 'omitnan');
                c75_db_gSDx    = mean(dBGazeStdX, 'omitnan');
                c75_db_gSDy    = mean(dBGazeStdY, 'omitnan');
                c75_db_bcea    = mean(dBBCEA, 'omitnan');
                c75_db_pups    = mean(dBPupilSize, 'omitnan');
                c75_db_msrate  = ms_tc_summary;
                c75_db_velHorz = mean(dBVelH, 'omitnan');
                c75_db_velVert = mean(dBVelV, 'omitnan');
                c75_db_vel2D   = mean(dBVel2D, 'omitnan');

                subj_data_gaze_trial_c75 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'dBGazeDeviation',dBGazeDeviation, ...
                    'GazeStdX',gazeSDx,        'BaselineGazeStdX',baselineGazeSDx,        'dBGazeStdX',      dBGazeStdX, ...
                    'GazeStdY',gazeSDy,        'BaselineGazeStdY',baselineGazeSDy,        'dBGazeStdY',      dBGazeStdY, ...
                    'BCEA',bcea,               'BaselineBCEA',baselineBcea,               'dBBCEA',          dBBCEA, ...
                    'PupilSize',pupilSize,     'BaselinePupilSize',baselinePupilSize,     'dBPupilSize',     dBPupilSize, ...
                    'MSRate',microsaccadeRate, 'BaselineMSRate',baselineMSRate,           'dBMSRate',        dBMSRate, ...
                    'VelH',velHorz,            'BaselineVelH',baselineVelH,               'dBVelH',          dBVelH, ...
                    'VelV',velVert,            'BaselineVelV',baselineVelV,               'dBVelV',          dBVelV, ...
                    'Vel2D',vel2D,             'BaselineVel2D',baselineVel2D,             'dBVel2D',         dBVel2D );

                velTS_c75        = velTS_noBL;
                velTS_c75_bl_db = velTS_BL_db;

                velTS_trials_c75      = velocityData;
                velTS_db_trials_c75  = velocityData_db;
                velOCC_trials_c75     = velOCCFT;

                pupTS_c75        = pupTS_noBL;
                pupTS_c75_bl_db = pupTS_BL_db;

                msTS_c75        = msTS_noBL;
                msTS_c75_bl_db = msTS_BL_db;

            case 'c100'
                c100_gdev      = mean(gazeDev,'omitnan');
                c100_bl_gdev   = mean(baselineGazeDev,'omitnan');
                c100_gSDx      = mean(gazeSDx,'omitnan');
                c100_bl_gSDx   = mean(baselineGazeSDx,'omitnan');
                c100_gSDy      = mean(gazeSDy,'omitnan');
                c100_bl_gSDy   = mean(baselineGazeSDy,'omitnan');
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

                % Condition mean of trial-level dB metrics
                c100_db_gdev    = mean(dBGazeDeviation, 'omitnan');
                c100_db_gSDx    = mean(dBGazeStdX, 'omitnan');
                c100_db_gSDy    = mean(dBGazeStdY, 'omitnan');
                c100_db_bcea    = mean(dBBCEA, 'omitnan');
                c100_db_pups    = mean(dBPupilSize, 'omitnan');
                c100_db_msrate  = ms_tc_summary;
                c100_db_velHorz = mean(dBVelH, 'omitnan');
                c100_db_velVert = mean(dBVelV, 'omitnan');
                c100_db_vel2D   = mean(dBVel2D, 'omitnan');

                subj_data_gaze_trial_c100 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'dBGazeDeviation',dBGazeDeviation, ...
                    'GazeStdX',gazeSDx,        'BaselineGazeStdX',baselineGazeSDx,        'dBGazeStdX',      dBGazeStdX, ...
                    'GazeStdY',gazeSDy,        'BaselineGazeStdY',baselineGazeSDy,        'dBGazeStdY',      dBGazeStdY, ...
                    'BCEA',bcea,               'BaselineBCEA',baselineBcea,               'dBBCEA',          dBBCEA, ...
                    'PupilSize',pupilSize,     'BaselinePupilSize',baselinePupilSize,     'dBPupilSize',     dBPupilSize, ...
                    'MSRate',microsaccadeRate, 'BaselineMSRate',baselineMSRate,           'dBMSRate',        dBMSRate, ...
                    'VelH',velHorz,            'BaselineVelH',baselineVelH,               'dBVelH',          dBVelH, ...
                    'VelV',velVert,            'BaselineVelV',baselineVelV,               'dBVelV',          dBVelV, ...
                    'Vel2D',vel2D,             'BaselineVel2D',baselineVel2D,             'dBVel2D',         dBVel2D );

                velTS_c100        = velTS_noBL;
                velTS_c100_bl_db = velTS_BL_db;

                velTS_trials_c100      = velocityData;
                velTS_db_trials_c100  = velocityData_db;
                velOCC_trials_c100     = velOCCFT;

                pupTS_c100        = pupTS_noBL;
                pupTS_c100_bl_db = pupTS_BL_db;

                msTS_c100        = msTS_noBL;
                msTS_c100_bl_db = msTS_BL_db;
        end
    end

    %% CREATE SUBJECT‐LEVEL STRUCTS ACROSS CONDITIONS
    savepath = fullfile(paths.features, subjects{subj}, 'gaze');
    if ~exist(savepath,'dir'); mkdir(savepath); end

    % load the pre‐extracted GCP_preprocessing metrics
    load(fullfile(savepath,'gaze_metrics'),'c25_blinks','c50_blinks','c75_blinks','c100_blinks', ...
        'c25_fixations','c50_fixations','c75_fixations','c100_fixations', ...
        'c25_saccades','c50_saccades','c75_saccades','c100_saccades', ...
        'c25_bl_blinks','c50_bl_blinks','c75_bl_blinks','c100_bl_blinks', ...
        'c25_bl_fixations','c50_bl_fixations','c75_bl_fixations','c100_bl_fixations', ...
        'c25_bl_saccades','c50_bl_saccades','c75_bl_saccades','c100_bl_saccades');

    c25_db_blinks     = compute_db_baseline(c25_blinks, c25_bl_blinks);
    c50_db_blinks     = compute_db_baseline(c50_blinks, c50_bl_blinks);
    c75_db_blinks     = compute_db_baseline(c75_blinks, c75_bl_blinks);
    c100_db_blinks    = compute_db_baseline(c100_blinks, c100_bl_blinks);
    c25_db_fixations  = compute_db_baseline(c25_fixations, c25_bl_fixations);
    c50_db_fixations  = compute_db_baseline(c50_fixations, c50_bl_fixations);
    c75_db_fixations  = compute_db_baseline(c75_fixations, c75_bl_fixations);
    c100_db_fixations = compute_db_baseline(c100_fixations, c100_bl_fixations);
    c25_db_saccades   = compute_db_baseline(c25_saccades, c25_bl_saccades);
    c50_db_saccades   = compute_db_baseline(c50_saccades, c50_bl_saccades);
    c75_db_saccades   = compute_db_baseline(c75_saccades, c75_bl_saccades);
    c100_db_saccades  = compute_db_baseline(c100_saccades, c100_bl_saccades);

    subj_data_gaze = struct( ...
        'ID',        num2cell(subject_id(1:4))', ...
        'Condition', num2cell([1;2;3;4]), ...
        'GazeDeviation', num2cell([c25_gdev;   c50_gdev;   c75_gdev;   c100_gdev]), ...
        'GazeStdX',      num2cell([c25_gSDx;   c50_gSDx;   c75_gSDx;   c100_gSDx]), ...
        'GazeStdY',      num2cell([c25_gSDy;   c50_gSDy;   c75_gSDy;   c100_gSDy]), ...
        'BCEA',          num2cell([c25_bcea;   c50_bcea;   c75_bcea;   c100_bcea]), ...
        'PupilSize',     num2cell([c25_pups;   c50_pups;   c75_pups;   c100_pups]), ...
        'MSRate',        num2cell([c25_msrate; c50_msrate; c75_msrate; c100_msrate]), ...
        'VelH',          num2cell([c25_velHorz;   c50_velHorz;   c75_velHorz;   c100_velHorz]), ...
        'VelV',          num2cell([c25_velVert;   c50_velVert;   c75_velVert;   c100_velVert]), ...
        'Vel2D',         num2cell([c25_vel2D;  c50_vel2D;  c75_vel2D;  c100_vel2D]), ...
        'Blinks',        num2cell([c25_blinks; c50_blinks; c75_blinks; c100_blinks]), ...
        'Fixations',     num2cell([c25_fixations;c50_fixations;c75_fixations;c100_fixations]), ...
        'Saccades',      num2cell([c25_saccades; c50_saccades; c75_saccades; c100_saccades]), ...
        'BaselineGazeDeviation', num2cell([c25_bl_gdev;   c50_bl_gdev;   c75_bl_gdev;   c100_bl_gdev]), ...
        'BaselineGazeStdX',      num2cell([c25_bl_gSDx;   c50_bl_gSDx;   c75_bl_gSDx;   c100_bl_gSDx]), ...
        'BaselineGazeStdY',      num2cell([c25_bl_gSDy;   c50_bl_gSDy;   c75_bl_gSDy;   c100_bl_gSDy]), ...
        'BaselineBCEA',          num2cell([c25_bl_bcea;   c50_bl_bcea;   c75_bl_bcea;   c100_bl_bcea]), ...
        'BaselinePupilSize',     num2cell([c25_bl_pups;   c50_bl_pups;   c75_bl_pups;   c100_bl_pups]), ...
        'BaselineMSRate',        num2cell([c25_bl_msrate; c50_bl_msrate; c75_bl_msrate; c100_bl_msrate]), ...
        'BaselineVelH',          num2cell([c25_bl_velHorz;   c50_bl_velHorz;   c75_bl_velHorz;   c100_bl_velHorz]), ...
        'BaselineVelV',          num2cell([c25_bl_velVert;   c50_bl_velVert;   c75_bl_velVert;   c100_bl_velVert]), ...
        'BaselineVel2D',         num2cell([c25_bl_vel2D;  c50_bl_vel2D;  c75_bl_vel2D;  c100_bl_vel2D]), ...
        'BaselineBlinks',        num2cell([c25_bl_blinks; c50_bl_blinks; c75_bl_blinks; c100_bl_blinks]), ...
        'BaselineFixations',     num2cell([c25_bl_fixations; c50_bl_fixations; c75_bl_fixations; c100_bl_fixations]), ...
        'BaselineSaccades',      num2cell([c25_bl_saccades; c50_bl_saccades; c75_bl_saccades; c100_bl_saccades]), ...
        'dBGazeDeviation', num2cell([c25_db_gdev;   c50_db_gdev;   c75_db_gdev;   c100_db_gdev]), ...
        'dBGazeStdX',      num2cell([c25_db_gSDx;   c50_db_gSDx;   c75_db_gSDx;   c100_db_gSDx]), ...
        'dBGazeStdY',      num2cell([c25_db_gSDy;   c50_db_gSDy;   c75_db_gSDy;   c100_db_gSDy]), ...
        'dBBCEA',          num2cell([c25_db_bcea;   c50_db_bcea;   c75_db_bcea;   c100_db_bcea]), ...
        'dBPupilSize',     num2cell([c25_db_pups;   c50_db_pups;   c75_db_pups;   c100_db_pups]), ...
        'dBMSRate',        num2cell([c25_db_msrate; c50_db_msrate; c75_db_msrate; c100_db_msrate]), ...
        'dBVelH',          num2cell([c25_db_velHorz;   c50_db_velHorz;   c75_db_velHorz;   c100_db_velHorz]), ...
        'dBVelV',          num2cell([c25_db_velVert;   c50_db_velVert;   c75_db_velVert;   c100_db_velVert]), ...
        'dBVel2D',         num2cell([c25_db_vel2D;  c50_db_vel2D;  c75_db_vel2D;  c100_db_vel2D]), ...
        'dBBlinks',        num2cell([c25_db_blinks; c50_db_blinks; c75_db_blinks; c100_db_blinks]), ...
        'dBFixations',     num2cell([c25_db_fixations; c50_db_fixations; c75_db_fixations; c100_db_fixations]), ...
        'dBSaccades',      num2cell([c25_db_saccades; c50_db_saccades; c75_db_saccades; c100_db_saccades]) );

    subj_data_gaze_baseline = struct( ...
        'ID',        num2cell(subject_id(1:4))', ...
        'Condition', num2cell([1;2;3;4]), ...
        'BaselineGazeDeviation', num2cell([c25_bl_gdev;   c50_bl_gdev;   c75_bl_gdev;   c100_bl_gdev]), ...
        'BaselineGazeStdX',      num2cell([c25_bl_gSDx;   c50_bl_gSDx;   c75_bl_gSDx;   c100_bl_gSDx]), ...
        'BaselineGazeStdY',      num2cell([c25_bl_gSDy;   c50_bl_gSDy;   c75_bl_gSDy;   c100_bl_gSDy]), ...
        'BaselineBCEA',          num2cell([c25_bl_bcea;   c50_bl_bcea;   c75_bl_bcea;   c100_bl_bcea]), ...
        'BaselinePupilSize',     num2cell([c25_bl_pups;   c50_bl_pups;   c75_bl_pups;   c100_bl_pups]), ...
        'BaselineMSRate',        num2cell([c25_bl_msrate; c50_bl_msrate; c75_bl_msrate; c100_bl_msrate]), ...
        'BaselineVelH',          num2cell([c25_bl_velHorz;   c50_bl_velHorz;   c75_bl_velHorz;   c100_bl_velHorz]), ...
        'BaselineVelV',          num2cell([c25_bl_velVert;   c50_bl_velVert;   c75_bl_velVert;   c100_bl_velVert]), ...
        'BaselineVel2D',         num2cell([c25_bl_vel2D;  c50_bl_vel2D;  c75_bl_vel2D;  c100_bl_vel2D]), ...
        'BaselineBlinks',        num2cell([c25_bl_blinks; c50_bl_blinks; c75_bl_blinks; c100_bl_blinks]), ...
        'BaselineFixations',     num2cell([c25_bl_fixations; c50_bl_fixations; c75_bl_fixations; c100_bl_fixations]), ...
        'BaselineSaccades',      num2cell([c25_bl_saccades; c50_bl_saccades; c75_bl_saccades; c100_bl_saccades]) );

    subj_data_gaze_dbchange = struct( ...
        'ID',        num2cell(subject_id(1:4))', ...
        'Condition', num2cell([1;2;3;4]), ...
        'dBGazeDeviation', num2cell([c25_db_gdev;   c50_db_gdev;   c75_db_gdev;   c100_db_gdev]), ...
        'dBGazeStdX',      num2cell([c25_db_gSDx;   c50_db_gSDx;   c75_db_gSDx;   c100_db_gSDx]), ...
        'dBGazeStdY',      num2cell([c25_db_gSDy;   c50_db_gSDy;   c75_db_gSDy;   c100_db_gSDy]), ...
        'dBBCEA',          num2cell([c25_db_bcea;   c50_db_bcea;   c75_db_bcea;   c100_db_bcea]), ...
        'dBPupilSize',     num2cell([c25_db_pups;   c50_db_pups;   c75_db_pups;   c100_db_pups]), ...
        'dBMSRate',        num2cell([c25_db_msrate; c50_db_msrate; c75_db_msrate; c100_db_msrate]), ...
        'dBVelH',          num2cell([c25_db_velHorz;   c50_db_velHorz;   c75_db_velHorz;   c100_db_velHorz]), ...
        'dBVelV',          num2cell([c25_db_velVert;   c50_db_velVert;   c75_db_velVert;   c100_db_velVert]), ...
        'dBVel2D',         num2cell([c25_db_vel2D;  c50_db_vel2D;  c75_db_vel2D;  c100_db_vel2D]), ...
        'dBBlinks',        num2cell([c25_db_blinks; c50_db_blinks; c75_db_blinks; c100_db_blinks]), ...
        'dBFixations',     num2cell([c25_db_fixations; c50_db_fixations; c75_db_fixations; c100_db_fixations]), ...
        'dBSaccades',      num2cell([c25_db_saccades; c50_db_saccades; c75_db_saccades; c100_db_saccades]) );


    %% Save
    save(fullfile(savepath,'gaze_matrix_trial'),  ...
        'subj_data_gaze_trial_c25',  'subj_data_gaze_trial_c50',  ...
        'subj_data_gaze_trial_c75',  'subj_data_gaze_trial_c100');
    save(fullfile(savepath,'gaze_matrix_subj'),   'subj_data_gaze');
    save(fullfile(savepath,'gaze_matrix_baseline'),'subj_data_gaze_baseline');
    save(fullfile(savepath,'gaze_matrix_dbchange'),'subj_data_gaze_dbchange');

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

    % Append to across‐subjects raw struct
    gaze_data = [gaze_data; subj_data_gaze];
end

%% Save files
save(fullfile(paths.features, 'GCP_gaze_raw.mat'), 'gaze_x_c25','gaze_y_c25','gaze_x_c50','gaze_y_c50','gaze_x_c75','gaze_y_c75','gaze_x_c100','gaze_y_c100');
save(fullfile(paths.features, 'GCP_gaze_matrix.mat'), 'gaze_data');
clc;
fprintf('[GCP] Gaze Fex done! %d/%d Subjects\n', subj, numel(subjects));

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