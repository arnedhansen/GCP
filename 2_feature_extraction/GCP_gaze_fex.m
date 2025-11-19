%% GCP Gaze Feature Extraction Sternberg
%
% Extracted features:
%   Gaze deviation (Euclidean distances)
%   Gaze standard deviation
%   Pupil size
%   Microsaccades
%   Eye Velocity
%
% Gaze metrics labelled by eye-tracker (saccades, blinks and
% fixations) are extracted already in GCP_preprocessing.m

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');
% path = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/';
% dirs  = dir(path);
% folders = dirs([dirs.isdir] & ~ismember({dirs.name},{'.','..'}));
% subjects = {folders.name};

% preallocate across‐all‐subjects matrix of raw data
gaze_data = struct('ID',{},'Condition',{},'GazeDeviation',{},...
    'GazeStdX',{},'GazeStdY',{},'PupilSize',{},...
    'MSRate',{}, ...
    'VelH',{},'VelV',{},'Vel2D',{}, ...
    'Blinks',{},'Fixations',{},'Saccades',{});

gaze_data_bl = struct('ID',{},'Condition',{},'PctGazeDeviation',{},...
    'PctGazeStdX',{},'PctGazeStdY',{},'PctPupilSize',{},...
    'PctMSRate',{}, ...
    'PctVelH',{},'PctVelV',{},'PctVel2D',{});

% time‐windows
baseline_period = [-1.5, -0.25];
analysis_period = [ 0.3,  2   ];
win_size        = 50;    % blink‐removal window (samples)
fsample         = 500;   % eye‐tracker sampling rate

% prepare raw gaze storage across all subjects
gaze_x_c25   = {};  gaze_y_c25   = {};
gaze_x_c50   = {};  gaze_y_c50   = {};
gaze_x_c75   = {};  gaze_y_c75   = {};
gaze_x_c100  = {};  gaze_y_c100  = {};

%% Loop over subjects
clc
for subj = 1:numel(subjects)

    % load preprocessed eye‐tracker data
    datapath = fullfile(path, subjects{subj}, 'gaze');
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
        pupilSize          = [];  baselinePupilSize = [];
        microsaccadeRate   = [];  baselineMSRate    = [];

        % Prepare velocity time series container for this condition
        velocityData     = dataET;    % copy meta-info
        velocityData.label = {'VelH','VelV','Vel2D'};

        velocityData_pct = dataET;    % baseline-normalised (% change)
        velocityData_pct.label = {'PctVelH','PctVelV','PctVel2D'};

        fsample          = dataET.fsample;
        vel_win_size     = round(0.1 * fsample);   % 100 ms window
        nTrials          = numel(dataET.trial);

        % For linking baseline and trial-wise means
        velH        = nan(1, nTrials);
        velV        = nan(1, nTrials);
        vel2D       = nan(1, nTrials);
        baselineVelH  = nan(1, nTrials);
        baselineVelV  = nan(1, nTrials);
        baselineVel2D = nan(1, nTrials);


        %% Trial loop
        for trl = 1:numel(dataET.trialinfo)
            raw    = dataET.trial{trl};
            tVec   = dataET.time{trl};

            % BASELINE WINDOW
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
            baseline_pupil    = mean(bl_dat(3,:),'omitnan')/1000;
            [baseline_msrate, ~] = detect_microsaccades(fsample, [bl_x; bl_y], numel(bl_x));

            % Baseline eye velocity (full baseline window)
            if numel(bl_x) > 1
                bl_x_sm = movmean(bl_x, vel_win_size);
                bl_y_sm = movmean(bl_y, vel_win_size);

                vx_bl    = diff(bl_x_sm) * fsample;        % px/s
                vy_bl    = diff(bl_y_sm) * fsample;        % px/s
                speed_bl = sqrt(vx_bl.^2 + vy_bl.^2);

                baselineVelH(trl)  = mean(abs(vx_bl), 'omitnan');
                baselineVelV(trl)  = mean(abs(vy_bl), 'omitnan');
                baselineVel2D(trl) = mean(speed_bl,   'omitnan');
            else
                baselineVelH(trl)  = NaN;
                baselineVelV(trl)  = NaN;
                baselineVel2D(trl) = NaN;
            end

            % ANALYSIS WINDOW
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
            pupil       = mean(an_dat(3,:),'omitnan')/1000;
            [msrate, ~] = detect_microsaccades(fsample, [x; y], numel(x));

            % Velocity time series in analysis window [0.3, 2.0] s
            if numel(x) > 1
                x_sm = movmean(x, vel_win_size);
                y_sm = movmean(y, vel_win_size);

                vx = diff(x_sm) * fsample;      % px/s
                vy = diff(y_sm) * fsample;      % px/s
                speed = sqrt(vx.^2 + vy.^2);    % 2D speed

                % Time axis for velocity (one sample shorter than position)
                t_an   = tVec(an_idx);          % analysis window times
                t_vel  = t_an(1:numel(speed));  % trim to match diff()

                % Store raw velocity time series
                velocityData.trial{trl} = zeros(3, numel(speed));
                velocityData.trial{trl}(1,:) = abs(vx);
                velocityData.trial{trl}(2,:) = abs(vy);
                velocityData.trial{trl}(3,:) = speed;
                velocityData.time{trl}       = t_vel;

                % Trial-wise mean velocities (for scalar features)
                velH(trl)   = mean(abs(vx), 'omitnan');
                velV(trl)   = mean(abs(vy), 'omitnan');
                vel2D(trl)  = mean(speed,   'omitnan');

                % Baseline-normalised velocity time series (% change)
                if ~isnan(baselineVelH(trl))
                    pct_vx    = (abs(vx) ./ baselineVelH(trl)  - 1) * 100;
                else
                    pct_vx    = nan(size(vx));
                end
                if ~isnan(baselineVelV(trl))
                    pct_vy    = (abs(vy) ./ baselineVelV(trl)  - 1) * 100;
                else
                    pct_vy    = nan(size(vy));
                end
                if ~isnan(baselineVel2D(trl))
                    pct_speed = (speed   ./ baselineVel2D(trl) - 1) * 100;
                else
                    pct_speed = nan(size(speed));
                end

                velocityData_pct.trial{trl} = zeros(3, numel(speed));
                velocityData_pct.trial{trl}(1,:) = pct_vx;
                velocityData_pct.trial{trl}(2,:) = pct_vy;
                velocityData_pct.trial{trl}(3,:) = pct_speed;
                velocityData_pct.time{trl}       = t_vel;

            else
                % Too few points → fill with NaNs
                velocityData.trial{trl}      = nan(3,0);
                velocityData.time{trl}       = [];
                velocityData_pct.trial{trl}  = nan(3,0);
                velocityData_pct.time{trl}   = [];
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
            pupilSize(end+1)        = pupil;
            baselinePupilSize(end+1)= baseline_pupil;
            microsaccadeRate(end+1) = msrate;
            baselineMSRate(end+1)   = baseline_msrate;
        end

        %% SUBJECT‐BY‐CONDITION AVERAGES
        switch cond
            case 'c25'
                c25_gdev      = mean(gazeDev,'omitnan');
                c25_bl_gdev   = mean(baselineGazeDev,'omitnan');
                c25_gSDx      = mean(gazeSDx,'omitnan');
                c25_bl_gSDx   = mean(baselineGazeSDx,'omitnan');
                c25_gSDy      = mean(gazeSDy,'omitnan');
                c25_bl_gSDy   = mean(baselineGazeSDy,'omitnan');
                c25_pups      = mean(pupilSize,'omitnan');
                c25_bl_pups   = mean(baselinePupilSize,'omitnan');
                c25_msrate    = mean(microsaccadeRate,'omitnan');
                c25_bl_msrate = mean(baselineMSRate,'omitnan');

                % Velocity: raw and baseline means
                c25_velH      = mean(velH,'omitnan');
                c25_bl_velH   = mean(baselineVelH,'omitnan');
                c25_velV      = mean(velV,'omitnan');
                c25_bl_velV   = mean(baselineVelV,'omitnan');
                c25_vel2D     = mean(vel2D,'omitnan');
                c25_bl_vel2D  = mean(baselineVel2D,'omitnan');

                % Percentage change relative to baseline
                c25_pct_gdev    = (c25_gdev    - c25_bl_gdev   ) / c25_bl_gdev   * 100;
                c25_pct_gSDx    = (c25_gSDx    - c25_bl_gSDx   ) / c25_bl_gSDx   * 100;
                c25_pct_gSDy    = (c25_gSDy    - c25_bl_gSDy   ) / c25_bl_gSDy   * 100;
                c25_pct_pups    = (c25_pups    - c25_bl_pups   ) / c25_bl_pups   * 100;
                c25_pct_msrate  = (c25_msrate  - c25_bl_msrate ) / c25_bl_msrate * 100;

                c25_pct_velH    = (c25_velH    - c25_bl_velH   ) / c25_bl_velH   * 100;
                c25_pct_velV    = (c25_velV    - c25_bl_velV   ) / c25_bl_velV   * 100;
                c25_pct_vel2D   = (c25_vel2D   - c25_bl_vel2D  ) / c25_bl_vel2D  * 100;

                subj_data_gaze_trial_c25 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'PctGazeDeviation',(gazeDev./baselineGazeDev-1)*100, ...
                    'GazeStdX',gazeSDx,        'BaselineGazeStdX',baselineGazeSDx,        'PctGazeStdX',      (gazeSDx./baselineGazeSDx-1)*100, ...
                    'GazeStdY',gazeSDy,        'BaselineGazeStdY',baselineGazeSDy,        'PctGazeStdY',      (gazeSDy./baselineGazeSDy-1)*100, ...
                    'PupilSize',pupilSize,     'BaselinePupilSize',baselinePupilSize,     'PctPupilSize',     (pupilSize./baselinePupilSize-1)*100, ...
                    'MSRate',microsaccadeRate, 'BaselineMSRate',baselineMSRate,           'PctMSRate',        (microsaccadeRate./baselineMSRate-1)*100, ...
                    'VelH',velH,               'BaselineVelH',baselineVelH,               'PctVelH',          (velH./baselineVelH-1)*100, ...
                    'VelV',velV,               'BaselineVelV',baselineVelV,               'PctVelV',          (velV./baselineVelV-1)*100, ...
                    'Vel2D',vel2D,             'BaselineVel2D',baselineVel2D,             'PctVel2D',         (vel2D./baselineVel2D-1)*100 );

            case 'c50'
                c50_gdev      = mean(gazeDev,'omitnan');
                c50_bl_gdev   = mean(baselineGazeDev,'omitnan');
                c50_gSDx      = mean(gazeSDx,'omitnan');
                c50_bl_gSDx   = mean(baselineGazeSDx,'omitnan');
                c50_gSDy      = mean(gazeSDy,'omitnan');
                c50_bl_gSDy   = mean(baselineGazeSDy,'omitnan');
                c50_pups      = mean(pupilSize,'omitnan');
                c50_bl_pups   = mean(baselinePupilSize,'omitnan');
                c50_msrate    = mean(microsaccadeRate,'omitnan');
                c50_bl_msrate = mean(baselineMSRate,'omitnan');

                % Velocity: raw and baseline means
                c50_velH      = mean(velH,'omitnan');
                c50_bl_velH   = mean(baselineVelH,'omitnan');
                c50_velV      = mean(velV,'omitnan');
                c50_bl_velV   = mean(baselineVelV,'omitnan');
                c50_vel2D     = mean(vel2D,'omitnan');
                c50_bl_vel2D  = mean(baselineVel2D,'omitnan');

                % Percentage change relative to baseline
                c50_pct_gdev    = (c50_gdev    - c50_bl_gdev   ) / c50_bl_gdev   * 100;
                c50_pct_gSDx    = (c50_gSDx    - c50_bl_gSDx   ) / c50_bl_gSDx   * 100;
                c50_pct_gSDy    = (c50_gSDy    - c50_bl_gSDy   ) / c50_bl_gSDy   * 100;
                c50_pct_pups    = (c50_pups    - c50_bl_pups   ) / c50_bl_pups   * 100;
                c50_pct_msrate  = (c50_msrate  - c50_bl_msrate ) / c50_bl_msrate * 100;

                c50_pct_velH    = (c50_velH    - c50_bl_velH   ) / c50_bl_velH   * 100;
                c50_pct_velV    = (c50_velV    - c50_bl_velV   ) / c50_bl_velV   * 100;
                c50_pct_vel2D   = (c50_vel2D   - c50_bl_vel2D  ) / c50_bl_vel2D  * 100;

                subj_data_gaze_trial_c50 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'PctGazeDeviation',(gazeDev./baselineGazeDev-1)*100, ...
                    'GazeStdX',gazeSDx,        'BaselineGazeStdX',baselineGazeSDx,        'PctGazeStdX',      (gazeSDx./baselineGazeSDx-1)*100, ...
                    'GazeStdY',gazeSDy,        'BaselineGazeStdY',baselineGazeSDy,        'PctGazeStdY',      (gazeSDy./baselineGazeSDy-1)*100, ...
                    'PupilSize',pupilSize,     'BaselinePupilSize',baselinePupilSize,     'PctPupilSize',     (pupilSize./baselinePupilSize-1)*100, ...
                    'MSRate',microsaccadeRate, 'BaselineMSRate',baselineMSRate,           'PctMSRate',        (microsaccadeRate./baselineMSRate-1)*100, ...
                    'VelH',velH,               'BaselineVelH',baselineVelH,               'PctVelH',          (velH./baselineVelH-1)*100, ...
                    'VelV',velV,               'BaselineVelV',baselineVelV,               'PctVelV',          (velV./baselineVelV-1)*100, ...
                    'Vel2D',vel2D,             'BaselineVel2D',baselineVel2D,             'PctVel2D',         (vel2D./baselineVel2D-1)*100 );

            case 'c75'
                c75_gdev      = mean(gazeDev,'omitnan');
                c75_bl_gdev   = mean(baselineGazeDev,'omitnan');
                c75_gSDx      = mean(gazeSDx,'omitnan');
                c75_bl_gSDx   = mean(baselineGazeSDx,'omitnan');
                c75_gSDy      = mean(gazeSDy,'omitnan');
                c75_bl_gSDy   = mean(baselineGazeSDy,'omitnan');
                c75_pups      = mean(pupilSize,'omitnan');
                c75_bl_pups   = mean(baselinePupilSize,'omitnan');
                c75_msrate    = mean(microsaccadeRate,'omitnan');
                c75_bl_msrate = mean(baselineMSRate,'omitnan');

                % Velocity: raw and baseline means
                c75_velH      = mean(velH,'omitnan');
                c75_bl_velH   = mean(baselineVelH,'omitnan');
                c75_velV      = mean(velV,'omitnan');
                c75_bl_velV   = mean(baselineVelV,'omitnan');
                c75_vel2D     = mean(vel2D,'omitnan');
                c75_bl_vel2D  = mean(baselineVel2D,'omitnan');

                % Percentage change relative to baseline
                c75_pct_gdev    = (c75_gdev    - c75_bl_gdev   ) / c75_bl_gdev   * 100;
                c75_pct_gSDx    = (c75_gSDx    - c75_bl_gSDx   ) / c75_bl_gSDx   * 100;
                c75_pct_gSDy    = (c75_gSDy    - c75_bl_gSDy   ) / c75_bl_gSDy   * 100;
                c75_pct_pups    = (c75_pups    - c75_bl_pups   ) / c75_bl_pups   * 100;
                c75_pct_msrate  = (c75_msrate  - c75_bl_msrate ) / c75_bl_msrate * 100;

                c75_pct_velH    = (c75_velH    - c75_bl_velH   ) / c75_bl_velH   * 100;
                c75_pct_velV    = (c75_velV    - c75_bl_velV   ) / c75_bl_velV   * 100;
                c75_pct_vel2D   = (c75_vel2D   - c75_bl_vel2D  ) / c75_bl_vel2D  * 100;

                subj_data_gaze_trial_c75 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'PctGazeDeviation',(gazeDev./baselineGazeDev-1)*100, ...
                    'GazeStdX',gazeSDx,        'BaselineGazeStdX',baselineGazeSDx,        'PctGazeStdX',      (gazeSDx./baselineGazeSDx-1)*100, ...
                    'GazeStdY',gazeSDy,        'BaselineGazeStdY',baselineGazeSDy,        'PctGazeStdY',      (gazeSDy./baselineGazeSDy-1)*100, ...
                    'PupilSize',pupilSize,     'BaselinePupilSize',baselinePupilSize,     'PctPupilSize',     (pupilSize./baselinePupilSize-1)*100, ...
                    'MSRate',microsaccadeRate, 'BaselineMSRate',baselineMSRate,           'PctMSRate',        (microsaccadeRate./baselineMSRate-1)*100, ...
                    'VelH',velH,               'BaselineVelH',baselineVelH,               'PctVelH',          (velH./baselineVelH-1)*100, ...
                    'VelV',velV,               'BaselineVelV',baselineVelV,               'PctVelV',          (velV./baselineVelV-1)*100, ...
                    'Vel2D',vel2D,             'BaselineVel2D',baselineVel2D,             'PctVel2D',         (vel2D./baselineVel2D-1)*100 );

            case 'c100'
                c100_gdev      = mean(gazeDev,'omitnan');
                c100_bl_gdev   = mean(baselineGazeDev,'omitnan');
                c100_gSDx      = mean(gazeSDx,'omitnan');
                c100_bl_gSDx   = mean(baselineGazeSDx,'omitnan');
                c100_gSDy      = mean(gazeSDy,'omitnan');
                c100_bl_gSDy   = mean(baselineGazeSDy,'omitnan');
                c100_pups      = mean(pupilSize,'omitnan');
                c100_bl_pups   = mean(baselinePupilSize,'omitnan');
                c100_msrate    = mean(microsaccadeRate,'omitnan');
                c100_bl_msrate = mean(baselineMSRate,'omitnan');

                % Velocity: raw and baseline means
                c100_velH      = mean(velH,'omitnan');
                c100_bl_velH   = mean(baselineVelH,'omitnan');
                c100_velV      = mean(velV,'omitnan');
                c100_bl_velV   = mean(baselineVelV,'omitnan');
                c100_vel2D     = mean(vel2D,'omitnan');
                c100_bl_vel2D  = mean(baselineVel2D,'omitnan');

                % Percentage change relative to baseline
                c100_pct_gdev    = (c100_gdev    - c100_bl_gdev   ) / c100_bl_gdev   * 100;
                c100_pct_gSDx    = (c100_gSDx    - c100_bl_gSDx   ) / c100_bl_gSDx   * 100;
                c100_pct_gSDy    = (c100_gSDy    - c100_bl_gSDy   ) / c100_bl_gSDy   * 100;
                c100_pct_pups    = (c100_pups    - c100_bl_pups   ) / c100_bl_pups   * 100;
                c100_pct_msrate  = (c100_msrate  - c100_bl_msrate ) / c100_bl_msrate * 100;

                c100_pct_velH    = (c100_velH    - c100_bl_velH   ) / c100_bl_velH   * 100;
                c100_pct_velV    = (c100_velV    - c100_bl_velV   ) / c100_bl_velV   * 100;
                c100_pct_vel2D   = (c100_vel2D   - c100_bl_vel2D  ) / c100_bl_vel2D  * 100;

                subj_data_gaze_trial_c100 = struct( ...
                    'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                    'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'PctGazeDeviation',(gazeDev./baselineGazeDev-1)*100, ...
                    'GazeStdX',gazeSDx,        'BaselineGazeStdX',baselineGazeSDx,        'PctGazeStdX',      (gazeSDx./baselineGazeSDx-1)*100, ...
                    'GazeStdY',gazeSDy,        'BaselineGazeStdY',baselineGazeSDy,        'PctGazeStdY',      (gazeSDy./baselineGazeSDy-1)*100, ...
                    'PupilSize',pupilSize,     'BaselinePupilSize',baselinePupilSize,     'PctPupilSize',     (pupilSize./baselinePupilSize-1)*100, ...
                    'MSRate',microsaccadeRate, 'BaselineMSRate',baselineMSRate,           'PctMSRate',        (microsaccadeRate./baselineMSRate-1)*100, ...
                    'VelH',velH,               'BaselineVelH',baselineVelH,               'PctVelH',          (velH./baselineVelH-1)*100, ...
                    'VelV',velV,               'BaselineVelV',baselineVelV,               'PctVelV',          (velV./baselineVelV-1)*100, ...
                    'Vel2D',vel2D,             'BaselineVel2D',baselineVel2D,             'PctVel2D',         (vel2D./baselineVel2D-1)*100 );
        end

        % Store eye velocity time series data
        switch cond
            case 'c25'
                velTS_c25      = velocityData;
                velTS_pct_c25  = velocityData_pct;
            case 'c50'
                velTS_c50      = velocityData;
                velTS_pct_c50  = velocityData_pct;
            case 'c75'
                velTS_c75      = velocityData;
                velTS_pct_c75  = velocityData_pct;
            case 'c100'
                velTS_c100     = velocityData;
                velTS_pct_c100 = velocityData_pct;
        end
    end

    %% CREATE SUBJECT‐LEVEL STRUCTS ACROSS CONDITIONS
    savepath = fullfile(path, subjects{subj}, 'gaze');
    if ~exist(savepath,'dir'); mkdir(savepath); end

    % load the pre‐extracted GCP_preprocessing metrics
    load(fullfile(savepath,'gaze_metrics'),'c25_blinks','c50_blinks','c75_blinks','c100_blinks', ...
        'c25_fixations','c50_fixations','c75_fixations','c100_fixations', ...
        'c25_saccades','c50_saccades','c75_saccades','c100_saccades');

    subj_data_gaze = struct( ...
        'ID',        num2cell(subject_id(1:4))', ...
        'Condition', num2cell([1;2;3;4]), ...
        'GazeDeviation', num2cell([c25_gdev;   c50_gdev;   c75_gdev;   c100_gdev]), ...
        'GazeStdX',      num2cell([c25_gSDx;   c50_gSDx;   c75_gSDx;   c100_gSDx]), ...
        'GazeStdY',      num2cell([c25_gSDy;   c50_gSDy;   c75_gSDy;   c100_gSDy]), ...
        'PupilSize',     num2cell([c25_pups;   c50_pups;   c75_pups;   c100_pups]), ...
        'MSRate',        num2cell([c25_msrate; c50_msrate; c75_msrate; c100_msrate]), ...
        'VelH',          num2cell([c25_velH;   c50_velH;   c75_velH;   c100_velH]), ...
        'VelV',          num2cell([c25_velV;   c50_velV;   c75_velV;   c100_velV]), ...
        'Vel2D',         num2cell([c25_vel2D;  c50_vel2D;  c75_vel2D;  c100_vel2D]), ...
        'Blinks',        num2cell([c25_blinks; c50_blinks; c75_blinks; c100_blinks]), ...
        'Fixations',     num2cell([c25_fixations;c50_fixations;c75_fixations;c100_fixations]), ...
        'Saccades',      num2cell([c25_saccades; c50_saccades; c75_saccades; c100_saccades]) );

    subj_data_gaze_baseline = struct( ...
        'ID',        num2cell(subject_id(1:4))', ...
        'Condition', num2cell([1;2;3;4]), ...
        'BaselineGazeDeviation', num2cell([c25_bl_gdev;   c50_bl_gdev;   c75_bl_gdev;   c100_bl_gdev]), ...
        'BaselineGazeStdX',      num2cell([c25_bl_gSDx;   c50_bl_gSDx;   c75_bl_gSDx;   c100_bl_gSDx]), ...
        'BaselineGazeStdY',      num2cell([c25_bl_gSDy;   c50_bl_gSDy;   c75_bl_gSDy;   c100_bl_gSDy]), ...
        'BaselinePupilSize',     num2cell([c25_bl_pups;   c50_bl_pups;   c75_bl_pups;   c100_bl_pups]), ...
        'BaselineMSRate',        num2cell([c25_bl_msrate; c50_bl_msrate; c75_bl_msrate; c100_bl_msrate]), ...
        'BaselineVelH',          num2cell([c25_bl_velH;   c50_bl_velH;   c75_bl_velH;   c100_bl_velH]), ...
        'BaselineVelV',          num2cell([c25_bl_velV;   c50_bl_velV;   c75_bl_velV;   c100_bl_velV]), ...
        'BaselineVel2D',         num2cell([c25_bl_vel2D;  c50_bl_vel2D;  c75_bl_vel2D;  c100_bl_vel2D]) );

    subj_data_gaze_pctchange = struct( ...
        'ID',        num2cell(subject_id(1:4))', ...
        'Condition', num2cell([1;2;3;4]), ...
        'PctGazeDeviation', num2cell([c25_pct_gdev;   c50_pct_gdev;   c75_pct_gdev;   c100_pct_gdev]), ...
        'PctGazeStdX',      num2cell([c25_pct_gSDx;   c50_pct_gSDx;   c75_pct_gSDx;   c100_pct_gSDx]), ...
        'PctGazeStdY',      num2cell([c25_pct_gSDy;   c50_pct_gSDy;   c75_pct_gSDy;   c100_pct_gSDy]), ...
        'PctPupilSize',     num2cell([c25_pct_pups;   c50_pct_pups;   c75_pct_pups;   c100_pct_pups]), ...
        'PctMSRate',        num2cell([c25_pct_msrate; c50_pct_msrate; c75_pct_msrate; c100_pct_msrate]), ...
        'PctVelH',          num2cell([c25_pct_velH;   c50_pct_velH;   c75_pct_velH;   c100_pct_velH]), ...
        'PctVelV',          num2cell([c25_pct_velV;   c50_pct_velV;   c75_pct_velV;   c100_pct_velV]), ...
        'PctVel2D',         num2cell([c25_pct_vel2D;  c50_pct_vel2D;  c75_pct_vel2D;  c100_pct_vel2D]) );


    %% Save
    save(fullfile(savepath,'gaze_matrix_trial'),  ...
        'subj_data_gaze_trial_c25',  'subj_data_gaze_trial_c50',  ...
        'subj_data_gaze_trial_c75',  'subj_data_gaze_trial_c100');
    save(fullfile(savepath,'gaze_matrix_subj'),   'subj_data_gaze');
    save(fullfile(savepath,'gaze_matrix_baseline'),'subj_data_gaze_baseline');
    save(fullfile(savepath,'gaze_matrix_pctchange'),'subj_data_gaze_pctchange');
    save(fullfile(savepath, 'gaze_velocity_timeseries'), ...
     'velTS_c25','velTS_c50','velTS_c75','velTS_c100', ...
     'velTS_pct_c25','velTS_pct_c50','velTS_pct_c75','velTS_pct_c100');

    % Append to across‐subjects raw struct
    gaze_data = [gaze_data; subj_data_gaze];
    gaze_data_bl = [gaze_data_bl; subj_data_gaze_pctchange];

    fprintf('Subject %d/%d done.\n', subj, numel(subjects));
end

%% Save files
save(fullfile(path,'gaze_raw'),    'gaze_x_c25','gaze_y_c25','gaze_x_c50','gaze_y_c50','gaze_x_c75','gaze_y_c75','gaze_x_c100','gaze_y_c100');
save(fullfile(path,'gaze_matrix'), 'gaze_data');
save(fullfile(path,'gaze_matrix_bl'), 'gaze_data_bl');