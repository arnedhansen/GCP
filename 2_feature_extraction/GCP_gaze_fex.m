%% GCP Gaze Feature Extraction Sternberg
%
% Extracted features:
%   Gaze deviation (Euclidean distances)
%   Gaze standard deviation
%   Pupil size
%   Microsaccades
%
% Gaze metrics labelled by eye-tracker (saccades, blinks and
% fixations) are extracted already in GCP_preprocessing.m

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/features/';
dirs  = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name},{'.','..'}));
subjects = {folders.name};

% preallocate across‐all‐subjects matrix of raw data
gaze_data = struct('ID',{},'Condition',{},'GazeDeviation',{},...
                   'GazeStdX',{},'GazeStdY',{},'PupilSize',{},...
                   'MSRate',{},'Blinks',{},'Fixations',{},'Saccades',{});
gaze_data_bl = struct('ID',{},'Condition',{},'PctGazeDeviation',{},...
                   'PctGazeStdX',{},'PctGazeStdY',{},'PctPupilSize',{},...
                   'PctMSRate',{});

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
for subj = 1:numel(subjects)
    
    % load preprocessed eye‐tracker data
    datapath = fullfile(path, subjects{subj}, 'gaze');
    load(fullfile(datapath,'dataET'));
    
    %% Loop over conditions
    for conds = {'c25','c50','c75','c100'}
        cond = conds{1};
        
        % pick the right dataET struct
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
        
        %% Trial loop
        for trl = 1:numel(dataET.trialinfo)
            raw    = dataET.trial{trl};
            tVec   = dataET.time{trl};
            
            %—— BASELINE WINDOW ——%
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
            
            %—— ANALYSIS WINDOW ——%
            an_idx = tVec >= analysis_period(1) & tVec <= analysis_period(2);
            an_dat = raw(:,an_idx);
            valid  = an_dat(1,:)>=0 & an_dat(1,:)<=800 & an_dat(2,:)>=0 & an_dat(2,:)<=600;
            an_dat = an_dat(1:3, valid);
            an_dat(2,:) = 600 - an_dat(2,:);
            an_dat = remove_blinks(an_dat, win_size);
            
            % extract the cleaned gaze trace
            x = an_dat(1,:);  
            y = an_dat(2,:);
            
            % store raw gaze traces for this subject/condition/trial
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
            
            % compute analysis‐window metrics
            dx = x - 400 - nanmean(x-400);
            dy = y - 300 - nanmean(y-300);
            
            mean_eucdev = nanmean( sqrt(dx.^2 + dy.^2) );
            std_x       = nanstd(x);
            std_y       = nanstd(y);
            pupil       = mean(an_dat(3,:),'omitnan')/1000;
            [msrate, ~] = detect_microsaccades(fsample, [x; y], numel(x));
            
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
            baselineMSRate(end+1)    = baseline_msrate;
        end % trial loop
        
        %% 4) SUBJECT‐BY‐CONDITION AVERAGES
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
            
            c25_pct_gdev    = (c25_gdev    - c25_bl_gdev   ) / c25_bl_gdev   * 100;
            c25_pct_gSDx    = (c25_gSDx    - c25_bl_gSDx   ) / c25_bl_gSDx   * 100;
            c25_pct_gSDy    = (c25_gSDy    - c25_bl_gSDy   ) / c25_bl_gSDy   * 100;
            c25_pct_pups    = (c25_pups    - c25_bl_pups   ) / c25_bl_pups   * 100;
            c25_pct_msrate  = (c25_msrate  - c25_bl_msrate ) / c25_bl_msrate * 100;
            
            subj_data_gaze_trial_c25      = struct( ...
                'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'PctGazeDeviation',   (gazeDev./baselineGazeDev-1)*100, ...
                'GazeStdX',gazeSDx,         'BaselineGazeStdX',baselineGazeSDx,         'PctGazeStdX',         (gazeSDx./baselineGazeSDx-1)*100, ...
                'GazeStdY',gazeSDy,         'BaselineGazeStdY',baselineGazeSDy,         'PctGazeStdY',         (gazeSDy./baselineGazeSDy-1)*100, ...
                'PupilSize',pupilSize,      'BaselinePupilSize',baselinePupilSize,      'PctPupilSize',        (pupilSize./baselinePupilSize-1)*100, ...
                'MSRate',microsaccadeRate,  'BaselineMSRate',baselineMSRate,            'PctMSRate',           (microsaccadeRate./baselineMSRate-1)*100);
            
          case 'c50'
            % … exactly the same pattern for c50 …
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
            
            c50_pct_gdev    = (c50_gdev    - c50_bl_gdev   ) / c50_bl_gdev   * 100;
            c50_pct_gSDx    = (c50_gSDx    - c50_bl_gSDx   ) / c50_bl_gSDx   * 100;
            c50_pct_gSDy    = (c50_gSDy    - c50_bl_gSDy   ) / c50_bl_gSDy   * 100;
            c50_pct_pups    = (c50_pups    - c50_bl_pups   ) / c50_bl_pups   * 100;
            c50_pct_msrate  = (c50_msrate  - c50_bl_msrate ) / c50_bl_msrate * 100;
            
            subj_data_gaze_trial_c50 = struct( ...
                'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'PctGazeDeviation',(gazeDev./baselineGazeDev-1)*100, ...
                'GazeStdX',gazeSDx,         'BaselineGazeStdX',baselineGazeSDx,         'PctGazeStdX',        (gazeSDx./baselineGazeSDx-1)*100, ...
                'GazeStdY',gazeSDy,         'BaselineGazeStdY',baselineGazeSDy,         'PctGazeStdY',        (gazeSDy./baselineGazeSDy-1)*100, ...
                'PupilSize',pupilSize,      'BaselinePupilSize',baselinePupilSize,      'PctPupilSize',       (pupilSize./baselinePupilSize-1)*100, ...
                'MSRate',microsaccadeRate,  'BaselineMSRate',baselineMSRate,            'PctMSRate',          (microsaccadeRate./baselineMSRate-1)*100);
            
          case 'c75'
            % … same again …
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
            
            c75_pct_gdev    = (c75_gdev    - c75_bl_gdev   ) / c75_bl_gdev   * 100;
            c75_pct_gSDx    = (c75_gSDx    - c75_bl_gSDx   ) / c75_bl_gSDx   * 100;
            c75_pct_gSDy    = (c75_gSDy    - c75_bl_gSDy   ) / c75_bl_gSDy   * 100;
            c75_pct_pups    = (c75_pups    - c75_bl_pups   ) / c75_bl_pups   * 100;
            c75_pct_msrate  = (c75_msrate  - c75_bl_msrate ) / c75_bl_msrate * 100;
            
            subj_data_gaze_trial_c75 = struct( ...
                'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'PctGazeDeviation',(gazeDev./baselineGazeDev-1)*100, ...
                'GazeStdX',gazeSDx,         'BaselineGazeStdX',baselineGazeSDx,         'PctGazeStdX',        (gazeSDx./baselineGazeSDx-1)*100, ...
                'GazeStdY',gazeSDy,         'BaselineGazeStdY',baselineGazeSDy,         'PctGazeStdY',        (gazeSDy./baselineGazeSDy-1)*100, ...
                'PupilSize',pupilSize,      'BaselinePupilSize',baselinePupilSize,      'PctPupilSize',       (pupilSize./baselinePupilSize-1)*100, ...
                'MSRate',microsaccadeRate,  'BaselineMSRate',baselineMSRate,            'PctMSRate',          (microsaccadeRate./baselineMSRate-1)*100);
            
          case 'c100'
            % … and once more for c100 …
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
            
            c100_pct_gdev    = (c100_gdev    - c100_bl_gdev   ) / c100_bl_gdev   * 100;
            c100_pct_gSDx    = (c100_gSDx    - c100_bl_gSDx   ) / c100_bl_gSDx   * 100;
            c100_pct_gSDy    = (c100_gSDy    - c100_bl_gSDy   ) / c100_bl_gSDy   * 100;
            c100_pct_pups    = (c100_pups    - c100_bl_pups   ) / c100_bl_pups   * 100;
            c100_pct_msrate  = (c100_msrate  - c100_bl_msrate ) / c100_bl_msrate * 100;
            
            subj_data_gaze_trial_c100 = struct( ...
                'ID',subject_id,'Trial',trial_num,'Condition',condition, ...
                'GazeDeviation',gazeDev,   'BaselineGazeDeviation',baselineGazeDev,   'PctGazeDeviation',(gazeDev./baselineGazeDev-1)*100, ...
                'GazeStdX',gazeSDx,         'BaselineGazeStdX',baselineGazeSDx,         'PctGazeStdX',        (gazeSDx./baselineGazeSDx-1)*100, ...
                'GazeStdY',gazeSDy,         'BaselineGazeStdY',baselineGazeSDy,         'PctGazeStdY',        (gazeSDy./baselineGazeSDy-1)*100, ...
                'PupilSize',pupilSize,      'BaselinePupilSize',baselinePupilSize,      'PctPupilSize',       (pupilSize./baselinePupilSize-1)*100, ...
                'MSRate',microsaccadeRate,  'BaselineMSRate',baselineMSRate,            'PctMSRate',          (microsaccadeRate./baselineMSRate-1)*100);
        end
    end  % end conds loop
    
    %% 5) CREATE SUBJECT‐LEVEL STRUCTS ACROSS CONDITIONS
    savepath = fullfile(path, subjects{subj}, 'gaze');
    if ~exist(savepath,'dir'); mkdir(savepath); end
    
    % load the pre‐extracted GCP_preprocessing metrics
    load(fullfile(savepath,'gaze_metrics'),'c25_blinks','c50_blinks','c75_blinks','c100_blinks', ...
                                                'c25_fixations','c50_fixations','c75_fixations','c100_fixations', ...
                                                'c25_saccades','c50_saccades','c75_saccades','c100_saccades');
    
    subj_data_gaze = struct( ...
      'ID',        num2cell(subject_id(1:4))', ...
      'Condition', num2cell([1;2;3;4]), ...
      'GazeDeviation', num2cell([c25_gdev;  c50_gdev;  c75_gdev;  c100_gdev]), ...
      'GazeStdX',      num2cell([c25_gSDx;  c50_gSDx;  c75_gSDx;  c100_gSDx]), ...
      'GazeStdY',      num2cell([c25_gSDy;  c50_gSDy;  c75_gSDy;  c100_gSDy]), ...
      'PupilSize',     num2cell([c25_pups;  c50_pups;  c75_pups;  c100_pups]), ...
      'MSRate',        num2cell([c25_msrate;c50_msrate;c75_msrate;c100_msrate]), ...
      'Blinks',        num2cell([c25_blinks;c50_blinks;c75_blinks;c100_blinks]), ...
      'Fixations',     num2cell([c25_fixations;c50_fixations;c75_fixations;c100_fixations]), ...
      'Saccades',      num2cell([c25_saccades; c50_saccades; c75_saccades; c100_saccades]) );

    subj_data_gaze_baseline = struct( ...
      'ID',        num2cell(subject_id(1:4))', ...
      'Condition', num2cell([1;2;3;4]), ...
      'BaselineGazeDeviation', num2cell([c25_bl_gdev;  c50_bl_gdev;  c75_bl_gdev;  c100_bl_gdev]), ...
      'BaselineGazeStdX',      num2cell([c25_bl_gSDx;  c50_bl_gSDx;  c75_bl_gSDx;  c100_bl_gSDx]), ...
      'BaselineGazeStdY',      num2cell([c25_bl_gSDy;  c50_bl_gSDy;  c75_bl_gSDy;  c100_bl_gSDy]), ...
      'BaselinePupilSize',     num2cell([c25_bl_pups;  c50_bl_pups;  c75_bl_pups;  c100_bl_pups]), ...
      'BaselineMSRate',        num2cell([c25_bl_msrate;c50_bl_msrate;c75_bl_msrate;c100_bl_msrate]) );

    subj_data_gaze_pctchange = struct( ...
      'ID',        num2cell(subject_id(1:4))', ...
      'Condition', num2cell([1;2;3;4]), ...
      'PctGazeDeviation', num2cell([c25_pct_gdev;  c50_pct_gdev;  c75_pct_gdev;  c100_pct_gdev]), ...
      'PctGazeStdX',      num2cell([c25_pct_gSDx;  c50_pct_gSDx;  c75_pct_gSDx;  c100_pct_gSDx]), ...
      'PctGazeStdY',      num2cell([c25_pct_gSDy;  c50_pct_gSDy;  c75_pct_gSDy;  c100_pct_gSDy]), ...
      'PctPupilSize',     num2cell([c25_pct_pups;  c50_pct_pups;  c75_pct_pups;  c100_pct_pups]), ...
      'PctMSRate',        num2cell([c25_pct_msrate;c50_pct_msrate;c75_pct_msrate;c100_pct_msrate]) );

    %% 6) SAVE EVERYTHING
    save(fullfile(savepath,'gaze_matrix_trial'),  ...
         'subj_data_gaze_trial_c25',  'subj_data_gaze_trial_c50',  ...
         'subj_data_gaze_trial_c75',  'subj_data_gaze_trial_c100');
    save(fullfile(savepath,'gaze_matrix_subj'),   'subj_data_gaze');
    save(fullfile(savepath,'gaze_matrix_baseline'),'subj_data_gaze_baseline');
    save(fullfile(savepath,'gaze_matrix_pctchange'),'subj_data_gaze_pctchange');
    
    % optionally append to your across‐subjects raw struct
    gaze_data = [gaze_data; subj_data_gaze];
    gaze_data_bl = [gaze_data_bl; subj_data_gaze_pctchange];
    
    fprintf('Subject %d/%d done.\n', subj, numel(subjects));
end

%% Save master files
save(fullfile(path,'gaze_raw'),    'gaze_x_c25','gaze_y_c25','gaze_x_c50','gaze_y_c50','gaze_x_c75','gaze_y_c75','gaze_x_c100','gaze_y_c100');
save(fullfile(path,'gaze_matrix'), 'gaze_data');
save(fullfile(path,'gaze_matrix_bl'), 'gaze_data_bl');
