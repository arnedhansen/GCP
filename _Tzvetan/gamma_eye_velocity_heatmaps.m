%% gamma do eye velocity and heatmaps
clear all
close all
addpath('/Users/tpopov/Documents/matlabtools/eeglab2022.1')
eeglab
close all hidden
%%
subj = {'601';'602';'603';'604';'605';'606';'607';'608';'609';'610'};

for s = 1:length(subj)
    
    path = strcat('/Users/tpopov/Documents/DATA4FT/gamma_arne/',subj{s});
    cd (path)
    
    
    
    
    fileList = dir(fullfile(path, [subj{s}, '*.mat']));
    for b=1:numel(fileList)
        load(fileList(b).name)
        data25 = pop_epoch( EEG, {'61'}, [-2 3]);
        data50 = pop_epoch( EEG, {'62'}, [-2 3]);
        data75  = pop_epoch( EEG, {'63'}, [-2 3]);
        data100  = pop_epoch( EEG, {'64'}, [-2 3]);
        %%
        dat25{b}=eeglab2fieldtrip(data25,'raw');
        dat50{b}=eeglab2fieldtrip(data50,'raw');
        dat75{b}=eeglab2fieldtrip(data75,'raw');
        dat100{b}=eeglab2fieldtrip(data100,'raw');
    end
    %%
    data25 = ft_appenddata([],dat25{:});
    data50 = ft_appenddata([],dat50{:});
    data75 = ft_appenddata([],dat75{:});
    data100 = ft_appenddata([],dat100{:});
    %% select et data
    cfg = [];
    cfg.channel = {'L-GAZE-X';'L-GAZE-Y';'L-AREA'};
    data25et = ft_selectdata(cfg,data25);
    data50et = ft_selectdata(cfg,data50);
    data75et = ft_selectdata(cfg,data75);
    data100et = ft_selectdata(cfg,data100);
    %%
    cfg = [];
    cfg.channel = {'EEG' ; '-L-GAZE-X';'-L-GAZE-Y';'-L-AREA'};
    data25 = ft_selectdata(cfg,data25);
    data50 = ft_selectdata(cfg,data50);
    data75 = ft_selectdata(cfg,data75);
    data100 = ft_selectdata(cfg,data100);
    %%
    cfg = [];
    cfg.reref = 'yes';
    cfg.refchannel = 'all';
    data25 = ft_preprocessing(cfg,data25);
    data50 = ft_preprocessing(cfg,data50);
    data75 = ft_preprocessing(cfg,data75);
    data100 = ft_preprocessing(cfg,data100);
    %%
    save data25 data25
    save data50 data50
    save data75 data75
    save data100 data100
    
    
    save data25et data25et
    save data50et data50et
    save data75et data75et
    save data100et data100et
    %% compute velocity
    et=data25et;
    et_sm=et;
    et_sm=rmfield(et_sm,'trial');
    et_sm=rmfield(et_sm,'time');
    
    sampling_rate= et.fsample;
    window_size = 25;% corresponds to 100 ms window
    for trl=1:length(et.trial)
        
        
        % Compute the moving average
        smoothed_data = movmean(et.trial{trl}(1,:), window_size);
        velocity = diff(smoothed_data) * sampling_rate;
        smoothed_velocity = movmean(abs(velocity), window_size);
        
        et_sm.trial{trl}(1,:) = abs(velocity);
        et_sm.time{trl}=et.time{trl}(1:end-1);
        smoothed_data = movmean(et.trial{trl}(2,:), window_size);
        velocity = diff(smoothed_data) * sampling_rate;
        smoothed_velocity = movmean(abs(velocity), window_size);
        et_sm.trial{trl}(2,:) = abs(velocity);
        et_sm.trial{trl}(3,:)=et.trial{trl}(3,1:end-1);
    end
    et_sm25=et_sm;
    %% velocity 50
    et=data50et;
    et_sm=et;
    et_sm=rmfield(et_sm,'trial');
    et_sm=rmfield(et_sm,'time');
    
    sampling_rate= et.fsample;
    window_size = 25;% corresponds to 100 ms window
    for trl=1:length(et.trial)
        
        
        % Compute the moving average
        smoothed_data = movmean(et.trial{trl}(1,:), window_size);
        velocity = diff(smoothed_data) * sampling_rate;
        smoothed_velocity = movmean(abs(velocity), window_size);
        
        et_sm.trial{trl}(1,:) = abs(velocity);
        et_sm.time{trl}=et.time{trl}(1:end-1);
        smoothed_data = movmean(et.trial{trl}(2,:), window_size);
        velocity = diff(smoothed_data) * sampling_rate;
        smoothed_velocity = movmean(abs(velocity), window_size);
        et_sm.trial{trl}(2,:) = abs(velocity);
        et_sm.trial{trl}(3,:)=et.trial{trl}(3,1:end-1);
    end
    et_sm50=et_sm;
    %% velocity 75
    et=data75et;
    et_sm=et;
    et_sm=rmfield(et_sm,'trial');
    et_sm=rmfield(et_sm,'time');
    
    sampling_rate= et.fsample;
    window_size = 25;% corresponds to 100 ms window
    for trl=1:length(et.trial)
        
        % Compute the moving average
        smoothed_data = movmean(et.trial{trl}(1,:), window_size);
        velocity = diff(smoothed_data) * sampling_rate;
        smoothed_velocity = movmean(abs(velocity), window_size);
        
        et_sm.trial{trl}(1,:) = abs(velocity);
        et_sm.time{trl}=et.time{trl}(1:end-1);
        smoothed_data = movmean(et.trial{trl}(2,:), window_size);
        velocity = diff(smoothed_data) * sampling_rate;
        smoothed_velocity = movmean(abs(velocity), window_size);
        et_sm.trial{trl}(2,:) = abs(velocity);
        et_sm.trial{trl}(3,:)=et.trial{trl}(3,1:end-1);
    end
    et_sm75=et_sm;
    %% velocity 50
    et=data100et;
    et_sm=et;
    et_sm=rmfield(et_sm,'trial');
    et_sm=rmfield(et_sm,'time');
    
    sampling_rate= et.fsample;
    window_size = 25;% corresponds to 100 ms window
    for trl=1:length(et.trial)
        
        
        % Compute the moving average
        smoothed_data = movmean(et.trial{trl}(1,:), window_size);
        velocity = diff(smoothed_data) * sampling_rate;
        smoothed_velocity = movmean(abs(velocity), window_size);
        
        et_sm.trial{trl}(1,:) = abs(velocity);
        et_sm.time{trl}=et.time{trl}(1:end-1);
        smoothed_data = movmean(et.trial{trl}(2,:), window_size);
        velocity = diff(smoothed_data) * sampling_rate;
        smoothed_velocity = movmean(abs(velocity), window_size);
        et_sm.trial{trl}(2,:) = abs(velocity);
        et_sm.trial{trl}(3,:)=et.trial{trl}(3,1:end-1);
    end
    et_sm100=et_sm;
    %% apply baseline and average over trials
    cfg = [];
    cfg.baseline = [-2 0];
    et_sm25 = ft_timelockbaseline(cfg,et_sm25);
    et_sm50 = ft_timelockbaseline(cfg,et_sm50);
    et_sm75 = ft_timelockbaseline(cfg,et_sm75);
    et_sm100 = ft_timelockbaseline(cfg,et_sm100);
    cfg = [];
    cfg.latency = [-1.5 2];
    cfg.keeptrials = 'no';
    tlk_et_sm25=ft_timelockanalysis(cfg,et_sm25);
    tlk_et_sm50=ft_timelockanalysis(cfg,et_sm50);
    tlk_et_sm75=ft_timelockanalysis(cfg,et_sm75);
    tlk_et_sm100=ft_timelockanalysis(cfg,et_sm100);
    %%
    save et_sm25 et_sm25
    save et_sm50 et_sm50
    save et_sm75 et_sm75
    save et_sm100 et_sm100
    
end
%%
close all
clear all
subj = {'601';'602';'603';'604';'605';'606';'607';'608';'609';'610'};

for s = 1:length(subj)
    
    path = strcat('/Users/tpopov/Documents/DATA4FT/gamma_arne/',subj{s});
    cd (path)
    load et_sm25
    load et_sm50
    load et_sm75
    load et_sm100
%            
        cfg = [];
    cfg.latency = [-1.5 2];
    cfg.keeptrials = 'no';
    tlk_et_sm25=ft_timelockanalysis(cfg,et_sm25);
    tlk_et_sm50=ft_timelockanalysis(cfg,et_sm50);
    tlk_et_sm75=ft_timelockanalysis(cfg,et_sm75);
    tlk_et_sm100=ft_timelockanalysis(cfg,et_sm100);
    
     cfg = [];
    cfg.baseline = [-2 0];
    tlk_et_sm25 = ft_timelockbaseline(cfg,tlk_et_sm25);
    tlk_et_sm50 = ft_timelockbaseline(cfg,tlk_et_sm50);
    tlk_et_sm75 = ft_timelockbaseline(cfg,tlk_et_sm75);
    tlk_et_sm100 = ft_timelockbaseline(cfg,tlk_et_sm100);
    %%
    alltlk25et{s}=tlk_et_sm25;
    alltlk50et{s}=tlk_et_sm50;
    alltlk75et{s}=tlk_et_sm75;
    alltlk100et{s}=tlk_et_sm100;
    
end
%%

cfg = [];
cfg.keepindividual='yes';
ga25et = ft_timelockgrandaverage(cfg,alltlk25et{:});
ga50et = ft_timelockgrandaverage(cfg,alltlk50et{:});
ga75et = ft_timelockgrandaverage(cfg,alltlk75et{:});
ga100et = ft_timelockgrandaverage(cfg,alltlk100et{:});
%% Plot eye velocity for 25/50/75/100% contrast (mean ± SEM)
close all
channels = {'L-GAZE-X','L-GAZE-Y','L-AREA'};
channeltitles = {'horizontal velocity','vertical velocity','pupil'};
ylabs    = {'eye velocity X [a.u.]','eye velocity Y [a.u.]','pupil area [a.u.]'};

labels  = {'25% contrast','50% contrast','75% contrast','100% contrast'};
lineCs  = [1 0 0; 1 .5 0; .6 0.2 .8; 0 0 1];   % red, orange, purple, blue

% assume equal n_subj across contrasts
n_subj = size(ga25et.individual, 1);

figure('Color','w');

for c = 1:numel(channels)
    subplot(2,2,c); hold on;

    % Select the channel for each contrast
    cfg = [];
    cfg.figure = 'gcf';
    cfg.channel     = channels{c};
    cfg.avgoverchan = 'yes';
    cfg.latency     = [-.2 2];

    et25  = ft_selectdata(cfg, ga25et);
    et50  = ft_selectdata(cfg, ga50et);
    et75  = ft_selectdata(cfg, ga75et);
    et100 = ft_selectdata(cfg, ga100et);

    ets = {et25, et50, et75, et100};

    hl = gobjects(1, numel(ets));
    for k = 1:numel(ets)
        x  = ets{k}.time(:);
        yi = squeeze(ets{k}.individual);        % [subj x time]
        if size(yi,1) ~= n_subj
            warning('n_subj mismatch for %s', labels{k});
        end

        y   = mean(yi, 1)';                      % mean over subjects
        e   = std(yi, [], 1)' ./ sqrt(n_subj);   % SEM
        low = y - e;
        high = y + e;

        % Light face color for band
        faceC = 0.8*lineCs(k,:) + 0.2*[1 1 1];

        hp = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], lineCs(k,:));
        set(hp, 'facecolor', faceC, 'edgecolor', 'none', 'facealpha', 0.6);

        hl(k) = plot(x, y, 'LineWidth', 2, 'Color', lineCs(k,:));
    end

    set(gca, 'FontSize', 18);
    xlabel('Time [sec]');
    ylabel(ylabs{c});
    title(strrep(channeltitles{c}, '_', '\_'));
    xlim([-.2 2]); box on; grid on;
    xticks([-.2 0 .2 .5 1 2]);

    if c == 1
        lgd = legend(hl, labels, 'Location', 'northeast', 'FontSize', 12);
        set(lgd, 'Color', 'none');
    end
end

% Optional overall title
% sgtitle('Gaze (X/Y) and Pupil AREA across contrasts');

%% do Ftest on eye velocity
cfg                  = [];
cfg.method           = 'montecarlo'; 
cfg.statistic        = 'ft_statfun_depsamplesFunivariate';                                  
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;       
cfg.clusterstatistic = 'maxsum';   
cfg.neighbours       = []; 
cfg.tail             = 1;        
cfg.clustertail      = cfg.tail;  
cfg.alpha            = 0.05;      
cfg.numrandomization = 1000;        % number of draws from the permutation distribution

n_25 = size(ga25et.individual, 1);
n_50  = size(ga50et.individual, 1);
n_75 = size(ga75et.individual, 1);
n_100 = size(ga100et.individual, 1);

clear design
design = zeros(2,4*n_100);
cfg.design(1,:)           = [ones(1,n_25), ones(1,n_50)*2,ones(1,n_75)*3,ones(1,n_100)*4]; 
cfg.design(2,:)           = [1:n_25,1:n_50, 1:n_75,1:n_100]; 
cfg.ivar                  = 1; 
cfg.uvar                  = 2;
[statF_eyevelocity] = ft_timelockstatistics(cfg, ga25et,ga50et, ga75et,ga100et);
%% insert in GA figure

ax4 = subplot(2,2,4); pos4 = get(ax4,'Position'); delete(ax4);
p4 = uipanel('Parent', gcf, 'Position', pos4, ...
    'BorderType','none', 'BackgroundColor','w');

tl = tiledlayout(p4, 3, 1, 'TileSpacing','compact', 'Padding','compact');

% --- base cfg: draw into current axes, NO internal mask fill
cfg = [];
cfg.parameter = 'stat';
cfg.figure    = 'gca';   % draw in current axes
cfg.channel   = [];      % set per loop
cfg.linecolor = 'k';
% IMPORTANT: don't let FT draw the mask (prevents seams)
cfg = rmfield(cfg, intersect(fieldnames(cfg), {'maskparameter','maskstyle'}));

statChans  = {'L-GAZE-X','L-GAZE-Y','L-AREA'};
statTitles = {'F-test: horizontal velocity','F-test: vertical velocity','F-test: pupil'};
for i = 1:3
    ax = nexttile(tl,i); axes(ax); cla; hold on; box on;
    set(ax, 'Color','w', 'FontSize',12);

    % plot the stat trace
    cfg.channel = statChans{i};
    ft_singleplotER(cfg, statF_eyevelocity);
    title(statTitles{i}, 'Interpreter','none');
    xlabel('Time [sec]'); ylabel('F-value');   % change label if needed
    xlim([-.2 2]); xticks([-.2 0 .2 .5 1 2]);
    grid on
    % ---- add clean mask shading (no vertical edges) ----
    % find the channel index in the stat struct
    chIdx = find(strcmp(statF_eyevelocity.label, statChans{i}), 1);
    if ~isempty(chIdx) && isfield(statF_eyevelocity,'mask')
        t = statF_eyevelocity.time(:)';                 % [1 x time]
        m = logical(statF_eyevelocity.mask(chIdx, :));  % [1 x time]

        % merge contiguous "true" runs
        dm = diff([false, m, false]);
        starts = find(dm == 1);
        stops  = find(dm == -1) - 1;

        yl = ylim(ax);
        for r = 1:numel(starts)
            t1 = t(starts(r));
            t2 = t(stops(r));
            % single seamless rectangle per run
            patch('XData',[t1 t2 t2 t1], 'YData',[yl(1) yl(1) yl(2) yl(2)], ...
                  'FaceColor',[0 0 0], 'FaceAlpha',0.12, 'EdgeColor','none', ...
                  'Parent', ax);
        end
        uistack(findobj(ax,'Type','line'), 'top');   % keep the line above patches
    end
end
%% Make gaze heatmaps and grand-average across subjects
% clear; 
% close all;


subjects = {'601','602','603','604','605','606','607','608','609','610'};
base_dir = '/Users/tpopov/Documents/DATA4FT/gamma_arne/';
addpath('/Volumes/Homestore/OCC/arne/funcs'); 

% Gaze heatmap parameters
sampling_rate = 500;
threshold = 20;
pad_ms = 150;
num_bins = 1000;
smooth_val = 5;
x_grid = linspace(0, 800, num_bins);
y_grid = linspace(0, 600, num_bins);

% initiate subjects
nS = numel(subjects);

% raw
allgazebase25      = cell(1,nS); allgazetaskearly25      = cell(1,nS); allgazetasklate25      = cell(1,nS);
allgazebase50      = cell(1,nS); allgazetaskearly50      = cell(1,nS); allgazetasklate50      = cell(1,nS);
allgazebase75      = cell(1,nS); allgazetaskearly75      = cell(1,nS); allgazetasklate75      = cell(1,nS);
allgazebase100     = cell(1,nS); allgazetaskearly100     = cell(1,nS); allgazetasklate100     = cell(1,nS);

% normalized
allgazebase25_norm  = cell(1,nS); allgazetaskearly25_norm  = cell(1,nS); allgazetasklate25_norm  = cell(1,nS);
allgazebase50_norm  = cell(1,nS); allgazetaskearly50_norm  = cell(1,nS); allgazetasklate50_norm  = cell(1,nS);
allgazebase75_norm  = cell(1,nS); allgazetaskearly75_norm  = cell(1,nS); allgazetasklate75_norm  = cell(1,nS);
allgazebase100_norm = cell(1,nS); allgazetaskearly100_norm = cell(1,nS); allgazetasklate100_norm = cell(1,nS);


results = struct([]);

%% Loop over subjects
for s = 1:nS
    subj = subjects{s};
    subj_dir = fullfile(base_dir, subj);
    cd(subj_dir);

    load data25et
    load data50et
    load data75et
    load data100et

    labels = {'25','50','75','100'};
    vars   = {'data25et','data50et','data75et','data100et'};

    for k = 1:numel(vars)
        dsName = labels{k};
        dataET = eval(vars{k});  

        
        cfg = []; 
        cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
        dataet = ft_selectdata(cfg, dataET);
        nTrials = numel(dataet.trial);

        % --- Blink correction & trial QC ---
        dataetnan = dataet;
        valid_trials = true(nTrials,1);

        for i = 1:nTrials
            x = dataet.trial{i}(1,:);
            y = dataet.trial{i}(2,:);
            t = dataet.time{i};

            [x_nan, y_nan, ~, ~, ~, is_valid] = ...
                removeAndInterpolateBlinks_checktrials(x, y, t, sampling_rate, threshold, pad_ms);

            if ~is_valid
                valid_trials(i) = false;
                continue;
            end
            dataetnan.trial{i}(1,:) = x_nan;
            dataetnan.trial{i}(2,:) = y_nan;
        end

        % Keep only valid trials (IMPORTANT so subjects contribute comparable data)
        cfg = [];
        cfg.trials = find(valid_trials);
        dataetnan = ft_selectdata(cfg, dataetnan);

        % --- Time windows ---
        cfg = []; 
        cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
        cfg.latency = [0 1];        
        dat_early = ft_selectdata(cfg, dataetnan);
        cfg.latency = [1 2];        
        dat_late  = ft_selectdata(cfg, dataetnan);
        cfg.latency = [-1.25 -.25];  
        dat_base  = ft_selectdata(cfg, dataetnan);

        % compute heamap (raw + normalized)
        [freq_early, freq_early_norm] = computeGazeHeatmap(dat_early, x_grid, y_grid, sampling_rate, smooth_val);
        [freq_late,  freq_late_norm ] = computeGazeHeatmap(dat_late,  x_grid, y_grid, sampling_rate, smooth_val);
        [freq_base,  freq_base_norm ] = computeGazeHeatmap(dat_base,  x_grid, y_grid, sampling_rate, smooth_val);

        % subject-level cell arrays for grand average ---
        switch dsName
            case '25'
                allgazebase25{s}       = freq_base;
                allgazetaskearly25{s}  = freq_early;
                allgazetasklate25{s}   = freq_late;
                allgazebase25_norm{s}      = freq_base_norm;
                allgazetaskearly25_norm{s} = freq_early_norm;
                allgazetasklate25_norm{s}  = freq_late_norm;

            case '50'
                allgazebase50{s}       = freq_base;
                allgazetaskearly50{s}  = freq_early;
                allgazetasklate50{s}   = freq_late;
                allgazebase50_norm{s}      = freq_base_norm;
                allgazetaskearly50_norm{s} = freq_early_norm;
                allgazetasklate50_norm{s}  = freq_late_norm;

            case '75'
                allgazebase75{s}       = freq_base;
                allgazetaskearly75{s}  = freq_early;
                allgazetasklate75{s}   = freq_late;
                allgazebase75_norm{s}      = freq_base_norm;
                allgazetaskearly75_norm{s} = freq_early_norm;
                allgazetasklate75_norm{s}  = freq_late_norm;

            case '100'
                allgazebase100{s}       = freq_base;
                allgazetaskearly100{s}  = freq_early;
                allgazetasklate100{s}   = freq_late;
                allgazebase100_norm{s}      = freq_base_norm;
                allgazetaskearly100_norm{s} = freq_early_norm;
                allgazetasklate100_norm{s}  = freq_late_norm;
        end

        % also keep per-subject in a struct
        results(s).subject = subj;
        results(s).(sprintf('valid_%s', dsName)) = find(valid_trials);
        results(s).(sprintf('n_valid_%s', dsName)) = numel(find(valid_trials));
        results(s).(sprintf('n_total_%s', dsName)) = nTrials;
    end
end

%% compute grand average raw
gagaze_b25 = ft_freqgrandaverage([], allgazebase25{:});
gagaze_e25 = ft_freqgrandaverage([], allgazetaskearly25{:});
gagaze_l25 = ft_freqgrandaverage([], allgazetasklate25{:});

gagaze_b50 = ft_freqgrandaverage([], allgazebase50{:});
gagaze_e50 = ft_freqgrandaverage([], allgazetaskearly50{:});
gagaze_l50 = ft_freqgrandaverage([], allgazetasklate50{:});

gagaze_b75 = ft_freqgrandaverage([], allgazebase75{:});
gagaze_e75 = ft_freqgrandaverage([], allgazetaskearly75{:});
gagaze_l75 = ft_freqgrandaverage([], allgazetasklate75{:});

gagaze_b100 = ft_freqgrandaverage([], allgazebase100{:});
gagaze_e100 = ft_freqgrandaverage([], allgazetaskearly100{:});
gagaze_l100 = ft_freqgrandaverage([], allgazetasklate100{:});

%% compute grand average normalized
gagaze_b25_norm = ft_freqgrandaverage([], allgazebase25_norm{:});
gagaze_e25_norm = ft_freqgrandaverage([], allgazetaskearly25_norm{:});
gagaze_l25_norm = ft_freqgrandaverage([], allgazetasklate25_norm{:});

gagaze_b50_norm = ft_freqgrandaverage([], allgazebase50_norm{:});
gagaze_e50_norm = ft_freqgrandaverage([], allgazetaskearly50_norm{:});
gagaze_l50_norm = ft_freqgrandaverage([], allgazetasklate50_norm{:});

gagaze_b75_norm = ft_freqgrandaverage([], allgazebase75_norm{:});
gagaze_e75_norm = ft_freqgrandaverage([], allgazetaskearly75_norm{:});
gagaze_l75_norm = ft_freqgrandaverage([], allgazetasklate75_norm{:});

gagaze_b100_norm = ft_freqgrandaverage([], allgazebase100_norm{:});
gagaze_e100_norm = ft_freqgrandaverage([], allgazetaskearly100_norm{:});
gagaze_l100_norm = ft_freqgrandaverage([], allgazetasklate100_norm{:});
%% make difference maps: (baseline - early) and (baseline - late)
contrasts = {'25','50','75','100'};
diffmaps = struct();   

for ic = 1:numel(contrasts)
    L   = contrasts{ic};         % '25','50','75','100'
    fld = ['c' L];               % valid struct field: 'c25','c50','c75','c100'

    base  = eval(['allgazebase' L]);         % cell 1×N
    early = eval(['allgazetaskearly' L]);    % cell 1×N
    late  = eval(['allgazetasklate' L]);     % cell 1×N
    nS = numel(base);

    % Preallocate outputs 
    eval(sprintf('cont%searly = cell(1,%d);', L, nS));
    eval(sprintf('cont%slate  = cell(1,%d);', L, nS));

    % preallocate inside diffmaps
    diffmaps.(fld).early = cell(1, nS);
    diffmaps.(fld).late  = cell(1, nS);

    for s = 1:nS
        A = base{s};  B = early{s};  C = late{s};
        if ~isequal(size(A.powspctrm), size(B.powspctrm)) || ~isequal(size(A.powspctrm), size(C.powspctrm))
            error('powspctrm size mismatch for contrast %s, subject %d', L, s);
        end

        % baseline - early
        D = A;  D.powspctrm = A.powspctrm - B.powspctrm;
        % baseline - late
        E = A;  E.powspctrm = A.powspctrm - C.powspctrm;

        % Save 
        eval(sprintf('cont%searly{%d} = D;', L, s));
        eval(sprintf('cont%slate{%d}  = E;', L, s));

        % Save into struct with valid field name
        diffmaps.(fld).early{s} = D;
        diffmaps.(fld).late{s}  = E;
    end
end

%% Grand-averages of difference maps 
ga_cont25_early  = ft_freqgrandaverage([], diffmaps.c25.early{:});
ga_cont25_late   = ft_freqgrandaverage([], diffmaps.c25.late{:});
ga_cont50_early  = ft_freqgrandaverage([], diffmaps.c50.early{:});
ga_cont50_late   = ft_freqgrandaverage([], diffmaps.c50.late{:});
ga_cont75_early  = ft_freqgrandaverage([], diffmaps.c75.early{:});
ga_cont75_late   = ft_freqgrandaverage([], diffmaps.c75.late{:});
ga_cont100_early = ft_freqgrandaverage([], diffmaps.c100.early{:});
ga_cont100_late  = ft_freqgrandaverage([], diffmaps.c100.late{:});
%% plot
figure;
cfg = [];
cfg.zlim = [-.01 .01];
% cfg.ylim = [0 1];
cfg.figure= 'gcf';
subplot(2,2,1);
ft_singleplotTFR(cfg,ga_cont25_early);
title('25%');
subplot(2,2,2);
ft_singleplotTFR(cfg,ga_cont50_early);
title('50%');
subplot(2,2,3);
ft_singleplotTFR(cfg,ga_cont75_early);
title('75%');
subplot(2,2,4);
ft_singleplotTFR(cfg,ga_cont100_early);
title('100%');
%% plot late
figure;
cfg = [];
cfg.zlim = [-.01 .01];
% cfg.ylim = [0 1];
cfg.figure= 'gcf';
subplot(2,2,1);
ft_singleplotTFR(cfg,ga_cont25_late);
title('25%');
subplot(2,2,2);
ft_singleplotTFR(cfg,ga_cont50_late);
title('50%');
subplot(2,2,3);
ft_singleplotTFR(cfg,ga_cont75_late);
title('75%');
subplot(2,2,4);
ft_singleplotTFR(cfg,ga_cont100_late);
title('100%');
%% F test
cfg                  = [];
cfg.method           = 'montecarlo'; 
cfg.statistic        = 'ft_statfun_depsamplesFunivariate';                                  
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;      
cfg.clusterstatistic = 'maxsum';  
cfg.neighbours       = []; 
cfg.tail             = 1;          
cfg.clustertail      = cfg.tail;  
cfg.alpha            = 0.05;      
cfg.numrandomization = 1000;        

n_25  = numel(allgazetasklate25);
n_50  = numel(allgazetasklate50);
n_75 = numel(allgazetasklate75);
n_100 = numel(allgazetasklate100);

clear design
design = zeros(2,4*n_100);
cfg.design(1,:)           = [ones(1,n_25), ones(1,n_50)*2,ones(1,n_75)*3,ones(1,n_100)*4]; 
cfg.design(2,:)           = [1:n_25,1:n_50, 1:n_75,1:n_100]; 
cfg.ivar                  = 1; 
cfg.uvar                  = 2;
[statFlate] = ft_freqstatistics(cfg, allgazetasklate25{:},allgazetasklate50{:},allgazetasklate75{:},allgazetasklate100{:});
statFlate.stat(statFlate.mask==0)=0;
[statFearly] = ft_freqstatistics(cfg, allgazetaskearly25{:},allgazetaskearly50{:},allgazetaskearly75{:},allgazetaskearly100{:});
statFearly.stat(statFearly.mask==0)=0;
% do also for normalized data, that is % of time spend on location
% expressed in color
[statFlate_norm] = ft_freqstatistics(cfg, allgazetasklate25_norm{:},allgazetasklate50_norm{:},allgazetasklate75_norm{:},allgazetasklate100_norm{:});
statFlate_norm.stat(statFlate_norm.mask==0)=0;
[statFearly_norm] = ft_freqstatistics(cfg, allgazetaskearly25_norm{:},allgazetaskearly50_norm{:},allgazetaskearly75_norm{:},allgazetaskearly100_norm{:});
statFearly_norm.stat(statFearly_norm.mask==0)=0;

%%
close all
cfg = [];
cfg.parameter = 'stat';
cfg.figure = 'gcf';
cfg.maskparameter ='mask';
cfg.maskstyle = 'outline';
cfg.zlim = [-8 8];
figure; 
subplot(2,2,1);ft_singleplotTFR(cfg,statFearly);
subplot(2,2,2);ft_singleplotTFR(cfg,statFlate);
%% plot bar graph of single values
vals_late = mask_means_per_subject(statFlate, ...
    {allgazetasklate25_norm, allgazetasklate50_norm, allgazetasklate75_norm, allgazetasklate100_norm});

% group summaries
m_late = cellfun(@(v) mean(v,'omitnan'), vals_late);
e_late = cellfun(@(v) std(v,[],'omitnan') ./ sqrt(numel(v)), vals_late);

labels = {'25%','50%','75%','100%'};
lineCs = [1 0 0; 1 .5 0; .6 0.2 .8; 0 0 1];
faceCs = 0.8*lineCs + 0.2;                 

figure('Color','w');
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');


ax1 = nexttile(tl,1); axes(ax1); hold on; box on;
cfg = [];
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle     = 'outline';     
cfg.colormap = 'YlOrRd';
cfg.figure        = 'gca';
ft_singleplotTFR(cfg, statFlate);
title('F-test (1-2 sec post stim)','Fontsize',20);

xlabel('x pos [px]'); ylabel('y pos [px]');
box on
grid on
% add colorbar
caxis(ax1,[0 8]);
cb = colorbar(ax1,'Location','eastoutside');
cb.Label.String = 'F-value';
cb.Label.FontSize = 20;
cb.TickDirection = 'out';
% cb.Position = cbpos;
set(gcf,'Units','pixels','Position',[2223 275 482 739]);

ax2 = nexttile(tl,2); axes(ax2); hold on; box on;

% Bars (group means)
for k = 1:4
    bar(k, m_late(k), 'FaceColor', faceCs(k,:), ...
        'EdgeColor', lineCs(k,:), 'LineWidth', 1.2);
end

% Error bars
errorbar(1:4, m_late, e_late, 'k', 'LineStyle','none', 'LineWidth',1.2);

% Stack to [nSubj x 4]; allow missing with NaN
nSubj = max(cellfun(@numel, vals_late));
subMat = nan(nSubj,4);
for k = 1:4
    v = vals_late{k};
    subMat(1:numel(v),k) = v(:);
end

% Subject-constant jitter so lines pass through dots
jitWidth = 0.18;                     % horizontal spread of dots
sjit = (rand(nSubj,1)-0.5)*2*jitWidth;

% Draw faint connecting lines first (behind bars/dots)
for s = 1:nSubj
    y = subMat(s,:);
    if all(isnan(y)), continue; end
    x = (1:4) + sjit(s);
    plot(x, y, '-', 'Color', [0 0 0 0.3], 'LineWidth', 1);
end
title(ax2, { ...
    'Higher contrast → more distributed gaze', ...
    'Lower contrast → more focused gaze' }, ...
    'FontSize',18,'FontWeight','normal');


% Plot the dots on top, colored by condition
for k = 1:4
    y = subMat(:,k);
    x = k + sjit;
    scatter(x(~isnan(y)), y(~isnan(y)), 28, ...
        'MarkerFaceColor', lineCs(k,:), ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
end

set(ax2,'XTick',1:4,'XTickLabel',{'25%','50%','75%','100%'},'FontSize',20);
ylabel('mean within mask [a.u.]'); xlim([0.5 4.5]); grid on;

% Optional: re-outline bars to keep them crisp on top
for k = 1:4
    rectangle('Position',[k-0.4, 0, 0.8, m_late(k)], ...
        'EdgeColor', lineCs(k,:), 'LineWidth', 1.2, 'FaceColor', 'none');
end

%% do post hoc ttest
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.latency =[300 500];
% cfg.frequency =[200 400];
% cfg.statistic        = 'ft_statfun_diff';
% cfg.clusterthreshold ='nonparametric_common';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

cfg.neighbours=[];
clear design
subj = numel(subjects);
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;
%%
[stat25early] = ft_freqstatistics(cfg, allgazetaskearly25{:},allgazebase25{:});
cohensd=((stat25early.stat)./sqrt(numel(subjects)));
stat25early.stat=cohensd;
stat25early.stat(stat25early.mask==0)=0;% set everything not relevant to zero

[stat50early] = ft_freqstatistics(cfg, allgazetaskearly50{:},allgazebase50{:});
cohensd=((stat50early.stat)./sqrt(numel(subjects)));
stat50early.stat=cohensd;
stat50early.stat(stat50early.mask==0)=0;% set everything not relevant to zero


[stat75early] = ft_freqstatistics(cfg, allgazetaskearly75{:},allgazebase75{:});
cohensd=((stat75early.stat)./sqrt(numel(subjects)));
stat75early.stat=cohensd;
stat75early.stat(stat75early.mask==0)=0;% set everything not relevant to zero

[stat100early] = ft_freqstatistics(cfg,allgazetaskearly100{:},allgazebase100{:});
cohensd=((stat100early.stat)./sqrt(numel(subjects)));
stat100early.stat=cohensd;
stat100early.stat(stat100early.mask==0)=0;
%
[stat25late] = ft_freqstatistics(cfg, allgazetasklate25{:},allgazebase25{:});
cohensd=((stat25late.stat)./sqrt(numel(subjects)));
stat25late.stat=cohensd;
stat25late.stat(stat25late.mask==0)=0;% set everything not relevant to zero

[stat50late] = ft_freqstatistics(cfg, allgazetasklate50{:},allgazebase50{:});
cohensd=((stat50late.stat)./sqrt(numel(subjects)));
stat50late.stat=cohensd;
stat50late.stat(stat50late.mask==0)=0;% set everything not relevant to zero


[stat75late] = ft_freqstatistics(cfg, allgazetasklate75{:},allgazebase75{:});
cohensd=((stat75late.stat)./sqrt(numel(subjects)));
stat75late.stat=cohensd;
stat75late.stat(stat75late.mask==0)=0;% set everything not relevant to zero

[stat100late] = ft_freqstatistics(cfg,allgazetasklate100{:},allgazebase100{:});
cohensd=((stat100late.stat)./sqrt(numel(subjects)));
stat100late.stat=cohensd;
stat100late.stat(stat100late.mask==0)=0;
%%
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';
cfg.zlim = [-.5 .5];
% cfg.xlim =[300 500];
% cfg.ylim =[200 400];
cfg.figure = 'gcf';
figure;
subplot(3,2,1);
ft_singleplotTFR(cfg,stat25early);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');

grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.5 0 .5];
title(c,'Effect size \it d')
title('25% early')

subplot(3,2,2);
ft_singleplotTFR(cfg,stat25late);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');

grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.5 0 .5];
title(c,'Effect size \it d')
title('25% late')

subplot(3,2,3);
ft_singleplotTFR(cfg,stat50early);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');

grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.5 0 .5];
title(c,'Effect size \it d')
title('50% early')

subplot(3,2,4);
ft_singleplotTFR(cfg,stat50late);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');

grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.5 0 .5];
title(c,'Effect size \it d')
title('50% late')

subplot(3,2,5);
ft_singleplotTFR(cfg,stat100early);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');

grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.5 0 .5];
title(c,'Effect size \it d')
title('100% early')

subplot(3,2,6);
ft_singleplotTFR(cfg,stat100late);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');

grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.5 0 .5];
title(c,'Effect size \it d')
title('100% late')
% hgsave(gcf, 'gaze_sternberg.fig', '-v7.3');   % legacy equivalent

%% Helper functions

function [freq_raw, freq_norm] = computeGazeHeatmap(data, x_grid, y_grid, fs, smoothing)
    pos = horzcat(data.trial{:});
    binned = histcounts2(pos(1,:), pos(2,:), x_grid, y_grid);
    dwell_time = binned / fs;
    smoothed = imgaussfilt(dwell_time, smoothing);

    freq_raw = [];
    freq_raw.powspctrm(1,:,:) = flipud(smoothed');
    freq_raw.time = x_grid(2:end);
    freq_raw.freq = y_grid(2:end);
    freq_raw.label = {'et'};
    freq_raw.dimord = 'chan_freq_time';

    norm_time = dwell_time / sum(dwell_time(:));
    norm_smooth = imgaussfilt(norm_time, smoothing);

    freq_norm = [];
    freq_norm.powspctrm(1,:,:) = flipud(norm_smooth');
    freq_norm.time = x_grid(2:end);
    freq_norm.freq = y_grid(2:end);
    freq_norm.label = {'et'};
    freq_norm.dimord = 'chan_freq_time';
end
%% Function to interpolate power spectrum over frequencies
function smoothed_data = smooth_tfr(data, orig_freq, new_freq)
% Get dimensions
[n_channels, n_freqs, n_time] = size(data.powspctrm); % [125 × 19 × 71]

% Preallocate the new power spectrum array
powspctrm_interp = nan(n_channels, length(new_freq), n_time);

% Loop over channels and time points to interpolate each frequency spectrum
for ch = 1:n_channels
    for t = 1:n_time
        % Interpolate across the frequency dimension
        powspctrm_interp(ch, :, t) = interp1(orig_freq, squeeze(data.powspctrm(ch, :, t)), new_freq, 'spline');
    end
end

% Update the data structure
smoothed_data = data;
smoothed_data.freq = new_freq;
smoothed_data.powspctrm = powspctrm_interp;
end
% extract single subj values from mask
function vals_per_cond = mask_means_per_subject(statTFR, cond_cells)
% statTFR: struct from ft_freqstatistics with .mask [1 x F x T]
% cond_cells: {cell25, cell50, cell75, cell100}, each 1xN cell of TFR structs

    mask = squeeze(statTFR.mask(1,:,:));   % [F x T] logical
    ft   = statTFR.time(:)';               % [1 x T]
    ff   = statTFR.freq(:)';               % [1 x F]

    [Tgrid, Fgrid] = meshgrid(ft, ff);     % for interp2

    nC = numel(cond_cells);
    vals_per_cond = cell(1,nC);
    for c = 1:nC
        nc = numel(cond_cells{c});
        vals = nan(nc,1);
        for s = 1:nc
            D = cond_cells{c}{s};          % subject TFR
            P = squeeze(D.powspctrm(1,:,:));     % [F x T] or [T x F]?
            % make sure it's [F x T]
            if size(P,1) ~= numel(D.freq) || size(P,2) ~= numel(D.time)
                P = P'; % fallback if stored transposed
            end

            % Align to stat grid if needed
            if ~isequal(D.time(:)', ft) || ~isequal(D.freq(:)', ff)
                [Tsub, Fsub] = meshgrid(D.time(:)', D.freq(:)');
                P = interp2(Tsub, Fsub, P, Tgrid, Fgrid, 'linear', NaN);
            end

            vals(s) = mean(P(mask), 'omitnan');
        end
        vals_per_cond{c} = vals;           % vector length = #subjects in that contrast
    end
end
