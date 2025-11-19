%% GCP Eye Velocity

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');

%% Load data
clc
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/gaze');
    disp(['Loading Subject ', num2str(subjects{subj})])
    load([datapath, filesep 'gaze_velocity_timeseries'])
    % velTS_cXX = timelocked average velocity (no baseline), FieldTrip timelock structs.
    % velTS_cXX_bl = timelocked average with subtractive baseline (baseline_period).
    % velTS_cXX_bl_pct = timelocked average with FieldTrip relative-change baseline, scaled to percent.
    % velTS_trials_cXX and velTS_pct_trials_cXX = your original per-trial velocity time series (cleaned,
    % analysis-only) and their scalar-baseline % change, for heatmaps or trial-wise analyses.

    alltlk25et{subj}  = velTS_c25;
    alltlk50et{subj}  = velTS_c50;
    alltlk75et{subj}  = velTS_c75;
    alltlk100et{subj} = velTS_c100;

        % alltlk25et{subj}  = velTS_c25_bl; % OR PERCENTAGE: velTS_c25_bl_pct
    % alltlk50et{subj}  = velTS_c50_bl;
    % alltlk75et{subj}  = velTS_c75_bl;
    % alltlk100et{subj} = velTS_c100_bl;
end

%% Grand average
cfg = [];
cfg.keepindividual='yes';
ga25et  = ft_timelockgrandaverage(cfg,alltlk25et{:});
ga50et  = ft_timelockgrandaverage(cfg,alltlk50et{:});
ga75et  = ft_timelockgrandaverage(cfg,alltlk75et{:});
ga100et = ft_timelockgrandaverage(cfg,alltlk100et{:});

%% Plot eye velocity for 25/50/75/100% contrast (mean ± SEM)
close all
channels      = {'VelH','VelV','Vel2D'};
channeltitles = {'horizontal velocity','vertical velocity','2D velocity'};
ylabs         = {'eye velocity X [px/s]', ...
                 'eye velocity Y [px/s]', ...
                 '2D eye velocity [px/s]'};

labels  = {'25% contrast','50% contrast','75% contrast','100% contrast'};

% assume equal n_subj across contrasts
n_subj = size(ga25et.individual, 1);

figure('Color','w');
set(gcf, "Position", [0 0 1512 982])

for c = 1:numel(channels)
    subplot(2,2,c); hold on;

    cfg = [];
    cfg.figure     = 'gcf';
    cfg.channel    = channels{c};
    cfg.avgoverchan = 'yes';
    cfg.latency    = [-1 2];

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

        faceC = 0.8*colors(k,:) + 0.2*[1 1 1];

        hp = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], colors(k,:));
        set(hp, 'facecolor', faceC, 'edgecolor', 'none', 'facealpha', 0.6);

        hl(k) = plot(x, y, 'LineWidth', 2, 'Color', colors(k,:));
    end

    set(gca, 'FontSize', 18);
    xlabel('Time [sec]');
    ylabel(ylabs{c});
    title(strrep(channeltitles{c}, '_', '\_'));
    %xticks([-.2 0 .2 .5 1 2]);

    if c == 1
        lgd = legend(hl, labels, 'Location', 'northeast', 'FontSize', 12);
        set(lgd, 'Color', 'none');
    end
end

saveName = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/gaze/velocity/GCP_gaze_velocity_raw.png';
saveas(gcf, saveName)