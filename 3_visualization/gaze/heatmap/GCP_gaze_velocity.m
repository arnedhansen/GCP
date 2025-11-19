%% GCP Eye Velocity

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');

% representational modes to plot, in order
modes = {'raw','bl','pct'};

% y-labels
ylabs_px = {'Eye Velocity X [px/s]', ...
            'Eye Velocity Y [px/s]', ...
            'Combined Eye Velocity [px/s]'};

ylabs_pct = {'Eye Velocity X [%]', ...
             'Eye Velocity Y [%]', ...
             'Combined Eye Velocity [%]'};

channels      = {'VelH','VelV','Vel2D'};
channeltitles = {'horizontal velocity','vertical velocity','2D velocity'};
labels        = {'25% contrast','50% contrast','75% contrast','100% contrast'};

%% Loop over data modes: raw → baselined → baselined %
for m = 1:numel(modes)

    data_mode = modes{m};
    fprintf('\nGCP Eye Velocity: plotting %s data\n', data_mode);

    %% Load data (per subject, pick correct representation)
    clc
    alltlk25et  = cell(1, numel(subjects));
    alltlk50et  = cell(1, numel(subjects));
    alltlk75et  = cell(1, numel(subjects));
    alltlk100et = cell(1, numel(subjects));

    for subj = 1:length(subjects)
        datapath = strcat(path, subjects{subj}, '/gaze');
        disp(['Loading Subject ', num2str(subjects{subj})])
        load([datapath, filesep 'gaze_velocity_timeseries'])

        % velTS_cXX           = timelocked average velocity (no baseline)
        % velTS_cXX_bl        = timelocked average with subtractive baseline
        % velTS_cXX_bl_pct    = timelocked average with relative-change baseline (%)

        switch data_mode
            case 'raw'
                alltlk25et{subj}  = velTS_c25;
                alltlk50et{subj}  = velTS_c50;
                alltlk75et{subj}  = velTS_c75;
                alltlk100et{subj} = velTS_c100;

            case 'bl'
                alltlk25et{subj}  = velTS_c25_bl;
                alltlk50et{subj}  = velTS_c50_bl;
                alltlk75et{subj}  = velTS_c75_bl;
                alltlk100et{subj} = velTS_c100_bl;

            case 'pct'
                alltlk25et{subj}  = velTS_c25_bl_pct;
                alltlk50et{subj}  = velTS_c50_bl_pct;
                alltlk75et{subj}  = velTS_c75_bl_pct;
                alltlk100et{subj} = velTS_c100_bl_pct;

            otherwise
                error('Unknown data_mode: %s', data_mode);
        end
    end

    %% Grand average
    cfg = [];
    cfg.keepindividual = 'yes';
    ga25et  = ft_timelockgrandaverage(cfg, alltlk25et{:});
    ga50et  = ft_timelockgrandaverage(cfg, alltlk50et{:});
    ga75et  = ft_timelockgrandaverage(cfg, alltlk75et{:});
    ga100et = ft_timelockgrandaverage(cfg, alltlk100et{:});

    %% Plot Eye Velocity for 25/50/75/100% contrast (mean ± SEM)
    close all

    % choose y-labels depending on representation
    switch data_mode
        case {'raw','bl'}
            ylabs = ylabs_px;
        case 'pct'
            ylabs = ylabs_pct;
    end

    % assume equal n_subj across contrasts
    n_subj = size(ga25et.individual, 1);

    figure('Color','w');
    set(gcf, "Position", [0 0 1512/2 982])

    for c = 1:numel(channels)
        subplot(3,1,c); hold on;

        cfg = [];
        cfg.figure      = 'gcf';
        cfg.channel     = channels{c};
        cfg.avgoverchan = 'yes';
        cfg.latency     = [-1 2];

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

            hp = patch([x; x(end:-1:1); x(1)], ...
                       [low; high(end:-1:1); low(1)], ...
                       colors(k,:));
            set(hp, 'facecolor', faceC, 'edgecolor', 'none', 'facealpha', 0.6);

            hl(k) = plot(x, y, 'LineWidth', 2, 'Color', colors(k,:));
        end

        set(gca, 'FontSize', 18);
        xlabel('Time [sec]');
        ylabel(ylabs{c});
        yline(0, '--');
        title(strrep(channeltitles{c}, '_', '\_'));

        if c == 1
            lgd = legend(hl, labels, 'Location', 'northeast', 'FontSize', 12);
            set(lgd, 'Color', 'none');
        end
    end

    % Add title and save with appropriate suffix
    switch data_mode
        case 'raw'
            suffix = 'raw';
            sgtitle('GCP Gaze Velocity RAW');
        case 'bl'
            suffix = 'bl';
            sgtitle('GCP Gaze Velocity BASELINED');
        case 'pct'
            suffix = 'blpct';
            sgtitle('GCP Gaze Velocity PERCENTAGE CHANGE');
    end

    saveName = sprintf(['/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/' ...
                        'gaze/velocity/GCP_gaze_velocity_%s.png'], suffix);
    saveas(gcf, saveName);

end
