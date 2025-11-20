%% GCP Eye Velocity Stats

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');

% channels for stats and plotting
channels      = {'VelH','VelV','Vel2D'};
channeltitles = {'Horizontal Velocity','Vertical Velocity','Combined Velocity'};

modes = {'raw','bl','pct'};

for m = 1:numel(modes)
    clc
    data_mode = modes{m};
    fprintf('\nGCP Eye Velocity: stats & plotting %s data\n', data_mode);

    %% Load data (per subject, pick correct representation)
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

    %% STATS: F-test on eye velocity
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

    % Design
    nSubj = size(ga25et.individual, 1);
    design      = zeros(2, 4*nSubj);
    design(1,:) = [ones(1,nSubj), 2*ones(1,nSubj), 3*ones(1,nSubj), 4*ones(1,nSubj)];
    design(2,:) = repmat(1:nSubj, 1, 4);

    cfg.design = design;
    cfg.ivar   = 1;
    cfg.uvar   = 2;

    % Stats
    statF_eyevelocity = ft_timelockstatistics(cfg, ga25et, ga50et, ga75et, ga100et);

    %% Decide suffix + title for this mode
    switch data_mode
        case 'raw'
            suffix     = 'raw';
            mode_title = 'RAW';
        case 'bl'
            suffix     = 'bl';
            mode_title = 'BASELINED';
        case 'pct'
            suffix     = 'blpct';
            mode_title = 'PERCENTAGE CHANGE';
    end

    %% Visualise F-statistics for eye velocity
    close all

    figF = figure('Color', 'w');
    set(figF, 'Position', [0 0 1512 982]);

    tl = tiledlayout(figF, 3, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    sgtitle(tl, sprintf('GCP Eye Velocity F-test (%s)', mode_title));

    statTitles = {'F-test: Horizontal Velocity', ...
                  'F-test: Vertical Velocity', ...
                  'F-test: Combined Velocity'};

    for i = 1:3
        ax = nexttile(tl, i);
        hold(ax, 'on');
        box(ax, 'on');
        set(ax, 'Color', 'w', 'FontSize', 12);

        % plot F-values over time
        cfg = [];
        cfg.parameter = 'stat';
        cfg.figure    = 'gca';
        cfg.channel   = channels{i};
        cfg.linecolor = 'k';

        ft_singleplotER(cfg, statF_eyevelocity);

        title(ax, statTitles{i});
        xlabel(ax, 'Time [s]');
        ylabel(ax, 'F-value');
        xlim(ax, [-1 2]);
        xline(0)

        % add cluster shading based on stat.mask
        if isfield(statF_eyevelocity, 'mask')
            chIdx = find(strcmp(statF_eyevelocity.label, channels{i}), 1);

            if ~isempty(chIdx)
                t = statF_eyevelocity.time(:)';                 % [1 x time]
                m = logical(statF_eyevelocity.mask(chIdx, :));  % [1 x time]

                if any(m)
                    dm     = diff([false, m, false]);
                    starts = find(dm == 1);
                    stops  = find(dm == -1) - 1;

                    yl = ylim(ax);
                    for r = 1:numel(starts)
                        t1 = t(starts(r));
                        t2 = t(stops(r));
                        patch('XData', [t1 t2 t2 t1], ...
                              'YData', [yl(1) yl(1) yl(2) yl(2)], ...
                              'FaceColor', [0 0 0], ...
                              'FaceAlpha', 0.12, ...
                              'EdgeColor', 'none', ...
                              'Parent', ax);
                    end

                    % keep F-line on top
                    lines = findobj(ax, 'Type', 'line');
                    uistack(lines, 'top');
                end
            end
        end
    end

    %% Save figure + stats for this mode
    outdir_stats = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/gaze/velocity/';

    saveas(figF, fullfile(outdir_stats, ...
        sprintf('GCP_gaze_velocity_Ftest_%s.png', suffix)));

    save(fullfile(outdir_stats, ...
        sprintf('GCP_gaze_velocity_Ftest_%s.mat', suffix)), ...
        'statF_eyevelocity');
end
