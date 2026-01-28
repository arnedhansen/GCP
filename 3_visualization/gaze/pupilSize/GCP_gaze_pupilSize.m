%% GCP Pupil Size

%% Setup
startup
[subjects, path, colors, headmodel] = setup('GCP');

% representational modes to plot, in order
modes = {'raw','bl','pct'};

% y-labels
channels      = {'Pupil'};
channeltitles = {'Pupil Size'};
ylabs_px      = {'Pupil Size [a.u.]'};
ylabs_pct     = {'Pupil Size [%]'};
labels        = {'25% contrast','50% contrast','75% contrast','100% contrast'};

%% Loop over data modes: raw → baselined → baselined %
for m = 1:numel(modes)
    clc
    data_mode = modes{m};
    fprintf('\nGCP Pupil Size: plotting %s data\n', data_mode);

    %% Load data (per subject, pick correct representation)
    alltlk25et  = cell(1, numel(subjects));
    alltlk50et  = cell(1, numel(subjects));
    alltlk75et  = cell(1, numel(subjects));
    alltlk100et = cell(1, numel(subjects));

    for subj = 1:length(subjects)
        datapath = strcat(path, subjects{subj}, '/gaze');
        disp(['Loading Subject ', num2str(subjects{subj})])
        load([datapath, filesep 'gaze_pupil_timeseries'])

        % velTS_cXX           = timelocked average Pupil Size (no baseline)
        % velTS_cXX_bl        = timelocked average with subtractive baseline
        % velTS_cXX_bl_pct    = timelocked average with relative-change baseline (%)

        switch data_mode
            case 'raw'
                alltlk25et{subj}  = pupTS_c25;
                alltlk50et{subj}  = pupTS_c50;
                alltlk75et{subj}  = pupTS_c75;
                alltlk100et{subj} = pupTS_c100;

            case 'bl'
                alltlk25et{subj}  = pupTS_c25_bl;
                alltlk50et{subj}  = pupTS_c50_bl;
                alltlk75et{subj}  = pupTS_c75_bl;
                alltlk100et{subj} = pupTS_c100_bl;

            case 'pct'
                alltlk25et{subj}  = pupTS_c25_bl_pct;
                alltlk50et{subj}  = pupTS_c50_bl_pct;
                alltlk75et{subj}  = pupTS_c75_bl_pct;
                alltlk100et{subj} = pupTS_c100_bl_pct;
        end

    end

    %% Grand average
    cfg = [];
    cfg.keepindividual = 'yes';
    ga25et  = ft_timelockgrandaverage(cfg, alltlk25et{:});
    ga50et  = ft_timelockgrandaverage(cfg, alltlk50et{:});
    ga75et  = ft_timelockgrandaverage(cfg, alltlk75et{:});
    ga100et = ft_timelockgrandaverage(cfg, alltlk100et{:});

    %% Plot Pupil Size for 25/50/75/100% contrast (mean ± SEM)
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

    fig_ga = figure('Color','w');
    set(fig_ga, "Position", [0 0 1512 982])

    ax = gobjects(1, numel(channels));  % store axes handles

    for c = 1:numel(channels)
        ax(c) = subplot(1,1,c);
        hold(ax(c), 'on');

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
            set(hp, 'facecolor', faceC, 'edgecolor', 'none', 'facealpha', 0.4);

            hl(k) = plot(x, y, 'LineWidth', 2, 'Color', colors(k,:));
        end

        set(gca, 'FontSize', 18);
        xlabel('Time [sec]');
        ylabel(ylabs{c});
        if ismember(data_mode, {'bl','pct'})
            yline(0, '--');
        end
        title(strrep(channeltitles{c}, '_', '\_'));

        lgd = legend(hl, labels, 'Location', 'northeast', 'FontSize', 12);
        set(lgd, 'Color', 'none');
    end

    % Add title and save with appropriate suffix
    switch data_mode
        case 'raw'
            suffix = 'raw';
            sgtitle('GCP Gaze Pupil Size RAW');
        case 'bl'
            suffix = 'bl';
            sgtitle('GCP Gaze Pupil Size BASELINED');
        case 'pct'
            suffix = 'blpct';
            sgtitle('GCP Gaze Pupil Size PERCENTAGE CHANGE');
    end

    outdir_overview = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/gaze/pupil/';
    saveName = sprintf('%sGCP_gaze_pupil_overview_%s.png', outdir_overview, suffix);
    saveas(gcf, saveName);

    % Also save individual channel plots for this mode
    for c = 1:numel(channels)
        fig_single = figure('Color', 'w');
        set(fig_single, 'Position', [0 0 1512 982]);

        % copy this subplot's axes into its own figure
        ax_new = copyobj(ax(c), fig_single);
        set(ax_new, 'Position', [0.13 0.15 0.8 0.75]);  % nice full-figure layout

        % set legend
        cond_lines = findobj(ax_new, 'Type', 'line', '-not', 'LineStyle', '--');
        cond_lines = flipud(cond_lines);
        legend(ax_new, cond_lines, labels, 'Location', 'northeast', 'FontSize', 15);

        % optional: add a figure-level title indicating mode + channel
        switch data_mode
            case 'raw'
                mode_title = 'RAW';
            case 'bl'
                mode_title = 'BASELINED';
            case 'pct'
                mode_title = 'PERCENTAGE CHANGE';
        end

        sgtitle(fig_single, sprintf('GCP Gaze Pupil Size %s — %s', ...
            mode_title, channeltitles{c}));

        % build single plots filename
        outdir_single = outdir_overview;
        saveName_single = sprintf('%sGCP_gaze_pupil_%s_%s.png', outdir_single, channels{c}, suffix);
        saveas(fig_single, saveName_single);
        close(fig_single);
    end
end
