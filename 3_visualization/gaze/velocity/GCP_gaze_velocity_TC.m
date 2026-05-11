%% GCP Eye Velocity

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');

% Plot labels (three distinct velocity measures)
velocityYLabels = {'Eye Velocity X [%]', ...
    'Eye Velocity Y [%]', ...
    'Combined Eye Velocity [%]'};
channels      = {'VelH', 'VelV', 'Vel2D'};
channeltitles = {'Horizontal Velocity', 'Vertical Velocity', 'Combined Velocity'};
labels        = {'25% contrast', '50% contrast', '75% contrast', '100% contrast'};
if numel(velocityYLabels) ~= numel(channels)
    error('Velocity y-labels must match number of channels.');
end

% Plotting style aligned to microsaccade TC template
fontSize = 20;
t_store  = [-1.5 2.5];
t_win    = [-1 2];
lineW    = 2.5;

outdir = fullfile(paths.figures, 'gaze', 'velocity');
mkdir(outdir);

%% Load percentage-change data (per subject)
clc
alltlk25et  = cell(1, numel(subjects));
alltlk50et  = cell(1, numel(subjects));
alltlk75et  = cell(1, numel(subjects));
alltlk100et = cell(1, numel(subjects));

for subj = 1:numel(subjects)
    datapath = fullfile(paths.features, subjects{subj}, 'gaze');
    disp(['[GCP Eye Velocity] Loading Subject ', num2str(subjects{subj})])
    dat = load(fullfile(datapath, 'gaze_velocity_timeseries'));
    alltlk25et{subj}  = dat.velTS_c25_bl_pct;
    alltlk50et{subj}  = dat.velTS_c50_bl_pct;
    alltlk75et{subj}  = dat.velTS_c75_bl_pct;
    alltlk100et{subj} = dat.velTS_c100_bl_pct;
end

%% Grand average
cfg = [];
cfg.keepindividual = 'yes';
ga25et  = ft_timelockgrandaverage(cfg, alltlk25et{:});
ga50et  = ft_timelockgrandaverage(cfg, alltlk50et{:});
ga75et  = ft_timelockgrandaverage(cfg, alltlk75et{:});
ga100et = ft_timelockgrandaverage(cfg, alltlk100et{:});

conditionData   = {ga25et, ga50et, ga75et, ga100et};
conditionColors = colors(1:4, :);
conditionLabels = labels(1:4);

%% Overview figure (all channels)
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
for c = 1:numel(channels)
    subplot(3, 1, c);
    hold on

    cfg = [];
    cfg.channel     = channels{c};
    cfg.avgoverchan = 'yes';
    cfg.latency     = t_store;
    ets = cell(1, numel(conditionData));
    for k = 1:numel(conditionData)
        ets{k} = ft_selectdata(cfg, conditionData{k});
    end

    leg_p = gobjects(numel(ets), 1);
    for k = 1:numel(ets)
        x_store = ets{k}.time(:);
        disp_idx = x_store >= t_win(1) & x_store <= t_win(2);
        x = x_store(disp_idx);

        yi = squeeze(ets{k}.individual); % [subj x time]
        if isvector(yi)
            yi = yi(:)';
        end
        yi_disp = yi(:, disp_idx);

        mu = nanmean(yi_disp, 1)';
        nValid = sum(~isnan(yi_disp(:, 1)));
        sem = nanstd(yi_disp, 0, 1)' ./ sqrt(max(nValid, 1));

        fill([x; flipud(x)], ...
            [mu + sem; flipud(mu - sem)], ...
            conditionColors(k, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        plot(x, mu, '-', 'Color', conditionColors(k, :), 'LineWidth', lineW);
        leg_p(k) = patch(NaN, NaN, conditionColors(k, :), 'EdgeColor', 'none');
    end

    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlim(t_win);
    xlabel('Time [s]');
    ylabel(velocityYLabels{c});
    legend(leg_p, conditionLabels, 'Location', 'northeast', 'FontSize', fontSize - 4, 'Box', 'off');
    set(gca, 'FontSize', fontSize);
    hold off
end

exportgraphics(gcf, fullfile(outdir, 'GCP_gaze_velocity_overview_TC.png'), ...
    'Resolution', 600);

%% Single-channel figures
for c = 1:numel(channels)
    figure('Position', [0 0 1512 982], 'Color', 'w');
    hold on

    cfg = [];
    cfg.channel     = channels{c};
    cfg.avgoverchan = 'yes';
    cfg.latency     = t_store;
    ets = cell(1, numel(conditionData));
    for k = 1:numel(conditionData)
        ets{k} = ft_selectdata(cfg, conditionData{k});
    end

    leg_p = gobjects(numel(ets), 1);
    for k = 1:numel(ets)
        x_store = ets{k}.time(:);
        disp_idx = x_store >= t_win(1) & x_store <= t_win(2);
        x = x_store(disp_idx);

        yi = squeeze(ets{k}.individual); % [subj x time]
        if isvector(yi)
            yi = yi(:)';
        end
        yi_disp = yi(:, disp_idx);

        mu = nanmean(yi_disp, 1)';
        nValid = sum(~isnan(yi_disp(:, 1)));
        sem = nanstd(yi_disp, 0, 1)' ./ sqrt(max(nValid, 1));

        fill([x; flipud(x)], ...
            [mu + sem; flipud(mu - sem)], ...
            conditionColors(k, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        plot(x, mu, '-', 'Color', conditionColors(k, :), 'LineWidth', lineW);
        leg_p(k) = patch(NaN, NaN, conditionColors(k, :), 'EdgeColor', 'none');
    end

    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlim(t_win);
    xlabel('Time [s]');
    ylabel(velocityYLabels{c});
    legend(leg_p, conditionLabels, 'Location', 'northeast', 'FontSize', fontSize - 4, 'Box', 'off');
    set(gca, 'FontSize', fontSize);
    hold off

    exportgraphics(gcf, fullfile(outdir, sprintf('GCP_gaze_velocity_%s_TC.png', channels{c})), ...
        'Resolution', 600);
    close(gcf);
end
