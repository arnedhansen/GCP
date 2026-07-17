%% GCP Eye Velocity Time Course
% Loads dB-baselined velocity time courses from
% 2_feature_extraction/GCP_gaze_fex.m (velTS_cXX_bl_db).

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar')

% Plot labels (three distinct velocity measures)
velocityYLabels = {'Eye Velocity X [dB]', ...
    'Eye Velocity Y [dB]', ...
    'Eye Velocity [dB]'};
channels      = {'VelH', 'VelV', 'Vel2D'};
channeltitles = {'Horizontal Velocity', 'Vertical Velocity', 'Combined Velocity'};
    labels        = {' 25% Contrast', ' 50% Contrast', ' 75% Contrast', ' 100% Contrast'};
if numel(velocityYLabels) ~= numel(channels)
    error('Velocity y-labels must match number of channels.');
end

% Plotting style aligned to microsaccade TC template
fontSize = 50;
t_store  = [-1.5 2.5];
t_win    = [-0.5 2];
lineW    = 4;

outdir = fullfile(paths.figures, 'gaze', 'velocity');
mkdir(outdir);

%% Load dB-baselined data (per subject)
clc
alltlk25et  = cell(1, numel(subjects));
alltlk50et  = cell(1, numel(subjects));
alltlk75et  = cell(1, numel(subjects));
alltlk100et = cell(1, numel(subjects));

for subj = 1:numel(subjects)
    datapath = fullfile(paths.features, subjects{subj}, 'gaze');
    disp(['[GCP Eye Velocity] Loading Subject ', num2str(subjects{subj})])
    dat = load(fullfile(datapath, 'gaze_velocity_timeseries'));
    alltlk25et{subj}  = dat.velTS_c25_bl_db;
    alltlk50et{subj}  = dat.velTS_c50_bl_db;
    alltlk75et{subj}  = dat.velTS_c75_bl_db;
    alltlk100et{subj} = dat.velTS_c100_bl_db;
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

%% Plot
close all
for c = 3%%%%%1:numel(channels)
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

        eb = shadedErrorBar(x, mu, sem, 'lineProps', {'-'}, 'transparent', true);
        set(eb.mainLine, 'Color', conditionColors(k, :), 'LineWidth', lineW);
        set(eb.patch, 'FaceColor', conditionColors(k, :), 'FaceAlpha', 0.125);
        set(eb.edge(1), 'Color', 'none');
        set(eb.edge(2), 'Color', 'none');
    end

    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
    yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
    xlim(t_win);
    xlabel('Time [s]', 'FontSize', fontSize*0.8);
    ylabel(velocityYLabels{c}, 'FontSize', fontSize*0.8);
    leg_p = gobjects(numel(ets), 1);
    for k = 1:numel(ets)
        leg_p(k) = patch(nan, nan, conditionColors(k, :), 'FaceAlpha', 0.25, ...
            'EdgeColor', conditionColors(k, :), 'LineWidth', 1.5);
    end
    set(gca, 'FontSize', fontSize*0.8);
    legend(leg_p, conditionLabels, 'Location', 'best', 'FontSize', fontSize*0.65, 'Box', 'off');
    box off
    hold off

    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, fullfile(outdir, sprintf('GCP_gaze_velocity_%s_TC_db.png', channels{c})), '-dpng', '-r600');
end
