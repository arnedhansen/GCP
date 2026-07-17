%% GCP Pupil Size Time Course
% Loads dB-baselined pupil time courses from
% 2_feature_extraction/GCP_gaze_fex.m (pupTS_cXX_bl_db).

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar')

% Plot labels
channels      = {'Pupil'};
channeltitles = {'Pupil Size'};
ylabs_db      = {'Pupil Size [dB]'};
labels        = {' 25% Contrast', ' 50% Contrast', ' 75% Contrast', ' 100% Contrast'};
fontSize = 50;
t_comp   = [-2 3];
t_win    = [-0.5 2];
lineW    = 4;

% Gaussian smoothing kernel
fsample    = 500;
sigma_ms   = 50;
sigma_samp = round(sigma_ms / (1000 / fsample));
kHalf      = 3 * sigma_samp;
x_kern     = -kHalf : kHalf;
gKernel    = exp(-x_kern.^2 / (2 * sigma_samp^2));
gKernel    = gKernel / sum(gKernel);

outdir = fullfile(paths.figures, 'gaze', 'pupil');
mkdir(outdir);

%% Load data
clc
alltlk25et  = cell(1, numel(subjects));
alltlk50et  = cell(1, numel(subjects));
alltlk75et  = cell(1, numel(subjects));
alltlk100et = cell(1, numel(subjects));

for subj = 1:numel(subjects)
    datapath = fullfile(paths.features, subjects{subj}, 'gaze');
    disp(['[GCP Pupil Size] Loading Subject ', num2str(subjects{subj})])
    dat = load(fullfile(datapath, 'gaze_pupil_timeseries'));
    alltlk25et{subj}  = dat.pupTS_c25_bl_db;
    alltlk50et{subj}  = dat.pupTS_c50_bl_db;
    alltlk75et{subj}  = dat.pupTS_c75_bl_db;
    alltlk100et{subj} = dat.pupTS_c100_bl_db;
end

%% Grand average
cfg = [];
cfg.keepindividual = 'yes';
ga25et  = ft_timelockgrandaverage(cfg, alltlk25et{:});
ga50et  = ft_timelockgrandaverage(cfg, alltlk50et{:});
ga75et  = ft_timelockgrandaverage(cfg, alltlk75et{:});
ga100et = ft_timelockgrandaverage(cfg, alltlk100et{:});
ylabs = ylabs_db;

%% Figure Pupil Time Course
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on

cfg = [];
cfg.channel     = channels{1};
cfg.avgoverchan = 'yes';
cfg.latency     = t_comp;
et25  = ft_selectdata(cfg, ga25et);
et50  = ft_selectdata(cfg, ga50et);
et75  = ft_selectdata(cfg, ga75et);
et100 = ft_selectdata(cfg, ga100et);
ets = {et25, et50, et75, et100};

for k = 1:numel(ets)
    x_store = ets{k}.time(:);
    disp_idx = x_store >= t_win(1) & x_store <= t_win(2);
    x = x_store(disp_idx);

    yi = squeeze(ets{k}.individual); % [subj x time]
    if isvector(yi)
        yi = yi(:)';
    end

    yi_sm = nan(size(yi));
    for s = 1:size(yi, 1)
        yi_sm(s, :) = conv(yi(s, :), gKernel, 'same');
    end
    yi_disp = yi_sm(:, disp_idx);

    mu = nanmean(yi_disp, 1)';
    nValid = sum(~isnan(yi_disp), 1)';
    sem = nanstd(yi_disp, 0, 1)' ./ sqrt(max(nValid, 1));

    eb = shadedErrorBar(x, mu, sem, 'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(k, :), 'LineWidth', lineW);
    set(eb.patch, 'FaceColor', colors(k, :), 'FaceAlpha', 0.125);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end

xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
xlim(t_win);
ylim([-1.25 0.4])
xlabel('Time [s]', 'FontSize', fontSize*0.8);
ylabel(ylabs{1}, 'FontSize', fontSize*0.8);
leg_p = gobjects(numel(ets), 1);
for k = 1:numel(ets)
    leg_p(k) = patch(nan, nan, colors(k, :), 'FaceAlpha', 0.25, ...
        'EdgeColor', colors(k, :), 'LineWidth', 1.5);
end
set(gca, 'FontSize', fontSize*0.8);
legend(leg_p, labels, 'Location', 'best', 'FontSize', fontSize*0.65, 'Box', 'off');
box off
hold off
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(outdir, 'GCP_gaze_pupil_size_TC_db.png'), '-dpng', '-r600');