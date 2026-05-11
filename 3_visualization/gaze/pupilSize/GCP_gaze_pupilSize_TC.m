%% GCP Pupil Size

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');

% Plot labels
channels      = {'Pupil'};
channeltitles = {'Pupil Size'};
ylabs_px      = {'Pupil Size [a.u.]'};
ylabs_pct     = {'Pupil Size [%]'};
labels        = {'25% contrast', '50% contrast', '75% contrast', '100% contrast'};
fontSize = 20;
t_comp   = [-2 3];
t_win    = [-1 2];
lineW    = 2.5;

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
    alltlk25et{subj}  = dat.pupTS_c25_bl_pct;
    alltlk50et{subj}  = dat.pupTS_c50_bl_pct;
    alltlk75et{subj}  = dat.pupTS_c75_bl_pct;
    alltlk100et{subj} = dat.pupTS_c100_bl_pct;
end

%% Grand average
cfg = [];
cfg.keepindividual = 'yes';
ga25et  = ft_timelockgrandaverage(cfg, alltlk25et{:});
ga50et  = ft_timelockgrandaverage(cfg, alltlk50et{:});
ga75et  = ft_timelockgrandaverage(cfg, alltlk75et{:});
ga100et = ft_timelockgrandaverage(cfg, alltlk100et{:});
ylabs = ylabs_pct;

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

leg_p = gobjects(numel(ets), 1);
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

    fill([x; flipud(x)], ...
        [mu + sem; flipud(mu - sem)], ...
        colors(k, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    plot(x, mu, '-', 'Color', colors(k, :), 'LineWidth', lineW);
    leg_p(k) = patch(NaN, NaN, colors(k, :), 'EdgeColor', 'none');
end

xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlim(t_win);
xlabel('Time [s]');
ylabel(ylabs{1});
legend(leg_p, labels, 'Location', 'southeast', 'FontSize', fontSize - 4, 'Box', 'off');
set(gca, 'FontSize', fontSize);
hold off
exportgraphics(gcf, fullfile(outdir, sprintf('GCP_gaze_pupil_size_TC.png')), 'Resolution', 600);