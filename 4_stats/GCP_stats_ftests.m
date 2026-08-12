%% GCP Gaze Time-Course F-tests
% Repeated-measures (dependent-samples) univariate F-tests over time for
% three gaze variables, with FieldTrip cluster permutation correction.
%
% Panels (3 rows):
%   A  F-test on Pupil Size
%   B  F-test on Microsaccade Rate
%   C  F-test on Combined Eye Velocity
%
% Null: condition means equal at each time point (25/50/75/100% contrast).
% Data: baselined % change time courses (pupil, MS, Vel2D).

%% Setup
startup
[subjects, paths, ~, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
nSubj = numel(subjects);
fprintf('[STATS FTEST] Included cohort: N = %d (%s)\n', nSubj, strjoin(subjects, ', '));

tWin = [-0.5 2];
xTicks = [-0.5 0 0.5 1 1.5 2];
lineW = 2.5;
fontSize = 18;
titleFontSize = 20;
outDir = fullfile(paths.figures, 'stats', 'f-tests');
panelLetters = {'A', 'B', 'C'};
panelTitles = { ...
    'F-test on Pupil Size', ...
    'F-test on Microsaccade Rate', ...
    'F-test on Eye Velocity'};

%% Load pupil / microsaccade / velocity subject timelocks
allPup25 = cell(1, nSubj); allPup50 = cell(1, nSubj);
allPup75 = cell(1, nSubj); allPup100 = cell(1, nSubj);
allMs25 = cell(1, nSubj);  allMs50 = cell(1, nSubj);
allMs75 = cell(1, nSubj);  allMs100 = cell(1, nSubj);
allVel25 = cell(1, nSubj); allVel50 = cell(1, nSubj);
allVel75 = cell(1, nSubj); allVel100 = cell(1, nSubj);

for subj = 1:nSubj
    gazePath = fullfile(paths.features, subjects{subj}, 'gaze');
    clc
    fprintf('[STATS FTEST] Loading %d/%d (%s)\n', subj, nSubj, subjects{subj});

    pup = load(fullfile(gazePath, 'gaze_pupil_timeseries.mat'), ...
        'pupTS_c25_bl', 'pupTS_c50_bl', 'pupTS_c75_bl', 'pupTS_c100_bl');
    allPup25{subj} = pup.pupTS_c25_bl;
    allPup50{subj} = pup.pupTS_c50_bl;
    allPup75{subj} = pup.pupTS_c75_bl;
    allPup100{subj} = pup.pupTS_c100_bl;

    ms = load(fullfile(gazePath, 'gaze_microsaccade_timeseries.mat'), ...
        'msTS_c25_bl', 'msTS_c50_bl', 'msTS_c75_bl', 'msTS_c100_bl');
    allMs25{subj} = ms.msTS_c25_bl;
    allMs50{subj} = ms.msTS_c50_bl;
    allMs75{subj} = ms.msTS_c75_bl;
    allMs100{subj} = ms.msTS_c100_bl;

    vel = load(fullfile(gazePath, 'gaze_velocity_timeseries.mat'), ...
        'velTS_c25_bl', 'velTS_c50_bl', 'velTS_c75_bl', 'velTS_c100_bl');
    allVel25{subj} = select_channel_timelock(vel.velTS_c25_bl, 'Vel2D');
    allVel50{subj} = select_channel_timelock(vel.velTS_c50_bl, 'Vel2D');
    allVel75{subj} = select_channel_timelock(vel.velTS_c75_bl, 'Vel2D');
    allVel100{subj} = select_channel_timelock(vel.velTS_c100_bl, 'Vel2D');
end

%% Grand averages (keep individuals)
cfg = [];
cfg.keepindividual = 'yes';
gaPup = { ...
    ft_timelockgrandaverage(cfg, allPup25{:}), ...
    ft_timelockgrandaverage(cfg, allPup50{:}), ...
    ft_timelockgrandaverage(cfg, allPup75{:}), ...
    ft_timelockgrandaverage(cfg, allPup100{:})};
gaMs = { ...
    ft_timelockgrandaverage(cfg, allMs25{:}), ...
    ft_timelockgrandaverage(cfg, allMs50{:}), ...
    ft_timelockgrandaverage(cfg, allMs75{:}), ...
    ft_timelockgrandaverage(cfg, allMs100{:})};
gaVel = { ...
    ft_timelockgrandaverage(cfg, allVel25{:}), ...
    ft_timelockgrandaverage(cfg, allVel50{:}), ...
    ft_timelockgrandaverage(cfg, allVel75{:}), ...
    ft_timelockgrandaverage(cfg, allVel100{:})};

%% Restrict to analysis window
cfg = [];
cfg.latency = tWin;
for k = 1:4
    gaPup{k} = ft_selectdata(cfg, gaPup{k});
    gaMs{k} = ft_selectdata(cfg, gaMs{k});
    gaVel{k} = ft_selectdata(cfg, gaVel{k});
end

%% F-tests (FieldTrip depsamples Funivariate + cluster)
statPup = run_depsamples_ftest(gaPup);
statMs = run_depsamples_ftest(gaMs);
statVel = run_depsamples_ftest(gaVel);
stats = {statPup, statMs, statVel};

%% Figure
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:3
    ax = nexttile(tl, i);
    hold(ax, 'on');
    set(ax, 'Color', 'w', 'FontSize', fontSize, 'Box', 'off');
    grid(ax, 'off');

    st = stats{i};
    t = st.time(:)';
    f = st.stat;
    if size(f, 1) > 1
        f = f(1, :);
    end
    f = f(:);

    % Cluster mask shading
    if isfield(st, 'mask') && ~isempty(st.mask)
        m = logical(st.mask);
        if size(m, 1) > 1
            m = m(1, :);
        end
        m = m(:)'; % column
        [starts, stops] = mask_runs(m);
        ylTmp = [0, max(1, max(f, [], 'omitnan') * 1.08)];
        for r = 1:numel(starts)
            patch(ax, ...
                'XData', [t(starts(r)) t(stops(r)) t(stops(r)) t(starts(r))], ...
                'YData', [ylTmp(1) ylTmp(1) ylTmp(2) ylTmp(2)], ...
                'FaceColor', [0 0 0], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
        end
    end

    plot(ax, t, f, 'k-', 'LineWidth', lineW);
    xline(ax, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'LineStyle', '--');

    xlim(ax, tWin);
    xticks(ax, xTicks);
    ymax = max(1, max(f, [], 'omitnan') * 1.08);
    ylim(ax, [0 ymax]);
    ylabel(ax, 'F-value');
    title(ax, sprintf('\\bf{%s} | %s', panelLetters{i}, panelTitles{i}), ...
        'Interpreter', 'tex', 'FontWeight', 'bold', 'FontSize', titleFontSize, ...
        'HorizontalAlignment', 'center');
    set(ax, 'TickDir', 'out');
    uistack(findobj(ax, 'Type', 'line'), 'top');
end
xlabel(tl, 'Time [s]', 'FontSize', fontSize);

outFig = fullfile(outDir, 'GCP_stats_gaze_ftests.png');
exportgraphics(gcf, outFig, 'Resolution', 600, 'BackgroundColor', 'white');
fprintf('[STATS FTEST] Saved %s\n', outFig);

%% Print cluster summary
for i = 1:3
    st = stats{i};
    fprintf('\n%s | %s\n', panelLetters{i}, panelTitles{i});
    fprintf('  max F = %.3f\n', max(st.stat(:), [], 'omitnan'));
    if ~isfield(st, 'mask') || ~any(st.mask(:))
        fprintf('  no significant clusters\n');
        continue
    end
    m = logical(st.mask);
    if size(m, 1) > 1, m = m(1, :); end
    m = m(:);
    t = st.time(:);
    [starts, stops] = mask_runs(m);
    for r = 1:numel(starts)
        fprintf('  cluster [%.3f, %.3f] s\n', t(starts(r)), t(stops(r)));
    end
end
fprintf('\n[STATS FTEST] Done.\n');

%% Local functions
function [starts, stops] = mask_runs(maskCol)
maskCol = logical(maskCol(:));
dm = diff([false; maskCol; false]);
starts = find(dm == 1);
stops = find(dm == -1) - 1;
end

function stat = run_depsamples_ftest(gaCell)
% gaCell: {ga25, ga50, ga75, ga100} with .individual
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesFunivariate';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.neighbours = [];
cfg.tail = 1;
cfg.clustertail = cfg.tail;
cfg.alpha = 0.05;
cfg.numrandomization = 10000;

nSubj = size(gaCell{1}.individual, 1);
cfg.design = zeros(2, 4 * nSubj);
cfg.design(1, :) = [ones(1, nSubj), 2 * ones(1, nSubj), ...
    3 * ones(1, nSubj), 4 * ones(1, nSubj)];
cfg.design(2, :) = repmat(1:nSubj, 1, 4);
cfg.ivar = 1;
cfg.uvar = 2;

stat = ft_timelockstatistics(cfg, gaCell{:});
end

function tlkOut = select_channel_timelock(tlkIn, channel)
cfg = [];
cfg.channel = channel;
tlkOut = ft_selectdata(cfg, tlkIn);
end
