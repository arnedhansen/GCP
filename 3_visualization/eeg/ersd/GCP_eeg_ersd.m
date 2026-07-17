%% GCP EEG GED ERS/ERD time course
% Loads subject-condition GED TFR averages from GCP_eeg_GED_TFR.mat
% (already ERS/ERD dB-baselined [-1.5 -0.25]), averages over gamma
% (30-90 Hz), and plots mean +/- SEM time courses across subjects.
%
% GED data are a single virtual channel, so scalp topoplots are omitted.

%% Setup
startup
[subjects, paths, colors] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
figDir = fullfile(paths.figures, 'eeg', 'ersd');
if ~isfolder(figDir), mkdir(figDir); end
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar')

in_path = fullfile(paths.features, 'GCP_eeg_GED_TFR.mat');
if ~isfile(in_path)
    error('Missing GED TFR file: %s', in_path);
end
dat = load(in_path, 'tfr_cond_avg', 'subjects', 'condLabels');

condLabelsPlot = {' 25% Contrast', ' 50% Contrast', ' 75% Contrast', ' 100% Contrast'};
nCond = numel(dat.condLabels);
subj_idx = arrayfun(@(s) find(strcmp(dat.subjects, subjects{s}), 1), 1:numel(subjects));
nSubj = numel(subj_idx);

%% Collect per-subject TFR and grand-average with individuals kept
tfr_all = cell(1, nCond);
for c = 1:nCond
    tfr_all{c} = {};
    for s = 1:nSubj
        si = subj_idx(s);
        if c <= size(dat.tfr_cond_avg, 1) && si <= size(dat.tfr_cond_avg, 2)
            cur = dat.tfr_cond_avg{c, si};
            if ~isempty(cur)
                tfr_all{c}{end + 1} = cur;
            end
        end
    end
    fprintf('Condition %s: %d subjects with GED TFR\n', dat.condLabels{c}, numel(tfr_all{c}));
end

cfg = [];
cfg.keepindividual = 'yes';
ga_cond = cell(1, nCond);
for c = 1:nCond
    if isempty(tfr_all{c})
        error('No GED TFR data for condition %s', dat.condLabels{c});
    end
    ga_cond{c} = ft_freqgrandaverage(cfg, tfr_all{c}{:});
end

%% Average over gamma band
cfg = [];
cfg.frequency = [30 90];
cfg.avgoverfreq = 'yes';
ga_erd = cell(1, nCond);
for c = 1:nCond
    ga_erd{c} = ft_selectdata(cfg, ga_cond{c});
end

%% Plot ERSD time course (GED virtual channel, mean +/- SEM)
cfg = [];
cfg.channel = 'GED';
cfg.latency = [-.5 2];
tlk = cell(1, nCond);
for c = 1:nCond
    tlk{c} = ft_selectdata(cfg, ga_erd{c});
end
timeVec = tlk{1}.time(:);
nSubjPlot = size(tlk{1}.powspctrm, 1);

close all
figure('Position', [0 0 1512 982], 'Color', 'w');
fontSize = 50;
legendFontSize = fontSize * 0.666;
mask = timeVec >= -.5 & timeVec <= 2;
x = timeVec(mask);

y = cell(1, nCond);
e = cell(1, nCond);
for c = 1:nCond
    tc = squeeze(tlk{c}.powspctrm);
    if isvector(tc)
        tc = tc(:)';
    end
    y{c} = mean(tc(:, mask), 1, 'omitnan');
    e{c} = std(tc(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubjPlot);
end

hold on;
for c = 1:nCond
    eb = shadedErrorBar(x, y{c}, e{c}, 'lineProps', {'-', 'Color', colors(c, :)});
    set(eb.mainLine, 'LineWidth', 5, 'Color', colors(c, :));
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.10);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end
xline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
yline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);

set(gca, 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Power Change [dB]');
xlim([-0.5 2]);
yAll = [y{:}];
eAll = [e{:}];
yPad = 0.15 * max(abs([yAll + eAll, yAll - eAll]), [], 'omitnan');
if ~isfinite(yPad) || yPad == 0
    yPad = 0.5;
end
ylim([min(yAll - eAll, [], 'omitnan') - yPad, max(yAll + eAll, [], 'omitnan') + yPad]);

leg_p = gobjects(1, nCond);
for c = 1:nCond
    leg_p(c) = patch(NaN, NaN, colors(c, :), 'FaceAlpha', 0.25, ...
        'EdgeColor', colors(c, :), 'LineWidth', 1.5);
end
legend(leg_p, condLabelsPlot, 'Location', 'northeast', 'FontSize', legendFontSize, 'Box', 'off');
drawnow; pause(0.05);
saveas(gcf, fullfile(figDir, 'GCP_eeg_ersd_GED_timecourse.png'));
