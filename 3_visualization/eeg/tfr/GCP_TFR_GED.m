%% GED-based TFR Visualization and pairwise CBPT for GCP
%
% Loads GED-projected TFR outputs from GCP_eeg_GED_TFR.mat and produces:
%   1) Grand-average TFR per condition (2x2; display latency -0.5..2 s)
%   2) Pairwise CBPT difference TFR: 100% vs 25% (cluster outline)
%   3) All six pairwise CBPT difference TFRs (cluster outlines)
% CBPT latency is 0..2 s; difference maps still show -0.5..2 s.
%
% Input (from GCP_eeg_fex_GED.m when do_tfr is true):
%   tfr_cond_avg{c, s}  FieldTrip freq struct, label = {'GED'}, baseline-corrected dB
%
% Outputs (paths.figures/eeg/tfr):
%   GCP_eeg_tfr_GED.png
%   GCP_eeg_tfr_GED_cbpt_100vs25.png
%   GCP_eeg_tfr_GED_cbpt_pairwise.png

%% Setup
startup
[subjects, paths, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);

in_path = fullfile(paths.features, 'GCP_eeg_GED_TFR.mat');
fig_dir = fullfile(paths.figures, 'eeg', 'tfr');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

dat = load(in_path, 'tfr_cond_avg', 'subjects', 'condLabels');
condLabels = dat.condLabels;
nCond = numel(condLabels);
if nCond < 4
    error('GCP_TFR_GED:NeedFourConditions', ...
        'Expected four contrast conditions in %s.', in_path);
end

subj_idx = arrayfun(@(s) find(strcmp(dat.subjects, subjects{s}), 1), 1:numel(subjects));
nSubj = numel(subj_idx);
fprintf('[VIZ EEG TFR GED] Included cohort: N = %d\n', nSubj);

%% Collect subject TFRs (matched inclusion order)
tfr_all = cell(nCond, nSubj);
for c = 1:nCond
    nValid = 0;
    for s = 1:nSubj
        si = subj_idx(s);
        if c <= size(dat.tfr_cond_avg, 1) && si <= size(dat.tfr_cond_avg, 2)
            cur = dat.tfr_cond_avg{c, si};
            if ~isempty(cur)
                tfr_all{c, s} = cur;
                nValid = nValid + 1;
            end
        end
    end
    fprintf('[VIZ EEG TFR GED] Condition %s: %d subjects with TFR\n', ...
        condLabels{c}, nValid);
end

complete_mask = true(1, nSubj);
for c = 1:nCond
    complete_mask = complete_mask & ~cellfun(@isempty, tfr_all(c, :));
end
complete_idx = find(complete_mask);
nComplete = numel(complete_idx);
if nComplete < 2
    error('GCP_TFR_GED:TooFewSubjects', ...
        ['Need at least 2 subjects with TFR in all four conditions ' ...
        '(found %d).'], nComplete);
end
fprintf('[VIZ EEG TFR GED] Subjects with all four conditions: N = %d\n', nComplete);

%% Grand averages (all available subjects per condition; for 2x2 display)
grand = cell(1, nCond);
for c = 1:nCond
    valid = tfr_all(c, ~cellfun(@isempty, tfr_all(c, :)));
    if isempty(valid)
        grand{c} = [];
    else
        grand{c} = ft_freqgrandaverage([], valid{:});
    end
end

%% Common color limits for condition TFRs (30..90 Hz, 0..2 s)
max_spctrm = 0;
for c = 1:nCond
    if isempty(grand{c})
        continue;
    end
    freq_idx = grand{c}.freq >= 30 & grand{c}.freq <= 90;
    time_idx = grand{c}.time >= 0 & grand{c}.time <= 2;
    if any(freq_idx) && any(time_idx)
        subpow = grand{c}.powspctrm(:, freq_idx, time_idx);
        max_spctrm = max(max_spctrm, max(abs(subpow), [], 'all'));
    end
end
if ~isfinite(max_spctrm) || max_spctrm <= 0
    max_spctrm = 1;
end
clim = double([-0.9 * max_spctrm, 0.9 * max_spctrm]);

if exist('cbrewer', 'file') == 2
    color_map = flipud(cbrewer('div', 'RdBu', 64));
else
    color_map = flipud(parula(64));
end

%% Figure 1: 4 conditions (unchanged consumer path for manuscript assemble)
close all
fontSize = 20;
cfg = [];
cfg.channel = 'GED';
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = [-0.5 2];
cfg.ylim = [30 90];
cfg.shading = 'interp';

figure('Position', [0 0 1512 982], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
cfg.figure = gcf;

for c = 1:min(4, nCond)
    subplot(2, 2, c);
    if isempty(grand{c})
        axis off
        title(sprintf('%s (no data)', condLabels{c}));
        continue;
    end
    ft_singleplotTFR(cfg, grand{c});
    colormap(color_map);
    set(gca, 'CLim', clim);
    cb = colorbar;
    ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    yticks([30 40 50 60 70 80 90]);
    xticks([-0.5 0 1 2]);
    xline(0, '--');
    set(gca, 'FontSize', fontSize);
    title(sprintf('%s Contrast', condLabels{c}));
end

set(gcf, 'Renderer', 'painters');
drawnow;
print(gcf, fullfile(fig_dir, 'GCP_eeg_tfr_GED.png'), '-dpng', '-r600');
fprintf('[VIZ EEG TFR GED] Saved %s\n', fullfile(fig_dir, 'GCP_eeg_tfr_GED.png'));

%% Pairwise CBPT (higher contrast minus lower contrast)
% Condition order from feature extraction: 25, 50, 75, 100
% Plot order is 2x3: top row among lower contrasts, bottom row all 100% pairs
pair_defs = {
    2, 1, '50% vs 25%'
    3, 1, '75% vs 25%'
    3, 2, '75% vs 50%'
    4, 1, '100% vs 25%'
    4, 2, '100% vs 50%'
    4, 3, '100% vs 75%'
    };

plot_latency = [-0.5 2];
stat_latency = [0 2];
stat_frequency = [30 90];
nRand = 10000;

tfr_complete = cell(nCond, 1);
for c = 1:nCond
    tfr_complete{c} = tfr_all(c, complete_idx);
end

pair_stats = cell(1, size(pair_defs, 1));
pair_diffs = cell(1, size(pair_defs, 1));
pair_titles = cell(1, size(pair_defs, 1));

for p = 1:size(pair_defs, 1)
    i_hi = pair_defs{p, 1};
    i_lo = pair_defs{p, 2};
    pair_titles{p} = pair_defs{p, 3};
    fprintf('[VIZ EEG TFR GED] CBPT %s (N = %d, %d randomizations)\n', ...
        pair_titles{p}, nComplete, nRand);

    [pair_stats{p}, pair_diffs{p}] = run_ged_tfr_pairwise_cbpt( ...
        tfr_complete{i_hi}, tfr_complete{i_lo}, ...
        plot_latency, stat_latency, stat_frequency, nRand);

    report_tfr_cbpt(pair_stats{p}, pair_titles{p});
end

%% Shared color limits for difference maps
max_diff = 0;
for p = 1:numel(pair_diffs)
    freq_idx = pair_diffs{p}.freq >= stat_frequency(1) & ...
        pair_diffs{p}.freq <= stat_frequency(2);
    time_idx = pair_diffs{p}.time >= plot_latency(1) & ...
        pair_diffs{p}.time <= plot_latency(2);
    subpow = pair_diffs{p}.powspctrm(:, freq_idx, time_idx);
    max_diff = max(max_diff, max(abs(subpow), [], 'all'));
end
if ~isfinite(max_diff) || max_diff <= 0
    max_diff = 1;
end
clim_diff = double([-max_diff, max_diff]);

cfg_diff = [];
cfg_diff.channel = 'GED';
cfg_diff.parameter = 'powspctrm';
cfg_diff.maskparameter = 'mask';
cfg_diff.maskstyle = 'outline';
cfg_diff.colorbar = 'yes';
cfg_diff.zlim = 'maxabs';
cfg_diff.xlim = plot_latency;
cfg_diff.ylim = stat_frequency;
cfg_diff.shading = 'interp';

%% Figure 2: planned secondary 100% vs 25%
idx_100_25 = find(strcmp(pair_titles, '100% vs 25%'), 1);
figure('Position', [0 0 1512 982], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
cfg_diff.figure = gcf;
ft_singleplotTFR(cfg_diff, pair_diffs{idx_100_25});
colormap(color_map);
set(gca, 'CLim', clim_diff);
cb = colorbar;
ylabel(cb, 'Power difference [dB]', 'FontSize', 25);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
yticks([30 40 50 60 70 80 90]);
xticks([-0.5 0 1 2]);
xline(0, '--');
set(gca, 'FontSize', 25);
title(sprintf('GED TFR CBPT: %s (N = %d)', pair_titles{idx_100_25}, nComplete));
set(gcf, 'Renderer', 'painters');
drawnow;
out_100_25 = fullfile(fig_dir, 'GCP_eeg_tfr_GED_cbpt_100vs25.png');
print(gcf, out_100_25, '-dpng', '-r600');
fprintf('[VIZ EEG TFR GED] Saved %s\n', out_100_25);

%% Figure 3: all six pairwise CBPT difference TFRs
nRow = 2;
nCol = 3;
figure('Position', [0 0 1512 982], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
ax_handles = gobjects(numel(pair_diffs), 1);
cb_handles = gobjects(numel(pair_diffs), 1);
for p = 1:numel(pair_diffs)
    row = ceil(p / nCol);
    col = mod(p - 1, nCol) + 1;
    is_left = col == 1;
    is_bottom = row == nRow;
    is_right = col == nCol;

    subplot(nRow, nCol, p);
    cfg_p = cfg_diff;
    cfg_p.figure = gcf;
    cfg_p.colorbar = 'no';
    ft_singleplotTFR(cfg_p, pair_diffs{p});
    colormap(color_map);
    ax = gca;
    set(ax, 'CLim', clim_diff);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'Frequency [Hz]');
    yticks(ax, [30 50 70 90]);
    xticks(ax, [-0.5 0 1 2]);
    xline(ax, 0, '--');
    set(ax, 'FontSize', fontSize);
    title(ax, pair_titles{p});

    if ~is_left
        ylabel(ax, '');
        yticklabels(ax, []);
    end
    if ~is_bottom
        xlabel(ax, '');
        xticklabels(ax, []);
    end
    ax_pos = ax.Position;
    cb = colorbar(ax);
    ax.Position = ax_pos;
    cb.Position = [ax_pos(1) + ax_pos(3) + 0.008, ax_pos(2), 0.018, ax_pos(4)];
    if is_right
        ylabel(cb, 'dB', 'FontSize', fontSize);
    end
    ax_handles(p) = ax;
    cb_handles(p) = cb;
end
% Shrink all axes from the top so sgtitle does not overlap subplot titles
top_pad = 0.05;
for p = 1:numel(pair_diffs)
    pos = ax_handles(p).Position;
    ax_handles(p).Position = [pos(1), pos(2), pos(3), max(0.05, pos(4) - top_pad)];
    ax_pos = ax_handles(p).Position;
    cb_handles(p).Position = [ax_pos(1) + ax_pos(3) + 0.008, ax_pos(2), 0.018, ax_pos(4)];
end
sgtitle(sprintf('GED TFR pairwise CBPT (N = %d)', nComplete), ...
    'FontSize', fontSize + 4, 'FontWeight', 'bold');
set(gcf, 'Renderer', 'painters');
drawnow;
out_pairs = fullfile(fig_dir, 'GCP_eeg_tfr_GED_cbpt_pairwise.png');
print(gcf, out_pairs, '-dpng', '-r600');
fprintf('[VIZ EEG TFR GED] Saved %s\n', out_pairs);

fprintf('[VIZ EEG TFR GED] Done.\n');

%% Local helpers
function [stat, tfr_diff] = run_ged_tfr_pairwise_cbpt( ...
    tfr_hi, tfr_lo, plot_latency, stat_latency, frequency, nRand)
% Pairwise within-subject CBPT on GED TFR power (higher minus lower).
% Difference maps keep plot_latency; CBPT is restricted to stat_latency.

nSubj = numel(tfr_hi);
if nSubj ~= numel(tfr_lo)
    error('GCP_TFR_GED:SubjectMismatch', ...
        'Higher and lower contrast subject lists must match.');
end

cfg_plot = [];
cfg_plot.channel = 'GED';
cfg_plot.latency = plot_latency;
cfg_plot.frequency = frequency;
tfr_hi_plot = cell(1, nSubj);
tfr_lo_plot = cell(1, nSubj);
for s = 1:nSubj
    tfr_hi_plot{s} = ft_selectdata(cfg_plot, tfr_hi{s});
    tfr_lo_plot{s} = ft_selectdata(cfg_plot, tfr_lo{s});
end

cfg_stat = [];
cfg_stat.latency = stat_latency;
tfr_hi_sel = cell(1, nSubj);
tfr_lo_sel = cell(1, nSubj);
for s = 1:nSubj
    tfr_hi_sel{s} = ft_selectdata(cfg_stat, tfr_hi_plot{s});
    tfr_lo_sel{s} = ft_selectdata(cfg_stat, tfr_lo_plot{s});
end

cfg = [];
cfg.channel = 'GED';
cfg.latency = stat_latency;
cfg.frequency = frequency;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.neighbours = [];
cfg.minnbchan = 0;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = nRand;
cfg.avgoverchan = 'no';
cfg.avgoverfreq = 'no';
cfg.avgovertime = 'no';

cfg.design = zeros(2, 2 * nSubj);
cfg.design(1, :) = [1:nSubj, 1:nSubj];
cfg.design(2, :) = [ones(1, nSubj), 2 * ones(1, nSubj)];
cfg.uvar = 1;
cfg.ivar = 2;

stat = ft_freqstatistics(cfg, tfr_hi_sel{:}, tfr_lo_sel{:});

cfg_ga = [];
ga_hi = ft_freqgrandaverage(cfg_ga, tfr_hi_plot{:});
ga_lo = ft_freqgrandaverage(cfg_ga, tfr_lo_plot{:});
tfr_diff = ga_hi;
tfr_diff.powspctrm = ga_hi.powspctrm - ga_lo.powspctrm;
tfr_diff.mask = false(size(tfr_diff.powspctrm));
if isfield(stat, 'mask') && ~isempty(stat.mask)
    [~, iFreq] = ismember(stat.freq, tfr_diff.freq);
    [~, iTime] = ismember(stat.time, tfr_diff.time);
    if any(iFreq == 0) || any(iTime == 0)
        error('GCP_TFR_GED:MaskAlign', ...
            'CBPT mask frequency/time grid does not match plot grid.');
    end
    tfr_diff.mask(:, iFreq, iTime) = logical(stat.mask);
end
end

function report_tfr_cbpt(stat, label)
fprintf('  [%s] ', label);
if ~isfield(stat, 'mask') || ~any(stat.mask(:))
    fprintf('no significant clusters\n');
    return;
end
if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
    ppos = [stat.posclusters.prob];
    fprintf('pos clusters: %s; ', mat2str(ppos, 3));
end
if isfield(stat, 'negclusters') && ~isempty(stat.negclusters)
    pneg = [stat.negclusters.prob];
    fprintf('neg clusters: %s; ', mat2str(pneg, 3));
end
fprintf('mask voxels: %d\n', nnz(stat.mask));
end
