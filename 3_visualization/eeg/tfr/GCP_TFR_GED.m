%% GED-based TFR Visualization for GCP
%
% Loads GED-projected TFR outputs from GCP_eeg_GED_TFR.mat and plots:
%   1) Grand-average TFR per condition (2x2 panel)
%   2) Difference TFR (100% - 25%)
%
% Input structure (from GCP_eeg_fex_GED_TFR.m):
%   tfr_cond_avg{c, s}  FieldTrip freq struct, label = {'GED'}, baseline-corrected dB

%% Setup
startup
[subjects, paths, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);

%% Load data
in_path = fullfile(paths.features, 'GCP_eeg_GED_TFR.mat');
fig_dir = fullfile(paths.figures, 'eeg', 'tfr');
dat = load(in_path);
condLabels = dat.condLabels;
nCond = numel(condLabels);
subj_idx = arrayfun(@(s) find(strcmp(dat.subjects, subjects{s}), 1), 1:numel(subjects));
nSubj = numel(subj_idx);

%% Grand averages across subjects
grand = cell(1, nCond);
for c = 1:nCond
    valid_subj_tfr = {};
    for s = 1:nSubj
        si = subj_idx(s);
        if c <= size(dat.tfr_cond_avg, 1) && si <= size(dat.tfr_cond_avg, 2)
            cur = dat.tfr_cond_avg{c, si};
            if ~isempty(cur)
                valid_subj_tfr{end+1} = cur;
            end
        end
    end
    if isempty(valid_subj_tfr)
        grand{c} = [];
    else
        grand{c} = ft_freqgrandaverage([], valid_subj_tfr{:});
    end
end

%% Common color limits (30..90 Hz, 0..2 s)
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

%% Colormap
if exist('cbrewer', 'file') == 2
    color_map = flipud(cbrewer('div', 'RdBu', 64));
else
    color_map = flipud(parula(64));
end

%% Figure: 4 conditions
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
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(fig_dir, 'GCP_eeg_tfr_GED.png'), '-dpng', '-r600');

%% Difference: 100% - 25%
if nCond < 4 || isempty(grand{4}) || isempty(grand{1})
    warning('Cannot compute 100%% - 25%% difference: missing condition data.');
else
    tfr_diff = grand{4};
    tfr_diff.powspctrm = grand{4}.powspctrm - grand{1}.powspctrm;

    freq_idx = tfr_diff.freq >= 30 & tfr_diff.freq <= 90;
    time_idx = tfr_diff.time >= -0.5 & tfr_diff.time <= 2;
    max_spctrm_diff = max(abs(tfr_diff.powspctrm(:, freq_idx, time_idx)), [], 'all');
    if ~isfinite(max_spctrm_diff) || max_spctrm_diff <= 0
        max_spctrm_diff = 1;
    end
    clim_diff = double([-max_spctrm_diff, max_spctrm_diff]);

    cfg_diff = cfg;
    cfg_diff.figure = [];

    figure('Position', [0 0 1512 982], 'Color', 'w');
    set(gcf, 'PaperPositionMode', 'auto');
    ft_singleplotTFR(cfg_diff, tfr_diff);
    colormap(color_map);
    set(gca, 'CLim', clim_diff);
    cb = colorbar;
    ylabel(cb, 'Power [dB]', 'FontSize', 25);
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    yticks([30 40 50 60 70 80 90]);
    xticks([-0.5 0 1 2]);
    xline(0, '--');
    set(gca, 'FontSize', 25);
    title('GED TFR Difference: 100% Contrast - 25% Contrast');

    set(gcf, 'Renderer', 'painters');
    drawnow;
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, fullfile(fig_dir, 'GCP_eeg_tfr_GED_diff.png'), '-dpng', '-r600');
end
