%% GED-based TFR Visualization for GCP
%
% Loads GED-projected TFR outputs and plots:
%   1) Grand-average TFR per condition (2x2 panel)
%   2) Difference TFR (100% - 50%)

%% Setup
startup
[subjects, paths, ~] = setup('GCP');

%% Paths
in_path = fullfile(paths.features, 'GCP_eeg_GED_TFR.mat');
fig_dir = fullfile(paths.figures, 'eeg', 'tfr');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

if ~isfile(in_path)
    error('Missing GED-TFR input file: %s', in_path);
end

S = load(in_path);
required_vars = {'tfr_cond_avg', 'condLabels'};
for vi = 1:numel(required_vars)
    if ~isfield(S, required_vars{vi})
        error('Variable "%s" missing in %s', required_vars{vi}, in_path);
    end
end

condLabels = S.condLabels;
nCond = numel(condLabels);
nSubj = numel(subjects);

%% Grand averages across subjects
grand = cell(1, nCond);
for c = 1:nCond
    valid_subj_tfr = {};
    for s = 1:nSubj
        if c <= size(S.tfr_cond_avg, 1) && s <= size(S.tfr_cond_avg, 2)
            cur = S.tfr_cond_avg{c, s};
            if ~isempty(cur)
                valid_subj_tfr{end+1} = cur; %#ok<AGROW>
            end
        end
    end
    if isempty(valid_subj_tfr)
        grand{c} = [];
    else
        grand{c} = ft_freqgrandaverage([], valid_subj_tfr{:});
    end
end

% Dynamic color limits from [30..90] Hz and [0..2] s
max_spctrm = 0;
for c = 1:nCond
    if isempty(grand{c})
        continue;
    end
    fmask = grand{c}.freq >= 30 & grand{c}.freq <= 90;
    tmask = grand{c}.time >= 0 & grand{c}.time <= 2;
    if any(fmask) && any(tmask)
        subpow = grand{c}.powspctrm(:, fmask, tmask);
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
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
sgtitle('Grand Average GED Time-Frequency Representations', 'FontSize', 30, 'FontWeight', 'bold');

for c = 1:min(4, nCond)
    subplot(2, 2, c);
    if isempty(grand{c})
        axis off
        title(sprintf('%s (no data)', condLabels{c}));
        continue;
    end
    p = squeeze(grand{c}.powspctrm(1, :, :)); % freq x time
    imagesc(grand{c}.time, grand{c}.freq, p);
    axis xy;
    xlim([-0.5 2]);
    ylim([30 90]);
    colormap(color_map);
    set(gca, 'CLim', clim);
    colb = colorbar;
    colb.Label.String = 'Power [dB]';
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    yticks([30 40 50 60 70 80 90]);
    xline(0, '--');
    set(gca, 'FontSize', 20);
    title(sprintf('%s Contrast', condLabels{c}));
end

set(gcf, 'Renderer', 'painters');
exportgraphics(gcf, fullfile(fig_dir, 'GCP_eeg_tfr_GED.png'), 'Resolution', 600);

%% Difference: 100% - 50%
if nCond >= 4 && ~isempty(grand{4}) && ~isempty(grand{2})
    tfr_diff = grand{4};
    tfr_diff.powspctrm = grand{4}.powspctrm - grand{2}.powspctrm;

    figure('Position', [0 0 1512 982]);
    set(gcf, 'Color', 'w');
    p = squeeze(tfr_diff.powspctrm(1, :, :)); % freq x time
    imagesc(tfr_diff.time, tfr_diff.freq, p);
    axis xy;
    xlim([-0.5 2]);
    ylim([30 90]);
    colormap(color_map);
    set(gca, 'CLim', clim);
    colb = colorbar;
    colb.Label.String = 'Power [dB]';
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    yticks([30 40 50 60 70 80 90]);
    xline(0, '--');
    set(gca, 'FontSize', 25);
    title('GED TFR Difference: 100% Contrast - 50% Contrast');

    set(gcf, 'Renderer', 'painters');
    exportgraphics(gcf, fullfile(fig_dir, 'GCP_eeg_tfr_GED_diff.png'), 'Resolution', 600);
end
