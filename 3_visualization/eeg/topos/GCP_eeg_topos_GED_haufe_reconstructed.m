%% GCP condition-specific GED topographies
%
% Figure 1 shows condition-specific Haufe activation patterns for the
% combined full-window GED signal. Each subject pattern is RMS-normalized
% before averaging, so these maps describe spatial shape rather than power.
%
% Figure 2 shows 30-90 Hz stimulus-to-baseline power change after
% multicomponent GED subspace reconstruction at the sensors. Reconstruction
% uses one pooled activation matrix per subject, allowing spatial power
% differences between conditions to be compared.

%% Setup
clear
startup
[subjects, paths] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);

fig_dir = fullfile(paths.figures, 'eeg', 'topos');
if ~isfolder(fig_dir), mkdir(fig_dir); end

ged = load(fullfile(paths.features, 'GCP_eeg_GED.mat'), ...
    'subjects', 'condLabels', 'all_topo_labels', ...
    'all_haufe_pattern_full', 'freq_reconstructed_multicomp_full');
[subject_found, subject_index] = ismember(subjects, ged.subjects);
if ~all(subject_found)
    warning('Some included subjects are absent from GCP_eeg_GED.mat and will be omitted.');
end

n_cond = numel(ged.condLabels);
haufe_subject = cell(1, n_cond);
reconstructed_subject = cell(1, n_cond);

%% Prepare subject-level FieldTrip structures
for cond = 1:n_cond
    haufe_subject{cond} = {};
    reconstructed_subject{cond} = {};

    for subj = 1:numel(subjects)
        if ~subject_found(subj)
            continue;
        end
        source_index = subject_index(subj);

        pattern = ged.all_haufe_pattern_full{cond, source_index};
        labels = ged.all_topo_labels{source_index};
        if ~isempty(pattern) && ~isempty(labels) && numel(pattern) == numel(labels)
            pattern = double(pattern(:));
            pattern_rms = sqrt(mean(pattern .^ 2, 'omitnan'));
            if isfinite(pattern_rms) && pattern_rms > eps
                pattern = pattern / pattern_rms;
                freq_pattern = [];
                freq_pattern.label = labels(:);
                freq_pattern.freq = 1;
                freq_pattern.powspctrm = pattern;
                freq_pattern.dimord = 'chan_freq';
                haufe_subject{cond}{end + 1} = freq_pattern; %#ok<SAGROW>
            end
        end

        freq_reconstructed = ged.freq_reconstructed_multicomp_full{cond, source_index};
        if ~isempty(freq_reconstructed)
            reconstructed_subject{cond}{end + 1} = freq_reconstructed; %#ok<SAGROW>
        end
    end
end

%% Grand averages
haufe_grand = cell(1, n_cond);
reconstructed_grand = cell(1, n_cond);
for cond = 1:n_cond
    if ~isempty(haufe_subject{cond})
        haufe_grand{cond} = ft_freqgrandaverage([], haufe_subject{cond}{:});
    end
    if ~isempty(reconstructed_subject{cond})
        reconstructed_grand{cond} = ft_freqgrandaverage([], reconstructed_subject{cond}{:});
    end
end

%% Plot configuration
load(fullfile(paths.base_students, 'toolboxes', ...
    'headmodel', 'layANThead.mat'), 'layANThead')

color_map = flipud(cbrewer('div', 'RdBu', 64));
cfg_common = [];
cfg_common.layout = layANThead;
cfg_common.channel = layANThead.label(1:end-2);
cfg_common.channel = cfg_common.channel(~strcmp(cfg_common.channel, 'M2'));
cfg_common.channel = cfg_common.channel(~strcmp(cfg_common.channel, 'M1'));
cfg_common.marker = 'off';
cfg_common.comment = 'no';
cfg_common.colormap = color_map;
cfg_common.gridscale = 300;

%% Condition-specific Haufe activation patterns
haufe_limit = 0;
for cond = 1:n_cond
    if ~isempty(haufe_grand{cond})
        haufe_limit = max(haufe_limit, ...
            max(abs(haufe_grand{cond}.powspctrm), [], 'all', 'omitnan'));
    end
end
if ~isfinite(haufe_limit) || haufe_limit <= 0
    haufe_limit = 1;
end

figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('Condition Specific GED Haufe Activation Patterns', ...
    'FontSize', 30, 'FontWeight', 'bold');
for cond = 1:n_cond
    subplot(2, 2, cond);
    if isempty(haufe_grand{cond})
        axis off
        title(sprintf('%s Contrast: no data', ged.condLabels{cond}));
        continue;
    end
    cfg = cfg_common;
    cfg.ylim = [1 1];
    cfg.zlim = double([-haufe_limit, haufe_limit]);
    ft_topoplotER(cfg, haufe_grand{cond});
    title(sprintf('%s Contrast', ged.condLabels{cond}), 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'RMS-normalized activation [a.u.]', 'FontSize', 20);
end
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(fig_dir, ...
    'GCP_eeg_topos_GED_haufe_condition_specific.png'), '-dpng', '-r600');

%% Reconstructed multicomponent sensor power
power_limit = 0;
for cond = 1:n_cond
    if isempty(reconstructed_grand{cond})
        continue;
    end
    freq_mask = reconstructed_grand{cond}.freq >= 30 & ...
        reconstructed_grand{cond}.freq <= 90;
    gamma_map = mean(reconstructed_grand{cond}.powspctrm(:, freq_mask), ...
        2, 'omitnan');
    power_limit = max(power_limit, max(abs(gamma_map), [], 'omitnan'));
end
if ~isfinite(power_limit) || power_limit <= 0
    power_limit = 1;
end

figure('Position', [0 0 1512 982], 'Color', 'w');
sgtitle('GED Reconstructed Sensor Gamma Power', ...
    'FontSize', 30, 'FontWeight', 'bold');
for cond = 1:n_cond
    subplot(2, 2, cond);
    if isempty(reconstructed_grand{cond})
        axis off
        title(sprintf('%s Contrast: no data', ged.condLabels{cond}));
        continue;
    end
    cfg = cfg_common;
    cfg.ylim = [30 90];
    cfg.zlim = double([-power_limit, power_limit]);
    ft_topoplotER(cfg, reconstructed_grand{cond});
    title(sprintf('%s Contrast', ged.condLabels{cond}), 'FontSize', 25);
    cb = colorbar;
    cb.FontSize = 20;
    ylabel(cb, 'Power change [dB]', 'FontSize', 20);
end
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(fig_dir, ...
    'GCP_eeg_topos_GED_reconstructed_multicomp_power.png'), '-dpng', '-r600');
