%% GCP Oculo-Cortical Coherence (OCC) in the Gamma Band
%
% Computes magnitude-squared coherence between saccade-cleaned eye velocity
% (velOCC) and occipital EEG across contrast conditions (25/50/75/100%).
%
% Method (following the SNF proposal):
%   1. Load per-subject velOCC (saccade-zeroed, rectified) and EEG data.
%   2. Combine ET velocity and EEG channels into a single FieldTrip struct.
%   3. Compute cross-spectral density (multitaper, dpss, 3 Hz smoothing).
%   4. Derive magnitude-squared coherence between velocity and EEG channels.
%   5. Compare coherence spectra across contrast levels.
%
% Outputs:
%   - Grand-average coherence spectra per condition (gamma range)
%   - Cluster-based F-test across conditions
%   - Topography of gamma-band coherence
%   - Scatter: gamma peak frequency vs. peak OCC coherence

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('GCP');
nSubj = numel(subjects);

condNames  = {'c25','c50','c75','c100'};
condLabels = {'25%','50%','75%','100%'};
condCodes  = [61 62 63 64];
nCond = numel(condNames);

% Frequency parameters
foi_range   = [30 90];
tapsmofrq   = 3;        % spectral smoothing (Hz) for dpss tapers
foi_step    = 1;

% Time window for coherence (post-stimulus)
coh_latency = [0.3 2.0]; % avoid early evoked transient

% Occipital channel selection (for scalp-level summary)
occ_pattern = '^(O|PO|I)';

% Figure output
fig_dir = fullfile(paths.figures, 'coherence');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); endy

%% Preallocate
coh_spectra = cell(nCond, nSubj); % each: struct with .cohspctrm [chan x freq]
coh_gamma_peak = nan(nCond, nSubj);
coh_gamma_mean = nan(nCond, nSubj);

%% Subject loop
for subj = 1:nSubj
    fprintf('Subject %d/%d: %s\n', subj, nSubj, subjects{subj});

    % Load EEG
    eeg_path = fullfile(paths.features, subjects{subj}, 'eeg', 'dataEEG.mat');
    if ~isfile(eeg_path)
        warning('No EEG data for %s, skipping.', subjects{subj});
        continue
    end
    eeg_data = load(eeg_path, 'dataEEG_c25','dataEEG_c50','dataEEG_c75','dataEEG_c100');

    % Load velOCC
    gaze_path = fullfile(paths.features, subjects{subj}, 'gaze', 'gaze_velocity_timeseries.mat');
    if ~isfile(gaze_path)
        warning('No gaze velocity data for %s, skipping.', subjects{subj});
        continue
    end
    vel_data = load(gaze_path, 'velOCC_trials_c25','velOCC_trials_c50','velOCC_trials_c75','velOCC_trials_c100');

    eeg_structs = {eeg_data.dataEEG_c25, eeg_data.dataEEG_c50, ...
                   eeg_data.dataEEG_c75, eeg_data.dataEEG_c100};
    vel_structs = {vel_data.velOCC_trials_c25, vel_data.velOCC_trials_c50, ...
                   vel_data.velOCC_trials_c75, vel_data.velOCC_trials_c100};

    for ci = 1:nCond
        dataEEG = eeg_structs{ci};
        dataVel = vel_structs{ci};

        % Use only the 2D speed channel from velOCC
        cfg_sel = [];
        cfg_sel.channel = {'Vel2D'};
        dataVel = ft_selectdata(cfg_sel, dataVel);

        % Cut both to common time window
        cfg_lat = [];
        cfg_lat.latency = coh_latency;
        dataEEG_cut = ft_selectdata(cfg_lat, dataEEG);
        dataVel_cut = ft_selectdata(cfg_lat, dataVel);

        % Match trial counts (use intersection of valid trials)
        nTrlEEG = numel(dataEEG_cut.trial);
        nTrlVel = numel(dataVel_cut.trial);
        nTrl = min(nTrlEEG, nTrlVel);
        if nTrl < 10
            warning('%s %s: only %d trials, skipping.', subjects{subj}, condNames{ci}, nTrl);
            continue
        end

        cfg_trl = [];
        cfg_trl.trials = 1:nTrl;
        dataEEG_cut = ft_selectdata(cfg_trl, dataEEG_cut);
        dataVel_cut = ft_selectdata(cfg_trl, dataVel_cut);

        % Append ET velocity as an extra channel to the EEG struct
        dataCombined = dataEEG_cut;
        vel_label = dataVel_cut.label{1};
        dataCombined.label{end+1} = vel_label;
        for trl = 1:nTrl
            nSampEEG = size(dataCombined.trial{trl}, 2);
            nSampVel = numel(dataVel_cut.trial{trl});
            nSamp = min(nSampEEG, nSampVel);
            dataCombined.trial{trl} = [dataCombined.trial{trl}(:,1:nSamp); ...
                                       dataVel_cut.trial{trl}(1:nSamp)];
            dataCombined.time{trl} = dataCombined.time{trl}(1:nSamp);
        end

        % Cross-spectral density via multitaper
        cfg_freq = [];
        cfg_freq.method    = 'mtmfft';
        cfg_freq.output    = 'fourier';
        cfg_freq.foilim    = foi_range;
        cfg_freq.tapsmofrq = tapsmofrq;
        cfg_freq.keeptrials = 'yes';
        cfg_freq.channel   = 'all';
        cfg_freq.pad       = 'nextpow2';
        freq = ft_freqanalysis(cfg_freq, dataCombined);

        % Coherence (full channel x channel x freq matrix)
        cfg_coh = [];
        cfg_coh.method = 'coh';
        coh_full = ft_connectivityanalysis(cfg_coh, freq);

        % Extract coherence between Vel2D and each EEG channel
        vel_idx = find(strcmp(coh_full.label, vel_label));
        eeg_idx = setdiff(1:numel(coh_full.label), vel_idx);

        coh = [];
        coh.label = coh_full.label(eeg_idx);
        coh.freq  = coh_full.freq;
        coh.dimord = 'chan_freq';
        coh.cohspctrm = squeeze(coh_full.cohspctrm(eeg_idx, vel_idx, :));
        if size(coh.cohspctrm, 2) ~= numel(coh.freq)
            coh.cohspctrm = coh.cohspctrm';
        end

        coh_spectra{ci, subj} = coh;

        % Extract occipital mean coherence in gamma
        occ_mask = cellfun(@(l) ~isempty(regexp(l, occ_pattern, 'once')), coh.label);
        if any(occ_mask)
            occ_coh = mean(coh.cohspctrm(occ_mask, :), 1, 'omitnan');
            coh_gamma_mean(ci, subj) = mean(occ_coh, 'omitnan');
            [~, pk_idx] = max(occ_coh);
            coh_gamma_peak(ci, subj) = coh.freq(pk_idx);
        end
    end
end

fprintf('Coherence computation complete.\n');

%% Grand-average coherence spectra
% Compute mean coherence spectrum across subjects for occipital channels

ga_coh_spectra = nan(nCond, numel(coh_spectra{1,find(~cellfun(@isempty, coh_spectra(1,:)),1)}.freq));
ga_coh_sem     = nan(size(ga_coh_spectra));
freq_axis = [];

for ci = 1:nCond
    valid_subj = find(~cellfun(@isempty, coh_spectra(ci,:)));
    if isempty(valid_subj), continue; end
    freq_axis = coh_spectra{ci, valid_subj(1)}.freq;
    nFreq = numel(freq_axis);
    all_occ_coh = nan(numel(valid_subj), nFreq);
    for si = 1:numel(valid_subj)
        s = valid_subj(si);
        occ_mask = cellfun(@(l) ~isempty(regexp(l, occ_pattern, 'once')), ...
            coh_spectra{ci,s}.label);
        all_occ_coh(si,:) = mean(coh_spectra{ci,s}.cohspctrm(occ_mask,:), 1, 'omitnan');
    end
    ga_coh_spectra(ci,:) = mean(all_occ_coh, 1, 'omitnan');
    ga_coh_sem(ci,:) = std(all_occ_coh, [], 1, 'omitnan') ./ sqrt(numel(valid_subj));
end

%% Figure 1: Grand-average coherence spectra across conditions
lineColors = [1 0 0; 1 .5 0; .6 0.2 .8; 0 0 1];

fig1 = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;
hl = gobjects(1, nCond);
for ci = 1:nCond
    y = ga_coh_spectra(ci,:);
    e = ga_coh_sem(ci,:);
    faceC = 0.8*lineColors(ci,:) + 0.2;
    patch([freq_axis fliplr(freq_axis)], [y-e fliplr(y+e)], ...
        lineColors(ci,:), 'FaceColor', faceC, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hl(ci) = plot(freq_axis, y, 'Color', lineColors(ci,:), 'LineWidth', 2.5);
end
xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('Coherence (Vel2D \leftrightarrow Occipital EEG)', 'FontSize', 16);
title('Oculo-Cortical Coherence: Gamma Band', 'FontSize', 20);
legend(hl, condLabels, 'Location', 'northeast', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim(foi_range); grid on; box on;
saveas(fig1, fullfile(fig_dir, 'GCP_OCC_gamma_spectra.png'));

%% Figure 2: Topography of gamma-band coherence per condition
fig2 = figure('Position', [0 0 1512 982], 'Color', 'w');
tiledlayout(1, nCond, 'TileSpacing', 'compact', 'Padding', 'compact');

for ci = 1:nCond
    nexttile;
    valid_subj = find(~cellfun(@isempty, coh_spectra(ci,:)));
    if isempty(valid_subj), continue; end
    ref_coh = coh_spectra{ci, valid_subj(1)};

    % Build a FieldTrip-compatible freq struct for topoplot
    topo_data = [];
    topo_data.label = ref_coh.label;
    topo_data.freq = mean(freq_axis);
    topo_data.dimord = 'chan_freq';

    all_coh_topo = nan(numel(valid_subj), numel(ref_coh.label));
    for si = 1:numel(valid_subj)
        s = valid_subj(si);
        all_coh_topo(si,:) = mean(coh_spectra{ci,s}.cohspctrm, 2, 'omitnan')';
    end
    topo_data.powspctrm = mean(all_coh_topo, 1, 'omitnan')';

    cfg_topo = [];
    cfg_topo.layout = 'EEG1010.lay';
    cfg_topo.parameter = 'powspctrm';
    cfg_topo.comment = 'no';
    cfg_topo.colorbar = 'no';
    cfg_topo.figure = 'gca';
    cfg_topo.style = 'straight';
    cfg_topo.marker = 'off';
    ft_topoplotER(cfg_topo, topo_data);
    title(condLabels{ci}, 'FontSize', 16);
end
cb = colorbar; cb.Label.String = 'Coherence'; cb.FontSize = 14;
sgtitle('Scalp Topography of Gamma OCC', 'FontSize', 20);
saveas(fig2, fullfile(fig_dir, 'GCP_OCC_gamma_topography.png'));

%% Figure 3: Bar plot of mean occipital gamma coherence per condition
fig3 = figure('Position', [0 0 1512 982], 'Color', 'w');
hold on;

m_coh = mean(coh_gamma_mean, 2, 'omitnan');
e_coh = std(coh_gamma_mean, [], 2, 'omitnan') ./ sqrt(sum(isfinite(coh_gamma_mean), 2));

for ci = 1:nCond
    bar(ci, m_coh(ci), 'FaceColor', 0.8*lineColors(ci,:)+0.2, ...
        'EdgeColor', lineColors(ci,:), 'LineWidth', 1.5);
end
errorbar(1:nCond, m_coh, e_coh, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Individual subject lines
jit = (rand(nSubj,1)-0.5)*0.3;
for s = 1:nSubj
    y = coh_gamma_mean(:,s);
    if all(isnan(y)), continue; end
    plot((1:nCond)+jit(s), y, '-', 'Color', [0 0 0 0.25], 'LineWidth', 0.8);
    for ci = 1:nCond
        if isfinite(y(ci))
            scatter(ci+jit(s), y(ci), 20, lineColors(ci,:), 'filled', ...
                'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'w');
        end
    end
end

set(gca, 'XTick', 1:nCond, 'XTickLabel', condLabels, 'FontSize', 16);
ylabel('Mean Occipital Coherence (30-90 Hz)', 'FontSize', 16);
title('Gamma OCC by Contrast', 'FontSize', 20);
xlim([0.5 nCond+0.5]); grid on; box on;
saveas(fig3, fullfile(fig_dir, 'GCP_OCC_gamma_bar.png'));

%% Figure 4: Peak OCC frequency vs. gamma peak frequency (from GED)
% Load gamma peak frequencies from GED feature extraction
ged_path = fullfile(paths.data, 'features', 'GCP_eeg_GED.mat');
if isfile(ged_path)
    ged = load(ged_path);

    fig4 = figure('Position', [0 0 1512 982], 'Color', 'w');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Panel A: OCC peak frequency vs. gamma peak frequency
    nexttile; hold on;
    for ci = 1:nCond
        if isfield(ged, 'peakFreq_full')
            gf = ged.peakFreq_full(ci,:);
        elseif isfield(ged, 'gamma_freq')
            gf = ged.gamma_freq(ci,:);
        else
            warning('Cannot find gamma peak frequency field in GED output.');
            break
        end
        scatter(gf, coh_gamma_peak(ci,:), 50, lineColors(ci,:), 'filled', ...
            'MarkerFaceAlpha', 0.7);
    end
    xlabel('Gamma Peak Frequency (Hz)', 'FontSize', 14);
    ylabel('OCC Peak Frequency (Hz)', 'FontSize', 14);
    title('Gamma Freq vs. OCC Peak Freq', 'FontSize', 16);
    legend(condLabels, 'Location', 'northwest', 'FontSize', 12);
    axis equal; grid on; box on;
    plot(xlim, xlim, '--k', 'LineWidth', 1);

    % Panel B: Gamma peak frequency vs. mean OCC magnitude
    nexttile; hold on;
    for ci = 1:nCond
        if isfield(ged, 'peakFreq_full')
            gf = ged.peakFreq_full(ci,:);
        elseif isfield(ged, 'gamma_freq')
            gf = ged.gamma_freq(ci,:);
        else
            break
        end
        scatter(gf, coh_gamma_mean(ci,:), 50, lineColors(ci,:), 'filled', ...
            'MarkerFaceAlpha', 0.7);
    end
    xlabel('Gamma Peak Frequency (Hz)', 'FontSize', 14);
    ylabel('Mean Occipital Coherence', 'FontSize', 14);
    title('Gamma Freq vs. OCC Strength', 'FontSize', 16);
    legend(condLabels, 'Location', 'northwest', 'FontSize', 12);
    grid on; box on;

    sgtitle('Relationship: Gamma Oscillations and Oculo-Cortical Coherence', 'FontSize', 20);
    saveas(fig4, fullfile(fig_dir, 'GCP_OCC_vs_gamma_scatter.png'));
end

%% Statistics: Cluster-based F-test on coherence spectra across conditions
% Build FieldTrip-compatible freq structs for ft_freqstatistics

ga_freq_structs = cell(nCond, 1);
for ci = 1:nCond
    valid_subj = find(~cellfun(@isempty, coh_spectra(ci,:)));
    nValid = numel(valid_subj);
    ref = coh_spectra{ci, valid_subj(1)};

    ga_struct = [];
    ga_struct.label = ref.label;
    ga_struct.freq = ref.freq;
    ga_struct.dimord = 'subj_chan_freq';
    ga_struct.powspctrm = nan(nValid, numel(ref.label), numel(ref.freq));
    for si = 1:nValid
        ga_struct.powspctrm(si,:,:) = coh_spectra{ci, valid_subj(si)}.cohspctrm;
    end
    ga_freq_structs{ci} = ga_struct;
end

% F-test (within-subjects, 4 levels)
nValid = size(ga_freq_structs{1}.powspctrm, 1);
cfg_stat = [];
cfg_stat.method           = 'montecarlo';
cfg_stat.statistic        = 'ft_statfun_depsamplesFunivariate';
cfg_stat.correctm         = 'cluster';
cfg_stat.clusteralpha     = 0.05;
cfg_stat.clusterstatistic = 'maxsum';
cfg_stat.neighbours       = [];
cfg_stat.tail             = 1;
cfg_stat.clustertail      = 1;
cfg_stat.alpha            = 0.05;
cfg_stat.numrandomization = 1000;
cfg_stat.design(1,:) = [ones(1,nValid) 2*ones(1,nValid) 3*ones(1,nValid) 4*ones(1,nValid)];
cfg_stat.design(2,:) = [1:nValid 1:nValid 1:nValid 1:nValid];
cfg_stat.ivar = 1;
cfg_stat.uvar = 2;

statF_coh = ft_freqstatistics(cfg_stat, ga_freq_structs{1}, ga_freq_structs{2}, ...
    ga_freq_structs{3}, ga_freq_structs{4});

%% Figure 5: F-statistic spectrum with significant clusters highlighted
fig5 = figure('Position', [0 0 1512 982], 'Color', 'w');

occ_mask_stat = cellfun(@(l) ~isempty(regexp(l, occ_pattern, 'once')), statF_coh.label);
f_vals = mean(statF_coh.stat(occ_mask_stat,:), 1);
mask_sig = any(statF_coh.mask(occ_mask_stat,:), 1);

subplot(2,1,1); hold on;
plot(statF_coh.freq, f_vals, 'k', 'LineWidth', 2);
yl = ylim;
sig_runs = diff([false mask_sig false]);
starts = find(sig_runs == 1);
stops = find(sig_runs == -1) - 1;
for r = 1:numel(starts)
    patch([statF_coh.freq(starts(r)) statF_coh.freq(stops(r)) ...
           statF_coh.freq(stops(r)) statF_coh.freq(starts(r))], ...
        [yl(1) yl(1) yl(2) yl(2)], [0 0 0], ...
        'FaceAlpha', 0.12, 'EdgeColor', 'none');
end
xlabel('Frequency (Hz)', 'FontSize', 14);
ylabel('F-value', 'FontSize', 14);
title('Cluster F-test: OCC across Contrast Levels', 'FontSize', 18);
set(gca, 'FontSize', 14); grid on; box on;

% Post-hoc: pairwise differences in mean gamma coherence (Wilcoxon)
subplot(2,1,2); hold on;
pairs = nchoosek(1:nCond, 2);
pvals = nan(size(pairs,1), 1);
for p = 1:size(pairs,1)
    a = coh_gamma_mean(pairs(p,1),:)';
    b = coh_gamma_mean(pairs(p,2),:)';
    valid = isfinite(a) & isfinite(b);
    if sum(valid) >= 5
        pvals(p) = signrank(a(valid), b(valid));
    end
end
pair_labels = arrayfun(@(i) sprintf('%s vs %s', condLabels{pairs(i,1)}, condLabels{pairs(i,2)}), ...
    1:size(pairs,1), 'UniformOutput', false);
barh(1:size(pairs,1), -log10(pvals), 'FaceColor', [0.4 0.6 0.9]);
hold on; xline(-log10(0.05), 'r--', 'LineWidth', 1.5);
set(gca, 'YTick', 1:size(pairs,1), 'YTickLabel', pair_labels, 'FontSize', 13);
xlabel('-log_{10}(p)', 'FontSize', 14);
title('Pairwise Wilcoxon signed-rank tests (mean gamma OCC)', 'FontSize', 16);
grid on; box on;

saveas(fig5, fullfile(fig_dir, 'GCP_OCC_gamma_statistics.png'));

%% Save results
save(fullfile(paths.data, 'features', 'GCP_OCC_coherence.mat'), ...
    'coh_spectra', 'coh_gamma_mean', 'coh_gamma_peak', ...
    'ga_coh_spectra', 'ga_coh_sem', 'freq_axis', ...
    'statF_coh', 'coh_latency', 'foi_range', 'tapsmofrq', '-v7.3');

fprintf('All figures saved to: %s\n', fig_dir);
fprintf('Data saved to: %s\n', fullfile(paths.data, 'features', 'GCP_OCC_coherence.mat'));
