%% GCP Gaze Event Time Courses (Saccades, Blinks, Fixations)
% Builds event-rate time courses per contrast condition from EyeLink events.
% Rate is reconstructed from onset events, Gaussian-smoothed, dB baseline-
% normalised (10*log10(stim/baseline)), and plotted with SEM shading.

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar')

outdir = fullfile(paths.figures, 'gaze', 'events');
mkdir(outdir);

nSubj = length(subjects);

%% Parameters
fsample = 500;
merge_tol_ms = 20;
saccade_blink_exclusion_ms = 100;

% Time windows
t_comp     = [-2.0 2.5];
n_comp     = round(diff(t_comp) * fsample) + 1;
t_comp_vec = linspace(t_comp(1), t_comp(2), n_comp);

t_store = [-1.5 2.5];
[~, store_start] = min(abs(t_comp_vec - t_store(1)));
[~, store_end]   = min(abs(t_comp_vec - t_store(2)));
store_idx   = store_start : store_end;
n_store     = length(store_idx);
t_store_vec = t_comp_vec(store_idx);

t_win = [-0.5 2];
disp_idx = t_store_vec >= t_win(1) & t_store_vec <= t_win(2);
t_vec    = t_store_vec(disp_idx);

bl_win = [-1.5 -0.5];
bl_idx = t_store_vec >= bl_win(1) & t_store_vec <= bl_win(2);

% Plotting style (matching velocity TC)
fontSize = 50;
lineW    = 4;

% Conditions
condCodes  = {'61', '62', '63', '64'};
condLabels = {' 25% Contrast', ' 50% Contrast', ' 75% Contrast', ' 100% Contrast'};
nConds     = length(condCodes);

% Event types to plot
eventDefs = struct( ...
    'name',      {'Saccades',    'Blinks',     'Fixations'}, ...
    'labels',    {{'L_saccade', 'R_saccade'}, {'L_blink', 'R_blink'}, {'L_fixation', 'R_fixation'}}, ...
    'yLabel',    {'Saccade Rate [dB]', 'Blink Rate [dB]', 'Fixation Rate [dB]'}, ...
    'fileTag',   {'saccades',   'blinks',     'fixations'}, ...
    'sigma_ms',  {50, 140, 50});

%% Process each event type
for ev = 1:numel(eventDefs)
    fprintf('\n=== %s ===\n', eventDefs(ev).name);

    % Event specific smoothing
    sigma_samp = max(1, round(eventDefs(ev).sigma_ms / (1000 / fsample)));
    kHalf      = 3 * sigma_samp;
    x_kern     = -kHalf : kHalf;
    gKernel    = exp(-x_kern.^2 / (2 * sigma_samp^2));
    gKernel    = gKernel / sum(gKernel);

    subjRate = nan(nSubj, n_store, nConds);

    for subj = 1:nSubj
        subjMergedPath = fullfile(paths.merged, subjects{subj});

        for c = 1:nConds
            condOnset = [];

            for block = 1:4
                clc
                fprintf('[VIZ GAZE ETEVENTS] %s | Subject %d/%d (%s) Cond %s Block %d/4\n', ...
                    eventDefs(ev).name, subj, nSubj, subjects{subj}, condCodes{c}, block);

                mergedFile = fullfile(subjMergedPath, ...
                    sprintf('%s_EEG_ET_GCP_block%d_merged.mat', subjects{subj}, block));
                if ~isfile(mergedFile)
                    continue
                end

                try
                    B = load(mergedFile, 'EEG');
                catch
                    warning('Could not load %s, skipping.', mergedFile);
                    continue
                end
                if ~isfield(B, 'EEG'), continue; end
                EEG = B.EEG;
                if ~isfield(EEG, 'event') || isempty(EEG.event), continue; end

                EEG_ep = pop_epoch(EEG, {condCodes{c}}, t_comp);
                if EEG_ep.trials < 1, continue; end

                fsample = EEG_ep.srate;
                nTrialSamples = EEG_ep.pnts;
                nTrialsBlock  = EEG_ep.trials;
                merge_tol_samp = max(0, round(merge_tol_ms / 1000 * fsample));
                sacc_blink_excl_samp = max(0, round(saccade_blink_exclusion_ms / 1000 * fsample));

                if ~all(isfield(EEG_ep.event, {'type', 'latency', 'epoch', 'duration'}))
                    continue
                end

                eventTypes = cell(1, numel(EEG_ep.event));
                for e = 1:numel(EEG_ep.event)
                    thisType = EEG_ep.event(e).type;
                    if ischar(thisType)
                        eventTypes{e} = thisType;
                    elseif isstring(thisType)
                        eventTypes{e} = char(thisType);
                    elseif iscell(thisType) && ~isempty(thisType)
                        eventTypes{e} = char(string(thisType{1}));
                    else
                        eventTypes{e} = '';
                    end
                end

                isTarget = false(1, numel(eventTypes));
                for lab = 1:numel(eventDefs(ev).labels)
                    isTarget = isTarget | strcmp(eventTypes, eventDefs(ev).labels{lab});
                end
                isBlink = strcmp(eventTypes, 'L_blink') | strcmp(eventTypes, 'R_blink');

                for trl = 1:nTrialsBlock
                    onsetVec = zeros(1, nTrialSamples);
                    onsetCandidates = [];

                    trlMask = [EEG_ep.event.epoch] == trl;
                    trlEvents = find(isTarget & trlMask);
                    trlBlinkEvents = find(isBlink & trlMask);
                    trlBlinkLat = double([EEG_ep.event(trlBlinkEvents).latency]);

                    for ei = 1:numel(trlEvents)
                        latencyVal = double(EEG_ep.event(trlEvents(ei)).latency);
                        if ~isfinite(latencyVal), continue; end
                        % Match feature extraction logic: remove saccades near blinks.
                        if strcmp(eventDefs(ev).name, 'Saccades') && ~isempty(trlBlinkLat)
                            if any(abs(latencyVal - trlBlinkLat) <= sacc_blink_excl_samp)
                                continue
                            end
                        end
                        onsetTrial = mod(round(latencyVal) - 1, nTrialSamples) + 1;
                        onsetCandidates(end+1) = onsetTrial; %#ok<AGROW>
                    end

                    if ~isempty(onsetCandidates)
                        onsetCandidates = sort(onsetCandidates);
                        keepOnsets = onsetCandidates(1);
                        for k = 2:numel(onsetCandidates)
                            if onsetCandidates(k) - keepOnsets(end) > merge_tol_samp
                                keepOnsets(end+1) = onsetCandidates(k); %#ok<AGROW>
                            end
                        end
                        onsetVec(keepOnsets) = 1;
                    end

                    if nTrialSamples >= n_comp
                        onsetVec = onsetVec(1:n_comp);
                    else
                        onsetVec(end+1 : n_comp) = 0;
                    end

                    condOnset(end+1, :) = onsetVec; %#ok<AGROW>
                end
            end

            if size(condOnset, 1) >= 3
                rateComp = mean(condOnset, 1) * fsample;
                smoothed = conv(rateComp, gKernel, 'same');
                subjRate(subj, :, c) = smoothed(store_idx);
            end
        end
    end

    %% dB baseline normalisation
    subjRate_db = nan(size(subjRate));
    for subj = 1:nSubj
        for c = 1:nConds
            ts = subjRate(subj, :, c);
            if all(isnan(ts)), continue; end
            bl_mean = nanmean(ts(bl_idx));
            if bl_mean <= 0 || isnan(bl_mean), continue; end
            ratio = ts ./ bl_mean;
            db_ts = 10 * log10(ratio);
            db_ts(~isfinite(ts) | ~isfinite(ratio) | ratio <= 0) = NaN;
            subjRate_db(subj, :, c) = db_ts;
        end
    end

    %% Grand averages (display portion)
    subjRate_db_disp = subjRate_db(:, disp_idx, :);
    grandMean = squeeze(nanmean(subjRate_db_disp, 1));
    nValid    = squeeze(sum(~isnan(subjRate_db_disp), 1));
    grandSEM  = squeeze(nanstd(subjRate_db_disp, 0, 1)) ./ sqrt(max(nValid, 1));
    grandSEM(nValid < 2) = NaN;

    %% Plot
    close all
    figure('Position', [0 0 1512 982], 'Color', 'w');
    hold on

    for c = 1:nConds
        mu  = grandMean(:, c);
        sem = grandSEM(:, c);

        eb = shadedErrorBar(t_vec, mu, sem, 'lineProps', {'-'}, 'transparent', true);
        set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', lineW);
        set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.2);
        set(eb.edge(1), 'Color', 'none');
        set(eb.edge(2), 'Color', 'none');
    end

    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
    yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
    xlim(t_win);
    xlabel('Time [s]', 'FontSize', fontSize*0.8);
    ylabel(eventDefs(ev).yLabel, 'FontSize', fontSize*0.8);

    leg_p = gobjects(nConds, 1);
    for c = 1:nConds
        leg_p(c) = patch(nan, nan, colors(c, :), 'FaceAlpha', 0.33, ...
            'EdgeColor', colors(c, :), 'LineWidth', 1.5);
    end
    set(gca, 'FontSize', fontSize*0.8);
    legend(leg_p, condLabels, 'Location', 'best', 'FontSize', fontSize*0.65, 'Box', 'off');
    box off
    hold off

    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, fullfile(outdir, sprintf('GCP_gaze_%s_rate_db_TC.png', eventDefs(ev).fileTag)), '-dpng', '-r600');
end

fprintf('\n[VIZ GAZE ETEVENTS] All event TC figures saved to %s\n', outdir);
