%% GCP Gaze Fixation Time Course Contrast Conditions
% Builds fixation-rate time courses per contrast condition from
% EyeLink fixation events. Fixation rate is reconstructed from fixation
% onset events (L_fixation, R_fixation) as fixations/s, then baseline-
% normalised as dB change (10*log10(stim/baseline)) and plotted with SEM.

%% Setup
startup
[subjects, paths, colors, ~] = setup('GCP');
subjects = gcp_subject_inclusion(subjects, paths);
addpath('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar')

figpath  = fullfile(paths.figures, 'gaze', 'fixations');
mkdir(figpath);

nSubj = length(subjects);

%% Parameters
fsample = 500;          % Hz (updated from data when available)
merge_tol_ms = 20;      % merge L/R onsets within +/-20 ms to avoid double counting

% Gaussian smoothing kernel for continuous rate estimation
sigma_ms   = 50;                                     % kernel SD in ms
sigma_samp = round(sigma_ms / (1000 / fsample));     % SD in samples
kHalf      = 3 * sigma_samp;
x_kern     = -kHalf : kHalf;
gKernel    = exp(-x_kern.^2 / (2 * sigma_samp^2));
gKernel    = gKernel / sum(gKernel);                 % normalise

% Computation window (wider than display to avoid edge artefacts)
t_comp     = [-2.0 2.5];                                % seconds
n_comp     = round(diff(t_comp) * fsample) + 1;         % samples in computation window
t_comp_vec = linspace(t_comp(1), t_comp(2), n_comp);    % computation time axis

% Storage window (includes baseline + display)
t_store = [-1.5 2];                                      % seconds
[~, store_start] = min(abs(t_comp_vec - t_store(1)));
[~, store_end]   = min(abs(t_comp_vec - t_store(2)));
store_idx   = store_start : store_end;
n_store     = length(store_idx);
t_store_vec = t_comp_vec(store_idx);

% Display time window (for plotting)
t_win = [-1 2];                                          % seconds
disp_idx = t_store_vec >= t_win(1) & t_store_vec <= t_win(2);
n_disp   = sum(disp_idx);
t_vec    = t_store_vec(disp_idx);

% Baseline period for dB normalisation
bl_win = [-1 -0.25];                                     % seconds
bl_idx = t_store_vec >= bl_win(1) & t_store_vec <= bl_win(2);

% Plotting
fontSize = 20;

% Condition definitions
condCodes   = {'61', '62', '63', '64'};
condLabels  = {'25% contrast', '50% contrast', '75% contrast', '100% contrast'};
nConds      = length(condCodes);

%% Process all conditions
fprintf('Processing GCP Fixation Rate Conditions');

% Preallocate: subjects x timepoints (storage window) x conditions
subjFix = nan(nSubj, n_store, nConds);

for subj = 1:nSubj
    subjMergedPath = fullfile(paths.merged, subjects{subj});

    for c = 1:nConds
        condCode = condCodes{c};

        condOnset = [];

        for block = 1:4
            clc
            fprintf('[GCP Fixations Time Course] Subject %d/%d (%s) Block %d/4 \n', subj, nSubj, subjects{subj}, block);
            mergedFile = fullfile(subjMergedPath, ...
                sprintf('%s_EEG_ET_GCP_block%d_merged.mat', subjects{subj}, block));
            if ~isfile(mergedFile)
                continue
            end

            B = load(mergedFile, 'EEG');
            if ~isfield(B, 'EEG')
                error('Merged file does not contain EEG: %s', mergedFile);
            end

            EEG = B.EEG;
            if ~isfield(EEG, 'event') || isempty(EEG.event)
                error('No event array in merged EEG for %s block %d.', subjects{subj}, block);
            end

            EEG_ep = pop_epoch(EEG, {condCode}, t_comp);
            if EEG_ep.trials < 1
                continue
            end

            fsample = EEG_ep.srate;
            nTrialSamples = EEG_ep.pnts;
            nTrialsBlock = EEG_ep.trials;
            merge_tol_samp = max(0, round(merge_tol_ms / 1000 * fsample));

            if ~all(isfield(EEG_ep.event, {'type', 'latency', 'epoch', 'duration'}))
                error(['Missing required event fields in %s block %d, cond %s. ' ...
                       'Expected: type, latency, epoch, duration.'], ...
                      subjects{subj}, block, condCode);
            end

            eventTypes = cell(1, numel(EEG_ep.event));
            for ev = 1:numel(EEG_ep.event)
                thisType = EEG_ep.event(ev).type;
                if ischar(thisType)
                    eventTypes{ev} = thisType;
                elseif isstring(thisType)
                    eventTypes{ev} = char(thisType);
                elseif iscell(thisType) && ~isempty(thisType)
                    if ischar(thisType{1})
                        eventTypes{ev} = thisType{1};
                    elseif isstring(thisType{1})
                        eventTypes{ev} = char(thisType{1});
                    else
                        error('Unsupported event.type cell content in %s block %d.', subjects{subj}, block);
                    end
                else
                    error('Unsupported event.type format in %s block %d.', subjects{subj}, block);
                end
            end

            isFixEvent = strcmp(eventTypes, 'L_fixation') | strcmp(eventTypes, 'R_fixation');

            for trl = 1:nTrialsBlock
                onsetVec = zeros(1, nTrialSamples);
                onsetCandidates = [];

                trlMask = [EEG_ep.event.epoch] == trl;
                trlFixEvents = find(isFixEvent & trlMask);

                for ev = 1:numel(trlFixEvents)
                    e = EEG_ep.event(trlFixEvents(ev));

                    latencyVal = double(e.latency);
                    durationVal = double(e.duration);
                    if ~isfinite(latencyVal) || ~isfinite(durationVal) || durationVal <= 0
                        error('Invalid fixation latency/duration in %s block %d, cond %s, trial %d.', ...
                            subjects{subj}, block, condCode, trl);
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

                condOnset(end+1, :) = onsetVec;
            end
        end

        if size(condOnset, 1) >= 3
            % Mean onset probability per sample -> fixation rate (Hz),
            % smooth over full computation window, then crop to storage.
            rateComp = mean(condOnset, 1) * fsample;
            smoothed = conv(rateComp, gKernel, 'same');
            subjFix(subj, :, c) = smoothed(store_idx);
        end
    end
end

%% dB baseline normalisation (per subject, per condition)
subjFix_db = nan(size(subjFix));
for subj = 1:nSubj
    for c = 1:nConds
        fix_ts = subjFix(subj, :, c);
        if all(isnan(fix_ts)); continue; end
        bl_mean = nanmean(fix_ts(bl_idx));
        if bl_mean <= 0 || isnan(bl_mean); continue; end
        ratio = fix_ts ./ bl_mean;
        db_ts = 10 * log10(ratio);
        db_ts(~isfinite(fix_ts) | ~isfinite(ratio) | ratio <= 0) = NaN;
        subjFix_db(subj, :, c) = db_ts;
    end
end

%% Grand averages display portion only (dB change)
subjFix_db_disp = subjFix_db(:, disp_idx, :);
grandMean_db = squeeze(nanmean(subjFix_db_disp, 1));                          % n_disp x nConds
nValid_db    = squeeze(sum(~isnan(subjFix_db_disp(:, 1, :)), 1));             % nConds x 1
grandSEM_db  = squeeze(nanstd(subjFix_db_disp, 0, 1)) ./ sqrt(nValid_db');   % n_disp x nConds

%% FIGURE Fixation Rate (dB change)
close all; figure('Position', [0 0 1512 982], 'Color', 'w');
hold on

for c = 1:nConds
    mu  = grandMean_db(:, c);
    sem = grandSEM_db(:, c);

    eb = shadedErrorBar(t_vec, mu, sem, 'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', 2.5);
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.20);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end

xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlim(t_win);
xlabel('Time [s]');
ylabel('Fixations [dB]');
leg_p_db = gobjects(nConds, 1);
for c = 1:nConds
    leg_p_db(c) = patch(nan, nan, colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.60);
end
legend(leg_p_db, condLabels, 'Location', 'northeast', 'FontSize', fontSize - 4, 'Box', 'off');
set(gca, 'FontSize', fontSize);
box off
hold off

exportgraphics(gcf, fullfile(figpath, 'GCP_gaze_fixations_rate_db.png'), 'Resolution', 600);

fprintf('\n=== All fixation figures saved to %s ===\n', figpath);
