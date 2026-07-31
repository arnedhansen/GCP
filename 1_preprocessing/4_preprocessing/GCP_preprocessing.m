%% Preprocessing GCP

%% Setup
startup
[subjects, paths] = setup('GCP');
mergedPath = paths.merged;

% EyeLink event-rate windows (must match GCP_gaze_fex)
baseline_window = [-1.5 -0.5];
analysis_full   = [0 2];
analysis_early  = [0 1];
analysis_late   = [1 2];
eyelink_wins    = {analysis_full, analysis_early, analysis_late};
epoch_window    = [-2 3.5];

%% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    clearvars -except subjects subj mergedPath paths ...
        baseline_window analysis_full analysis_early analysis_late eyelink_wins epoch_window
    clc; fprintf('[PREP] Subject %s (%d/%d)\n', subjects{subj}, subj, length(subjects));
    datapath = fullfile(mergedPath, subjects{subj});
    cd(datapath)

    eegFile = fullfile(paths.features, subjects{subj}, 'eeg', 'dataEEG.mat');
    needEEG = isempty(dir(eegFile));

    %% Read blocks (needed for EEG epoching and/or EyeLink rates)
    alleeg = cell(1, 4);
    for block = 1:4
        try
            load(sprintf('%s_EEG_ET_GCP_block%d_merged.mat', subjects{subj}, block))
            alleeg{block} = EEG;
            clear EEG
            fprintf('[PREP] Subject %s (%d/%d): Block %d loaded\n', subjects{subj}, subj, length(subjects), block)
        catch ME
            ME.message
            fprintf('[PREP] ERROR loading Block %d!\n', block)
        end
    end

    if needEEG
        %% Segment data into epochs -2s before and 3.5s after stim onset and
        %  convert to FieldTrip data structure
        % 51 = PRESENTATION_C25_TASK    (Trigger for presentation of 25% contrast concentric dynamic inward grating WITH button press response)
        % 52 = PRESENTATION_C50_TASK    (Trigger for presentation of 50% contrast concentric dynamic inward grating WITH button press response)
        % 53 = PRESENTATION_C75_TASK    (Trigger for presentation of 75% contrast concentric dynamic inward grating WITH button press response)
        % 54 = PRESENTATION_C100_TASK   (Trigger for presentation of 100% contrast concentric dynamic inward grating WITH button press response)
        % 61 = PRESENTATION_C25_NOTASK  (Trigger for presentation of 25% contrast concentric dynamic inward grating WITHOUT button press response)
        % 62 = PRESENTATION_C50_NOTASK  (Trigger for presentation of 50% contrast concentric dynamic inward grating WITHOUT button press response)
        % 63 = PRESENTATION_C75_NOTASK  (Trigger for presentation of 75% contrast concentric dynamic inward grating WITHOUT button press response)
        % 64 = PRESENTATION_C100_NOTASK (Trigger for presentation of 100% contrast concentric dynamic inward grating WITHOUT button press response)
        data_c25 = cell(1, 4);
        data_c50 = cell(1, 4);
        data_c75 = cell(1, 4);
        data_c100 = cell(1, 4);
        for block = 1:4
            if isempty(alleeg{block}), continue; end
            try
                EEG_c25 = pop_epoch(alleeg{block}, {'61'}, epoch_window);
                data_c25{block} = eeglab2fieldtrip(EEG_c25, 'raw');

                EEG_c50 = pop_epoch(alleeg{block}, {'62'}, epoch_window);
                data_c50{block} = eeglab2fieldtrip(EEG_c50, 'raw');

                EEG_c75 = pop_epoch(alleeg{block}, {'63'}, epoch_window);
                data_c75{block} = eeglab2fieldtrip(EEG_c75, 'raw');

                EEG_c100 = pop_epoch(alleeg{block}, {'64'}, epoch_window);
                data_c100{block} = eeglab2fieldtrip(EEG_c100, 'raw');
            catch ME
                ME.message
                fprintf('[PREP] ERROR segmenting Block %d!\n', block)
            end
        end

        %% Remove empty blocks
        data_c25 = data_c25(~cellfun('isempty', data_c25));
        data_c50 = data_c50(~cellfun('isempty', data_c50));
        data_c75 = data_c75(~cellfun('isempty', data_c75));
        data_c100 = data_c100(~cellfun('isempty', data_c100));

        %% Equalize labels
        data_c25 = update_labels(data_c25);
        data_c50 = update_labels(data_c50);
        data_c75 = update_labels(data_c75);
        data_c100 = update_labels(data_c100);

        %% Add trialinfo
        for block = 1:numel(data_c25)
            try
                data_c25{block}.trialinfo = zeros(numel(data_c25{block}.trial), 1) + 61;
            catch ME
                ME.message
                fprintf('[PREP] ERROR adding trialinfo (c25) in Block %d!\n', block)
            end
        end
        for block = 1:numel(data_c50)
            try
                data_c50{block}.trialinfo = zeros(numel(data_c50{block}.trial), 1) + 62;
            catch ME
                ME.message
                fprintf('[PREP] ERROR adding trialinfo (c50) in Block %d!\n', block)
            end
        end
        for block = 1:numel(data_c75)
            try
                data_c75{block}.trialinfo = zeros(numel(data_c75{block}.trial), 1) + 63;
            catch ME
                ME.message
                fprintf('[PREP] ERROR adding trialinfo (c75) in Block %d!\n', block)
            end
        end
        for block = 1:numel(data_c100)
            try
                data_c100{block}.trialinfo = zeros(numel(data_c100{block}.trial), 1) + 64;
            catch ME
                ME.message
                fprintf('[PREP] ERROR adding trialinfo (c100) in Block %d!\n', block)
            end
        end

        %% Append data for conditions
        cfg = [];
        cfg.keepsampleinfo = 'yes';
        data_c25 = ft_appenddata(cfg, data_c25{:});
        data_c50 = ft_appenddata(cfg, data_c50{:});
        data_c75 = ft_appenddata(cfg, data_c75{:});
        data_c100 = ft_appenddata(cfg, data_c100{:});

        %% Select EyeTracking data
        cfg = [];
        cfg.channel = {'L-GAZE-X'  'L-GAZE-Y' 'L-AREA', 'R-GAZE-X'  'R-GAZE-Y' 'R-AREA'};
        dataET_c25 = ft_selectdata(cfg, data_c25);
        dataET_c50 = ft_selectdata(cfg, data_c50);
        dataET_c75 = ft_selectdata(cfg, data_c75);
        dataET_c100 = ft_selectdata(cfg, data_c100);

        %% Select EEG data (excl. ET and EOG data)
        cfg = [];
        cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' , '-L-AREA', '-R-GAZE-X'  '-R-GAZE-Y' '-R-AREA'};
        dataEEG_c25 = ft_selectdata(cfg, data_c25);
        dataEEG_c50 = ft_selectdata(cfg, data_c50);
        dataEEG_c75 = ft_selectdata(cfg, data_c75);
        dataEEG_c100 = ft_selectdata(cfg, data_c100);

        %% Re-reference data to average/common reference
        cfg = [];
        cfg.reref   = 'yes';
        cfg.refchannel = 'all';
        dataEEG_c25 = ft_preprocessing(cfg, dataEEG_c25);
        dataEEG_c50 = ft_preprocessing(cfg, dataEEG_c50);
        dataEEG_c75 = ft_preprocessing(cfg, dataEEG_c75);
        dataEEG_c100 = ft_preprocessing(cfg, dataEEG_c100);

        %% Save EEG / ET
        savepath = fullfile(paths.features, subjects{subj}, 'eeg');
        mkdir(savepath)
        cd(savepath)
        save dataEEG dataEEG_c25 dataEEG_c50 dataEEG_c75 dataEEG_c100
        savepathET = fullfile(paths.features, subjects{subj}, 'gaze');
        mkdir(savepathET)
        cd(savepathET)
        save dataET dataET_c25 dataET_c50 dataET_c75 dataET_c100
    end

    %% Multi-window EyeLink event rates (blinks, fixations, saccades)
    % Windows: full [0 2], early [0 1], late [1 2], baseline [-1.5 -0.5]
    % Rates in Hz. Saccades within 100 ms of a blink are excluded.
    ev = eyelink_event_rates_from_alleeg(alleeg, eyelink_wins, baseline_window);

    c25_blinks = ev.blinks(1,1); c50_blinks = ev.blinks(2,1);
    c75_blinks = ev.blinks(3,1); c100_blinks = ev.blinks(4,1);
    c25_fixations = ev.fixations(1,1); c50_fixations = ev.fixations(2,1);
    c75_fixations = ev.fixations(3,1); c100_fixations = ev.fixations(4,1);
    c25_saccades = ev.saccades(1,1); c50_saccades = ev.saccades(2,1);
    c75_saccades = ev.saccades(3,1); c100_saccades = ev.saccades(4,1);

    c25_blinks_early = ev.blinks(1,2); c50_blinks_early = ev.blinks(2,2);
    c75_blinks_early = ev.blinks(3,2); c100_blinks_early = ev.blinks(4,2);
    c25_fixations_early = ev.fixations(1,2); c50_fixations_early = ev.fixations(2,2);
    c75_fixations_early = ev.fixations(3,2); c100_fixations_early = ev.fixations(4,2);
    c25_saccades_early = ev.saccades(1,2); c50_saccades_early = ev.saccades(2,2);
    c75_saccades_early = ev.saccades(3,2); c100_saccades_early = ev.saccades(4,2);

    c25_blinks_late = ev.blinks(1,3); c50_blinks_late = ev.blinks(2,3);
    c75_blinks_late = ev.blinks(3,3); c100_blinks_late = ev.blinks(4,3);
    c25_fixations_late = ev.fixations(1,3); c50_fixations_late = ev.fixations(2,3);
    c75_fixations_late = ev.fixations(3,3); c100_fixations_late = ev.fixations(4,3);
    c25_saccades_late = ev.saccades(1,3); c50_saccades_late = ev.saccades(2,3);
    c75_saccades_late = ev.saccades(3,3); c100_saccades_late = ev.saccades(4,3);

    c25_bl_blinks = ev.bl_blinks(1); c50_bl_blinks = ev.bl_blinks(2);
    c75_bl_blinks = ev.bl_blinks(3); c100_bl_blinks = ev.bl_blinks(4);
    c25_bl_fixations = ev.bl_fixations(1); c50_bl_fixations = ev.bl_fixations(2);
    c75_bl_fixations = ev.bl_fixations(3); c100_bl_fixations = ev.bl_fixations(4);
    c25_bl_saccades = ev.bl_saccades(1); c50_bl_saccades = ev.bl_saccades(2);
    c75_bl_saccades = ev.bl_saccades(3); c100_bl_saccades = ev.bl_saccades(4);

    c25_pct_blinks = compute_pct_baseline(c25_blinks, c25_bl_blinks);
    c50_pct_blinks = compute_pct_baseline(c50_blinks, c50_bl_blinks);
    c75_pct_blinks = compute_pct_baseline(c75_blinks, c75_bl_blinks);
    c100_pct_blinks = compute_pct_baseline(c100_blinks, c100_bl_blinks);
    c25_pct_blinks_early = compute_pct_baseline(c25_blinks_early, c25_bl_blinks);
    c50_pct_blinks_early = compute_pct_baseline(c50_blinks_early, c50_bl_blinks);
    c75_pct_blinks_early = compute_pct_baseline(c75_blinks_early, c75_bl_blinks);
    c100_pct_blinks_early = compute_pct_baseline(c100_blinks_early, c100_bl_blinks);
    c25_pct_blinks_late = compute_pct_baseline(c25_blinks_late, c25_bl_blinks);
    c50_pct_blinks_late = compute_pct_baseline(c50_blinks_late, c50_bl_blinks);
    c75_pct_blinks_late = compute_pct_baseline(c75_blinks_late, c75_bl_blinks);
    c100_pct_blinks_late = compute_pct_baseline(c100_blinks_late, c100_bl_blinks);

    c25_pct_fixations = compute_pct_baseline(c25_fixations, c25_bl_fixations);
    c50_pct_fixations = compute_pct_baseline(c50_fixations, c50_bl_fixations);
    c75_pct_fixations = compute_pct_baseline(c75_fixations, c75_bl_fixations);
    c100_pct_fixations = compute_pct_baseline(c100_fixations, c100_bl_fixations);
    c25_pct_fixations_early = compute_pct_baseline(c25_fixations_early, c25_bl_fixations);
    c50_pct_fixations_early = compute_pct_baseline(c50_fixations_early, c50_bl_fixations);
    c75_pct_fixations_early = compute_pct_baseline(c75_fixations_early, c75_bl_fixations);
    c100_pct_fixations_early = compute_pct_baseline(c100_fixations_early, c100_bl_fixations);
    c25_pct_fixations_late = compute_pct_baseline(c25_fixations_late, c25_bl_fixations);
    c50_pct_fixations_late = compute_pct_baseline(c50_fixations_late, c50_bl_fixations);
    c75_pct_fixations_late = compute_pct_baseline(c75_fixations_late, c75_bl_fixations);
    c100_pct_fixations_late = compute_pct_baseline(c100_fixations_late, c100_bl_fixations);

    c25_pct_saccades = compute_pct_baseline(c25_saccades, c25_bl_saccades);
    c50_pct_saccades = compute_pct_baseline(c50_saccades, c50_bl_saccades);
    c75_pct_saccades = compute_pct_baseline(c75_saccades, c75_bl_saccades);
    c100_pct_saccades = compute_pct_baseline(c100_saccades, c100_bl_saccades);
    c25_pct_saccades_early = compute_pct_baseline(c25_saccades_early, c25_bl_saccades);
    c50_pct_saccades_early = compute_pct_baseline(c50_saccades_early, c50_bl_saccades);
    c75_pct_saccades_early = compute_pct_baseline(c75_saccades_early, c75_bl_saccades);
    c100_pct_saccades_early = compute_pct_baseline(c100_saccades_early, c100_bl_saccades);
    c25_pct_saccades_late = compute_pct_baseline(c25_saccades_late, c25_bl_saccades);
    c50_pct_saccades_late = compute_pct_baseline(c50_saccades_late, c50_bl_saccades);
    c75_pct_saccades_late = compute_pct_baseline(c75_saccades_late, c75_bl_saccades);
    c100_pct_saccades_late = compute_pct_baseline(c100_saccades_late, c100_bl_saccades);

    savepathET = fullfile(paths.features, subjects{subj}, 'gaze');
    mkdir(savepathET)
    save(fullfile(savepathET, 'gaze_metrics'), ...
        'c25_blinks', 'c50_blinks', 'c75_blinks', 'c100_blinks', ...
        'c25_fixations', 'c50_fixations', 'c75_fixations', 'c100_fixations', ...
        'c25_saccades', 'c50_saccades', 'c75_saccades', 'c100_saccades', ...
        'c25_blinks_early', 'c50_blinks_early', 'c75_blinks_early', 'c100_blinks_early', ...
        'c25_fixations_early', 'c50_fixations_early', 'c75_fixations_early', 'c100_fixations_early', ...
        'c25_saccades_early', 'c50_saccades_early', 'c75_saccades_early', 'c100_saccades_early', ...
        'c25_blinks_late', 'c50_blinks_late', 'c75_blinks_late', 'c100_blinks_late', ...
        'c25_fixations_late', 'c50_fixations_late', 'c75_fixations_late', 'c100_fixations_late', ...
        'c25_saccades_late', 'c50_saccades_late', 'c75_saccades_late', 'c100_saccades_late', ...
        'c25_bl_blinks', 'c50_bl_blinks', 'c75_bl_blinks', 'c100_bl_blinks', ...
        'c25_bl_fixations', 'c50_bl_fixations', 'c75_bl_fixations', 'c100_bl_fixations', ...
        'c25_bl_saccades', 'c50_bl_saccades', 'c75_bl_saccades', 'c100_bl_saccades', ...
        'c25_pct_blinks', 'c50_pct_blinks', 'c75_pct_blinks', 'c100_pct_blinks', ...
        'c25_pct_blinks_early', 'c50_pct_blinks_early', 'c75_pct_blinks_early', 'c100_pct_blinks_early', ...
        'c25_pct_blinks_late', 'c50_pct_blinks_late', 'c75_pct_blinks_late', 'c100_pct_blinks_late', ...
        'c25_pct_fixations', 'c50_pct_fixations', 'c75_pct_fixations', 'c100_pct_fixations', ...
        'c25_pct_fixations_early', 'c50_pct_fixations_early', 'c75_pct_fixations_early', 'c100_pct_fixations_early', ...
        'c25_pct_fixations_late', 'c50_pct_fixations_late', 'c75_pct_fixations_late', 'c100_pct_fixations_late', ...
        'c25_pct_saccades', 'c50_pct_saccades', 'c75_pct_saccades', 'c100_pct_saccades', ...
        'c25_pct_saccades_early', 'c50_pct_saccades_early', 'c75_pct_saccades_early', 'c100_pct_saccades_early', ...
        'c25_pct_saccades_late', 'c50_pct_saccades_late', 'c75_pct_saccades_late', 'c100_pct_saccades_late', ...
        'analysis_full', 'analysis_early', 'analysis_late', 'baseline_window');

    clc
    if subj == length(subjects)
        fprintf('[PREP] Subject %s (%d/%d) done. PREPROCESSING FINALIZED.\n', subjects{subj}, subj, length(subjects))
    else
        fprintf('[PREP] Subject %s (%d/%d) done. Loading next subject...\n', subjects{subj}, subj, length(subjects))
    end
end

function ev = eyelink_event_rates_from_alleeg(alleeg, winList, baselineWin)
% Rates [Hz] per condition (rows 1..4) x analysis window (cols), plus baseline.
nWin = numel(winList);
nCond = 4;
condCodes = {'61','62','63','64'};
ev = struct();
ev.blinks = nan(nCond, nWin);
ev.fixations = nan(nCond, nWin);
ev.saccades = nan(nCond, nWin);
ev.bl_blinks = nan(nCond, 1);
ev.bl_fixations = nan(nCond, 1);
ev.bl_saccades = nan(nCond, 1);

for c = 1:nCond
    counts = zeros(nWin + 1, 3);
    nTrials = zeros(nWin + 1, 1);
    for block = 1:numel(alleeg)
        if isempty(alleeg{block}) || ~isfield(alleeg{block}, 'event') || isempty(alleeg{block}.event)
            continue
        end
        EEG = alleeg{block};
        for wi = 1:(nWin + 1)
            if wi <= nWin
                tw = winList{wi};
            else
                tw = baselineWin;
            end
            try
                EEG_ep = pop_epoch(EEG, condCodes(c), tw);
            catch
                continue
            end
            if EEG_ep.trials < 1, continue; end
            nTrials(wi) = nTrials(wi) + EEG_ep.trials;
            [nb, nf, ns] = count_eyelink_events(EEG_ep);
            counts(wi, :) = counts(wi, :) + [nb, nf, ns];
        end
    end
    for wi = 1:nWin
        dur = diff(winList{wi});
        if nTrials(wi) > 0 && dur > 0
            rates = counts(wi, :) ./ (nTrials(wi) * dur);
            ev.blinks(c, wi) = rates(1);
            ev.fixations(c, wi) = rates(2);
            ev.saccades(c, wi) = rates(3);
        end
    end
    durBl = diff(baselineWin);
    if nTrials(end) > 0 && durBl > 0
        ratesBl = counts(end, :) ./ (nTrials(end) * durBl);
        ev.bl_blinks(c) = ratesBl(1);
        ev.bl_fixations(c) = ratesBl(2);
        ev.bl_saccades(c) = ratesBl(3);
    end
end
end

function [nBlink, nFix, nSacc] = count_eyelink_events(EEG_ep)
types = cell(1, numel(EEG_ep.event));
for ev = 1:numel(EEG_ep.event)
    t = EEG_ep.event(ev).type;
    if ischar(t)
        types{ev} = t;
    elseif isstring(t)
        types{ev} = char(t);
    elseif iscell(t) && ~isempty(t)
        types{ev} = char(string(t{1}));
    else
        types{ev} = '';
    end
end
isBlink = strcmp(types, 'L_blink') | strcmp(types, 'R_blink');
isFix = strcmp(types, 'L_fixation') | strcmp(types, 'R_fixation');
isSacc = strcmp(types, 'L_saccade') | strcmp(types, 'R_saccade');
nBlink = sum(isBlink);
nFix = sum(isFix);
blinkLat = [EEG_ep.event(isBlink).latency];
nSacc = 0;
saccIdx = find(isSacc);
for k = 1:numel(saccIdx)
    saccLat = EEG_ep.event(saccIdx(k)).latency;
    if isempty(blinkLat) || ~any(abs(saccLat - blinkLat) <= 50)
        nSacc = nSacc + 1;
    end
end
end

function pct = compute_pct_baseline(stim, baseline)
% Percentage change: 100*(stim-baseline)/baseline. Non-positive baselines -> NaN.
pct = 100 * (stim - baseline) ./ baseline;
pct(~isfinite(stim) | ~isfinite(baseline) | ~isfinite(pct) | baseline <= 0) = NaN;
end
