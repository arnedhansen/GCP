%% Preprocessing GCP

%% Setup
startup
[subjects, paths] = setup('GCP');
mergedPath = paths.merged;
runMode = askRunMode();

% EyeLink event-rate windows (must match GCP_gaze_fex)
baseline_window = [-1.5 -0.5];
analysis_full   = [0 2];
epoch_window    = [-2 3.5];

et_channels = {'L-GAZE-X' 'L-GAZE-Y' 'L-AREA' 'R-GAZE-X' 'R-GAZE-Y' 'R-AREA'};
eeg_channels = {'all' '-B*' '-HEOGR' '-HEOGL' '-VEOGU' '-VEOGL' ...
    '-L-GAZE-X' '-L-GAZE-Y' '-L-AREA' '-R-GAZE-X' '-R-GAZE-Y' '-R-AREA'};

%% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    clearvars -except subjects subj mergedPath paths runMode ...
        baseline_window analysis_full epoch_window et_channels eeg_channels
    clc; fprintf('[PREP] Subject %s (%d/%d)\n', subjects{subj}, subj, length(subjects));
    datapath = fullfile(mergedPath, subjects{subj});
    cd(datapath)

    eegFile = fullfile(paths.features, subjects{subj}, 'eeg', 'dataEEG.mat');
    processThis = strcmp(runMode, 'all') || isempty(dir(eegFile));
    if ~processThis
        fprintf('[PREP] Subject %s already has dataEEG.mat. Skipping (NEW mode).\n', subjects{subj});
        continue
    end

    %% Read blocks
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
    dataEEG_c25 = cell(1, 4);
    dataEEG_c50 = cell(1, 4);
    dataEEG_c75 = cell(1, 4);
    dataEEG_c100 = cell(1, 4);
    dataET_c25 = cell(1, 4);
    dataET_c50 = cell(1, 4);
    dataET_c75 = cell(1, 4);
    dataET_c100 = cell(1, 4);
    eeglab_c25 = cell(1, 4);
    eeglab_c50 = cell(1, 4);
    eeglab_c75 = cell(1, 4);
    eeglab_c100 = cell(1, 4);

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
            data_c25{block} = []; data_c50{block} = [];
            data_c75{block} = []; data_c100{block} = [];
            continue
        end

        [data_c25{block}, data_c50{block}, data_c75{block}, data_c100{block}, ...
            keep25, keep50, keep75, keep100] = apply_fixation_mask_block( ...
            data_c25{block}, data_c50{block}, data_c75{block}, data_c100{block}, ...
            paths.raw_occ, subjects{subj}, block);

        EEG_c25 = select_eeglab_trials(EEG_c25, keep25);
        EEG_c50 = select_eeglab_trials(EEG_c50, keep50);
        EEG_c75 = select_eeglab_trials(EEG_c75, keep75);
        EEG_c100 = select_eeglab_trials(EEG_c100, keep100);

        cfg_et = [];
        cfg_et.channel = et_channels;
        cfg_eeg = [];
        cfg_eeg.channel = eeg_channels;
        cfg_ref = [];
        cfg_ref.reref = 'yes';
        cfg_ref.refchannel = 'all';

        dataET_b = {[], [], [], []};
        dataEEG_b = {[], [], [], []};
        EEG_b = {EEG_c25, EEG_c50, EEG_c75, EEG_c100};
        data_b = {data_c25{block}, data_c50{block}, data_c75{block}, data_c100{block}};
        for ci = 1:4
            if isempty(data_b{ci}), continue; end
            dataET_b{ci} = ft_selectdata(cfg_et, data_b{ci});
            dataEEG_b{ci} = ft_selectdata(cfg_eeg, data_b{ci});
            dataEEG_b{ci} = ft_preprocessing(cfg_ref, dataEEG_b{ci});
        end

        n_initial = n_ft_trials(dataEEG_b{1}) + n_ft_trials(dataEEG_b{2}) + ...
            n_ft_trials(dataEEG_b{3}) + n_ft_trials(dataEEG_b{4});
        [dataEEG_b, keep90] = reject_90uv_block(dataEEG_b, n_initial, subjects{subj}, block);
        for ci = 1:4
            dataET_b{ci} = select_ft_trials(dataET_b{ci}, keep90{ci});
            EEG_b{ci} = select_eeglab_trials(EEG_b{ci}, keep90{ci});
        end

        dataEEG_c25{block} = dataEEG_b{1};
        dataEEG_c50{block} = dataEEG_b{2};
        dataEEG_c75{block} = dataEEG_b{3};
        dataEEG_c100{block} = dataEEG_b{4};
        dataET_c25{block} = dataET_b{1};
        dataET_c50{block} = dataET_b{2};
        dataET_c75{block} = dataET_b{3};
        dataET_c100{block} = dataET_b{4};
        eeglab_c25{block} = EEG_b{1};
        eeglab_c50{block} = EEG_b{2};
        eeglab_c75{block} = EEG_b{3};
        eeglab_c100{block} = EEG_b{4};
    end

    %% Remove empty blocks
    dataEEG_c25 = dataEEG_c25(~cellfun('isempty', dataEEG_c25));
    dataEEG_c50 = dataEEG_c50(~cellfun('isempty', dataEEG_c50));
    dataEEG_c75 = dataEEG_c75(~cellfun('isempty', dataEEG_c75));
    dataEEG_c100 = dataEEG_c100(~cellfun('isempty', dataEEG_c100));
    dataET_c25 = dataET_c25(~cellfun('isempty', dataET_c25));
    dataET_c50 = dataET_c50(~cellfun('isempty', dataET_c50));
    dataET_c75 = dataET_c75(~cellfun('isempty', dataET_c75));
    dataET_c100 = dataET_c100(~cellfun('isempty', dataET_c100));
    eeglab_c25 = eeglab_c25(~cellfun('isempty', eeglab_c25));
    eeglab_c50 = eeglab_c50(~cellfun('isempty', eeglab_c50));
    eeglab_c75 = eeglab_c75(~cellfun('isempty', eeglab_c75));
    eeglab_c100 = eeglab_c100(~cellfun('isempty', eeglab_c100));

    if isempty(dataEEG_c25) || isempty(dataEEG_c50) || isempty(dataEEG_c75) || isempty(dataEEG_c100)
        fprintf(['[PREP] Subject %s: at least one condition has no epochs ' ...
            'after fixation / 90 uV rejection. Skipping EEG/ET save.\n'], subjects{subj});
    else
        %% Equalize labels
        dataEEG_c25 = update_labels(dataEEG_c25);
        dataEEG_c50 = update_labels(dataEEG_c50);
        dataEEG_c75 = update_labels(dataEEG_c75);
        dataEEG_c100 = update_labels(dataEEG_c100);
        dataET_c25 = update_labels(dataET_c25);
        dataET_c50 = update_labels(dataET_c50);
        dataET_c75 = update_labels(dataET_c75);
        dataET_c100 = update_labels(dataET_c100);

        %% Add trialinfo
        dataEEG_c25 = set_trialinfo_cells(dataEEG_c25, 61, 'c25');
        dataEEG_c50 = set_trialinfo_cells(dataEEG_c50, 62, 'c50');
        dataEEG_c75 = set_trialinfo_cells(dataEEG_c75, 63, 'c75');
        dataEEG_c100 = set_trialinfo_cells(dataEEG_c100, 64, 'c100');
        dataET_c25 = set_trialinfo_cells(dataET_c25, 61, 'c25');
        dataET_c50 = set_trialinfo_cells(dataET_c50, 62, 'c50');
        dataET_c75 = set_trialinfo_cells(dataET_c75, 63, 'c75');
        dataET_c100 = set_trialinfo_cells(dataET_c100, 64, 'c100');

        %% Append data for conditions
        cfg = [];
        cfg.keepsampleinfo = 'yes';
        dataEEG_c25 = ft_appenddata(cfg, dataEEG_c25{:});
        dataEEG_c50 = ft_appenddata(cfg, dataEEG_c50{:});
        dataEEG_c75 = ft_appenddata(cfg, dataEEG_c75{:});
        dataEEG_c100 = ft_appenddata(cfg, dataEEG_c100{:});
        dataET_c25 = ft_appenddata(cfg, dataET_c25{:});
        dataET_c50 = ft_appenddata(cfg, dataET_c50{:});
        dataET_c75 = ft_appenddata(cfg, dataET_c75{:});
        dataET_c100 = ft_appenddata(cfg, dataET_c100{:});

        %% Save EEG / ET
        savepath = fullfile(paths.features, subjects{subj}, 'eeg');
        mkdir(savepath)
        cd(savepath)
        save dataEEG dataEEG_c25 dataEEG_c50 dataEEG_c75 dataEEG_c100
        savepathET = fullfile(paths.features, subjects{subj}, 'gaze');
        mkdir(savepathET)
        cd(savepathET)
        save dataET dataET_c25 dataET_c50 dataET_c75 dataET_c100

        %% EyeLink event rates from remaining (fixation- and 90 uV-filtered) epochs
        ev = eyelink_event_rates_from_eeglab( ...
            eeglab_c25, eeglab_c50, eeglab_c75, eeglab_c100, analysis_full, baseline_window);

        c25_blinks = ev.blinks(1); c50_blinks = ev.blinks(2);
        c75_blinks = ev.blinks(3); c100_blinks = ev.blinks(4);
        c25_fixations = ev.fixations(1); c50_fixations = ev.fixations(2);
        c75_fixations = ev.fixations(3); c100_fixations = ev.fixations(4);
        c25_saccades = ev.saccades(1); c50_saccades = ev.saccades(2);
        c75_saccades = ev.saccades(3); c100_saccades = ev.saccades(4);

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

        c25_pct_fixations = compute_pct_baseline(c25_fixations, c25_bl_fixations);
        c50_pct_fixations = compute_pct_baseline(c50_fixations, c50_bl_fixations);
        c75_pct_fixations = compute_pct_baseline(c75_fixations, c75_bl_fixations);
        c100_pct_fixations = compute_pct_baseline(c100_fixations, c100_bl_fixations);

        c25_pct_saccades = compute_pct_baseline(c25_saccades, c25_bl_saccades);
        c50_pct_saccades = compute_pct_baseline(c50_saccades, c50_bl_saccades);
        c75_pct_saccades = compute_pct_baseline(c75_saccades, c75_bl_saccades);
        c100_pct_saccades = compute_pct_baseline(c100_saccades, c100_bl_saccades);

        save(fullfile(savepathET, 'gaze_metrics'), ...
            'c25_blinks', 'c50_blinks', 'c75_blinks', 'c100_blinks', ...
            'c25_fixations', 'c50_fixations', 'c75_fixations', 'c100_fixations', ...
            'c25_saccades', 'c50_saccades', 'c75_saccades', 'c100_saccades', ...
            'c25_bl_blinks', 'c50_bl_blinks', 'c75_bl_blinks', 'c100_bl_blinks', ...
            'c25_bl_fixations', 'c50_bl_fixations', 'c75_bl_fixations', 'c100_bl_fixations', ...
            'c25_bl_saccades', 'c50_bl_saccades', 'c75_bl_saccades', 'c100_bl_saccades', ...
            'c25_pct_blinks', 'c50_pct_blinks', 'c75_pct_blinks', 'c100_pct_blinks', ...
            'c25_pct_fixations', 'c50_pct_fixations', 'c75_pct_fixations', 'c100_pct_fixations', ...
            'c25_pct_saccades', 'c50_pct_saccades', 'c75_pct_saccades', 'c100_pct_saccades', ...
            'analysis_full', 'baseline_window');
    end

    clc
    if subj == length(subjects)
        fprintf('[PREP] Subject %s (%d/%d) done. PREPROCESSING FINALIZED.\n', subjects{subj}, subj, length(subjects))
    else
        fprintf('[PREP] Subject %s (%d/%d) done. Loading next subject...\n', subjects{subj}, subj, length(subjects))
    end
end

function n = n_ft_trials(data)
if isempty(data) || ~isstruct(data) || ~isfield(data, 'trial')
    n = 0;
else
    n = numel(data.trial);
end
end

function data = select_ft_trials(data, keep)
if isempty(data) || isempty(keep)
    data = [];
    return
end
keep = keep(:) ~= 0;
if ~any(keep)
    data = [];
    return
end
if all(keep)
    return
end
cfg = [];
cfg.trials = find(keep);
data = ft_selectdata(cfg, data);
end

function EEG = select_eeglab_trials(EEG, keep)
if isempty(EEG) || isempty(keep)
    EEG = [];
    return
end
keep = keep(:) ~= 0;
if ~any(keep)
    EEG = [];
    return
end
if all(keep)
    return
end
EEG = pop_select(EEG, 'trial', find(keep));
end

function cells = set_trialinfo_cells(cells, code, tag)
for block = 1:numel(cells)
    try
        cells{block}.trialinfo = zeros(numel(cells{block}.trial), 1) + code;
    catch ME
        ME.message
        fprintf('[PREP] ERROR adding trialinfo (%s) in Block %d!\n', tag, block)
    end
end
end

function [d25, d50, d75, d100, k25, k50, k75, k100] = apply_fixation_mask_block(d25, d50, d75, d100, raw_occ, subj, block)
% Among non-catch trials, drop fixation == 0; keep fixation == 1 (including replacements).
behav_file = fullfile(raw_occ, subj, sprintf('%s_GCP_block%d.mat', subj, block));
if ~isfile(behav_file)
    error('GCP_preprocessing:MissingBehav', ...
        'Behavioural file not found for fixation mask: %s', behav_file);
end
S = load(behav_file);
if ~isfield(S, 'saves') || ~isfield(S.saves, 'data')
    error('GCP_preprocessing:BehavFormat', 'No saves.data in %s', behav_file);
end
D = S.saves.data;
need = {'whiteCross', 'fixation', 'grating'};
for k = 1:numel(need)
    if ~isfield(D, need{k})
        error('GCP_preprocessing:BehavField', 'Missing %s in %s', need{k}, behav_file);
    end
end
whiteCross = D.whiteCross(:);
fixation = D.fixation(:);
grating = D.grating(:);
if numel(whiteCross) ~= numel(fixation) || numel(fixation) ~= numel(grating)
    error('GCP_preprocessing:BehavLength', ...
        'whiteCross/fixation/grating length mismatch in %s', behav_file);
end
noncatch = whiteCross == 0;
fix_nc = fixation(noncatch);
grat_nc = grating(noncatch);

n_ep = n_ft_trials(d25) + n_ft_trials(d50) + n_ft_trials(d75) + n_ft_trials(d100);
n_nc = numel(fix_nc);
if n_ep ~= n_nc
    error('GCP_preprocessing:EpochBehavMismatch', ...
        ['Subject %s block %d: %d epochs (61-64) vs %d non-catch behavioural trials. ' ...
        'Refusing silent misalignment.'], subj, block, n_ep, n_nc);
end

[d25, k25] = mask_fixation_condition(d25, fix_nc, grat_nc, 1, subj, block, 61);
[d50, k50] = mask_fixation_condition(d50, fix_nc, grat_nc, 2, subj, block, 62);
[d75, k75] = mask_fixation_condition(d75, fix_nc, grat_nc, 3, subj, block, 63);
[d100, k100] = mask_fixation_condition(d100, fix_nc, grat_nc, 4, subj, block, 64);
end

function [data, keep_behav] = mask_fixation_condition(data, fix_nc, grat_nc, cond_code, subj, block, trig)
keep_behav = fix_nc(grat_nc == cond_code) == 1;
n_ep = n_ft_trials(data);
if n_ep ~= numel(keep_behav)
    error('GCP_preprocessing:CondMismatch', ...
        ['Subject %s block %d trigger %d: %d epochs vs %d non-catch trials of that contrast. ' ...
        'Refusing silent misalignment.'], subj, block, trig, n_ep, numel(keep_behav));
end
if n_ep == 0 || ~any(keep_behav)
    data = [];
    return
end
if all(keep_behav)
    return
end
cfg = [];
cfg.trials = find(keep_behav);
data = ft_selectdata(cfg, data);
end

function [sets, keep] = reject_90uv_block(sets, n_initial, subj, block)
% Flag after average rereference. Denominator is n_initial (post-fixation
% epochs entering 90 uV), not the count after flagged trials are dropped.
keep = {[], [], [], []};
n_flagged = 0;
for i = 1:4
    keep{i} = eeg_90uv_keep_mask(sets{i});
    if isempty(keep{i})
        n_flagged = n_flagged + n_ft_trials(sets{i});
    else
        n_flagged = n_flagged + sum(~keep{i});
    end
end
if n_initial == 0
    return
end
if (n_flagged / n_initial) > 0.5
    fprintf(['[PREP] Subject %s block %d: %d/%d initial epochs exceed +/-90 uV ' ...
        '(>50%%). Discarding block.\n'], subj, block, n_flagged, n_initial);
    sets = {[], [], [], []};
    keep = {false(0, 1), false(0, 1), false(0, 1), false(0, 1)};
    return
end
for i = 1:4
    sets{i} = select_ft_trials(sets{i}, keep{i});
end
if n_flagged > 0
    fprintf('[PREP] Subject %s block %d: dropped %d/%d initial epochs exceeding +/-90 uV.\n', ...
        subj, block, n_flagged, n_initial);
end
end

function keep = eeg_90uv_keep_mask(data)
n_t = n_ft_trials(data);
keep = true(n_t, 1);
if n_t == 0
    keep = false(0, 1);
    return
end
cfg = [];
cfg.continuous = 'no';
cfg.feedback = 'no';
cfg.artfctdef.threshold.channel = data.label;
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.min = -90;
cfg.artfctdef.threshold.max = 90;
[~, artifact] = ft_artifact_threshold(cfg, data);
if isempty(artifact)
    return
end
cfg_r = [];
cfg_r.artfctdef.reject = 'complete';
cfg_r.artfctdef.threshold.artifact = artifact;
d_keep = ft_rejectartifact(cfg_r, data);
n_kept = n_ft_trials(d_keep);
if n_kept == 0
    keep = false(n_t, 1);
    return
end
if isfield(data, 'sampleinfo') && isfield(d_keep, 'sampleinfo') && ~isempty(d_keep.sampleinfo)
    keep = ismember(data.sampleinfo, d_keep.sampleinfo, 'rows');
else
    keep = true(n_t, 1);
    keep((n_kept + 1):end) = false;
end
end

function ev = eyelink_event_rates_from_eeglab(eeg_c25, eeg_c50, eeg_c75, eeg_c100, stimWin, baselineWin)
% Rates [Hz] from remaining EEGLAB epochs (same trials as saved dataET).
sets = {eeg_c25, eeg_c50, eeg_c75, eeg_c100};
ev = struct();
ev.blinks = nan(4, 1);
ev.fixations = nan(4, 1);
ev.saccades = nan(4, 1);
ev.bl_blinks = nan(4, 1);
ev.bl_fixations = nan(4, 1);
ev.bl_saccades = nan(4, 1);
durStim = diff(stimWin);
durBl = diff(baselineWin);
for c = 1:4
    nTrials = 0;
    countsStim = [0 0 0];
    countsBl = [0 0 0];
    blocks = sets{c};
    for bi = 1:numel(blocks)
        EEG = blocks{bi};
        if isempty(EEG) || ~isfield(EEG, 'event') || EEG.trials < 1
            continue
        end
        nTrials = nTrials + EEG.trials;
        countsStim = countsStim + count_eyelink_events_window(EEG, stimWin);
        countsBl = countsBl + count_eyelink_events_window(EEG, baselineWin);
    end
    if nTrials > 0 && durStim > 0
        rates = countsStim ./ (nTrials * durStim);
        ev.blinks(c) = rates(1);
        ev.fixations(c) = rates(2);
        ev.saccades(c) = rates(3);
    end
    if nTrials > 0 && durBl > 0
        ratesBl = countsBl ./ (nTrials * durBl);
        ev.bl_blinks(c) = ratesBl(1);
        ev.bl_fixations(c) = ratesBl(2);
        ev.bl_saccades(c) = ratesBl(3);
    end
end
end

function counts = count_eyelink_events_window(EEG, tWin)
counts = [0 0 0];
if isempty(EEG) || ~isfield(EEG, 'event') || isempty(EEG.event)
    return
end
blink_t = [];
sacc_t = [];
n_sacc = 0;
for k = 1:numel(EEG.event)
    t_ev = eeglab_event_time_s(EEG, k);
    if ~isfinite(t_ev) || t_ev < tWin(1) || t_ev > tWin(2)
        continue
    end
    typ = eeglab_event_type(EEG.event(k).type);
    if strcmp(typ, 'L_blink') || strcmp(typ, 'R_blink')
        counts(1) = counts(1) + 1;
        blink_t(end+1) = t_ev; %#ok<AGROW>
    elseif strcmp(typ, 'L_fixation') || strcmp(typ, 'R_fixation')
        counts(2) = counts(2) + 1;
    elseif strcmp(typ, 'L_saccade') || strcmp(typ, 'R_saccade')
        n_sacc = n_sacc + 1;
        sacc_t(end+1) = t_ev; %#ok<AGROW>
    end
end
nSaccKeep = 0;
for k = 1:numel(sacc_t)
    if isempty(blink_t) || ~any(abs(sacc_t(k) - blink_t) <= 0.1)
        nSaccKeep = nSaccKeep + 1;
    end
end
counts(3) = nSaccKeep;
end

function t_ev = eeglab_event_time_s(EEG, evIdx)
lat = EEG.event(evIdx).latency;
if isfield(EEG.event, 'epoch') && ~isempty(EEG.event(evIdx).epoch)
    ep = EEG.event(evIdx).epoch;
else
    ep = 1;
end
if lat >= 1 && lat <= EEG.pnts
    samp = lat;
else
    samp = lat - (ep - 1) * EEG.pnts;
end
t_ev = EEG.xmin + (samp - 1) / EEG.srate;
end

function typ = eeglab_event_type(t)
if ischar(t)
    typ = t;
elseif isstring(t)
    typ = char(t);
elseif iscell(t) && ~isempty(t)
    typ = char(string(t{1}));
else
    typ = '';
end
end

function pct = compute_pct_baseline(stim, baseline)
% Percentage change: 100*(stim-baseline)/baseline. Non-positive baselines -> NaN.
pct = 100 * (stim - baseline) ./ baseline;
pct(~isfinite(stim) | ~isfinite(baseline) | ~isfinite(pct) | baseline <= 0) = NaN;
end
