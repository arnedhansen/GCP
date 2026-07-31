%% GCP pilot recording durations (CW output only; nothing saved)
%
% Pilot durations summary (subjects 601 to 610; behavioral timing.mat fields;
% data on /Volumes/g_psyplafor_methlab_data$/OCC/GCP/data):
%   Per task block: ~15.1 min (range ~15.0 to 15.8)
%   Sum of 4 task blocks: ~60.6 min (60.3 to 61.5)
%   Training: ~1.1 min
%   Resting EEG: ~5.05 min (602/610 empty pnts; 607 ~15 min)
%   Wall clock training start to block 4 end: ~80 min mean (72 to 99;
%     inter block gaps ~18 min mean)
%   Flags: 604 missing block4 task EEG; 602/610 resting pnts=0; 607 resting ~15 min
%
% Timing source: saves.timing.startTime / endTime / duration in
%   {ID}_training.mat and {ID}_GCP_block{N}.mat (see paradigm/GCP_gratingsTask.m)

clear
clc
[subjects, paths] = setup('GCP', 0);

% Raw behavioral mats live on the methlab_data share (setup raw_occ may point elsewhere)
dataPath = fullfile('/Volumes/g_psyplafor_methlab_data$', 'OCC', 'GCP', 'data');
if ~isfolder(dataPath) && isfolder(paths.raw_occ)
    dataPath = paths.raw_occ;
end

dirs = dir(dataPath);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
subjects = subjects(~startsWith(subjects, '_') & ~startsWith(subjects, '.'));
subjects = sort(subjects);

fmt = 'dd/MM/yy-HH:mm:ss';
sessionSec = nan(numel(subjects), 1);
taskSec = nan(numel(subjects), 1);

fprintf('\n[CTRL RECORDING] Recording durations (training start -> last block end)\n');
fprintf('%-6s  %-19s  %-19s  %10s  %10s\n', 'ID', 'Start', 'End', 'Session', 'TaskOnly');
fprintf('%s\n', repmat('-', 1, 72));

for s = 1:numel(subjects)
    sid = subjects{s};
    subjDir = fullfile(dataPath, sid);

    cand = [ ...
        dir(fullfile(subjDir, [sid '_training.mat'])); ...
        dir(fullfile(subjDir, [sid '_GCP_block*.mat']))];
    % Drop EEG/ET companions
    keep = true(numel(cand), 1);
    for k = 1:numel(cand)
        nm = cand(k).name;
        if contains(nm, {'_task_', '_EEG', '_ET'})
            keep(k) = false;
        end
    end
    cand = cand(keep);

    tStart = NaT;
    tEnd = NaT;
    taskSum = 0;
    nTaskBlocks = 0;

    for k = 1:numel(cand)
        f = fullfile(cand(k).folder, cand(k).name);
        try
            S = load(f, 'saves');
        catch
            continue
        end
        if ~isfield(S, 'saves') || ~isfield(S.saves, 'timing')
            continue
        end
        tm = S.saves.timing;
        if ~isfield(tm, 'startTime') || ~isfield(tm, 'endTime')
            continue
        end
        st = datetime(tm.startTime, 'InputFormat', fmt);
        en = datetime(tm.endTime, 'InputFormat', fmt);
        if isnat(tStart) || st < tStart
            tStart = st;
        end
        if isnat(tEnd) || en > tEnd
            tEnd = en;
        end
        if contains(cand(k).name, '_GCP_block') && isfield(tm, 'duration')
            taskSum = taskSum + double(tm.duration);
            nTaskBlocks = nTaskBlocks + 1;
        end
    end

    if isnat(tStart) || isnat(tEnd)
        fprintf('%-6s  %-19s  %-19s  %10s  %10s\n', sid, 'NA', 'NA', 'NA', 'NA');
        continue
    end

    sess = seconds(tEnd - tStart);
    sessionSec(s) = sess;
    taskSec(s) = taskSum;

    fprintf('%-6s  %-19s  %-19s  %7.1f min  %7.1f min (%d blk)\n', ...
        sid, datestr(tStart, 'dd/mm/yy-HH:MM:SS'), datestr(tEnd, 'dd/mm/yy-HH:MM:SS'), ...
        sess/60, taskSum/60, nTaskBlocks);
end

ok = ~isnan(sessionSec);
fprintf('%s\n', repmat('-', 1, 72));
fprintf('Subjects with timing: %d / %d\n', sum(ok), numel(subjects));
fprintf('Mean session (train start to last block end): %.1f min\n', mean(sessionSec(ok))/60);
fprintf('Mean task only (sum of block durations):      %.1f min\n', mean(taskSec(ok))/60);
fprintf('\nTOTAL session time across subjects: %.1f min (%.2f h)\n', ...
    nansum(sessionSec)/60, nansum(sessionSec)/3600);
fprintf('TOTAL task only across subjects:    %.1f min (%.2f h)\n\n', ...
    nansum(taskSec)/60, nansum(taskSec)/3600);
