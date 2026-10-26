function data_upload_MethLab()
% Upload GCP .edf, .asc, and .mat files that are on this computer and missing on methlab_data.
%
% A filename counts as already uploaded when it is in the subject folder or in
% archive on the server. It is not copied a second time into the other folder.
% Nothing on methlab_data is overwritten or deleted.
% The share is unmounted before this function returns.

%% Paths
% Share name matches the macOS volume g_psyplafor_methlab_data$.
% Host matches startup.m (smb://idnas37.d.uzh.ch/g_psyplafor_methlab_data$).
% On this Ubuntu machine the same GCP folder is mountPoint/OCC/GCP/data,
% which is /Volumes/g_psyplafor_methlab_data$/OCC/GCP/data on macOS.
shareUNC = '//idnas37.d.uzh.ch/g_psyplafor_methlab_data$';
mountPoint = '/home/methlab/methlab_data';
localRoot = '/home/methlab/Desktop/GCP_data';
serverRoot = fullfile(mountPoint, 'OCC', 'GCP', 'data');

%% State for the unmount, which also runs if the function exits early
shouldUnmount = false;
shareUnmounted = false;
unmountOK = false;
unmountDetail = '';
finished = false;
wb = [];
sudoKeepPid = '';
cleanupUnmount = onCleanup(@() finalizeUnmount()); %#ok<NASGU>
cleanupWaitbar = onCleanup(@() closeWaitbar()); %#ok<NASGU>

try
    if ~(isunix && ~ismac)
        waitfor(msgbox('This script is for the Ubuntu stimulus computer in the lab.'));
        return
    end

    %% make sure the internet conenction is established
    waitfor(msgbox('Make sure the Stimulus Computer is connected to the internet (not to EYELINK), then press ok'));
    connected=0;
    cnt=1;
    while connected==0

        try
            url ='https://www.google.de';
            code=webread(url); %#ok<NASGU>
            connected=1;
        catch
            connected=0;
        end

        if connected==0
            waitfor(msgbox(['Still no connection, try connecting again, then press ok. ('...
                num2str(3-cnt) ' more attempts)' ]));
        else
            break;
        end

        if cnt==3
            break;
        end
        cnt=cnt+1;

    end

    %% Mount methlab_data
    if ~isfolder(mountPoint)
        [mkOK, mkMsg] = mkdir(mountPoint);
        if ~mkOK
            finish('Could not prepare the mount point', { ...
                ['Mount point: ' mountPoint], ...
                mkMsg, ...
                'Nothing was copied.'});
            return
        end
    end

    mountedSource = mountSource(mountPoint);
    if ~isempty(mountedSource)
        if contains(mountedSource, 'g_psyplafor_methlab_data')
            shouldUnmount = true;
            disp('methlab_data is already mounted. It will be unmounted at the end.');
        else
            finish('Mount point is already in use', { ...
                ['Mount point: ' mountPoint], ...
                ['Currently mounted source: ' mountedSource], ...
                'Nothing was copied. That drive was left mounted.'});
            return
        end
    elseif folderHasEntries(mountPoint)
        finish('Mount point is not empty', { ...
            ['Mount point: ' mountPoint], ...
            'The folder is not a network drive and already contains files.', ...
            'Nothing was copied, so those local files stay visible.'});
        return
    else
        mounted = false;
        attempt = 0;
        while attempt < 3
            shortAnswer = inputdlg( ...
                {'Please enter your UZH shortname:'}, ...
                'Connection to methlab_data', [1 50]);
            if isempty(shortAnswer) || isempty(strtrim(shortAnswer{1}))
                finish('Connection cancelled', {'Nothing was copied.'});
                return
            end
            shortname = strtrim(shortAnswer{1});
            if isempty(regexp(shortname, '^[A-Za-z0-9]+$', 'once'))
                waitfor(msgbox('The UZH shortname can only contain letters and digits. Try again.'));
                continue
            end

            attempt = attempt + 1;
            waitfor(msgbox('Enter the Ubuntu password first, then the MethLab server password, in the MATLAB command window.'));
            mountCmd = sprintf(['sudo mount -t cifs %s %s -o "username=%s,domain=UZH,uid=$(id -u),gid=$(id -g),file_mode=0664,dir_mode=0775"'], ...
                shellQuote(shareUNC), shellQuote(mountPoint), shortname);
            [~, mountOut] = system(mountCmd, '-echo');
            mountedSource = mountSource(mountPoint);
            if contains(mountedSource, 'g_psyplafor_methlab_data')
                mounted = true;
                shouldUnmount = true;
                startSudoKeepAlive();
                break
            end
            if attempt < 3
                detail = strtrim(mountOut);
                if isempty(detail)
                    detail = 'No details were returned.';
                end
                waitfor(msgbox(['Connection failed. Check the Ubuntu password, the UZH shortname, and the MethLab server password.' newline newline detail]));
            end
        end
        if ~mounted
            finish('Connection could not be established', { ...
                'The MethLab server was not mounted. Wrong password is a common cause.', ...
                'Nothing was copied. Please run data_upload_MethLab again.'});
            return
        end
    end

    %% GCP .edf, .asc, and .mat files missing on the server
    if ~isfolder(localRoot)
        finish('Local GCP folder not found', { ...
            ['Expected local folder: ' localRoot], ...
            'Nothing was copied.'});
        return
    end
    if ~isfolder(serverRoot)
        finish('Server GCP folder not found', { ...
            ['Expected folder on methlab_data: ' serverRoot], ...
            'This is OCC/GCP/data on g_psyplafor_methlab_data$.', ...
            'Nothing was copied. The folder was not created, so data cannot land in the wrong place.'});
        return
    end

    localIDs = listSubjectFolders(localRoot);
    serverBefore = listSubjectFolders(serverRoot);
    [missingIDs, missingFiles] = missingDataFiles(localRoot, serverRoot, localIDs);

    if isempty(missingIDs)
        finish('Nothing to upload', { ...
            serverRoot, ...
            '', ...
            'Every .edf, .asc, and .mat filename is already on methlab_data.', ...
            '', ...
            sprintf('On methlab_data (%d)', numel(serverBefore)), ...
            formatIDBlock(serverBefore)});
        return
    end

    nFiles = 0;
    for iMissing = 1:numel(missingIDs)
        nFiles = nFiles + numel(missingFiles{iMissing});
    end
    if ~confirmUpload(missingIDs, nFiles, serverRoot)
        finish('Upload cancelled', { ...
            serverRoot, ...
            '', ...
            'Nothing was copied.', ...
            '', ...
            sprintf('On methlab_data (%d)', numel(serverBefore)), ...
            formatIDBlock(serverBefore)});
        return
    end

    %% Copy missing .edf, .asc, and .mat files only
    uploadedIDs = {};
    failedIDs = {};
    failedDetail = {};
    wb = waitbar(0, 'Uploading GCP files. Please wait.');
    for iSub = 1:numel(missingIDs)
        thisID = missingIDs{iSub};
        relFiles = missingFiles{iSub};
        waitbar((iSub - 1) / numel(missingIDs), wb, ...
            sprintf('Uploading %s (%d of %d)', thisID, iSub, numel(missingIDs)));
        srcDir = fullfile(localRoot, thisID);
        dstDir = fullfile(serverRoot, thisID);
        fprintf('Uploading files for %s\n  from %s\n  to   %s\n', thisID, srcDir, dstDir);

        if exist(dstDir, 'file') == 2
            failedIDs{end + 1} = thisID; %#ok<AGROW>
            failedDetail{end + 1} = 'A file already occupies this subject path on methlab_data. It was left unchanged.'; %#ok<AGROW>
            continue
        end
        if ~isfolder(dstDir)
            [mkOK, mkMsg] = mkdir(dstDir);
            if ~mkOK
                failedIDs{end + 1} = thisID; %#ok<AGROW>
                failedDetail{end + 1} = mkMsg; %#ok<AGROW>
                continue
            end
        end

        copiedThis = {};
        failedThis = {};
        for iFile = 1:numel(relFiles)
            fileName = relFiles{iFile};
            if fileNameOnServer(dstDir, fileName)
                continue
            end
            [src, relPath] = localFileSource(srcDir, fileName);
            if ~isAllowedRelativePath(relPath)
                failedThis{end + 1} = [fileName ' was skipped because it is outside the subject folder or archive.']; %#ok<AGROW>
                continue
            end
            dst = fullfile(dstDir, relPath);
            destParent = fileparts(dst);
            if ~isfolder(destParent)
                if exist(destParent, 'file') == 2
                    failedThis{end + 1} = [fileName ' was skipped because a file already occupies the destination folder path.']; %#ok<AGROW>
                    continue
                end
                [mkOK, mkMsg] = mkdir(destParent);
                if ~mkOK
                    failedThis{end + 1} = [fileName ': ' mkMsg]; %#ok<AGROW>
                    continue
                end
            end
            if fileNameOnServer(dstDir, fileName) || exist(dst, 'file') ~= 0
                continue
            end
            if exist(src, 'file') ~= 2
                failedThis{end + 1} = [fileName ' was not found in the local subject folder or archive.']; %#ok<AGROW>
                continue
            end
            copyCmd = ['cp -n -- ' shellQuote(src) ' ' shellQuote(dst)];
            [copyStatus, copyOut] = system(copyCmd, '-echo');
            if copyStatus == 0 && exist(dst, 'file') == 2
                copiedThis{end + 1} = relPath; %#ok<AGROW>
            else
                detail = strtrim(copyOut);
                if isempty(detail)
                    detail = sprintf('cp returned status %d.', copyStatus);
                end
                if fileNameOnServer(dstDir, fileName)
                    detail = [detail ' A partial file may now be on methlab_data. It was left in place. The next run will treat this filename as already present.'];
                end
                failedThis{end + 1} = [fileName ': ' detail]; %#ok<AGROW>
            end
        end
        if ~isempty(copiedThis)
            uploadedIDs{end + 1} = thisID; %#ok<AGROW>
        end
        if ~isempty(failedThis)
            failedIDs{end + 1} = thisID; %#ok<AGROW>
            failedDetail{end + 1} = strjoin(failedThis, ' '); %#ok<AGROW>
        end
    end
    closeWaitbar();

    serverAfter = listSubjectFolders(serverRoot);
    lostIDs = sortSubjectIDs(setdiff(serverBefore, serverAfter, 'stable'));

    report = { ...
        serverRoot, ...
        '', ...
        sprintf('Uploaded (%d)', numel(uploadedIDs)), ...
        formatIDBlock(uploadedIDs)};
    if ~isempty(failedIDs)
        report = [report, {''}, {sprintf('Not fully uploaded (%d):', numel(failedIDs))}]; %#ok<AGROW>
        for iFail = 1:numel(failedIDs)
            report{end + 1} = ['    ' failedIDs{iFail} ': ' failedDetail{iFail}]; %#ok<AGROW>
        end
    end
    if ~isempty(lostIDs)
        report = [report, {''}, { ...
            'WARNING: subject folders that were already on the server are no longer listed.'}, ...
            {'This script has no delete step. Check the share before recording again.'}, ...
            {formatIDBlock(lostIDs)}]; %#ok<AGROW>
    end
    report = [report, {''}, { ...
        sprintf('On methlab_data (%d)', numel(serverAfter))}, ...
        {formatIDBlock(serverAfter)}]; %#ok<AGROW>

    if isempty(failedIDs) && isempty(lostIDs)
        titleText = 'Upload finished';
    else
        titleText = 'Upload finished with problems';
    end
    finish(titleText, report);

catch ME
    finish('Upload stopped', { ...
        'Nothing further was copied after this error:', ...
        ME.message});
end

    function startSudoKeepAlive()
        % The Ubuntu password was just entered for the mount. Refresh that
        % sudo ticket so unmount does not ask for it again.
        pidFile = [tempname '.pid'];
        cmd = sprintf(['nohup bash -c ''echo $$ > %s; while true; do sudo -n -v || exit 0; sleep 50; done'' >/dev/null 2>&1 &'], pidFile);
        system(cmd);
        pause(0.2);
        if isfile(pidFile)
            sudoKeepPid = strtok(strtrim(fileread(pidFile)));
            delete(pidFile);
        end
    end

    function stopSudoKeepAlive()
        if isempty(sudoKeepPid)
            return
        end
        system(['kill ' sudoKeepPid ' >/dev/null 2>&1']);
        sudoKeepPid = '';
    end

    function mountedNow = methlabDataIsMounted()
        mountedNow = contains(mountSource(mountPoint), 'g_psyplafor_methlab_data');
    end

    function reportMountState(mountedNow)
        % Printed after the window closes, so this is the last command window line.
        fprintf('\n');
        if mountedNow
            fprintf('METHLAB_DATA IS STILL MOUNTED\n');
            fprintf('%s\n', strtrim(mountSource(mountPoint)));
        else
            fprintf('METHLAB_DATA IS UNMOUNTED\n');
        end
        fprintf('\n');
    end

    function closeWaitbar()
        if ~isempty(wb) && isgraphics(wb)
            delete(wb);
            wb = [];
        end
    end

    function finalizeUnmount()
        if shareUnmounted || ~shouldUnmount
            shareUnmounted = true;
            return
        end
        shareUnmounted = true;
        if ~contains(mountSource(mountPoint), 'g_psyplafor_methlab_data')
            unmountOK = true;
            unmountDetail = '';
            return
        end
        stopSudoKeepAlive();
        [~, umountOut] = system(['sudo -n umount ' shellQuote(mountPoint)]);
        if ~contains(mountSource(mountPoint), 'g_psyplafor_methlab_data')
            unmountOK = true;
            unmountDetail = '';
        else
            unmountOK = false;
            unmountDetail = strtrim(umountOut);
            if isempty(unmountDetail)
                unmountDetail = 'umount failed. The methlab_data share is still mounted.';
            end
        end
    end

    function finish(titleText, bodyLines)
        if finished
            return
        end
        finished = true;
        closeWaitbar();
        finalizeUnmount();
        mountedNow = methlabDataIsMounted();
        bodyLines = bodyLines(:);
        bodyLines{end + 1} = '';
        if mountedNow
            bodyLines{end + 1} = 'methlab_data is still mounted.';
            if ~isempty(unmountDetail)
                bodyLines{end + 1} = unmountDetail;
            end
        else
            bodyLines{end + 1} = 'methlab_data is unmounted.';
        end
        showTextDialog(titleText, bodyLines);
        reportMountState(mountedNow);
    end
end

function [missingIDs, missingFiles] = missingDataFiles(localRoot, serverRoot, localIDs)
missingIDs = {};
missingFiles = {};
for iID = 1:numel(localIDs)
    thisID = localIDs{iID};
    localFiles = listDataFiles(fullfile(localRoot, thisID));
    if isempty(localFiles)
        continue
    end
    serverFiles = listDataFiles(fullfile(serverRoot, thisID));
    need = sort(setdiff(localFiles, serverFiles, 'stable'));
    if isempty(need)
        continue
    end
    missingIDs{end + 1} = thisID; %#ok<AGROW>
    missingFiles{end + 1} = need; %#ok<AGROW>
end
missingIDs = missingIDs(:);
missingFiles = missingFiles(:);
end

function names = listDataFiles(subjectDir)
names = {};
if ~isfolder(subjectDir)
    return
end
names = listDataFileNamesIn(subjectDir);
archiveDir = fullfile(subjectDir, 'archive');
if isfolder(archiveDir)
    archiveNames = listDataFileNamesIn(archiveDir);
    if isempty(names)
        names = archiveNames(:);
    elseif ~isempty(archiveNames)
        names = [names(:); archiveNames(:)];
    end
end
if isempty(names)
    return
end
names = unique(names);
names = sort(names);
end

function names = listDataFileNamesIn(folderPath)
names = {};
if ~isfolder(folderPath)
    return
end
d = dir(folderPath);
d = d(~[d.isdir]);
names = {d.name};
names = names(endsWith(lower(names), {'.edf', '.asc', '.mat'}));
names = names(:);
end

function tf = fileNameOnServer(subjectDir, fileName)
tf = exist(fullfile(subjectDir, fileName), 'file') == 2 ...
    || exist(fullfile(subjectDir, 'archive', fileName), 'file') == 2;
end

function [src, relPath] = localFileSource(subjectDir, fileName)
archiveFile = fullfile(subjectDir, 'archive', fileName);
rootFile = fullfile(subjectDir, fileName);
if exist(archiveFile, 'file') == 2
    src = archiveFile;
    relPath = fullfile('archive', fileName);
elseif exist(rootFile, 'file') == 2
    src = rootFile;
    relPath = fileName;
else
    src = '';
    relPath = fileName;
end
end

function tf = isAllowedRelativePath(relPath)
[parent, base, ext] = fileparts(relPath);
tf = ~contains(relPath, '..') ...
    && ~isempty(base) ...
    && ismember(lower(ext), {'.edf', '.asc', '.mat'}) ...
    && (isempty(parent) || strcmp(parent, 'archive'));
end

function ids = listSubjectFolders(rootDir)
d = dir(rootDir);
ids = {};
if isempty(d)
    return
end
keep = [d.isdir];
names = {d(keep).name};
junk = {'.', '..', 'archive', '$RECYCLE.BIN', 'System Volume Information', '.Trash-1000', 'lost+found'};
names = names(~ismember(names, junk));
names = names(~startsWith(names, '.'));
ids = sortSubjectIDs(names);
end

function ids = sortSubjectIDs(ids)
ids = ids(:);
if isempty(ids)
    ids = {};
    return
end
nums = str2double(ids);
if all(~isnan(nums))
    [~, order] = sort(nums);
    ids = ids(order);
else
    ids = sort(ids);
end
end

function tf = folderHasEntries(rootDir)
d = dir(rootDir);
names = {d.name};
tf = any(~ismember(names, {'.', '..'}));
end

function src = mountSource(mountPoint)
src = '';
fid = fopen('/proc/mounts', 'r');
if fid < 0
    return
end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
txt = fread(fid, '*char')';
rows = strsplit(txt, newline);
for iRow = 1:numel(rows)
    parts = strsplit(strtrim(rows{iRow}));
    if numel(parts) >= 2 && strcmp(unescapeMountPath(parts{2}), mountPoint)
        src = unescapeMountPath(parts{1});
        return
    end
end
end

function s = unescapeMountPath(s)
s = regexprep(s, '\\([0-7]{3})', '${char(base2dec($1, 8))}');
end

function q = shellQuote(s)
q = ['''' strrep(s, '''', '''\''''') ''''];
end

function block = formatIDBlock(ids)
if isempty(ids)
    block = '    (none)';
    return
end
block = strjoin(strcat({'    '}, ids(:)'), newline);
end

function ok = confirmUpload(subjectIDs, nFiles, serverRoot)
ok = false;
nSubjects = numel(subjectIDs);
figW = 440;
figH = 560;
screenSize = get(0, 'ScreenSize');
figX = max(0, round((screenSize(3) - figW) / 2));
figY = max(0, round((screenSize(4) - figH) / 2));
bg = [0.965 0.965 0.960];
ink = [0.13 0.13 0.13];
muted = [0.45 0.45 0.45];
card = [1 1 1];

f = figure( ...
    'Name', 'GCP', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'Resize', 'off', ...
    'WindowStyle', 'modal', ...
    'Position', [figX figY figW figH], ...
    'Color', bg, ...
    'CloseRequestFcn', @cancelUpload);

uicontrol(f, ...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'Position', [0 figH - 3 figW 3], ...
    'BackgroundColor', ink, ...
    'String', '');
uicontrol(f, ...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'Position', [32 figH - 48 376 14], ...
    'BackgroundColor', bg, ...
    'ForegroundColor', muted, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 9, ...
    'String', 'Study');
uicontrol(f, ...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'Position', [32 figH - 78 376 26], ...
    'BackgroundColor', bg, ...
    'ForegroundColor', ink, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 18, ...
    'String', 'GCP');
uicontrol(f, ...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'Position', [32 figH - 104 376 14], ...
    'BackgroundColor', bg, ...
    'ForegroundColor', muted, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 9, ...
    'String', 'Server folder');
uicontrol(f, ...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'Position', [32 figH - 128 376 20], ...
    'BackgroundColor', bg, ...
    'ForegroundColor', ink, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 11, ...
    'String', serverRoot);
if nSubjects == 1
    subjectWord = 'subject';
else
    subjectWord = 'subjects';
end
if nFiles == 1
    fileWord = 'filename';
else
    fileWord = 'filenames';
end
uicontrol(f, ...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'Position', [32 figH - 188 376 48], ...
    'BackgroundColor', bg, ...
    'ForegroundColor', ink, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 12, ...
    'String', { ...
    sprintf('These %d %s have %d %s', nSubjects, subjectWord, nFiles, fileWord), ...
    'that are not yet on methlab_data.'});
uicontrol(f, ...
    'Style', 'edit', ...
    'Units', 'pixels', ...
    'Position', [32 76 figW - 64 figH - 280], ...
    'BackgroundColor', card, ...
    'ForegroundColor', ink, ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 15, ...
    'Max', 2, ...
    'Min', 0, ...
    'Enable', 'inactive', ...
    'String', strjoin(subjectIDs(:)', newline));
uicontrol(f, ...
    'Style', 'pushbutton', ...
    'String', 'Cancel', ...
    'Units', 'pixels', ...
    'Position', [32 24 100 36], ...
    'BackgroundColor', bg, ...
    'ForegroundColor', muted, ...
    'FontSize', 12, ...
    'Callback', @cancelUpload);
uicontrol(f, ...
    'Style', 'pushbutton', ...
    'String', 'Upload', ...
    'Units', 'pixels', ...
    'Position', [figW - 152 24 120 36], ...
    'BackgroundColor', ink, ...
    'ForegroundColor', card, ...
    'FontSize', 12, ...
    'Callback', @acceptUpload);
uiwait(f);
if isgraphics(f)
    delete(f);
end

    function acceptUpload(~, ~)
        ok = true;
        uiresume(f);
    end

    function cancelUpload(~, ~)
        ok = false;
        uiresume(f);
    end
end

function showTextDialog(titleText, bodyLines)
if ischar(bodyLines) || isstring(bodyLines)
    txt = char(bodyLines);
else
    parts = cellfun(@(c) char(c), bodyLines(:), 'UniformOutput', false);
    txt = strjoin(parts, newline);
end
figW = 440;
figH = 620;
screenSize = get(0, 'ScreenSize');
figX = max(0, round((screenSize(3) - figW) / 2));
figY = max(0, round((screenSize(4) - figH) / 2));
bg = [0.965 0.965 0.960];
ink = [0.13 0.13 0.13];
card = [1 1 1];

f = figure( ...
    'Name', 'GCP', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'Resize', 'off', ...
    'WindowStyle', 'modal', ...
    'Position', [figX figY figW figH], ...
    'Color', bg, ...
    'CloseRequestFcn', @(src, ~) uiresume(src));
uicontrol(f, ...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'Position', [0 figH - 3 figW 3], ...
    'BackgroundColor', ink, ...
    'String', '');
uicontrol(f, ...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'Position', [32 figH - 72 376 36], ...
    'BackgroundColor', bg, ...
    'ForegroundColor', ink, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 16, ...
    'String', titleText);
uicontrol(f, ...
    'Style', 'edit', ...
    'Units', 'pixels', ...
    'Position', [32 76 figW - 64 figH - 168], ...
    'BackgroundColor', card, ...
    'ForegroundColor', ink, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 12, ...
    'Max', 2, ...
    'Min', 0, ...
    'Enable', 'inactive', ...
    'String', txt);
uicontrol(f, ...
    'Style', 'pushbutton', ...
    'String', 'OK', ...
    'Units', 'pixels', ...
    'Position', [figW - 152 24 120 36], ...
    'BackgroundColor', ink, ...
    'ForegroundColor', card, ...
    'FontSize', 12, ...
    'Callback', @(src, ~) uiresume(ancestor(src, 'figure')));
uiwait(f);
if isgraphics(f)
    delete(f);
end
end
