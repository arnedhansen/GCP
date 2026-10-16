%% GCP microsaccade contrast test
% Standalone grating task for testing contrast effects on microsaccades.
% Eye tracking only. No EEG, no master script.
%
% Run this file directly. It sets paths, opens the screen, calibrates the
% EyeLink, and saves an EDF plus a .mat file.
%
% Contrasts: 12, 25, 50, 100% Michelson.
%   Lower contrast raises the average microsaccade rate because
%   inhibition is weaker and later (Bonneh et al., 2015, JOV;
%   Rolfs et al., 2008, JOV). 4% was dropped: on this large grating it
%   can fall near threshold and add search saccades. 12% should remain
%   clearly visible while still sitting below the flat 25 to 100% range
%   seen in the GCP pilots and in Kuo and Juan (2026).
%
% This script uses constant mean luminance (true Michelson contrast) and
% the same spatial envelope at every level. GCP_gratingsTask does not:
%   its formula multiplies the whole image by the contrast factor, so mean
%   luminance falls as nominal contrast falls, and the 50% grating is the
%   only one without the taper. Those differences will move microsaccades
%   on their own. Results here are therefore a contrast test, not a
%   pixel matched replay of the main task.
%
% Viewing is passive. The white cross button task is omitted because that
% task is itself a strong microsaccade modulator.
%
% Playback of the drift matches GCP_gratingsTask (one precomputed texture
% cycle spread across the 2 s presentation). The header of that script
% states 2 Hz; the playback loop does not advance one texture per refresh.

%% General settings, screen and paths
clear all;
close all;
clc;

DATA_PATH = '/home/methlab/Desktop/GCP_data';
FUNS_PATH = '/home/methlab/Desktop/GCP';

addpath(FUNS_PATH);
addpath(fileparts(mfilename('fullpath'))); % helpers next to this script
screenSettings;

subject.ID = 991; % Set to 99X for tests

% Variables required by calibrateET (EDF name: <ID>_5.edf, not blocks 1-4)
TASK     = 'GCP';
TRAINING = 0;
BLOCK    = 5;

ListenChar(2);

try

%% Contrast and trial settings
contrastLevels = [0.12, 0.25, 0.50, 1.00]; % Michelson, constant mean
exp.nPerContrast = 150; % one observer: ~2 SE for a 0.2 Hz low vs high rate difference
exp.nContrasts   = numel(contrastLevels);
exp.nTrials      = exp.nPerContrast * exp.nContrasts;

enableScreenshots = 0; % Writing images inside the flip loop destroys timing

%% EyeLink
disp('CALIBRATING ET...');
calibrateET;

HideCursor(whichScreen);

%% EyeLink messages (no EEG triggers)
TASK_START         = 10;
BLOCK_START        = 15;
FIXCROSS           = 17;
FIXATION_REMINDER  = 91;
BLOCK_END          = 75;
TASK_END           = 90;
% Presentation messages are 81..84 in contrastLevels order (12, 25, 50, 100%)

%% Instructions
startExperimentText = [
    'You will see a series of gratings. \n\n' ...
    'Please always look at the center of the screen. \n\n' ...
    'No button press is required. \n\n' ...
    '\n\n' ...
    'Press any key to continue...'];
loadingText = 'Loading CONTRAST TEST...';

%% Psychtoolbox
global GL;
AssertOpenGL;

global ptb_drawformattedtext_disableClipping;
ptb_drawformattedtext_disableClipping = 1;

Screen('Preference', 'Verbosity', 0);
KbName('UnifyKeyNames');

%% Window
screen.ID = whichScreen;

if gray == white
    gray = white / 2;
end

inc = white - gray;
backgroundColor = gray;

[ptbWindow, winRect] = Screen('OpenWindow', screen.ID, backgroundColor);

[screen.centerX, screen.centerY] = RectCenter(winRect);
screen.width       = 48;
screen.height      = 29.89;
screen.resolutionX = 800;
screen.resolutionY = 600;
screen.viewDist    = 80;

screen.totVisDeg = 2 * atan(screen.width / (2 * screen.viewDist)) * (180 / pi);
screen.ppd       = screen.resolutionX / screen.totVisDeg;
% Pixel sizes matched to the prior ppd = 50 setup (500 px grating, 18 px cross)
gratingSizePixTarget  = 500;
fixationSizePixTarget = 18;

ifi       = Screen('GetFlipInterval', ptbWindow);
frameRate = Screen('FrameRate', screen.ID);

Screen('BlendFunction', ptbWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

Screen('TextSize', ptbWindow, 25);
DrawFormattedText(ptbWindow, loadingText, 'center', 'center', black);
Screen('Flip', ptbWindow);

%% Fixation cross (black only)
fixationSize_dva  = fixationSizePixTarget / screen.ppd;
fixationLineWidth = 1.5;
fixationColor     = [0 0 0];

fixationSize_pix = round(fixationSize_dva * screen.ppd);
fixHorizontal    = [round(-fixationSize_pix / 2) round(fixationSize_pix / 2) 0 0];
fixVertical      = [0 0 round(-fixationSize_pix / 2) round(fixationSize_pix / 2)];
fixCoords        = [fixHorizontal; fixVertical];

timing.cfilower         = 2000;
timing.cfiupper         = 3000;
timing.fixCheckDuration = 0.5;

%% Grating geometry (same size, spatial frequency and duration as the main task)
visualAngleGrating = gratingSizePixTarget / screen.ppd;
gratingSize        = visualAngleGrating * screen.ppd;
gratingRadius      = round(gratingSize / 2);
gratingSize        = 2 * gratingRadius;

targetCPD      = 3;
driftFreq      = 2; % Nominal Hz. Playback below follows GCP_gratingsTask, not one frame per refresh.
nFramesInCycle = round((1 / driftFreq) / ifi);

movieDurationSecs = 2;

gratingDim      = [0 0 2 * gratingRadius 2 * gratingRadius];
gratingPosition = CenterRectOnPointd(gratingDim, screen.centerX, screen.centerY);

priorityLevel = MaxPriority(ptbWindow);
Priority(priorityLevel);

%% Textures: Michelson contrast, mean held at gray, same envelope for every level
[x, y] = meshgrid(-gratingRadius:gratingRadius, -gratingRadius:gratingRadius);
f = screen.ppd / (2 * pi * targetCPD);

L = 2 * gratingRadius + 1;
w1D = hann(L);
xx  = linspace(-gratingRadius, gratingRadius, L);
[X, Y] = meshgrid(xx);
r   = sqrt(X.^2 + Y.^2);
w2D = zeros(L);
w2D(r <= gratingRadius) = interp1(xx, w1D, r(r <= gratingRadius));

radius    = sqrt(x.^2 + y.^2);
maxRadius = gratingRadius;
taperStart = maxRadius * 0.25;
taperMask  = 0.5 * (1 + cos(pi * (radius - taperStart) / (maxRadius - taperStart)));
taperMask(radius <= taperStart) = 1;
taperMask(radius > maxRadius)   = 0;

tex = zeros(exp.nContrasts, nFramesInCycle);
for iContrast = 1:exp.nContrasts
    for jFrame = 1:nFramesInCycle
        phase = (jFrame / nFramesInCycle) * 2 * pi;
        m = sin(sqrt(x.^2 + y.^2) / f + phase);
        % gray +/- contrast * windowed sinusoid. Mean stays at gray.
        grating = gray + contrastLevels(iContrast) * (w2D .* (inc * m));
        grating = grating .* taperMask + gray * (1 - taperMask);
        grating = min(max(grating, black), white);
        tex(iContrast, jFrame) = Screen('MakeTexture', ptbWindow, grating);
    end
end

%% Trial sequence (balanced, no run longer than 2)
seq = repmat(1:exp.nContrasts, 1, exp.nPerContrast);
gratingSequence = seq(randperm(numel(seq)));
maxRun = 3;
for iShuffle = 1:1000
    runLen = 1;
    maxRun = 1;
    for k = 2:numel(gratingSequence)
        if gratingSequence(k) == gratingSequence(k - 1)
            runLen = runLen + 1;
            if runLen > maxRun
                maxRun = runLen;
            end
        else
            runLen = 1;
        end
    end
    if maxRun <= 2
        break;
    end
    gratingSequence = seq(randperm(numel(seq)));
end

data.contrast     = nan(1, exp.nTrials);
data.contrastPct  = nan(1, exp.nTrials);
data.fixation     = nan(1, exp.nTrials);
data.cfi          = nan(1, exp.nTrials);
data.trlDuration  = nan(1, exp.nTrials);

%% Start
DrawFormattedText(ptbWindow, startExperimentText, 'center', 'center', black);
Screen('Flip', ptbWindow);
disp(upper('Participant is reading the instructions.'));
KbReleaseWait;
KbWait(-1, 2);

Eyelink('Message', num2str(TASK_START));
Eyelink('command', 'record_status_message "MS CONTRAST START"');
Eyelink('Message', num2str(BLOCK_START));
Eyelink('command', 'record_status_message "MS CONTRAST BLOCK"');

HideCursor(whichScreen);
timing.startTime = datestr(now, 'dd/mm/yy-HH:MM:SS');
breakAfter = floor(exp.nTrials / 2);

clc;
disp('GCP MICROSACCADE CONTRAST TEST...');
disp(['Contrasts: ' sprintf('%d%% ', round(contrastLevels * 100))]);

%% Trial loop
for trl = 1:exp.nTrials
    tic;
    iContrast = gratingSequence(trl);
    data.contrast(trl)    = contrastLevels(iContrast);
    data.contrastPct(trl) = round(contrastLevels(iContrast) * 100);

    if trl == breakAfter + 1
        DrawFormattedText(ptbWindow, ...
            'Break. Rest for a moment. \n\n Press any key to continue.', ...
            'center', 'center', black);
        Screen('Flip', ptbWindow);
        KbReleaseWait;
        KbWait(-1, 2);
    end

    %% Fixation
    Screen('FillRect', ptbWindow, backgroundColor);
    Screen('Flip', ptbWindow);
    timing.cfi(trl) = (randsample(timing.cfilower:timing.cfiupper, 1)) / 1000;
    data.cfi(trl) = timing.cfi(trl);

    Screen('DrawLines', ptbWindow, fixCoords, fixationLineWidth, fixationColor, ...
        [screen.centerX screen.centerY], 2);
    Screen('Flip', ptbWindow);
    Eyelink('Message', num2str(FIXCROSS));
    Eyelink('command', 'record_status_message "FIXCROSS"');

    WaitSecs(timing.cfi(trl) - timing.fixCheckDuration);
    noFixation = checkFixation(screen.centerX, screen.centerY, timing.fixCheckDuration);
    data.fixation(trl) = ~noFixation;

    if noFixation
        DrawFormattedText(ptbWindow, 'Please fixate on the cross!', 'center', 'center', black);
        Screen('Flip', ptbWindow);
        Eyelink('Message', num2str(FIXATION_REMINDER));
        Eyelink('command', 'record_status_message "FIX_REMINDER"');
        WaitSecs(1);
    end

    %% Grating (same playback timing as GCP_gratingsTask)
    Screen('Flip', ptbWindow);
    Eyelink('Message', num2str(80 + iContrast));
    Eyelink('Message', sprintf('CONTRAST %d', data.contrastPct(trl)));
    Eyelink('command', 'record_status_message "PRESENTATION"');

    probeStartTime = GetSecs;
    frameDuration  = movieDurationSecs / nFramesInCycle;
    whileCount     = 1;
    while (GetSecs - probeStartTime) < movieDurationSecs && whileCount <= nFramesInCycle
        Screen('DrawTexture', ptbWindow, tex(iContrast, whileCount), [], gratingPosition);
        Screen('Flip', ptbWindow);
        if enableScreenshots == 1 && whileCount == 1
            screenshot(sprintf('GCP_ms_contrast_%d.png', data.contrastPct(trl)), ptbWindow, 1);
        end

        elapsedTime = GetSecs - probeStartTime;
        waitTime = whileCount * frameDuration - elapsedTime;
        if waitTime > 0
            WaitSecs(waitTime);
        end
        whileCount = whileCount + 1;
    end

    data.trlDuration(trl) = toc;
    if data.fixation(trl)
        fixLabel = 'Y';
    else
        fixLabel = 'N';
    end
    fprintf('Trial %3d/%d | %3d%% contrast | Fix: %s\n', ...
        trl, exp.nTrials, data.contrastPct(trl), fixLabel);
end

nFixFails = sum(data.fixation == 0);
disp(['Block summary: ' num2str(exp.nTrials) ' trials. Fixation failures: ' num2str(nFixFails) '.']);

%% End recording
Screen('Flip', ptbWindow);
Eyelink('Message', num2str(BLOCK_END));
Eyelink('command', 'record_status_message "END BLOCK"');
Eyelink('Message', num2str(TASK_END));
Eyelink('command', 'record_status_message "TASK_END"');

timing.endTime = datestr(now, 'dd/mm/yy-HH:MM:SS');
startTime = datetime(timing.startTime, 'InputFormat', 'dd/MM/yy-HH:mm:ss');
endTime   = datetime(timing.endTime, 'InputFormat', 'dd/MM/yy-HH:mm:ss');
timing.duration = seconds(endTime - startTime);

subjectID = num2str(subject.ID);
filePath  = fullfile(DATA_PATH, subjectID);
if ~exist(filePath, 'dir')
    mkdir(filePath);
end
fileName = [subjectID '_GCP_ms_contrast.mat'];

saves.data            = data;
saves.gratingSequence = gratingSequence;
saves.contrastLevels  = contrastLevels;
saves.experiment      = exp;
saves.screen          = screen;
saves.subjectID       = subjectID;
saves.subject         = subject;
saves.timing          = timing;
saves.ifi             = ifi;
saves.frameRate       = frameRate;
saves.targetCPD       = targetCPD;
saves.driftFreqNominal = driftFreq;
saves.contrastMode    = 'michelson_constant_mean';

trigger.TASK_START        = TASK_START;
trigger.BLOCK_START       = BLOCK_START;
trigger.FIXCROSS          = FIXCROSS;
trigger.FIXATION_REMINDER = FIXATION_REMINDER;
trigger.BLOCK_END         = BLOCK_END;
trigger.TASK_END          = TASK_END;
trigger.PRESENTATION      = 80 + (1:exp.nContrasts);
trigger.contrastPct       = round(contrastLevels * 100);

disp('SAVING DATA...');
save(fullfile(filePath, fileName), 'saves', 'trigger');

fprintf('Stop Recording Track\n');
Eyelink('StopRecording');
Eyelink('CloseFile');
fprintf('Downloading File\n');
EL_DownloadDataFile;
EL_Cleanup;

pathEdf2Asc = '/usr/bin/edf2asc';
if isfile(pathEdf2Asc)
    disp('CONVERTING EDF to ASCII...');
    system([pathEdf2Asc ' "' fullfile(filePath, edfFile) '" -y']);
else
    disp('edf2asc not found. EDF left unconverted.');
end

DrawFormattedText(ptbWindow, ...
    'End of the contrast test. \n\n Thank you. \n\n Press any key to close.', ...
    'center', 'center', black);
Screen('Flip', ptbWindow);
KbReleaseWait;
KbWait(-1, 2);

Priority(0);
Screen('Close');
Screen('CloseAll');
ListenChar(0);
disp('CONTRAST TEST FINISHED');

catch ME
    ListenChar(0);
    Priority(0);
    Screen('CloseAll');
    rethrow(ME);
end
