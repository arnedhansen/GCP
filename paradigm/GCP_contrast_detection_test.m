%% GCP Contrast Detection Test
% Screening task run before training. Participants judge whether contrast
% increased or decreased relative to the previous grating (arrow keys).
% Adjacent gratings differ by exactly one contrast step.
% No EEG or eye tracking.
%
% Template: GCP_gratingsTask.m

%% Start message
disp('START OF CONTRAST DETECTION TEST');

%% Hide cursor on participant screen
HideCursor(whichScreen);

%% Set up experiment parameters
exp.nTrials = 21; % First trial unscored; 20 comparisons
exp.passThreshold = 80; % Minimum overall accuracy (%) to continue
contrastPct = [25, 50, 75, 100];

%% Set up text parameters
startExperimentText = [
    'You will see a series of gratings. \n\n' ...
    'After the first grating, press the UP arrow \n\n' ...
    'if contrast increased, or the DOWN arrow \n\n' ...
    'if contrast decreased, compared to the previous grating. \n\n' ...
    'Respond while the grating is on screen. \n\n' ...
    'Always look at the center of the screen. \n\n' ...
    '\n\n' ...
    'Press any key to continue...'];
loadingText = 'Loading CONTRAST DETECTION...';

%% Set up standard Psychtoolbox Settings
global GL;
AssertOpenGL;

global ptb_drawformattedtext_disableClipping;
ptb_drawformattedtext_disableClipping = 1;

Screen('Preference', 'Verbosity', 0);

upKeyCode = KbName('UpArrow');
downKeyCode = KbName('DownArrow');

%% Imaging set up
screen.ID = whichScreen;

if gray == white
    gray = white / 2;
end

inc = white - gray;
backgroundColorGray = 256;

[ptbWindow, winRect] = Screen('OpenWindow', screen.ID, backgroundColorGray);

[screen.centerX, screen.centerY] = RectCenter(winRect);
screen.width       = 48;
screen.height      = 29.89;
screen.resolutionX = 800;
screen.resolutionY = 600;
screen.viewDist    = 80;

screen.totVisDeg = 2 * atan(screen.width / (2 * screen.viewDist)) * (180 / pi);
screen.ppd       = screen.resolutionX / screen.totVisDeg;
screen.ppd       = 50;

ifi       = Screen('GetFlipInterval', ptbWindow);
frameRate = Screen('FrameRate', screen.ID);

Screen('BlendFunction', ptbWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% Text parameters
Screen('TextSize', ptbWindow, 25);

DrawFormattedText(ptbWindow, loadingText, 'center', 'center', black);
Screen('Flip', ptbWindow);

%% Fixation cross parameters (black only)
fixationSize_dva   = .35;
fixationLineWidth  = 1.5;
fixationColorBlack = [0 0 0];

fixationSize_pix = round(fixationSize_dva * screen.ppd);
fixHorizontal    = [round(-fixationSize_pix / 2) round(fixationSize_pix / 2) 0 0];
fixVertical      = [0 0 round(-fixationSize_pix / 2) round(fixationSize_pix / 2)];
fixCoords        = [fixHorizontal; fixVertical];

timing.cfilower = 1000; % 1-2 s jittered fixation
timing.cfiupper = 2000;

%% Settings for inward moving circular grating
visualAngleGrating = 10;
gratingSize        = visualAngleGrating * screen.ppd;
gratingRadius      = round(gratingSize / 2);
gratingSize        = 2 * gratingRadius;

driftFreq      = 2;
nFramesInCycle = round((1 / driftFreq) / ifi);

movieDurationSecs = 2;
nFramesTotal      = round(movieDurationSecs * frameRate);

gratingDim      = [0 0 2 * gratingRadius 2 * gratingRadius];
gratingYpos     = screen.centerY;
gratingXpos     = screen.centerX;
gratingPosition = CenterRectOnPointd(gratingDim, gratingXpos, gratingYpos);

priorityLevel = MaxPriority(ptbWindow);
Priority(priorityLevel);

%% Generate grating textures
[x, y] = meshgrid(-gratingRadius:gratingRadius, -gratingRadius:gratingRadius);
f      = 0.55 * 2 * pi;

L                     = 2 * gratingRadius + 1;
w1D                   = hann(L);
xx                    = linspace(-gratingRadius, gratingRadius, L);
[X, Y]                = meshgrid(xx);
r                     = sqrt(X.^2 + Y.^2);
w2D                   = zeros(L);
w2D(r <= gratingRadius) = interp1(xx, w1D, r(r <= gratingRadius));

radius                          = sqrt(x.^2 + y.^2);
maxRadius                       = gratingRadius;
taperStart                      = maxRadius * 0.25;
taperMask                       = 0.5 * (1 + cos(pi * (radius - taperStart) / (maxRadius - taperStart)));
taperMask(radius <= taperStart) = 1;
taperMask(radius > maxRadius)   = 0;

contrastLevels = [0.25, 0.5, 0.75, 1];
tex_c25  = zeros(1, nFramesInCycle);
tex_c50  = zeros(1, nFramesInCycle);
tex_c75  = zeros(1, nFramesInCycle);
tex_c100 = zeros(1, nFramesInCycle);

for jFrame = 1:nFramesInCycle
    phase = (jFrame / nFramesInCycle) * 2 * pi;
    m     = sin(sqrt(x.^2 + y.^2) / f + phase);

    grating_c25  = (w2D .* (inc * m) + gray) * contrastLevels(1);
    grating_c25  = grating_c25 .* taperMask + (gray / 2) * (1 - taperMask);

    grating_c50  = (w2D .* (inc * m) + gray) * contrastLevels(2);

    grating_c75  = (w2D .* (inc * m) + gray) * contrastLevels(3);
    grating_c75  = grating_c75 .* taperMask + (gray / 2) * (1 - taperMask);

    grating_c100 = (w2D .* (inc * m) + gray) * contrastLevels(4);
    grating_c100 = grating_c100 .* taperMask + (gray / 2) * (1 - taperMask);

    tex_c25(jFrame)  = Screen('MakeTexture', ptbWindow, grating_c25);
    tex_c50(jFrame)  = Screen('MakeTexture', ptbWindow, grating_c50);
    tex_c75(jFrame)  = Screen('MakeTexture', ptbWindow, grating_c75);
    tex_c100(jFrame) = Screen('MakeTexture', ptbWindow, grating_c100);
end

%% Instruction screen
DrawFormattedText(ptbWindow, startExperimentText, 'center', 'center', black);
Screen('Flip', ptbWindow);
clc;
disp(upper('Participant is reading the contrast detection instructions.'));
waitResponse = 1;
while waitResponse
    [~, ~] = KbWait(-1, 2);
    waitResponse = 0;
end

HideCursor(whichScreen);

%% Attempt loop until pass criterion
passed = false;
attempt = 0;

while ~passed
    attempt = attempt + 1;
    timing.startTime = datestr(now, 'dd/mm/yy-HH:MM:SS');

    % Build one-step-adjacent random sequence
    gratingSequence = zeros(1, exp.nTrials);
    gratingSequence(1) = randi(4);
    for iSeq = 2:exp.nTrials
        prev = gratingSequence(iSeq - 1);
        candidates = [];
        if prev > 1
            candidates(end + 1) = prev - 1; %#ok<AGROW>
        end
        if prev < 4
            candidates(end + 1) = prev + 1; %#ok<AGROW>
        end
        gratingSequence(iSeq) = candidates(randi(numel(candidates)));
    end

    data                  = struct;
    data.grating          = gratingSequence;
    data.contrastChange   = [NaN, diff(gratingSequence)]; % >0 increase, <0 decrease
    data.responses        = nan(1, exp.nTrials); % 1 = up, -1 = down, 0 = none
    data.correct          = nan(1, exp.nTrials);
    data.reactionTime     = nan(1, exp.nTrials);
    data.trlDuration      = nan(1, exp.nTrials);
    data.cfi              = nan(1, exp.nTrials);
    data.tooSlow          = false(1, exp.nTrials);

    clc;
    disp(['CONTRAST DETECTION TEST (attempt ' num2str(attempt) ')...']);

    for trl = 1:exp.nTrials
        tic;
        currentGrating = gratingSequence(trl);

        if currentGrating == 1
            gratingForm = ' 25% contrast';
        elseif currentGrating == 2
            gratingForm = ' 50% contrast';
        elseif currentGrating == 3
            gratingForm = ' 75% contrast';
        else
            gratingForm = '100% contrast';
        end

        %% Fixation cross (black only, 1-2 s jitter)
        Screen('FillRect', ptbWindow, backgroundColorGray);
        Screen('Flip', ptbWindow);
        data.cfi(trl) = (randsample(timing.cfilower:timing.cfiupper, 1)) / 1000;

        Screen('DrawLines', ptbWindow, fixCoords, fixationLineWidth, ...
            fixationColorBlack, [screen.centerX screen.centerY], 2);
        Screen('Flip', ptbWindow);
        WaitSecs(data.cfi(trl));

        %% Present grating and collect response
        Screen('Flip', ptbWindow);
        responseGiven = false;
        maxProbeDuration = movieDurationSecs;
        frameDuration = maxProbeDuration / nFramesInCycle;

        probeStartTime = GetSecs;
        whileCount = 1;
        while (GetSecs - probeStartTime) < maxProbeDuration
            if currentGrating == 1
                Screen('DrawTexture', ptbWindow, tex_c25(whileCount), [], gratingPosition);
            elseif currentGrating == 2
                Screen('DrawTexture', ptbWindow, tex_c50(whileCount), [], gratingPosition);
            elseif currentGrating == 3
                Screen('DrawTexture', ptbWindow, tex_c75(whileCount), [], gratingPosition);
            else
                Screen('DrawTexture', ptbWindow, tex_c100(whileCount), [], gratingPosition);
            end
            Screen('Flip', ptbWindow);

            if trl > 1 && ~responseGiven
                [keyIsDown, responseTime, keyCode] = KbCheck(-1);
                if keyIsDown
                    if keyCode(upKeyCode)
                        responseGiven = true;
                        data.responses(trl) = 1;
                        data.reactionTime(trl) = responseTime - probeStartTime;
                    elseif keyCode(downKeyCode)
                        responseGiven = true;
                        data.responses(trl) = -1;
                        data.reactionTime(trl) = responseTime - probeStartTime;
                    end
                end
            end

            elapsedTime = GetSecs - probeStartTime;
            expectedTimeCurrentFrame = whileCount * frameDuration;
            waitTime = expectedTimeCurrentFrame - elapsedTime;
            if waitTime > 0
                WaitSecs(waitTime);
            end

            whileCount = whileCount + 1;
            if whileCount > nFramesInCycle
                whileCount = 1;
            end
        end

        %% Score response (trial 1 unscored)
        if trl == 1
            data.correct(trl) = NaN;
            feedbackText = 'n/a (first)';
        else
            expectedDir = sign(data.contrastChange(trl)); % 1 increase, -1 decrease
            if ~responseGiven
                data.responses(trl) = 0;
                data.tooSlow(trl) = true;
                data.correct(trl) = 0;
                feedbackText = 'TOO SLOW!';
                DrawFormattedText(ptbWindow, 'TOO SLOW!', 'center', 'center', black);
                Screen('Flip', ptbWindow);
                WaitSecs(1);
            elseif data.responses(trl) == expectedDir
                data.correct(trl) = 1;
                feedbackText = 'Correct';
            else
                data.correct(trl) = 0;
                feedbackText = 'Incorrect';
            end
        end

        %% CW trial output
        scored = data.correct(2:trl);
        scored = scored(~isnan(scored));
        if isempty(scored)
            overall_accuracy = NaN;
            accStr = 'n/a';
        else
            overall_accuracy = round((sum(scored) / numel(scored)) * 100);
            accStr = sprintf('%3d%%', overall_accuracy);
        end

        if isnan(data.reactionTime(trl))
            rtStr = ' n/a';
        else
            rtStr = sprintf('%4.2fs', data.reactionTime(trl));
        end

        if trl == 1
            changeStr = '---';
        elseif data.contrastChange(trl) > 0
            changeStr = 'UP ';
        else
            changeStr = 'DOWN';
        end

        fprintf('Contrast trial %2d/%d | %s | change: %s | %s | Acc: %s | RT: %s\n', ...
            trl, exp.nTrials, gratingForm, changeStr, feedbackText, accStr, rtStr);

        data.trlDuration(trl) = toc;
    end

    %% Overall accuracy (exclude first trial)
    scoredCorrect = data.correct(2:end);
    nScored = numel(scoredCorrect);
    nCorrect = sum(scoredCorrect == 1);
    overall_accuracy = round((nCorrect / nScored) * 100);
    data.overallAccuracy = overall_accuracy;
    data.nCorrect = nCorrect;
    data.nScored = nScored;
    data.attempt = attempt;

    timing.endTime = datestr(now, 'dd/mm/yy-HH:MM:SS');
    startTime = datetime(timing.startTime, 'InputFormat', 'dd/MM/yy-HH:mm:ss');
    endTime = datetime(timing.endTime, 'InputFormat', 'dd/MM/yy-HH:mm:ss');
    timing.duration = seconds(endTime - startTime);

    fprintf('CONTRAST DETECTION overall accuracy: %d%% (%d/%d). Attempt %d.\n', ...
        overall_accuracy, nCorrect, nScored, attempt);

    %% End-of-attempt screen
    if overall_accuracy >= exp.passThreshold
        passed = true;
        endText = sprintf([
            'Contrast detection finished.\n\n' ...
            'Accuracy: %d%% (%d/%d)\n\n' ...
            'Press any key to continue...'], ...
            overall_accuracy, nCorrect, nScored);
    else
        endText = sprintf([
            'Accuracy: %d%% (%d/%d)\n\n' ...
            'At least %d%% correct is required.\n\n' ...
            'Press any key to retry...'], ...
            overall_accuracy, nCorrect, nScored, exp.passThreshold);
    end

    DrawFormattedText(ptbWindow, endText, 'center', 'center', black);
    Screen('Flip', ptbWindow);
    waitResponse = 1;
    while waitResponse
        [~, ~] = KbWait(-1, 2);
        waitResponse = 0;
    end

    %% Save after each completed attempt (overwrite primary file on pass)
    subjectID = num2str(subject.ID);
    filePath = fullfile(DATA_PATH, subjectID);
    mkdir(filePath);

    saves                    = struct;
    saves.data               = data;
    saves.data.upKeyCode     = upKeyCode;
    saves.data.downKeyCode   = downKeyCode;
    saves.experiment         = exp;
    saves.screen             = screen;
    saves.startExperimentText = startExperimentText;
    saves.subjectID          = subjectID;
    saves.subject            = subject;
    saves.timing             = timing;
    saves.passed             = passed;
    saves.gratingSequence    = gratingSequence;

    if passed
        fileName = [subjectID, '_contrast_detection.mat'];
    else
        fileName = [subjectID, '_contrast_detection_attempt' num2str(attempt) '.mat'];
    end
    save(fullfile(filePath, fileName), 'saves');
    disp(['Saved: ' fullfile(filePath, fileName)]);
end

%% Close Psychtoolbox window
Priority(0);
Screen('Close');
Screen('CloseAll');
disp('CONTRAST DETECTION TEST PASSED');
