%% GCP Gratings Task
% This code requires PsychToolbox. https://psychtoolbox.org
% This was tested with PsychToolbox version 3.0.15, and with MATLAB R2023b.
%
% The code for the grating stimulus was copied from DriftDemo and modified
% to display a (masked) animated concentric grating moving inward. Adapted
% from van Es and Schoffelen, 2019.
% https://github.com/Donders-Institute/dyncon_erfosc/blob/master/concentric_grating_experiment.m

%% Initialize EEG and ET

% Start of block message in CW
if TRAINING == 1
    disp('START OF BLOCK 0 (TRAINING)');
else
    disp(['START OF BLOCK ' num2str(BLOCK)]);
end

% Calibrate ET (Tobii Pro Fusion)
disp('CALIBRATING ET...');
calibrateET;

if TRAINING == 0
    % Start recording EEG
    disp('STARTING EEG RECORDING...');
    initEEG;

    % Wait ten seconds to initialize EEG
    clc
    disp('INITIALIZING EEG... PLEASE WAIT 10 SECONDS')
    for i=1:10
        if i > 1
            wbar = findall(0,'type','figure','tag','TMWWaitbar');
            delete(wbar)
        end
        waitbar(i/10, 'INITIALIZING EEG');
        pause(1);
    end
    wbar = findall(0,'type','figure','tag','TMWWaitbar');
    delete(wbar)
    clc
    disp('EEG INITIALIZED!')
end

% Hide cursor on participant screen
HideCursor(whichScreen);

%% Define TRIGGERS
TASK_START = 10; % trigger for ET cutting

BLOCK1                   = 11; % Trigger for start of block 1
BLOCK2                   = 12; % Trigger for start of block 2
BLOCK3                   = 13; % Trigger for start of block 3
BLOCK4                   = 14; % Trigger for start of block 4
BLOCK0                   = 15; % Trigger for start of training block (block 0)

FIXCROSSR                = 16; % Trigger for white (task) fixation cross
FIXCROSSB                = 17; % Trigger for black fixation cross

PRESENTATION_C25_TASK    = 51; % Trigger for presentation of 25% contrast concentric dynamic inward grating WITH button press response
PRESENTATION_C50_TASK    = 52; % Trigger for presentation of 50% contrast concentric dynamic inward grating WITH button press response
PRESENTATION_C75_TASK    = 53; % Trigger for presentation of 75% contrast concentric dynamic inward grating WITH button press response
PRESENTATION_C100_TASK   = 54; % Trigger for presentation of 100% contrast concentric dynamic inward grating WITH button press response
PRESENTATION_C25_NOTASK  = 61; % Trigger for presentation of 25% contrast concentric dynamic inward grating WITHOUT button press response
PRESENTATION_C50_NOTASK  = 62; % Trigger for presentation of 50% contrast concentric dynamic inward grating WITHOUT button press response
PRESENTATION_C75_NOTASK  = 63; % Trigger for presentation of 75% contrast concentric dynamic inward grating WITHOUT button press response
PRESENTATION_C100_NOTASK = 64; % Trigger for presentation of 100% contrast concentric dynamic inward grating WITHOUT button press response

BLOCK1_END               = 71; % End of block 1
BLOCK2_END               = 72; % End of block 2
BLOCK3_END               = 73; % End of block 3
BLOCK4_END               = 74; % End of block 4
BLOCK0_END               = 75; % End of block 0

RESP_YES                 = 87; % Trigger for response yes (spacebar)
RESP_NO                  = 88; % Trigger for response no (no input)

TASK_END                 = 90; % Trigger for ET cutting

%% Set up experiment parameters
% Block and Trial Number
exp.nTrlTrain = 10; % n gratings per training block
exp.nTrlTask = 176; % n gratings per task block

if TRAINING == 1
    exp.nTrials = exp.nTrlTrain;
else
    exp.nTrials = exp.nTrlTask;
end

% Enable (= 1) or disable (= 0) screenshots
enableScreenshots = 1;

%% Set up text parameters
% Define startExperimentText
startExperimentText = [
    'You will see a series of gratings. \n\n' ...
    'Between gratings, a fixation cross \n\n' ...
    'will appear on the screen. If the \n\n' ...
    'fixation cross appers in WHITE, it is \n\n' ...
    'your task to press SPACE as soon as the next \n\n' ...
    'grating appears. Use your right hand. Please \n\n' ...
    'always look at the center of the screen. \n\n' ...
    '\n\n' ...
    'Press any key to continue...'];
if TRAINING == 1
    loadingText = 'Loading TRAINING...';
elseif TRAINING == 0
    loadingText = 'Loading TASK...';
end

%% Set up standard Psychtoolbox Settings
global GL;
AssertOpenGL; % Check OpenGL Psychtoolbox

% Disable clipping of text
global ptb_drawformattedtext_disableClipping;
ptb_drawformattedtext_disableClipping = 1;

% Set verbosity to disallow CW output
Screen('Preference','Verbosity', 0);
%Screen('Preference', 'SkipSyncTests', 0); % For linux (can be 0)

% Retrieve key code for spacebar
spaceKeyCode = KbName('Space');

%% Imaging set up
screen.ID = whichScreen; % Get index for stimulus presentation screen

if gray == white % Ensure well defined gray, even on floating point framebuffers
    gray = white / 2;
end

% Contrast 'inc'rement range for given white and gray values
inc = white - gray;

% Background color
backgroundColorGray = 256; % Needs to be fixed to 256 for gratings blending into background

% Open a double buffered fullscreen window and select a gray background color:
[ptbWindow, winRect] = Screen('OpenWindow', screen.ID, backgroundColorGray);

% Get screen size and center coordinates
[screen.centerX, screen.centerY] = RectCenter(winRect); % Screen center in pixels
screen.width                     = 48; % Screen width in cm
screen.height                    = 29.89; % Screen height in cm
screen.resolutionX               = 800; % Screen resolution width in pixels
screen.resolutionY               = 600; % Screen resolution height in pixels
screen.viewDist                  = 80; % Viewing distance in cm from participant on head rest to screen center

% Calculate visual parameters
screen.totVisDeg = 2*atan(screen.width / (2*screen.viewDist))*(180/pi); % Calculate degrees of visual angle
screen.ppd       = screen.resolutionX / screen.totVisDeg; % Pixels per degree
% MethLab 20.5761 ppd; estimated with MeasureDpi function: 20
screen.ppd       = 50;

% Get frame duration
ifi              = Screen('GetFlipInterval', ptbWindow);
frameRate        = Screen('FrameRate', screen.ID); % MethLab 100 Hz

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', ptbWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% Text parameters
Screen('TextSize', ptbWindow, 25); % Font size for instructions and stimuli

% Show loading text
DrawFormattedText(ptbWindow,loadingText,'center','center',black);
Screen('Flip',ptbWindow);

%% Fixation cross parameters
% Size
fixationSize_dva      = .35;             % Size of fixation cross in degress of visual angle
fixationLineWidth     = 1.5;            % Line width of fixation cross

% Color
blackLevel            = BlackIndex(ptbWindow);   % Get RGB values for black
fixationColorBlack    = [0 0 0];                 % Black fixation cross [0 0 0]
whiteLevel            = WhiteIndex(ptbWindow);   % Get RGB values for white
fixationColorWhite    = [1023 1023 1023];        % White fixation cross [1023 1023 1023]

% Location
fixationSize_pix      = round(fixationSize_dva*screen.ppd);
fixHorizontal         = [round(-fixationSize_pix/2) round(fixationSize_pix/2) 0 0];
fixVertical           = [0 0 round(-fixationSize_pix/2) round(fixationSize_pix/2)];
fixCoords             = [fixHorizontal; fixVertical];
fixPos                = [screen.centerX, screen.centerY];

% Temporal parameters
timing.cfilower       = 2000; % Lower limit of CFI duration
timing.cfiupper       = 3000; % Upper limit of CFI duration
timing.cfi_task       = 0.5;  % Duration of white fixation cross

%% Settings for inward moving circular grating
% Size
visualAngleGrating    = 10; %7.1
visualAngleLocation   = 15;
gratingSize           = visualAngleGrating*screen.ppd; % Grating stimulus size in pixels
gratingRadius         = round(gratingSize/2); % Grating can only exist of integers -> round
gratingSize           = 2*gratingRadius; % To prevent consistency errors, redifine gratingSize

% Frequency
driftFreq             = 2; % Every pixel of the grating completes two cycles per second (black-white-black)
nFramesInCycle        = round((1/driftFreq)/ifi); % Temporal period, in frames, of the drifting grating

% Set duration
movieDurationSecs = 2;
nFramesTotal = round(movieDurationSecs * frameRate); % Convert movieDuration in seconds to duration in frames

% Location
gratingDim      = [0 0 2*gratingRadius 2*gratingRadius];
gratingYpos     = screen.centerY;
gratingXpos     = screen.centerX;
gratingPosition = CenterRectOnPointd(gratingDim, gratingXpos, gratingYpos); % Move the object to those coordinates

% Use realtime priority for better timing precision
priorityLevel = MaxPriority(ptbWindow);
Priority(priorityLevel);

%% Generate grating textures
% Generate stimulus
[x,y]                           = meshgrid(-gratingRadius:gratingRadius,-gratingRadius:gratingRadius);
f                               = 0.55*2*pi; % Period of the grating

% Circular hanning mask
L                               = 2*gratingRadius+1;
w1D                             = hann(L); % 1D hann window
xx                              = linspace(-gratingRadius,gratingRadius,L);
[X,Y]                           = meshgrid(xx);
r                               = sqrt( X.^2 + Y.^2 );
w2D                             = zeros(L);
w2D(r<=gratingRadius)           = interp1(xx,w1D,r(r<=gratingRadius)); % 2D hanning window

% Tapering mask for the high contrast condition
[rows, cols]                    = size(x);
radius                          = sqrt(x.^2 + y.^2); % Distance from the center
maxRadius                       = gratingRadius; % Maximum radius of the grating
taperStart                      = maxRadius * 0.25; % Start tapering at 50% of the grating radius
taperMask                       = 0.5 * (1 + cos(pi * (radius - taperStart) / (maxRadius - taperStart)));
taperMask(radius <= taperStart) = 1; % Flat region in the center
taperMask(radius > maxRadius)   = 0; % Fully tapered outside the grating
% The mask is created based on the radial distance (radius) from the center
% of the grating. The cosine taper smoothly declines the intensity starting
% at 50% of the grating size (taperStart) and reaches zero at the maximum radius.

% Compute each frame of the movie and convert those frames stored in
% MATLAB matrices, into Psychtoolbox OpenGL textures using 'MakeTexture'
contrastLevels                  = [0.25, 0.5, 0.75, 1];
tex                             = zeros(nFramesInCycle,1);
for jFrame = 1:nFramesInCycle
    phase                       = (jFrame / nFramesInCycle) * 2 * pi; % Change the phase of the grating according to frame number
    m                           = sin(sqrt(x.^2 + y.^2) / f + phase); % Formula for sinusoidal grating

    % inc*m fluctuates from [-gray, gray]. Multiply this with the hanning mask to let the grating die off at 0
    % 25% contrast grating
    grating_c25                  = (w2D .* (inc * m) + gray) * contrastLevels(1);
    % Multiply by taperMask to gradually fade grating towards gray background color (64)
    grating_c25                 = grating_c25 .* taperMask + (gray/2) * (1 - taperMask);

    % 50% contrast grating
    grating_c50                  = (w2D .* (inc * m) + gray) * contrastLevels(2);

    % 75% contrast grating
    grating_c75                  = (w2D .* (inc * m) + gray) * contrastLevels(3);
    % Multiply by taperMask to gradually fade grating towards gray background color (64)
    grating_c75                 = grating_c75 .* taperMask + (gray/2) * (1 - taperMask);

    % 100% contrast grating
    grating_c100                 = (w2D .* (inc * m) + gray) * contrastLevels(4);
    % Multiply by taperMask to gradually fade grating towards gray background color (64)
    grating_c100                 = grating_c100 .* taperMask + (gray/2) * (1 - taperMask);

    % Create textures for low and high contrast gratings
    tex_c25(jFrame)              = Screen('MakeTexture', ptbWindow, grating_c25);
    tex_c50(jFrame)              = Screen('MakeTexture', ptbWindow, grating_c50);
    tex_c75(jFrame)              = Screen('MakeTexture', ptbWindow, grating_c75);
    tex_c100(jFrame)             = Screen('MakeTexture', ptbWindow, grating_c100);
end

%% Create data structure for preallocating data
data                             = struct;
if TRAINING == 1
    gratingSequence              = [2 3 1 4 3 1 2 4 1 4];
else
    nums                         = repmat(1:4, 1, floor(exp.nTrials/4));
    gratingSequence              = nums(randperm(length(nums), exp.nTrials)); % Define grating sequence
end
%countSequence                   = histcounts(gratingSequence, 1:5); % Check equal occurence
data.grating(1, exp.nTrials)     = NaN; % Saves grating form (see below)
% grating = 1 is 25% contrast concentric dynamic inward
% grating = 2 is 50% contrast concentric dynamic inward
% grating = 3 is 75% contrast concentric dynamic inward
% grating = 4 is 100% contrast concentric dynamic inward
data.whiteCross(1, exp.nTrials)  = NaN; % Binary measure for task condition
data.responses(1, exp.nTrials)   = NaN; % Binary measure for (no) response
data.correct(1, exp.nTrials)     = NaN; % Binary measure for correct responses
data.reactionTime(1:exp.nTrials) = NaN; % Reaction time
data.fixation(1:exp.nTrials)     = NaN; % Fixation check info
data.trlDuration(1:exp.nTrials)  = NaN; % Trial duration in seconds
count5trials                     = NaN; % Initialize accuracy reminder loop variable

%% Show task instruction text
DrawFormattedText(ptbWindow, startExperimentText, 'center', 'center', black);
Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1); % black background for photo diode
Screen('Flip',ptbWindow);
clc;
disp(upper('Participant is reading the instructions.'));
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

% Send triggers for start of task (ET cutting)
if TRAINING == 1
    Eyelink('Message', num2str(TASK_START));
    Eyelink('command', 'record_status_message "START"');
else
    Eyelink('Message', num2str(TASK_START));
    Eyelink('command', 'record_status_message "START"');
    sendtrigger(TASK_START,port,SITE,stayup);
end

% Send triggers for block and output
if BLOCK == 1
    TRIGGER = BLOCK1;
elseif BLOCK == 2
    TRIGGER = BLOCK2;
elseif BLOCK == 3
    TRIGGER = BLOCK3;
elseif BLOCK == 4
    TRIGGER = BLOCK4;
end

if TRAINING == 1
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "START BLOCK"');
else
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "START BLOCK"');
    sendtrigger(TRIGGER,port,SITE,stayup);
end

% Experiment prep
HideCursor(whichScreen); % Make sure to hide cursor from participant screen
timing.startTime = datestr(now, 'dd/mm/yy-HH:MM:SS'); % Measure duration

%% Experiment Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
disp('GCP GRATING TASK...');
for trl = 1:exp.nTrials
    tic;
    % Add trial grating info
    if gratingSequence(trl) == 1
        gratingForm = ' 25% contrast';
    elseif gratingSequence(trl) == 2
        gratingForm = ' 50% contrast';
    elseif gratingSequence(trl) == 3
        gratingForm = ' 75% contrast';
    elseif gratingSequence(trl) == 4
        gratingForm = '100% contrast';
    end

    % Randomized selection of task (white fication cross) trials (10%)
    if TRAINING == 0
        if randi(10) == 1
            data.whiteCross(trl) = 1;
        else
            data.whiteCross(trl) = 0;
        end
    elseif TRAINING == 1
        if randi(3) == 1
            data.whiteCross(trl) = 1;
        else
            data.whiteCross(trl) = 0;
        end
    end

    %% Present fixation cross (white for task condition)
    % Fill gray screen
    Screen('FillRect', ptbWindow, backgroundColorGray);
    Screen('Flip', ptbWindow);
    % Set jittered trial-specific durations for CFIs
    timing.cfi(trl) = (randsample(timing.cfilower:timing.cfiupper, 1))/1000; % Randomize the jittered central fixation interval on trial
    start_time = GetSecs;
    while (GetSecs - start_time) < timing.cfi(trl)
        if data.whiteCross(trl) == 0 % No task condition
            Screen('DrawLines', ptbWindow, fixCoords,fixationLineWidth,fixationColorBlack,[screen.centerX screen.centerY],2);
            Screen('Flip', ptbWindow);
            screenshot('GCP_screenshot_blackcross.png', ptbWindow, enableScreenshots);
            TRIGGER = FIXCROSSB;
            if TRAINING == 1
                Eyelink('Message', num2str(TRIGGER));
                Eyelink('command', 'record_status_message "FIXCROSS"');
            else
                Eyelink('Message', num2str(TRIGGER));
                Eyelink('command', 'record_status_message "FIXCROSS"');
                sendtrigger(TRIGGER,port,SITE,stayup);
            end
            WaitSecs(timing.cfi(trl));
        elseif data.whiteCross(trl) == 1 % Task condition
            Screen('DrawLines', ptbWindow, fixCoords,fixationLineWidth,fixationColorWhite,[screen.centerX screen.centerY],2);
            Screen('Flip', ptbWindow);
            TRIGGER = FIXCROSSR;
            screenshot('GCP_screenshot_whiteCross.png', ptbWindow, enableScreenshots);
            if TRAINING == 1
                Eyelink('Message', num2str(TRIGGER));
                Eyelink('command', 'record_status_message "FIXCROSS"');
            else
                Eyelink('Message', num2str(TRIGGER));
                Eyelink('command', 'record_status_message "FIXCROSS"');
                sendtrigger(TRIGGER,port,SITE,stayup);
            end
            WaitSecs(timing.cfi_task); % Show white cross for 500 ms
            Screen('DrawLines', ptbWindow, fixCoords,fixationLineWidth,fixationColorBlack,[screen.centerX screen.centerY],2);
            Screen('Flip', ptbWindow);
            TRIGGER = FIXCROSSB;
            if TRAINING == 1
                Eyelink('Message', num2str(TRIGGER));
                Eyelink('command', 'record_status_message "FIXCROSS"');
            else
                Eyelink('Message', num2str(TRIGGER));
                Eyelink('command', 'record_status_message "FIXCROSS"');
                sendtrigger(TRIGGER,port,SITE,stayup);
            end
            WaitSecs(timing.cfi(trl)-timing.cfi_task); % Show black cross for the rest of the CFI time
        end
    end

    %% Define grating depending on sequence number
    data.grating(trl) = gratingSequence(trl);
    % grating = 1 is 25% contrast concentric dynamic inward
    % grating = 2 is 50% contrast concentric dynamic inward
    % grating = 3 is 75% contrast concentric dynamic inward
    % grating = 4 is 100% contrast concentric dynamic inward

    %% Present grating and get response
    Screen('Flip', ptbWindow); % Preparatory flip
    responseGiven = false;
    maxProbeDuration = 2; % Maximum time to show the grating
    frameDuration = maxProbeDuration / length(tex);

    % Send presentation triggers
    if gratingSequence(trl) == 1 && data.whiteCross(trl) == 1
        TRIGGER = PRESENTATION_C25_TASK;
    elseif gratingSequence(trl) == 2 && data.whiteCross(trl) == 1
        TRIGGER = PRESENTATION_C50_TASK;
    elseif gratingSequence(trl) == 3 && data.whiteCross(trl) == 1
        TRIGGER = PRESENTATION_C75_TASK;
    elseif gratingSequence(trl) == 4 && data.whiteCross(trl) == 1
        TRIGGER = PRESENTATION_C100_TASK;
    elseif gratingSequence(trl) == 1 && data.whiteCross(trl) == 0
        TRIGGER = PRESENTATION_C25_NOTASK;
    elseif gratingSequence(trl) == 2 && data.whiteCross(trl) == 0
        TRIGGER = PRESENTATION_C50_NOTASK;
    elseif gratingSequence(trl) == 3 && data.whiteCross(trl) == 0
        TRIGGER = PRESENTATION_C75_NOTASK;
    elseif gratingSequence(trl) == 4 && data.whiteCross(trl) == 0
        TRIGGER = PRESENTATION_C100_NOTASK;
    end

    if TRAINING == 1
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "PRESENTATION"');
    else
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "PRESENTATION"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end

    probeStartTime = GetSecs;
    whileCount = 1;
    % Draw gratings depending on gratingSequence
    while (GetSecs - probeStartTime) < maxProbeDuration
        if gratingSequence(trl) == 1 % 25% contrast grating
            Screen('DrawTexture', ptbWindow, tex_c25(whileCount), [], gratingPosition);
            Screen('Flip', ptbWindow);
        elseif gratingSequence(trl) == 2 % 50% contrast grating
            Screen('DrawTexture', ptbWindow, tex_c50(whileCount), [], gratingPosition);
            Screen('Flip', ptbWindow);
        elseif gratingSequence(trl) == 3 % 75% contrast grating
            Screen('DrawTexture', ptbWindow, tex_c75(whileCount), [], gratingPosition);
            Screen('Flip', ptbWindow);
        elseif gratingSequence(trl) == 4 % 100% contrast grating
            Screen('DrawTexture', ptbWindow, tex_c100(whileCount), [], gratingPosition);
            Screen('Flip', ptbWindow);
        end
        screenshot(sprintf('GCP_screenshot_%s.png', gratingForm), ptbWindow, enableScreenshots);

        % Check for participant response
        if ~responseGiven
            [keyIsDown, responseTime, keyCode] = KbCheck;
            if keyIsDown
                responseGiven = true;
                data.reactionTime(trl) = responseTime - probeStartTime;
                data.responses(trl) = 1; % Response made
            end
        end

        % Calculate elapsed time for each while loop and wait
        elapsedTime = GetSecs - probeStartTime;
        expectedTimeCurrentFrame = whileCount * frameDuration;
        waitTime = expectedTimeCurrentFrame - elapsedTime;
        if waitTime > 0
            WaitSecs(waitTime);
        end

        % Increment while loop count for textures
        whileCount = whileCount+1;
    end

    % If no response is given, record default
    if ~responseGiven
        data.responses(trl) = 0; % No response
        data.reactionTime(trl) = NaN;
    end

    %% Check if response was correct
    if data.whiteCross(trl) == 1 && data.responses(trl) == 1 % White fixation cross + button press = correct
        data.correct(trl) = 1;
        feedbackText = 'Correct!  ';
    elseif data.whiteCross(trl) == 0 && data.responses(trl) == 0 % No white fixation cross + no button press = correct
        data.correct(trl) = 1;
        feedbackText = 'Correct!  ';
    else % Anything else is wrong response
        data.correct(trl) = 0;
        feedbackText = 'Incorrect!';
    end

    %% Feedback for training block and CW output
    % Give feedback in training block
    if TRAINING == 1
        DrawFormattedText(ptbWindow,feedbackText,'center','center',black);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('Flip',ptbWindow);
        WaitSecs(2);
        % Give feedback for no response (too slow)
    elseif TRAINING == 0 && data.correct(trl) == 0 && data.responses(trl) == 0
        feedbackText = 'TOO SLOW! ';
        DrawFormattedText(ptbWindow,feedbackText,'center','center',black);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('Flip',ptbWindow);
        WaitSecs(2);
    end

    %% Dynamically compute accuracy for past 10 trials and remind participant if accuracy drops below threshhold of 90%
    responsesLastTrials = 0;
    if trl >= 10
        responsesLastTrials = data.correct(trl-9 : trl);
        percentLastTrialsCorrect = sum(responsesLastTrials)*10;
        if percentLastTrialsCorrect < 90 && count5trials <= trl-5
            count5trials = trl;
            feedbackLastTrials = ['Your accuracy has declined!'...
                '\n\n ' ...
                '\n\n Please stay focused on the task!'];
            disp(['Participant was made aware of low accuracy in the last 10 trials: ' num2str(percentLastTrialsCorrect) ' %. [' num2str(responsesLastTrials) ']']);
            DrawFormattedText(ptbWindow,feedbackLastTrials,'center','center',black);
            Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
            Screen('Flip',ptbWindow);
            WaitSecs(3);
        end
    end

    %% Trial Info CW output
    overall_accuracy = round((sum(data.correct(1:trl))/trl)*100);
    reactionTime = num2str(round(data.reactionTime(trl), 2), '%.2f');
    if trl < 10
        disp(['Response to Trial ' num2str(trl) '/' num2str(exp.nTrials) ...
            ' in Block ' num2str(BLOCK) ' is ' feedbackText '  (White FixCross: ' ...
            '' num2str(data.whiteCross(trl)) ' | Acc: ' num2str(overall_accuracy) ...
            '% | RT: ' reactionTime 's | ' gratingForm ')']);
    else
        disp(['Response to Trial ' num2str(trl) '/' num2str(exp.nTrials) ...
            ' in Block ' num2str(BLOCK) ' is ' feedbackText ' (White FixCross: ' ...
            num2str(data.whiteCross(trl)) ' | Acc: ' num2str(overall_accuracy) ...
            '% | RT: ' reactionTime 's | ' gratingForm ')']);
    end

    % Save trial duration in seconds
    data.trlDuration(trl) = toc;
end

%% End task and save data

% Send triggers to end task
Screen('Flip',ptbWindow);

% Send triggers for block and output
if BLOCK == 1 && TRAINING == 1
    TRIGGER = BLOCK0_END; % Training block
elseif BLOCK == 1 && TRAINING == 0
    TRIGGER = BLOCK1_END;
elseif BLOCK == 2 && TRAINING == 0
    TRIGGER = BLOCK2_END;
elseif BLOCK == 3 && TRAINING == 0
    TRIGGER = BLOCK3_END;
elseif BLOCK == 4 && TRAINING == 0
    TRIGGER = BLOCK4_END;
end

if TRAINING == 1
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "END BLOCK"');
    disp('End of Training Block.');
else
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "END BLOCK"');
    sendtrigger(TRIGGER,port,SITE,stayup);
    disp(['End of Block ' num2str(BLOCK)]);
end

% Send triggers for end of task (ET cutting)
if TRAINING == 1
    Eyelink('Message', num2str(TASK_END));
    Eyelink('command', 'record_status_message "TASK_END"');
else
    Eyelink('Message', num2str(TASK_END));
    Eyelink('command', 'record_status_message "TASK_END"');
    sendtrigger(TASK_END,port,SITE,stayup);
end

%% Record block duration
timing.endTime = datestr(now, 'dd/mm/yy-HH:MM:SS');
% Convert to datetime objects
startTime = datetime(timing.startTime, 'InputFormat', 'dd/MM/yy-HH:mm:ss');
endTime = datetime(timing.endTime, 'InputFormat', 'dd/MM/yy-HH:mm:ss');
% Calculate block duration in seconds
timing.duration = seconds(endTime - startTime);

%% Save data
subjectID = num2str(subject.ID);
filePath = fullfile(DATA_PATH, subjectID);
mkdir(filePath)

if TRAINING == 1
    fileName = [subjectID, '_training.mat'];
elseif TRAINING == 0
    fileName = [subjectID '_', TASK, '_block' num2str(BLOCK), '.mat'];
end

% Save data
saves                           = struct;
saves.data                      = data;
saves.data.spaceKeyCode         = spaceKeyCode;
saves.data.reactionTime         = data.reactionTime;
saves.experiment                = exp;
saves.screen                    = screen;
saves.startExperimentText       = startExperimentText;
saves.subjectID                 = subjectID;
saves.subject                   = subject;
saves.timing                    = timing;

% Save triggers
trigger                         = struct;
trigger.TASK_START              = TASK_START;
trigger.BLOCK1                  = BLOCK1;
trigger.BLOCK2                  = BLOCK2;
trigger.BLOCK3                  = BLOCK3;
trigger.BLOCK4                  = BLOCK4;
trigger.BLOCK0                  = BLOCK0;

trigger.FIXCROSSR               = FIXCROSSR;
trigger.FIXCROSSB               = FIXCROSSB;

trigger.PRESENTATION_C25_TASK    = PRESENTATION_C25_TASK;
trigger.PRESENTATION_C50_TASK    = PRESENTATION_C50_TASK;
trigger.PRESENTATION_C75_TASK    = PRESENTATION_C75_TASK;
trigger.PRESENTATION_C100_TASK   = PRESENTATION_C100_TASK;
trigger.PRESENTATION_C25_NOTASK  = PRESENTATION_C25_NOTASK;
trigger.PRESENTATION_C50_NOTASK  = PRESENTATION_C50_NOTASK;
trigger.PRESENTATION_C75_NOTASK  = PRESENTATION_C75_NOTASK;
trigger.PRESENTATION_C100_NOTASK = PRESENTATION_C100_NOTASK;

trigger.BLOCK1_END              = BLOCK1_END;
trigger.BLOCK2_END              = BLOCK2_END;
trigger.BLOCK3_END              = BLOCK3_END;
trigger.BLOCK4_END              = BLOCK4_END;
trigger.BLOCK0_END              = BLOCK0_END;

trigger.RESP_YES                = RESP_YES;
trigger.RESP_NO                 = RESP_NO;

trigger.TASK_END                = TASK_END;

%% Stop and close EEG and ET recordings
if TRAINING == 1
    disp('TRAINING FINISHED...');
else
    disp(['BLOCK ' num2str(BLOCK) ' FINISHED...']);
end
disp('SAVING DATA...');
save(fullfile(filePath, fileName), 'saves', 'trigger');
closeEEGandET;

try
    PsychPortAudio('Close');
catch
end

%% Show break instruction text
if TRAINING == 1 && BLOCK == 1
    breakInstructionText = 'Well done! \n\n Press any key to finalize the training block.';

elseif TRAINING == 0 && BLOCK == 4
    breakInstructionText = ['End of the Task! ' ...
        '\n\n Thank you very much for your participation.'...
        '\n\n Please press any key to finalize the exp.'];
else
    breakInstructionText = ['Break! Rest for a while... ' ...
        '\n\n Press any key to start the break.'];
end

DrawFormattedText(ptbWindow,breakInstructionText,'center','center',black);
Screen('Flip',ptbWindow);
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

% Show final screen
if BLOCK == 4 && TRAINING == 0
    FinalText = ['You are done.' ...
        '\n\n Have a great day!'];
    DrawFormattedText(ptbWindow, FinalText, 'center', 'center', black);
elseif BLOCK == 1 && TRAINING == 0 || BLOCK == 2 && TRAINING == 0 || BLOCK == 3 && TRAINING == 0
    breakText = 'Enjoy your break...';
    DrawFormattedText(ptbWindow, breakText, 'center', 'center', black);
end
Screen('Flip',ptbWindow);

%% Close Psychtoolbox window
Priority(0);
Screen('Close');
Screen('CloseAll');