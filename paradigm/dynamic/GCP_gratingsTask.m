%% GCP Gratings Task
%
% This code requires PsychToolbox. https://psychtoolbox.org
% This was tested with PsychToolbox version 3.0.15, and with MATLAB R2023b.

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
    disp('INITIALIZING EEG... PLEASE WAIT 10 SECONDS')
    for i=1:10
        waitbar(i/10, 'INITIALIZING EEG');
        pause(1);
    end
    wbar = findall(0,'type','figure','tag','TMWWaitbar');
    delete(wbar)
    disp('EEG INITIALIZED!')
end

% Hide cursor on participant screen
HideCursor(whichScreen);

%% Define TRIGGERS
TASK_START = 10; % trigger for ET cutting

BLOCK1                 = 11; % Trigger for start of block 1
BLOCK2                 = 12; % Trigger for start of block 2
BLOCK3                 = 13; % Trigger for start of block 3
BLOCK4                 = 14; % Trigger for start of block 4
BLOCK0                 = 15; % Trigger for start of training block (block 0)

FIXCROSSR              = 16; % Trigger for red (task) fixation cross
FIXCROSSB              = 17; % Trigger for black fixation cross

PRESENTATION_LC_TASK   = 51; % Trigger for presentation of low contrast concentric dynamic inward grating WITH button press response
PRESENTATION_HC_TASK   = 52; % Trigger for presentation of high contrast concentric dynamic inward grating WITH button press response
PRESENTATION_LC_NOTASK = 61; % Trigger for presentation of low contrast concentric dynamic inward grating WITHOUT button press response
PRESENTATION_HC_NOTASK = 62; % Trigger for presentation of high contrast concentric dynamic inward grating WITHOUT button press response

BLOCK1_END             = 71; % End of block 1
BLOCK2_END             = 72; % End of block 2
BLOCK3_END             = 73; % End of block 3
BLOCK4_END             = 74; % End of block 4
BLOCK0_END             = 75; % End of block 0

RESP_YES               = 87; % Trigger for response yes (spacebar)
RESP_NO                = 88; % Trigger for response no (no input)

TASK_END               = 90; % Trigger for ET cutting

%% Set up experiment parameters
% Block and Trial Number
exp.nTrlTrain = 10; % n gratings per training block
exp.nTrlTask = 125; % n gratings per task block

if TRAINING == 1
    exp.nTrials = exp.nTrlTrain;
else
    exp.nTrials = exp.nTrlTask;
end

% Set up equipment parameters
equipment.viewDist = 800;               % Viewing distance in millimetres
equipment.ppm = 3.6;                    % Pixels per millimetre !NEEDS TO BE SET using the MeasureDpi function!
equipment.greyVal = .5;
equipment.blackVal = 0;
equipment.whiteVal = 1;
equipment.gammaVals = [1 1 1];          % The gamma values for color calibration of the monitor

% Set up stimulus parameters Fixation
stimulus.fixationOn = 1;                  % Toggle fixation on (1) or off (0)
stimulus.fixationSize_dva = .3;           % Size of fixation cross in degress of visual orientationAngle
stimulus.fixationColor0 = [0 0 0];        % Color of fixation cross (1 = white, 0 = black, [1 0 0] = red)
stimulus.fixationColor1 = [1 0 0];        % Red fixation cross
stimulus.fixationLineWidth = 1.3;         % Line width of fixation cross

% Location
stimulus.regionHeight_dva = 7.3;         % Height of the region
stimulus.regionWidth_dva = 4;            % Width of the region
stimulus.regionEccentricity_dva = 3;     % Eccentricity of regions from central fixation

% Set up color parameters
color.textVal = 0;                      % Color of text (0 = black)
color.Black = color.textVal;
color.White = 1;

%  Retrieve key code for spacebar
spaceKeyCode = KbName('Space');

%% Set up text parameters
% Define startExperimentText
startExperimentText = [
    'You will see a series of gratings. \n\n' ...
    'Between gratings, a fixation cross \n\n' ...
    'will appear on the screen. If the \n\n' ...
    'fixation cross appers in RED, it is \n\n' ...
    'your task to press SPACE during the next \n\n' ...
    'grating. Use your right hand. Please \n\n' ...
    'always look at the center of the screen \n\n' ...
    '\n\n' ...
    'Press any key to continue.'];
if TRAINING == 1
    loadingText = 'Loading TRAINING...';
elseif TRAINING == 0
    loadingText = 'Loading TASK...';
end

%% Set up temporal parameters (all in seconds)
timing.blank = 1;                               % Duration of blank screen

% For Fixation Cross
timing.cfilower = 2000; % lower limit of CFI duration
timing.cfiupper = 3000; % upper limit of CFI duration
timing.cfi_task = 0.5;

% For Stimulus Trial
timing.stimuliDuration = 2;

% Shuffle rng for random elements
rng('default');
rng('shuffle');                     % Use MATLAB twister for rng

%% Set up Psychtoolbox Pipeline
global GL;
AssertOpenGL;

% Imaging set up
screenID = whichScreen;
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
Screen('Preference', 'SkipSyncTests', 0); % For linux (can be 0)

% Set verbosity to disallow CW output
Screen('Preference','Verbosity', 0);

% Window setup
[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, equipment.greyVal);
PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
[screenWidth, screenHeight] = RectSize(winRect);
screenCenterX = round(screenWidth/2);
screenCenterY = round(screenHeight/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);
Screen('BlendFunction', ptbWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
exp.runPriority = MaxPriority(ptbWindow);

% Set font size for instructions and stimuli
Screen('TextSize', ptbWindow, 36);

global psych_default_colormode;                     % Sets colormode to be unclamped   0-1 range.
psych_default_colormode = 1;

global ptb_drawformattedtext_disableClipping;       % Disable clipping of text
ptb_drawformattedtext_disableClipping = 1;

% Show loading text
DrawFormattedText(ptbWindow,loadingText,'center','center',color.textVal);
Screen('Flip',ptbWindow);

%% Calculate equipment parameters
equipment.mpd = (equipment.viewDist/2)*tan(deg2rad(2*stimulus.regionEccentricity_dva))/stimulus.regionEccentricity_dva; % Millimetres per degree
equipment.ppd = equipment.ppm*equipment.mpd;    % Pixels per degree

% Fix coordinates for fixation cross
stimulus.fixationSize_pix = round(stimulus.fixationSize_dva*equipment.ppd);
fixHorizontal = [round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2) 0 0];
fixVertical = [0 0 round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2)];
fixCoords = [fixHorizontal; fixVertical];
fixPos = [screenCenterX, screenCenterY];

%% Settings for inward moving circular grating
% Copied from DriftDemo and modified to display a (masked) animated concentric
% grating moving inward. Adapted from van Es and Schoffelen, 2019.
% https://github.com/Donders-Institute/dyncon_erfosc/blob/master/concentric_grating_experiment.m

% Find the color values which correspond to white and black(black = 0; white = 255)
white = WhiteIndex(ptbWindow);
black = BlackIndex(ptbWindow);
grey  = round((white+black)/2);

% Contrast 'inc'rement range for given white and grey values:
inc = white-grey;

% Open a double buffered fullscreen window and select a black background color:
[ptbWindow, windowRect]=Screen('OpenWindow',1, grey);

% Query the frame duration
ifi = Screen('GetFlipInterval', ptbWindow);
frameRate = Screen('FrameRate', 1); %1/ifi

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', ptbWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Grating size
equipment.ppd = screen.resolutionX/screen.totVisDeg;
visualAngleGrating = 7.1;
visualAngleLocation = 15;
gratingSize = visualAngleGrating*equipment.ppd;
gratingRadius = round(gratingSize/2); % Grating can only exist of integers -> round
gratingSize = 2*gratingRadius; % To prevent consistency errors, redifine gratingSize
% rLocation = round(visualAngleLocation*equipment.ppd/2);

% Frequency
driftFreq = 2; % Every pixel of the grating completes two cycles per second (black-white-black)
nFramesInCycle = round((1/driftFreq)/ifi); % Temporal period, in frames, of the drifting grating

% Generate stimulus
[x,y] = meshgrid(-gratingRadius:gratingRadius,-gratingRadius:gratingRadius);
f     = 0.55*2*pi; % period of the grating.

% Circular hanning mask
L                     = 2*gratingRadius+1;
w1D                   = hann(L); % 1D hann window
xx                    = linspace(-gratingRadius,gratingRadius,L);
[X,Y]                 = meshgrid(xx);
r                     = sqrt( X.^2 + Y.^2 );
w2D                   = zeros(L);
w2D(r<=gratingRadius) = interp1(xx,w1D,r(r<=gratingRadius)); % 2D hanning window

% Generate grating texture
% Compute each frame of the movie and convert those frames stored in
% MATLAB matrices, into Psychtoolbox OpenGL textures using 'MakeTexture'
tex             = zeros(nFramesInCycle,1);
for jFrame      = 1:nFramesInCycle
    phase       = (jFrame / nFramesInCycle) * 2 * pi; % Change the phase of the grating according to the framenumber
    m           = sin(sqrt(x.^2+y.^2) / f + phase); % Formula sinusoidal
    grating     = (w2D.*(inc*m)+grey);
    % inc*m fluctuates from [-grey, grey]. Multiply this with the
    % hanning mask to let the grating die off at 0
    tex(jFrame) = Screen('MakeTexture', ptbWindow, grating);
end

% Set location
gratingDim = [0 0 2*gratingRadius 2*gratingRadius];
gratingYpos = screenCenterY;
gratingXpos = screenCenterX;
frameTexId = mod(0:(nFramesTotal-1), nFramesInCycle) + 1; % Assign the right texture index to each frame
position = CenterRectOnPointd(gratingDim, gratingXpos, gratingYpos); % Move the object to those coordinates

% Set duration
movieDurationSecs = 2;
nFramesTotal = round(movieDurationSecs * frameRate); % Convert movieDuration in seconds to duration in frames

% Use realtime priority for better timing precision
priorityLevel = MaxPriority(ptbWindow);
Priority(priorityLevel);

disp(char('INWARD MOVING CONCENTRIC GRATING SETUP COMPLETED'));

%% Create data structure for preallocating data
data                             = struct;
nums                             = repmat(1:2, 1, 100);
gratingSequence                  = nums(randperm(length(nums), exp.nTrials)); % Define grating sequence
data.grating(1, exp.nTrials)     = NaN; % Saves grating form (1 or 2, see below)
% grating = 1 is low contrast concentric dynamic inward
% grating = 2 is high contrast concentric dynamic inward
data.contrast(1, exp.nTrials)    = NaN; % Binary measure for low/high contrast
data.redCross(1, exp.nTrials)    = NaN; % Binary measure for task condition
data.responses(1, exp.nTrials)   = NaN; % Binary measure for (no) response
data.correct(1, exp.nTrials)     = NaN; % Binary measure for correct responses
data.reactionTime(1:exp.nTrials) = NaN; % Reaction time
data.fixation(1:exp.nTrials)     = NaN; % Fixation check info
data.trlDuration(1:exp.nTrials)  = NaN; % Trial duration in seconds

%% Show task instruction text
DrawFormattedText(ptbWindow,startExperimentText,'center','center',color.textVal);
Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1); % black background for photo diode
Screen('Flip',ptbWindow);
clc;
disp('Participant is reading the instructions.');
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

HideCursor(whichScreen);
timing.startTime = datestr(now, 'dd/mm/yy-HH:MM:SS');
count5trials = 0;

%% Experiment Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
disp('GCP GRATING TASK...');
for trl = 1:exp.nTrials
    tic;
    % Add trial grating info
    if gratingSequence(trl) == 1
        gratingForm = 'low contrast concentric dynamic inward';
    elseif gratingSequence(trl) == 2
        gratingForm = 'high contrast concentric dynamic inward';
    end

    % Randomized selection of task (red fication cross) trials (25%)
    if randi(4) == 1
        data.redCross(trl) = 1;
    else
        data.redCross(trl) = 0;
    end

    %% Present fixation cross (red for task condition)
    % Fill gray screen
    Screen('FillRect', ptbWindow, gray);
    Screen('Flip', ptbWindow);
    % Set jittered trial-specific durations for CFIs
    timing.cfi(trl) = (randsample(timing.cfilower:timing.cfiupper, 1))/1000; % Randomize the jittered central fixation interval on trial
    start_time = GetSecs;
    while (GetSecs - start_time) < timing.cfi(trl)
        if data.redCross(trl) == 0 % No task condition
            Screen('DrawLines', ptbWindow, fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor0,[screenCenterX screenCenterY],2);
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
            WaitSecs(timing.cfi(trl));
        elseif data.redCross(trl) == 1 % Task condition
            Screen('DrawLines', ptbWindow, fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor1,[screenCenterX screenCenterY],2);
            Screen('Flip', ptbWindow);
            TRIGGER = FIXCROSSR;
            if TRAINING == 1
                Eyelink('Message', num2str(TRIGGER));
                Eyelink('command', 'record_status_message "FIXCROSS"');
            else
                Eyelink('Message', num2str(TRIGGER));
                Eyelink('command', 'record_status_message "FIXCROSS"');
                sendtrigger(TRIGGER,port,SITE,stayup);
            end
            WaitSecs(timing.cfi_task); % Show red cross for 500 ms
            Screen('DrawLines', ptbWindow, fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor0,[screenCenterX screenCenterY],2);
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
    % grating =  1 is low contrast concentric dynamic inward
    % grating =  2 is high contrast concentric dynamic inward

    %% Present grating and get response
    Screen('Flip', ptbWindow); % Preparatory flip
    responseGiven = false;
    maxProbeDuration = 2; % Maximum time to show the grating

    % Send presentation triggers
    if gratingSequence(trl) == 1 && data.redCross(trl) == 1
        TRIGGER = PRESENTATION_LC_TASK;
    elseif gratingSequence(trl) == 2 && data.redCross(trl) == 1
        TRIGGER = PRESENTATION_HC_TASK;
    elseif gratingSequence(trl) == 1 && data.redCross(trl) == 0
        TRIGGER = PRESENTATION_LC_NOTASK;
    elseif gratingSequence(trl) == 2 && data.redCross(trl) == 0
        TRIGGER = PRESENTATION_HC_NOTASK;
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
    % Draw gratings depending on gratingSequence
    while (GetSecs - probeStartTime) < maxProbeDuration
        % Calculate the current phase based on the elapsed time
        elapsedTime = GetSecs - startTime;
        phase = (elapsedTime * driftFreq) * 2 * pi; % Adjust frequency

        % Generate the updated grating texture based on the current phase
        grating = (w2D .* (inc * sin(sqrt(x.^2 + y.^2) / f + phase)) + grey);

        % Create the new texture for the current frame
        tex = Screen('MakeTexture', ptbWindow, grating);

        if gratingSequence(trl) == 1 % low contrast concentric dynamic inward
            Screen('DrawTexture', ptbWindow, tex, [], position);
            Screen('Flip', ptbWindow);
        elseif gratingSequence(trl) == 2 % high contrast concentric dynamic inward
            Screen('DrawTexture', ptbWindow, tex, [], position);
            Screen('Flip', ptbWindow);
        end

        % Take screenshot of current screen
        screenshotFilename = sprintf('GCP_screenshot_%s.png', gratingForm);
        imageArray = Screen('GetImage', ptbWindow);
        imwrite(imageArray, screenshotFilename);

        % Check for participant response
        if ~responseGiven
            [keyIsDown, responseTime, keyCode] = KbCheck;
            if keyIsDown
                responseGiven = true;
                data.reactionTime(trl) = responseTime - probeStartTime;
                data.responses(trl) = 1; % Response made
            end
        end
    end

    % If no response is given, record default
    if ~responseGiven
        data.responses(trl) = 0; % No response
        data.reactionTime(trl) = NaN;
    end

    %% Check if response was correct
    if data.redCross(trl) == 1 && data.responses(trl) == 1 % Red fixation cross + button press = correct
        data.correct(trl) = 1;
        feedbackText = 'Correct!  ';
    elseif data.redCross(trl) == 0 && data.responses(trl) == 0 % No red fixation cross + no button press = correct
        data.correct(trl) = 1;
        feedbackText = 'Correct!  ';
    else % Anything else is wrong response
        data.correct(trl) = 0;
        feedbackText = 'Incorrect!';
    end

    %% Feedback for training block and CW output
    % Give feedback in training block
    if TRAINING == 1
        DrawFormattedText(ptbWindow,feedbackText,'center','center',color.Black);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('Flip',ptbWindow);
        WaitSecs(2);
        % Give feedback for no response (too slow)
    elseif TRAINING == 0 && data.correct(trl) == 0 && data.responses(trl) == 0
        feedbackText = 'TOO SLOW! ';
        DrawFormattedText(ptbWindow,feedbackText,'center','center',color.Black);
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
            DrawFormattedText(ptbWindow,feedbackLastTrials,'center','center',color.White);
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
            ' in Block ' num2str(BLOCK) ' is ' feedbackText ' (Red FixCross: ' ...
            '' num2str(data.redCross(trl)) ' | Acc: ' num2str(overall_accuracy) ...
            '% | RT: ' reactionTime 's | ' gratingForm ')']);
    else
        disp(['Response to Trial ' num2str(trl) '/' num2str(exp.nTrials) ...
            ' in Block ' num2str(BLOCK) ' is ' feedbackText '(Red FixCross: ' ...
            num2str(data.redCross(trl)) ' | Acc: ' num2str(overall_accuracy) ...
            '% | RT: ' reactionTime 's | ' gratingForm ')']);
    end
    % Save duration in seconds
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
saves.screenWidth               = screenWidth;
saves.screenHeight              = screenHeight;
saves.screenCenterX             = screenCenterX;
saves.screenCenterY             = screenCenterY;
saves.startExperimentText       = startExperimentText;
saves.stimulus                  = stimulus;
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

trigger.PRESENTATION_LC_TASK    = PRESENTATION_LC_TASK;
trigger.PRESENTATION_HC_TASK    = PRESENTATION_HC_TASK;
trigger.PRESENTATION_LC_NOTASK  = PRESENTATION_LC_NOTASK;
trigger.PRESENTATION_HC_NOTASK  = PRESENTATION_HC_NOTASK;

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

DrawFormattedText(ptbWindow,breakInstructionText,'center','center',color.textVal);
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
    DrawFormattedText(ptbWindow, FinalText, 'center', 'center', color.textVal);
elseif BLOCK == 1 && TRAINING == 0 || BLOCK == 2 && TRAINING == 0 || BLOCK == 3 && TRAINING == 0
    breakText = 'Enjoy your break...';
    DrawFormattedText(ptbWindow, breakText, 'center', 'center', color.textVal);
end

Screen('Flip',ptbWindow);

%% Close Psychtoolbox window
Priority(0);
Screen('Close');
Screen('CloseAll');