%% GCP Gratings Task
%
% This code requires PsychToolbox. https://psychtoolbox.org
% This was tested with PsychToolbox version 3.0.15, and with MATLAB R2023b.

%% Initialize EEG and ET
% Calibrate ET (Tobii Pro Fusion)
disp('CALIBRATING ET...');
calibrateET;

if TRAINING == 0
    % Start recording EEG
    disp('STARTING EEG RECORDING...');
    initEEG;

    %     % Wait ten seconds to initialize EEG
    %     disp('INITIALIZING EEG... PLEASE WAIT 10 SECONDS')
    %     for i=1:15
    %         waitbar(i/15, 'INITIALIZING EEG');
    %         pause(1);
    %     end
    %     wbar = findall(0,'type','figure','tag','TMWWaitbar');
    %     delete(wbar)
    %     disp('EEG INITIALIZED!')
end
TRAINING == 1

% Hide cursor on participant screen
HideCursor(whichScreen);

%% Define TRIGGERS
TASK_START = 10; % trigger for ET cutting

BLOCK1 = 11; % trigger for start of block 1
BLOCK2 = 12;
BLOCK3 = 13;
BLOCK4 = 14;
BLOCK0 = 15; % trigger for start of training block (block 0)

FIXATION = 15; % trigger for fixation cross

PRESENTATION1 = 21; % trigger for presentation of low contrast horizontal
PRESENTATION2 = 22; % trigger for presentation of high contrast horizontal
PRESENTATION3 = 31; % trigger for presentation of low contrast vertical
PRESENTATION4 = 32; % trigger for presentation of high contrast vertical
PRESENTATION5 = 41; % trigger for presentation of low contrast concentric
PRESENTATION6 = 42; % trigger for presentation of high contrast concentric

BLOCK1_END = 71; % end of block 1
BLOCK2_END = 72;
BLOCK3_END = 73;
BLOCK4_END = 74;
BLOCK0_END = 75; % end of block 0

RESP_YES = 87; % trigger for response yes (spacebar)
RESP_NO = 88; % trigger for response no (no input)

TASK_END = 90; % trigger for ET cutting

%% Set up experiment parameters
% Block and Trial Number
blocks = 4;
exp.nGratingsTrain = 8; % n gratings per training block
exp.nGratingsTask = 30; % n gratings per test block 4 blocks

if TRAINING == 1
    exp.nTrials = exp.nGratingsTrain;
else
    exp.nTrials = exp.nGratingsTask;
end

% Set up equipment parameters
equipment.viewDist = 800;               % Viewing distance in millimetres
equipment.ppm = 3.6;                    % Pixels per millimetre !! NEEDS TO BE SET. USE THE MeasureDpi FUNCTION !!
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
if TRAINING == 1
    loadingText = 'Loading training task...';
    startExperimentText = ['Training task: \n\n' ...
        'You will see a series of gratings. \n\n' ...
        'Your task is to press SPACE if you see \n\n' ...
        'a red fixation cross appearing for a short time. \n\n' ...
        'Otherwise, do not press any key. \n\n' ...
        'Please always look at the screen center \n\n' ...
        'and use your right hand. \n\n' ...
        'Press any key to continue.'];
elseif TRAINING == 0
    loadingText = 'Loading test task...';
    startExperimentText = [
        'You will see a series of gratings. \n\n' ...
        'Your task is to press SPACE if you see \n\n' ...
        'a red fixation cross appearing for a short time. \n\n' ...
        'Otherwise, do not press any key. \n\n' ...
        'Please always look at the screen center \n\n' ...
        'and use your right hand. \n\n' ...
        'Press any key to continue.'];
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
screenCentreX = round(screenWidth/2);
screenCentreY = round(screenHeight/2);
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

% Fix coordiantes for fixation cross
stimulus.fixationSize_pix = round(stimulus.fixationSize_dva*equipment.ppd);
fixHorizontal = [round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2) 0 0];
fixVertical = [0 0 round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2)];
fixCoords = [fixHorizontal; fixVertical];
fixPos = [screenCentreX, screenCentreY];
virtualSize = 400;  % Default x + y size

%% Set circular grating settings
rshader = [PsychtoolboxRoot 'PsychDemos/ExpandingRingsShader.vert.txt'];
expandingRingShader = LoadGLSLProgramFromFiles({ rshader, [PsychtoolboxRoot 'PsychDemos/ExpandingRingsShader.frag.txt'] }, 1);
ringwidth = 8; % Width of a single ring (radius) in pixels
ringtex = Screen('SetOpenGLTexture', ptbWindow, [], 0, GL.TEXTURE_RECTANGLE_EXT, virtualSize, virtualSize, 1, expandingRingShader);
glUseProgram(expandingRingShader);
glUniform2f(glGetUniformLocation(expandingRingShader, 'RingCenter'), virtualSize/2, virtualSize/2); % Center ring
glUseProgram(0); % Done with setup, disable shader

%% Create data structure for preallocating data
data = struct;
% Define grating sequence
nums = repmat(1:6, 1, 10);
gratingSequence = nums(randperm(length(nums), exp.nTrials));
data.grating(1, exp.nTrials) = NaN;
% grating = 1 is low contrast horizontal
% grating = 2 is high contrast horizontal
% grating = 3 is low contrast vertical
% grating = 4 is high contrast vertical
% grating = 5 is low contrast concentric
% grating = 6 is high contrast concentric
data.contrast(1, exp.nTrials) = NaN;
data.redCross(1, exp.nTrials) = NaN;
data.responses(1, exp.nTrials) = NaN;
data.correct(1, exp.nTrials) = NaN;
data.reactionTime(1:exp.nTrials) = NaN;

%% Show task instruction text
DrawFormattedText(ptbWindow,startExperimentText,'center','center',color.textVal);
Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1); % black background for photo diode
Screen('Flip',ptbWindow);
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
    disp('Start of Block 0 (Training)');
else
    disp(['Start of Block ' num2str(BLOCK)]);
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
tic;
timing.startTime =  datestr(now, 'dd/mm/yy-HH:MM');
count5trials = 0;

%% Experiment Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('GRATING TASK...');
for trl = 1:exp.nTrials

    if randi(2) == 1
        gratingSequence(trl) = 5
    else
        gratingSequence(trl) = 6
    end










    if gratingSequence(trl) == 1
        gratingForm = 'low contrast horizontal';
    elseif gratingSequence(trl) == 2
        gratingForm = 'high contrast horizontal';
    elseif gratingSequence(trl) == 3
        gratingForm = 'low contrast vertical';
    elseif gratingSequence(trl) == 4
        gratingForm = 'high contrast vertical';
    elseif gratingSequence(trl) == 5
        gratingForm = 'low contrast concentric';
    elseif gratingSequence(trl) == 6
        gratingForm = 'high contrast concentric';
    end

    % Randomized selection of task trials (33%)
    if randi(3) == 1
        data.redCross(trl) = 1;
    else
        data.redCross(trl) = 0;
    end

    % Start of trial CW output
    disp(['Start of Trial ' num2str(trl) '/' num2str(exp.nTrials) ' in Block ' num2str(BLOCK) ' (' gratingForm ')']);

    %% Present fixation cross
    % Set jittered trial-specific durations for CFIs
    timing.cfi(trl) = (randsample(timing.cfilower:timing.cfiupper, 1))/1000; % Randomize the jittered central fixation interval on trial
    start_time = GetSecs;
    while (GetSecs - start_time) < timing.cfi(trl)
        Screen('DrawLines', ptbWindow, fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor0,[screenCentreX screenCentreY],2);
        Screen('Flip', ptbWindow)
        TRIGGER = FIXATION;
        if TRAINING == 1
            Eyelink('Message', num2str(TRIGGER));
            Eyelink('command', 'record_status_message "FIXCROSS"');
        else
            Eyelink('Message', num2str(TRIGGER));
            Eyelink('command', 'record_status_message "FIXCROSS"');
            sendtrigger(TRIGGER,port,SITE,stayup);
        end
        WaitSecs(timing.cfi(trl))
    end

    %% Define grating depending on sequence number
    data.grating(trl) = gratingSequence(trl);
    % Grating = 1 is low contrast horizontal
    % Grating = 2 is high contrast horizontal
    % Grating = 3 is low contrast vertical
    % Grating = 4 is high contrast vertical
    % Grating = 5 is low contrast concentric
    % Grating = 6 is high contrast concentric
    % Common Settings for parameters in DrawTexture
    phase = 0; % Start phase
    frequency = 0.05; % Radial frequency (adjust for ring spacing)
    contrast = 0.5; % High contrast for concentric pattern
    sigma = 1.0; % < 0 is sinusiodal; > 0 is square-wave grating
    radius = floor(virtualSize / 2); % Radius of the disc edge
    PsychDefaultSetup(2); % Setup defaults and unit color range

    %% Linear grating (horizontal, vertical) setup
    if gratingSequence(trl) < 5

        % Set orientationAngle
        if gratingSequence(trl) == 1 || gratingSequence(trl) == 2
            orientationAngle = 90;   % 90° = horizontal
        elseif gratingSequence(trl) == 3 || gratingSequence(trl) == 4
            orientationAngle = 0;    % 0° = vertical
        end

        % Color Settings
        if gratingSequence(trl) == 1 || gratingSequence(trl) == 3 % Low contrast
            color1 = [0.25 0.25 0.25 1];
            color2 = [0.75 0.75 0.75 1];
            contrastPreMultiplicator = 0.5;
            data.contrast(trl) = 0;
        elseif gratingSequence(trl) == 2 || gratingSequence(trl) == 4 % High contrast
            color1 = [0 0 0 1];
            color2 = [1 1 1 1];
            contrastPreMultiplicator = 1;
            data.contrast(trl) = 1;
        end

        % Query frame duration for later use to time 'Flips' properly for an
        % animation with constant framerate:
        ifi = Screen('GetFlipInterval', ptbWindow);

        % Create linear grating texture
        evalc(['[texture, gratingRect] = CreateProceduralSquareWaveGrating(ptbWindow, ' ...
            'virtualSize, virtualSize, [0 0 0 1], radius, contrastPreMultiplicator)']);
    end

    %% Concentric grating setup
    if gratingSequence(trl) == 5 || gratingSequence(trl) == 6
        % Retrieve monitor refresh duration
        ifi = Screen('GetFlipInterval', ptbWindow);

        % Concentric colour settings
        if gratingSequence(trl) == 5 % Low contrast (grey scale)
            firstColor   = [0.75 0.75 0.75 1]; % Dark grey
            secondColor  = [0.25 0.25 0.25 1]; % Light grey
            data.contrast(trl) = 0; % Low contrast
        elseif gratingSequence(trl) == 6 % High contrast (black and white)
            firstColor   = [1 1 1 1]; % Black
            secondColor  = [0 0 0 1]; % White
            data.contrast(trl) = 1; % High contrast
        end
        shiftvalue = 0;
    end


    %% Present grating and get response
    % Suppress output from shader creation
    vbl = Screen('Flip', ptbWindow); % Preparatory flip
    responseGiven = false;
    probeStartTime = GetSecs;
    maxProbeDuration = 2; % Maximum time to show the grating

    while (GetSecs - probeStartTime) < maxProbeDuration
        % Draw grating
        if gratingSequence(trl) < 5
                dstRect = CenterRectOnPointd(gratingRect, screenCentreX, screenCentreY);
            Screen('DrawTexture', ptbWindow, texture, [], dstRect, orientationAngle, [], [], [], [], [], [phase, frequency, contrast, sigma]);
        else
%             Screen('DrawTexture', ptbWindow, texture, [], dstRect, 0, [], [], [], [], [], [phase, frequency, contrast, sigma]);
        Screen('DrawTexture', ptbWindow, ringtex, [], [], [], [], [], ...
            firstColor, [], [], [secondColor(1), secondColor(2), secondColor(3), ...
            secondColor(4), shiftvalue, ringwidth, radius, 0]);
        
        end

        % Send presentation triggers
        if gratingSequence(trl) == 1
            TRIGGER = PRESENTATION1;
        elseif gratingSequence(trl) == 2
            TRIGGER = PRESENTATION2;
        elseif gratingSequence(trl) == 3
            TRIGGER = PRESENTATION3;
        elseif gratingSequence(trl) == 4
            TRIGGER = PRESENTATION4;
        elseif gratingSequence(trl) == 5
            TRIGGER = PRESENTATION5;
        elseif gratingSequence(trl) == 6
            TRIGGER = PRESENTATION6;
        end
        if TRAINING == 1
            Eyelink('Message', num2str(TRIGGER));
            Eyelink('command', 'record_status_message "PRESENTATION"');
        else
            Eyelink('Message', num2str(TRIGGER));
            Eyelink('command', 'record_status_message "PRESENTATION"');
            sendtrigger(TRIGGER,port,SITE,stayup);
        end

        % Draw fixation cross on top of the grating
        if data.redCross(trl) == 1 % Red fixation cross (task condition)
            Screen('DrawLines', ptbWindow, fixCoords, stimulus.fixationLineWidth, stimulus.fixationColor1, [screenCentreX screenCentreY], 2);
        else % Normal fixation cross
            Screen('DrawLines', ptbWindow, fixCoords, stimulus.fixationLineWidth, stimulus.fixationColor0, [screenCentreX screenCentreY], 2);
        end
        vbl = Screen('Flip', ptbWindow, vbl + ifi);

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
        %     elseif TRAINING == 0 && data.correct(trl) == 0 && data.responses(trl) == 0
        %         feedbackText = 'TOO SLOW! ';
        %         DrawFormattedText(ptbWindow,feedbackText,'center','center',color.Black);
        %         Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        %         Screen('Flip',ptbWindow);
        %         WaitSecs(2);
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
        disp(['Response to Trial ' num2str(trl) '/' num2str(exp.nTrials) ' in Block ' num2str(BLOCK) ' is ' feedbackText '  (Acc: ' num2str(overall_accuracy) '% | RT: ' reactionTime 's)']);
    else
        disp(['Response to Trial ' num2str(trl) '/' num2str(exp.nTrials) ' in Block ' num2str(BLOCK) ' is ' feedbackText ' (Acc: ' num2str(overall_accuracy) '% | RT: ' reactionTime 's)']);
    end
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

% Record block duration
timing.duration = toc;
timing.endTime = datestr(now, 'dd/mm/yy-HH:MM');

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
saves = struct;
saves.data = data;
saves.data.spaceKeyCode = spaceKeyCode;
saves.data.reactionTime = reactionTime;
saves.experiment = exp;
saves.screenWidth = screenWidth;
saves.screenHeight = screenHeight;
saves.screenCentreX = screenCentreX;
saves.screenCentreY = screenCentreY;
saves.startExperimentText = startExperimentText;
saves.stimulus = stimulus;
saves.subjectID = subjectID;
saves.subject = subject;
saves.timing = timing;

% Save triggers
trigger = struct;
trigger.TASK_START = 10;
trigger.BLOCK1 = BLOCK1;
trigger.BLOCK2 = BLOCK2;
trigger.BLOCK3 = BLOCK3;
trigger.BLOCK4 = BLOCK4;
trigger.BLOCK0 = BLOCK0;

trigger.FIXATION = FIXATION;

trigger.PRESENTATION1 = PRESENTATION1;
trigger.PRESENTATION2 = PRESENTATION2;
trigger.PRESENTATION3 = PRESENTATION3;
trigger.PRESENTATION4 = PRESENTATION4;
trigger.PRESENTATION5 = PRESENTATION5;
trigger.PRESENTATION6 = PRESENTATION6;

trigger.BLOCK1_END = BLOCK1_END;
trigger.BLOCK2_END = BLOCK2_END;
trigger.BLOCK3_END = BLOCK3_END;
trigger.BLOCK4_END = BLOCK4_END;
trigger.BLOCK0_END = BLOCK0_END;

trigger.RESP_YES = RESP_YES;
trigger.RESP_NO = RESP_NO;

trigger.TASK_END = 90;

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
        '\n\n Press any key to start the mandatory break of at least 15 seconds.'];
end

DrawFormattedText(ptbWindow,breakInstructionText,'center','center',color.textVal);
Screen('Flip',ptbWindow);
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

% Create final screen
if BLOCK == 4 && TRAINING == 0
    FinalText = ['You are done.' ...
        '\n\n Have a great day!'];
    DrawFormattedText(ptbWindow, FinalText, 'center', 'center', color.textVal);
elseif BLOCK == 1 && TRAINING == 0 || BLOCK == 2 && TRAINING == 0 || BLOCK == 3 && TRAINING == 0
    BreakText = 'Enjoy your break...';
    DrawFormattedText(ptbWindow, BreakText, 'center', 'center', color.textVal);
end

Screen('Flip',ptbWindow);

%% Close Psychtoolbox window
Screen('CloseAll');