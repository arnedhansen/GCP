%% Gratings
%
% This script executes the paradigm with gatings of Lea's MA. 
%
% This code requires PsychToolbox. https://psychtoolbox.org
% This was tested with PsychToolbox version 3.0.15, and with MATLAB R2022a.

%% Summary
% 1. Definition of task, trial number / durations
% 2. Block for loop: definition of triggers, text properties and introductory text
% 3. Trial for loop: display fixation cross (20% red) and stimuli
% 4. End of trial for loop: save trial data
% 5. End of block for loop: save block data


%% TASK

% Check if mat files exist, in order to avoid overwriting these and start
% at correct block in case of breakdown
if TRAINING == 0
    if isfile([DATA_PATH, '/', num2str(subject.ID), '/', [num2str(subject.ID), '_G_block4.mat']])
        return
    elseif isfile([DATA_PATH, '/', num2str(subject.ID), '/', [num2str(subject.ID), '_G_block3.mat']])
        disp('BLOCK 3 DATA ALREADY EXISTS');
        start = 4;
    elseif isfile([DATA_PATH, '/', num2str(subject.ID), '/', [num2str(subject.ID), '_G_block2.mat']])
        disp('BLOCK 2 DATA ALREADY EXISTS');
        start = 3;
    elseif isfile([DATA_PATH, '/', num2str(subject.ID), '/', [num2str(subject.ID), '_G_block1.mat']])
        disp('BLOCK 1 DATA ALREADY EXISTS');
        start = 2;
    else
        start = 1;
    end
elseif TRAINING == 1
    start = 1;
end


%% Define trial number and duration

% Block and Trial Number
NumberOfBlocks = 4;
experiment.nGratingsTrain = 8; % n gratings per train block, normally 8
experiment.nGratingsTest = 125; % n gratings per test block, normally 125

if TRAINING == 1
    experiment.nGratings = experiment.nGratingsTrain;
elseif TRAINING == 0
    experiment.nGratings = experiment.nGratingsTest;
end

% For Fixation Cross
timing.cfilower = 2000; % lower limit of CFI duration
timing.cfiupper = 3000; % upper limit of CFI duration
timing.cfi_task = 0.5;

% For Stimulus Trial
timing.StimuliDuration = 2; 


%% For loop

for BLOCK = start : NumberOfBlocks
    % Start the actual task (EEG recording will start here, if TRAINING = 0)
    disp('GRATING TASK...');
    WaitSecs(10);


    %% EEG and ET
    if TRAINING == 0
        % Start recording EEG
        disp('STARTING EEG RECORDING...');
        initEEG;
        WaitSecs(10);
    end

    % Calibrate ET (Tobii Pro Fusion)
    disp('CALIBRATING ET...');
    calibrateET

    %% Task
    HideCursor(whichScreen);

    % define triggers
    BLOCK1 = 11; % trigger for start of task (block)
    BLOCK2 = 12;
    BLOCK3 = 13;
    BLOCK4 = 14;
    BLOCK0 = 15; % trigger for start of training block
    
    FIXATION0 = 16; % trigger for fixation cross
    FIXATION1 = 17; % trigger for fixation cross
    PRESENTATION1 = 21; % trigger for presentation of high horizontal
    PRESENTATION2 = 22; % trigger for presentation of high vertical
    PRESENTATION3 = 23; % trigger for presentation of high 45
    PRESENTATION4 = 24; % trigger for presentation of high 115
    PRESENTATION5 = 25; % trigger for presentation of low horizontal
    PRESENTATION6 = 26; % trigger for presentation of low vertical
    PRESENTATION7 = 27; % trigger for presentation of low 45
    PRESENTATION8 = 28; % trigger for presentation of low 115
    
    RESP = 87; % trigger for response yes (spacebar)
    NO_RESP = 88; % trigger for response no (no input)
    
    BLOCK1_END = 91;
    BLOCK2_END = 92;
    BLOCK3_END = 93;
    BLOCK4_END = 94;
    BLOCK0_END = 95;

    % Set up experiment parameters
    % Number of trials for the experiment
    if TRAINING == 1
        experiment.nTrials = experiment.nGratingsTrain * 2;   % for each grating video, there should be a fixation cross => hence nTrial should be nGratings times two
    else
        experiment.nTrials = experiment.nGratingsTest * 2;   % for each grating video, there should be a fixation cross => hence nTrial should be nGratings times two
    end

    % Set up equipment parameters
    equipment.viewDist = 800;               % Viewing distance in millimetres
    equipment.ppm = 3.6;                    % Pixels per millimetre !! NEEDS TO BE SET. USE THE MeasureDpi FUNCTION !!
    equipment.greyVal = .5;
    equipment.blackVal = 0;
    equipment.whiteVal = 1;
    equipment.gammaVals = [1 1 1];          % The gamma values for color calibration of the monitor

    % Set up stimulus parameters Fixation
    stimulus.fixationOn = 1;                % Toggle fixation on (1) or off (0)
    stimulus.fixationSize_dva = .3;         % Size of fixation cross in degress of visual angle
    stimulus.fixationColor0 = [0 0 0];       % Color of fixation cross (1 = white, 0 = black, [1 0 0] = red)
    stimulus.fixationColor1 = [1 0 0];
    stimulus.fixationLineWidth = 1.3;         % Line width of fixation cross

    % Location
    stimulus.regionHeight_dva = 7.3;         % Height of the region
    stimulus.regionWidth_dva = 4;            % Width of the region
    stimulus.regionEccentricity_dva = 3;     % Eccentricity of regions from central fixation

    % Set up color parameters
    color.textVal = 0;                      % Color of text (0 = black)


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

    % Set up temporal parameters (all in seconds)
    timing.blank = 1;                               % Duration of blank screen

    % Shuffle rng for random elements
    rng('default');
    rng('shuffle');                     % Use MATLAB twister for rng

    % Set up Psychtoolbox Pipeline
    AssertOpenGL;

    % Imaging set up
    screenID = whichScreen;
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
    Screen('Preference', 'SkipSyncTests', 0); % For lin.ux

    % Window set-up
    [ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, equipment.greyVal);
    PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
    [screenWidth, screenHeight] = RectSize(winRect);
    screenCentreX = round(screenWidth/2);
    screenCentreY = round(screenHeight/2);
    flipInterval = Screen('GetFlipInterval', ptbWindow);
    Screen('BlendFunction', ptbWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    experiment.runPriority = MaxPriority(ptbWindow);

    % Set font size for instructions and stimuli
    Screen('TextSize', ptbWindow, 36);

    global psych_default_colormode;                     % Sets colormode to be unclamped   0-1 range.
    psych_default_colormode = 1;

    global ptb_drawformattedtext_disableClipping;       % Disable clipping of text
    ptb_drawformattedtext_disableClipping = 1;

    % Show loading text
    DrawFormattedText(ptbWindow,loadingText,'center','center',color.textVal);
    Screen('Flip',ptbWindow);

    % Retrieve response key
    spaceKeyCode = KbName('Space'); % Retrieve key code for spacebar

    % Calculate equipment parameters
    equipment.mpd = (equipment.viewDist/2)*tan(deg2rad(2*stimulus.regionEccentricity_dva))/stimulus.regionEccentricity_dva; % Millimetres per degree
    equipment.ppd = equipment.ppm*equipment.mpd;    % Pixels per degree

    % Fix coordiantes for fixation cross
    stimulus.fixationSize_pix = round(stimulus.fixationSize_dva*equipment.ppd);
    fixHorizontal = [round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2) 0 0];
    fixVertical = [0 0 round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2)];
    fixCoords = [fixHorizontal; fixVertical];

    % Create data structure for preallocating data
    data = struct;
    data.trialMatch(1:experiment.nGratings) = NaN;
    data.allResponses(1:experiment.nGratings) = NaN;
    data.allCorrect(1:experiment.nGratings) = NaN;

    % Preallocate reaction time variable
    reactionTime(1:experiment.nGratings) = 0;

    % Preallocate cfi timing variable
    timing.cfi(1:experiment.nGratings) = 0;

    % Define grating and cross sequences
    Shuffle;
    if TRAINING == 1 && BLOCK == 1
        gratingSequence = gratingSequence0;
    elseif TRAINING == 0 && BLOCK == 1
            gratingSequence = gratingSequence1;
    elseif TRAINING == 0 && BLOCK == 2
        gratingSequence = gratingSequence2;
    elseif TRAINING == 0 && BLOCK == 3
        gratingSequence = gratingSequence3;
    elseif TRAINING == 0 && BLOCK == 4
        gratingSequence = gratingSequence4;
    end

    if TRAINING == 1 && BLOCK == 1
        crossSequence = crossSequence0;
    elseif TRAINING == 0 && BLOCK == 1
        crossSequence = crossSequence1;
    elseif TRAINING == 0 && BLOCK == 2
        crossSequence = crossSequence2;
    elseif TRAINING == 0 && BLOCK == 3
        crossSequence = crossSequence3;
    elseif TRAINING == 0 && BLOCK == 4
        crossSequence = crossSequence4;
    end

    % Save grating sequence in mat file
    data.GratingSequence = gratingSequence;
    data.CrossSequence = crossSequence;


    %% Show task instruction text
    DrawFormattedText(ptbWindow,startExperimentText,'center','center',color.textVal);
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1); % black background for photo diode
    startExperimentTime = Screen('Flip',ptbWindow);
    disp('Participant is reading the instructions.');
    waitResponse = 1;
    while waitResponse
        [time, keyCode] = KbWait(-1,2);
        waitResponse = 0;
    end

%     % Send triggers: task starts. If training, send only ET triggers
%     if TRAINING == 1
%         %     EThndl.sendMessage(TASK_START); % ET
%         Eyelink('Message', num2str(BLOCK0));
%         Eyelink('command', 'record_status_message "TASK_START"');
%     else
%         %     EThndl.sendMessage(TASK_START); % ET
%         Eyelink('Message', num2str(TASK_START));
%         Eyelink('command', 'record_status_message "TASK_START"');
%         sendtrigger(TASK_START,port,SITE,stayup); % EEG
%     end

    % Send triggers for block and output
    if TRAINING == 1 % Training condition
        TRIGGER = BLOCK0;
    elseif TRAINING == 0 && BLOCK == 1
        TRIGGER = BLOCK1;
    elseif TRAINING == 0 && BLOCK == 2
        TRIGGER = BLOCK2;
    elseif TRAINING == 0 && BLOCK == 3
        TRIGGER = BLOCK3;
    elseif TRAINING == 0 && BLOCK == 4
        TRIGGER = BLOCK4;
    end

    if TRAINING == 1
        %     EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "START BLOCK"');
    else
        %     EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "START BLOCK"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end

    if TRAINING == 1
        disp('Start of Training Block.');
    else
        disp(['Start of Block ' num2str(BLOCK)]);
    end
    HideCursor(whichScreen);

    %% Experiment Loop
    noFixation = 0;

    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');


    [ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, equipment.greyVal);
    PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
    Screen('BlendFunction', ptbWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

    %  ptbWindow = PsychImaging('OpenWindow', screenID, equipment.greyVal);
    for trl = 1:experiment.nTrials
        %% Define thisGrating and thisCFI for usage as index
        thisGrating = trl/2;
        thisCFI = round(trl/2);


        %% Define stimulus or CFI trial
        moviename = '/home/methlab/Desktop/LEAMA/videoCFI.mp4';


        %% Start trial

        if mod(trl,2) == 0
            disp(['Start of Trial ' num2str(thisGrating)]); % Output of current trial iteration
        end


        %% Wait for release of possible keyboard presses
        KbReleaseWait;

        %% Prepare for trial execution

        % Select fixation cross or grid videos
        if mod(trl,2) == 1
            % CFI trial
            if crossSequence(thisCFI) == 0
                TRIGGER = FIXATION0;
            elseif crossSequence(thisCFI) == 1
                TRIGGER = FIXATION1;
            end
        elseif mod(trl,2) == 0
            % Stimulus trial
            if gratingSequence(thisGrating) == 1
                TRIGGER = PRESENTATION1;
            elseif gratingSequence(thisGrating) == 2
                TRIGGER = PRESENTATION2;
            elseif gratingSequence(thisGrating) == 3
                TRIGGER = PRESENTATION3;
            elseif gratingSequence(thisGrating) == 4
                TRIGGER = PRESENTATION4;
            elseif gratingSequence(thisGrating) == 5
                TRIGGER = PRESENTATION5;
            elseif gratingSequence(thisGrating) == 6
                TRIGGER = PRESENTATION6;
            elseif gratingSequence(thisGrating) == 7
                TRIGGER = PRESENTATION7;
            elseif gratingSequence(thisGrating) == 8
                TRIGGER = PRESENTATION8;
            end
        end

        if TRAINING == 1
            %         EThndl.sendMessage(TRIGGER);
            Eyelink('Message', num2str(TRIGGER));
            Eyelink('command', 'record_status_message "PRESENTATION"');
        else
            %         EThndl.sendMessage(TRIGGER);
            Eyelink('Message', num2str(TRIGGER));
            Eyelink('command', 'record_status_message "PRESENTATION"');
            sendtrigger(TRIGGER,port,SITE,stayup);
        end

        % Create keyboard monitoring queue with target buttons
        keyFlags = zeros(1,256); % An array of zeros
        keyFlags(spaceKeyCode) = 1; % Monitor only spaces
        % GetKeyboardIndices; % For checking keyboard number
        %  gki = 9; % has to be checked every time, that's why we created a dialog window now
        KbQueueCreate(gki, keyFlags); % Initialize the Queue


        % Set durations for CFIs
         if mod(trl,2) == 1
            % CFI trial
            timing.cfi(thisCFI) = (randsample(timing.cfilower:timing.cfiupper, 1))/1000; % Randomize the jittered central fixation interval on trial
%             maxTime = GetSecs + timing.cfi;
%         elseif mod(trl,2) == 0
%             % Stimulus trial
%             maxTime = GetSecs + timing.StimuliDuration; % Set maxTime to max. 2 seconds from start of video
        end

        % Start keyboard monitoring
        KbQueueStart(gki);
        queueStartTime = GetSecs;

        %% Display grating or fixation cross

        if mod(trl,2) == 1

            %% Present fixation cross

            start_time = GetSecs;

            while (GetSecs - start_time) < timing.cfi(thisCFI)

                if crossSequence(thisCFI) == 0 % No task condition
                    Screen('DrawLines', ptbWindow, fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor0,[screenCentreX screenCentreY],2);
                    Screen('Flip', ptbWindow)
                    WaitSecs(timing.cfi(thisCFI))
                elseif crossSequence(thisCFI) == 1 % Task condition
                    Screen('DrawLines', ptbWindow, fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor1,[screenCentreX screenCentreY],2);
                    Screen('Flip', ptbWindow)
                    WaitSecs(timing.cfi_task) % Show red cross for a short time
                    Screen('DrawLines', ptbWindow, fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor0,[screenCentreX screenCentreY],2);
                    Screen('Flip', ptbWindow)
                    WaitSecs(timing.cfi(thisCFI)-timing.cfi_task) % Show black cross for the rest of the CFI time
                end

            end


        else

            %% Present grating depending on sequence number

            % 4 different angles 
            if gratingSequence(thisGrating) == 1
                angle = 90;
            elseif gratingSequence(thisGrating) == 2
                angle = 0;      % 35 times
            elseif gratingSequence(thisGrating) == 3
                angle = 45;     % 35 times
            elseif gratingSequence(thisGrating) == 4
                angle = 115;    % 35 times
            elseif gratingSequence(thisGrating) == 5
                angle = 90;     % 35 times
            elseif gratingSequence(thisGrating) == 6
                angle = 0;      % 35 times
            elseif gratingSequence(thisGrating) == 7
                angle = 45;     % 35 times
            elseif gratingSequence(thisGrating) == 8
                angle = 115;    % 35 times
            end

            % Black and white grating
            color1 = [0 0 0 1];
            color2 = [1 1 1 1];
            baseColor = [0.5 0.5 0.5 1];

            % Defining contrast 
            if gratingSequence(thisGrating) == 1 || gratingSequence(thisGrating) == 2 || gratingSequence(thisGrating) == 3 || gratingSequence(thisGrating) == 4
                contrast = 1;
            elseif gratingSequence(thisGrating) == 5 || gratingSequence(thisGrating) == 6 || gratingSequence(thisGrating) == 7 || gratingSequence(thisGrating) == 8
                contrast = 0.5;
            end
contrast = 1;
            % default x + y size
            virtualSize = 400;

            % radius of the disc edge
            radius = floor(virtualSize / 2);

            % These settings are the parameters passed in directly to DrawTexture
%           angle % set above
            % phase
            phase = 0;
            % spatial frequency
            frequency = 0.05;
%           contrast % set above
            % sigma < 0 is a sinusoid.
            sigma = -1.0;

            % Setup defaults and unit color range:
            PsychDefaultSetup(2);

            % Query frame duration: We use it later on to time 'Flips' properly for an
            % animation with constant framerate:
            ifi = Screen('GetFlipInterval', ptbWindow);

            % Build a procedural texture, we also keep the shader as we will show how to
            % modify it (though not as efficient as using parameters in drawtexture)
            texture = CreateProceduralColorGrating(ptbWindow, virtualSize, virtualSize,...
                color1, color2, radius);

            % Preperatory flip
            % showTime = 3;
            vbl = Screen('Flip', ptbWindow);
            % tstart = vbl + ifi; %start is on the next frame
            start_time = GetSecs;


                % Draw the shader texture with parameters
                Screen('DrawTexture', ptbWindow, texture, [], [],...
                    angle, [], [], baseColor, [], [],...
                    [phase, frequency, contrast, sigma]);

                vbl = Screen('Flip', ptbWindow, vbl * ifi);
%               phase = phase - 15; This would be for moving the grating
            
screenshot = Screen('GetImage', ptbWindow)
 
        end % Presentation of grating or cfi


        % Stop keyboard monitoring
        KbQueueStop(gki);

        % Get infos for key presses during KbQueue (only for Gratings)
        [pressed, firstPress, ~, ~, ~] = KbQueueCheck(gki);

        if mod(trl,2) == 1
            % Get and save reaction time for each trial
            if firstPress(spaceKeyCode) > 0
                reactionTime(thisCFI) = firstPress(spaceKeyCode) - queueStartTime;
            elseif firstPress(spaceKeyCode) == 0
                reactionTime(thisCFI) = NaN;
            end
        end

        % Remove all unprocessed events from the queue and zeros out any already scored events
        KbQueueFlush(gki);

        % Release queue-associated resources
        KbQueueRelease(gki);


        % Save responses as triggers (only for Gratings, always sent out after video is stopped)
        if mod(trl,2) == 1
            % Save responses (RESPONSE TIMESTAMPS ARE WRONG --> use reactionTime) and send triggers
            if pressed > 0
                data.allResponses(thisCFI) = spaceKeyCode;
                TRIGGER = RESP;
            elseif pressed == 0
                data.allResponses(thisCFI) = 0;
                TRIGGER = NO_RESP;
            end
        end

        %     if TRAINING == 1
        %         %     EThndl.sendMessage(TRIGGER);
        %         Eyelink('Message', num2str(TRIGGER));
        %         Eyelink('command', 'record_status_message "END BLOCK"');
        %         disp('End of Training Block.');
        %     else
        %         %     EThndl.sendMessage(TRIGGER);
        %         Eyelink('Message', num2str(TRIGGER));
        %         Eyelink('command', 'record_status_message "END BLOCK"');
        %         sendtrigger(TRIGGER,port,SITE,stayup);
        %         disp(['End of Block ' num2str(BLOCK)]);
        %     end
        %     if mod(trl,2) == 0
        %         % Stop playback:
        %         Screen('PlayMovie', movie, 0);
        %         % Close movie:
        %         Screen('CloseMovie', movie);

        %     end



        if mod(trl,2) == 1
            % Check if response was correct (only for Gratings)
            if crossSequence(thisCFI) == 0 && data.allResponses(thisCFI) == 0 % Skewed + NO button press = correct answer
                data.allCorrect(thisCFI) = 1;
            elseif crossSequence(thisCFI) == 1 && data.allResponses(thisCFI) == spaceKeyCode % Vert + button press
                data.allCorrect(thisCFI) = 1;
            else
                data.allCorrect(thisCFI) = 0; % Wrong response
            end

            % Display (in-)correct response in CW (only for Gratings)
            if data.allCorrect(thisCFI) == 1
                feedbackText = 'Correct!';
            elseif data.allCorrect(thisCFI) == 0
                feedbackText = 'Incorrect!';
            end
            disp(['Response to Trial ' num2str(thisCFI) ' is ' feedbackText]);
        end



        %% Fixation Check

        %     % Check if subject fixate at center, give warning if not
        %     checkFixation;
        %     if noFixation > 2
        %         disp('Insufficient fixation!')
        %         noFixation = 0; % reset
        %     end

        %     Screen('DrawLines',ptbWindow,fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor,[screenCentreX screenCentreY],2); % Draw fixation cross
        %     Screen('Flip', ptbWindow, [], 1);



        %% End for loop over all trials
    end % trial loop

    % Wait a few seconds before continuing in order to record the last stimuli presentation
    WaitSecs(4)


    %% End task and save data

    % Send triggers to end task
    endT = Screen('Flip',ptbWindow);
%     if TRAINING == 1
%         %     EThndl.sendMessage(TASK_END,endT);
%         Eyelink('Message', num2str(TASK_END));
%         Eyelink('command', 'record_status_message "TASK_END"');
%     else
%         %     EThndl.sendMessage(TASK_END,endT);
%         Eyelink('Message', num2str(TASK_END));
%         Eyelink('command', 'record_status_message "TASK_END"');
%         sendtrigger(TASK_END,port,SITE,stayup)
%     end

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
        %     EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "END BLOCK"');
        disp('End of Training Block.');
    else
        %     EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "END BLOCK"');
        sendtrigger(TRIGGER,port,SITE,stayup);
        disp(['End of Block ' num2str(BLOCK)]);
    end

    % Save data
    subjectID = num2str(subject.ID);
    filePath = fullfile(DATA_PATH, subjectID);
    
    if exist(filePath,'dir') == 0 % if matlab broke down and this is a retry, path could already exist
        mkdir(filePath)
    end

    if TRAINING == 1
        fileName = [subjectID, '_training.mat'];
    elseif TRAINING == 0
        fileName = [subjectID '_', TASK, '_block' num2str(BLOCK), '.mat'];
    end

    % % Compute accuracy and report after each block (no additional cash for training task)
    % if TRAINING == 1
    %     % Get sum of correct responses, but ignore first and last data point
    %     totalCorrect = sum(data.allCorrect(1, 2:end-1));
    %     totalTrials = trl-2;
    %     percentTotalCorrect = totalCorrect / totalTrials * 100;
    %     format bank % Change format for display
    %     feedbackBlockText = ['Your accuracy in the training task was ' num2str(percentTotalCorrect) ' %. '];
    %     disp(['Participant ' subjectID ' had an accuracy of ' num2str(percentTotalCorrect) ' % in the training task.'])
    %     DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.textVal);
    %     format default % Change format back to default
    %     Screen('Flip',ptbWindow);
    %     WaitSecs(5);
    % elseif BLOCK == 1
    %     % Get sum of correct responses, but ignore first and last data point
    %     totalCorrect = sum(data.allCorrect(1, 2:end-1));
    %     totalTrials = trl-2;
    %     percentTotalCorrect(BLOCK) = totalCorrect / totalTrials * 100;
    %     format bank % Change format for display
    %     if percentTotalCorrect(BLOCK) > 80
    %         amountCHFextra(BLOCK) = percentTotalCorrect(BLOCK)*0.02;
    %         feedbackBlockText = ['Your accuracy was ' num2str(percentTotalCorrect(BLOCK)) ' %. ' ...
    %             '\n\n Because of your accuracy you have been awarded an additional ' num2str(amountCHFextra(BLOCK)) ' CHF.' ...
    %             '\n\n Keep it up!'];
    %     elseif percentTotalCorrect(BLOCK) < 80 && BLOCK == 1
    %         amountCHFextra(BLOCK) = 0;
    %         feedbackBlockText = ['Your accuracy was ' num2str(percentTotalCorrect(BLOCK)) ' %. ' ...
    %             '\n\n Your accuracy was very low in this block. Please stay focused!'];
    %         disp(['Low accuracy in Block ' num2str(BLOCK) '.']);
    %     end
    %     DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.textVal);
    %     disp(['Participant ' subjectID ' was awarded CHF ' num2str(amountCHFextra(BLOCK)) ' for an accuracy of ' num2str(percentTotalCorrect(BLOCK)) ' % in Block ' num2str(BLOCK) '.'])
    %     format default % Change format back to default
    %     Screen('Flip',ptbWindow);
    %     WaitSecs(5);
    % end

    % Save data
    saves = struct;
    saves.data = data;
    saves.data.spaceKeyCode = spaceKeyCode;
    saves.data.reactionTime = reactionTime;
    saves.experiment = experiment;
    saves.screenWidth = screenWidth;
    saves.screenHeight = screenHeight;
    saves.screenCentreX = screenCentreX;
    saves.screenCentreY = screenCentreY;
    saves.startExperimentTime = startExperimentTime;
    saves.startExperimentText = startExperimentText;
    saves.stimulus = stimulus;
    saves.subjectID = subjectID;
    saves.subject = subject;
    saves.timing = timing;
    saves.waitResponse = waitResponse;
    saves.flipInterval = flipInterval;

    % Save triggers
    trigger = struct;
    trigger.RESTING_START = 10;
    trigger.BLOCK1 = BLOCK1;
    trigger.BLOCK2 = BLOCK2;
    trigger.BLOCK3 = BLOCK3;
    trigger.BLOCK4 = BLOCK4;
    trigger.BLOCK0 = BLOCK0;

    trigger.FIXATION0 = FIXATION0;
    trigger.FIXATION1 = FIXATION1;
    trigger.PRESENTATION1 = PRESENTATION1;
    trigger.PRESENTATION2 = PRESENTATION2;
    trigger.PRESENTATION3 = PRESENTATION3;
    trigger.PRESENTATION4 = PRESENTATION4;
    trigger.PRESENTATION5 = PRESENTATION5;
    trigger.PRESENTATION6 = PRESENTATION6;
    trigger.PRESENTATION7 = PRESENTATION7;
    trigger.PRESENTATION8 = PRESENTATION8;

    trigger.RESP_YES = RESP;
    trigger.RESP_NO = NO_RESP;

    trigger.RESTING_END = 90;
    trigger.BLOCK1_END = BLOCK1_END;
    trigger.BLOCK2_END = BLOCK2_END;
    trigger.BLOCK3_END = BLOCK3_END;
    trigger.BLOCK4_END = BLOCK4_END;
    trigger.BLOCK0_END = BLOCK0_END;
    

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
      %  if percentTotalCorrect >= THRESH
            breakInstructionText = 'Well done! \n\n Press any key to finalize the training block.';
     %   else
     %      breakInstructionText = ['Score too low! ' num2str(percentTotalCorrect) ' % correct. ' ...
     %           '\n\n Press any key to repeat the training task.'];
     %   end
    elseif TRAINING == 0 && BLOCK == 4
        breakInstructionText = ['End of the Task! ' ...
            '\n\n Thank you very much for your participation.'...
            '\n\n Please press any key to finalize the experiment.'];
    else
        breakInstructionText = ['Break! Rest for a while... ' ...
            '\n\n Press any key to start the mandatory break of at least 30 seconds.'];
    end

    DrawFormattedText(ptbWindow,breakInstructionText,'center','center',color.textVal);
    Screen('Flip',ptbWindow);
    waitResponse = 1;
    while waitResponse
        [time, keyCode] = KbWait(-1,2);
        waitResponse = 0;
    end

    % Break for loop if it was the training block
    if TRAINING == 1 && BLOCK == 1
        break
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


    % Wait at least 30 Seconds between Blocks (only after Block 1 has finished, not after Block 2)
    % if TRAINING == 1 && percentTotalCorrect < THRESH
    %     waitTime = 30;
    %     intervalTime = 1;
    %     timePassed = 0;
    %     printTime = 30;
    %
    %     waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
    %         ' \n\n ' ...
    %         ' \n\n You can repeat the training task afterwards.'];
    %
    %     DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
    %     Screen('Flip',ptbWindow);
    %
    %     while timePassed < waitTime
    %         pause(intervalTime);
    %         timePassed = timePassed + intervalTime;
    %         printTime = waitTime - timePassed;
    %         waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
    %             ' \n\n ' ...
    %             ' \n\n You can repeat the training task afterwards.'];
    %         DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
    %         Screen('Flip',ptbWindow);
    %     end
    % elseif BLOCK == 1 && TRAINING == 1
    %     waitTime = 30;
    %     intervalTime = 1;
    %     timePassed = 0;
    %     printTime = 30;
    %
    %     waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
    %         ' \n\n ' ...
    %         ' \n\n The Grating task will start afterwards.'];
    %
    %     DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
    %     Screen('Flip',ptbWindow);
    %
    %     while timePassed < waitTime
    %         pause(intervalTime);
    %         timePassed = timePassed + intervalTime;
    %         printTime = waitTime - timePassed;
    %         waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
    %             ' \n\n ' ...
    %             ' \n\n The Grating task will start afterwards.'];
    %         DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
    %         Screen('Flip',ptbWindow);
    %     end
    % end

end 

%% Close Psychtoolbox window
Screen('CloseAll');