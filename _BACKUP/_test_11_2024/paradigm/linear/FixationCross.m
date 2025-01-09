%% Running

% Set up equipment parameters
equipment.viewDist = 800;               % Viewing distance in millimetres
equipment.ppm = 3.6;                    % Pixels per millimetre !! NEEDS TO BE SET. USE THE MeasureDpi FUNCTION !!
equipment.greyVal = .47;
equipment.blackVal = 0;
equipment.whiteVal = 1;
equipment.gammaVals = [1 1 1];

% Set up stimulus parameters Fixation
stimulus.fixationOn = 1;                % Toggle fixation on (1) or off (0)
stimulus.fixationSize_dva = 0.3;          % Size of fixation cross in degress of visual angle (smaller version: 1.5)
stimulus.fixationColor = [0 0 0];       % Color of fixation cross (1 = white, 0 = black, [1 0 0] = red, [1 0.05 0.5] = pink)
stimulus.fixationLineWidth = 1.3;         % Line width of fixation cross (thicker version: 1.5)

% Location
stimulus.regionHeight_dva = 7.3;         % Height of the region
stimulus.regionWidth_dva = 4;            % Width of the region
stimulus.regionEccentricity_dva = 3;     % Eccentricity of regions from central fixation

try
% Imaging set up
screenNumber=max(Screen('Screens'));
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
%     Screen('Preference', 'SkipSyncTests', 0); % For linux

% Window set-up
[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, equipment.greyVal);
PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
[screenWidth, screenHeight] = RectSize(winRect);
screenCentreX = round(screenWidth/2);
screenCentreY = round(screenHeight/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);
Screen('BlendFunction', ptbWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
experiment.runPriority = MaxPriority(ptbWindow);

% Calculate equipment parameters
equipment.mpd = (equipment.viewDist/2)*tan(deg2rad(2*stimulus.regionEccentricity_dva))/stimulus.regionEccentricity_dva; % Millimetres per degree
equipment.ppd = equipment.ppm*equipment.mpd;  

% Fix coordiantes for fixation cross
stimulus.fixationSize_pix = round(stimulus.fixationSize_dva*equipment.ppd); 
fixHorizontal = [round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2) 0 0];
fixVertical = [0 0 round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2)];
fixCoords = [fixHorizontal; fixVertical];


% Trial run
Screen('DrawLines', ptbWindow, fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor,[screenCentreX screenCentreY],2);
Screen('Flip', ptbWindow)
WaitSecs(timing.cfi)
% sca;
catch
%      sca;

end

