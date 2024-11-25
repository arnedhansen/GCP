try
    % Setup for Calibration
    % TODO ADD MORE INFO FOR SETTINGS!
    try
        Snd('Close');
        PsychPortAudio('Close');
    catch
    end
    fprintf('Starting Calibration \n' );
    Eyelink('Command', 'saccade_velocity_threshold = 35');
    Eyelink('Command', 'saccade_acceleration_threshold = 9500');
    Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA');
    Eyelink('Command', 'active_eye = BOTH');
    Eyelink('Command', 'calibration_type = HV9');
    Eyelink('Command', 'enable_automatic_calibration = YES');
    Eyelink('Command', 'automatic_calibration_pacing = 500');
    Eyelink('Command', 'set_idle_mode');

    % Run, Calibration, and Validation need to be accessed over the Stim computer
    HideCursor(whichScreen);
    Eyelink('StartSetup',1);
    disp('CALIBRATION DONE');
    try
        Snd('Close');
        PsychPortAudio('Close');
    catch
    end
catch
    disp('ERROR running the calibration');
end