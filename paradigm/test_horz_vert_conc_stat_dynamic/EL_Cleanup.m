%% Cleanup Routine for The Eylink Eyetracker
IOPort('CloseAll');
Eyelink('Command', 'clear_screen 0');

% Shutdown Eyelink:
Eyelink('StopRecording');
Eyelink('CloseFile');
Eyelink('Shutdown');
fprintf('Stopped the Eyetracker\n');