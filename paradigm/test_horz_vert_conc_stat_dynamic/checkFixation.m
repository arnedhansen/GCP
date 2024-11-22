function noFixation = checkFixation(screenWidth, screenHeight, screenCentreX, screenCentreY)
    numSamples = 125; % 500ms 250; % 1 second 
    fixThresh = numSamples * 0.75; % 75% threshold
    distOK = 45 + 45; % 1 degree from the center (+ 0.5 deg of ET error)

    % Initialize noFixation counter
    noFixation = 0;

    % Collect gaze data for a specified number of samples
    samples = zeros(numSamples, 2); % Initialize matrix for gaze samples

    for i = 1:numSamples
        % Fetch gaze data sample
        evt = Eyelink('NewestFloatSample');
        gaze_x = evt.gx(1);
        gaze_y = evt.gy(1);
        samples(i, :) = [gaze_x, gaze_y]; % Store gaze sample
        WaitSecs(0.002) % 2 ms for 125 Hz            (0.004); % Wait for 4 ms to get approximately 250 Hz sampling rate
    end

    % Check fixation
    validSamples = sum(samples(:, 1) > 0 & samples(:, 2) > 0); % Count valid samples
    if validSamples > fixThresh % Ensure enough valid samples before checking fixation
        xFix = sum(samples(:, 1) > screenCentreX - distOK & samples(:, 1) < screenCentreX + distOK);
        yFix = sum(samples(:, 2) > screenCentreY - distOK & samples(:, 2) < screenCentreY + distOK);

        if xFix <= fixThresh || yFix <= fixThresh
            % No fixation detected
            noFixation = 1;
        end
    else
        disp('Not enough valid samples for fixation check.');
    end
end