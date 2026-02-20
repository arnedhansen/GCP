function noFixation = checkFixation(screenCentreX, screenCentreY, fixCheckDuration)
    % checkFixation  Online pre-stimulus fixation check via EyeLink.
    %
    %   noFixation = checkFixation(screenCentreX, screenCentreY, fixCheckDuration)
    %
    %   Polls EyeLink gaze samples for fixCheckDuration seconds and evaluates
    %   whether the participant maintained fixation within a window around screen
    %   center. Call this while a fixation cross is already drawn on screen —
    %   the function does not perform any Screen operations.
    %
    %   Returns:
    %       noFixation = 0  — participant was fixating (pass)
    %       noFixation = 1  — participant was NOT fixating (fail)

    sampleInterval = 0.004; % 4 ms between samples (~250 Hz effective polling)
    numSamples = ceil(fixCheckDuration / sampleInterval);
    fixThresh = 0.80; % 80% of valid samples must be within fixation window
    distOK = 90; % Fixation window radius in pixels (~1.8 dva at 50 ppd)

    noFixation = 0;
    samples = NaN(numSamples, 2);
    nCollected = 0;

    startTime = GetSecs;
    for i = 1:numSamples
        if (GetSecs - startTime) >= fixCheckDuration
            break;
        end

        evt = Eyelink('NewestFloatSample');
        gx = evt.gx(1);
        gy = evt.gy(1);

        % EyeLink returns large negative values (e.g. -32768) for missing data
        if gx >= 0 && gy >= 0
            nCollected = nCollected + 1;
            samples(nCollected, :) = [gx, gy];
        end

        WaitSecs(sampleInterval);
    end

    minValidSamples = round(numSamples * 0.50);
    if nCollected < minValidSamples
        noFixation = 1;
        return;
    end

    validSamples = samples(1:nCollected, :);
    withinWindow = abs(validSamples(:,1) - screenCentreX) < distOK & ...
                   abs(validSamples(:,2) - screenCentreY) < distOK;
    fixRatio = sum(withinWindow) / nCollected;

    if fixRatio < fixThresh
        noFixation = 1;
    end
end
