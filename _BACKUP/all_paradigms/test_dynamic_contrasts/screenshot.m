%% Screenshot function for Psychtoolbox (PTB) presentations
function screenshot(screenshotFilename, ptbWindow, enableScreenshots)
    if enableScreenshots == 1
        imageArray = Screen('GetImage', ptbWindow);
        imwrite(imageArray, screenshotFilename);
    end
end