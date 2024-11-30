%% Screenshot function for Psychtoolbox (PTB) presentations
function screenshot(screenshotFilename, enableScreenshots)
    if enableScreenshots == 1
        imageArray = Screen('GetImage', ptbWindow);
        imwrite(imageArray, screenshotFilename);
    end
end