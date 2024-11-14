%% Reminder to turn power on again

% Create dialog box for GKI

defAns = {''};

while true
    prompt = {'\bf \fontsize{19} \color[rgb]{1,0,0} HAVE YOU TURNED THE POWER STRIP BACK ON? (if so, press Ctrl+C and type "yes")'};
    options.Interpreter = 'tex';
    dims = [1 120];
    box = inputdlg(prompt, 'POWER', dims, defAns,options);
    answer = char(box(1));
    if length(answer) == 3 % Ensure that yes response is made
        break
    end
end
    