%% dialogGKI

% Create dialog box for GKI (Keyboard input needs to be checked every time)
defAns = {''};
while true
    prompt = {'GKI'};
    box = inputdlg(prompt, 'Enter Subject Information', 1,defAns);
    gki = char(box(1));
    if length(gki) == 1 || length(gki) == 2 % Ensure response made in subject ID
        break
    end
end

% Convert subject data to numeric
gki = str2num(gki);
