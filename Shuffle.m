%% CREATE GRATING SEQUENCES

% Define proportions 
N_type = 50; % 70%

% create list of N_look * 8 and N_task * 8
gratingselection = [
    zeros(N_type + 50,1) + 1; % + 25 because need double amount of horizontal gratings
    zeros(N_type,1) + 2;
    zeros(N_type,1) + 3;
    zeros(N_type,1) + 4;
    zeros(N_type + 50,1) + 5; % + 25 because need double amount of horizontal gratings
    zeros(N_type,1) + 6;
    zeros(N_type,1) + 7;
    zeros(N_type,1) + 8];

% Shuffle the initial list and save it
rng(10) % 'set seed' so randperm always shuffles in the same way
gratingSelection_shuffled = gratingselection(randperm(length(gratingselection)));
gratingSeq = transpose(gratingSelection_shuffled);


% create four lists out of the one
gratingSequence1 = gratingSeq(:,1:125);
gratingSequence2 = gratingSeq(:,126:250);
gratingSequence3 = gratingSeq(:,251:375);
gratingSequence4 = gratingSeq(:,376:500);

%% CREATE CROSS SEQUENCES

% Define proportions 
N_look = 100; % 80%
N_cross = 25; % 20%

% create list of N_look * 8 and N_task * 8
crossSelection = [
    zeros(N_look,1);
    zeros(N_cross,1) + 1;];

% Shuffle the initial list and save it
rng(10) % 'set seed' so randperm always shuffles in the same way
crossSequence1 = transpose(crossSelection(randperm(length(crossSelection))));
rng(11) % 'set seed' so randperm always shuffles in the same way
crossSequence2 = transpose(crossSelection(randperm(length(crossSelection))));
rng(12) % 'set seed' so randperm always shuffles in the same way
crossSequence3 = transpose(crossSelection(randperm(length(crossSelection))));
rng(13) % 'set seed' so randperm always shuffles in the same way
crossSequence4 = transpose(crossSelection(randperm(length(crossSelection))));


%% Create Training Sequences

% Define proportions 
N_type = 1; % 70%

% create list of N_look * 8 and N_task * 8
gratingselection = [
    zeros(N_type,1) + 1; % + 25 because need double amount of horizontal gratings
    zeros(N_type,1) + 2;
    zeros(N_type,1) + 3;
    zeros(N_type,1) + 4;
    zeros(N_type,1) + 5; % + 25 because need double amount of horizontal gratings
    zeros(N_type,1) + 6;
    zeros(N_type,1) + 7;
    zeros(N_type,1) + 8];

% Shuffle the initial list and save it
rng(10) % 'set seed' so randperm always shuffles in the same way
gratingSelection_shuffled = gratingselection(randperm(length(gratingselection)));
gratingSeqTraining = transpose(gratingSelection_shuffled);

gratingSequence0 = gratingSeqTraining;


% Define proportions 
N_look = 7; % 80%
N_cross = 1; % 20%

% create list of N_look * 8 and N_task * 8
crossSelection = [
    zeros(N_look,1);
    zeros(N_cross,1) + 1;];

% Shuffle the initial list and save it
rng(10) % 'set seed' so randperm always shuffles in the same way
crossSequenceTraining = transpose(crossSelection(randperm(length(crossSelection))));

crossSequence0 = crossSequenceTraining;









