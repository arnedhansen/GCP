%% Setup
startup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.2');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/merged/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%%
for subj = 1%:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)

    % Read blocks
    EEG_allblocks{subj} = [];
    for block = 1:4
        try % Do not load empty blocks
            load(strcat(subjects{subj}, '_EEG_ET_GCP_block', num2str(block), '_merged.mat'))

            % Identify the indices of the EEG channels
            chanlocs = {EEG.chanlocs.labels}; % Channel labels
            exclude_channels = {'B*', 'HEOGR', 'HEOGL', 'VEOGU', 'VEOGL', ...
                'L-GAZE-X', 'L-GAZE-Y', 'L-AREA', ...
                'R-GAZE-X', 'R-GAZE-Y', 'R-AREA'};
            include_idx = find(~cellfun(@(x) any(cellfun(@(y) ~isempty(regexp(x, y, 'once')), exclude_channels)), chanlocs));

            % Select only EEG channels using indices
            EEG = pop_select(EEG, 'channel', include_idx);

            % Filter EEG data
            EEG = pop_eegfiltnew(EEG, 30, 90);
            if isempty(EEG_allblocks{subj})
                EEG_allblocks{subj} = EEG; % Initialise with the first block
            else
                EEG_allblocks{subj} = pop_mergeset(EEG_allblocks{subj}, EEG); % Append subsequent blocks
            end
            clear EEG
            fprintf('Subject GCP %.3s (%.3d/%.3d): Block %.1d loaded and appended\n', subjects{subj}, subj, length(subjects), block)
        catch ME
            fprintf('Block %.1d for Subject GCP %.3s could not be loaded: %s\n', block, subjects{subj}, ME.message)
        end
    end
end

%%
EEG = alleeg{1};
stim_idx = ismember({EEG.event.type}, {'61', '62'});
latencies = cell2mat({EEG.event.latency});
stim_latency = latencies(stim_idx);

pre_timepoints = 500
post_timpoints = 1000


%pre covariance
covPre = 0;
for ti=1:length(stim_latency)
    tmpdat = EEG.data(:, stim_latency(ti)-pre_timepoints : stim_latency(ti)); %indices if pre
    tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
    nwin = length(stim_latency(ti)-pre_timepoints : stim_latency(ti));
    covPre   = covPre + (tmpdat*tmpdat')/nwin; % nwin = length in TIME POINTS von pre für jedes trial
end
covPre = covPre./ti;

figure;
imagesc(covPre, [-5, 5])

%post covariance
covPost = 0;
for ti=1:length(stim_latency)
    tmpdat = EEG.data(:, stim_latency(ti) : stim_latency(ti)+post_timpoints); %indices if pre
    tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
    nwin = length(stim_latency(ti) : stim_latency(ti)+post_timpoints);
    covPost   = covPost + (tmpdat*tmpdat')/nwin; % nwin = length in TIME POINTS von pre für jedes trial
end
covPost = covPost./ti;

figure;
imagesc(covPost, [-5, 5])


%% GED

% GED
[evecsT,evals] = eig(covPost,covPre);

% find best component and compute filter projection
[~,maxcomp] = sort(diag(evals));
post_map_orig    = covPost*evecsT(:,maxcomp(end));

% theta time series component
%%% nicht EEG.data, sondern einfach gamma-gefilterte daten
gammacomp = EEG.data' * evecsT(:,maxcomp(end)); %%%% VERY IMPORTANT: select approprioate component (end - XXX)
%maxcopm(end) = eigenvector with highest eigenvalue

figure
topoplot(gammacomp, EEG.chanlocs)

figure;
plot(gammacomp(1:1000))


hilb = hilbert(gammacomp')';
pwr = abs(hilb).^2;

figure;
plot(pwr(1:100))
%pwr ist zeitreihe
%power in segmente prä vs. post als zeitreihe plotten


%powerscptrm von gammacomp


% für jedes subj plot


