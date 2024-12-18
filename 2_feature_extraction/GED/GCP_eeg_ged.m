%% GCP General Eigendecomposition for Gamma Band (30 - 90 Hz) Activations

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

%% General Eigendecomposition
for subj = 10%:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    
    % Initialise a structure to hold all blocks
    EEG_all = [];
    
    % Read blocks
    for block = 1:1
        try % Do not load empty blocks
            load(strcat(subjects{subj}, '_EEG_ET_GCP_block', num2str(block), '_merged.mat'))
            
            EEG.data = double(EEG.data);
            
            % Identify the indices of the EEG channels
            chanlocs = {EEG.chanlocs.labels}; % Channel labels
            exclude_channels = {'B*', 'HEOGR', 'HEOGL', 'VEOGU', 'VEOGL', ...
                                'L-GAZE-X', 'L-GAZE-Y', 'L-AREA', ...
                                'R-GAZE-X', 'R-GAZE-Y', 'R-AREA'};
            include_idx = find(~cellfun(@(x) any(cellfun(@(y) ~isempty(regexp(x, y, 'once')), exclude_channels)), chanlocs));
            
            % Select only EEG channels using indices
            EEG = pop_select(EEG, 'channel', include_idx);

            % Rereference to average
            EEG = pop_reref(EEG, []);
            

            % Filter EEG data
            EEG = pop_eegfiltnew(EEG, 0.2, 30);
            
            % Append data
            if isempty(EEG_all)
                EEG_all = EEG; % Initialise with the first block
            else
                EEG_all = pop_mergeset(EEG_all, EEG); % Append subsequent blocks
            end
            clear EEG
            fprintf('Subject GCP %.3s (%.3d/%.3d): Block %.1d loaded and appended\n', subjects{subj}, subj, length(subjects), block)
        catch ME
            fprintf('Block %.1d for Subject GCP %.3s could not be loaded: %s\n', block, subjects{subj}, ME.message)
        end
    end

%% Find timepoints for LC and HC stimulus presentation
stim_idx = ismember({EEG_all.event.type}, {'61', '62'});
latencies = cell2mat({EEG_all.event.latency});
stim_latency = latencies(stim_idx);

% Define how many timepoints (at 500Hz) should be uses for pre and post segments
pre_timepoints = 250; % -1000 ms until 0 (stimuls presentation)
post_timpoints = 250; % 0 (stimulus presentation) until +2000ms 

%% Pre-stimulus covariance
addpath('/Users/Arne/Documents/matlabtools/colormaps/')
% Get colormap from cbrewer
cmap = cbrewer('div', 'RdBu', 256);
% Ensure values are in the range [0, 1]
cmap = cmap - min(cmap(:));  % Normalize the colormap
cmap = cmap / max(cmap(:));  % Scale to [0, 1]
cmap = flipud(cmap);  % Flip the colormap upside down
cmapRdBu = cmap;

covPre = 0;
tmpdat = [];
for ti=1:length(stim_latency)
    tmpdat = EEG_all.data(:, stim_latency(ti)-pre_timepoints : stim_latency(ti) + post_timpoints); %indices if pre
    tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
    nwin = length(stim_latency(ti)-pre_timepoints : stim_latency(ti) + post_timpoints);
    covPre   = covPre + (tmpdat*tmpdat')/nwin; % nwin = length in TIME POINTS von pre für jedes trial
end
covPre = covPre./ti;

% Sanity check
close all
figure;
set(gcf, 'Position', [0, 0, 800, 800], 'Color', 'w'); % Set figure size and background
maxabs = max(abs(covPre(:))); % Maximum absolute value in the covariance matrix
imagesc(covPre, [-maxabs*0.5, maxabs*0.5]); 
colormap(cmapRdBu);
colorbar;
title('Pre-stimulus Covariance', 'FontSize', 20, 'FontWeight', 'bold');
xlabel('Channels', 'FontSize', 16);
ylabel('Channels', 'FontSize', 16);
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'Box', 'on');

%% Post-stimulus covariance
covPost = 0;
tmpdat = [];
for ti=1:length(stim_latency)
    tmpdat = EEG_all.data(:, stim_latency(ti) : stim_latency(ti)+post_timpoints); %indices if pre
    tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
    nwin = length(stim_latency(ti) : stim_latency(ti)+post_timpoints);
    covPost   = covPost + (tmpdat*tmpdat')/nwin; % nwin = length in TIME POINTS von pre für jedes trial
end
covPost = covPost./ti;

% Sanity check
close all
figure;
set(gcf, 'Position', [0, 0, 800, 800], 'Color', 'w'); % Set figure size and background
maxabs = max(abs(covPost(:))); % Maximum absolute value in the covariance matrix
imagesc(covPost, [-maxabs*0.5, maxabs*0.5]); 
colormap(cmapRdBu);
colorbar;
title('Pre-stimulus Covariance', 'FontSize', 20, 'FontWeight', 'bold');
xlabel('Channels', 'FontSize', 16);
ylabel('Channels', 'FontSize', 16);
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'Box', 'on');

%% GED
% Conduct general eigendecomposition
[eigenvecs,eigenvals] = eig(covPost,covPre);

% Find best component 
[~,maxcomp] = sort(diag(eigenvals));

% Compute filter projection by multiplying filtered data with best component
% maxcopm(end) = eigenvector with highest eigenvalue
topo_map_prepost = covPost * eigenvecs(:, maxcomp(end));

% Gamma time series component (raw data * eigenvector)
gammacomp_pre = EEG_all.data' * eigenvecs(:, maxcomp(end));


% Pre-stimulus gamma component


% Post-stimulus gamma component
figure;
% plot(EEG_all.times(1:1000), EEG_all.data(, 1:1000))

figure;
imagesc(EEG_all.data, [-5, 5])


%% Plot hilbert-transformed power timeseries
close all
figure
gcomp = gammacomp_pre;
hilb = hilbert(gcomp')';
pow = abs(hilb).^2;
plot(pow(1:100))
end

%%
figure
topoplot(topo_map_prepost, EEG_all.chanlocs)
% 
% figure;
% plot(gammacomp(1:1000))
% 
% 
% hilb = hilbert(gammacomp')';
% pwr = abs(hilb).^2;
% 
% figure;
% plot(pwr(1:100))
%pwr ist zeitreihe
%power in segmente prä vs. post als zeitreihe plotten


%powerscptrm von gammacomp


% für jedes subj plot


