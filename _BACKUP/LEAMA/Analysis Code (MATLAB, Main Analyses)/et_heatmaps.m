%% Eye tracker: Analysis %%

%% Description %%
% This script analyses the eye tracker data of the ten subjects descriptively and statistically.
% It includes the performance of cluster-based permutation t-tests on the heatmap difference between the contrasts
% (and the task conditions vs. baseline for the full eye tracker data).


%% Script Structure
% 0. Preparation
% 1. Prepare the data
% 2. Load the data
% 3. Average over trials (and subjects)
% 4. First Level Analysis (descriptive)
% 5. Second Level Analysis (descriptive and statistical)



%% 0. Preparation

clear all
close all

%% 0.1 Path and subjects
subjects = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '11'};
path = '/Volumes/methlab/Students/Lea Baechlin/data/';

%% 0.2 Add clean version of Fieldtrip
cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin')
run startup

%% 0.3 Settings

% Smoothing: 'Granular', 'Coarse'
smoothing = 'Granular';

% Variable to investigate: 'none' (whole ET data), 'L_fixation', 'L_saccade'
variable = 'none';

% Threshold for blink exclusion
BLINK_THRESHOLD = 0.1;

% Definition of kernel standard deviation
if strcmp(smoothing, 'Granular')
    sig = 10;
elseif strcmp(smoothing, 'Coarse')
    sig = 80;
end

%% 0.4 Set parameteres for circle (to add to figures)

set(gcf, 'Renderer', 'opengl');
center = [400, 300];    % Screen center
radius = 200;           % Same radius as stimuli had
circlecolor1 = [0.9 0.9 0.9];
circlecolor2 = [0.7 0.7 0.7];

% Number of points to draw the circle
numPoints = 10000; % for a smooth circle

% Parametric equations for a circle
theta = linspace(0, 2*pi, numPoints);
x = center(1) + radius * cos(theta);
y = center(2) + radius * sin(theta);


%% 0.5 Prepare the data

% The three main hypotheses of this study are summarized and schematically visualized in Figure 3.
for subj =1:length( subjects)

    datapath = strcat(path, subjects{subj});
    cd(datapath)

    % Load the data
    load dataet


    %% 0.5.1 Identify conditions

    % High contrast
    indhc_hor = find(dataet.trialinfo == 0);
    indhc_ver = find(dataet.trialinfo == 90);
    indhc_45 = find(dataet.trialinfo == 45);
    indhc_115 = find(dataet.trialinfo == 115);

    % Low contrast
    indlc_hor = find(dataet.trialinfo == 1);
    indlc_ver = find(dataet.trialinfo == 91);
    indlc_45 = find(dataet.trialinfo == 46);
    indlc_115 = find(dataet.trialinfo == 116);

    %% 0.5.1 Split into high and low contrast and baseline

    cfg             = [];
    cfg.latency     = [0 2];
    cfg.trials      = [indlc_hor; indlc_ver; indlc_45; indlc_115];
    dataetlow = ft_selectdata(cfg, dataet);
    cfg.trials      = [indhc_hor; indhc_ver; indhc_45; indhc_115];
    dataethigh = ft_selectdata(cfg, dataet);

    cfg             = [];
    cfg.latency     = [-2 0];
    cfg.trials      = [indlc_hor; indlc_ver; indlc_45; indlc_115];
    dataetlowbl = ft_selectdata(cfg, dataet);
    cfg.trials      = [indhc_hor; indhc_ver; indhc_45; indhc_115];
    dataethighbl = ft_selectdata(cfg, dataet);


    %% 1. Prepare the data

    %% 1.1 Low Contrast

    %% 1.1.1 Create heatmaps for stimulus presentation period

    clear taskgaze basegaze gaze gazebl

    for trl = 1:length(dataetlow.trial)

        % Extract x position, y position, and pupil data for the first trial in datatask
        x_coords = dataetlow.trial{trl}(1,:);
        y_coords = dataetlow.trial{trl}(2,:);
        pupil_data = dataetlow.trial{trl}(3,:);

        % Remove blink indices
        valid_indices = pupil_data > BLINK_THRESHOLD;

        % Construct the heatmap using the filtered x and y position data for the first trial
        heatmap = histcounts2(x_coords(valid_indices), ...
            y_coords(valid_indices), ...
            [0:800], [0:600]);

        % Smooth the heatmap
        sigma = sig;
        smoothed_heatmap = imgaussfilt(heatmap, sigma);
        taskgaze(trl,:,:) = smoothed_heatmap';

    end

    %% 1.1.2. Reshape
    gaze.powspctrm = reshape(taskgaze, [numel(dataetlow.trial) 1 600 800]); % 4D data
    gaze.time = [1:800];
    gaze.freq = [1:600];
    gaze.dimord = 'rpt_chan_freq_time';
    gaze.label = {'et'};


    %% 1.1.3. Create heatmaps for baseline
    for trl = 1:length(dataetlowbl.trial)

        % Extract x position, y position, and pupil data for the first trial in datatask
        x_coords = dataetlowbl.trial{trl}(1,:);
        y_coords = dataetlowbl.trial{trl}(2,:);
        pupil_data = dataetlowbl.trial{trl}(3,:);

        % Remove blink indices
        valid_indices = pupil_data > BLINK_THRESHOLD;

        % Construct the heatmap using the filtered x and y position data for the first trial
        heatmap = histcounts2(x_coords(valid_indices), ...
            y_coords(valid_indices), ...
            [0:800], [0:600]);

        % Smooth the heatmap
        sigma = sig;
        smoothed_heatmap = imgaussfilt(heatmap, sigma);
        basegaze(trl,:,:)=smoothed_heatmap';

    end

    %% 1.1.4. Reshape
    gazebl.powspctrm = reshape(basegaze, [numel(dataetlowbl.trial) 1 600 800]); % 4D data
    gazebl.time = [1:800];
    gazebl.freq = [1:600];
    gazebl.dimord = 'rpt_chan_freq_time';
    gazebl.label = {'et'};
    figure; ft_singleplotTFR([],gazebl);


    %% 1.1.5. Rename the gaze variables
    if strcmp(smoothing, 'Granular')
        gazelow = gaze;
        gazelow.trialinfo = dataetlow.trialinfo;
        gazelowbl = gazebl;
        gazelowbl.trialinfo = dataetlow.trialinfo;
    elseif strcmp(smoothing, 'Coarse')
        gazelow_coarse = gaze;
        gazelow_coarse.trialinfo = dataetlow.trialinfo;
        gazelowbl_coarse = gazebl;
        gazelowbl_coarse.trialinfo = dataetlow.trialinfo;
    end


    %% 1.2. High Contrast

    %% 1.2.1. Create heatmaps for stimulus presentation period

    clear taskgaze basegaze gaze gazebl

    for trl=1:length(dataethigh.trial)

        % Extract x position, y position, and pupil data for the first trial in datatask
        x_coords = dataethigh.trial{trl}(1,:);
        y_coords = dataethigh.trial{trl}(2,:);
        pupil_data = dataethigh.trial{trl}(3,:);

        % Remove blink trials
        valid_indices = pupil_data > BLINK_THRESHOLD;

        % Construct the heatmap using the filtered x and y position data for the first trial
        heatmap = histcounts2(x_coords(valid_indices), ...
            y_coords(valid_indices), ...
            [0:800], [0:600]);

        % Smooth the heatmap
        sigma = sig;
        smoothed_heatmap = imgaussfilt(heatmap, sigma);
        taskgaze(trl,:,:) = smoothed_heatmap';

    end

    %% 1.2.2 Reshape

    gaze.powspctrm=reshape(taskgaze, [numel(dataethigh.trial) 1 600 800]);
    gaze.time = [1:800];
    gaze.freq = [1:600];
    gaze.dimord = 'rpt_chan_freq_time';
    gaze.label = {'et'};


    %% 1.2.3 Create heatmaps for baseline

    for trl = 1:length(dataethighbl.trial)

        % Extract x position, y position, and pupil data for the first trial in datatask
        x_coords = dataethighbl.trial{trl}(1,:);
        y_coords = dataethighbl.trial{trl}(2,:);
        pupil_data = dataethighbl.trial{trl}(3,:);

        % Remove blink trials
        valid_indices = pupil_data > BLINK_THRESHOLD;

        % Construct the heatmap using the filtered x and y position data for the first trial
        heatmap = histcounts2(x_coords(valid_indices), ...
            y_coords(valid_indices), ...
            [0:800], [0:600]);

        % Smooth the heatmap
        sigma = sig;
        smoothed_heatmap = imgaussfilt(heatmap, sigma);
        basegaze(trl,:,:) = smoothed_heatmap';

    end

    %% 1.2.4 Reshape

    gazebl.powspctrm=reshape(basegaze, [numel(dataethighbl.trial) 1 600 800]);
    gazebl.time = [1:800];
    gazebl.freq = [1:600];
    gazebl.dimord = 'rpt_chan_freq_time';
    gazebl.label = {'et'};
    figure; ft_singleplotTFR([],gazebl);


    %% 1.2.5 Rename the gaze variables

    if strcmp(smoothing, 'Granular')
        gazehigh_fine = gaze;
        gazehigh_fine.trialinfo = dataethigh.trialinfo;
        gazehighbl_fine = gazebl;
        gazehighbl_fine.trialinfo = dataethigh.trialinfo;
    elseif strcmp(smoothing, 'Coarse')
        gazehigh_coarse = gaze;
        gazehigh_coarse.trialinfo = dataethigh.trialinfo;
        gazehighbl_coarse = gazebl;
        gazehighbl_coarse.trialinfo = dataethigh.trialinfo;
    end


    %% 1.3. Save the new data frames

    if strcmp(smoothing, 'Granular')

        save gazehigh_fine gazehigh_fine gazehighbl_fine
        save gazelow_fine gazelow_fine gazelowbl_fine

    elseif strcmp(smoothing, 'Coarse')

        save gazehigh_coarse gazehigh_coarse gazehighbl_coarse
        save gazelow_coarse gazelow_coarse gazelowbl_coarse

    end

end % Loop over subjects



%% 2. Load the data

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)

    close all

    if strcmp(smoothing, 'Granular') && strcmp(variable, 'none')
        load gazehigh
        load gazelow
    elseif strcmp(smoothing, 'Coarse') && strcmp(variable, 'none')
        load gazehigh_coarse
        load gazelow_coarse
        gazehigh = gazehigh_coarse;
        gazelow = gazelow_coarse;
    elseif strcmp(variable, 'L_saccade')
        load gazehighET_gran_sacc
        load gazelowET_gran_sacc
        gazelow = gazelowET_gran;
        gazehigh = gazehighET_gran;
    elseif strcmp(variable, 'L_fixation')
        load gazehighET_gran_fix
        load gazelowET_gran_fix
        gazelow = gazelowET_gran;
        gazehigh = gazehighET_gran;
    end

    % Save in a structure
    gaLCl{subj} = ft_freqdescriptives([], gazelow);
    gaHCl{subj} = ft_freqdescriptives([], gazehigh);

    if strcmp(variable, 'none') % Only available for full eye tracker data
        gaLCbl{subj} = ft_freqdescriptives([],gazelowbl);
        gaHCbl{subj} = ft_freqdescriptives([],gazehighbl);
    end

end


%% 3. Calculate condition differences...

%% 3.1 ... per subject

diffsub = ga_gazeh;

for subj = 1:length(subjects)
    diffsub = gaLCl{subj};
    diffsub.powspctrm = gaHCl{subj}.powspctrm - gaLCl{subj}.powspctrm;
    differences{subj} = diffsub;
end

%% 3.1 ... across all subjects

% Grand averages
ga_gazel = ft_freqgrandaverage([], gaLCl{:});
ga_gazeh = ft_freqgrandaverage([], gaHCl{:});

% Calculate overall difference [high - low]
diff = ga_gazeh;
diff.powspctrm = ga_gazeh.powspctrm - ga_gazel.powspctrm;



%% 4. First Level Analysis

%% 4.1 High and Low contrast

close all

% Preparation
cfg             = [];
cfg.figure      = 'gcf';
cfg.zlim        = [-.5 .5];
cfg.figure      = 'gcf';
cfg.colorbartext = 'Average Gaze Frequency';

% Low Contrast

figure;
for sub = 1:length(subjects)

    s = str2num(subjects{sub});
    subplot(2,5,sub); ft_singleplotTFR(cfg, gaLCl{sub});
    set(gcf, 'Position', [0, 0, 1500, 465]);
    set(gca, 'YDir','reverse')
    set(gcf,'color','w');
    set(gca,'Fontsize', 15);
    xlabel('x [pix]');
    ylabel('y [pix]');
    if ismember(sub, [2 3 4 5 7 8 9 10])
        ylabel('')
        set(gca, 'YTickLabel', [])
    end
    %     if ismember(sub, [1 2 3 4 5])
    %         set(gca, 'XTickLabel', [])
    %         xlabel('')
    %     end
    set(gca, 'XTickLabel', [])
    xlabel('')
    ax = gca;
    ax.FontSize = 15;
    ax.TitleFontSizeMultiplier = 1.3;
    cfg.colormap = '*RdBu';
    colorbar('off')
    title(strcat('S', num2str(sub)))

    % Plot circle
    hold on
    plot(x, y, 'Color', circlecolor1, 'LineWidth', 1);

end

% Save
saveas(gcf, '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_all\Granular\Low_All', 'png');


% High contrast

figure;
for sub = 1:length(subjects)

    s = str2num(subjects{sub});
    subplot(2,5,sub); ft_singleplotTFR(cfg, gaHCl{sub});
    set(gcf, 'Position', [0, 0, 1500, 465]);
    set(gca, 'YDir','reverse')
    set(gcf,'color','w');
    set(gca,'Fontsize', 15);
    xlabel('x [pix]');
    ylabel('y [pix]');
    if ismember(sub, [2 3 4 5 7 8 9 10])
        ylabel('')
        set(gca, 'YTickLabel', [])
    end
    %     if ismember(sub, [1 2 3 4 5])
    %         set(gca, 'XTickLabel', [])
    %         xlabel('')
    %     end
    set(gca, 'XTickLabel', [])
    xlabel('')
    ax = gca;
    ax.FontSize = 15;
    ax.TitleFontSizeMultiplier = 1.3;
    cfg.colormap = '*RdBu';
    colorbar('off')
    title(strcat('S', num2str(sub)))

    % Plot circle
    hold on
    plot(x, y, 'Color', circlecolor1, 'LineWidth', 1);

end

% Save
saveas(gcf, '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_all\Granular\High_All', 'png');

%% 4.2. Plot Differences

close all

% Preparation
cfg             = [];
cfg.figure      = 'gcf';
cfg.zlim        = [-.1 .1];
cfg.colormap    = '*RdBu';
cfg.figure      = 'gcf';
cfg.colorbartext = 'Average Gaze Frequency';


% Plot and save as single figures

s = 0;

for sub = 1:length(subjects)

    s = s + 1;
    figure;
    set(gcf, 'Position', [0, 0, 920, 630]); % Specify the figure size
    ft_singleplotTFR(cfg, differences{sub});
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.FontSize = 28;
    ax.TitleFontSizeMultiplier = 1.2;
    set(gcf,'color','w');
    xlabel('x [pix]');
    ylabel('y [pix]');
    hold on
    plot(x, y, 'Color', circlecolor1, 'LineWidth', 1);
    title(strcat('S', ' ', num2str(sub)))

    % Save
    saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_all\Granular\Subjects\Sub', num2str(s)), 'png');

end



%% 5. Second Level Analyss

%% 5.1 Descriptive Analysis

%% 5.1.1 Plot heatmaps for contrast and the heatmap for the difference

% Preparation
cfg             = [];
cfg.figure      = 'gcf';
cfg.zlim        = [0 .7];
cfg.colormap    = '*RdBu';
cfg.colorbartext = 'Average Gaze Frequency';


% Low contrast
figure; ft_singleplotTFR(cfg, ga_gazel);

set(gcf, 'Position', [0, 0, 900, 610]);
set(gca, 'YDir','reverse')
ax = gca;
ax.FontSize = 28;
ax.TitleFontSizeMultiplier = 1.2;
set(gcf,'color','w');
xlabel('x [pix]')
ylabel('y [pix]')
title('')

% Plot circle
hold on
plot(x, y, 'Color', circlecolor1, 'LineWidth', 1);

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_all\', smoothing, '\Low'), 'png');



% High contrast
figure;ft_singleplotTFR(cfg, ga_gazeh);

set(gcf, 'Position', [0, 0, 900, 610]);
set(gca, 'YDir','reverse')
ax = gca;
ax.FontSize = 28;
ax.TitleFontSizeMultiplier = 1.2;
set(gcf,'color','w');
xlabel('x [pix]')
ylabel('y [pix]')
title('')

% Plot circle
hold on
plot(x, y, 'Color', circlecolor1, 'LineWidth', 1);

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_all\', smoothing, '\High'), 'png');


% Differencee
cfg.zlim        = [-.025 .025];
figure;ft_singleplotTFR(cfg, diff);
set(gcf, 'Position', [0, 0, 900, 610]);
set(gca, 'YDir','reverse')
ax = gca;
ax.FontSize = 28;
ax.TitleFontSizeMultiplier = 1.2;
set(gcf,'color','w');
xlabel('x [pix]')
ylabel('y [pix]')
title('')

% Plot circle
hold on
plot(x, y, 'Color', circlecolor1, 'LineWidth', 1);

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_all\', smoothing, '\Difference'), 'png');




%% 5.2 Statistical analysis

%% 5.2.1 Calculate difference for each subject

if strcmp(variable, 'none') % only possible if full eye tracker data

    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'subtract';
    for subj=1:length(subjects)
        gazel{subj}=ft_math(cfg,gaLCl{subj},gaLCbl{subj});
        gazeh{subj}=ft_math(cfg,gaHCl{subj},gaHCbl{subj});
    end

end


%% 5.2.2 Monte Carlo estimation

% Design preparationn
clear design
subj = numel(subjects);
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;


% Other Preparation
cfg                  = [];
cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;
cfg.spmversion       = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';

cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

cfg.neighbours       = [];


% Estimate

if strcmp(variable, 'none') % Only available for full eye tracker data

    [statLC] = ft_freqstatistics(cfg, gaLCl{:}, gaLCbl{:});
    [statHC] = ft_freqstatistics(cfg, gaHCl{:}, gaHCbl{:});

end

[stat] = ft_freqstatistics(cfg, gaHCl{:}, gaLCl{:});



%% 5.2.3 Visualization of significant cluster

% Mask out all non significant
stat.stat(isnan(stat.stat)) = 0;
if strcmp(variable, 'none') % Only available for full eye tracker data
    statLC.stat(isnan(statLC.stat)) = 0;
    statHC.stat(isnan(statHC.stat)) = 0;
end

% Preparation
cfg                 = [];
cfg.parameter       = 'stat';
cfg.maskparameter   = 'mask';
cfg.maskstyle       = 'outline';
cfg.figure          = 'gcf';
cfg.zlim            = [-4 4];
% cfg.zlim            = 'absmax';
cfg.colormap = '*RdBu';
cfg.colorbartext = 't-value';

% Plot
figure;
set(gcf, 'Position', [0, 0, 900, 610]);
ft_singleplotTFR(cfg, stat);
set(gca, 'YDir','reverse')
ax = gca;
ax.FontSize = 28;
ax.TitleFontSizeMultiplier = 1.5;
title('')
yticks([0 200 400 600])
ylabel('y [pix]')
xlabel('x [pix]')

% Add circle
hold on
plot(x, y, 'Color', circlecolor2, 'LineWidth', 1);

% Save
if strcmp(variable, 'none')
    saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_all\', smoothing, '\TTest_Difference'), 'png');
elseif strcmp(variable, 'L_saccade')
    saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_ET\Saccade\', smoothing, '\TTest_Difference'), 'png');
elseif strcmp(variable, 'L_fixation')
    saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_ET\Fixation\', smoothing, '\TTest_Difference'), 'png');
end


% Also plot cluster for task vs. baseline, if available (only for full eye tracker data)

if strcmp(variable, 'none')

    % Low contrast
    figure;
    ft_singleplotTFR(cfg, statLC);
    set(gca, 'YDir','reverse')
    set(gcf, 'Position', [0, 0, 900, 610]); % Specify the figure size
    ax = gca;
    ax.FontSize = 28;
    ax.TitleFontSizeMultiplier = 1.2;
    title('')
    ylabel('y [pix]')
    xlabel('x [pix]')

    % Add circle
    hold on
    plot(x, y, 'Color', circlecolor2, 'LineWidth', 1);

    % Save
    if strcmp(variable, 'none')
        saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_all\', smoothing, '\TTest_low'), 'png');
    elseif strcmp(variable, 'L_saccade')
        saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_ET\Saccade\', smoothing, '\TTest_low'), 'png');
    elseif strcmp(variable, 'L_fixation')
        saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_ET\Fixation\', smoothing, '\TTest_low'), 'png');
    end


    % High contrast
    figure;
    ft_singleplotTFR(cfg, statHC);
    set(gca, 'YDir','reverse')
    set(gcf, 'Position', [0, 0, 900, 610]); % Specify the figure size
    ax = gca;
    ax.FontSize = 28;
    ax.TitleFontSizeMultiplier = 1.2;
    title('')
    ylabel('y [pix]')
    xlabel('x [pix]')

    % Add circle
    hold on
    plot(x, y, 'Color', circlecolor2, 'LineWidth', 1);

    % Save
    if strcmp(variable, 'none')
        saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_all\', smoothing, '\TTest_high'), 'png');
    elseif strcmp(variable, 'L_saccade')
        saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_ET\Saccade\', smoothing, '\TTest_high'), 'png');
    elseif strcmp(variable, 'L_fixation')
        saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\Gaze_ET\Fixation\', smoothing, '\TTest_high'), 'png');
    end

end
