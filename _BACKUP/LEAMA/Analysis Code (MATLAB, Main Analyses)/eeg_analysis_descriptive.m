%% EEG: Descriptive Analysis %%

%% Description %%
% This script analyses the EEG data of the ten subjects descriptively. It includes the generation of
% TFR spectra, topographic plots and power spectra across subjects and for individual subjects.
% The type of data that should be analysed can be specified.


%% Script Structure
% 0. Preparation
% 1. FFT and average over all trials of each condition
% 2. Baseline correction (if necessary)
% 3. Grand average (over participants)
% 4. Multiplots for exploration
% 5. Select electrodes
% 6. First Level Analysis
% 7. Second Level Analysis



%% 0. Preparation

clear all
close all

%% 0.1 Add clean version of Fieldtrip
cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin')
run startup

%% 0.2 Define data types
% Data type: 'manual_preprocessing', 'automagic', 'automagic_opticat'
datatype = 'automagic_opticat';

% Taper type: 'singletaper', 'multitaper'
tapertype = 'singletaper';

% Time window: 'baseline', 'stimuli'
timewindow = 'baseline';

% Exclusion criteria: 'blk' 'blkET' 'deviation' 'no'
exclusion = 'no';

% Specify folder according to selected analyses settings
if strcmp(datatype, 'automagic_opticat') && strcmp(tapertype, 'singletaper') && strcmp(exclusion, 'no') && strcmp(timewindow, 'stimuli')
    folder = 'Singletaper';
elseif strcmp(datatype, 'automagic') && strcmp(tapertype, 'singletaper') && strcmp(exclusion, 'no') && strcmp(timewindow, 'stimuli')
    folder = 'Singletaper\WithoutOPTICAT';
elseif strcmp(datatype, 'automagic_opticat') && strcmp(tapertype, 'multitaper') && strcmp(exclusion, 'no') && strcmp(timewindow, 'stimuli')
    folder = 'Multitaper';
elseif strcmp(datatype, 'automagic') && strcmp(tapertype, 'multitaper') && strcmp(exclusion, 'no') && strcmp(timewindow, 'stimuli')
    folder = 'Multitaper\WithoutOPTICAT';
elseif strcmp(exclusion, 'blk')
    folder = 'BlinkExclusion';
elseif strcmp(exclusion, 'blkET')
    folder = 'BlinkExclusionET';
elseif strcmp(exclusion, 'deviation')
    folder = 'DeviationExclusion';
elseif strcmp(timewindow, 'baseline')
    folder = 'Baseline';
end

%% 0.3 Participants and path

% Participants
if strcmp(exclusion, 'blkET')
    subjects = {'1' '2' '3' '6' '8' '9' '11'}; % Only select participants with existing data after exclusion
elseif strcmp(exclusion, 'deviation')
    subjects = {'1' '2' '3' '4' '5' '6' '8' '9' '11'}; % Only select participants with existing data after exclusion
else
    subjects = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '11'};
end

% Path
if strcmp(datatype, 'manual_preprocessing')
    path = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\data\';
elseif strcmp(datatype,'automagic')
    path = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_results\';
elseif strcmp(datatype,'automagic_opticat')
    path = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results\';
end


%% 1. FFT and average over all trials of each condition

for subj = 1:length(subjects)

    %% 1.1 Load data
    %     keep subj subjects path
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    load data


    %% 1.2 Identify indices of trials belonging to conditions

    % High contrast
    indhc_hor = find(data.trialinfo == 0);
    indhc_ver = find(data.trialinfo == 90);
    indhc_45 = find(data.trialinfo == 45);
    indhc_115 = find(data.trialinfo == 115);

    % Low contrast
    indlc_hor = find(data.trialinfo == 1);
    indlc_ver = find(data.trialinfo ==  91);
    indlc_45 = find(data.trialinfo == 46);
    indlc_115 = find(data.trialinfo == 116);


    %% 1.3 Fast Fourier Transform (automatically averaging over trials)

    if strcmp(tapertype, 'singletaper')

        % With single taper
        cfg             = [];
        cfg.output      = 'pow';
        cfg.method      = 'mtmconvol';
        cfg.taper       = 'hanning';
        cfg.foi         = 4:1:120;                        % Analysis 2 to 30 Hz in steps of 1 Hz
        cfg.t_ftimwin   = ones(length(cfg.foi), 1).*0.5;  % Length of time window = 0.5 sec
        cfg.toi         = -1:0.05:3;                      % What time window we are interested in, rsp. to stim presentation
        cfg.keeptrials  = 'no';

        % FFT
        cfg.trials = [indhc_hor; indhc_ver; indhc_45; indhc_115];
        tfr_hc = ft_freqanalysis(cfg, data);
        cfg.trials = [indlc_hor; indlc_ver; indlc_45; indlc_115];
        tfr_lc = ft_freqanalysis(cfg, data);
        tfr_lc.cfg = [];
        tfr_hc.cfg = [];

        % Save
        save tfr_lc tfr_lc
        save tfr_hc tfr_hc

    elseif strcmp(tapertype, 'multitaper')
        % With multiple tapers
        cfg             = [];
        cfg.output      = 'pow';
        cfg.method      = 'mtmconvol';
        cfg.taper       = 'dpss';
        cfg.foi         = 30:1:90;                          % Analysis 30 to 90 Hz in steps of 1 Hz
        cfg.tapsmofrq   = 10;                               % 10 tapers
        cfg.t_ftimwin   = ones(length(cfg.foi),1).*0.5;     % Length of time window = 0.5 sec
        cfg.toi         = -1:0.05:3;                        % What time window we are interested in, rsp. to stim presentation
        cfg.keeptrials  = 'no';

        % FFT
        cfg.trials      = [indhc_hor; indhc_ver; indhc_45; indhc_115];
        tfr_hc_mul = ft_freqanalysis(cfg, data);
        cfg.trials      = [indlc_hor; indlc_ver; indlc_45; indlc_115];
        tfr_lc_mul= ft_freqanalysis(cfg, data);

        % Save
        tfr_lc_mul.cfg  = [];
        tfr_hc_mul.cfg  = [];
        save tfr_lc_mul tfr_lc_mul
        save tfr_hc_mul tfr_hc_mul
    end

end


%% 2. Baseline correction

for subj = 1:length(subjects)

    %% 2.1 Load correct data
    datapath = strcat(path, subjects{subj});
    disp(datapath)
    cd(datapath)

    if strcmp(tapertype, 'singletaper') && strcmp(exclusion, 'no')
        load tfr_lc
        load tfr_hc
    elseif strcmp(tapertype, 'multitaper') && strcmp(exclusion, 'no')
        load tfr_lc_mul
        load tfr_hc_mul
        tfr_lc = tfr_lc_mul;
        tfr_hc = tfr_hc_mul;
    elseif strcmp(exclusion, 'blk')
        load tfr_lc_nb
        load tfr_hc_nb
    elseif strcmp(exclusion, 'blkET')
        load tfr_lc_nbET
        load tfr_hc_nbET
        tfr_lc = tfr_lc_nbET;
        tfr_hc = tfr_hc_nbET;
    elseif strcmp(exclusion, 'deviation')
        load tfr_lc_nodev
        load tfr_hc_nodev
        tfr_lc = tfr_lc_nodev;
        tfr_hc = tfr_hc_nodev;
    end

    disp(subj)

    %% 2.2 Select baseline and calculate difference
    cfg             = [];
    cfg.baseline    = [-Inf -.25];
    cfg.baselinetype = 'db';

    if strcmp(timewindow, 'stimuli') % if 'baseline', do not correct for baseline
        if strcmp(tapertype, 'singletaper')
            tfr_hc = ft_freqbaseline(cfg, tfr_hc);
            tfr_lc = ft_freqbaseline(cfg, tfr_lc);
        elseif strcmp(tapertype, 'multitaper')
            tfr_hc = ft_freqbaseline(cfg, tfr_hc_mul);
            tfr_lc = ft_freqbaseline(cfg, tfr_lc_mul);
        end
    end

    %% 2.3 Save baseline-corrected data for each participant
    lc{subj} = tfr_lc;
    hc{subj} = tfr_hc;

end % End participant loop


%% 3. Compute grand average (over participants)
cfg             = [];
cfg.keepindividual = 'yes';
gahc = ft_freqgrandaverage(cfg, hc{:});
galc = ft_freqgrandaverage(cfg, lc{:});


%% 4. Multiplot baseline differences for both contrasts

% Load layout
cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\headmodel')
load layout

cfg             = [];
cfg.layout      = layout;
cfg.zlim = [-3 3];
cfg.figure      = 'gcf';
cfg.colormap = '*RdBu';

% Figures without line noise
figure; ft_multiplotTFR(cfg, gahc);
figure; ft_multiplotTFR(cfg, galc);


%% 5. Select electrodes

% Only execute this for the data that was computed using OPTICAT and single tapers for FFT without exclusions
if strcmp(tapertype, 'singletaper') && strcmp(datatype, 'automagic_opticat') && strcmp(timewindow, 'stimuli')  && strcmp(exclusion, 'no')

    %% 5.1 Contrast conditions

    % Low Contrast: Average the power over the gamma range for each electrode
    counter = 0;

    for chnnl = 1:length(galc.label) % For all channels
        cfg = [];
        cfg.channel = char(galc.label(chnnl));
        cfg.frequency = [30 90];        % Gamma range
        cfg.time = [0 2];               % Stimulus presentation period
        powerspectrum = ft_freqdescriptives(cfg, galc);
        meanPower = mean(powerspectrum.powspctrm(:)); % Calculate mean power
        if meanPower > 0
            counter = counter + 1;
            chnnls(counter, 1) = galc.label(chnnl);
            chnnls_means(counter, 2) = meanPower;
        end
    end

    % Select parietal and occipital electrodes with positive gamma power, on average
    indicesPlow = find(cellfun(@(str) startsWith(str, 'P'), chnnls));
    indicesOlow = find(cellfun(@(str) startsWith(str, 'O'), chnnls));
    indicesPOlow = [indicesPlow; indicesOlow];
    electrodeslow = chnnls(indicesPOlow);


    % High Contrast: Average the power over the specified frequency range for each electrode
    clear powerspectrum meanPower chnnls chnnls_means
    counter = 0;

    for chnnl = 1:length(gahc.label)
        cfg = [];
        cfg.channel = char(gahc.label(chnnl));
        cfg.frequency = [30 90];        % Gamma range
        cfg.time = [0 2];               % Stimulus presentation period
        powerspectrum = ft_freqdescriptives(cfg, gahc);
        meanPower = mean(powerspectrum.powspctrm(:));
        if meanPower > 0
            counter = counter + 1;
            chnnls(counter, 1) = gahc.label(chnnl);
            chnnls_means(counter, 2) = meanPower; % Calculate mean power
        end
    end

    % Select parietal and occipital electrodes with positive gamma power, on average
    indicesPhigh = find(cellfun(@(str) startsWith(str, 'P'), chnnls));
    indicesOhigh = find(cellfun(@(str) startsWith(str, 'O'), chnnls));
    indicesPOhigh = [indicesPhigh; indicesOhigh];
    electrodeshigh = chnnls(indicesPOhigh);

    %% 5.2 Common electrodes

    % Select electrodes with positive gamma power in both conditions
    electrodes = intersect(electrodeslow, electrodeshigh);
    save '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\scripts\electrodes' electrodes

else % For all other analyses, load electrodes determined with the singletaper opticat data
    electrodesstruct = load('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\scripts\electrodes');
    electrodes = electrodesstruct.electrodes;
end


%% 6. First Level Analysis

%% 6.1 Time-frequency spectra for all subjects

% Load layout
cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\headmodel')
load layout

% Preparation
cfg             = [];
cfg.channel     = electrodes; % Selected electrodes
cfg.layout      = layout;
cfg.figure      = 'gcf';
cfg.colormap    = '*RdBu';
cfg.zlim        = [-3 3];
if strcmp(timewindow, 'stimuli')
    cfg.colorbartext = 'dB';
elseif strcmp(timewindow, 'baseline')
    cfg.colorbartext = '\muV{^2}/Hz';
end

if strcmp(timewindow, 'stimuli')
    cfg.xlim        = [-.25 2];
elseif strcmp(timewindow, 'baseline')
    cfg.xlim        = [-.1 1];
end

% Low contrast
figure;
for sub = 1:length(subjects)
    s = str2num(subjects{sub});
    subplot(2,5,sub); ft_singleplotTFR(cfg, lc{sub});
    set(gcf, 'Position', [0, 0, 1500, 500]); % Specify the figure size
    set(gcf,'color','w');
    set(gca,'Fontsize', 15);
    xlabel('Time [sec]');
    ylabel('Frequency [Hz]');
    if ismember(sub, [2 3 4 5 7 8 9 10])
        ylabel('')
        set(gca, 'YTickLabel', [])
    end
    if ismember(sub, [1 2 3 4 5])
        set(gca, 'XTickLabel', [])
        xlabel('')
    end
    ax = gca;
    ax.FontSize = 15;
    ax.TitleFontSizeMultiplier = 1.3;
    cfg.colormap = '*RdBu';
    colorbar('off')
    title(strcat('S', num2str(s)))
end

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\TFR_Low_All'), 'png');


% High Contrast
figure;
for sub = 1:length(subjects)
    s = str2num(subjects{sub});
    subplot(2,5,sub); ft_singleplotTFR(cfg, hc{sub});
    set(gcf, 'Position', [0, 0, 1500, 500]); % Specify the figure size
    set(gcf,'color','w');
    set(gca,'Fontsize', 15);
    xlabel('Time [sec]');
    ylabel('Frequency [Hz]');
    if ismember(sub, [2 3 4 5 7 8 9 10])
        ylabel('')
        set(gca, 'YTickLabel', [])
    end
    if ismember(sub, [1 2 3 4 5])
        set(gca, 'XTickLabel', [])
        xlabel('')
    end
    colorbar('off')
    ax = gca;
    ax.FontSize = 15;
    ax.TitleFontSizeMultiplier = 1.3;
    cfg.colormap = '*RdBu';
    title(strcat('S', num2str(s)))
end

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\TFR_High_All'), 'png');


%% 6.2 Power spectra for all participants

%% 6.2.1 Calculate power for all participants

cfg             = [];
cfg.latency     = [0 2];
cfg.frequency   = [30 120];
cfg.avgovertime = 'yes';

for sub = 1:length(subjects)

    rpow_lc{sub} = ft_selectdata(cfg, lc{sub});
    rpow_hc{sub} = ft_selectdata(cfg, hc{sub});
    rpow_lc{sub} = rmfield(rpow_lc{sub}, 'time');
    rpow_hc{sub} = rmfield(rpow_hc{sub}, 'time');

    rpow_lc{sub}.dimord = 'chan_freq';
    rpow_hc{sub}.dimord = 'chan_freq';
end

%% 6.2.2 Plot

% Load layout
cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\headmodel')
load layout

% Preparation
cfg                     = [];
cfg.figure              = 'gcf';
cfg.layout              = layout;
cfg.channel             = electrodes;
figure;

% Plot
for sub = 1:length(subjects)

    subplot(2,5,sub);
    ft_singleplotER(cfg, rpow_lc{1,sub}, rpow_hc{1,sub});
    hold on;

    % Add shadedErrorBar path
    addpath('C:\Users\dummy\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\raacampbell_shadedErrorBar')

    % Add the standard error bars
    channels = ismember(rpow_lc{subj}.label, cfg.channel);
    lowbar = shadedErrorBar(rpow_lc{1,sub}.freq, mean(rpow_lc{1,sub}.powspctrm(channels, :), 1), ...
        std(rpow_lc{1,sub}.powspctrm(channels, :))/sqrt(size(rpow_lc{1,sub}.powspctrm(channels, :), 1)), 'lineProps', {'b', 'MarkerFaceColor', 'b'});
    highbar = shadedErrorBar(rpow_hc{1,sub}.freq, mean(rpow_hc{1,sub}.powspctrm(channels, :), 1), ...
        std(rpow_hc{1,sub}.powspctrm(channels, :))/sqrt(size(rpow_hc{1,sub}.powspctrm(channels, :), 1)), 'lineProps', {'r', 'MarkerFaceColor', 'r'});

    % Settings
    set(gcf, 'Position', [0, 0, 1950, 480]); % Specify the figure size
    set(gcf, 'color', 'none');
    ax = gca;
    ax.FontSize = 15;
    ax.TitleFontSizeMultiplier = 1.3;

    % Xlim
    xlim([30 120])
    xlabel('Frequency [Hz]');
    ylabel('Power [dB]');
    xticks([40 60 80 100 120]);
    xticklabels({'40','60','80','100','120'});

    % Ylim
    max_high = max(mean(rpow_hc{1,sub}.powspctrm(channels, :), 1));
    max_low = max(mean(rpow_lc{1,sub}.powspctrm(channels, :), 1));
    min_high = min(mean(rpow_hc{1,sub}.powspctrm(channels, :), 1));
    min_low = min(mean(rpow_lc{1,sub}.powspctrm(channels, :), 1));
    max_val = max(max_high, max_low);
    min_val = min(min_high, min_low);
    set(gca, 'YLim', [min_val+0.6*min_val max_val+0.3*max_val]);

    if ismember(sub, [2 3 4 5 7 8 9 10])
        ylabel('')
    end
    if ismember(sub, [1 2 3 4 5])
        set(gca, 'XTickLabel', [])
        xlabel('')
    end

    title(strcat('S', num2str(sub)))
    box on

end

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\Power_All'), 'png');


%% 6.3 Time-frequency spectra of contrast effect for all participants

%% 6.3.1 Calculate difference for all participants

for subj = 1:length(subjects)
    diffall{subj} = hc{subj};
    diffall{subj}.powspctrm = (hc{subj}.powspctrm - lc{subj}.powspctrm);
end

%% 6.3.2 Plot

% Preparation
cfg             = [];
cfg.channel     = electrodes;
cfg.layout      = layout;
cfg.figure      = 'gcf';
cfg.colormap    = '*RdBu';
cfg.zlim        = [-3 3];
if strcmp(timewindow, 'stimuli') % Adapt z axis text
    cfg.colorbartext = 'dB';
elseif strcmp(timewindow, 'baseline')
    cfg.colorbartext = '\muV{^2}/Hz';
end
if strcmp(timewindow, 'stimuli') % Adapt x axis limits
    cfg.xlim        = [-.25 2];
elseif strcmp(timewindow, 'baseline')
    cfg.xlim        = [-.1 1];
end

figure;
for sub = 1:length(subjects)
    s = str2num(subjects{sub});
    subplot(2,5,sub); ft_singleplotTFR(cfg, diffall{sub});
    set(gcf, 'Position', [0, 0, 1500, 500]);
    set(gcf,'color','w');
    set(gca,'Fontsize', 15);
    xlabel('Time [sec]');
    ylabel('Frequency [Hz]');
    if ismember(sub, [2 3 4 5 7 8 9 10])
        ylabel('')
        set(gca, 'YTickLabel', [])
    end
    if ismember(sub, [1 2 3 4 5])
        set(gca, 'XTickLabel', [])
        xlabel('')
    end
    ax = gca;
    ax.FontSize = 15;
    ax.TitleFontSizeMultiplier = 1.3;
    cfg.colormap = '*RdBu';
    colorbar('off')
    title(strcat('S', num2str(s)))
end

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\Diff_All'), 'png');



%% 7. Second Level Analysis

%% 7.1 Time-frequency spectra relative to baseline

% Load layout
cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\headmodel')
load layout

% Preparation
close all
cfg             = [];
cfg.layout      = layout;
cfg.channel     = electrodes; % Selected channels
cfg.colormap = '*RdBu';

if strcmp(tapertype, 'singletaper') && strcmp(timewindow, 'stimuli')
    cfg.zlim        = [-2 2];
elseif strcmp(tapertype, 'singletaper') && strcmp(timewindow, 'stimuli') && strcmp(exclusion, 'deviation')
    cfg.zlim        = [-2 2];
elseif strcmp(tapertype, 'multitaper') && strcmp(timewindow, 'stimuli')
    cfg.zlim        = [-1.5 1.5];
elseif strcmp(timewindow, 'baseline')
end

if strcmp(timewindow, 'stimuli')
    cfg.xlim        = [-.25 2];
elseif strcmp(timewindow, 'baseline')
    cfg.xlim        = [-1 1];
    cfg.zlim        = [0 1];
end

cfg.figure      = 'gcf';
if strcmp(timewindow, 'stimuli')
    cfg.colorbartext = 'dB';
elseif strcmp(timewindow, 'baseline')
    cfg.colorbartext = '\muV{^2}/Hz';
end


% Low contrast
figure; ft_singleplotTFR(cfg, galc);
set(gcf, 'Position', [0, 0, 810, 600]); % Specify the figure size
set(gcf,'color','w');
set(gca,'Fontsize', 20);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
ax = gca;
ax.FontSize = 25;
ax.TitleFontSizeMultiplier = 1.3;
title('');

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\TFR_Low'), 'png');


% High contrast
figure; ft_singleplotTFR(cfg, gahc);
set(gcf, 'Position', [0, 0, 810, 600]); % Specify the figure size
set(gca,'Fontsize', 20);
xlabel('Time [sec]');
ylabel('Frequency [Hz]', 'Color', 'white');
ax = gca;
ax.FontSize = 25;
ax.TitleFontSizeMultiplier = 1.3;
title('');

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\TFR_High'), 'png');


%% 7.2 Time-frequency spectrum of differences between contrasts

%% 7.2.1 Calculate difference

diff = gahc;
diff.powspctrm = (gahc.powspctrm - galc.powspctrm);

%% 7.2.2 Plot

% Preparation
cfg             = [];
cfg.layout      = layout;
cfg.channel     = electrodes; % Selected channels
cfg.figure      = 'gcf';
cfg.colormap = '*RdBu';
cfg.zlim        = [-1 1];

if strcmp(timewindow, 'stimuli')
    cfg.xlim        = [-.25 2];
elseif strcmp(timewindow, 'baseline')
    cfg.xlim        = [-1 1];
end

if strcmp(tapertype, 'multitaper') && strcmp(timewindow, 'stimuli')
    cfg.zlim        = [-1.5 1.5];
elseif strcmp(timewindow, 'baseline')
    cfg.zlim        = [-0.01 0.01];
end

if strcmp(timewindow, 'stimuli')
    cfg.colorbartext = 'dB';
elseif strcmp(timewindow, 'baseline')
    cfg.colorbartext = '\muV{^2}/Hz';
end

% Plot
figure;ft_singleplotTFR(cfg, diff);
set(gcf, 'Position', [0, 0, 810, 600]);
set(gcf, 'color', 'w');
set(gca, 'Fontsize', 20);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
if strcmp(timewindow, 'baseline')
    zticks([-0.3, 0, 0.03])
end
ax = gca;
ax.FontSize = 25;
ax.TitleFontSizeMultiplier = 1.3;
title('');

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\TFR_Diff'), 'png');


%% 7.3 Power spectra relative to baseline for the two conditions

%% 7.3.1 Calculate the power spectra information

cfg             = [];
cfg.latency     = [0 2];
cfg.frequency   = [30 120];
cfg.avgovertime = 'yes'; % Average across stimulus presentation period

relpow_lc = ft_selectdata(cfg, galc);
relpow_hc = ft_selectdata(cfg, gahc);
relpow_lc = rmfield(relpow_lc, 'time');
relpow_hc = rmfield(relpow_hc, 'time');

relpow_lc.dimord = 'subj_chan_freq';
relpow_hc.dimord = 'subj_chan_freq';

% Difference (not necessary)
diff = relpow_hc;
diff.powspctrm = relpow_hc.powspctrm - relpow_lc.powspctrm;


%% 7.3.2 Plot

% Run startup again
cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin')
run startup

% Prepare data averaged across electrodes for shaded error bars
cfg             = [];
cfg.channel     = electrodes;
cfg.avgoverchan = 'yes'; % Average over electrodes defined in cfg.channel
tmplow = ft_selectdata(cfg, relpow_lc);
tmphigh = ft_selectdata(cfg, relpow_hc);

% Load layout
cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\headmodel')
load layout

cfg             = [];
cfg.layout      = layout;
cfg.channel     = electrodes;
cfg.avgoverchan = 'yes';
cfg.figure      = 'gcf';
figure; ft_singleplotER(cfg, relpow_lc, relpow_hc);

% Settings
set(gcf, 'color', 'w');
set(gca, 'Fontsize', 30);
title('')
hold on;

% Add shadedErrorBar path
addpath('C:\Users\dummy\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\raacampbell_shadedErrorBar')

% Add errorbar
clf
lowbar = shadedErrorBar(tmplow.freq, mean(tmplow.powspctrm, 1), ...
    std(tmplow.powspctrm)/sqrt(size(tmplow.powspctrm, 1)), 'lineprops', {'b'}, 'patchSaturation', 0.05);
set(lowbar.edge,'LineStyle','none');
lowbar.mainLine.LineWidth = 1;
hold on;
highbar = shadedErrorBar(tmphigh.freq, mean(tmphigh.powspctrm, 1), ...
    std(tmphigh.powspctrm)/sqrt(size(tmphigh.powspctrm, 1)), 'lineprops', {'r', 'markerfacecolor', 'r'}, 'patchSaturation', 0.05);
set(highbar.edge,'LineStyle','none')
highbar.mainLine.LineWidth = 1;

% More settings
set(gcf, 'Position', [0, 0, 820, 600]); % Specify the figure size
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
xticks([30 40 50 60 70 80 90]);
% xticklabels({'30', '40', '50','60','70', '80', '90'});
xlim([30 90]);
ylim([-2 1.4]);
if strcmp(exclusion, 'deviation')
    ylim([-2 1.4]);
elseif strcmp(exclusion, 'blkET')
    ylim([-2 2]);
elseif strcmp(exclusion, 'blk')
    ylim([-2 1.4]);
else
    ylim([-1 1.4]);
end

ax = gca;
ax.FontSize = 25;
ax.TitleFontSizeMultiplier = 1.3;

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\Power'), 'png');


%% 7.3.3 Calculate average GFP and GPP for both condition averages

meanPowerh = mean(tmphigh.powspctrm, 1);
meanPowerh = squeeze(meanPowerh);
[maxPowerh, maxFreqIndexh] = max(meanPowerh);
maxFrequencyh = tmphigh.freq(maxFreqIndexh);

meanPowerl = mean(tmplow.powspctrm, 1);
meanPowerl = squeeze(meanPowerl);
[maxPowerl, maxFreqIndexl] = max(meanPowerl);
maxFrequencyl = tmplow.freq(maxFreqIndexl);


%% 7.4 Topoplots

% Preparation
cfg             = [];
cfg.layout      = layout;
cfg.comment     = 'no';
cfg.gridscale   = 300;
cfg.figure      = 'gcf';
cfg.zlim        = [-0.8 0.8];
cfg.marker      = 'off';
cfg.highlight   = 'labels';
cfg.highlightchannel    = electrodes;
cfg.highlightsymbol     = '.';
cfg.highlightsize       = 20;
cfg.highlightfontsize  = 10;
cfg.colorbar = 'yes';
cfg.colormap = '*RdBu';

if strcmp(timewindow, 'stimuli')
    cfg.colorbartext = 'dB';
elseif strcmp(timewindow, 'baseline')
    cfg.colorbartext = '\muV{^2}/Hz';
end


%  Low contrast
cfg.xlim        = [30 90];
figure; ft_topoplotER(cfg, relpow_lc);
set(gcf, 'Position', [0, 0, 800, 600]); % Specify the figure size
set(gcf, 'color','w');
set(gca, 'Fontsize', 25);
title('')
cfg.highlightfontsize  = 10;
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\Topo_Low'), 'png');

% High contrast
cfg.xlim        = [30 90];
figure; ft_topoplotER(cfg, relpow_hc);
set(gcf, 'Position', [0, 0, 800, 600]); % Specify the figure size
set(gca, 'Fontsize', 25);
title('')
cfg.highlightfontsize  = 10;
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\Topo_High'), 'png');

