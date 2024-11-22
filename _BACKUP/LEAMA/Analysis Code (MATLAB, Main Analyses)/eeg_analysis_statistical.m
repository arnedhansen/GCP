%% EEG: Statistical Analysis %%

%% Description %%
% This script analyses the EEG data of the ten subjects statistically. It includes the performance of cluster-based 
% permutation t-tests on the TFR difference between the contrasts and the power spectra, as well as two-sided t-tests 
% for dependent samples for the individual GPP and GPF values.


%% Script Structure
% 0. Preparation
% 1. Baseline correction (if necessary)
% 2. Grand average (over participants)
% 3. Average over stimulus presentation period
% 4. Calculate difference of power spectra
% 5. Statistical Analysis: TFR and ER
% 6. Statistical Analysis: Timelocked data
% 7. Statistical Analysis: GPP and GPF


%% 0. Preparation

clear all
close all

%% 0.1 Add clean version of Fieldtrip
cd('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin')
run startup

%% 0.2 Define data type
% 'blk', 'automagic', 'automagic_opticat'
datatype = 'automagic_opticat';

% 'singletaper', 'multitaper'
tapertype = 'singletaper';

% 'baseline', 'stimuli'
timewindow = 'stimuli';

% 'blk' 'blkET' 'no' 'deviation'
exclusion = 'no';

% Specify folder
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

%% 0.3 Path and subjects
if strcmp(exclusion, 'blkET')
    subjects = {'1' '2' '3' '6' '8' '9' '11'};
elseif strcmp(exclusion, 'deviation')
    subjects = {'1' '2' '3' '4' '5' '6' '8' '9' '11'};
else
  subjects = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '11'};  
end


if strcmp(datatype, 'blk')
    path = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\data\';
elseif strcmp(datatype,'automagic')
    path = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_results\';
elseif strcmp(datatype,'automagic_opticat')
    path = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results\';
end

%% 0.4 Define electrodes
electrodesstruct = load('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\scripts\electrodes');
electrodes = electrodesstruct.electrodes;


%% 1. Baseline correction
for subj = 1:length(subjects)
    
    % Load data
    datapath = strcat(path, subjects{subj});
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
    
    % Baseline correction
    cfg                 = [];
    cfg.baseline        = [-Inf -.25];
    cfg.baselinetype    = 'db';
    
    if strcmp(timewindow, 'stimuli')
            tfr_hc = ft_freqbaseline(cfg, tfr_hc);
            tfr_lc = ft_freqbaseline(cfg, tfr_lc);
    end
    
    % Select 2.5 seconds
    cfg                 = [];
    
    if strcmp(timewindow, 'stimuli')
        cfg.latency         = [-.5 2];
    elseif strcmp(timewindow, 'baseline')
        cfg.latency         = [-1 1];
    end
    
    tfr_hc = ft_selectdata(cfg, tfr_hc);
    tfr_lc = ft_selectdata(cfg, tfr_lc);
    
    % Save each subject structure to a common structure
    HC{subj} = tfr_hc;
    LC{subj} = tfr_lc;
end


%% 2. Grand average

gahc = ft_freqgrandaverage([],HC{:});
galc  = ft_freqgrandaverage([],LC{:});


%% 3. Average over stimulus presentation period

cfg                     = [];
cfg.latency             = [0 2]; 
cfg.avgovertime         = 'yes';
cfg.frequency           = [30 120];

% Create extra datasets
for subj=1:length(subjects)
    tmp = ft_selectdata(cfg, HC{subj});
    tmp.dimord = 'chan_freq';
    tmp = rmfield(tmp,'time');
    HCtlk{subj} = tmp;
    
    tmp = ft_selectdata(cfg, LC{subj});
    tmp.dimord = 'chan_freq';
    tmp = rmfield(tmp,'time');
    LCtlk{subj} = tmp;
end

% Average frequency powers for time and frequencies specified above
gahc_tlk = ft_freqgrandaverage([], HCtlk{:});
galc_tlk  = ft_freqgrandaverage([], LCtlk{:});


%% 4. Calculate difference between contrast conditions

diff = gahc;
diff.powspctrm = (gahc.powspctrm - galc.powspctrm);



%% 5. Statistical analysis

%% 5.1 Preparation: Average across all selected electrodes for each individual subject

% TFR
cfg             = [];
cfg.channel     = electrodes;
cfg.avgoverchan = 'yes'; % Average over electrodes defined in cfg.channel
allchan_lowTFR = cell(1, length(subjects));  % Initialize cell arrays
allchan_highTFR = cell(1, length(subjects));
for sub = 1:length(subjects)
    allchan_lowTFR{1, sub} = ft_selectdata(cfg, LC{1,sub});
    allchan_highTFR{1, sub} = ft_selectdata(cfg, HC{1,sub});
end

% ER 
cfg             = [];
cfg.channel     = electrodes;
cfg.avgoverchan = 'yes'; % Average over electrodes defined in cfg.channel
allchan_lowER = cell(1, length(subjects));  % Initialize cell arrays
allchan_highER = cell(1, length(subjects));
for sub = 1:length(subjects)
    allchan_lowER{1, sub} = ft_selectdata(cfg, LCtlk{1,sub});
    allchan_highER{1, sub} = ft_selectdata(cfg, HCtlk{1,sub});
end


%% 5.2 Statistical analysis

% Define type of data that is analysed: 'TFR', 'ER'
type  = 'TFR';

% Define if averaged across electrodes: 'True', 'False'
avg = 'True';


% Design matrix
clear design
subj = numel(subjects);
design = zeros(2,2*subj);

for i = 1:subj
    design(1, i)        = i;
end
for i = 1:subj
    design(1, subj+i)   = i;
end

design(2, 1:subj)       = 1;
design(2, subj + 1:2*subj) = 2;


% Prepare cfg structure for Monte Carlo cluster analysis with frequency range, electrodes, and number of iterations
cfg                     = [];
cfg.design              = design;
cfg.uvar                = 1;
cfg.ivar                = 2;
cfg.spmversion          = 'spm12';
cfg.method              = 'montecarlo'; % Monte Carlo method
cfg.statistic           = 'ft_statfun_depsamplesT'; % Dependent samples t-test
cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum'; % Sum of t-values is the test value for a cluster
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.05;
cfg.numrandomization    = 1000;
cfg.frequency           = [30 120];  % For gamma frequency range
cfg.neighbours          = [];
if strcmp(avg, 'True') % Otherwise, select electrode
    cfg.channel  = electrodes;
    cfg.avgoverchan='yes';
end


% Calculate statistical parameters according to the data
if strcmp(tapertype,'multitaper') || strcmp(datatype, 'automagic')
    
    cfg.channel  = 'PPO2';
    chnl = 1;
    
    if strcmp(type, 'TFR') && strcmp(avg, 'False')
        [stat(chnl)] = ft_freqstatistics(cfg,  HC{:}, LC{:});
    elseif strcmp(type, 'ER') && strcmp(avg, 'False')
        [stat(chnl)] = ft_freqstatistics(cfg, HCtlk{:}, LCtlk{:});
    elseif strcmp(type, 'TFR') && strcmp(avg, 'True')
        [stat] = ft_freqstatistics(cfg, allchan_highTFR{:}, allchan_lowTFR{:});
    elseif strcmp(type, 'ER') && strcmp(avg, 'True')
        [stat] = ft_freqstatistics(cfg, allchan_highER{:}, allchan_lowER{:});
    end
    
elseif strcmp(timewindow,'baseline')
    
    cfg.channel  = electrodes;
    cfg.avgoverchan='yes';
    if strcmp(type, 'TFR') && strcmp(avg, 'True')
        [stat] = ft_freqstatistics(cfg,  HC{:}, LC{:});
    elseif strcmp(type, 'ER') && strcmp(avg, 'True')
        [stat] = ft_freqstatistics(cfg, HCtlk{:}, LCtlk{:});
    end
    
else % If singletaper OPTICAT data, estimate for all selected electrodes
    
    for chnl = 1:size(electrodes)
        cfg.channel  = electrodes(chnl);
        
        if strcmp(type, 'TFR') && strcmp(avg, 'False')
            [stat(chnl)] = ft_freqstatistics(cfg,  HC{:}, LC{:});
        elseif strcmp(type, 'ER') && strcmp(avg, 'False')
            [stat(chnl)] = ft_freqstatistics(cfg, HCtlk{:}, LCtlk{:});
        elseif strcmp(type, 'TFR') && strcmp(avg, 'True')
            [stat] = ft_freqstatistics(cfg, allchan_highTFR{:}, allchan_lowTFR{:});
        elseif strcmp(type, 'ER') && strcmp(avg, 'True')
            [stat] = ft_freqstatistics(cfg, allchan_highER{:}, allchan_lowER{:});
        end
        
    end
end


%% 5.3 Figure of the results of all channels (if applicable)

figure;
for chan=1:length(stat)
    
    % Preparation
    cfg                     = [];
    cfg.channel = stat(chan).label;
    cfg.parameter           = 'stat';
    cfg.maskparameter       = 'mask';
    cfg.maskstyle           = 'outline';
    cfg.colormap = '*RdBu';
    cfg.figure = 'gcf';
    cfg.colorbartext = 't-value';
    cfg.zlim        = [-6 6];
    
    figure;
    set(gcf, 'Position', [0, 0, 650, 1200]); % Specify the figure size
    
    % Plot
    if strcmp(type, 'TFR')
        ft_singleplotTFR(cfg, stat(chan));
    elseif strcmp(type, 'ER')
        ft_singleplotER(cfg, stat(chan));
    end

    xlabel('Time [sec]');
    ylabel('Frequency [Hz]');
       ax = gca;
    ax.FontSize = 28;
    ax.TitleFontSizeMultiplier = 1.3;
    if ismember(chan, [2 3 5 6 8 9 11 12])
        ylabel('')
        set(gca, 'YTickLabel', [])
    end
    if ismember(chan, [1 2 3 4 5 6 7 8 9 10])
        set(gca, 'XTickLabel', [])
        xlabel('')
    end
    title('')
    title(stat(chan).label)
    hold on
end

% Save
saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\Cluster_All'), 'png');



%% 5.4 Singleplots for electrodes with significant effects

sign_electrodes = {'PPO2', 'Pz'}; % Determined visually

% Preparation
cfg                     = [];
cfg.parameter           = 'stat';
cfg.maskparameter       = 'mask';
cfg.maskstyle           = 'outline';
cfg.xlim                = [0 2];
cfg.layout              = layout;
cfg.colormap            = '*RdBu';
cfg.figure              = 'gcf';
cfg.colorbartext        = 't-value';


% Plot for all selected channels
for chan = 1:length(sign_electrodes) 
    
    % Select correct statistical result
    statindex = find(strcmp([stat.label], char(sign_electrodes(chan))));
    
    % Plot
    figure;
    if strcmp(type, 'TFR')
        cfg.zlim            = [-6 6];
        ft_singleplotTFR(cfg, stat(statindex));
    elseif strcmp(type, 'ER')
        ft_singleplotER(cfg, stat(statindex));
    end
    
    set(gcf, 'Position', [0, 0, 1300, 800]); % Specify the figure size
    ax = gca;
    ax.FontSize = 28;
    ax.TitleFontSizeMultiplier = 1.3;
    ylabel('Frequency [Hz]');
    xlabel('Time [sec]');
    title('')
   
    % Save
    saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\Cluster_', char(sign_electrodes(chan))), 'png');
    
    hold on

end


%% 6. Statistical analysis: timelocked data (averaged over stimulus presentation period)

%% 6.1 Estimation

% Cluster-based permutation test of the effect of contrast on power
for chan = 1:length(sign_electrodes)
    cfg                     = [];
    cfg.spmversion          = 'spm12';
    cfg.method              = 'montecarlo'; % Monte Carlo method
    cfg.statistic           = 'ft_statfun_depsamplesT'; % Dependent samples t-test
    cfg.correctm            = 'cluster';
    cfg.clusteralpha        = 0.05;
    cfg.clusterstatistic    = 'maxsum'; % Sum of t-values is the test value for a cluster
    cfg.channel             = char(sign_electrodes(chan)); % Perform the test for the selected channel
    cfg.tail                = 0;
    cfg.clustertail         = 0;
    cfg.alpha               = 0.05;
    cfg.numrandomization    = 1000;
    cfg.frequency           = [30 120];  % For the gamma frequency range
    cfg.neighbours          = [];
    
    clear design
    subj = numel(subjects);
    design = zeros(2,2*subj);
    for i = 1:subj
        design(1,i)         = i;
    end
    for i = 1:subj
        design(1, subj + i)    = i;
    end
    design(2, 1:subj)        = 1;
    design(2, subj + 1:2*subj) = 2;
    
    cfg.design              = design;
    cfg.uvar                = 1;
    cfg.ivar                = 2;
    
    [stattlk(chan)] = ft_freqstatistics(cfg, HCtlk{:}, LCtlk{:}); % over whole time
end


%% 6.2 Plot

color1 = [141/255, 102/255, 139/255];
color2 = [205/255 192/255 176/255];

for chan = 1:length(sign_electrodes)
    
    tvalues = stattlk(chan).stat;
    freq = stattlk(chan).freq;
    mask = double(stattlk(chan).mask);
    
    if ismember(1, mask) % Figure with box if there are significant results
        
        ylimlow = -4.5;
        ylimup = 4.5;
        xlim = [min(freq) max(freq)];
        ylim = [ylimlow ylimup];
        
        figure;
        box_width = 0.5; % Width of the box
        box_start = find(mask, 1, 'first'); % Find the starting point of consecutive 1s
        box_end = find(mask, 1, 'last'); % Find the ending point of consecutive 1s
        % Create a box for the frequency range with significant effects of contrast
        patch([freq(box_start)-box_width/2, freq(box_end)+box_width/2, freq(box_end)+box_width/2, freq(box_start)-box_width/2], [ylimlow, ylimlow, ylimup, ylimup], color2, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        
        set(gcf, 'Position', [0, 0, 1250, 800]); % Specify the figure size
        ax = gca;
        ax.FontSize = 28;
        
        hold on;
        plot(freq, tvalues, 'Color', color1, 'LineWidth', 3)
        xlabel('Frequency [Hz]');
        ylabel('t-value');
    else % Figure without box if there are no significant results
        figure;plot(freq, tvalues, 'Color', color1, 'LineWidth', 3)
        set(gcf, 'Position', [0, 0, 1250, 800]); % Specify the figure size
        ax = gca;
        ax.FontSize = 28;
        xlabel('Frequency [Hz]');
        ylabel('t-value');
        xlim = [min(freq) max(freq)];
    end
    
    % Extract the frequencies with significant effects of contrast 
    ind = find(stattlk(chan).negclusterslabelmat == 1);
    freq_sign_low{chan} = [stattlk(chan).freq(ind)];
    
    saveas(gcf, strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\figures\FinalForMA\EEG\', folder, '\Tvalues_', char(sign_electrodes(chan))), 'png');
end

close all




%% 7. Statistical Analysis: GPP and GPF

%% 7.1 Extract GPP and GPF values

cfg                     = [];
cfg.frequency           = [30 90]; % Define frequency range

for chnnl = 1:length(electrodes)
    cfg.channel             = char(electrodes(chnnl));
    
    for subj = 1:length(subjects) % Low contrast
        tmp = ft_selectdata(cfg, HCtlk{subj});
        PeaksChannel(chnnl).peakfreqHC(subj) = tmp.freq(find(tmp.powspctrm == max(tmp.powspctrm))); % Find frequency peak
        PeaksChannel(chnnl).peakpowHC(subj) = max(tmp.powspctrm); % Find power peak
    end
    
    for subj = 1:length(subjects) % High contrast
        tmp = ft_selectdata(cfg, LCtlk{subj});
        PeaksChannel(chnnl).peakfreqLC(subj) = tmp.freq(find(tmp.powspctrm == max(tmp.powspctrm))); % Find frequency peak
        PeaksChannel(chnnl).peakpowLC(subj) = max(tmp.powspctrm); % Find power peak
    end
    
    % Write table for later plotting (first two lines for low, third and fourth line for high contrast)
    t = table();
    t(1,:) = table(round((PeaksChannel(chnnl).peakfreqLC), 2));
    t(2,:) = table(round((PeaksChannel(chnnl).peakpowLC), 2));
    t(3,:) = table(round((PeaksChannel(chnnl).peakfreqHC), 2));
    t(4,:) = table(round((PeaksChannel(chnnl).peakpowHC), 2));
    
    % Define filename
    if strcmp(datatype, 'blk') && strcmp(exclusion, 'no')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\PeakMeans\SubMeans_blk';
    elseif strcmp(datatype,'automagic') && strcmp(tapertype, 'singletaper') && strcmp(exclusion, 'no')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\PeakMeans\SubMeans_NoOPTICAT_single';
    elseif strcmp(datatype,'automagic') && strcmp(tapertype, 'multitaper') && strcmp(exclusion, 'no')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\PeakMeans\SubMeans_NoOPTICAT_multi';
    elseif strcmp(datatype,'automagic_opticat') && strcmp(tapertype, 'singletaper') && strcmp(exclusion, 'no')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\PeakMeans\SubMeans_Singletaper';
    elseif strcmp(datatype,'automagic_opticat') && strcmp(tapertype, 'multitaper') && strcmp(exclusion, 'no')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\PeakMeans\SubMeans_Multitaper';
    elseif strcmp(exclusion, 'blk')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\PeakMeans\SubMeans_blk';
    elseif strcmp(exclusion, 'blkET')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\PeakMeans\SubMeans_blkET';
    elseif strcmp(exclusion, 'deviation')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\PeakMeans\SubMeans_deviation';
    end
    
    % Save
    writetable(t, strcat(filename, '_', char(electrodes(chnnl)),'.xlsx'),'Sheet',1,'Range','D1')
end



%% 7.2 T-tests

for chnnl = 1:length(electrodes)
    
    %% 7.1 T-test for frequency
    [H_freq,P_freq,CI_freq,STATS_freq] = ttest(PeaksChannel(chnnl).peakfreqHC, PeaksChannel(chnnl).peakfreqLC);
    
    %% 7.2 T-test for power
    [H_pow,P_pow,CI_pow,STATS_pow] = ttest(PeaksChannel(chnnl).peakpowHC, PeaksChannel(chnnl).peakpowLC);
    
    %% 7.3 Write table
    % First line for frequency, second line for power
    t = table();
    t(1,1) = table(round((H_freq), 3));
    t(1,2) = table(round((P_freq), 3));
    t(1,3) = table(CI_freq(1));
    t(1,4) = table(CI_freq(2));
    t(1,5) = table(round((STATS_freq.tstat), 3));
    t(1,6) = table(round((STATS_freq.sd), 3));
    t(2,1) = table(round((H_pow), 3));
    t(2,2) = table(round((P_pow), 3));
    t(2,3) = table(CI_pow(1));
    t(2,4) = table(CI_pow(2));
    t(2,5) = table(round((STATS_pow.tstat), 3));
    t(2,6) = table(round((STATS_pow.sd), 3));
    
    % Define file name
    if strcmp(datatype, 'blk') && strcmp(exclusion, 'no')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\TTestsEEG\TTestresults_blk';
    elseif strcmp(datatype,'automagic') && strcmp(tapertype, 'singletaper') && strcmp(exclusion, 'no')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\TTestsEEG\TTestresults_NoOPTICAT_single';
    elseif strcmp(datatype,'automagic') && strcmp(tapertype, 'multitaper') && strcmp(exclusion, 'no')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\TTestsEEG\TTestresults_NoOPTICAT_multi';
    elseif strcmp(datatype,'automagic_opticat') && strcmp(tapertype, 'singletaper') && strcmp(exclusion, 'no')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\TTestsEEG\TTestresults_Singletaper';
    elseif strcmp(datatype,'automagic_opticat') && strcmp(tapertype, 'multitaper') && strcmp(exclusion, 'no')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\TTestsEEG\TTestresults_Multitaper';
    elseif strcmp(exclusion, 'blk')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\TTestsEEG\TTestresults_blk';
    elseif strcmp(exclusion, 'blkET')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\TTestsEEG\TTestresults_blkET';
    elseif strcmp(exclusion, 'deviation')
        filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Input\TTestsEEG\TTestresults_deviation';
    end
    
    % Save
    writetable(t,strcat(filename, '_', char(electrodes(chnnl)),'.xlsx'),'Sheet',1,'Range','D1')
end






