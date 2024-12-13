%% Comparison of Eye Event Frequency %%

%% Description %%
% This script analyses the Eye tracker data of the ten subjects. It counts the amount of eye events for
% each contrast condition (absolute sums and averages) and then compares these numbers with t-tests.
% The script's output is tables with the descriptive differences.

%% Script Structure

% 0. Preparation – eeglab, paths, subject IDs
% 1. Counting eye events for the different conditions
% 2. Averages across subjects / trials and subjects
% 3. T-tests
% 4. Create and save tables with descriptive differences



%% 0. Preparation

%% 0.1  Prepare eeglab and FieldTrip
addpath('\\psyger-stor02.d.uzh.ch\methlab\Students\Arne\MA\eeglab2022.1');
eeglab
ft_defaults;
close all hidden

%% 0.2 Path and subjects
subjects = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '11'};
path = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results\';

%% 0.3 Trigger information

% High Contrast: 21-24
% 21 = horizontal (0 deg), 22 = vertical (90  deg), 23 = 45 degrees, 24 = 115 degrees

% Low Contrast: 25-28
% 25 = horizontal (0 deg), 26 = vertical (90  deg), 27 = 45 degrees, 28 = 115 degrees



%% 1. Count saccades, blinks and fixations

for subj = 1:length(subjects)

    %% 1.1 Load data

    datapath = strcat(path, subjects{subj});
    cd(datapath)

    %% 1.2 Find merged data, sort out bad blocks, and load

    filter = strcat('\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\automagic_opticat_results\', num2str(subj), '\*ip_EEGmerged*.mat');
    MergedData = dir(filter(:));

    counter = 0;

    for block = 1:length(MergedData)

        counter = counter + 1;

        if startsWith(MergedData(block).name, 'bip') && block < 4
            block = block + 1;
        elseif startsWith(MergedData(block).name, 'bip') && block == 4
            break
        end

        load(MergedData(block).name)


        %% 1.3 Segment data into epochs 0 before and 2 after stim onset

        %% High Contrast Triggers: 21 hor, 22 vert, 23 45deg , 24 115deg
        hc = pop_epoch(EEG,{'21', '22', '23', '24'}, [0 2]); % Stimulus presentation time

        %% Score number of saccades etc.
        trl_hc(block) = hc.trials;
        sacc_hc(block) = sum(ismember({hc.event.type},'L_saccade'));
        fix_hc(block) = sum(ismember({hc.event.type},'L_fixation'));
        blink_hc(block) = sum(ismember({hc.event.type},'L_blink'));

        %% Low Contrast Triggers: 25 hor, 26 vert, 27 45deg , 28 115deg
        lc = pop_epoch(EEG,{'25', '26', '27', '28'}, [0 2]); % Stimulus presentation time

        %% Score number of saccades etc.
        trl_lc(block) = lc.trials;
        sacc_lc(block) = sum(ismember({lc.event.type},'L_saccade'));
        fix_lc(block) = sum(ismember({lc.event.type},'L_fixation'));
        blink_lc(block) = sum(ismember({lc.event.type},'L_blink'));


    end % Block loop


    %% 1.2 Calculate sums and save per subject

    % High Contrast: absolute amount
    saccades_hc(subj) = sum(sacc_hc(:));
    fixations_hc(subj) = sum(fix_hc(:));
    blinks_hc(subj) = sum(blink_hc(:));

    % High Contrast: average across trials
    saccades_hc_mean(subj) = sum(sacc_hc(:))/sum(trl_hc(:));
    fixations_hc_mean(subj) = sum(fix_hc(:))/sum(trl_hc(:));
    blinks_hc_mean(subj) = sum(blink_hc(:))/sum(trl_hc(:));


    % Low Contrast: absolute amount
    saccades_lc(subj) = sum(sacc_lc(:));
    fixations_lc(subj) = sum(fix_lc(:));
    blinks_lc(subj) = sum(blink_lc(:));

    % Low Contrast: average across trials
    saccades_lc_mean(subj) = sum(sacc_lc(:))/sum(trl_lc(:));
    fixations_lc_mean(subj) = sum(fix_lc(:))/sum(trl_lc(:));
    blinks_lc_mean(subj) = sum(blink_lc(:))/sum(trl_lc(:));


end % Subject loop



%% 2. Averages of...

% ...Absolute Amount
mean_sacc_lc = mean(saccades_lc);
mean_fix_lc = mean(fixations_lc);
mean_blink_lc = mean(blinks_lc);

mean_sacc_hc = mean(saccades_hc);
mean_fix_hc = mean(fixations_hc);
mean_blink_hc = mean(blinks_hc);

% Difference
diff_sacc = mean_sacc_hc - mean_sacc_lc;
diff_fix = mean_fix_hc - mean_fix_lc;
diff_blink = mean_blink_hc - mean_blink_lc;


% ...Averages across trials
mean_sacc_lcm = mean(saccades_lc_mean);
mean_fix_lcm = mean(fixations_lc_mean);
mean_blink_lcm = mean(blinks_lc_mean);

mean_sacc_hcm = mean(saccades_hc_mean);
mean_fix_hcm = mean(fixations_hc_mean);
mean_blink_hcm = mean(blinks_hc_mean);

% Difference
diff_sacc_mean = mean_sacc_hcm - mean_sacc_lcm;
diff_fix_mean = mean_fix_hcm - mean_fix_lcm;
diff_blink_mean = mean_blink_hcm - mean_blink_lcm;



%% 3. T-tests for averages

% Absolute Amounts
[H,P_sacc,CI,STATS] = ttest(saccades_lc, saccades_hc);
[H,P_fix,CI,STATS] = ttest(fixations_lc, fixations_hc);
[H,P_blink,CI,STATS] = ttest(blinks_lc, blinks_hc);

% Averages across trials
[H,P_sacc,CI,STATS] = ttest(saccades_lc_mean, saccades_hc_mean);
[H,P_fix,CI,STATS] = ttest(fixations_lc_mean, fixations_hc_mean);
[H,P_blink,CI,STATS] = ttest(blinks_lc_mean, blinks_hc_mean);


%% 4. Create tables with descriptive values

% Absolute Amounts
t = table(saccades_lc .', saccades_hc .', fixations_lc .', fixations_hc .', blinks_lc .', blinks_hc.');
filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Numbers_Eyeevents.xlsx';

writetable(t,filename,'Sheet',1,'Range','D1')


% Averages across trials
t = table(saccades_lc_mean .', saccades_hc_mean .', fixations_lc_mean .', fixations_hc_mean .', blinks_lc_mean .', blinks_hc_mean .');
filename = '\\psyger-stor02.d.uzh.ch\methlab\Students\Lea Baechlin\Results\Numbers_Eyeevents_Means.xlsx';

writetable(t,filename,'Sheet',1,'Range','D1')
