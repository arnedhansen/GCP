%% Preprocessing GCP

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.2');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/GCP/data/merged/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    clearvars -except subjects subj path
    datapath = strcat(path,subjects{subj});
    cd(datapath)

    %% Read blocks
    for block = 1:2
        try % Do not load emtpy blocks
            load(strcat(subjects{subj}, '_EEG_ET_GCP_block',num2str(block),'_merged.mat'))
            alleeg{block} = EEG;
            clear EEG
            fprintf('Subject GCP %.3s (%.3d/%.3d): Block %.1d loaded \n', subjects{subj}, subj, length(subjects), block)
        end
    end

    %% Segment data into epochs -2s before and 3.5s after stim onset and
    %  convert to Fieldtrip data structure
    epoch_window = [-2 3.5];
    for block = 1:2
        % 21 = trigger for presentation of low contrast horizontal
        % 22 = trigger for presentation of high contrast horizontal
        % 31 = trigger for presentation of low contrast vertical
        % 32 = trigger for presentation of high contrast vertical
        % 41 = trigger for presentation of low contrast concentric static
        % 42 = trigger for presentation of high contrast concentric static
        % 51 = trigger for presentation of low contrast concentric dynamic inward
        % 52 = trigger for presentation of high contrast concentric dynamic inward
        % 61 = trigger for presentation of low contrast concentric dynamic outward
        % 62 = trigger for presentation of high contrast concentric dynamic outward
        try
            % Horizontal
            EEG_horz_lc = pop_epoch(alleeg{block}, {'21'}, epoch_window);
            EEG_horz_lc = exclude_epochs(EEG_horz_lc, '16'); % exclude task condition (red fixcross) trials
            data_horz_lc{block} = eeglab2fieldtrip(EEG_horz_lc, 'raw');

            EEG_horz_hc = pop_epoch(alleeg{block}, {'22'}, epoch_window);
            EEG_horz_hc = exclude_epochs(EEG_horz_hc, '16'); % exclude task condition (red fixcross) trials
            data_horz_hc{block} = eeglab2fieldtrip(EEG_horz_hc, 'raw');

            % Vertical
            EEG_vert_lc = pop_epoch(alleeg{block}, {'31'}, epoch_window);
            EEG_vert_lc = exclude_epochs(EEG_vert_lc, '16'); % exclude task condition (red fixcross) trials
            data_vert_lc{block} = eeglab2fieldtrip(EEG_vert_lc, 'raw');

            EEG_vert_hc = pop_epoch(alleeg{block}, {'32'}, epoch_window);
            EEG_vert_hc = exclude_epochs(EEG_vert_hc, '16'); % exclude task condition (red fixcross) trials
            data_vert_hc{block} = eeglab2fieldtrip(EEG_vert_hc, 'raw');

            % Concentric Static
            EEG_concentric_static_lc = pop_epoch(alleeg{block}, {'41'}, epoch_window);
            EEG_concentric_static_lc = exclude_epochs(EEG_concentric_static_lc, '16'); % exclude task condition (red fixcross) trials
            data_concentric_static_lc{block} = eeglab2fieldtrip(EEG_concentric_static_lc, 'raw');

            EEG_concentric_static_hc = pop_epoch(alleeg{block}, {'42'}, epoch_window);
            EEG_concentric_static_hc = exclude_epochs(EEG_concentric_static_hc, '16'); % exclude task condition (red fixcross) trials
            data_concentric_static_hc{block} = eeglab2fieldtrip(EEG_concentric_static_hc, 'raw');

            % Concentric Dynamic Inward
            EEG_concentric_dynamic_inward_lc = pop_epoch(alleeg{block}, {'51'}, epoch_window);
            EEG_concentric_dynamic_inward_lc = exclude_epochs(EEG_concentric_dynamic_inward_lc, '16'); % exclude task condition (red fixcross) trials
            data_concentric_dynamic_inward_lc{block} = eeglab2fieldtrip(EEG_concentric_dynamic_inward_lc, 'raw');

            EEG_concentric_dynamic_inward_hc = pop_epoch(alleeg{block}, {'52'}, epoch_window);
            EEG_concentric_dynamic_inward_hc = exclude_epochs(EEG_concentric_dynamic_inward_hc, '16'); % exclude task condition (red fixcross) trials
            data_concentric_dynamic_inward_hc{block} = eeglab2fieldtrip(EEG_concentric_dynamic_inward_hc, 'raw');

            % Concentric Dynamic Outward
            EEG_concentric_dynamic_outward_lc = pop_epoch(alleeg{block}, {'61'}, epoch_window);
            EEG_concentric_dynamic_outward_lc = exclude_epochs(EEG_concentric_dynamic_outward_lc, '16'); % exclude task condition (red fixcross) trials
            data_concentric_dynamic_outward_lc{block} = eeglab2fieldtrip(EEG_concentric_dynamic_outward_lc, 'raw');

            EEG_concentric_dynamic_outward_hc = pop_epoch(alleeg{block}, {'62'}, epoch_window);
            EEG_concentric_dynamic_outward_hc = exclude_epochs(EEG_concentric_dynamic_outward_hc, '16'); % exclude task condition (red fixcross) trials
            data_concentric_dynamic_outward_hc{block} = eeglab2fieldtrip(EEG_concentric_dynamic_outward_hc, 'raw');
        end
    end

    %% Remove empty blocks
    data_horz_lc = data_horz_lc(~cellfun('isempty', data_horz_lc));
    data_horz_hc = data_horz_hc(~cellfun('isempty', data_horz_hc));
    data_vert_lc = data_vert_lc(~cellfun('isempty', data_vert_lc));
    data_vert_hc = data_vert_hc(~cellfun('isempty', data_vert_hc));
    data_concentric_static_lc = data_concentric_static_lc(~cellfun('isempty', data_concentric_static_lc));
    data_concentric_static_hc = data_concentric_static_hc(~cellfun('isempty', data_concentric_static_hc));
    data_concentric_dynamic_inward_lc = data_concentric_dynamic_inward_lc(~cellfun('isempty', data_concentric_dynamic_inward_lc));
    data_concentric_dynamic_inward_hc = data_concentric_dynamic_inward_hc(~cellfun('isempty', data_concentric_dynamic_inward_hc));
    data_concentric_dynamic_outward_lc = data_concentric_dynamic_outward_lc(~cellfun('isempty', data_concentric_dynamic_outward_lc));
    data_concentric_dynamic_outward_hc = data_concentric_dynamic_outward_hc(~cellfun('isempty', data_concentric_dynamic_outward_hc));

    %% Equalize labels
    update_labels(data_horz_lc);
    update_labels(data_horz_hc);
    update_labels(data_vert_lc);
    update_labels(data_vert_hc);
    update_labels(data_concentric_static_lc);
    update_labels(data_concentric_static_hc);
    update_labels(data_concentric_dynamic_inward_lc);
    update_labels(data_concentric_dynamic_inward_hc);
    update_labels(data_concentric_dynamic_outward_lc);
    update_labels(data_concentric_dynamic_outward_hc);

    %% Append data for conditions
    cfg = [];
    cfg.keepsampleinfo = 'no';
    data_horz_lc = ft_appenddata(cfg, data_horz_lc{:});
    data_horz_hc = ft_appenddata(cfg, data_horz_hc{:});
    data_vert_lc = ft_appenddata(cfg, data_vert_lc{:});
    data_vert_hc = ft_appenddata(cfg, data_vert_hc{:});
    data_concentric_static_lc = ft_appenddata(cfg, data_concentric_static_lc{:});
    data_concentric_static_hc = ft_appenddata(cfg, data_concentric_static_hc{:});
    data_concentric_dynamic_inward_lc = ft_appenddata(cfg, data_concentric_dynamic_inward_lc{:});
    data_concentric_dynamic_inward_hc = ft_appenddata(cfg, data_concentric_dynamic_inward_hc{:});
    data_concentric_dynamic_outward_lc = ft_appenddata(cfg, data_concentric_dynamic_outward_lc{:});
    data_concentric_dynamic_outward_hc = ft_appenddata(cfg, data_concentric_dynamic_outward_hc{:});

    %% Add trialinfo
    % Horizontal
    data_horz_lc.trialinfo = zeros(numel(data_horz_lc.trial), 1) + 21;
    data_horz_hc.trialinfo = zeros(numel(data_horz_hc.trial), 1) + 22;

    % Vertical
    data_vert_lc.trialinfo = zeros(numel(data_vert_lc.trial), 1) + 31;
    data_vert_hc.trialinfo = zeros(numel(data_vert_hc.trial), 1) + 32;

    % Concentric Static
    data_concentric_static_lc.trialinfo = zeros(numel(data_concentric_static_lc.trial), 1) + 41;
    data_concentric_static_hc.trialinfo = zeros(numel(data_concentric_static_hc.trial), 1) + 42;

    % Concentric Dynamic Inward
    data_concentric_dynamic_inward_lc.trialinfo = zeros(numel(data_concentric_dynamic_inward_lc.trial), 1) + 51;
    data_concentric_dynamic_inward_hc.trialinfo = zeros(numel(data_concentric_dynamic_inward_hc.trial), 1) + 52;

    % Concentric Dynamic Outward
    data_concentric_dynamic_outward_lc.trialinfo = zeros(numel(data_concentric_dynamic_outward_lc.trial), 1) + 61;
    data_concentric_dynamic_outward_hc.trialinfo = zeros(numel(data_concentric_dynamic_outward_hc.trial), 1) + 62;

    %% Get EyeTracking data
    cfg = [];
    cfg.channel = {'L-GAZE-X'  'L-GAZE-Y' 'L-AREA'};
    dataET_horz_lc = ft_selectdata(cfg, data_horz_lc);
    dataET_horz_hc = ft_selectdata(cfg, data_horz_hc);
    dataET_vert_lc = ft_selectdata(cfg, data_vert_lc);
    dataET_vert_hc = ft_selectdata(cfg, data_vert_hc);
    dataET_concentric_static_lc = ft_selectdata(cfg, data_concentric_static_lc);
    dataET_concentric_static_hc = ft_selectdata(cfg, data_concentric_static_hc);
    dataET_concentric_dynamic_inward_lc = ft_selectdata(cfg, data_concentric_dynamic_inward_lc);
    dataET_concentric_dynamic_inward_hc = ft_selectdata(cfg, data_concentric_dynamic_inward_hc);
    dataET_concentric_dynamic_outward_lc = ft_selectdata(cfg, data_concentric_dynamic_outward_lc);
    dataET_concentric_dynamic_outward_hc = ft_selectdata(cfg, data_concentric_dynamic_outward_hc);

    %% Get EEG data (excl. ET and EOG data)
    cfg = [];
    cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' , '-L-AREA'};
    dataEEG_horz_lc = ft_selectdata(cfg, data_horz_lc);
    dataEEG_horz_hc = ft_selectdata(cfg, data_horz_hc);
    dataEEG_vert_lc = ft_selectdata(cfg, data_vert_lc);
    dataEEG_vert_hc = ft_selectdata(cfg, data_vert_hc);
    dataEEG_concentric_static_lc = ft_selectdata(cfg, data_concentric_static_lc);
    dataEEG_concentric_static_hc = ft_selectdata(cfg, data_concentric_static_hc);
    dataEEG_concentric_dynamic_inward_lc = ft_selectdata(cfg, data_concentric_dynamic_inward_lc);
    dataEEG_concentric_dynamic_inward_hc = ft_selectdata(cfg, data_concentric_dynamic_inward_hc);
    dataEEG_concentric_dynamic_outward_lc = ft_selectdata(cfg, data_concentric_dynamic_outward_lc);
    dataEEG_concentric_dynamic_outward_hc = ft_selectdata(cfg, data_concentric_dynamic_outward_hc);

    %% Re-reference data to average or common reference
    cfg = [];
    cfg.reref   = 'yes';
    cfg.refchannel = 'all';
    dataEEG_horz_lc = ft_preprocessing(cfg, dataEEG_horz_lc);
    dataEEG_horz_hc = ft_preprocessing(cfg, dataEEG_horz_hc);
    dataEEG_vert_lc = ft_preprocessing(cfg, dataEEG_vert_lc);
    dataEEG_vert_hc = ft_preprocessing(cfg, dataEEG_vert_hc);
    dataEEG_concentric_static_lc = ft_preprocessing(cfg, dataEEG_concentric_static_lc);
    dataEEG_concentric_static_hc = ft_preprocessing(cfg, dataEEG_concentric_static_hc);
    dataEEG_concentric_dynamic_inward_lc = ft_preprocessing(cfg, dataEEG_concentric_dynamic_inward_lc);
    dataEEG_concentric_dynamic_inward_hc = ft_preprocessing(cfg, dataEEG_concentric_dynamic_inward_hc);
    dataEEG_concentric_dynamic_outward_lc = ft_preprocessing(cfg, dataEEG_concentric_dynamic_outward_lc);
    dataEEG_concentric_dynamic_outward_hc = ft_preprocessing(cfg, dataEEG_concentric_dynamic_outward_hc);

    %% Save data
    savepath = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save dataEEG dataEEG_horz_lc dataEEG_horz_hc dataEEG_vert_lc dataEEG_vert_hc dataEEG_concentric_static_lc ...
        dataEEG_concentric_static_hc dataEEG_concentric_dynamic_inward_lc dataEEG_concentric_dynamic_inward_hc ...
        dataEEG_concentric_dynamic_outward_lc dataEEG_concentric_dynamic_outward_hc
    savepathET = strcat('/Volumes/methlab/Students/Arne/GCP/data/features/',subjects{subj}, '/gaze/');
    mkdir(savepathET)
    cd(savepathET)
    save dataET dataET_horz_lc dataET_horz_hc dataET_vert_lc dataET_vert_hc dataET_concentric_static_lc ...
        dataET_concentric_static_hc dataET_concentric_dynamic_inward_lc dataET_concentric_dynamic_inward_hc ...
        dataET_concentric_dynamic_outward_lc dataET_concentric_dynamic_outward_hc
    clc
    if subj == length(subjects)
        disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. PREPROCESSING FINALIZED.'])
    else
        disp(['Subject GCP ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. Loading next subject...'])
    end
end

%% Function to update labels
function update_labels(data)
blocks = size(data);
for block = 1:blocks
    if isempty(data{block})
        break;
    else
        try
            for i = 1:blocks
                if ~isempty(data{i}.label)
                    data{block}.label = data{i}.label;
                    break;
                end
            end
        catch
            warning('Error occurred while processing block %d in data structure.', block);
        end
    end
end
end