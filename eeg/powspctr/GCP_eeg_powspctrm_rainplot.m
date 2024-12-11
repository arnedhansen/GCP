%% 

%% Load the data
clear
clc
close all
data = readtable('/Volumes/methlab/Students/Arne/GCP/data/features/merged_data.csv');
%summary(data)

%%

% Extract relevant columns
IDs = data.ID; % Participant IDs
conditions = data.Condition; % Experimental conditions
gamma_power = data.GammaPower; % Gamma Power
gamma_freq = data.GammaFreq; % Gamma Frequency

% Split data into two conditions
cond1_power = gamma_power(conditions == 1);
cond2_power = gamma_power(conditions == 2);
cond1_freq = gamma_freq(conditions == 1);
cond2_freq = gamma_freq(conditions == 2);

% Prepare data for rm_raincloud
data_power = {cond1_power, cond2_power};
data_freq = {cond1_freq, cond2_freq};

%% make figure
close all
figure;
set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');

% Generate the repeated measures raincloud plot
colours = [0.5 0.5 0.5; 1 1 1]; % Use MATLAB's default line colours (or adjust as needed)
rcp = rm_raincloud(data_power, colours);

% Adjust the Y-axis limits
set(gca, 'YLim', [-0.3 1.6]);

% Add title to the plot
title(['Figure M10' newline 'Repeated measures raincloud plot - some aesthetic options']);

% Define a new colour for a specific subset
new_cl = [0.2 0.2 0.2]; % Grey colour

% Change one subset to the new colour and adjust dot size
% Assuming subset indices are structured in rcp as rcp.p and rcp.s
rcp.p{2, 2}.FaceColor         = new_cl; % Change patch colour
rcp.s{2, 2}.MarkerFaceColor   = new_cl; % Change scatter face colour
% rcp.m(2, 2).MarkerEdgeColor   = 'none'; % Remove edge colour
% rcp.m(2, 2).MarkerFaceColor   = new_cl; % Change marker face colour
rcp.s{2, 2}.SizeData          = 300;    % Adjust scatter dot size