% Load the data
clear
clc
close all
data = readtable('/Volumes/methlab/Students/Arne/GCP/data/features/merged_data.csv');

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

% Set colours for the conditions
colours = [0.2, 0.6, 1; 1, 0.6, 0.2]; % Example colours: blue and orange

% Create Gamma Power plot
figure;
subplot(1,2,1); % Gamma Power
h_power = rm_raincloud(data_power, colours, 0, 'ks', 0.3);
title('Gamma Power Across Conditions');
ylabel('Gamma Power');
xlabel('Condition');
hold on;

% Add connection lines for Gamma Power
for i = 1:length(cond1_power)
    plot([1, 2], [cond1_power(i), cond2_power(i)], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
end
hold off;

% Create Gamma Frequency plot
subplot(1,2,2); % Gamma Frequency
h_freq = rm_raincloud(data_freq, colours, 0, 'ks', 0.3);
title('Gamma Frequency Across Conditions');
ylabel('Gamma Frequency (Hz)');
xlabel('Condition');
hold on;

% Add connection lines for Gamma Frequency
for i = 1:length(cond1_freq)
    plot([1, 2], [cond1_freq(i), cond2_freq(i)], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
end
hold off;

% Adjust layout
set(gcf, 'Position', [100, 100, 1200, 600]); % Set figure size
sgtitle('Gamma Power and Frequency Across Conditions'); % Add a super title
