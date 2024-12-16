%% Gratings RGB
clc
clear
close all
load('/Volumes/methlab/Students/Arne/GCP/figures/paradigm/grating_lc.mat')
load('/Volumes/methlab/Students/Arne/GCP/figures/paradigm/grating_hc.mat')

% Define colours for the GCP project
colors = color_def('GCP'); % Get colours: Beige and Purple

% Create a figure
figure;
set(gcf, 'Position', [100, 200, 1000, 800], 'Color', 'w'); % Adjust size for each individual plot
set(gca, 'FontSize', 20)
hold on;

% Iterate through both gratings
for i = 1:2
    if i == 2
        grating = grating_lc;
        displayName = 'Low Contrast'; 
        color = colors(1, :); % Beige
        alphaValue = 0.85; 
    elseif i == 1
        grating = grating_hc;
        displayName = 'High Contrast';
        color = colors(2, :); % Purple
        alphaValue = 0.5;
    end
    
    % Ensure consistent RGB processing
    if size(grating, 3) == 1
        grating = repmat(grating, [1, 1, 3]);
    end

    % Flatten all channels for the histogram
    R_values = grating(:, :, 1);
    G_values = grating(:, :, 2);
    B_values = grating(:, :, 3);

    R_values = R_values(:);
    G_values = G_values(:);
    B_values = B_values(:);

    % Plot the histogram for intensity
    histogram(R_values, 256, 'FaceColor', color, 'EdgeColor', 'none', ...
        'Normalization', 'probability', 'FaceAlpha', alphaValue, 'DisplayName', displayName);
end

% Set y-axis to logarithmic scale
set(gca, 'YScale', 'log');

% Finalise the plot
xline(128, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
xline(256, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.25);
xlabel('Pixel Intensity [RGB 1-256]');
ylabel('Density');
title('RGB Distribution for Low and High Contrast Gratings');
legend({'Low Contrast', 'High Contrast'}); % Show legend
% grid on;
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/GCP/figures/paradigm/GCP_gratings_rgb_distribution.png')