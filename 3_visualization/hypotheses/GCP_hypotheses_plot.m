%% GCP Hypotheses Visualization
clear; close all; clc

% Freqs
fMin = 30;
fMax = 90;
nPts = 2000;
f = linspace(fMin, fMax, nPts);

% Curves
peakFreq1 = 50; peakAmp1 = 4; 
peakFreq2 = 55; peakAmp2 = 7;
peakFreq3 = 60; peakAmp3 = 10;
peakFreq4 = 65; peakAmp4 = 5;
peakWidth = 6;
offset = 0;
curves = [peakFreq1, peakAmp1, peakWidth, offset
          peakFreq2, peakAmp2, peakWidth, offset
          peakFreq3, peakAmp3, peakWidth, offset
          peakFreq4, peakAmp4, peakWidth, offset];

% Visual style
colors = color_def('GCP');
lineW = 10;
lineWdashed = 3;

% Preallocate
P = zeros(size(curves,1), numel(f));

% Build spectra
for k = 1:size(curves,1)
    peakFreq  = curves(k,1);
    peakAmp   = curves(k,2);
    peakWidth = curves(k,3);
    offset    = curves(k,4);

    % Gaussian peak
    P(k,:) = offset + peakAmp .* exp(-0.5 .* ((f - peakFreq)./peakWidth).^2);
end

% Plot
fig = figure('Color','w', 'Position', [0 0 1512 982]);
ax = axes('Position',[0.12 0.14 0.78 0.76]); hold(ax,'on')

% Dashed vertical lines at peaks
for k = 1:size(curves,1)
    plot(ax, [curves(k,1), curves(k,1)], [-.5, curves(k,2)], '--', 'LineWidth', lineWdashed, 'Color', colors(k,:), 'HandleVisibility', 'off')
end

for k = 1:size(P,1)
    plot(f, P(k,:), 'LineWidth', lineW, 'Color', colors(k,:))
end

% Axis labels and formatting
xlabel('Frequency [Hz]')
ylabel('Power [dB]')
xlim([fMin fMax])
yMax = max(P(:))*1.025;
ylim([-.5 yMax])
set(ax, 'Box','off', ...
        'TickDir','out', ...
        'LineWidth',1.5, ...
        'FontName','Helvetica', ...
        'FontSize', 30)
legend({'  25% Contrast', '  50% Contrast', '  75% Contrast', '100% Contrast'}, "Box", "off", "FontSize", 30)

% Save
print(set(fig, 'Renderer', 'painters'), '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/hypotheses/GCP_hypotheses.png', '-dpng', '-r600')