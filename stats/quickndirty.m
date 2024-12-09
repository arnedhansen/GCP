% Load your data into MATLAB, assuming you have converted it to a table format
close all
data = merged_table

% Split the data by contrast condition
low_contrast = data(data.Condition == 1, :);
high_contrast = data(data.Condition == 2, :);

% Calculate descriptive statistics
disp('Descriptive Statistics for Low Contrast:');
summary(low_contrast)
disp('Descriptive Statistics for High Contrast:');
summary(high_contrast)

% Plot comparisons
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
for i = 1:length(variables)
    subplot(2, 3, i);
    boxplot(data.(variables{i}), data.Condition, 'Labels', {'Low Contrast', 'High Contrast'});
    title(variables{i});
    ylabel(variables{i});
end
