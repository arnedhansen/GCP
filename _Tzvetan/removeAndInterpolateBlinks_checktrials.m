function [x_nan, y_nan, x_interp, y_interp, blink_mask, is_valid] = ...
    removeAndInterpolateBlinks(x, y, time, sampling_rate, threshold, pad_ms)

% Convert pad duration to samples
pad_samples = round(pad_ms * sampling_rate / 1000);

% Initial blink detection
blink_mask_raw = (x < threshold) & (y < threshold);
blink_idx = find(blink_mask_raw);

blink_mask = false(size(x));  % default no-blinks

if ~isempty(blink_idx)
    if numel(blink_idx) == 1
        % Only one blink sample
        start_idx = max(blink_idx(1) - pad_samples, 1);
        end_idx   = min(blink_idx(1) + pad_samples, length(x));
        blink_mask(start_idx:end_idx) = true;
    else
        % Group contiguous blink samples
        diff_idx = diff(blink_idx);
        split_points = [0, find(diff_idx > 1), numel(blink_idx)];
        for i = 1:(length(split_points)-1)
            segment = blink_idx(split_points(i)+1 : split_points(i+1));
            start_idx = max(segment(1) - pad_samples, 1);
            end_idx   = min(segment(end) + pad_samples, length(x));
            blink_mask(start_idx:end_idx) = true;
        end
    end
end

% Apply NaNs
x_nan = x;
y_nan = y;
x_nan(blink_mask) = NaN;
y_nan(blink_mask) = NaN;

% Interpolate if possible
valid_x = ~isnan(x_nan);
valid_y = ~isnan(y_nan);

if sum(valid_x) < 2 || sum(valid_y) < 2
    % Not enough data to interpolate
    x_interp = nan(size(x));
    y_interp = nan(size(y));
    is_valid = false;
else
    x_interp = interp1(time(valid_x), x_nan(valid_x), time, 'linear', 'extrap');
    y_interp = interp1(time(valid_y), y_nan(valid_y), time, 'linear', 'extrap');
    is_valid = true;
end
end
