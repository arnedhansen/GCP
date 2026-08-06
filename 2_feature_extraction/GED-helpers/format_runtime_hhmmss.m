function runtime_str = format_runtime_hhmmss(runtime_seconds)
% Format a duration in seconds as HH:MM:SS.
if ~isfinite(runtime_seconds) || runtime_seconds < 0
    runtime_str = 'n/a';
    return;
end
runtime_seconds = round(runtime_seconds);
h = floor(runtime_seconds / 3600);
m = floor(mod(runtime_seconds, 3600) / 60);
s = mod(runtime_seconds, 60);
runtime_str = sprintf('%02d:%02d:%02d', h, m, s);
end
