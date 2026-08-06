function format_power_change_db_axis(axh)
% Format y-tick labels on a power-change (dB) axis.
if isempty(axh) || ~isgraphics(axh, 'axes')
    axh = gca;
end
yt = get(axh, 'YTick');
if isempty(yt)
    return;
end
yt_lbl = arrayfun(@(v) sprintf('%g', v), yt, 'UniformOutput', false);
set(axh, 'YTickLabel', yt_lbl);
end
