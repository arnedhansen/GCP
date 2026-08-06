function save_figure_png(fig_handle, out_path)
% Save a figure as PNG at 300 dpi.
if isempty(fig_handle) || ~ishandle(fig_handle)
    return;
end
if isempty(out_path)
    return;
end
out_dir = fileparts(out_path);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
drawnow;
pause(0.05);
set(fig_handle, 'PaperPositionMode', 'auto');
print(fig_handle, out_path, '-dpng', '-r300');
end
