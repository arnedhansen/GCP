function plot_rejection_info_text_columns(info_lines, info_viol)
% Draw gate metric text on a spectrum subplot (failed gates in red).
y0 = 1.42;
dy = 0.10;
criteria_font_size = 6;
if isempty(info_lines)
    return;
end
n_lines = numel(info_lines);
for li = 1:n_lines
    if info_viol(li)
        txt_col = [0.82 0.10 0.10];
    else
        txt_col = [0.10 0.10 0.10];
    end
    text(0.02, y0 - (li - 1) * dy, info_lines{li}, ...
        'Units', 'normalized', 'Clipping', 'off', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'FontSize', criteria_font_size, 'Color', txt_col, 'Interpreter', 'none');
end
end
