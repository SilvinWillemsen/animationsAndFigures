clear all;
close all;
figure('Position', [144 541 856 257])

numModes = 10;
range = 0:1/100:1;
offset = 0;
for m = 1:numModes
    plot(range, -offset - 0.25 * sin(pi * m * (range-range(1))), 'color', [0, 0, 0, 0.2], 'Linewidth', 1.5);%, 'color', repmat((m-1) * 1/numModes, 1, 3))
    hold on;
    plot(range, -offset + 0.25 * sin(pi * m * (range-range(1))), 'k', 'Linewidth', 1.5);%, 'color', repmat((m-1) * 1/numModes, 1, 3))

    plot ([range(1), range(end)], [-offset, -offset], '--', 'color', [0.5, 0.5, 0.5])
    text(range(end)-0.05, 0.5 - offset, "$" + m + "$", 'interpreter', 'latex', ...
        'Fontsize', 16)
    if m == 1
        text(range(1), -offset - 0.15, "$0$", 'horizontalAlignment', 'center', 'interpreter', 'latex', 'Fontsize', 14)
        text(range(floor(length(range) / 2)), -offset - 0.35, "$x$", 'horizontalAlignment', 'center', 'interpreter', 'latex', 'Fontsize', 14)
        text(range(end), -offset - 0.15, "$L$", 'horizontalAlignment', 'center', 'interpreter', 'latex', 'Fontsize', 14)
    end
    
    range = range + 1.3;
    
    if m == numModes/2
        offset = offset + 1;
        range = 0:1/100:1;
    end
end
ylim([-offset - 0.5, 0.75])
xlim([-0.03, range(end)-1.22])
% ylim([-offset, 1.5])
axis off
set(gca, 'Position', [0 0 1 1])
set(gcf, 'color', 'w')