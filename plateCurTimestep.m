close all;
clear all;

x = 0:6;
y = 0:6;

figure(1)
hold on;
Xnames = {'$l-2$'; '$l-1$'; '$l$'; '$l+1$'; '$l+2$'};
Ynames = {'$m-2$'; '$m-1$'; '$m$'; '$m+1$'; '$m+2$'};
for l = 1:6
    plot([x(1) x(end)],[y(l) y(l)], 'Color', [1 1 1], 'LineWidth', 1.5)
end
for m = 1:6
    plot([x(m) x(m)],[y(1) y(end)], 'Color', [1 1 1], 'LineWidth', 1.5)
end
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'xtick',[1:5],'xticklabel', Xnames, 'FontSize', 18)
set(gca,'ytick',[1:5],'yticklabel', Ynames, 'FontSize', 18) 

gca = axes('NextPlot', 'add');


coords = [1 3; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 3 5;
            4 2; 4 3; 4 4; 5 3];
labels = {"$-\mu^2$"
          "$-2\mu^2$"
          "$(8\mu^2 + \phi)$"
          "$-2\mu^2$"
          "$-\mu^2$"
          "$(8\mu^2 + \phi)$"
          "$(2 - 4\phi + 20\mu^2)$"
          "$(8\mu^2 + \phi)$"
          "$-\mu^2$"
          "$-2\mu^2$"
          "$(8\mu^2 + \phi)$"
          "$-2\mu^2$"
          "$-\mu^2$"};
for l = 1:6
    plot([x(1) x(end)],[y(l) y(l)], 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5)
end
for m = 1:6
    plot([x(m) x(m)],[y(1) y(end)], 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5)
end
for i = 1:length(coords(:,1))
    scatter(coords(i, 1), coords(i, 2), 'filled', 'k');
%     if l
%         text(coords(i, 1), coords(i, 2) - 0.25, labels(i), 'FontSize', 18, 'interpreter','latex', 'horizontalAlignment', 'center')
%     else
        text(coords(i, 1), coords(i, 2) + 0.25, labels(i), 'FontSize', 18, 'interpreter','latex', 'horizontalAlignment', 'center')
%     end
end
set(gca, 'XColor', [1,1,1], 'YColor', [1,1,1], 'xtick', [], 'ytick', []);
title('$u_{l,m}^n$', 'interpreter', 'latex', 'FontSize', 30);
xlim([0.5 5.5])
ylim([0.5 5.5])
