close all;
clear all;

x = 0:6;
y = 0:6;
% figure(2)
figure('Renderer', 'painters', 'Position', [100 195 569 505])
hold on;
Xnames = {'$l-2$'; '$l-1$'; '$l$'; '$l+1$'; '$l+2$'};
Ynames = {'$m-2$'; '$m-1$'; '$m$'; '$m+1$'; '$m+2$'};
grayScale = 0.85
% for l = 1:6
%     plot([x(1) x(end)],[y(l) y(l)], 'Color', [grayScale, grayScale , grayScale], 'LineWidth', 1.5)
% end
% for m = 1:6
%     plot([x(m) x(m)],[y(1) y(end)], 'Color', [grayScale, grayScale, grayScale], 'LineWidth', 1.5)
% end
set(gca, 'Position', [0.0957 0.0643 0.8974 0.9297])
set(gca,'TickLabelInterpreter', 'latex', 'Linewidth', 1.5);
set(gca,'xtick',[1:5],'xticklabel', Xnames, 'FontSize', 18)
set(gca,'ytick',[1:5],'yticklabel', Ynames, 'FontSize', 18) 
xlim([0.5 5.5])
ylim([0.5 5.5])
grid on
% gca = axes('NextPlot', 'add');
% xlim([0.5 5.5])
% ylim([0.5 5.5])

timestep = "prev"; % cur or prev

if timestep == "cur"
    coords = [1 3; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 3 5;
                4 2; 4 3; 4 4; 5 3];
    labels = {"$-\mu^2$"
              "$-2\mu^2$"
              "$(8\mu^2 + S)$"
              "$-2\mu^2$"
              "$-\mu^2$"
              "$(8\mu^2 + S)$"
              "$(2 - 20\mu^2 - 4S)$"
              "$(8\mu^2 + S)$"
              "$-\mu^2$"
              "$-2\mu^2$"
              "$(8\mu^2 + S)$"
              "$-2\mu^2$"
              "$-\mu^2$"};
else
    coords = [2 3; 3 2; 3 3; 3 4; 4 3];
    labels = {"$-S$"
              "$-S$"
              "$(\sigma_0k - 1 + 4S)$"
              "$-S$"
              "$-S$"};
end
for l = 1:6
    plot([x(1) x(end)],[y(l) y(l)], 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5)
end
for m = 1:6
    plot([x(m) x(m)],[y(1) y(end)], 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5)
end

for i = 1:length(coords(:,1))
    if timestep == "cur"
        scatter(coords(i, 1), coords(i, 2), 100, 'filled', ...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[109/255, 207/255 , 246/255]);
    else
        scatter(coords(i, 1), coords(i, 2), 100, 'filled', ...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0, 84/255 , 166/255]);
    end
    text(coords(i, 1), coords(i, 2) - 0.25, labels(i), 'FontSize', 18, 'interpreter','latex', 'horizontalAlignment', 'center')
end

% set(gca, 'XColor', [1,1,1], 'YColor', [1,1,1], 'xtick', [], 'ytick', []);
% if timestep == "cur"
%     title('$u_{l,m}^n$', 'interpreter', 'latex', 'FontSize', 30);
% else
%     title('$u_{l,m}^{n-1}$', 'interpreter', 'latex', 'FontSize', 30);
% end