close all;
clear all;
x = 0:6;
y = 0:4;
figure('Renderer', 'painters', 'Position', [200 200 600 400])
hold on;

Xnames = {'$l-2$'; '$l-1$'; '$l$'; '$l+1$'; '$l+2$'};
Ynames = {'$n-1$'; '$n$'; '$n+1$'};
for l = 1:5
    plot([x(1) x(end)],[y(l) y(l)], 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5)
end
for m = 1:7
    plot([x(m) x(m)],[y(1) y(end)], 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5)

end
xlim([0.5 5.5])
ylim([0.5 3.5])


set(gca,'TickLabelInterpreter', 'latex');
set(gca,'xtick',[1:5],'xticklabel', Xnames, 'FontSize', 18)
set(gca,'ytick',[1:5],'yticklabel', Ynames, 'FontSize', 18) 
% grid on;
% view(30,12);

coordsCur = [1 2; 2 2; 3 2; 4 2; 5 2];
coordsPrev = [2 1; 3 1; 4 1];

% coordsCur = [2 2; 3 2; 4 2];
% coordsPrev = [3 1];

for i = 1:length(coordsCur(:,1))
    scatter(coordsCur(i, 1), coordsCur(i, 2), 200, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[109/255, 207/255 , 246/255]);
%     text(coords(i, 1), coordsCur(i, 2) - 0.25, labels(i), 'FontSize', 18, 'interpreter','latex', 'horizontalAlignment', 'center')
end

for i = 1:length(coordsPrev(:,1))
    scatter(coordsPrev(i, 1), coordsPrev(i, 2), 200, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0, 84/255 , 166/255]);
end
scatter(3, 3, 200,...
'MarkerEdgeColor','k',...
'MarkerFaceColor',[1, 1 ,0]);

annotation('doublearrow',[0.685 0.825],[0.93 0.93])
annotation('doublearrow',[0.905 0.905],[0.564 0.8])
text(4.5, 3.65, "$h$", 'interpreter','latex', 'horizontalAlignment', 'center', 'Fontsize', 18)
text(5.65, 2.5, "$k$", 'interpreter','latex', 'horizontalAlignment', 'center', 'Fontsize', 18)
text(3, 0.05, "space", 'horizontalAlignment', 'center', 'Fontsize', 18)
text(-0.4, 2, "time", 'horizontalAlignment', 'center', 'Fontsize', 18, 'Rotation',90)
annotation('doublearrow',[0.16 0.905],[0.12 0.12])
annotation('arrow',[0.05 0.05],[0.2 0.905])
set(gca, 'Position', [0.1600 0.2000 0.7450 0.7250])
b = copyobj(gca, gcf);
set(b, 'Xcolor', [1 1 1], 'YColor', [1 1 1], 'XTickLabel', [], 'YTickLabel', [],'linewidth',2)

% text(3, 3 - 0.25, '$1+\sigma_0k$', 'FontSize', 18, 'interpreter','latex', 'horizontalAlignment', 'center')
