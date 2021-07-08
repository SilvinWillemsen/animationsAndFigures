close all;
clear all;
x = 0:6;
y = 0:6;
z = 1:3;
figure('Renderer', 'painters', 'Position', [100 195 569 505])
hold on;
Xnames = {'$l-2$'; '$l-1$'; '$l$'; '$l+1$'; '$l+2$'};
Ynames = {'$m-2$'; '$m-1$'; '$m$'; '$m+1$'; '$m+2$'};
Znames = {'$n-1$'; '$n$'; '$n+1$'};
for zz = 1:3
    for l = 1:7
        plot3([x(1) x(end)],[y(l) y(l)], [z(zz) z(zz)], 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5)
    end
    for m = 1:7
        plot3([x(m) x(m)],[y(1) y(end)], [z(zz) z(zz)], 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5)

    end
end
xlim([0 6.1])
ylim([0 6.1])
zlim([0.9 3.5])

set(gca,'TickLabelInterpreter', 'latex');
set(gca,'xtick',[1:5],'xticklabel', [], 'FontSize', 18)

text([1:5], repmat(-0.5, 5, 1), repmat(0.8, 5, 1), Xnames, 'horizontalAlignment', 'center', 'Fontsize', 18, 'interpreter', 'latex')
set(gca,'ytick',[1:5],'yticklabel', Ynames, 'FontSize', 18) 
set(gca,'ztick',[1:3],'zticklabel', Znames, 'FontSize', 18) 
% grid on;
view(30,16);

coordsCur = [1 3; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 3 5;
            4 2; 4 3; 4 4; 5 3];

coordsPrev = [2 3; 3 2; 3 3; 3 4; 4 3];

for i = 1:length(coordsCur(:,1))
    scatter3(coordsCur(i, 1), coordsCur(i, 2), 2, 100, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[109/255, 207/255 , 246/255]);
%     text(coords(i, 1), coordsCur(i, 2) - 0.25, labels(i), 'FontSize', 18, 'interpreter','latex', 'horizontalAlignment', 'center')
end

for i = 1:length(coordsPrev(:,1))
    scatter3(coordsPrev(i, 1), coordsPrev(i, 2), 1, 100, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0, 84/255 , 166/255]);
end
scatter3(3, 3, 3, 100,...
'MarkerEdgeColor','k',...
'MarkerFaceColor',[1, 1 ,0]);
camproj('perspective')
text(4.38991979746008, 0.601020069211444, 3.18889040789108, '$(1+\sigma_0k)$', 'FontSize', 18, 'interpreter','latex', 'horizontalAlignment', 'center')

set(gca, 'Position', [0.1 0.05 0.85 1], 'Projection', 'perspective')