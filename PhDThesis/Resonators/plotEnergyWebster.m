figure('Position', [173 578 827 220])

energyRangeEnd = 101;

subp1 = subplot(1, 2, 1)
kin = plot(0:energyRangeEnd-1, kinEnergy(1:energyRangeEnd), 'b', 'Linewidth' , 1.5);
hold on;
pot = plot(0:energyRangeEnd-1, potEnergy(1:energyRangeEnd), 'r', 'Linewidth', 1.5);
str = '#00DB00';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

bound = plot(0:energyRangeEnd-1, boundaryEnergy(1:energyRangeEnd), 'color', color, 'Linewidth', 1.5);

tot = plot(0:energyRangeEnd-1, kinEnergy(1:energyRangeEnd) ...
    + potEnergy(1:energyRangeEnd) ...
    + boundaryEnergy(1:energyRangeEnd), 'k', 'Linewidth' , 1.5);

yLim = ylim;
ylim([yLim(1), yLim(2)*1.05])
xLab = xlabel("$n$", 'interpreter', 'latex');
xLab.Position(2) = yLim(1) -0.075 * (yLim(2) - yLim(1));

% legend({'$t$', '$v$', '$h$'}, 'interpreter', 'latex')
set(gca, 'Linewidth', 1.5, 'Fontsize', 16, ...
    'Position', [0.0532 0.1500 0.4115 0.7591], ...
    'TickLabelInterpreter', 'latex')

subplot(1, 2, 2)
plot([0, 1], [0, 0], 'w'); % plot first so that top and right axes appear
hold on;
scatter(0:energyRangeEnd-1, (totEnergy(1:energyRangeEnd) - totEnergy(1)) / totEnergy(1), 60, 'k', 'Marker', '.', 'Linewidth', 1.5);
set(gca, 'Linewidth', 1.5, 'Fontsize', 16,...
    'Position', [0.5532 0.1500 0.4115 0.7591], ...
    'TickLabelInterpreter', 'latex')
xLab2 = xlabel("$n$", 'interpreter', 'latex');
yLim2 = ylim;
ylim([yLim2(1) - 0.1 * (yLim2(2) - yLim2(1)), yLim2(2) + 0.1 * (yLim2(2) - yLim2(1))]);
yLim2 = ylim;

xLab2.Position(2) = yLim2(1) -0.075 * (yLim2(2) - yLim2(1));

set(gcf, 'color', 'w')