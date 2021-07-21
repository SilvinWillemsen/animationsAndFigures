figure('Position', [173 578 827 220])

energyRangeEnd = 101;
if ~exist('connEnergy', 'var')
    plotConnEnergy = false;
    connEnergy = zeros(energyRangeEnd, 1);
else
    plotConnEnergy = true;
end
subp1 = subplot(1, 2, 1)
uEnergy = plot(0:energyRangeEnd-1, totEnergyU(1:energyRangeEnd), 'r', 'Linewidth' , 1.5);
hold on;
wEnergy = plot(0:energyRangeEnd-1, totEnergyW(1:energyRangeEnd), 'b', 'Linewidth', 1.5);
if plotConnEnergy
    cEnergy = plot(0:energyRangeEnd-1, connEnergy(1:energyRangeEnd), 'color', [0, 0.85, 0], 'Linewidth', 1.5);
end
tot = plot(0:energyRangeEnd-1, totEnergyU(1:energyRangeEnd) + totEnergyW(1:energyRangeEnd) + connEnergy(1:energyRangeEnd) , 'k', 'Linewidth' , 1.5);

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