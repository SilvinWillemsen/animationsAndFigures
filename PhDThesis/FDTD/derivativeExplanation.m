close all;

figure("Position", [182 430 1001 313])
startXPlot = 0.08;
subplotWidth = 0.28;

syms x u y yDisc;
u = -(x-1).^6+1;

range = -1:0.01:-0;

pointsX = [0.25, 0.5, 0.75] * (range(end) - range(1)) + range(1);
pointsY = double(subs(u, pointsX));

% get function for true derivative
du = diff(u);
slope = subs(du, pointsX(2));
b = subs(u, pointsX(2)) - slope * pointsX(2);
y = slope * x + b;

titles = ["Forward difference", "Backward difference", "Centred difference"];
legendEntries = ["$\delta_{t+}", "$\delta_{t-}", "$\delta_{t\cdot}"];

for pl = 1:3
    subplot(1, 3, pl)
    % axis off;
    if pl == 1
        slopeDisc = (pointsY(3) - pointsY(2)) / (pointsX(3) - pointsX(2));
        bDisc = subs(u, pointsX(2)) - slopeDisc * pointsX(2);
        yDisc = slopeDisc * x + bDisc;
        title("$\delta_{t+}$", 'interpreter', 'latex');
    elseif pl == 2
        slopeDisc = (pointsY(2) - pointsY(1)) / (pointsX(2) - pointsX(1));
        bDisc = subs(u, pointsX(2)) - slopeDisc * pointsX(2);
        yDisc = slopeDisc * x + bDisc;
        title("$\delta_{t-}$", 'interpreter', 'latex');
    elseif pl == 3
        slopeDisc = (pointsY(3) - pointsY(1)) / (pointsX(3) - pointsX(1));
        bDisc = subs(u, pointsX(3)) - slopeDisc * pointsX(3);
        yDisc = slopeDisc * x + bDisc;

    end
    approx = plot(range, subs(yDisc, range), 'b', 'Linewidth', 2)
    hold on;
    trueDeriv = plot(range, subs(y, range), '-.', 'Linewidth', 2, 'color', [0, 0.7, 0])

    uPlot = plot(range, subs(u, range), 'k', 'Linewidth', 2)
%     hold on;
    scatter(pointsX, pointsY, 40, 'r', 'Linewidth', 2)
    
    grid on;
    if pl == 1
        set(gca, 'Fontsize', 20, 'Linewidth', 2, ...
            'xticklabel', ["$t-k$", "$t$", "$t+k$"], ...
            'yticklabel', ["$u(t-k)$", "$u(t)$", "$u(t+k)$"], ...
            'TickLabelInterpreter', 'latex', ...
            'Position', [startXPlot 0.11 subplotWidth 0.815]);
    else
       set(gca, 'Fontsize', 20, 'Linewidth', 2, ...
            'xticklabel', ["$t-k$", "$t$", "$t+k$"], ...
            'yticklabel', [], ...
            'TickLabelInterpreter', 'latex', ...
            'Position', [startXPlot+subplotWidth*1.1*(pl-1) 0.11 subplotWidth 0.815]); 
    end
    xticks(pointsX)
    yticks(pointsY)   
    legend([uPlot, approx, trueDeriv], {"$u(t)$",legendEntries(pl) + "u(t)$", "$\partial_tu(t)$"}, ...
        'interpreter', 'latex', ...
        'Position', [0.25+subplotWidth*(pl-1)*1.1 0.17 0.1 0.16]);
    ylim([pointsY(1) - (pointsY(2)-pointsY(1)), pointsY(3) + (pointsY(3)-pointsY(2))])
    yLim = ylim;
    curTitle = title(titles(pl), 'interpreter', 'latex');
    curTitle.Position(2) = curTitle.Position(2) - (yLim(2) - yLim(1)) * 0.01;

end
subs(du, pointsX(2))

set(gcf, 'color', 'w')