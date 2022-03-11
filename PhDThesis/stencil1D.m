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

coords = cell(3, 1);
spatialSteps = [1, 3, 1];
numSideCoords = (spatialSteps - 1) / 2;
for i = 1:3
    xPos = (3-numSideCoords(i):3+numSideCoords(i))';
    coords{i} = [xPos, repmat(4-i, length(xPos), 1)];
end
% coords{1} = [3 3];              % next
% coords{2} = [2 2; 3 2; 4 2];    % cur
% coords{3} = [3 1];              % prev
colors = {[1, 1, 0]; ...
        [109/255, 207/255 , 246/255]; ...
        [0, 84/255 , 166/255]};
% coordsCur = [2 2; 3 2; 4 2];
% coordsPrev = [3 1];
xOffset = 0.1;
yOffset = 0.15;

for n = 1:3
    lRange = -(size(coords{n},1)-1)/2:(size(coords{n},1)-1)/2;
    if n == 1
        uText = "n+1";
    elseif n == 2
        uText = "n";
    else
        uText = "n-1";
    end
    
    for i = 1:length(coords{n}(:,1))
        scatter(coords{n}(i, 1), coords{n}(i, 2), 200, ...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',colors{n});

        lText = makeLText(lRange, i);
        text(coords{n}(i, 1) + xOffset, coords{n}(i, 2) + yOffset, ...
            "$q^{" + uText + "}_{" + lText + "}$", ...
            'interpreter', 'latex', 'Fontsize', 18);

    end
end

annotation('doublearrow',[0.685 0.825],[0.935 0.935])
annotation('doublearrow',[0.915 0.915],[0.564 0.8])
text(4.5, 3.67, "$h$", 'interpreter','latex', 'horizontalAlignment', 'center', 'Fontsize', 18)
text(5.70, 2.5, "$k$", 'interpreter','latex', 'horizontalAlignment', 'center', 'Fontsize', 18)
text(3, 0.05, "space", 'horizontalAlignment', 'center', 'Fontsize', 18)
text(-0.4, 2, "time", 'horizontalAlignment', 'center', 'Fontsize', 18, 'Rotation',90)
annotation('doublearrow',[0.16 0.905],[0.12 0.12])
annotation('arrow',[0.05 0.05],[0.2 0.905])
set(gca, 'Position', [0.1600 0.2000 0.7450 0.7250])
b = copyobj(gca, gcf);
set(b, 'Xcolor', [1 1 1], 'YColor', [1 1 1], 'XTickLabel', [], 'YTickLabel', [],'linewidth',2)

set(gcf, 'color', 'w')

function [lText] = makeLText (lRange, i)
    if lRange(i) < 0
        lText = "l" + lRange(i);
    elseif lRange(i) == 0
        lText = "l";
    else
        lText = "l+" + lRange(i);
    end
end