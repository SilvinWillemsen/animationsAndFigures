%{
    Drawing arrows for the Trombone system figure. To be used with
    trombone.m
%}
if pInterp
    %% Arrows in the usual case for p
    color = [alphArr, alphArr, 1];
    if drawVnmh
        arrowX = [floor((0.5:0.5:3)') + 0.5, floor((0:0.5:2.5)') + 1];
        arrowY = repmat([-3, -2] * vScale, length(arrowX), 1);
        for i = 1:length(arrowX)
            myArrow(arrowX(i, :), arrowY(i,:), arrowLineWidth, arrowHeadWidth, arrowHeadLength, color);
        end
    end
    
    arrowX = [floor((1.5:0.5:3)'), floor((1:0.5:2.5)') + 0.5];
    arrowY = repmat([-2, -1] * vScale, length(arrowX), 1);
    for i = 1:length(arrowX)
        myArrow(arrowX(i, :), arrowY(i,:), arrowLineWidth, arrowHeadWidth, arrowHeadLength, color);
    end
    
    arrowX = [1.5, 2; 2.5, 2];
    arrowY = repmat([-1, -0] * vScale, length(arrowX), 1);
    for i = 1:length(arrowX)
        myArrow(arrowX(i, :), arrowY(i,:), arrowLineWidth, arrowHeadWidth, arrowHeadLength, color);
    end
end
%% Arrows in the usual case for q
% color = [1, alphArr, alphArr];
% arrowX = [xQLocs(end-1), xQLocs(end-1) - 0.5;...
%             xQLocs(end-1), xQLocs(end-1) + 0.5; ...
%             xQLocs(end-1) - 0.5, xQLocs(end-1); ...
%             xQLocs(end-1) + 0.5, xQLocs(end-1)];
% arrowY = [-2, -1; -2, -1; -1, -0; -1, -0];
%
% for i = 1:length(arrowX)
%     myArrow(arrowX(i, :), arrowY(i,:), arrowLineWidth, arrowHeadWidth, arrowHeadLength, color);
% end
if ~pInterp
    alf = alf-1;
end
%% Arrows calculating p^n (or v^{n-1/2})
color = [alphArr, alphArr, 1];
if pInterp
    arrowX = [floor(xPLocs(end-2):0.5:xPLocs(end) - 0.5)' + 1, floor(xPLocs(end-1)-0.5:0.5:xPLocs(end))' + 0.5];
    arrowY = repmat([-3, -2] * vScale, length(arrowX), 1);
else
    arrowX = [floor(xVLocs(end-2):0.5:xVLocs(end))' + 0.5, floor(xVLocs(end-1):0.5:xVLocs(end)+0.5)'];
    arrowY = repmat([-4, -3] * vScale, length(arrowX), 1);
end
% if includeStraightArrows
%     arrowX = [arrowX; xPLocs(end) + 0.5, xPLocs(end) + 0.5];
%     arrowY = [arrowY; -3, -1];
% end
%
if drawVnmh
    for i = 1:length(arrowX)
        myArrow(arrowX(i, :), arrowY(i,:), arrowLineWidth, arrowHeadWidth, arrowHeadLength, color);
    end
end

%% Arrows calculating q^n (or w^{n-1/2})
color = [1, alphArr, alphArr];
arrowX = arrowX + 1 + alf;

if drawVnmh
    for i = 1:length(arrowX)
%         if i == length(arrowX) && numQLeft == 2
%             lineStyle = '--';
%             arrowStartX = arrowX(i, 1) - abs(arrowX(i, 2) - arrowX(i, 1)) * 0.5;
%             arrowStartY = arrowY(i, 2) - abs(arrowY(i, 2) - arrowY(i, 1)) * 0.5;
%             myArrow([arrowStartX ;arrowX(i, 2)], [arrowStartY, arrowY(i, 2)], 1, arrowHeadWidth, arrowHeadLength, color, '--');
%         else
            myArrow(arrowX(i, :), arrowY(i, :), arrowLineWidth, arrowHeadWidth, arrowHeadLength, color);
%         end
    end
end

%% Arrows calculating velocities v (or p^n)
color = [alphArr, alphArr, 1];
arrowX = arrowX - 0.5 - alf;
arrowY = arrowY + 1 * vScale;
for i = 1:length(arrowX)
    myArrow(arrowX(i, :), arrowY(i,:), arrowLineWidth, arrowHeadWidth, arrowHeadLength, color);
end

%% Arrows calculating velocities w (or q^n) 
color = [1, alphArr, alphArr];
arrowX = arrowX + alf;

for i = 1:length(arrowX)
    myArrow(arrowX(i, :), arrowY(i,:), arrowLineWidth, arrowHeadWidth, arrowHeadLength, color);
end

%% Arrows calculating pressures at inner boundaries
color = [alphArr, alphArr, 1];
if pInterp
    arrowX = [xPLocs(end) + 0.5, xPLocs(end); ...
        xPLocs(end) - 0.5, xPLocs(end)];
else
    arrowX = [xVLocs(end) + 0.5, xVLocs(end); ...
        xVLocs(end) - 0.5, xVLocs(end)];
end
arrowY = arrowY(1:2, :) + vScale;

for i = 1:length(arrowX)
    myArrow(arrowX(i, :), arrowY(i,:), arrowLineWidth, arrowHeadWidth, arrowHeadLength, color);
end

color = [1, alphArr, alphArr];
arrowX = arrowX + alf;
% arrowY = [-1, 0; -1, 0] * vScale;
for i = 1:length(arrowX)
    myArrow(arrowX(i, :), arrowY(i,:), arrowLineWidth, arrowHeadWidth, arrowHeadLength, color);
end

%% Draw curved arrows for connection points
aHO = 0.02; % arrowHeadOffset
if pInterp
    arrowX = [xPLocs(end-1), xPLocs(end-1) + alf - aHO; ...
        xPLocs(end), xPLocs(end) + 1 + 0.5 * aHO; ...
        xPLocs(end) + alf, xPLocs(end) + 1 - aHO; ...
        xPLocs(end), xQLocs(1) - 1 + aHO; ...
        xQLocs(1), xQLocs(1) - 1 - 0.5 * aHO; ...
        xQLocs(2), xPLocs(end) + 1 + aHO];
else
    arrowX = [xVLocs(end-1), xVLocs(end-1) + alf - aHO; ...
        xVLocs(end), xVLocs(end) + 1 + 0.5 * aHO; ...
        xVLocs(end) + alf, xVLocs(end) + 1 - aHO; ...
        xVLocs(end), xWLocs(1) - 1 + aHO; ...
        xWLocs(1), xWLocs(1) - 1 - 0.5 * aHO; ...
        xWLocs(2), xVLocs(end) + 1 + aHO];
end
colors = [alphArrInterp, alphArrInterp, 1; ...
    alphArrInterp, alphArrInterp, 1; ...
    1, alphArrInterp, alphArrInterp; ...
    alphArrInterp, alphArrInterp, 1; ...
    1, alphArrInterp, alphArrInterp; ...
    1, alphArrInterp, alphArrInterp];
% color = [1, alphArrInterp, 1]; % purple
if pInterp
    arrowY = repmat([-2, -2] * vScale, length(arrowX), 1);
else
    arrowY = repmat([-3, -3] * vScale, length(arrowX), 1);
end
for i = 1:length(arrowX)
    if i == 1 || i == length(arrowX)
        bend = 9.5/8 * curvatureArcs;
    else
        bend = curvatureArcs;
    end
    if i < 4
        myArrow(arrowX(i, :), arrowY(i,:), arrowLineWidth * 2/3, arrowHeadWidth, arrowHeadLength, colors(i,:), '-.', bend, false);
    else
        myArrow(arrowX(i, :), arrowY(i,:), arrowLineWidth * 2/3, arrowHeadWidth, arrowHeadLength, colors(i,:), '-.', -bend, false);
    end
end

if ~pInterp
    alf = alf+1;
end
