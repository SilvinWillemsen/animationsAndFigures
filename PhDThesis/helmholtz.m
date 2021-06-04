clear all;
close all;

figure('Position', [440 515 560 283])
N = 100;
numFigs = 4;

phase = 0.75 * pi;
amp = 1;
plotAmp = 1.5;
figNum = [1, 3, 4, 2];
for i = 1:numFigs
    c = 0.5 * cos(phase) + 0.5;
    sgn = sign(sin(phase));
    length1 = round(N * c);
    length2 = round(N * (1-c));
    line1 = linspace(0, sgn * amp * sin(c * pi), length1);
    line2 = linspace(sgn* amp * sin(c * pi) * (length2-1)/length2, 0, length2);
    
    subplot(2, 2, figNum(i))
    
    plot(amp * sin([0:1/(N-1):1] * pi), '--', 'color', [0.5, 0.5, 0.5])
    hold on;
    plot(-amp * sin([0:1/(N-1):1] * pi), '--', 'color', [0.5, 0.5, 0.5])
    plot([line1, line2], 'k', 'Linewidth', 1.5)
    rectangle('Position', [10, -1, 3, 2], 'EdgeColor', [0, 0, 0, 0.2], 'FaceColor', [0, 0, 0, 0.2])
    if figNum(i) < 3
        yPos = 0.5;
    else
        yPos = 0;
    end
    
    if mod(figNum(i), 2) == 0
        xPos = 0.52;
    else
        xPos = 0;
    end
    
    set(gca, 'Position', [xPos yPos 0.48 0.5])
    axis off;
    ylim([-plotAmp*amp plotAmp*amp])
    myArrow([11.5, 11.5],[-0.8, 0.8], 1, 8, 10, [0.2, 0.2, 0.2])

    phase = phase + 2 * pi / numFigs;
end
set(gcf, 'color','none')