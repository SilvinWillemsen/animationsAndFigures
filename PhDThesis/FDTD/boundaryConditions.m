clear all;
close all;

figure('Position', [1 408 600 353])
fs = 441000;
k = 1/fs;

boundaryCond = "d"; %(d)irichlet or (n)eumann

drawSpeed = 300;
plotDiscrete = false;
lengthSound = 300;

c = 735;
h = c*k;
L = 1;
N = floor(L/h);
% h = L/N;

lambdaSq = (c*k/h)^2;

uNext = zeros(N+1, 1);
u = zeros(N+1, 1);

% widthU = N/3+1;
widthU = 100;
halfWidthU = floor(widthU/2);

u(floor(N/3)-halfWidthU : floor(N/3) + halfWidthU) = hann(widthU + 1); 
uPrev = u;

uSave = zeros(N+1, lengthSound);
range = 2:N;
j = 0;
numSubPlots = 4;
for n = 1:lengthSound
    uNext(range) = 2*u(range) - uPrev(range) + lambdaSq * (u(range+1) - 2*u(range) + u(range-1));
    if boundaryCond == "n"
        uNext(1) = 2*u(1) - uPrev(1) + lambdaSq * (2*u(2) - 2*u(1));
    end
    if mod(n, 52) == 0 && n>52
        j = j + 1;
        subplot(numSubPlots, 1, j);
        plot(0:N, uNext, 'k', 'Linewidth', 2)
        ylim([-1, 1])
        drawnow; 
        axis off;
        set(gca, 'Position', [0.07 0.7673 - 0.24*(j-1) 0.9250 0.22]);
        pause(0.1)
        if j == numSubPlots
            break;
        end
    end
    uPrev = u;
    u = uNext;

end
set(gcf, 'color', 'w')
annotation('arrow', [0.035, 0.035], [0.95, 0.1], 'Linewidth', 1)
text(-40, 0, "$t$", 'interpreter', 'latex', 'Fontsize', 26, 'horizontalAlignment', 'center')