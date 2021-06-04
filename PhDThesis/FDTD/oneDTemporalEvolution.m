clear all;
close all;

figure('Position', [0, 300, 600, 200])
fs = 44100;
k = 1/fs;

drawSpeed = 300;
plotDiscrete = false;
lengthSound = 300;

c = 1470;
h = c*k;
L = 1;
N = floor(L/h);
% h = L/N;

lambdaSq = (c*k/h)^2;

uNext = zeros(N+1, 1);
u = zeros(N+1, 1);

% widthU = N/3+1;
widthU = 5;
halfWidthU = floor(widthU/2);

u(floor(N/2)-halfWidthU : floor(N/2) + halfWidthU) = hann(widthU); 
uPrev = u;

uSave = zeros(N+1, lengthSound);
range = 2:N;
for n = 1:lengthSound
    uNext(range) = 2*u(range) - uPrev(range) + lambdaSq * (u(range+1) - 2*u(range) + u(range-1));

    %% drawthings
    uSave(:, n) = uNext;
    if mod(n, drawSpeed) == 0
        if plotDiscrete
            plot(0:N, u, 'k', 'Marker', '.', 'Markersize', 20);
            xlabel('$l$', 'interpreter', 'latex')
            ylabel("$u_l^n$", 'interpreter', 'latex')
            title("$n = " + num2str(n-1)+"$", 'interpreter', 'latex');
        else
            hold off;
            plot((0:N) / N, u, 'k', 'Linewidth', 2);
            xlabel('$x$ (m)', 'interpreter', 'latex')
            ylabel("$u$ (m)", 'interpreter', 'latex')

        end
        ylim([-1, 1])

    %     grid on;
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
        set(gcf, 'color', 'w')
    
        drawnow;
%         pause(0.1);
    end
    uPrev = u;
    u = uNext;

end
drawEnd = lengthSound / 30;
view([-36.66, 30.77]);
figure('Position', [1 300 600 328])
for n = 1:1:drawEnd
    scatter3((0:N)/(N+1), ones(N+1, 1) * n, uSave(:, n), 40, repmat(1-n/drawEnd, N+1, 3), 'Linewidth', 2);
%     plot3((0:N)/(N+1), -ones(N+1) * n, uSave(:, n), 'color', repmat([1-n/drawEnd], 3, 1), 'Linewidth', 2)
    hold on;
%     drawnow
end
xlabel("$x$", 'interpreter', 'latex')
ylabel("$t$", 'interpreter', 'latex')
zlabel("$u(x,t)$", 'interpreter', 'latex')

set(gca, 'Linewidth', 2, 'Fontsize', 16)