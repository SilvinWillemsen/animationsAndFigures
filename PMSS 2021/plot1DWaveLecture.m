clear all;
close all;

figure('Position', [0, 300, 600, 200])
fs = 44100;
k = 1/fs;

plotDiscrete = false;
lengthSound = 100;

c = 1470/10;
h = c*k;
L = 1;
N = floor(L/h);
% h = L/N;

lambdaSq = (c*k/h)^2;

uNext = zeros(N+1, 1);
u = zeros(N+1, 1);
wNext = zeros(N+1, 1);
w = zeros(N+1, 1);

widthU = 101;
widthW = 11;

halfWidthU = floor(widthU/2);
halfWidthW = floor(widthW/2);

u(floor(N/2)-halfWidthU : floor(N/2) + halfWidthU) = hann(widthU); 
uPrev = u;

w(floor(N/2)-halfWidthW : floor(N/2) + halfWidthW) = hann(widthW); 
wPrev = w;

range = 2:N;
for n = 1:lengthSound
    uNext(range) = 2*u(range) - uPrev(range) + lambdaSq * (u(range+1) - 2*u(range) + u(range-1));
    wNext(range) = 2*w(range) - wPrev(range) + lambdaSq * (w(range+1) - 2*w(range) + w(range-1));

%     uNext(end)= 2*u(end) - uPrev(end) + lambdaSq * (2 * u(end-1) - 2*u(end));
%     uNext(1)= 2*u(1) - uPrev(1) + lambdaSq * (2 * u(2) - 2*u(1));

    output(n) = u(end-3);
    %% drawthings
    if mod(n,1) == 0
        if plotDiscrete
            plot(0:N, u, 'k', 'Marker', '.', 'Markersize', 20);
            xlabel('$l$', 'interpreter', 'latex')
            ylabel("$u_l^n$", 'interpreter', 'latex')
            title("$n = " + num2str(n-1)+"$", 'interpreter', 'latex');
        else
            hold off;
            plot((0:N) / N, u, 'b', 'Linewidth', 2);
            hold on;
            plot((0:N) / N, w, 'r', 'Linewidth', 1);

            xlabel('$x$', 'interpreter', 'latex')
%             legend('$u$', '$w$', 'interpreter', 'latex')
            ylabel("$u$", 'interpreter', 'latex')
%             title("$t = " + num2str((n-1)/fs)+"$", 'interpreter', 'latex');

        end
        ylim([-1, 1])

    %     grid on;
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
        set(gcf, 'color', 'w')

        drawnow;
        pause(0.1);
    end
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;

end