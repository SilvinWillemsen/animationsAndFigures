close all;
clear all;

fs = 44100;
k = 1/fs;
lengthSound = fs*2;

figure('Position', [0, 300, 600, 300]);
drawStuff = false;
drawSpeed = 1;
plotContinuous = false;

%% Parameters bar 1
L1 = 1;
kappa1 = 10;
rho1 = 7850;
A1 = 0.0005;

% grid spacing
h1 = sqrt(2 * kappa1 * k);
N1 = floor(L1 / h1);
h1 = L1/N1;

%% Parameters bar 2
L2 = 1;
kappa2 = 15;
rho2 = 7850;
A2 = 0.0005;

% grid spacing
h2 = sqrt(2 * kappa2 * k);
N2 = floor(L2 / h2);
h2 = L2/N2;

%% Initialise (and excite) states
uNext = zeros(N1+1, 1);
u = zeros(N1+1, 1);
width = 5;

startHann = floor(N1/pi) - width;
% startHann = width + 1;
endHann = startHann + 2*width;
u(startHann:endHann) = 0.5 * hann(width*2+1);
uPrev = u;

wNext = zeros(N2+1, 1);
w = zeros(N2+1, 1);
wPrev = w;

%% Loop

% set clamped range
range1 = 3:N1-1;
range2 = 3:N2-1;
conn1 = floor(N1*4/5);
conn2 = floor(N2 * 6/7);
for n = 1:lengthSound
    
    % update equation
    uNext(range1) = 2 * u(range1) - uPrev(range1) - kappa1^2 * k^2 / h1^4 * (u(range1+2) - 4 * u(range1+1) + 6 * u(range1) - 4 * u(range1-1) + u(range1-2));
    wNext(range2) = 2 * w(range2) - wPrev(range2) - kappa2^2 * k^2 / h2^4 * (w(range2+2) - 4 * w(range2+1) + 6 * w(range2) - 4 * w(range2-1) + w(range2-2));

    dttu = -kappa1^2 / h1^4 * (u(conn1+2) - 4 * u(conn1+1) + 6 * u(conn1) - 4 * u(conn1-1) + u(conn1-2));
    dttw = -kappa2^2 / h2^4 * (w(conn2+2) - 4 * w(conn2+1) + 6 * w(conn2) - 4 * w(conn2-1) + w(conn2-2));
    %% calculate force (dttu = dttw)
    Ftest = (dttu - dttw) / (1/(h1 * rho1 * A1) + 1/(h2 * rho2 * A2));
    
    %% calculate force (uNext = wNext)
    F = (uNext(conn1) - wNext(conn2)) / (k^2/(h1 * rho1 * A1) + k^2/(h2 * rho2 * A2));
%     F - Ftest
    uNext(conn1) = uNext(conn1) - k^2/(h1*rho1*A1) * F;
    wNext(conn2) = wNext(conn2) + k^2/(h2*rho2*A2) * F;
    % boundary conditions
%     uNext(2) = 2 * u(2) - uPrev(2) - kappa1^2 * k^2 / h1^4 * (u(4) - 4 * u(3) + 5 * u(2) - 4 * u(1));
%     uNext(end-1) = 2 * u(end-1) - uPrev(end-1) - kappa1^2 * k^2 / h1^4 * (-4 * u(end) + 5 * u(end-1) - 4 * u(end-2) + u(end-3));

    
    % retrieve output
    out(n) = uNext(floor(N1/pi));
%     out(n) = uNext(end-1);
    %% plotting stuff
    if mod(n, drawSpeed) == 0 && drawStuff
        hold off;
        plot([0, 0], [u(conn1), w(conn2) + 1], 'k', 'Linewidth', 1)
        hold on;
        if plotContinuous
            plot([0:N1]-conn1+1, u, 'r', 'Linewidth', 2)
            plot([0:N2]-conn2+1, w+1, 'b', 'Linewidth', 2)
        else
            scatter([0:N1]-conn1+1, u, 200, 'r', '.');
            scatter([0:N2]-conn2+1, w+1, 200, 'b', '.');
        end
        xlim([-3, max(N1,N2)+3]+1-(max(conn1, conn2)))
        ylim([-1, 2])
        xticks([])
        yticks([])
%         xlabel('$l$', 'interpreter', 'latex')
%         ylabel('$u_l^n$', 'interpreter', 'latex')
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
        set(gcf, 'color', 'w')
        drawnow;
        pause(0.05)

    end 
    %% update states
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;
end

plot ([0:lengthSound-1]/fs, out)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel("$u^n_{" + num2str(floor(N1/pi)) + "}$", 'interpreter', 'latex')

set(gca, 'Linewidth', 2, 'Fontsize', 16)
set(gcf, 'color', 'w')

figure 
plot(20 * log10 (abs(fft(out))))
xlim([0, 1500])
xlabel('Frequency (Hz)')
ylabel("dB")

set(gca, 'Linewidth', 2, 'Fontsize', 16, 'FontName', 'times')
set(gcf, 'color', 'w')
