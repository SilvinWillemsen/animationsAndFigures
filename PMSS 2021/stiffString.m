% close all;
clear all;

fs = 44100;
k = 1/fs;
lengthSound = fs;

% figure('Position', [0, 300, 600, 300]);
drawStuff = false;
drawSpeed = 1;

%% Parameters
L = 1;
c = 0;
sig0 = 1;
sig1 = 589;
kappa = 0.03;

% grid spacing
% htest = c * k;
htest = sqrt((c^2 * k^2 + 4 * sig1 * k + sqrt((c^2*k^2 + 4 * sig1 * k)^2 + 16 * kappa^2 * k^2)) / 2);
h = sqrt((2 * k * (sig1 + sqrt(sig1^2 + kappa^2))));
h-htest
N = floor(L / h);
h = L/N;

%% Initialise (and excite) states
uNext = zeros(N+1, 1);
u = zeros(N+1, 1);
width = 5;
% startHann = floor(N/pi) - width;
startHann = width + 1;
endHann = startHann + 2*width;
u(startHann:endHann) = hann(width*2+1);
uPrev = u;

%% Loop

% set clamped range
range = 3:N-1;
for n = 1:lengthSound
    
    % update equation
    uNext(range) = (2 * u(range) - uPrev(range) + c^2 * k^2 / h^2 * (u(range+1) - 2 * u(range) + u(range-1)) ...
        - kappa^2 * k^2 / h^4 * (u(range+2) - 4 * u(range+1) + 6 * u(range) - 4 * u(range-1) + u(range-2))...
        + sig0 * k * uPrev(range) ...
        + 2 * sig1 * k / h^2 * (u(range+1) - 2 * u(range) + u(range-1) - uPrev(range+1) + 2 * uPrev(range) - uPrev(range-1))) ...
        / (1 + sig0 * k);
    
    % boundary conditions
    uNext(2) = (2 * u(2) - uPrev(2) + c^2 * k^2 / h^2 * (u(3) - 2 * u(2) + u(1)) ...
        - kappa^2 * k^2 / h^4 * (u(4) - 4 * u(3) + 5 * u(2) - 4 * u(1)) ...
        + sig0 * k * uPrev(2) ...
        + 2 * sig1 * k / h^2 * (u(3) - 2 * u(2) + u(1) - uPrev(3) + 2 * uPrev(2) - uPrev(1))) ...
        / (1 + sig0 * k);
    uNext(end-1) = (2 * u(end-1) - uPrev(end-1) + c^2 * k^2 / h^2 * (2 * u(end) - 2 * u(end-1) + u(end-2)) ...
       - kappa^2 * k^2 / h^4 * (-4 * u(end) + 5 * u(end-1) - 4 * u(end-2) + u(end-3))...
       + sig0 * k * uPrev(end-1) ...
        + 2 * sig1 * k / h^2 * (u(end) - 2 * u(end-1) + u(end-2) - uPrev(end) + 2 * uPrev(end-1) - uPrev(end-2))) ...
        / (1 + sig0 * k);

    % retrieve output
    out(n) = uNext(floor(N/pi));
%     out(n) = uNext(end-2);
    %% plotting stuff
    if mod(n, drawSpeed) == 0 && drawStuff
        plot([0:N], u, 'k', 'Linewidth', 2)
        xlim([0, N])
        ylim([-1, 1])
        yticks([])
        xlabel('$l$', 'interpreter', 'latex')
        ylabel('$u_l^n$', 'interpreter', 'latex')
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
        set(gcf, 'color', 'w')
        drawnow;
        pause(0.05)

    end 
    %% update states
    uPrev = u;
    u = uNext;
end

% plot ([0:lengthSound-1]/fs, out)
% xlabel('$t$ (s)', 'interpreter', 'latex')
% ylabel("$u^n_{" + num2str(floor(N/pi)) + "}$", 'interpreter', 'latex')
% 
% set(gca, 'Linewidth', 2, 'Fontsize', 16)
% set(gcf, 'color', 'w')
% 
% figure 
hold on;
plot([0:lengthSound-1], 20 * log10 (abs(fft(out))))
xlim([0, 1500])

xlabel('Frequency (Hz)')
ylabel("dB")
ylim([-80, 80])
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'FontName', 'times')
set(gcf, 'color', 'w')
