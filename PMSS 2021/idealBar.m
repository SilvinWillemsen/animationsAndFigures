close all;
clear all;

fs = 44100;
k = 1/fs;
lengthSound = fs*2;

figure('Position', [0, 300, 600, 300]);
drawStuff = false;
drawSpeed = 1;

%% Parameters
L = 1;
kappa = 0.006;

% grid spacing
h = sqrt(2 * kappa * k);
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
    uNext(range) = 2 * u(range) - uPrev(range) - kappa^2 * k^2 / h^4 * (u(range+2) - 4 * u(range+1) + 6 * u(range) - 4 * u(range-1) + u(range-2));
    
    % boundary conditions
    uNext(2) = 2 * u(2) - uPrev(2) - kappa^2 * k^2 / h^4 * (u(4) - 4 * u(3) + 5 * u(2) - 4 * u(1));
    uNext(end-1) = 2 * u(end-1) - uPrev(end-1) - kappa^2 * k^2 / h^4 * (-4 * u(end) + 5 * u(end-1) - 4 * u(end-2) + u(end-3));

    
    % retrieve output
%     out(n) = uNext(floor(N/pi));
    out(n) = uNext(end-1);
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

plot ([0:lengthSound-1]/fs, out)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel("$u^n_{" + num2str(floor(N/pi)) + "}$", 'interpreter', 'latex')

set(gca, 'Linewidth', 2, 'Fontsize', 16)
set(gcf, 'color', 'w')

figure 
plot(20 * log10 (abs(fft(out))))
xlim([0, 1500])
xlabel('Frequency (Hz)')
ylabel("dB")

set(gca, 'Linewidth', 2, 'Fontsize', 16, 'FontName', 'times')
set(gcf, 'color', 'w')
