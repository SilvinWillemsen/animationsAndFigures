% close all;
% clear all;
clc;

fs = 44100*100;
k = 1/fs;
lengthSound = fs;

%% Set physical parameters
kappa = 1000;
L = 1;

%% Calculate grid spacing
h = sqrt(2 * kappa * k);

N = floor(L / h);
h = L / N;

%% Initialise state vectors
uNext = zeros(N+1, 1);
u = zeros(N+1, 1);
u(10:16) = hann(7);
uPrev = u;

range = 3:N-1;
for n = 1:lengthSound
    % update equation
    uNext(range) = 2 * u(range) - uPrev(range) ...
        - kappa^2 * k^2 / h^4 * (u(range+2) - 4*u(range+1) + 6*u(range) - 4*u(range-1) + u(range-2));
    
    % (u(4) - 4*u(3) + 6*u(2) - 4*u(1) + u(0))
    % u_{-1} = -u_1 -> u(0) = -u(2)
    % (u(4) - 4*u(3) + 5*u(2) - 4*u(1))
    
    uNext(2) = 2 * u(2) - uPrev(2) ...
        - kappa^2 * k^2 / h^4 * (u(4) - 4*u(3) + 5*u(2) - 4*u(1));
    uNext(N) = 2 * u(N) - uPrev(N) ...
        - kappa^2 * k^2 / h^4 * (-4*u(N+1) + 5*u(N) - 4*u(N-1) + u(N-2));

    
    out(n) = u (floor(N/pi));
    % plot things
%     plot(u);
%     xlim([N-10, N+1])
%     pause(0.5)
%     drawnow;
    
    % update states
    uPrev = u;
    u = uNext;
end
hold on
plot([0:lengthSound-1],abs(fft(out)));
xlim([0, 1500])


