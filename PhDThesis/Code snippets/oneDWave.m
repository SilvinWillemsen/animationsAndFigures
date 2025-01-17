close all;
clear all;

exc = "rc"; % excitation type (rc, imp, or tri)
initVel = false; % give an initial velocity

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs;   % Length of the simulation (1 second) [samples]             

rho = 7850;
A = pi * 0.0005^2;
c = 735;            % Wave speed [m/s]
T = c^2 * rho * A; 

L = 1;              % Length [m]
h = c * k;          % Grid spacing [m] (from CFL condition)
N = floor(L/h);     % Number of intervals between grid points
h = L / N;          % Recalculation of grid spacing based on integer N

lambdaSq = c^2 * k^2 / h^2; % Courant number squared

% Boundary conditions ([D]irichlet or [N]eumann)
bc = "D";

% Prepare Dxx matrix
if bc == "D"
    Dxx = toeplitz([-2, 1, zeros(1, N-3)]);
    I = eye(N-1);
    Nu = N-2;
elseif bc == "N"
    Dxx = toeplitz([-2, 1, zeros(1, N-1)]);
    Dxx(1, 2) = 2;
    Dxx(end, end-1) = 2;
    I = eye(N+1);
    Nu = N;
end
Dxx = Dxx / h^2;
    
%% Initialise state vectors (one more grid point than the number of intervals)
uNext = zeros(Nu+1, 1); 
u = zeros(Nu+1, 1);

%% Initial conditions (raised cosine)
if exc == "rc"
    if bc == "D"
        loc = round(0.5 * N);       % Center location (is correct as boundaries are excluded)
    elseif bc == "N"
        loc = round(0.5 * N) + 1;   % Center location + 1 (1-based matlab) (doesn't dampen the same modes due to different modal shapes)
    end

    halfWidth = round(Nu/10 - 1);    % Half-width of raised cosine
    width = 2 * halfWidth;      % Full width
    rcX = 0:width-1;              % x-locations for raised cosine

    rc = 0.5 - 0.5 * cos(2 * pi * rcX / (width - 1)); % raised cosine
    u(loc-halfWidth : loc+halfWidth-1) = rc; % initialise current state  
elseif exc == "imp"
    u(loc) = 1;
elseif exc == "tri"
    l0 = floor(0.2 * N);
    idx = 1;
    eamp = 1;
    for ll = 1:N-1 % dirichlet
        if ll <= l0
            u(ll) = eamp / l0 * ll;
        else
            u(ll) = eamp / (l0 - N) * (ll - N);

        end
        idx = idx + 1;
    end
end

% Set initial velocity to zero
if initVel
    uPrev = zeros(size(u));
else
    uPrev = u;
end

% Output location and scaling of boundary points for energy
if bc == "D"
    outLoc = round(0.1 * (N-1));
    scaling = ones(N-1, 1);
elseif bc == "N"
    outLoc = round(0.1 * (N+1));
    scaling = [0.5; ones(N-1, 1); 0.5];
end
out = zeros(lengthSound, 1);

kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

fp = 1/(pi * k) * asin (c * k / 2 * sqrt(-eig(Dxx)));
noDamping = true;
Amat = I;
B = (2 * I + c^2 * k^2 * Dxx);
C = -I;
N = N - 2;
% plotModalAnalysis;
figure('Position', [440 576 357 222])

%% Simulation loop
for n = 1:lengthSound
    
    %% Update equation
    uNext = (2 * I + c^2 * k^2 * Dxx) * u - uPrev;
    
    %% Energy
    kinEnergy(n) = rho * A / 2 * h * sum(scaling .* (1/k * (u-uPrev)).^2);
    
    if bc == "D"
        potEnergy(n) = T/(2*h) * sum(([u; 0] - [0; u]) ...
            .* ([uPrev; 0] - [0; uPrev]));
    else
        potEnergy(n) = T/(2*h) * sum((u(2:end) - u(1:end-1))...
            .* (uPrev(2:end) - uPrev(1:end-1)));
    end
    totEnergy(n) = kinEnergy(n) + potEnergy(n);
    
    out(n) = u(floor(N/2));
    
    if n == 1 || n == 6 || n == 11
        plot((0:Nu+2) / (Nu+2), [0;u;0], 'k', 'Linewidth', 1.5)
        xticks([0 1])
        xticklabels({'$0$','$N$'})
        ylim([-1.1, 1.1])
        xLab = xlabel('$l$', 'interpreter', 'latex', 'Fontsize', 16);
%         ylim([-3, 3])
        xLim = xlim;
        yLim = ylim;
        xLab.Position(2) = yLim(1) - 0.05 * (yLim(2) - yLim(1));
        yLab = ylabel('$u_l^n$', 'interpreter', 'latex', 'Fontsize', 16);
        yLab.Position(1) = -0.05;

        set(gca, 'Fontsize', 16, 'tickLabelInterpreter', 'latex', ...
            'Position', [0.1036 0.1396 0.8750 0.8306], 'Linewidth', 2)

        drawnow;

    end
    % Update system states
    uPrev = u;
    u = uNext;
end

plotEnergy;

%% Plotting
figure('Position', [173 578 827 220])

t = 1:lengthSound;
% plot(t, zeros(length(t), 1), '--', 'Linewidth', 2, 'color', [0.5, 0.5, 0.5])
% hold on;
subp1 = subplot(1, 2, 1)
plot(t, out, 'k', 'Linewidth', 2)

xlim([0, 200])
ylim([-1.1, 1.1])
xLab = xlabel("$n$", 'interpreter', 'latex');
yLab = ylabel("$u^n_3$", 'interpreter', 'latex');

xLab.Position(2) = -1.35;
yLab.Position(1) = -15;
% grid on;

set(gca, 'Linewidth', 1.5, 'Fontsize', 16, ...
    'Position', [0.0532 0.2000 0.4115 0.7591], ...
    'TickLabelInterpreter', 'latex')
subplot(1, 2, 2)
dbFFT = 20 * log10(abs(fft(out)));
plot(0:lengthSound-1, dbFFT, 'k', 'Linewidth', 2);
xlim([0, lengthSound / 2])
ylim([-300, 100])

% grid on;
% xticks([c/2 : c/2 : 3000])
set(gca, 'Linewidth', 1.5, 'Fontsize', 16,...
    'Position', [0.5532 0.2000 0.4115 0.7591], ...
    'TickLabelInterpreter', 'latex')

xLab2 = xlabel("$f$ [Hz]", 'interpreter', 'latex');
yLab2 = ylabel("Magnitude [dB]", 'interpreter', 'latex');

% xLab2.Position(1) = 750;
xLab2.Position(2) = -350;
% yLab2.Position(1) = -70;

set(gcf, 'color', 'w')
