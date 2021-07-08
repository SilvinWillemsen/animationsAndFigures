% close all;
clear all;

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs;   % Length of the simulation (1 second) [samples]             

rho = 7850;
A = pi * 0.0005^2;
c = 1470;            % Wave speed [m/s]
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
if bc == "D"
    loc = round(0.2 * N);       % Center location (is correct as boundaries are excluded)
elseif bc == "N"
    loc = round(0.5 * N) + 1;   % Center location + 1 (1-based matlab) (doesn't dampen the same modes due to different modal shapes)
end

halfWidth = round(Nu/10 - 1);    % Half-width of raised cosine
width = 2 * halfWidth;      % Full width
rcX = 0:width;              % x-locations for raised cosine

rc = 0.5 - 0.5 * cos(2 * pi * rcX / width); % raised cosine
u(loc-halfWidth : loc+halfWidth) = rc; % initialise current state  
uVec(loc-halfWidth : loc+halfWidth) = rc;

% Set initial velocity to zero
uPrev = u;
uPrevVec = uVec;

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
    
    out(n) = u(outLoc);

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
