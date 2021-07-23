close all;
clear all;


drawThings = false;
drawEnergy = false;
drawSpeed = 30;
drawStart = 132300;
plotSubplots = true;

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs * 4;   % Length of the simulation (1 second) [samples]             

% Material properties and geometry
L = 1;              % Length [m]
r = 5e-4;           % Radius [m]
A = pi * r^2;       % Cross-sectional area [m^2] (circular cross-section)
rho = 7850;         % Material density [kg / m^3]
E = 2e11;           % Young's modulus [Pa]
I = pi * r^4 / 4;   % Area moment of inertia [m^4]
T = 1000;            % Tension [N]

% Damping coefficients
sig0 = 1;         % Frequency-independent damping [s^{-1}]
sig1 = 0.005;       % Frequency-dependent damping [m^2/s]

% Scheme coefficients
c = sqrt(T / (rho * A));            % Wave speed [m/s]
kappa = sqrt(E * I / (rho * A));    % Stiffness coefficient [m^2/s]

% Grid spacing and number of intervals
h = sqrt(1/2 * (c^2*k^2 + 4*sig1*k ...
    + sqrt((c^2*k^2 + 4*sig1*k)^2 + 16*kappa^2*k^2)));
N = floor(L/h);     % Number of intervals between grid points
h = L / N;          % Recalculation of grid spacing based on integer N

% Update coefficients
lambdaSq = c^2 * k^2 / h^2;
muSq = kappa^2 * k^2 / h^4;

% Change N to the usable range
Norig = N;

% simply supported boundary conditions
N = N - 2;

%% Bow parameters
a = 100;
BM = sqrt(2 * a) * exp(0.5);
fB = 1;
FB = 1/(rho * A);
vB = 0.2;

Iu = zeros(1, N+1);
% Iu(floor(N/10-3) : floor(N/10+3)) = hann(7);
% Iu = Iu / sum(Iu)
Iu(floor(N/8)) = 1;
Ju = 1/h * Iu';

%% Initialise state vectors (one more grid point than the number of intervals)
uNext = zeros(N+1, 1); 
u = zeros(N+1, 1);

%% Initialise scheme matrices
Id  = eye(N+1);         % identity matrix

Dxx = toeplitz([-2, 1, zeros(1, N-1)]) / h^2;
Dxxxx = Dxx * Dxx;

Amat = (1 + sig0 * k);
B = 2 * Id + c^2 * k^2 * Dxx - kappa^2 * k^2 * Dxxxx + 2 * sig1 * k * Dxx;
C = -(1 - sig0 * k) * Id - 2 * sig1 * k * Dxx;

%% Initial conditions (raised cosine)
ratio = 0.3;
loc = floor(ratio * N);       % Center location
halfWidth = round(N/20);    % Half-width of raised cosine
width = 2 * halfWidth;      % Full width
rcX = 0:width;              % x-locations for raised cosine

rc = 0.5 - 0.5 * cos(2 * pi * rcX / width); % raised cosine
% u(loc-halfWidth : loc+halfWidth) = 1.1 * rc;

% Set initial velocity to zero
uPrev = u;

% Output location
outLoc = round(0.3 * N);

%% Initialise matrices for energy
Dxp = sparse(1:Norig, 1:Norig, -ones(1, Norig), Norig, Norig) + ...
        sparse(1:Norig-1, 2:Norig, ones(1, Norig-1), Norig, Norig);
Dxp = Dxp / h;

kinEnergy = zeros(lengthSound, 1); % kinetic energy
potEnergy = zeros(lengthSound, 1); % potential energy
totEnergy = zeros(lengthSound, 1); % hamiltonian (total energy)
qTot = 0;
bowTot = 0;

if drawThings
    %initialise figure
%     figure('Position', [1037 552 357 222])
    figure('Position', [173 437 852 361])
    width = 0.28;
    inc = 0.33;
    start = 0.04
    height = 0.35
    fig = 1;
end

tol = 1e-7;
vRel = 0;
%% Simulation loop
for n = 1:lengthSound
        
    % Update equation     
    b = - 2/k^2 * (Iu * u - Iu * uPrev) - c^2 * Iu * Dxx * u ...
        + kappa^2 * Iu * Dxxxx * u + (2/k + 2 * sig0) * vB ...
        - 2 * sig1 / k * (Iu * Dxx * u - Iu * Dxx * uPrev);
    
    eps = 1;
    i = 0;
%     vRel = ;
    while eps > tol && i < 100
        g = (2/k + 2 * sig0) * vRel + FB * Iu * Ju ...
            * BM * vRel * exp(-a * vRel^2) + b;
        gDeriv = 2/k + 2 * sig0 + FB * Iu * Ju ...
            * BM * (1 - 2 * a * vRel^2) * exp(-a * vRel^2);
        
        vRelNext = vRel - g/gDeriv;
        eps = abs(vRelNext - vRel);
        vRel = vRelNext;
        i = i + 1;
    end
    uNext = (B * u + C * uPrev - Ju * k^2 * FB * BM * vRel * exp(-a * vRel^2)) / Amat;

    % Retrieve output
    out(n) = u(10);
    
    %% Energy
    % energy in the system
    kinEnergy(n) = rho * A * h / 2 * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / 2 * h * sum((Dxp * [0; u]) .* (Dxp * [0; uPrev])) ...
        + E * I * h / 2 * sum((Dxx * u) .* (Dxx * uPrev));

    % damping
    q0(n) = 2 * sig0 * rho * A * h * sum((1/(2*k) * (uNext - uPrev)).^2);
    q1(n) = - 2 * sig1 * rho * A * h * sum(1/(2*k) * (uNext - uPrev)...
        .* (1/k * (Dxx * u - Dxx * uPrev)));
    idx = n - (1 * (n~=1));
    qTot = qTot + k * (q0(idx) + q1(idx));

    % bow
    bowEnergy(n) = fB * sqrt(2 * a) * vRel * exp(-a * vRel^2 + 1/2) * (vRel + vB);
    bowTot = bowTot + k * bowEnergy(idx);
    bowTotSave(n) = bowTot;
    
    % total energy 
    totEnergy(n) = kinEnergy(n) + potEnergy(n) + qTot + bowTot;
    
    %% plot stuff
    if drawThings && mod(n, drawSpeed) == 4 && n > drawStart
        if ~drawEnergy
            if plotSubplots
            	subplot(2, 3, fig)
            end
            plot((0:Norig) / Norig, [0;u;0], 'k', 'Linewidth', 1.5)
            xticks([0, find(Iu ~= 0) / Norig, 1])
            xticklabels({'$0$','$x_\textrm{\fontsize{7}{7}\selectfont B}h$','$N$'})
            if plotSubplots
                ylim([-1.2e-3, 1.2e-3])
            end
            xLab = xlabel('$l$', 'interpreter', 'latex', 'Fontsize', 16);
            yLim = ylim;
            xLab.Position(2) = yLim(1) - (yLim(2) - yLim(1)) * 0.07;
            title("$n = " + (n-3*fs) + " + 3 f_\textrm{\fontsize{7}{7}\selectfont s}$", 'interpreter', 'latex', 'Fontsize', 16)
            yLab = ylabel('$u_l^n$', 'interpreter', 'latex', 'Fontsize', 16);
            yLab.Position(1) = -0.06;

%                 xlabel('$x$', 'interpreter', 'latex', 'Fontsize', 16)
%                 ylabel('$y$', 'interpreter', 'latex', 'Fontsize', 16)
%                 title(titles{fig}, 'interpreter', 'latex', 'Fontsize', 16)
%                 colormap(gray);
%             xticks([])
%             yticks([])
%             grid on
            if plotSubplots
                if fig < 4
                    set(gca, 'Position', [start+(fig-1)*inc 0.5800 width height], 'Linewidth', 2, ...
                        'ticklabelinterpreter', 'latex', 'Fontsize', 16, 'XGrid','on')
                else
                    set(gca, 'Position', [start+(fig-4)*inc 0.0800 width height], 'Linewidth', 2, ...
                        'ticklabelinterpreter', 'latex', 'Fontsize', 16, 'XGrid','on')
                end
                fig = fig + 1;

                if fig > 6
                    return;
                end
            end
%             set(gca, 'Fontsize', 16, 'tickLabelInterpreter', 'latex', ...
%                 'Position', [0.1036 0.1396 0.8750 0.7658], 'Linewidth', 2)
            drawnow;
        
        else
            plot(totEnergy(1:n) / totEnergy(1) - 1)         
            drawnow;
        end
%         drawnow;
    end
    % Update system states
    uPrev = u;
    u = uNext;
    
end

% plotEnergy;

figure('Position', [173 549 786 249])
outrange = drawStart:drawStart+1000;
plot(outrange, out(outrange), 'k', 'Linewidth', 2)

% xlim([0, 0.1])
% yLim = ylim;
% ylim([yLim * 1.1])


% ylim([-1.5e-4, 1.5e-4])
xlim([outrange(1), outrange(end)])
yLim = ylim;
xLim = xlim;
grid on
xLab2 = xlabel("$n$", 'interpreter', 'latex');
yLab2 = ylabel("$u_{10}^n$", 'interpreter', 'latex', 'Fontsize', 20);
yLab2.Position(1) = xLim(1) -0.02 * (xLim(2) - xLim(1));
xLab2.Position(2) = yLim(1) -0.15 * (yLim(2) - yLim(1));

% xticks(0:0.02:0.1)
set(gca, 'Linewidth', 1.5, 'Fontsize', 16,...
    'Position', [0.0471 0.1807 0.9223 0.7309], ...
    'TickLabelInterpreter', 'latex')
set(gcf, 'color', 'w')


% plot(iSave)
% %% Plotting
% figure('Position', [173 578 827 220])
% 
% t = 1:lengthSound;
% % plot(t, zeros(length(t), 1), '--', 'Linewidth', 2, 'color', [0.5, 0.5, 0.5])
% % hold on;
% subp1 = subplot(1, 2, 1);
% plot(t, out, 'k', 'Linewidth', 2)
% 
% xlim([0, 300])
% ylim([-1.1, 1.1])
% xLab = xlabel("$n$", 'interpreter', 'latex');
% yLab = ylabel("$u^n_3$", 'interpreter', 'latex');
% 
% xLab.Position(2) = -1.35;
% yLab.Position(1) = -15;
% % grid on;
% 
% set(gca, 'Linewidth', 1.5, 'Fontsize', 16, ...
%     'Position', [0.0532 0.2000 0.4115 0.7591], ...
%     'TickLabelInterpreter', 'latex')
% subplot(1, 2, 2)
% dbFFT = 20 * log10(abs(fft(out)));
% plot(0:lengthSound-1, dbFFT, 'k', 'Linewidth', 2);
% xlim([0, lengthSound / 2])
% ylim([-100, 100])
% 
% % grid on;
% % xticks([c/2 : c/2 : 3000])
% set(gca, 'Linewidth', 1.5, 'Fontsize', 16,...
%     'Position', [0.5532 0.2000 0.4115 0.7591], ...
%     'TickLabelInterpreter', 'latex')
% 
% xLab2 = xlabel("$f$ [Hz]", 'interpreter', 'latex');
% yLab2 = ylabel("Magnitude [dB]", 'interpreter', 'latex');
% 
% % xLab2.Position(1) = 750;
% yLim = ylim;
% xLab2.Position(2) = yLim(1) -0.125 * (yLim(2) - yLim(1));
% 
% % yLab2.Position(1) = -70;
% 
% set(gcf, 'color', 'w')
