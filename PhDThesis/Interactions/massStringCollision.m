close all;
clear all;
drawThings = false;
drawEnergy = false;
drawSubplots = false;
drawSpeed = 3;
drawStart = 1;
%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = 250;   % Length of the simulation (1 second) [samples]             

% Material properties and geometry
L = 1;              % Length [m]
r = 5e-4;           % Radius [m]
A = pi * r^2;       % Cross-sectional area [m^2] (circular cross-section)
rho = 7850;         % Material density [kg / m^3]
E = 2e11;           % Young's modulus [Pa]
I = pi * r^4 / 4;   % Area moment of inertia [m^4]
T = 554;            % Tension [N]

% Scheme coefficients
c = sqrt(T / (rho * A));            % Wave speed [m/s]
kappa = sqrt(E * I / (rho * A));    % Stiffness coefficient [m^2/s]

% Grid spacing and number of intervals
h = sqrt(1/2 * (c^2*k^2 ...
    + sqrt((c^2*k^2)^2 + 16*kappa^2*k^2)));
N = floor(L/h);     % Number of intervals between grid points
h = L / N;          % Recalculation of grid spacing based on integer N

% Update coefficients
lambdaSq = c^2 * k^2 / h^2;
muSq = kappa^2 * k^2 / h^4;
lambdaSq + 4 * muSq

% Boundary conditions ([c]lamped, [s]imply supported or [f]ree)
bc = "s";            

% Change N to the usable range
Norig = N;
if bc == "c"
    N = N - 4;
elseif bc == "s"
    N = N - 2;
elseif bc == "f"
    N = N;
end
    
%% Initialise state vectors (one more grid point than the number of intervals)
uNext = zeros(N+1, 1); 
u = zeros(N+1, 1);
% u(3:10) = hann(8);
%% Initialise scheme matrices
Id  = eye(N+1);         % identity matrix

Dxx = toeplitz([-2, 1, zeros(1, N-1)]) / h^2;
Dxxxx = Dxx * Dxx;

if bc == "c"
    Dxxxx(1, 1) = 6 / h^4;
    Dxxxx(end,end) = 6 / h^4;
elseif bc == "f"
    Dxx(1, 2) = 2 / h^2;
    Dxx(end, end-1) = 2 / h^2;
    
    Dxxxx(2, 1:4) = [-2, 5, -4, 1] / h^4;
    Dxxxx(1,1:3) = [2, -4, 2] / h^4;
    Dxxxx(end-1, end-3:end) = [1, -4, 5, -2] / h^4;
    Dxxxx(end,end-2:end) = [2, -4, 2] / h^4;
end

B = 2 * Id + c^2 * k^2 * Dxx - kappa^2 * k^2 * Dxxxx;
C = -1;
%% Mass
K = 1000000;
M = 0.01;
w = -0.1;
wPrev = -0.2;

%% Collision potential
Kc = 1e6;
alfC = 1.3;
psiPrev = 0;

%% Initial conditions (raised cosine)
ratio = 0.3;
massLoc = floor(ratio * N);       % Center location
halfWidth = round(N/20);    % Half-width of raised cosine
width = 2 * halfWidth;      % Full width
rcX = 0:width;              % x-locations for raised cosine

rc = 0.5 - 0.5 * cos(2 * pi * rcX / width); % raised cosine
% u(loc-halfWidth : loc+halfWidth) = rc;

% Set initial velocity to zero
uPrev = u;

% Output location
outLoc = round(0.3 * N);

%% Energy stuff
Dxp = sparse(1:Norig+1, 1:Norig+1, -ones(1, Norig+1), Norig+1, Norig+1) + ...
        sparse(1:Norig, 2:Norig+1, ones(1, Norig), Norig+1, Norig+1);
Dxp = Dxp / h;
Ip = zeros(1, N+1);
massLoc = floor(2 * N / 3);
startLoc = massLoc - 4;
Ip(startLoc:startLoc+8) = hann(9);
Ip = Ip / sum(Ip);
Jp = 1/h * Ip';
subp = 1;
if drawThings
    figure('Position', [173 578 827 220])
end

%% Simulation loop
for n = 1:lengthSound
    
    
    % interpolation and spreading
%     massLoc = 2 * N / 3 + N / 5 * sin(2 * pi * 10 * n / lengthSound);
%     Ip = zeros(1, N+1);
%     massLoc = floor(2 * N / 3);
%     Ip = rand(1, N+1);
%     Ip = Ip / sum(Ip);
%     alfIp = massLoc - floor(massLoc);
%     Ip(floor(massLoc)) = (1-alfIp);
%     Ip(floor(massLoc) + 1) = alfIp;
%     Jp = 1/h * Ip';


    eta = w - Ip * u;
    etaPrev = wPrev - Ip * uPrev;
    if psiPrev >= 0 
        kappaG = 1;
    else
        kappaG = -1;
    end
    
    %% uNext and wNext without the collision term
    uNext = B * u + C * uPrev;
    wNext = 2 * w - wPrev - k^2 * K/M * w;
    uStar = Ip * uNext;
    wStar = wNext;
    if eta < 0
        etaStar = wStar - uStar;
        if etaStar - etaPrev == 0
            g = 0;
        else
            % Update equation 
            g = -2 * psiPrev / (etaStar - etaPrev);
        end
    else
        g = kappaG * sqrt(Kc * (alfC + 1) / 2) * eta^((alfC - 1) / 2);
    end    
    Jterm = Ip * Jp * k^2 / (rho * A);
   
    Amat = [(1 + Jterm * g^2 / 4), -Jterm * g^2/4;
         -g^2 * k^2 / (4*M), (1 + g^2 * k^2 / (4*M))];
    v = [uStar + Jterm * (- g^2 / 4 * etaPrev + psiPrev * g);
         wStar - k^2 / M * (- g^2 / 4 * etaPrev + psiPrev * g)];
    
    solut = Amat \ v;
    etaNext = solut(2) - solut(1);
    uNext = uNext + k^2 / (rho * A) * Jp * (g^2 / 4 * (etaNext - etaPrev) + psiPrev * g);
    wNext = wNext - k^2 / M * (g^2 / 4 * (etaNext - etaPrev) + psiPrev * g);
    
    % check
%     etaNext - (wNext - Ip * uNext)
    psi = psiPrev + g / 2 * (etaNext - etaPrev);
%     psi=0;
    % Retrieve output
    out(n) = u(3);
    
%     hold off;
%     plot(u)
%     hold on
%     scatter(massLoc, w)
        
    % energy in the string
    kinEnergyU(n) = rho * A * h / 2 * sum((1/k * (u - uPrev)).^2);
    potEnergyU(n) = T / 2 * h * sum((Dxp * [0; u; 0]) .* (Dxp * [0; uPrev; 0])) ...
        + E * I * h / 2 * sum((1/h^2 * ([u(2:end); 0] - 2 * u + [0; u(1:end-1)])) ...
         .* (1/h^2 * ([uPrev(2:end); 0] - 2 * uPrev + [0; uPrev(1:end-1)])));
     % damping
    % total energy 
    totEnergyU(n) = kinEnergyU(n) + potEnergyU(n);

    kinEnergyW(n) = M / 2 * (1/k * (w - wPrev))^2;
    potEnergyW(n) = K / 2 * w * wPrev;
    totEnergyW(n) = kinEnergyW(n) + potEnergyW(n);

    connEnergy(n) = psiPrev^2 / 2; % should be col energy but for plotEnergy.m it's handier

    totEnergy(n) = totEnergyU(n) + totEnergyW(n) + connEnergy(n);
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        if drawSubplots
            subplot(1, 3, subp)
        end
        hold off;
        %...string
        plot(uNext, 'r', 'Linewidth', 2);    
        hold on;
        %...mass
        scatter(massLoc, wNext, 400, 'b', '.');    
        %...barrier
        offsetPlot = 2e-5;
        widthPlot = 5e-5;
        legend('$u_l^n$', '$w^n$', ...
        'interpreter', 'latex', 'Location', 'northwest')
        xticks([0, N])
        xticklabels(["$0$","$N$"])

        ylim([-0.5,1.5])
        yLim = ylim;
        xLab = xlabel("$l$", 'interpreter', 'latex');
        xLab.Position(2) = yLim(1) -0.075 * (yLim(2) - yLim(1));
        title("$n = " + n + "$", 'interpreter', 'latex')
        set(gca, 'Linewidth', 1.5, 'Fontsize', 16, ...
            'Position', [0.0508 + (subp-1) * 0.3250 0.1455 0.2669 0.7473], ...
            'TickLabelInterpreter', 'latex')
        if subp == 1
            yLab = ylabel("Displacement (m)", 'Fontname', 'times', 'Fontsize', 16)
            xLim = xlim;
            yLab.Position(1) = xLim(1) - 0.12 * (xLim(2) - xLim(1));
        end
%         subplot(234)
%         hold off;
%         plot(energy1u(1:n), 'b')
%         hold on;
%         plot(energy2u(1:n), 'r')
%         plot(energy3u(1:n), 'g')
%         plot(totEnergyu(1:n), 'k')
%         subplot(235)
%         hold off
%         plot(colEnergy1u(1:n), 'r')
%         hold on;
%         plot(colEnergy2u(1:n), 'b')
% 
%         subplot(236)
%         plot(totEnergyu(3:n) / totEnergyu(3) - 1);
        
        drawnow
        if drawSubplots
            subp = subp + 1;
            if subp > 3
                return;
            end
        end
    end

    % Update system states
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;
    psiPrev = psi;
    
end
plotEnergyConnected;