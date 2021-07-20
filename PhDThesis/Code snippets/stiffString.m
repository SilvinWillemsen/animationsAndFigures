% close all;
% clear all;
drawThings = true;
drawEnergy = false;
drawSpeed = 1;

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = 200;   % Length of the simulation (1 second) [samples]             
% tSave = zeros(1001,1);
% tSaveIdx = 1;
% for tens = 1000:2000

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
%     tSave(tSaveIdx) = lambdaSq + 4 * muSq;
%     tSaveIdx = tSaveIdx + 1;
% end
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

Amat = (1 + sig0 * k);
B = 2 * Id + c^2 * k^2 * Dxx - kappa^2 * k^2 * Dxxxx + 2 * sig1 * k * Dxx;
C = -(1 - sig0 * k) * Id - 2 * sig1 * k * Dxx;

Amat = (1 + sig0 * k) * speye(N+1) - sig1 * k * Dxx;
B = 2 * Id + c^2 * k^2 * Dxx - kappa^2 * k^2 * Dxxxx;
C = -(1 - sig0 * k) * Id - sig1 * k * Dxx;

%% Initial conditions (raised cosine)
ratio = 0.3;
loc = floor(ratio * N);       % Center location
halfWidth = round(N/20);    % Half-width of raised cosine
width = 2 * halfWidth;      % Full width
rcX = 0:width;              % x-locations for raised cosine

rc = 0.5 - 0.5 * cos(2 * pi * rcX / width); % raised cosine
u(loc-halfWidth : loc+halfWidth) = 1.1 * rc;

% Set initial velocity to zero
uPrev = u;

% Output location
outLoc = round(0.3 * N);

%% Initialise matrices for energy
if bc == 'c'
    Dxp = sparse(1:Norig-2, 1:Norig-2, -ones(1, Norig-2), Norig-2, Norig-2) + ...
            sparse(1:Norig-3, 2:Norig-2, ones(1, Norig-3), Norig-2, Norig-2);
    DxxE = sparse(2:N+3, 1:N+2, ones(1, N+2), N+3, N+3) + ...
        sparse(1:N+3, 1:N+3, -2 * ones(1, N+3), N+3, N+3) + ... 
        sparse(1:N+2, 2:N+3, ones(1, N+2), N+3, N+3);
    DxxE = DxxE / h^2;

elseif bc == 's'
    Dxp = sparse(1:Norig, 1:Norig, -ones(1, Norig), Norig, Norig) + ...
            sparse(1:Norig-1, 2:Norig, ones(1, Norig-1), Norig, Norig);
elseif bc == 'f'
    Dxp = sparse(1:Norig+1, 1:Norig+1, -ones(1, Norig+1), Norig+1, Norig+1) + ...
    sparse(1:Norig, 2:Norig+1, ones(1, Norig), Norig+1, Norig+1);
    Dxp(end, :) = [];
    
    DxxE = sparse(2:N+1, 1:N, ones(1, N), N+1, N+1) + ...
        sparse(1:N+1, 1:N+1, -2 * ones(1, N+1), N+1, N+1) + ... 
        sparse(1:N, 2:N+1, ones(1, N), N+1, N+1);
    DxxE(1, :) = [];
    DxxE(end, :) = [];

    DxxE = DxxE / h^2;
end
Dxp = Dxp / h;


kinEnergy = zeros(lengthSound, 1); % kinetic energy
potEnergy = zeros(lengthSound, 1); % potential energy
totEnergy = zeros(lengthSound, 1); % hamiltonian (total energy)
qTot = 0;

% plotModalAnalysis;

if drawThings
    %initialise figure
    figure('Position', [440 576 357 222])
end
%% Simulation loop
for n = 1:lengthSound
        
    % Update equation 
    uNext = Amat \ (B * u + C  * uPrev);
    
    % Retrieve output
    out(n) = u(3);
    
    %% Energy
    if bc == "c"
        % energy in the system
        kinEnergy(n) = rho * A * h / 2 * sum((1/k * (u - uPrev)).^2);
        potEnergy(n) = T / 2 * h * sum((Dxp * [0; u]) .* (Dxp * [0; uPrev]))...
                        + E * I * h / 2 * sum((DxxE * [0; u; 0]) .* (DxxE * [0; uPrev; 0]));
%             + E * I * h / 2 * sum((Dxx * u) .* (Dxx * uPrev));

        % damping
        q0(n) = 2 * sig0 * rho * A * h * sum((1/(2*k) * (uNext - uPrev)).^2);
        q1(n) = - 2 * sig1 * rho * A * h * sum(1/(2*k) * ([0; uNext; 0] - [0; uPrev; 0])...
            .* (1/k * (DxxE * [0; u; 0] - DxxE * [0; uPrev; 0])));
        idx = n - (1 * (n~=1));
        qTot = qTot + k * (q0(idx) + q1(idx));

        % total energy 
        totEnergy(n) = kinEnergy(n) + potEnergy(n) + qTot;

    elseif bc == "s"
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

        % total energy 
        totEnergy(n) = kinEnergy(n) + potEnergy(n) + qTot;

    elseif bc == "f"

        scaling = ones(N+1, 1);
        scaling(1) = 0.5;
        scaling(end) = 0.5;

        kinEnergy(n) = rho * A * h / 2 * sum(scaling .* (1/k * (u - uPrev)).^2);
        potEnergy(n) = T / 2 * h * sum((Dxp * u) .* (Dxp * uPrev)) ...
            + E * I * h / 2 * sum((DxxE * u) .* (DxxE * uPrev));

        % damping
        q0(n) = 2 * sig0 * rho * A * h * sum(scaling .* (1/(2*k) * (uNext - uPrev)).^2);
        q1(n) =  - 2 * sig1 * rho * A * h * sum(1/(2*k) * scaling .* (uNext - uPrev)...
            .* (1/k * (Dxx * u - Dxx * uPrev)));
        idx = n - (1 * (n~=1));
        qTot = qTot + k * (q0(idx) + q1(idx));
        totEnergy(n) = kinEnergy(n) + potEnergy(n) + qTot;

    end
    
    %% plot stuff
    if drawThings && mod(n, drawSpeed) == 0
%         if drawEnergy
%             subplot(311)
%         end
%         if n == 1 || n == 41 || n == 81
%             if bc == "c"
%                 uPlot = [0; 0; u; 0; 0];
%             elseif bc == "s"
%                 uPlot = [0; u; 0];
%             elseif bc == "f"
%                 uPlot = u;
%             end
%     %         hold off;
            plot((0:Norig) / Norig, [0;u;0], 'k', 'Linewidth', 1.5)
            xticks([0 1])
            xticklabels({'$0$','$L$'})
            ylim([-1.1, 1.1])
            xLab = xlabel('$x$', 'interpreter', 'latex', 'Fontsize', 16);
            xLab.Position(2) = -1.25;
            yLab = ylabel('$u$', 'interpreter', 'latex', 'Fontsize', 16);
            yLab.Position(1) = -0.05;
            set(gca, 'Fontsize', 16, 'tickLabelInterpreter', 'latex', ...
                'Position', [0.1036 0.1396 0.8750 0.8306], 'Linewidth', 2)
drawnow;
%             pause
%         end
        
        if drawEnergy
%             subplot(312)
%             plot(uPlot - uVec)
%             subplot(313)
            plot(totEnergy(1:n) / totEnergy(1) - 1)         
            drawnow;
        end
%         drawnow;
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
subp1 = subplot(1, 2, 1);
plot(t, out, 'k', 'Linewidth', 2)

xlim([0, 300])
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
ylim([-100, 100])

% grid on;
% xticks([c/2 : c/2 : 3000])
set(gca, 'Linewidth', 1.5, 'Fontsize', 16,...
    'Position', [0.5532 0.2000 0.4115 0.7591], ...
    'TickLabelInterpreter', 'latex')

xLab2 = xlabel("$f$ [Hz]", 'interpreter', 'latex');
yLab2 = ylabel("Magnitude [dB]", 'interpreter', 'latex');

% xLab2.Position(1) = 750;
yLim = ylim;
xLab2.Position(2) = yLim(1) -0.125 * (yLim(2) - yLim(1));

% yLab2.Position(1) = -70;

set(gcf, 'color', 'w')
