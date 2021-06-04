close all;
clear all;

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs;   % Length of the simulation (1 second) [samples]             

% Material properties and geometry
L = 1;              % Length [m]
r = 5e-4;           % Radius [m]
A = pi * r^2;       % Cross-sectional area [m^2] (circular cross-section)
rho = 7850;         % Material density [kg / m^3]
E = 2e11;           % Young's modulus [Pa]
I = pi * r^4 / 4;   % Area moment of inertia [m^4]
T = 554;            % Tension [N]

% Damping coefficients
sig0 = 1;         % Frequency-independent damping [s^{-1}]
sig1 = 0.000;       % Frequency-dependent damping [m^2/s]

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

% Boundary conditions ([c]lamped, [s]imply supported or [f]ree)
bc = "f";            

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

uNextVec = zeros(Norig+1, 1); 
uVec = zeros(Norig+1, 1);
range = 3:Norig-1;

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
C = (-1 + sig0 * k) * Id - 2 * sig1 * k * Dxx;

%% Initial conditions (raised cosine)
ratio = 0.3;
loc = floor(ratio * N);       % Center location
halfWidth = round(N/20);    % Half-width of raised cosine
width = 2 * halfWidth;      % Full width
rcX = 0:width;              % x-locations for raised cosine

rc = 0.5 - 0.5 * cos(2 * pi * rcX / width); % raised cosine
u(loc-halfWidth : loc+halfWidth) = rc; % initialise current state  
if bc == "s"
    vecLoc = loc + 1;
elseif bc == "c"
    vecLoc = loc + 2;
else
    vecLoc = loc;

end
u(loc-halfWidth : loc+halfWidth) = rc;
uVec(vecLoc-halfWidth : vecLoc+halfWidth) = rc;

% Set initial velocity to zero
uPrev = u;
uPrevVec = uVec;

% Output location
outLoc = round(0.3 * N);


%% Initialise matrices for energy
Dxp = sparse(1:Norig+1, 1:Norig+1, -ones(1, Norig+1), Norig+1, Norig+1) + ...
        sparse(1:Norig, 2:Norig+1, ones(1, Norig), Norig+1, Norig+1);
Dxp = Dxp / h;

DxxE = sparse(2:N+3, 1:N+2, ones(1, N+2), N+3, N+3) + ...
        sparse(1:N+3, 1:N+3, -2 * ones(1, N+3), N+3, N+3) + ... 
        sparse(1:N+2, 2:N+3, ones(1, N+2), N+3, N+3);
DxxE = DxxE / h^2;

tEgy = zeros(lengthSound, 1); % kinetic energy
vEgy = zeros(lengthSound, 1); % potential energy
hEgy = zeros(lengthSound, 1); % hamiltonian (total energy)
qTot = 0;
%% Simulation loop
for n = 1:lengthSound
    
    vectorForm; % don't forget to remove! (and all initialisation before this) 
    
    % Update equation 
    uNext = B / Amat * u + C / Amat * uPrev;
    
    % Retrieve output
    out(n) = u(outLoc);
    
    %% plot stuff
    subplot(311)
    if bc == "c"
        uPlot = [0; 0; u; 0; 0];
    elseif bc == "s"
        uPlot = [0; u; 0];
    elseif bc == "f"
        uPlot = u;
    end
    hold off;
    plot(uPlot)
    hold on;
    plot(uVec)
    ylim([-1, 1])
    
    subplot(312)
    plot(uPlot - uVec)
    drawnow;
    
    %% Energy
    if bc == "c"
        % energy in the system
        tEgy(n) = rho * A * h / 2 * sum((1/k * (u - uPrev)).^2);
        vEgy(n) = T / 2 * h * sum((Dxp * [0; 0; u; 0; 0]) .* (Dxp * [0; 0; uPrev; 0; 0]))...
            + E * I * h / 2 * sum((DxxE * [0; u; 0]) .* (DxxE * [0; uPrev; 0]));
        
        % damping
        q0(n) = 2 * sig0 * rho * A * h * sum((1/(2*k) * (uNext - uPrev)).^2);
        q1(n) = - 2 * sig1 * rho * A * h * sum(1/(2*k) * ([0; uNext; 0] - [0; uPrev; 0])...
            .* (1/k * (DxxE * [0; u; 0] - DxxE * [0; uPrev; 0])));
        idx = n - (1 * (n~=1));
        qTot = qTot + k * (q0(idx) + q1(idx));
        
        % total energy 
        hEgy(n) = tEgy(n) + vEgy(n) + qTot;

    elseif bc == "s"
        % energy in the system
        tEgy(n) = rho * A * h / 2 * sum((1/k * (u - uPrev)).^2);
        vEgy(n) = T / 2 * h * sum((Dxp * [0; u; 0]) .* (Dxp * [0; uPrev; 0])) ...
            + E * I * h / 2 * sum((1/h^2 * ([u(2:end); 0] - 2 * u + [0; u(1:end-1)])) ...
             .* (1/h^2 * ([uPrev(2:end); 0] - 2 * uPrev + [0; uPrev(1:end-1)])));
        
        % damping
        q0(n) = 2 * sig0 * rho * A * h * sum((1/(2*k) * (uNext - uPrev)).^2);
        q1(n) = - 2 * sig1 * rho * A * h * sum(1/(2*k) * ([0; uNext; 0] - [0; uPrev; 0])...
            .* (1/k * (DxxE * [0; u; 0] - DxxE * [0; uPrev; 0])));
        idx = n - (1 * (n~=1));
        qTot = qTot + k * (q0(idx) + q1(idx));
        
        % total energy 
        hEgy(n) = tEgy(n) + vEgy(n) + qTot;
        
    elseif bc == "f"

        scaling = ones(N+1, 1);
        scaling(1) = 0.5;
        scaling(end) = 0.5;
        
        tEgy(n) = rho * A * h / 2 * sum(scaling .* (1/k * (u - uPrev)).^2);
        vEgy(n) = T / 2 * h * sum((1/h * (u(2:end) - u(1:end-1))) .* (1/h * (uPrev(2:end) - uPrev(1:end-1)))) ...
            + E * I * h / 2 * sum((1/h^2 * (u(3:end) - 2 * u(2:end-1) + u(1:end-2))) ...
             .* (1/h^2 * (uPrev(3:end) - 2 * uPrev(2:end-1) + uPrev(1:end-2))));
         
        % damping
        q0(n) = 2 * sig0 * rho * A * h * sum(scaling .* (1/(2*k) * (uNext - uPrev)).^2);
        q1(n) =  - 2 * sig1 * rho * A * h * sum(1/(2*k) * scaling .* (uNext - uPrev)...
            .* (1/k * (Dxx * u - Dxx * uPrev)));
        idx = n - (1 * (n~=1));
        qTot = qTot + k * (q0(idx) + q1(idx));
        hEgy(n) = tEgy(n) + vEgy(n) + qTot;
        
    end
    subplot(313)
    plot(hEgy(1:n) / hEgy(1) - 1)
    drawnow;

    % Update system states
    uPrev = u;
    u = uNext;
    
    uPrevVec = uVec;
    uVec = uNextVec;
end

