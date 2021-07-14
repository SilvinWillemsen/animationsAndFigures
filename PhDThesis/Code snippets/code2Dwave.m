close all;
clear all;

drawThings = true;
energyCalc = false;
drawSpeed = 1;
%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs;   % Length of the simulation (1 second) [samples]             

rho = 7850;                     % Material density [kg/m^3]
H =  0.0005;                    % Thickness [m]
T = 1000000;                    % Tension per unit length [N/m]
c = sqrt(T / (rho * H));        % Wave speed [m/s]

Lx = 1;                         % Length in x direction [m]
Ly = 2;                         % Length in y direction [m]

h = sqrt(2) * c * k;            % Grid spacing [m]
Nx = floor(Lx/h);               % Number of intervals in x direction
Ny = floor(Ly/h);               % Number of intervals in y direction
h = min(Lx/Nx, Ly/Ny);          % Recalculation of grid spacing

lambdaSq = c^2 * k^2 / h^2;     % Courant number squared
h = min(Lx/Nx, Ly/Ny);          % Recalculation of grid spacing

%% Create scheme matrices with Dirichlet boundary conditions 
Nxu = Nx - 1;
Nyu = Ny - 1;
Dxx = toeplitz([-2, 1, zeros(1, Nxu-2)]);
Dyy = toeplitz([-2, 1, zeros(1, Nyu-2)]);

% Kronecker sum
D = kron(speye(Nxu), Dyy) + kron(Dxx, speye(Nyu));
D = D / h^2;

% Total amount of grid points
Nu = Nxu * Nyu;
    
%% Initialise state vectors (one more grid point than the number of intervals)
uNext = zeros(Nu, 1); 
u = zeros(Nu, 1);

%% Initial conditions (2D raised cosine)
halfWidth = floor(min(Nx, Ny) / 5);
width = 2 * halfWidth + 1;
xLoc = floor(0.3 * Nx);
yLoc = floor(0.6 * Ny);
xRange = xLoc-halfWidth : xLoc+halfWidth;
yRange = yLoc-halfWidth : yLoc+halfWidth;

rcMat = zeros(Nyu, Nxu);
rcMat(yRange, xRange) = hann(width) * hann(width)';

% initialise current state  
u = reshape (rcMat, Nu, 1); 

% Set initial velocity to zero
uPrev = u;

% Output location
xOut = 0.45;
yOut = 0.25;
outLoc = round((xOut + yOut * Nyu) * Nxu);
out = zeros(lengthSound, 1);

%% Simulation loop
for n = 1:lengthSound
    
    %% Update equation
    uNext = (2 * eye(Nu) + c^2 * k^2 * D) * u - uPrev;
   
    % Update system states
    uPrev = u;
    u = uNext;
    
end
plot(out)