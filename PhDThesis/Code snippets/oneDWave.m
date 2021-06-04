% close all;
clear all;

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs;   % Length of the simulation (1 second) [samples]             

rho = 7850;
A = pi * 0.0005^2;
c = 300;            % Wave speed [m/s]
T = c^2 * rho * A; 

L = 1;              % Length [m]
h = c * k;          % Grid spacing [m] (from CFL condition)
N = floor(L/h);     % Number of intervals between grid points
h = L / N;          % Recalculation of grid spacing based on integer N

lambdaSq = c^2 * k^2 / h^2; % Courant number squared

% Boundary conditions ([D]irichlet or [N]eumann)
bcLeft = "D";            
bcRight = "D"; 

if bcLeft == "D" && bcRight == "D"
    Dxx = toeplitz([-2, 1, zeros(1, N-3)]);
    I = eye(N-1);
    Nu = N-2;
elseif bcLeft == "N" && bcRight == "N"
    Dxx = toeplitz([-2, 1, zeros(1, N-1)]);
    Dxx(1, 2) = 2;
    Dxx(end, end-1) = 2;
    I = eye(N+1);
    Nu = N;
elseif bcLeft ~= bcRight
    Dxx = toeplitz([-2, 1, zeros(1, N-2)]);
    if bcLeft == "N"
        Dxx(1, 2) = 2;
    else
        Dxx(end, end-1) = 2;
    end
    I = eye(N);
    Nu = N-1;
end
Dxx = Dxx / h^2;
    
%% Initialise state vectors (one more grid point than the number of intervals)
uNext = zeros(Nu+1, 1); 
u = zeros(Nu+1, 1);

%% Initial conditions (raised cosine)
loc = round(0.5 * Nu);       % Center location
halfWidth = round(Nu/2 - 1);    % Half-width of raised cosine
width = 2 * halfWidth;      % Full width
rcX = 0:width;              % x-locations for raised cosine

rc = 0.5 - 0.5 * cos(2 * pi * rcX / width); % raised cosine
u(loc-halfWidth : loc+halfWidth) = rc; % initialise current state  

% Set initial velocity to zero
uPrev = u;

% Range of calculation
range = 2:N;

% Output location
outLoc = round(0.3 * N);

[~, z, phi] = eig (c^2 * Dxx, 'vector');

%% Simulation loop
for n = 1:lengthSound
    
    uNext = (2 * I + c^2 * k^2 * Dxx) * u - uPrev;
    plot(u)
%     drawnow;
    % Update equation 
%     uNext(range) = (2 - 2 * lambdaSq) * u(range) - uPrev(range) ...
%         + lambdaSq * (u(range+1) + u(range-1)); 
%     
%     % boundary updates
%     if bcLeft == "N"
%         uNext(1) = (2 - 2 * lambdaSq) * u(1) - uPrev(1) ...
%         + 2 * lambdaSq * u(2); 
%     end
%     
%     if bcRight == "N"
%         uNext(N+1) = (2 - 2 * lambdaSq) * u(N+1) - uPrev(N+1) ...
%         + 2 * lambdaSq * u(N); 
%     end
%     
    kinEnergy(n) = rho * A/2 * h * sum((1/k * (u-uPrev)).^2);
    potEnergy(n) = T/(2*h) * sum((u(2:end) - u(1:end-1)) .* (uPrev(2:end) - uPrev(1:end-1)));
    hold off;
    plot(kinEnergy(1:n))
    hold on;
    plot(potEnergy(1:n));
    drawnow;
    out(n) = u(outLoc);
    
    % Update system states
    uPrev = u;
    u = uNext;
end
fftOut = abs(fft(out));
plot([0:fs-1],fftOut)
