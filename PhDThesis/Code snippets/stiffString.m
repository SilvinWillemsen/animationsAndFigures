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
T = 555;            % Tension [N]

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

% Boundary conditions ([c]lamped, [s]imply supported or [f]ree)
bc = "s";            

% Change number of intervals to the usable range
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
Id  = eye(N+1);
Dxp = sparse(1:Norig+1, 1:Norig+1, -ones(1, Norig+1), Norig+1, Norig+1) + ...
     sparse(1:Norig, 2:Norig+1, ones(1, Norig), Norig+1, Norig+1);
 
Dxp = Dxp / h;
DxxE = -Dxp' * Dxp;

DxxE = sparse(2:N+3, 1:N+2, ones(1, N+2), N+3, N+3) + ...
        sparse(1:N+3, 1:N+3, -2 * ones(1, N+3), N+3, N+3) + ... 
        sparse(1:N+2, 2:N+3, ones(1, N+2), N+3, N+3) ;
DxxE = DxxE / h^2;

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
elseif bc == "s"
%     DxxE(end, end-1) = 2 / h^2;
%     DxxE(1, 2) = 2 / h^2;
end

Amat = (1 + sig0 * k);
B = 2 * Id + c^2 * k^2 * Dxx - kappa^2 * k^2 * Dxxxx + 2 * sig1 * k * Dxx;
C = (-1 + sig0 * k) * Id - 2 * sig1 * k * Dxx;

%% Initial conditions (raised cosine)
ratio = 0.3;
loc = 20;       % Center location
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

[~, z, phi] = eig (c^2 * Dxx, 'vector');



%% Simulation loop
for n = 1:lengthSound
    
    uNextVec(range) = ((2 - 2 * lambdaSq - 6 * muSq - 4 * sig1 * k / h^2) * uVec(range)...
        + (lambdaSq + 4 * muSq + 2 * sig1 * k / h^2) * (uVec(range+1) + uVec(range-1)) ...
        - muSq * (uVec(range+2) + uVec(range-2)) + (-1 + sig0 * k + 4 * sig1 * k / h^2) * uPrevVec(range) ...
        - 2 * sig1 * k / h^2 * (uPrevVec(range+1) + uPrevVec(range-1))) / (1+sig0 * k);
    
    if bc == "s"
        uVirt1 = 2 * uVec(1) - uVec(2);       
        uVirtN1 = 2 * uVec(end) - uVec(end-1);

        uNextVec(2) = ((2 - 2 * lambdaSq - 6 * muSq - 4 * sig1 * k / h^2) * uVec(2)...
            + (lambdaSq + 4 * muSq + 2 * sig1 * k / h^2) * (uVec(3) + uVec(1)) ...
            - muSq * (uVec(4) + uVirt1) + (-1 + sig0 * k + 4 * sig1 * k / h^2) * uPrevVec(2) ...
            - 2 * sig1 * k / h^2 * (uPrevVec(3) + uPrevVec(1))) / (1+sig0 * k);
        
        uNextVec(end-1) = ((2 - 2 * lambdaSq - 6 * muSq - 4 * sig1 * k / h^2) * uVec(end-1)...
            + (lambdaSq + 4 * muSq + 2 * sig1 * k / h^2) * (uVec(end-2) + uVec(end)) ...
            - muSq * (uVec(end-3) + uVirtN1) + (-1 + sig0 * k + 4 * sig1 * k / h^2) * uPrevVec(end-1) ...
            - 2 * sig1 * k / h^2 * (uPrevVec(end-2) + uPrevVec(end))) / (1+sig0 * k);
    elseif bc == "f"
%         uVirt2 = uVec(2) - 2 * uVec(1) + 2 * uVirt1;
%         uVirtN2 = uVec(end-1) - 2 * uVec(end) + 2 * uVirtN1;
% 
%         uPrevVirt1 = 2 * uPrevVec(1) - uPrevVec(2);       
%         uPrevVirtN1 = 2 * uPrevVec(end) - uPrevVec(end-1);

        uNextVec(2) = ((2 - 2 * lambdaSq - 5 * muSq - 4 * sig1 * k / h^2) * uVec(2)...
            + (lambdaSq + 4 * muSq + 2 * sig1 * k / h^2) * (uVec(3)) ...
            + (lambdaSq + 2 * muSq + 2 * sig1 * k / h^2) * (uVec(1)) ...
            - muSq * uVec(4) + (-1 + sig0 * k + 4 * sig1 * k / h^2) * uPrevVec(2) ...
            - 2 * sig1 * k / h^2 * (uPrevVec(3) + uPrevVec(1))) / (1+sig0 * k);
        
        uNextVec(end-1) = ((2 - 2 * lambdaSq - 5 * muSq - 4 * sig1 * k / h^2) * uVec(end-1)...
            + (lambdaSq + 4 * muSq + 2 * sig1 * k / h^2) * (uVec(end-2)) ...
            + (lambdaSq + 2 * muSq + 2 * sig1 * k / h^2) * (uVec(end)) ...
            - muSq * uVec(end-3) + (-1 + sig0 * k + 4 * sig1 * k / h^2) * uPrevVec(end-1) ...
            - 2 * sig1 * k / h^2 * (uPrevVec(end-2) + uPrevVec(end))) / (1+sig0 * k);
        
        
        uNextVec(1) = ((2 - 2 * lambdaSq - 2 * muSq - 4 * sig1 * k / h^2) * uVec(1)...
            + (2 * lambdaSq + 4 * muSq + 4 * sig1 * k / h^2) * uVec(2) ...
            - 2 * muSq * uVec(3) + (-1 + sig0 * k + 4 * sig1 * k / h^2) * uPrevVec(1) ...
            - 4 * sig1 * k / h^2 * uPrevVec(2)) / (1+sig0 * k);
        
        uNextVec(end) = ((2 - 2 * lambdaSq - 2 * muSq - 4 * sig1 * k / h^2) * uVec(end)...
            + (2 * lambdaSq + 4 * muSq + 4 * sig1 * k / h^2) * uVec(end-1) ...
            - 2 * muSq * uVec(end-2) + (-1 + sig0 * k + 4 * sig1 * k / h^2) * uPrevVec(end) ...
            - 4 * sig1 * k / h^2 * uPrevVec(end-1)) / (1+sig0 * k);
        
%         uNextVec(1) = ((2 - 2 * lambdaSq - 6 * muSq - 4 * sig1 * k / h^2) * uVec(1)...
%             + (2 * lambdaSq + 4 * muSq + 4 * sig1 * k / h^2) * uVec(2) ...
%             - 2 * muSq * uVec(3) + (-1 + sig0 * k + 4 * sig1 * k / h^2) * uPrevVec(1) ...
%             - 4 * sig1 * k / h^2 * uPrevVec(2)) / (1+sig0 * k);
%         uNextVec(end) = ((2 - 2 * lambdaSq - 6 * muSq - 4 * sig1 * k / h^2) * uVec(end)...
%             + (2 * lambdaSq + 4 * muSq + 4 * sig1 * k / h^2) * uVec(end-1) ...
%             - 2 * muSq * uVec(end-2) + (-1 + sig0 * k + 4 * sig1 * k / h^2) * uPrevVec(end) ...
%             - 4 * sig1 * k / h^2 * uPrevVec(end-1)) / (1+sig0 * k);
    end
    
    uNext = B / Amat * u + C / Amat * uPrev;
    % Update equation 
    out(n) = u(outLoc);
    subplot(421)
    hold off;
    if bc == "c"
        uPlot = [0; 0; u; 0; 0];
    elseif bc == "s"
        uPlot = [0; u; 0];
    elseif bc == "f"
        uPlot = u;
    end
    plot(uPlot)
    hold on;
    plot(uVec)
    ylim([-1, 1])
    
    subplot(422)
    plot(uPlot - uVec)
    drawnow;
    
    
    % energy
    if bc == "f"
        scaling = ones(N+1, 1);
        scaling(1) = 0.5;
        scaling(end) = 0.5;
        
        kinEnergy(n) = rho * A * h / 2 * sum(scaling .* (1/k * (u - uPrev)).^2);
        potEnergy(n) = T / 2 * h * sum((1/h * (u(2:end) - u(1:end-1))) .* (1/h * (uPrev(2:end) - uPrev(1:end-1)))) ...
            + E * I * h / 2 * sum((1/h^2 * (u(3:end) - 2 * u(2:end-1) + u(1:end-2))) ...
             .* (1/h^2 * (uPrev(3:end) - 2 * uPrev(2:end-1) + uPrev(1:end-2))));
         
        totEnergy(n) = kinEnergy(n) + potEnergy(n);

    elseif bc == "s"
        kinEnergy(n) = rho * A * h / 2 * sum((1/k * (u - uPrev)).^2);
        potEnergy(n) = T / 2 * h * sum((Dxp * [0; u; 0]) .* (Dxp * [0; uPrev; 0])) ...
            + E * I * h / 2 * sum((1/h^2 * ([u(2:end); 0] - 2 * u + [0; u(1:end-1)])) ...
             .* (1/h^2 * ([uPrev(2:end); 0] - 2 * uPrev + [0; uPrev(1:end-1)])));
%         matForm = (DxxE * [0; u; 0]) .* (DxxE * [0; uPrev; 0]);
%         vecForm = (1/h^2 * ([u(2:end); 0] - 2 * u + [0; u(1:end-1)])) ...
%              .* (1/h^2 * ([uPrev(2:end); 0] - 2 * uPrev + [0; uPrev(1:end-1)]));
         
%         sum(matForm) - sum(vecForm)
        totEnergy(n) = kinEnergy(n) + potEnergy(n);
        subplot(412)

        plot(kinEnergy(1:n))
%         hold on;

        subplot(413)
%         plot(totEnergy(1:n))
                plot(potEnergy(1:n))

%         hold off;
%         plot(p1Test(1:n))
%         hold on;
%         plot(p2Test(1:n))
        
    elseif bc == "c"
        kinEnergy(n) = rho * A * h / 2 * sum((1/k * (u - uPrev)).^2);
        potEnergy(n) = T / 2 * h * sum((Dxp * [0; u; 0]) .* (Dxp * [0; uPrev; 0]))...
            + E * I * h / 2 * sum((DxxE * [0; u; 0]) .* (DxxE * [0; uPrev; 0]));
        totEnergy(n) = kinEnergy(n) + potEnergy(n);
%         hold off;
        
    end
    subplot(414)
    plot(totEnergy(1:n) / totEnergy(1) - 1)
    drawnow;

    % Update system states
    uPrev = u;
    u = uNext;
    
    uPrevVec = uVec;
    uVec = uNextVec;
end
fftOut = abs(fft(out));
plot([0:fs-1],fftOut)
