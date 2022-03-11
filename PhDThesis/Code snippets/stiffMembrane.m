close all;
clear all;

drawThings = true;
energyCalc = true;
plotPropagation = false;
if plotPropagation 
    figure('Position', [440 585 815 213])
end
drawSpeed = 1;
%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = 500;   % Length of the simulation (1 second) [samples]             

rho = 7850;
H =  0.005;
T = 10000;
E = 2e11;
nu = 0.3;
sig0 = 1;
sig1 = 0.005;
cSq = T / (rho * H);
Dparam = E * H^3/ (12*(1-nu^2));
kappaSq = Dparam / (rho * H); % Stiffness coefficient [m/s]

Lx = 2;               % Length in x direction [m]
Ly = 2;                 % Length in y direction [m]

h = sqrt(cSq * k^2 + 4 * sig1 * k + sqrt((cSq * k^2 + 4 * sig1 * k)^2 + 16 * kappaSq * k^2));    % Grid spacing [m] (from CFL condition)
Nx = floor(Lx/h);       % Number of intervals between grid points
Ny = floor(Ly/h);       % Number of intervals between grid points
h = min(Lx/Nx, Ly/Ny);  % Recalculation of grid spacing based on integer N

lambdaSq = cSq * k^2 / h^2;
muSq = kappaSq * k^2 / h^4; % Courant number squared
h = min(Lx/Nx, Ly/Ny);  % Recalculation of grid spacing based on integer N
% Boundary conditions ([C]lamped or [S]imply supported)
bc = "S";

% Prepare Dxx matrix
if bc == "S"
    Nxu = Nx - 1;
    Nyu = Ny - 1;
    Dxx = toeplitz([-2, 1, zeros(1, Nxu-2)]);
    Dyy = toeplitz([-2, 1, zeros(1, Nyu-2)]);
    DxxEn = toeplitz([-2, 1, zeros(1, Nxu)]);
    DyyEn = toeplitz([-2, 1, zeros(1, Nyu)]);
    DEn = kron(speye(Nx+1), DyyEn) + kron(DxxEn, speye(Ny+1));
    DEn = DEn / h^2;
% elseif bc == "C"
%     Nxu = Nx - 3;
%     Nyu = Ny - 3;
end    

D = kron(speye(Nxu), Dyy) + kron(Dxx, speye(Nyu));
D = D / h^2;
DD = D*D;

Nu = Nxu * Nyu;
    
%% Initialise state vectors (one more grid point than the number of intervals)
uNext = zeros(Nu, 1); 
u = zeros(Nu, 1);

%% Initial conditions (raised cosine)
% halfWidth = floor(min(Nx, Ny) / 10);
halfWidth = 1;
width = 2 * halfWidth + 1;
% xLoc = floor(0.25 * Nx);
% yLoc = floor(0.5 * Ny);
xLoc = 40;
yLoc = 40;
xRange = xLoc-halfWidth : xLoc+halfWidth;
yRange = yLoc-halfWidth : yLoc+halfWidth;

rcMat = zeros(Nyu, Nxu);
rcMat(yRange, xRange) = 1 * hann(width) * hann(width)';
u = reshape(rcMat, Nu, 1); % initialise current state  

% Set initial velocity to zero
uPrev = u;

% Output location and scaling of boundary points for energy
xOut = 0.45;
yOut = 0.25;
outLoc = round((xOut + yOut * Nyu) * Nxu);

out = zeros(lengthSound, 1);

kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

reshapedUPrev = reshape(uPrev, Nyu, Nxu);

Amat = speye(Nu) * (1 + sig0 * k);
B = 2 * speye(Nu) + cSq * k^2 * D - kappaSq * k^2 * DD + 2 * sig1 * k * D;
C = -(1-sig0 * k) * speye(Nu) - 2 * sig1 * k * D;

firstPlot = true;
N = Nu-1;
% [phi, lamb] = eig(full(D), 'vector');
% 
plotDampingAgainstFrequency = true;
% plotModalAnalysis;
percentCounter = 0;
nCounter = 0;
zeroU = zeros(Ny+1, Nx+1);
zeroUPrev = zeros(Ny+1, Nx+1);

qTot = 0;
plotNum = 1
%% Simulation loop
for n = 1:lengthSound
    
    %% Update equation
    uNext = Amat \ B * u + Amat \ C * uPrev;
    
    if energyCalc || drawThings
        reshapedU = reshape(u, Nyu, Nxu);
    end
    
    if energyCalc
        %% Energy
        kinEnergy(n) = rho * H / 2 * h^2 * sum((1/k * (u-uPrev)).^2);
        
%         zeroU(2:end-1, 2:end-1) = reshapedU;
%         zeroUPrev(2:end-1, 2:end-1) = reshapedUPrev;
%         uEn = reshape(zeroU, (Nx+1) * (Ny+1), 1);
%         uPrevEn = reshape(zeroUPrev, (Nx+1) * (Ny+1), 1);
        potEnergy(n) = T/2 * (sum(sum( ...
                ([zeros(Nyu, 1), reshapedU] - [reshapedU, zeros(Nyu, 1)]) ...
                .* ([zeros(Nyu, 1), reshapedUPrev] - [reshapedUPrev, zeros(Nyu, 1)]))) ...
                + sum(sum(([zeros(1, Nxu); reshapedU] - [reshapedU; zeros(1, Nxu)]) ...
                .* ([zeros(1, Nxu); reshapedUPrev] - [reshapedUPrev; zeros(1, Nxu)])))) ...
                + Dparam * h^2 / 2 * sum((D * u) .* (D * uPrev));
        q0(n) = 2 * sig0 * rho * H * h^2 * sum((1/(2*k) * (uNext - uPrev)).^2);
        q1(n) = -2 * sig1 * rho * H * h^2 * sum(1/(2*k) * (uNext - uPrev)...
            .* (1/k * (D * u - D * uPrev)));
        idx = n - (1 * (n~=1));
        qTot = qTot + k * (q0(idx) + q1(idx));

        totEnergy(n) = kinEnergy(n) + potEnergy(n) + qTot;

        out(n) = u(outLoc);
         
    end
    if drawThings && (mod(n, drawSpeed) == 0 || (plotPropagation && (n == 1 || n == 101 || n == 201)))
        
        if energyCalc
            subplot(211)
        end
        
        if plotPropagation
            subplot(1,3,plotNum)
        end
        imagesc(reshapedU);
        ax = gca;
        colormap gray;
        ax.CLim = [-0.5, 0.5];
        if plotPropagation
            plotNum = plotNum + 1;
            pause
        end
       
        if energyCalc
            subplot(212)
            hold off;
%             plot(kinEnergy(1:n))
%             hold on;
%             plot(potEnergy(1:n))
            plot(totEnergy(1:n) / totEnergy(1)  - 1)
        end
        drawnow;
    end
    % Update system states
    uPrev = u;
    u = uNext;
    if energyCalc
        reshapedUPrev = reshapedU;
    end
    
    if n > lengthSound * nCounter / 100
        percentCounter = percentCounter + 1;
        nCounter = nCounter + 1;
        disp((percentCounter) + "% done")
    end
end

plotEnergy;