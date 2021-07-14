% close all;
clear all;

drawThings = false;
energyCalc = false;
plotPropagation = false;
if plotPropagation 
    figure('Position', [440 606 815 192])
    figWidth = 0.26;
    inc = 0.33;
    start = 0.04;
    height = 0.72;
    plotNum = 1;
end
drawSpeed = 1;
%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs;   % Length of the simulation (1 second) [samples]             

rho = 7850;
H =  0.005;
E = 2e11;
nu = 0.3;
sig0 = 1;
sig1 = 0.005;
Dparam = E * H^3/ (12*(1-nu^2));
kappaSq = Dparam / (rho * H); % Stiffness coefficient [m/s]

Lx = 1.5;               % Length in x direction [m]
Ly = 1;                 % Length in y direction [m]

h = 2 * sqrt(k * (sig1 + sqrt (sig1^2 + kappaSq)));    % Grid spacing [m] (from CFL condition)
Nx = floor(Lx/h);       % Number of intervals between grid points
Ny = floor(Ly/h);       % Number of intervals between grid points
h = min(Lx/Nx, Ly/Ny);  % Recalculation of grid spacing based on integer N

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
end    

D = kron(speye(Nxu), Dyy) + kron(Dxx, speye(Nyu));
D = D / h^2;
DD = D*D;

Nu = Nxu * Nyu;
    
%% Initialise state vectors (one more grid point than the number of intervals)
uNext = zeros(Nu, 1); 
u = zeros(Nu, 1);

%% Initial conditions (raised cosine)
halfWidth = floor(min(Nx, Ny) / 10);
width = 2 * halfWidth + 1;
xLoc = floor(0.25 * Nx);
yLoc = floor(0.5 * Ny);
xRange = xLoc-halfWidth : xLoc+halfWidth;
yRange = yLoc-halfWidth : yLoc+halfWidth;

rcMat = zeros(Nyu, Nxu);
rcMat(yRange, xRange) = 2.5 * hann(width) * hann(width)';
u = reshape(rcMat, Nu, 1); % initialise current state  

% Set initial velocity to zero
uPrev = u;

% Output location and scaling of boundary points for energy
xOut = 0.15;
yOut = 0.85;
xIdx = round(xOut * Nxu);
yIdx = round(yOut*Nyu);
outLoc = yIdx + xIdx * Nyu + 1;

out = zeros(lengthSound, 1);

kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

reshapedUPrev = reshape(uPrev, Nyu, Nxu);

Amat = speye(Nu) * (1 + sig0 * k);
B = 2 * speye(Nu) - kappaSq * k^2 * DD + 2 * sig1 * k * D;
C = -(1-sig0 * k) * speye(Nu) - 2 * sig1 * k * D;

firstPlot = true;
N = Nu-1;
% [phi, lamb] = eig(full(D), 'vector');
% 
% plotModalAnalysis;
percentCounter = 0;
nCounter = 0;
zeroU = zeros(Ny+1, Nx+1);
zeroUPrev = zeros(Ny+1, Nx+1);

qTot = 0;
if plotPropagation
    reshapedUPre = zeros(Ny+1, Nx+1);
end
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
        potEnergy(n) = Dparam * h^2 / 2 * sum((D * u) .* (D * uPrev));
        q0(n) = 2 * sig0 * rho * H * h^2 * sum((1/(2*k) * (uNext - uPrev)).^2);
        q1(n) = -2 * sig1 * rho * H * h^2 * sum(1/(2*k) * (uNext - uPrev)...
            .* (1/k * (D * u - D * uPrev)));
        idx = n - (1 * (n~=1));
        qTot = qTot + k * (q0(idx) + q1(idx));

        totEnergy(n) = kinEnergy(n) + potEnergy(n) + qTot;

        out(n) = u(outLoc);
         
    end
    if drawThings && (~plotPropagation && mod(n, drawSpeed) == 0) || (plotPropagation && (n == 1 || n == 101 || n == 201))
        
        if energyCalc
            subplot(211)
        end
        
        if plotPropagation
            subplot(1,3,plotNum)
            reshapedUPre(2:end-1, 2:end-1) = reshapedU;
            reshapedU = reshapedUPre;
        end
        imagesc(reshapedU);
        ax = gca;
        colormap gray;
        ax.CLim = [-0.45, 0.45];
        if plotPropagation
            xLab = xlabel('$l$', 'interpreter', 'latex', 'Fontsize', 16)
            yLab = ylabel('$m$', 'interpreter', 'latex', 'Fontsize', 16)
            xLab.Position(2) = Ny * 1.08;
            yLab.Position(1) = Nx * -0.05;
            title("$n = " + n + "$", 'interpreter', 'latex', 'Fontsize', 16)
            colormap(gray);
            
            set(gca,'xtick',[1, Nx+1],'xticklabel', ["$0$", "$N_x$"], ...
                'ytick',[1, Ny+1],'yticklabel', ["$0$", "$N_y$"], ...
                'ticklabelinterpreter', 'latex', 'FontSize', 16)
            set(gca, 'Position', [start+(plotNum-1)*inc 0.1500 figWidth height])           
            plotNum = plotNum + 1;
            if plotNum > 3
                return;
            end
%             pause
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
    
    out(n) = u(outLoc);
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

% plotEnergy;

%% Plotting
figure('Position', [173 578 827 220])

t = 1:lengthSound;
% plot(t, zeros(length(t), 1), '--', 'Linewidth', 2, 'color', [0.5, 0.5, 0.5])
% hold on;
subp1 = subplot(1, 2, 1);
plot(t, out, 'k', 'Linewidth', 2)

xlim([0, 6000])
ylim([-1.1, 1.1])
xLab = xlabel("$n$", 'interpreter', 'latex');
yLab = ylabel("$u^n_{" + xIdx + "," + yIdx + "}$", 'interpreter', 'latex');

xLab.Position(2) = -1.35;
% yLab.Position(1) = -15;
% grid on;
xLim = xlim;
yLab.Position(1)  = xLim(1) -0.07 * (xLim(2) - xLim(1));


set(gca, 'Linewidth', 1.5, 'Fontsize', 16, ...
    'Position', [0.0532 0.2000 0.4115 0.7591], ...
    'TickLabelInterpreter', 'latex')
subplot(1, 2, 2)
dbFFT = 20 * log10(abs(fft(out)));
plot(0:lengthSound-1, dbFFT, 'k', 'Linewidth', 2);
xlim([0, 5000])
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