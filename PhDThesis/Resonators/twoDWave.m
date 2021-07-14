close all;
clear all;

drawThings = false;
energyCalc = false;
drawSpeed = 1;

plotPropagation = false;
if plotPropagation 
    figure('Position', [440 606 815 192])
    figWidth = 0.26;
    inc = 0.33;
    start = 0.04;
    height = 0.72;
    plotNum = 1;

end

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs;   % Length of the simulation (1 second) [samples]             

rho = 7850;
H =  0.0005;
T = 500000;
c = sqrt(T / (rho * H));            % Wave speed [m/s]
c = 0.1/ (sqrt(2)*k);
Lx = 1.5;                 % Length in x direction [m]
Ly = 1;                 % Length in y direction [m]

h = sqrt(2) * c * k;    % Grid spacing [m] (from CFL condition)
Nx = floor(Lx/h);       % Number of intervals between grid points
Ny = floor(Ly/h);       % Number of intervals between grid points
h = min(Lx/Nx, Ly/Ny);  % Recalculation of grid spacing based on integer N

lambdaSq = c^2 * k^2 / h^2; % Courant number squared
h = max(Lx/Nx, Ly/Ny);  % Recalculation of grid spacing based on integer N

% Boundary conditions ([D]irichlet or [N]eumann)
bc = "D";

% Prepare Dxx matrix
if bc == "D"
    Nxu = Nx - 1;
    Nyu = Ny - 1;
    Dxx = toeplitz([-2, 1, zeros(1, Nxu-2)]);
    Dyy = toeplitz([-2, 1, zeros(1, Nyu-2)]);
end    

D = kron(speye(Nxu), Dyy) + kron(Dxx, speye(Nyu));

Nu = Nxu * Nyu;
D = D / h^2;
    
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
rcMat(yRange, xRange) = 4.5 * hann(width) * hann(width)';
u = reshape(rcMat, Nu, 1); % initialise current state  

% Set initial velocity to zero
uPrev = u;

% Output location and scaling of boundary points for energy
xOut = 0.15;
yOut = 0.85;
xIdx = round(xOut * Nxu);
yIdx = round(yOut*Nyu);
outLoc = yIdx + xIdx * Nyu + 1;
testVec = zeros(Nu, 1)
testVec(outLoc) = 1
testVecReshaped = reshape(testVec, Nyu, Nxu);
imagesc(testVecReshaped)
find(testVecReshaped == 1)
if bc == "D"
%     scaling = ones(N-1, 1);
elseif bc == "N"
%     scaling = [0.5; ones(N-1, 1); 0.5];
end
out = zeros(lengthSound, 1);

kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

reshapedUPrev = reshape(uPrev, Nyu, Nxu);

B = (2 * speye(Nu) + c^2 * k^2 * D);
Amat = speye(Nu);
C = -speye(Nu);

firstPlot = true;
N = Nu-1;
[phi, lamb] = eig(full(D), 'vector');
% fp = 1/(pi * k) * asin(c * k/ 2 * sqrt(-eig(D)));
plot2DmodeShapes;
noDamping = true;
plotModalAnalysis;
percentCounter = 0;
nCounter = 0;

if plotPropagation
    reshapedUPre = zeros(Ny+1, Nx+1);
end
%% Simulation loop
for n = 1:lengthSound
    
    %% Update equation
    uNext = B * u - uPrev;
    
    if energyCalc || drawThings || plotPropagation
        reshapedU = reshape(u, Nyu, Nxu);

    end
    if energyCalc
        %% Energy
        kinEnergy(n) = rho * H / 2 * h^2 * sum((1/k * (u-uPrev)).^2);


        if bc == "D"
            potEnergy(n) = T/2 * (sum(sum( ...
                ([zeros(Nyu, 1), reshapedU] - [reshapedU, zeros(Nyu, 1)]) ...
                .* ([zeros(Nyu, 1), reshapedUPrev] - [reshapedUPrev, zeros(Nyu, 1)]))) ...
                + sum(sum(([zeros(1, Nxu); reshapedU] - [reshapedU; zeros(1, Nxu)]) ...
                .* ([zeros(1, Nxu); reshapedUPrev] - [reshapedUPrev; zeros(1, Nxu)]))));
        else
            potEnergy(n) = T/(2*h) * sum((u(2:end) - u(1:end-1))...
                .* (uPrev(2:end) - uPrev(1:end-1)));
        end
        totEnergy(n) = kinEnergy(n) + potEnergy(n);

        out(n) = u(outLoc);
         
    end
    if drawThings && (~plotPropagation && mod(n, drawSpeed) == 0) || (plotPropagation && (n == 1 || n == 51 || n == 101))
        
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
    
    nCounter = nCounter + 1;
    if nCounter > lengthSound / 100
        percentCounter = percentCounter + 1;
        nCounter = 0;
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

xlim([0, 600])
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