close all;
clear all;
drawThings = true;
drawEnergy = false;
drawSubplots = false;
drawSpeed = 150;
drawStart = 1;
%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = 150;   % Length of the simulation (1 second) [samples]             

%% Mass
M = 1;
u = -0.95;
uPrev = -1;

b = 0;

%% Collision potential
Kc = 1e9;
alfC = 1.3;
psiPrev = 0;

subp = 1;
if drawThings
    figure('Position', [173 612 332 186])
end

%% Simulation loop
for n = 1:lengthSound
    eta = u - b;
    etaPrev = uPrev - b;
    if psiPrev >= 0 
        kappaG = 1;
    else
        kappaG = -1;
    end
    
    %% uNext and wNext without the collision term
    uStar = 2 * u - uPrev;
    if eta < 0
        etaStar = uStar - b;
        if etaStar - etaPrev == 0
            g = 0;
        else
            % Update equation 
            g = -2 * psiPrev / (etaStar - etaPrev);
        end
    else
        g = kappaG * sqrt(Kc * (alfC + 1) / 2) * eta^((alfC - 1) / 2);
    end 
    
    uNext(n) = (M / k^2 * uStar + g^2 / 4 * uPrev - psiPrev * g) / (M / k^2 + g^2 / 4);  
    etaNext = uNext(n) - b;

    psi = psiPrev + g / 2 * (etaNext - etaPrev);
%     psi=0;
    % Retrieve output
    
%     hold off;
%     plot(u)
%     hold on
%     scatter(massLoc, w)
        
    % energy in the string

    kinEnergy(n) = M / 2 * (1/k * (u - uPrev))^2;

    connEnergy(n) = psiPrev^2 / 2; % should be col energy but for plotEnergy.m it's handier
    totEnergy(n) = kinEnergy(n) + connEnergy(n);
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        if drawSubplots
            subplot(1, 3, subp)
        end
        hold off;
        %...mass
        xlim([1, n])
        xLim = xlim;
        ylim([-4, 3])
        yLim = ylim;
        yLimDiff = yLim(2) - yLim(1);
        patchOffsetX = 0.5;
        patchOffsetY = 0.007*yLimDiff;
        patch([xLim(1)+patchOffsetX, xLim(2)-patchOffsetX*0.5, xLim(2)-patchOffsetX*0.5, xLim(1)+patchOffsetX], ...
            [yLim(2)-patchOffsetY, yLim(2)-patchOffsetY, 0, 0], [0.75, 0.75, 0.75], ...
            'Edgecolor', 'none')
%         plot(xLim, [0,0], 'color', [0.5, 0.5, 0.5], 'Linewidth', 2)
        hold on;
        plot(uNext(1:n), 'b', 'Linewidth', 2);   

%         if Kc == 1e8
            legend('$b$', '$u^n$', ...
            'interpreter', 'latex', 'Location', 'northeast')
%         end
        str = extractAfter(num2str(Kc,2),"e+")
        if eval(str) < 10
            str = extractAfter(num2str(Kc,2),"0")
        end
        
%         title("$K_\textrm{\fontsize{7}{7}\selectfont c} = 10^" + str + "$", 'interpreter', 'latex')
        
        xLim = xlim;
        xLab = xlabel("$n$", 'interpreter', 'latex', 'Fontsize', 16)
        xLab.Position(2) = yLim(1) - 0.05 * (yLim(2) - yLim(1));
        
        yLab = ylabel("Displacement (m)", 'Fontname', 'times', 'Fontsize', 16)
        yLab.Position(1) = xLim(1) - 0.10 * (xLim(2) - xLim(1));
        box on;
        set(gca, 'Linewidth', 1.5, 'Fontsize', 16, ...
            'TickLabelInterpreter', 'latex', 'Position', [0.1416 0.1310 0.8072 0.8475]) 
        
        drawnow;
    end

    % Update system states
    uPrev = u;
    u = uNext(n);
    
    psiPrev = psi;
    
end
% plotEnergyMassBarrier;