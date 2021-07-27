% close all;
clear all;

fs = 44100;
lengthSound = fs;
k = 1/fs;
drawThings = false;
plotSubplots = true;
% initialise variables for u
rhou = 7850;
Au = pi * 0.0005^2;
cu = 300;
Tu = cu^2 * rhou * Au; 

Lu = 1;
hu = cu * k;
Nu = floor(Lu / hu);
hu = Lu / Nu;
lambdaSqu = cu^2 * k^2 / hu^2;

% initialise variables for w
rhow = 7850;
Aw = pi * 0.0005^2;
cw = 400;
Tw = cw^2 * rhow * Aw; 

Lw = 1;
hw = cw * k;
Nw = floor(Lw / hw);
hw = Lw / Nw;
lambdaSqw = cw^2 * k^2 / hw^2

% states (dirichlet boundary conditions)
u = zeros(Nu-1, 1);
halfWidth = 5;
width = halfWidth*2 + 1;
startLoc = floor(Nu * 0.5) - halfWidth + 1;
amp = 1;
u(startLoc:startLoc+width-1) = amp * hann(width);
uPrev = u;

w = zeros(Nw-1, 1);
startLoc = floor(Nw * 0.5) - halfWidth + 1;
% w(startLoc:startLoc+width-1) = amp * hann(width);
wPrev = w;

Dxxu = 1/hu^2 * toeplitz([-2, 1, zeros(1, Nu-3)]);
Dxxw = 1/hw^2 * toeplitz([-2, 1, zeros(1, Nw-3)]);
Bu = 2 * eye(Nu-1) + cu^2 * k^2 * Dxxu;
Bw = 2 * eye(Nw-1) + cw^2 * k^2 * Dxxw;

connLocU = 0.25;
xcu = floor(connLocU * Nu);
alphaU = connLocU * Nu - xcu;
xcu = xcu + 1; %matlab
Iu = zeros(1, Nu-1);
% cubInterp = [-alphaU * (alphaU - 1) * (alphaU - 2) / 6, ...
%     (alphaU - 1) * (alphaU + 1) * (alphaU - 2) / 2, ...
%     -alphaU * (alphaU + 1) * (alphaU - 2) / 2, ...
%     alphaU * (alphaU + 1) * (alphaU - 1) / 6]
% Iu(xcu-1:xcu+2) = cubInterp;
Iu(xcu) = 1 - alphaU;
Iu(xcu+1) = alphaU;
Ju = 1/hu * Iu';

connLocW = 0.75;
xcw = floor(connLocW * Nw);
alphaW = connLocW * Nw - xcw;
xcw = xcw + 1; %matlab
Iw = zeros(1, Nw-1);
Iw(xcw) = 1 - alphaW;
Iw(xcw+1) = alphaW;
Jw = 1/hw * Iw';

if drawThings
    plotNum = 1;
    figure('Position', [180 454 820 344])
end
offset = amp * 0.5;

uw = zeros(Nu + Nw - 2, 1);
uw(1:Nu-1) = u;
uwPrev = uw;

B = zeros(Nu + Nw - 2);
fDiv = Iu * Ju / (rhou * Au) + Iw * Jw / (rhow * Aw);

fVec = [k^2 / (rhou * Au) * cu^2 * Iu * Dxxu / fDiv, ...
        - k^2 / (rhow * Aw) * cw^2 * Iw * Dxxw / fDiv];
% puu = -Ju * k^2 / (rhou*Au) * cu^2 * Iu * Dxxu;
% puw = -Ju * k^2 / (rhou*Au) * -cw^2 * Iw * Dxxw;
% pww = Jw * k^2 / (rhow*Aw) * -cw^2 * Iw * Dxxw;
% pwu = Jw * k^2 / (rhow*Aw) * cu^2 * Iu * Dxxu;
B(1:Nu-1, 1:Nu-1) = Bu;
B(Nu:end, Nu:end) = Bw;
B = B + [-Ju ; Jw] * fVec;

Amat = speye(size(B));
C = -speye(size(B));

% plotModalAnalysis;
%% Main Loop
for n = 1:lengthSound
    f = (cu^2 * Iu * Dxxu * u - cw^2 * Iw * Dxxw * w) / (Iu * Ju / (rhou * Au) + Iw * Jw / (rhow * Aw));
    uNext = Bu * u - uPrev - Ju * k^2 / (rhou*Au) * f;
    wNext = Bw * w - wPrev + Jw * k^2 / (rhow*Aw) * f;
   
    uwNext = B * uw - uwPrev;
    
    % energy
    kinEnergyU(n) = rhou * Au / 2 * hu * sum((1/k * (u-uPrev)).^2);
    potEnergyU(n) = Tu / (2*hu) * sum(([0; u] - [u; 0]) .* ([0; uPrev] - [uPrev; 0]));
    totEnergyU(n) = kinEnergyU(n) + potEnergyU(n);
    
    kinEnergyW(n) = rhow * Aw /2 * hw * sum((1/k * (w-wPrev)).^2);
    potEnergyW(n) = Tw / (2*hw) * sum(([0; w] - [w; 0]) .* ([0; wPrev] - [wPrev; 0]));
    totEnergyW(n) = kinEnergyW(n) + potEnergyW(n);
    totEnergy(n) = totEnergyU(n) + totEnergyW(n);
    
    out(n) = uw(2);
    outU(n) = u(2);
    outW(n) = w(2);
    %% Plot stuff
    if drawThings && mod(n, 18) == 1 && plotSubplots
        % system states
        subplot(2, 3, plotNum);
        hold off;
        plot((0:Nu) / Nu, [0; u; 0], 'r', 'Linewidth', 2)
        hold on;
        plot((0:Nw)/ Nw + (xcu+alphaU)/Nu - (xcw+alphaW)/Nw, [0; w; 0] - offset, ...
            'b', 'Linewidth', 2)

        plot([(xcu+alphaU), (xcu+alphaU)]/Nu, [Iu * u, Iw * w - offset], ...
            'color', [0.5, 0.5, 0.5], 'Linewidth', 2)
        if plotNum == 1
            legend(["$u_l^n$", "$w_m^n$"], 'interpreter', 'latex', ...
                'location', 'northwest', 'Fontsize', 16)
        end
        ylim([-amp * 0.5-offset, amp * 1])
        xticks([])
        yticks([])
        title("$n = " + n + "$", 'interpreter', 'latex', 'Fontsize', 16)
        % energy
%         subplot(212)
%         plot(totEnergy(1:n) / totEnergy(1) - 1);
        
        
        if plotNum<4
            set(gca, 'Position', [0.009+0.33*(plotNum-1) 0.53 0.32 0.4], ...
                'Linewidth', 2)
        else
            set(gca, 'Position', [0.009+0.33*(plotNum-4) 0.0438 0.32 0.4], ...
                'Linewidth', 2)
        end
        plotNum = plotNum + 1;
        if plotNum > 6
            return;
        end
    elseif ~plotSubplots && drawThings
        subplot(211)
        hold off
        plot((0:Nu) / Nu, [0; uw(1:Nu-1); 0], 'r', 'Linewidth', 2)
        hold on;
        plot((0:Nu) / Nu, [0; u; 0], '--g', 'Linewidth', 2)
        plot((0:Nw)/ Nw + (xcu+alphaU)/Nu - (xcw+alphaW)/Nw, [0; uw(Nu:end); 0] - offset, ...   
            'b', 'Linewidth', 2)
        plot((0:Nw)/ Nw + (xcu+alphaU)/Nu - (xcw+alphaW)/Nw, [0; w; 0] - offset, ...   
            '--g', 'Linewidth', 2)
        subplot(212)
        hold off;
        plot(uw(1:Nu-1) - u)
        hold on;
        plot(uw(Nu:end)- w)
        drawnow;
    end
    % Update States
    uPrev = u;
    u = uNext;
        
    wPrev = w;
    w = wNext;
    
    uwPrev = uw;
    uw = uwNext;

end

plotEnergyConnected;
