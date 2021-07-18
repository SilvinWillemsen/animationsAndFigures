close all;
clear all;

fs = 44100;
lengthSound = fs;
k = 1/fs;
drawThings = true;

% initialise variables for u
rhou = 7850;
Au = pi * 0.0005^2;
cu = 300;
Tu = cu^2 * rhou * Au; 

Lu = 1;
hu = cu * k;
Nu = floor(Lu / hu);
hu = Lu / Nu;
lambdaSqu = cu^2 * k^2 / hu^2

% initialise variables for w
rhow = 1050;
Aw = pi * 0.0005^2;
cw = 500;
Tw = cw^2 * rhow * Aw; 

Lw = 1;
hw = cw * k;
Nw = floor(Lw / hw);
hw = Lw / Nw;
lambdaSqw = cw^2 * k^2 / hw^2

% states (dirichlet boundary conditions)
u = zeros(Nu-1, 1);
u(3:10) = hann(8);
uPrev = u;

w = zeros(Nw-1, 1);
wPrev = w;

Dxxu = 1/hu^2 * toeplitz([-2, 1, zeros(1, Nu-3)]);
Dxxw = 1/hw^2 * toeplitz([-2, 1, zeros(1, Nw-3)]);
Bu = 2 * eye(Nu-1) + cu^2 * k^2 * Dxxu;
Bw = 2 * eye(Nw-1) + cw^2 * k^2 * Dxxw;

connLocU = 0.25;
xcu = floor(connLocU * Nu);
alphaU = connLocU * Nu - xcu;
Iu = zeros(1, Nu-1);
Iu(xcu) = 1 - alphaU;
Iu(xcu+1) = alphaU;
Ju = 1/hu * Iu';

connLocW = 0.75;
xcw = floor(connLocW * Nw);
alphaW = connLocW * Nw - xcw;
Iw = zeros(1, Nw-1);
Iw(xcw) = 1 - alphaW;
Iw(xcw+1) = alphaW;
Jw = 1/hw * Iw';

%% Main Loop
for n = 1:lengthSound
    f = (cu^2 * Iu * Dxxu * u - cw^2 * Iw * Dxxw * w) / (Iu * Ju / (rhou * Au) + Iw * Jw / (rhow * Aw));
    uNext = Bu * u - uPrev - Ju * k^2 / (rhou*Au) * f;
    wNext = Bw * w - wPrev + Jw * k^2 / (rhow*Aw) * f;
   
    % energy
    kinEnergyU(n) = rhou * Au / 2 * hu * sum((1/k * (u-uPrev)).^2);
    potEnergyU(n) = Tu / (2*hu) * sum(([0; u] - [u; 0]) .* ([0; uPrev] - [uPrev; 0]));
    totEnergyU(n) = kinEnergyU(n) + potEnergyU(n);
    
    kinEnergyW(n) = rhow * Aw /2 * hw * sum((1/k * (w-wPrev)).^2);
    potEnergyW(n) = Tw / (2*hw) * sum(([0; w] - [w; 0]) .* ([0; wPrev] - [wPrev; 0]));
    totEnergyW(n) = kinEnergyW(n) + potEnergyW(n);

    totEnergy(n) = totEnergyU(n) + totEnergyW(n);
    
    %% Plot stuff
    if drawThings
        % system states
        subplot(211)
        hold off;
        plot([0; u; 0])
        hold on;
        plot([0; w; 0] + 0.5)

        % energy
        subplot(212)
        plot(totEnergy(1:n) / totEnergy(1) - 1);

        drawnow;
    end
    % Update States
    uPrev = u;
    u = uNext;
        
    wPrev = w;
    w = wNext;

end

plotEnergyConnected;
