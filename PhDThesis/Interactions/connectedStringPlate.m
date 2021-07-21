close all;
clear all;

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs;   % Length of the simulation (1 second) [samples]             

%% Plotting
drawThings = false;     % Plot or not
drawSpeed = 10;         % 
plotSubplots = false;
calcEnergy = true;

%% Initialise variables for string
rhoS = 7850;
rS = 0.0005;
AS = pi * rS^2;
ES = 2e11;
IS = pi * rS^4 / 4; 
TS = 555;
cS = sqrt(TS / (rhoS * AS));
sig0S = 1;
sig1S = 0.005;
kappaSqS = ES * IS / (rhoS * AS);

LS = 1.5;
hS = sqrt((cS^2 * k^2 + 4 * sig1S * k + sqrt((cS^2 * k^2 + 4 * sig1S * k)^2 + 16 * kappaSqS * k^2))/2);
NS = floor(LS / hS);
hS = LS / NS;

%% Initialise variables for plate
rhoP = 7850;
EP = 2e11;
HP = 0.0005;
nu = 0.3;
sig0P = 1.0;
sig1P = 0.005;
Dvar = EP * HP^3 / (12 * (1-nu^2));
kappaSqP = Dvar / (rhoP * HP);

Lx = 1.5;
Ly = 1;
hP = 2 * sqrt(k * (sig1P + sqrt(kappaSqP + sig1P^2)));
Nx = floor(Lx / hP);
Ny = floor(Ly / hP);
hP = min(Lx/Nx, Ly/Ny);
muSqP = kappaSqP * k^2 / hP^4;

%% Initialise string matrices
IdS  = eye(NS-1);         % identity matrix

DxxS = toeplitz([-2, 1, zeros(1, NS-3)]) / hS^2;
DxxxxS = DxxS * DxxS;

Dxp = sparse(1:NS, 1:NS, -ones(1, NS), NS, NS) + ...
        sparse(1:NS-1, 2:NS, ones(1, NS-1), NS, NS);
Dxp = Dxp / hS;
AmatS = (1 + sig0S * k);
BS = 2 * IdS + cS^2 * k^2 * DxxS - kappaSqS * k^2 * DxxxxS + 2 * sig1S * k * DxxS;
CS = -(1 - sig0S * k) * IdS - 2 * sig1S * k * DxxS;


%% Initialise plate matrices
Nxw = Nx - 1;
Nyw = Ny - 1;
DxxP = toeplitz([-2, 1, zeros(1, Nxw-2)]);
DyyP = toeplitz([-2, 1, zeros(1, Nyw-2)]);
DxxEnP = toeplitz([-2, 1, zeros(1, Nxw)]);
DyyEnP = toeplitz([-2, 1, zeros(1, Nyw)]);
DEnPP = kron(speye(Nx+1), DyyEnP) + kron(DxxEnP, speye(Ny+1));
DEnPP = DEnPP / hP^2;
    
D = kron(speye(Nxw), DyyP) + kron(DxxP, speye(Nyw));
D = D / hP^2;
DD = D*D;

Nw = Nxw * Nyw;

AmatP = speye(Nw) * (1 + sig0P * k);
BP = 2 * speye(Nw) - kappaSqP * k^2 * DD + 2 * sig1P * k * D;
CP = -(1-sig0P * k) * speye(Nw) - 2 * sig1P * k * D;


%% Initialise states (simply supported boundary conditions)
amp = 5;
offset = amp * 0.5;
u = zeros(NS-1, 1);
halfWidth = 5;
loc = 0.33333;
width = halfWidth*2 + 1;
startLoc = floor(NS * loc) - halfWidth + 1;
u(startLoc:startLoc+width-1) = u(startLoc:startLoc+width-1) + amp * hann(width);
uPrev = u;

wNext = zeros(Nw, 1); 
excitationMat = zeros(Nyw, Nxw);
startX = floor(Nxw/2) - 5;
startY = floor(Nyw/2) - 5;
excitationMat(startX:startX+10, startY:startY+10) = hann(11) * hann(11)';
% w = reshape(excitationMat, Nw, 1);
w = zeros(Nw, 1);
wPrev = w;

%% Connection Variables
K1 = 1e4;
K3 = 1e7;
R = 10;

%% connection string
connLocU = 0.25;
xcu = floor(connLocU * NS);
alphaU = connLocU * NS - xcu;
xcu = xcu + 1; %matlab
Iu = zeros(1, NS-1);
Iu(xcu) = 1 - alphaU;
Iu(xcu+1) = alphaU;
Ju = 1/hS * Iu';

%% connection plate
connLocWX = 0.25;
connLocWY = 0.75;

xcwX = floor(connLocWX * Nxw);
alphaWX = connLocWX * Nxw - xcwX;
xcwX = xcwX + 1; %matlab

xcwY = floor(connLocWY * Nyw);
alphaWY = connLocWY * Nyw - xcwY;
xcwY = xcwY + 1; %matlab

Imat = zeros(Nyw, Nxw);
Imat(xcwY, xcwX) = (1 - alphaWX) * (1 - alphaWY);
Imat(xcwY, xcwX+1) = alphaWX * (1 - alphaWY);
Imat(xcwY+1, xcwX) = (1 - alphaWX) * alphaWY;
Imat(xcwY+1, xcwX+1) = alphaWX * alphaWY;

Iw = reshape(Imat, 1, Nw);
Jw = 1/hP^2 * Iw';

if drawThings
    plotNum = 1;
    figure('Position', [180 454 820 344])
end

out = zeros(lengthSound, 1);

%% energy initialisation
qTotU = 0;
qTotW = 0;
qTotC = 0;

kinEnergyU = zeros(lengthSound, 1);
potEnergyU = zeros(lengthSound, 1);
totEnergyU = zeros(lengthSound, 1);

q0U = zeros(lengthSound, 1);
q1U = zeros(lengthSound, 1);

kinEnergyW = zeros(lengthSound, 1);
potEnergyW = zeros(lengthSound, 1);
totEnergyW = zeros(lengthSound, 1);

q0W = zeros(lengthSound, 1);
q1W = zeros(lengthSound, 1);

connEnergy = zeros(lengthSound, 1);
qC = zeros(lengthSound, 1);
totEnergyC = zeros(lengthSound, 1);

totEnergy = zeros(lengthSound, 1);

%% Main Loop
for n = 1:lengthSound
    
    %% Calculate system states at connection locations without force term 
    uPre = AmatS \ (BS * u + CS * uPrev);
    wPre = AmatP \ (BP * w + CP * wPrev);
    
    uStar = Iu * uPre;
    wStar = Iw * wPre;

    % Relative displacement at n and n-1
    eta = Iu * u - Iw * w;
    etaPrev = Iu * uPrev - Iw * wPrev;

    %% Calculate force
    rPlus = K1 / 4 + K3 * eta^2 / 2 + R / (2*k);
    rMinus = K1 / 4 + K3 * eta^2 / 2 - R / (2*k);
    
    f = (uStar - wStar + K1 / (2 * rPlus) * eta + rMinus / rPlus * etaPrev) ...
            / (1/rPlus + Iu * Ju * k^2 / (rhoS * AS * (1+sig0S * k)) + Iw * Jw * k^2 / (rhoP * HP * (1+sig0P * k)));
    
    %% Calculate scheme
    uNext = uPre - Ju * k^2 * f / (rhoS * AS * (1 + sig0S * k));
    wNext = wPre + Jw * k^2 * f / (rhoP * HP * (1 + sig0P * k));
       
    out(n) = Iu * uNext;
    etaNext = Iu * uNext - Iw * wNext;
    if calcEnergy 
        %% energy string
        kinEnergyU(n) = rhoS * AS * hS / 2 * sum((1/k * (u - uPrev)).^2);
        potEnergyU(n) = TS / 2 * hS * sum((Dxp * [0; u]) .* (Dxp * [0; uPrev])) ...
            + ES * IS * hS / 2 * sum((DxxS * u) .* (DxxS * uPrev));

        % damping
        q0U(n) = 2 * sig0S * rhoS * AS * hS * sum((1/(2*k) * (uNext - uPrev)).^2);
        q1U(n) = - 2 * sig1S * rhoS * AS * hS * sum(1/(2*k) * (uNext - uPrev)...
            .* (1/k * (DxxS * u - DxxS * uPrev)));
        idx = n - (1 * (n~=1));
        qTotU = qTotU + k * (q0U(idx) + q1U(idx));

        % total energy 
        totEnergyU(n) = kinEnergyU(n) + potEnergyU(n) + qTotU;

        kinEnergyW(n) = rhoP * HP / 2 * hP^2 * sum((1/k * (w-wPrev)).^2);

    %         zeroU(2:end-1, 2:end-1) = reshapedU;
    %         zeroUPrev(2:end-1, 2:end-1) = reshapedUPrev;
    %         uEn = reshape(zeroU, (Nx+1) * (Ny+1), 1);
    %         uPrevEn = reshape(zeroUPrev, (Nx+1) * (Ny+1), 1);
        potEnergyW(n) = Dvar * hP^2 / 2 * sum((D * w) .* (D * wPrev));
        q0W(n) = 2 * sig0P * rhoP * HP * hP^2 * sum((1/(2*k) * (wNext - wPrev)).^2);
        q1W(n) = -2 * sig1P * rhoP * HP * hP^2 * sum(1/(2*k) * (wNext - wPrev)...
            .* (1/k * (D * w - D * wPrev)));
        idx = n - (1 * (n~=1));
        qTotW = qTotW + k * (q0W(idx) + q1W(idx));

        totEnergyW(n) = kinEnergyW(n) + potEnergyW(n) + qTotW;

        connEnergy(n) = K1 / 8 * (eta + etaPrev)^2 + K3 / 4 * (eta^2 * etaPrev^2);
        qC(n) = R * (1/(2*k) * (etaNext - etaPrev))^2;

        qTotC = qTotC + k * qC(idx);
        totEnergyC(n) = connEnergy(n) + qTotC;
        totEnergy(n) = totEnergyU(n) + totEnergyW(n) + totEnergyC(n);
    end
    
    %% Plot stuff
    if drawThings && ((~plotSubplots && mod(n, drawSpeed) == 0)...
            || plotSubplots && mod(n, 18) == 1)
        % system states
        if plotSubplots
            subplot(2, 3, plotNum);
        end
        hold off;
        
        u_s = plot3([0:NS] + xcwX +1 - xcu, xcwY * ones(NS+1), [0;u;0] + offset, 'r', 'Linewidth', 2);
        view(20,35);
        hold on;
        mesh(reshape(w, Nyw, Nxw))
        colormap gray;
        caxis([-0.05, 0.05]);
        plot3([xcwX, xcwX]+1, [xcwY, xcwY], [Iu * u + offset, Iw * w], 'color', [0, 0.85, 0], 'Linewidth', 2)
%         if plotNum == 1
%             legend(["$u_{l_u}^n$", "$w_{l_w}^n$"], 'interpreter', 'latex', ...
%                 'location', 'northwest')
%         end
%         plot([(xcu+alphaU), (xcu+alphaU)]/NS, [Iu * u, Iw * w - offset], ...
%             'color', [0.5, 0.5, 0.5], 'Linewidth', 2)
        

        
        zlim([-0.2, amp + offset])
        xticks([])
        yticks([])
        zticks([])
        title("$n = " + n + "$", 'interpreter', 'latex', 'Fontsize', 16)
        % energy
%         subplot(212)
%         plot(totEnergy(1:n) / totEnergy(1) - 1);
        box on
        grid on
        if plotSubplots

%             if plotNum == 1
%                 legend(["$u_{l_u}^n$", "$w_{l_w}^n$"], 'interpreter', 'latex', ...
%                     'location', 'northwest', 'Fontsize', 16)
%             end

            if plotNum<4
                set(gca, 'Position', [0.009+0.34*(plotNum-1) 0.53 0.31 0.4], ...
                    'Linewidth', 1, 'Projection', 'perspective')
            else
                set(gca, 'Position', [0.009+0.34*(plotNum-4) 0.0438 0.31 0.4], ...
                    'Linewidth', 1, 'Projection', 'perspective')
            end
                    plotNum = plotNum + 1;
            if plotNum > 6
                return;
            end
        end
        drawnow;
    end
%     elseif ~plotSubplots && drawThings
%         subplot(311)
%         hold off
%         plot((0:NS) / NS, [0; u; 0], 'r', 'Linewidth', 2)
%         hold on;
%         subplot(312)
% %         mesh(reshape(w, Nyw, Nxw))
% %         plot((0:Nw)/ Nw + (xcu+alphaU)/NS - (xcwX+alphaWX)/Nw, [0; uw(NS:end); 0] - offset, ...   
% %             'b', 'Linewidth', 2)
% %         plot((0:Nw)/ Nw + (xcu+alphaU)/NS - (xcwX+alphaWX)/Nw, [0; w; 0] - offset, ...   
% %             '--g', 'Linewidth', 2)
% %         subplot(212)
% %         hold off;
% %         plot(uw(1:NS-1) - u)
% %         hold on;
% %         plot(uw(NS:end)- w)
%         subplot(313)
% %         hold off
% %         plot(totEnergyS(1:n) - totEnergyS(1))
% %         hold on;
% %         plot(totEnergyP(1:n) - totEnergyP(1))
% %         plot(totEnergyC(1:n) - totEnergyC(1))
%         if calcEnergy
%             plot(totEnergy(1:n) / totEnergy(1) - 1)
%         end
%         drawnow;
%     end
    % Update States
    uPrev = u;
    u = uNext;
        
    wPrev = w;
    w = wNext;
    
end

plotEnergyConnected;
