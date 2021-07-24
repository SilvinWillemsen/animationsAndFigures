%{
    First order brass with coupled lip model
%}

clear all;
close all;

% drawing variables
drawThings = true;
drawSpeed = 10;
drawStart = 0;
centered = true;
radiation = true;
plotEnergy = false;
plotSubplots = false;
fs = 44100;             % Sample rate (Hz)
k = 1/fs;               % Time step (s)
lengthSound = 1000;   % Duration (s)

%% viscothermal effects
T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants (T);

%% Tube variables
h = c * k;              % Grid spacing (m)
L = 2.99;                  % Length

N = floor(L/h);         % Number of points (-)
h = L/N;                % Recalculate gridspacing from number of points

lambda = c * k / h      % courant number

%% Lip Collision
Kcol = 10000;
alfCol = 1; 
LnonExtended = 2.593;
NnonExtended = LnonExtended / h;

%% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N+1, NnonExtended, true);

%% Lip variables
f0 = 600;                   % fundamental freq lips
M = 5.37e-5;                % mass lips
omega0 = 2 * pi * f0;   % angular freq
K = omega0^2 * M;
sig = 5;                % damping
R = sig * M;
H0 = 2.9e-4;                % equilibrium

y = -H0;                      % initial lip state
yPrev = -H0;                  % previous lip state

w = 1e-2;                   % lip width
Sr = 1.46e-5;               % lip area

%% Initialise states
pNext = zeros(N+1, 1);        % pressure
p = zeros(N+1, 1);
% p(N/2-5 : N/2+5) =100 * hann(11);
vNext = zeros(N, 1);      % velocity
v = zeros(N, 1);

amp = 2000;                 % input pressure (Pa)

in = zeros(lengthSound, 1);
% p(floor(2*N / 3) - 5 : floor(2*N / 3) + 5) = 10000 * hann(11);

% Initialise output
out = zeros (lengthSound, 1);
outputPos = floor(4/5 * N);

% Set ranges
pRange = 2:N;         % range without boundaries
vRange = 1:N;         % range from 1/2 - N-1/2

% Virtual points
SNph = 2 * SBar(N+1) - SHalf(end);
SOnemh = 2 * SBar(1) - SHalf(1);

%% Initialise energies
kinEnergyTube = zeros (lengthSound, 1);
potEnergyTube = zeros (lengthSound, 1);
totEnergy = zeros (lengthSound, 1);
hTube = zeros (lengthSound, 1);
hReed = zeros (lengthSound, 1);
hColl = zeros (lengthSound, 1);
hRad = zeros (lengthSound, 1);

qReed = zeros (lengthSound, 1);
qHReed = zeros (lengthSound, 1);
pReed = zeros (lengthSound, 1);
pHReed = zeros (lengthSound, 1);
qRad = zeros (lengthSound, 1);
qHRad = zeros (lengthSound, 1);


totEnergy = zeros (lengthSound, 1);
scaledTotEnergy = zeros (lengthSound, 1);

scaling = ones(N+1,1);
if centered
    scaling(1) = 1 / 2;
    scaling(N+1) = 1 / 2;
end

%
psiPrev = 0;
etaC = 0;

%% Radiation impedance
R1 = rho * c;
rL = sqrt(SBar(N+1)) / (2 * pi);
Lr = 0.613 * rho * rL;
R2 = 0.505 * rho * c;
Cr = 1.111 * rL / (rho * c^2); 

zDiv = 2 * R1 * R2 * Cr + k * (R1 + R2);
if zDiv == 0
    z1 = 0;
    z2 = 0;
else
    z1 = 2 * R2 * k / zDiv;
    z2 = (2 * R1 * R2 * Cr - k * (R1 + R2)) / zDiv;

end
z3 = k/(2*Lr) + z1 / (2 * R2) + Cr * z1 / k;
z4 = (z2 + 1) / (2 * R2) + (Cr * z2 - Cr) / k;
    
p1 = 0;
v1 = 0;

if plotSubplots
    subp = 1;
    figure('Position', [1 467 1048 390])
    inc = 0.33;
end
for n = 1:lengthSound
    %% Calculate velocities before lip model
    vNext = v - lambda / (rho * c) * (p(vRange+1) - p(vRange));
    
    %% Variable input force
%     ramp = 10;
%     if n < ramp
%         Pm = amp * (n-1) / ramp;
%     else
        Pm = amp;
%     end
%     Pm = 0;
    %% Collision
    barr = -H0;
    etaC = barr - y;
    g = 0;
    if alfCol == 1
        if etaC > 0
            g = sqrt(Kcol * (alfCol+1) / 2);
        end
    else
        g = sqrt(Kcol * (alfCol+1) / 2) * subplus(etaC)^((alfCol - 1)/2);
    end
    
    %% Obtain deltaP
    a1 = 2 / k + omega0^2 * k + sig + g^2 * k / (2 * M);
    a2 = Sr / M;
    a3 = 2/k * 1/k * (y - yPrev) - omega0^2 * yPrev + g / M * psiPrev;
    b1 = SHalf(1) * vNext(1) + h * SBar(1) / (rho * c^2 * k) * (Pm  - p(1));
    b2 = h * SBar(1) / (rho * c^2 * k);
    c1 = w * subplus(y + H0) * sqrt(2 / rho);
    c2 = b2 + a2 * Sr / a1;
    c3 = b1 - a3 * Sr / a1;
    
    deltaP = sign(c3) * ((-c1 + sqrt(c1^2 + 4 * c2 * abs(c3)))/ (2 * c2))^2;
    
    %% Update lip scheme
    gammaR = g * k^2 / (2 * M);
    alpha = 2 + omega0^2 * k^2 + sig * k + g * gammaR;
    beta = sig * k - 2 - omega0^2 * k^2 + g * gammaR;
    xi = 2 * Sr * k^2 / M;

    yNext(n) = 4 / alpha * y + beta / alpha * yPrev + xi / alpha * deltaP + 4 * gammaR * psiPrev / alpha;

    %% Update collision potential
    psi = psiPrev - 0.5 * g * (yNext(n) - yPrev);
    
    %% Calculate flow velocities
    Ub = w * subplus(y + H0) * sign(deltaP) * sqrt(2 * abs(deltaP)/rho);
    Ur = Sr * 1/(2*k) * (yNext(n) - yPrev);

    %% Calculate pressure
    pNext(pRange) = p(pRange) - rho * c * lambda ./ SBar(pRange) .* (SHalf(pRange) .* vNext(pRange) - SHalf(pRange-1) .* vNext(pRange-1));
    pNext(1) = p(1) - rho * c * lambda / SBar(1) .* (-2 * (Ub + Ur) + 2 * SHalf(1) * vNext(1));
    if radiation
        pNext(N+1) = ((1 - rho * c * lambda * z3) * p(N+1) - 2 * rho * c * lambda * (v1 + z4 * p1 - (SHalf(end) .* vNext(end))/SBar(N+1))) / (1 + rho * c * lambda * z3);
    %     pNext(N+1) = alphaR * p(N+1) + betaR * SHalf(end) * vNext(end) + epsilonR * v1 + nuR * p1;
        v1Next = v1 + k / (2 * Lr) * (pNext(N+1) + p(N+1));
        p1Next = z1 / 2 * (pNext(N+1) + p(N+1)) + z2 * p1;
    %     p1Next = taoR * p1 + xiR / 2 * (pNext(N+1) + p(N+1)) ;
    end
    %% Set output from output position
    out(n) = p(outputPos);
    
    %% Energies
    kinEnergyTube(n) = 1/(2 * rho * c^2) * h * sum(SBar .* scaling .* p.^2);
    potEnergyTube(n) = rho / 2 * h * sum(SHalf .* vNext .* v);
    hTube(n) = kinEnergyTube(n) + potEnergyTube(n);
    kinEnergy(n) = M / 2 * (1/k * (y - yPrev))^2;
    potEnergy(n) = K / 4 * (y^2 + yPrev^2);
    hReed(n) = kinEnergy(n) + potEnergy(n);
    hColl(n) = psiPrev^2 / 2;
    if radiation
        hRad(n) = SBar(N+1) / 2 * (Lr * v1^2 + Cr * p1^2);

        v3Next = p1Next / R2;
        v3 = p1 / R2;
        pBar = 0.5 * (pNext(N+1) + p(N+1));
        muTPv2 = (pBar - 0.5 * (p1Next + p1)) / R1;
    else
        hRad(n) = 0;
    end
    % summed forms (damping and power input)
    idx = n - (1 * (n~=1));
    qReed(n) = R * (1/(2*k) * (yNext(n) - yPrev))^2 + Ub * deltaP;
    qHReed(n) = k * qReed(n) + qHReed(idx);
    pReed(n) = -(Ub + Ur) * Pm;
    pHReed(n) = k * pReed(n) + pHReed(idx);
    if radiation
        qRad(n) = SBar(N+1) * (R1 * muTPv2^2 + R2 * (0.5 * (v3Next + v3))^2);
        qHRad(n) = k * qRad(n) + qHRad(idx);
    else
        qRad(n) = 0;
        qHRad(n) = 0;
    end
%     qRad(n
%     qRad(n) = SBar(N+1) * (R1 / 2 * (v))

    % total energies
    totEnergy(n) = hTube(n) + hReed(n) + hColl(n) + hRad(n) + qHReed(idx) + pHReed(idx) + qHRad(idx);
%     scaledTotEnergy(n) = (totEnergy(n) - hTube(2) - hReed(2) - hColl(2) - hRad(2)) / 2^floor(log2(hTube(1) + hReed(1) + hColl(1) + hRad(1)));
    %% Draw things
    if drawThings && mod (n, drawSpeed) == 1
        if ~plotEnergy
%             scatter(0, y, 'b')
            pScaling = 0.0005;
            if plotSubplots
                subplot(2, 3, subp)
            end
            cla
            plot(0,0);
            hold on
            plotPressurePotential (p * pScaling, sqrt(S)/pi);
            plot (p*0.002*pScaling, 'k--', 'Linewidth', 1)
            plot( -sqrt(S)/pi, 'k', 'Linewidth', 1.5)
            plot(sqrt(S)/pi, 'k', 'Linewidth', 1.5)
            ylim(1.5*[-sqrt(S(1))/pi, sqrt(S(1))/pi])
%             plot([0 0], [-sqrt(S(1))/pi, -H0], 'k', 'Linewidth', 15)
            %% draw lips
            lipWidth = 20;
            
            % pressure colour between lips
            patch([-lipWidth, 1, 1, -lipWidth], ...
                [y, y, -H0, -H0], [clamp(p(1)*pScaling + 1, 0, 1), clamp(1-abs(p(1)*pScaling), 0, 1), clamp(1-p(1)*pScaling, 0, 1)], ...
                'EdgeColor', 'none')
            
            % bottom lip
            patch([-lipWidth, 1, 1, -lipWidth], ...
                [-H0, -H0, -1.2 * sqrt(S(1))/pi, -1.2 * sqrt(S(1))/pi], 'w', ...
                'EdgeColor', 'k', 'Linewidth', 1.5);
            text(-0.5*lipWidth, -H0 - 0.0003, "$-H_0$", 'horizontalAlignment',...
                'center', 'verticalAlignment', 'top', ...
                'interpreter', 'latex', ...
                'Fontsize', 16, 'color', 'b')
            plot([-lipWidth, 1], [-H0, -H0], 'color', 'b', 'Linewidth', 2)

            % top lip
            patch([-lipWidth, 1, 1, -lipWidth], ...
                [1.2 * sqrt(S(1))/pi, 1.2 * sqrt(S(1))/pi, y, y], 'w', ...
                'EdgeColor', 'k', 'Linewidth', 1.5)

            
            text(-0.5*lipWidth, y+ 0.0003, "$y^n$", 'horizontalAlignment',...
                'center', 'verticalAlignment', 'bottom', ...
                'interpreter', 'latex', ...
                'Fontsize', 16, 'color', [0, 0.85, 0])
            plot([-lipWidth, 1], [y, y], 'color', [0, 0.85, 0], 'Linewidth', 2)
            xticks([0, 50, 100, 150])
            xticklabels({'$0$','$50$' '$100$','$150$'})
            
            xlim([-lipWidth - 5, 150])
            xLab = xlabel("$l$", 'interpreter', 'latex', 'Fontsize', 16);
            yLim = ylim;
            xLim = xlim;
            xLab.Position(1) = 75;
            xLab.Position(2) = yLim(1) - 0.05 * (yLim(2) - yLim(1));
            title("$n = " + n + "$", 'interpreter', 'latex')
            ylabel("Displacement [m]")
%             grid on
            if plotSubplots

                if subp < 4
                    set(gca, 'Fontsize', 16, 'Linewidth', 1, ...
                        'Position', [0.0401+(subp-1)*inc 0.5838 0.2815 0.35], ...
                        'Fontname', 'times', 'TickLabelInterpreter', 'latex')
                else
                    set(gca, 'Fontsize', 16, 'Linewidth', 1, ...
                        'Position', [0.0401+(subp-4)*inc 0.0738 0.2815 0.35], ...
                        'Fontname', 'times', 'TickLabelInterpreter', 'latex')

                end
                subp = subp + 1;
                if subp > 6
                    return;
                end
            else
                set(gca, 'Fontsize', 16, 'Linewidth', 1, ...
                        'Fontname', 'times', 'TickLabelInterpreter', 'latex')
            end
%         hLocs = (0:length(p)-1) * h;
%         % Plot the velocity
% %         subplot(4,1,1)
%         cla
%         hold on;
% %         plotPressurePotential (p / 10000, sqrt(S));
% %         plot(p / 100000)
% %         plot(sqrt(S), 'k');
% %         plot(-sqrt(S), 'k');
%         plot(hLocs * N / L, p, '-o');
% %         xlim([1 N]);
% %         scatter(1, (y + H0) * 10000)
% %         ylim([-max(sqrt(S)) max(sqrt(S))] * 1.1);
% %         title("Pressure");
%         
% %         % Plot the velocity
% %         subplot(4,1,2)
% %         cla;
%         plot(hLocs(1:end-1) * N / L + 0.5, vNext * 100, 'Marker', '.', 'MarkerSize', 10);
% %         hold on;
% %         plot(sqrt(S) * amp, 'k');
% %         plot(-sqrt(S) * amp, 'k');
% %         xlim([1 N]);
% %         title("Particle Velocity")
% % pause(0.2)
        % Plot the output
%         subplot(4,1,3)
%         plot(out(1:n))
%         
%         % Plot scaled energy
%         subplot(4,1,4)
%         plot(scaledTotEnergy(2:n))
%         hold off;
%         plot(hTube(1:n))
%         hold on;
%         plot(hReed(1:n))
%         plot(hColl(1:n))
% hold off;
% plot(kinEnergyReed(1:n))
% hold on;
% plot(potEnergyReed(1:n))
        else
            plot(totEnergy(2:n)/totEnergy(2) - 1)
        end
        drawnow;

    end

    %% Update states
    v = vNext;
    p = pNext;
    if radiation
        p1 = p1Next;
        v1 = v1Next;
    end
    yPrev = y;
    y = yNext(n);
    psiPrev = psi; 

end   
plotTromboneEnergy;
function [S, SHalf, SBar] = setTube(N, NnonExtended, setToOnes)

    if N < NnonExtended
        S = 0.01 * ones(N,1);
        disp("N < NnonExtended")
    else
        lengths = [0.708, 0.177, 0.711, 0.241, 0.254, 0.502];
        radii = [0.0069, 0.0072, 0.0069, 0.0071, 0.0075, 0.0107]; % two radii for tuning slide

        lengthN = round(NnonExtended * lengths ./ sum(lengths));
        addPointsAt = round(lengthN(1) + lengthN(2) * 0.5) + (N-NnonExtended) * 0.5; % indicate split of two connected schemes (including offset if N differs from NnonExtended

        mouthPiece = 0.013 * (0.45 * (1 + cos(pi * ((1:lengthN(1))'-1) / (lengthN(1)-1))) + 0.1);
        inner1 = ones(lengthN(1), 1) * radii(1);
        inner2 = ones(lengthN(3), 1) * radii(3);
        gooseneck = ones(lengthN(4), 1) * radii(4);
        tuning = linspace(radii(5), radii(6), lengthN(5))';

        x0 = 0.0174; 
        b = 0.0063;
        flare = 0.7;
        bellL = lengthN(end);

        bell = b * ((lengths(6):-lengths(6) / (bellL - 1):0) + x0).^(-flare);


    %     pointsLeft = N - length([mp, m2t, bell]);
    %     tube = linspace(m2t(end), m2t(end), pointsLeft);    % tube
        totLengthN = length(inner1) + length(inner2) + length(gooseneck) + length(tuning) + length(bell);
    %     lengthN(2) = lengthN(2) + (N - NnonExtended + 1);
        slide = ones(N - totLengthN, 1) * radii(2);
        totRadii = [inner1; slide; inner2; gooseneck; tuning; bell'];

        % True geometry
        S = totRadii.^2 * pi;
    %     S = flipud(S);
    %     addPointsAt = N - addPointsAt;
        if setToOnes
    %         S = exp((-N-1:0)'/(0.5*N))/2;
    %         S = 0.5-0.5 *cos(0.5*pi:pi/(N-1):1.5*pi)';
    %         S = rand(size(S)) * 0.5 + 1;
            S = 5e-05 * ones(size(S));

        end
    end
    % Calculate approximations to the geometry
    SHalf = (S(1:end-1) + S(2:end)) * 0.5;                  % mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                        % mu_{x-}S_{l+1/2}
end

function [val] = clamp (input, min, max)
   if input < min
       val = min;
   elseif input > max
       val = max;
   else
       val = input;
   end
end