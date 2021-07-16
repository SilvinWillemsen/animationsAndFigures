%{
    Webster's equation
%}

clear all;
close all;

% drawing variables
drawThings = false;
drawSpeed = 1;

impulse = false;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = fs; % Duration (s)

L = 3;              % Length
c = 343;            % Wave speed (m/s)
h = c * k;          % Grid spacing (m)

N = floor(L/h);             % Number of points (-)
L = N * h;
h = L/N;                    % Recalculate gridspacing from number of points
lambdaSq = (c * k / h)^2    % Courant number

% Set cross-sectional geometry
LnonExtended = 2.593;
NnonExtended = LnonExtended / h;
[S, SHalf, SBar] = setTube (N+1, NnonExtended, false);

a1 = 1 / (4 * (0.8216)^2 * c);              % loss term
a2 = L / (0.8216 * sqrt(S(1)*S(N+1)/pi));     % inertia coefficient
% a1 = 0;
% a2 = 0;
%Initialise states
uNext = zeros(N+1, 1);
u = zeros(N+1, 1);
% a2 = 0
amp = 10000;
if ~impulse  
    % input signal
    t = (0:lengthSound - 1) / fs;
    freq = 213;
    in = cos(2 * pi * freq * t) - 0.5;
%     in = chirp(t, 200, t(end), 800);
    in = (in + abs(in)) / 2; % subplus
%     in = in - sum(in) / length(in);
    in = in * amp;
    rampLength = 1000; 
    env = [linspace(0, 1, rampLength), ones(1, lengthSound - rampLength)];
    in = in .* env;
%     in = in - sum(in) / length(in);
else
    in = zeros(lengthSound, 1);
    u(floor(N / 3) - 4 : floor(N / 3) + 4) = amp*hann(9);
%     u(1) = 1;
end
uPrev = u;

% output
out = zeros(lengthSound, 1);
outputPos = floor(1/5 * N);

% Initialise energy vectors
kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
boundaryEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

% Set ranges
range = 2:N;          % range without boundaries
potEnergyRange = 1:N; % range for the potential energy

% set up figure
figure(1);
subplot(2,1,1)
plot(sqrt(S), 'k');
hold on;
plot(-sqrt(S), 'k');
xlim([1 N+1]);
ylim([-max(sqrt(S)) max(sqrt(S))] * 1.1);
title("Pressure potential")

subplot(2,1,2)
title("Normalised energy (should be within machine precision)") 

% Problem 9.5 (weighted boundary conditions)
epsilonL = SHalf(1)/SBar(1);
epsilonR = SHalf(end)/SBar(end)

scaling = ones(N+1,1);
scaling(1) = epsilonL / 2;
scaling(end) = epsilonR / 2;

% virtual geometry points
SNph = 2 * SBar(end) - SHalf(end);
Smh = 2 * SBar(1) - SHalf(1);
qTot = 0;

% Initialise matrices
DxxS = sparse(1:N, 2:N+1, [2; SHalf(2:end) ./ SBar(2:end-1)], N+1, N+1) +...
       sparse(1:N+1, 1:N+1, -2 * ones(N+1, 1), N+1, N+1) + ...
        sparse(2:N+1, 1:N, [SHalf(1:end-1) ./ SBar(2:end-1); 2], N+1, N+1);

DxxS = DxxS / h^2;
alfPlus = h * (a1/k + a2) * lambdaSq * SNph / SBar(end);
alfMinus = h * (a1/k - a2) * lambdaSq * SNph / SBar(end);
B = 2 * speye(N+1) + c^2 * k^2 * DxxS;
Amat = speye(N+1);
Amat(end, end) = (1 + alfPlus);
C = -speye(N+1);
C(end, end) = -1 + alfMinus;
% 
% plotDampingAgainstFrequency = true;
% plotModalAnalysis;
% figure
% scatter(1:length(s)-1, (imag(s(2:end))-imag(s(1:end-1)))/(2*pi))
% drawnow
for n = 1:lengthSound
    
    vIn = 2 * h * lambdaSq * Smh / SBar(1) * in(n);
    % calculate scheme
    uNext = Amat \ B * u + Amat \ C * uPrev + Amat \ [vIn; zeros(N,1)];

%     uNext(1) = uNext(1) + 2 * h * lambdaSq * Smh / SBar(1) * in(n);
    % set output from output position
    out(n) = u(end);
    
    % energies
    kinEnergy(n) = 1/2 * sum(h * SBar .* scaling .* (1/k * (u - uPrev)).^2);
    potEnergy(n) = h * c^2 / 2 * 1/h^2 * sum(SHalf(potEnergyRange)...
        .* (u(potEnergyRange+1) - u(potEnergyRange)) .* (uPrev(potEnergyRange+1) - uPrev(potEnergyRange)));

    boundaryEnergy(n) = (2 - epsilonR) * SHalf(end) * c^2 * a2 / 4 * (u(end)^2 + uPrev(end)^2);
    
    q0(n) = c^2 * SHalf(end) * (2- epsilonR) * a1 * (1/(2*k) * (uNext(end) - uPrev(end)))^2;
    idx = n - (1 * (n~=1));
    qTot = qTot + k * (q0(idx));
    
    totEnergy(n) = kinEnergy(n) + potEnergy(n) + boundaryEnergy(n) + qTot;
    
    %% Rate-of-Changes of energy
    rOCkinEnergy(n) = h / (2*k^3) * sum(SBar .* scaling .* (uNext - 2 * u + uPrev) .* (uNext - uPrev));
    rOCpotEnergy(n) = -c^2 / (2 * k * h) * sum(SHalf .* (uNext(potEnergyRange+1) - uNext(potEnergyRange) - uPrev(potEnergyRange+1) + uPrev(potEnergyRange)) .* (u(potEnergyRange+1) - u(potEnergyRange)));
    
    rOCboundaryEnergy(n) = -(2-epsilonL) * SHalf(1) * c^2 * 1/(2*k) * (uNext(1) - uPrev(1)) * (-in(n));
    rOCboundaryEnergy(n) = rOCboundaryEnergy(n) + (2-epsilonR) * SHalf(end) * c^2 * (-a1 * (1/(2*k) * (uNext(end) - uPrev(end)))^2 - a2/(4*k) * (uNext(end)^2 - uPrev(end)^2));

    rOCtotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCboundaryEnergy(n);
    
    % draw things
    if drawThings && mod (n, drawSpeed) == 0
        
        % plot state
        subplot(3,1,1)
        cla
        hold on;
        plot (uNext * 1/amp);
        plot(sqrt(S) * amp, 'k');
        plot(-sqrt(S) * amp, 'k');
        xlim([1 N+1]);
%         ylim([-max(sqrt(S)) max(sqrt(S))] * amp * 1.1);
        title("Pressure potential. n = " + num2str(n))

        % plot total energy
        subplot(3,1,2)
        plot(uVec - u)
        
%         if n > 10
%             if impulse% && a1 == 0
%                 plot(totEnergy(10:n) / totEnergy(10) - 1);
%                 title("Normalised energy (should be within machine precision)") 
%             else
%                 plot(totEnergy(10:n))
%                 title("Total energy") 
%             end
%         end
        
        subplot(3,1,3)
        hold off
        plot(kinEnergy(1:n))
        hold on;
        plot(potEnergy(1:n))
        plot(boundaryEnergy(1:n))

        title("Rate-of-change of energy (should be very close to 0)");
        drawnow;
    end
   
    % update states    
    uPrev = u;
    u = uNext;
end   
plotEnergyWebster;


%% Plotting
figure('Position', [173 578 827 220])

t = 1:lengthSound;
% plot(t, zeros(length(t), 1), '--', 'Linewidth', 2, 'color', [0.5, 0.5, 0.5])
% hold on;
subp1 = subplot(1, 2, 1);
plot(t, out, 'k', 'Linewidth', 2)

xlim([0, 5000])
% ylim([-1.1, 1.1])
xLab = xlabel("$n$", 'interpreter', 'latex');
yLab = ylabel("$u^n_N$", 'interpreter', 'latex');

yLim = ylim;
ylim([yLim * 1.1])
yLim = yLim;

xLab.Position(2) = yLim(1) - (yLim(2) - yLim(1)) * 0.18;
xLim = xlim;
yLab.Position(1)  = xLim(1) -0.07 * (xLim(2) - xLim(1));
% grid on;

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
xLim = xlim;
xLab2.Position(2) = yLim(1) -0.125 * (yLim(2) - yLim(1));
yLab2.Position(1) = xLim(1) -0.09 * (xLim(2) - xLim(1));

% yLab2.Position(1) = -70;

set(gcf, 'color', 'w')

    
%% Plotting
figure('Position', [173 578 422 220])
subp1 = subplot(1, 2, 1);
plot(sqrt(S/pi), 'k', 'Linewidth', 2)
hold on;
plot(-sqrt(S/pi), 'k', 'Linewidth', 2)
xlim([1, N+1])
ylabel("Cross-section [m]", 'Fontname', 'times')
xLab = xlabel("$x$ [m]", 'interpreter', 'latex')
yLim = ylim;
ylim([yLim * 1.1])
yLim = yLim;

xLab.Position(2) = yLim(1) - (yLim(2) - yLim(1)) * 0.12;

set(gca,'xtick',[1, N+1],'xticklabel', ["$0$", "$L$"], ...
        'ticklabelinterpreter', 'latex', 'FontSize', 16, ...
        'Linewidth', 1.5, 'Position', [0.1398 0.1500 0.8389 0.8091])

grid off
  
figure('Position', [173 578 786 220])
plot(0:1/fs:1-1/fs, in, 'k', 'Linewidth', 2)

xlim([0, 0.1])
yLim = ylim;
ylim([yLim * 1.1])

xLab2 = xlabel("$t$", 'interpreter', 'latex');
yLab2 = ylabel("$v_\textrm{\fontsize{7}{7}\selectfont in}$", 'interpreter', 'latex', 'Fontsize', 20);

yLim = ylim;
xLab2.Position(2) = yLim(1) -0.07 * (yLim(2) - yLim(1));

yticks([])
xticks(0:0.02:0.1)
set(gca, 'Linewidth', 1.5, 'Fontsize', 16,...
    'Position', [0.0356 0.1500 0.9338 0.8091], ...
    'TickLabelInterpreter', 'latex')
set(gcf, 'color', 'w')

    
% figure
% plot(abs(fft(out)))
% function [S, SHalf, SBar] = setTube(N)
%     mp = linspace(0.0005, 0.0005, floor(N/20));       % mouthpiece
%     m2t = linspace(0.0005, 0.0001, floor(N/20));     % mouthpiece to tube
%     alpha = 0.25;
%     b = m2t(end) * exp(alpha * (0:18));               % bell
%     pointsLeft = N - length([mp, m2t, b]);
%     tube = linspace(m2t(end), m2t(end), pointsLeft);        % tube
% 
%     S = [mp, m2t, tube, b]';                        % True geometry
%     S = 0.001 * ones(length(S), 1);
%     % Calculate approximations to the geometry
%     SHalf = (S(1:N-1) + S(2:N)) * 0.5;           	% mu_{x+}
%     SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
%     SBar = [S(1); SBar; S(end)];                    % mu_{x-}S_{l+1/2}
% 
% end

function [S, SHalf, SBar] = setTube(N, NnonExtended, setToOnes)

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
        S = 0.5-0.5 *cos(0.5*pi:pi/(N-1):1.5*pi)';
%         S = rand(size(S)) * 0.5 + 1;
%         S = ones(size(S));

    end
    
    % Calculate approximations to the geometry
    SHalf = (S(1:end-1) + S(2:end)) * 0.5;                  % mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                        % mu_{x-}S_{l+1/2}
end