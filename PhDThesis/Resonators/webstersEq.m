%{
    Webster's equation
%}

clear all;
close all;

% drawing variables
drawThings = true;
drawSpeed = 1;
centered = true;

impulse = true;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = 2000; % Duration (s)

c = 441;            % Wave speed (m/s)
h = c * k;          % Grid spacing (m)
L = 1;              % Length

N = floor(L/h);             % Number of points (-)
h = L/N;                    % Recalculate gridspacing from number of points
lambdaSq = (c * k / h)^2    % Courant number

% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N+1);

a1 = 1 / (2 * (0.8216)^2 * c);              % loss term
a2 = L / (0.8216 * sqrt(S(1)*S(N+1)/pi));     % inertia coefficient
% a1 = 0;
% a2 = 0;

%Initialise states
uNextVec = zeros(N+1, 1);
uVec = zeros(N+1, 1);

uNext = zeros(N+1, 1);
u = zeros(N+1, 1);

amp = 1e-5;
if ~impulse  
    % input signal
    t = (0:lengthSound - 1) / fs;
    freq = 446/2;
    in = cos(2 * pi * freq * t) - 0.5;
    in = (in + abs(in)) / 2; % subplus
    in = in - sum(in) / length(in);
    in = in * amp;
    rampLength = 0; 
    env = [linspace(0, 1, rampLength), ones(1, lengthSound - rampLength)];
    in = in .* env;
    in = in - sum(in) / length(in);
else
    in = zeros(lengthSound, 1);
    uVec(floor(N / 3) - 4 : floor(N / 3) + 4) = amp*hann(9);
    u(floor(N / 3) - 4 : floor(N / 3) + 4) = amp*hann(9);
end
uPrevVec = uVec;
uPrev = u;

% output
out = zeros(lengthSound, 1);
outputPos = floor(1/5 * N);

% Initialise energy vectors
kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
boundaryEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

rOCkinEnergy = zeros(lengthSound, 1);
rOCpotEnergy = zeros(lengthSound, 1);
rOCtotEnergy = zeros(lengthSound, 1);
rOCboundaryEnergy = zeros(lengthSound, 1);

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
epsilonR = SHalf(end)/SBar(end);

scaling = ones(N+1,1);
if centered
    scaling(1) = epsilonL / 2;
    scaling(end) = epsilonR / 2;
end

SNph = 2 * SBar(end) - SHalf(end);
Smh = 2 * SBar(1) - SHalf(1);
qTot = 0;

Dxx = toeplitz([-2, 1, zeros(1, N-1)]);
Dxx(end, end-1) = 2;

Dxx = Dxx;

test = sparse(1:10, 2:11, 1:10, 11, 11) + sparse(2:11, 1:10, 1:10, 11, 11)


Smat = sparse(1:N, 2:N+1, SHalf(1:end) ./ SBar(1:end-1), N+1, N+1) +...
       sparse(1:N+1, 1:N+1, ones(N+1, 1), N+1, N+1) + ...
        sparse(2:N+1, 1:N, [SHalf(1:end)] ./ SBar(2:end), N+1, N+1);

DxxSmat = (Dxx .* Smat);
DxxSmat(1, 2) = 2;
DxxSmat(end, end-1) = 2;
B = 2 * speye(N+1) + c^2 * k^2 / h^2 * DxxSmat;
Amat = speye(N+1);
Amat(end, end) = (1 + lambdaSq * SNph / SBar(end) * h * (a1/k + a2));
C = -speye(N+1);
C(end, end) = -1 + h * lambdaSq * SNph / SBar(end) * (a1/k - a2);

for n = 1:lengthSound
    
    % calculate scheme
    uNextVec(range) = 2 * (1 - lambdaSq) * uVec(range) - uPrevVec(range) + lambdaSq * ((SHalf(range) ./ SBar(range)) .* uVec(range+1) + (SHalf(range - 1) ./ SBar(range)) .* uVec(range-1));
    uNext = Amat \ B * u + Amat \ C * uPrev;
    
    uNextVec(1) = 2 * (1 - lambdaSq) * uVec(1) - uPrevVec(1) + lambdaSq * 2 * uVec(2) + 2 * h * lambdaSq * Smh / SBar(1) * in(n);
    uNextVec(end) = (2 * (1 - lambdaSq) * uVec(end) - uPrevVec(end) + lambdaSq * 2 * uVec(end-1) + h * lambdaSq * SNph / SBar(end) * (a1/k - a2) * uPrevVec(end)) / (1 + lambdaSq * SNph / SBar(end) * h * (a1/k + a2));

%     plot(uNextVec-uNext)
%     drawnow
    % set output from output position
    out(n) = uNextVec(outputPos);
    
    % energies
    kinEnergy(n) = 1/2 * sum(h * SBar .* scaling .* (1/k * (uVec - uPrevVec)).^2);
    potEnergy(n) = h * c^2 / 2 * 1/h^2 * sum(SHalf(potEnergyRange)...
        .* (uVec(potEnergyRange+1) - uVec(potEnergyRange)) .* (uPrevVec(potEnergyRange+1) - uPrevVec(potEnergyRange)));

    if centered
        boundaryEnergy(n) = (2 - epsilonR) * SHalf(end) * c^2 * a2 / 4 * (uVec(end)^2 + uPrevVec(end)^2);
    else
        boundaryEnergy(n) = c^2 * SNph * a2 / 4 * (uVec(end)^2 + uPrevVec(end)^2);
    end
    q0(n) = c^2 * SHalf(end) * (2- epsilonR) * a1 * (1/(2*k) * (uNextVec(end) - uPrevVec(end)))^2;
    idx = n - (1 * (n~=1));
    qTot = qTot + k * (q0(idx));
    
    totEnergy(n) = kinEnergy(n) + potEnergy(n) + boundaryEnergy(n) + qTot;
    
    %% Rate-of-Changes of energy
    rOCkinEnergy(n) = h / (2*k^3) * sum(SBar .* scaling .* (uNextVec - 2 * uVec + uPrevVec) .* (uNextVec - uPrevVec));
    rOCpotEnergy(n) = -c^2 / (2 * k * h) * sum(SHalf .* (uNextVec(potEnergyRange+1) - uNextVec(potEnergyRange) - uPrevVec(potEnergyRange+1) + uPrevVec(potEnergyRange)) .* (uVec(potEnergyRange+1) - uVec(potEnergyRange)));
    
    if centered
        rOCboundaryEnergy(n) = -(2-epsilonL) * SHalf(1) * c^2 * 1/(2*k) * (uNextVec(1) - uPrevVec(1)) * (-in(n));
        rOCboundaryEnergy(n) = rOCboundaryEnergy(n) + (2-epsilonR) * SHalf(end) * c^2 * (-a1 * (1/(2*k) * (uNextVec(end) - uPrevVec(end)))^2 - a2/(4*k) * (uNextVec(end)^2 - uPrevVec(end)^2));
    else
        % no boundaryEnergy at the left boundary as u(1) = u(2) -> ... * (u(1)-u(2)) = 0
        rOCboundaryEnergy(n) = c^2 * SNph * (-a1 * (1/(2*k) * (uNextVec(end) - uPrevVec(end)))^2 - a2 / (4*k) * (uNextVec(end)^2 - uPrevVec(end)^2));
    end
    
    rOCtotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCboundaryEnergy(n);
    
    % draw things
    if drawThings && mod (n, drawSpeed) == 0
        
        % plot state
        subplot(3,1,1)
        cla
        hold on;
        plot (uNextVec * 1/amp);
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
    uPrevVec = uVec;
    uVec = uNextVec;
    
    uPrev = u;
    u = uNext;
end   
plotEnergyWebster;

function [S, SHalf, SBar] = setTube(N)
    mp = linspace(0.0005, 0.0005, floor(N/20));       % mouthpiece
    m2t = linspace(0.0005, 0.0001, floor(N/20));     % mouthpiece to tube
    alpha = 0.25;
    b = m2t(end) * exp(alpha * (0:18));               % bell
    pointsLeft = N - length([mp, m2t, b]);
    tube = linspace(m2t(end), m2t(end), pointsLeft);        % tube

    S = [mp, m2t, tube, b]';                        % True geometry
%     S = 0.001 * ones(length(S), 1);
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;           	% mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                    % mu_{x-}S_{l+1/2}

end
