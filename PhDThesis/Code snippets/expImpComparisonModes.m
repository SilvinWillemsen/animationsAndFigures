

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]

% Material properties and geometry
L = 1;              % Length [m]
r = 5e-4;           % Radius [m]
A = pi * r^2;       % Cross-sectional area [m^2] (circular cross-section)
rho = 7850;         % Material density [kg / m^3]
E = 2e14;           % Young's modulus [Pa]
I = pi * r^4 / 4;   % Area moment of inertia [m^4]
T = 1885;            % Tension [N]

% Damping coefficients
sig0 = 1;         % Frequency-independent damping [s^{-1}]
sig1 = 1;       % Frequency-dependent damping [m^2/s]

% Scheme coefficients
c = sqrt(T / (rho * A));            % Wave speed [m/s]
kappa = sqrt(E * I / (rho * A));    % Stiffness coefficient [m^2/s]

%% Explicit matrices

% Grid spacing and number of intervals
h = sqrt(1/2 * (c^2*k^2 + 4*sig1*k ...
    + sqrt((c^2*k^2 + 4*sig1*k)^2 + 16*kappa^2*k^2)));
N = floor(L/h);     % Number of intervals between grid points
h = L / N;          % Recalculation of grid spacing based on integer N

% Update coefficients
lambdaSq = c^2 * k^2 / h^2;
muSq = kappa^2 * k^2 / h^4;
 lambdaSq + 4 * muSq
%     tSave(tSaveIdx) = lambdaSq + 4 * muSq;
%     tSaveIdx = tSaveIdx + 1;
% end
% Boundary conditions ([c]lamped, [s]imply supported or [f]ree)
bc = "s";            

% Change N to the usable range
Norig = N;
if bc == "c"
    N = N - 4;
elseif bc == "s"
    N = N - 2;
elseif bc == "f"
    N = N;
end

%% Initialise scheme matrices
Id  = eye(N+1);         % identity matrix

Dxx = toeplitz([-2, 1, zeros(1, N-1)]) / h^2;
Dxxxx = Dxx * Dxx;

if bc == "c"
    Dxxxx(1, 1) = 6 / h^4;
    Dxxxx(end,end) = 6 / h^4;
elseif bc == "f"
    Dxx(1, 2) = 2 / h^2;
    Dxx(end, end-1) = 2 / h^2;
    
    Dxxxx(2, 1:4) = [-2, 5, -4, 1] / h^4;
    Dxxxx(1,1:3) = [2, -4, 2] / h^4;
    Dxxxx(end-1, end-3:end) = [1, -4, 5, -2] / h^4;
    Dxxxx(end,end-2:end) = [2, -4, 2] / h^4;
end

AmatExp = (1 + sig0 * k);
BExp = 2 * Id + c^2 * k^2 * Dxx - kappa^2 * k^2 * Dxxxx + 2 * sig1 * k * Dxx;
CExp = -(1 - sig0 * k) * Id - 2 * sig1 * k * Dxx;

%% Implicit matrices

% Grid spacing and number of intervals
h = sqrt(1/2 * (c^2*k^2 + sqrt(c^4*k^4  + 16*kappa^2*k^2)));
N = floor(L/h);     % Number of intervals between grid points
h = L / N;          % Recalculation of grid spacing based on integer N

% Update coefficients
lambdaSq = c^2 * k^2 / h^2;
muSq = kappa^2 * k^2 / h^4;
lambdaSq + 4 * muSq
%     tSave(tSaveIdx) = lambdaSq + 4 * muSq;
%     tSaveIdx = tSaveIdx + 1;
% end
% Boundary conditions ([c]lamped, [s]imply supported or [f]ree)
bc = "s";            

% Change N to the usable range
Norig = N;
if bc == "c"
    N = N - 4;
elseif bc == "s"
    N = N - 2;
elseif bc == "f"
    N = N;
end

%% Initialise scheme matrices
Id  = eye(N+1);         % identity matrix

Dxx = toeplitz([-2, 1, zeros(1, N-1)]) / h^2;
Dxxxx = Dxx * Dxx;

if bc == "c"
    Dxxxx(1, 1) = 6 / h^4;
    Dxxxx(end,end) = 6 / h^4;
elseif bc == "f"
    Dxx(1, 2) = 2 / h^2;
    Dxx(end, end-1) = 2 / h^2;
    
    Dxxxx(2, 1:4) = [-2, 5, -4, 1] / h^4;
    Dxxxx(1,1:3) = [2, -4, 2] / h^4;
    Dxxxx(end-1, end-3:end) = [1, -4, 5, -2] / h^4;
    Dxxxx(end,end-2:end) = [2, -4, 2] / h^4;
end
AmatImp = (1 + sig0 * k) * Id -  sig1 * k * Dxx;
BImp = 2 * Id + c^2 * k^2 * Dxx - kappa^2 * k^2 * Dxxxx;
CImp = -(1 - sig0 * k) * Id - sig1 * k * Dxx;

%% PLOT MODAL ANALYSIS
figure('Position', [173 578 827 220])

% create Q matrix (one-step form)
QExp = [AmatExp \ BExp, AmatExp \ CExp;
     eye(size(BExp)), zeros(size(BExp))];
 
% obtain complex frequencies
sExp = 1/k * log(eig(QExp));

% obtain positive frequencies and sort them
sExp = sExp(imag(sExp) >= 0);
[~, order] = sort(imag(sExp));
sExp = sExp(order);

% create Q matrix (one-step form)
QImp = [AmatImp \ BImp, AmatImp \ CImp;
     eye(size(BImp)), zeros(size(BImp))];
 
% obtain complex frequencies
sImp = 1/k * log(eig(QImp));

% obtain positive frequencies and sort them
sImp = sImp(imag(sImp) >= 0);
[~, order] = sort(imag(sImp));
sImp = sImp(order);


p = 1:length(sImp);
fp = c * p / 2 .* sqrt(1 + (kappa^2 * pi^2) / c^2 * p.^2)


subplot(121)
plot(0, 0)
hold on
scatter(1:length(sExp), imag(sExp)/(2*pi), 'b', 'Linewidth', 2)
scatter(1:length(sImp), imag(sImp)/(2*pi), 'r', 'Linewidth', 2)
plot(fp)
grid on
xlabel("Mode number $p$", 'interpreter', 'latex')
ylabel("$f_p$ [Hz]", 'interpreter', 'latex')
xlim([1, length(sImp)])
title("Modal frequency")
set(gca, 'Fontsize', 16, 'Linewidth', 2, 'FontName', 'times', ...
    'Position', [0.0632 0.2091 0.4115 0.6909])

subplot(122)
plot(0, 0)
hold on
scatter(1:length(sExp), real(sExp), 'b', 'Linewidth', 2)
scatter(1:length(sImp), real(sImp), 'r', 'Linewidth', 2)

grid on
xlim([1, length(sImp)])
title("Damping per mode")
xlabel("Mode number $p$", 'interpreter', 'latex')
ylabel("$\sigma_p$ [s$^{-1}$]", 'interpreter', 'latex')
legend('', 'Explicit', 'Implicit')
set(gca, 'Fontsize', 16, 'Linewidth', 2, 'FontName', 'times', ...
    'Position', [0.5732 0.2091 0.4115 0.6909])
