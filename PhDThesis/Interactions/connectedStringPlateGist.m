%{
    Implementation of a string connected to a plate through a nonlinear
    damped spring. This implementation accompanies the PhD Thesis: "The
    Emulated Ensemble" by Silvin Willemsen (more specifically Section 10.5)

    CC 3.0 Silvin Willemsen 2021.
%}

close all;
clear all;

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
lengthSound = fs;   % Length of the simulation (1 second) [samples]             

%% Plotting
drawThings = true;     % To plot or not
drawSpeed = 1;         % How fast to plot
if drawThings
    figure('Position', [180 454 820 344])
end

%% String variables
LS = 0.75;                      % Length [m]
rhoS = 7850;                    % Material density [kg/m^3]
rS = 0.0005;                    % Radius [m]
A = pi * rS^2;                  % Cross-sectional area [m^2]
ES = 2e11;                      % Young's modulus [Pa]
IS = pi * rS^4 / 4;             % Area moment of inertia [m^4]
T = 555;                        % Tension [N]
sig0S = 1;                      % Frequency-independent damping [1/s]
sig1S = 0.005;                  % Frequency-dependent damping [m^2/s]

cS = sqrt(T / (rhoS * A));        % Wave speed 
kappaSqS = ES * IS / (rhoS * A);   % Stiffness term

% Calculate grid spacing and number of points
hS = sqrt((cS^2 * k^2 + 4 * sig1S * k ...
    + sqrt((cS^2 * k^2 + 4 * sig1S * k)^2 + 16 * kappaSqS * k^2))/2);
NS = floor(LS / hS);
hS = LS / NS;

%% Plate variables
Lx = 0.75;                  % Length in x-direction [m]
Ly = 0.5;                   % Length in y-direction [m]
rhoP = 7850;                % Material density [kg/m^3]
EP = 2e11;                  % Young's modulus [Pa]
HP = 0.0005;                % Thickness [m]
nu = 0.3;                   % Poisson's ratio [-]
sig0P = 1.0;                % Frequency-independent damping [1/s]
sig1P = 0.005;              % Frequency-ddependent damping [m^2/s]

Dvar = EP * HP^3 / (12 * (1-nu^2)); % Stiffness term [Joules]
kappaSqP = Dvar / (rhoP * HP);      % Stiffness coefficient

% Calculate grid spacing and number of points
hP = 2 * sqrt(k * (sig1P + sqrt(kappaSqP + sig1P^2)));
Nx = floor(Lx / hP);
Ny = floor(Ly / hP);
hP = min(Lx/Nx, Ly/Ny);

%% Connection Variables
K1 = 1e4;
K3 = 1e7;
R = 10;

% Connection locations 
connLocS = 0.25;    % along string
connLocPX = 0.25;   % along x-direction plate
connLocPY = 0.75;   % along y-direction plate
    

%% Initialise string matrices (simply supported boundaries)
IdS  = speye(NS-1);  % identity matrix

% Dxx and Dxxxx matrices
DxxS = toeplitz([-2, 1, zeros(1, NS-3)]) / hS^2;
DxxxxS = DxxS * DxxS;

% Matrices used for update equation
AS = (1 + sig0S * k);
BS = 2 * IdS + cS^2 * k^2 * DxxS - kappaSqS * k^2 * DxxxxS + 2 * sig1S * k * DxxS;
CS = -(1 - sig0S * k) * IdS - 2 * sig1S * k * DxxS;

% Include the division before main loop for faster computation
BSoverA = BS / AS;
CSoverA = CS / AS;


%% Initialise plate matrices (simply supported boundaries)
Nxw = Nx - 1;
Nyw = Ny - 1;

% Kronecker sum for stacked matrix form
DxxP = toeplitz([-2, 1, zeros(1, Nxw-2)]);
DyyP = toeplitz([-2, 1, zeros(1, Nyw-2)]);
    
D = kron(speye(Nxw), DyyP) + kron(DxxP, speye(Nyw));
D = D / hP^2;
DD = D*D;

% Total number of points
Nw = Nxw * Nyw;

% Matrices used for update equation
AP = (1 + sig0P * k);
BP = 2 * speye(Nw) - kappaSqP * k^2 * DD + 2 * sig1P * k * D;
CP = -(1-sig0P * k) * speye(Nw) - 2 * sig1P * k * D;

% Include the division before main loop for faster computation
BPoverA = BP / AP;
CPoverA = CP / AP;

%% Matrices used for energy analysis

% Dx+ (string)
Dxp = sparse(1:NS, 1:NS, -ones(1, NS), NS, NS) + ...
        sparse(1:NS-1, 2:NS, ones(1, NS-1), NS, NS);
Dxp = Dxp / hS;

% Discrete Laplacian (plate)
DxxEnP = toeplitz([-2, 1, zeros(1, Nxw)]);
DyyEnP = toeplitz([-2, 1, zeros(1, Nyw)]);
DEnPP = kron(speye(Nxw+2), DyyEnP) + kron(DxxEnP, speye(Nyw+2));
DEnPP = DEnPP / hP^2;

%% Excitation (raised cosine)
amp = 5;                    % amplitude
offset = amp * 0.5;         % plot offset
halfWidth = 5;              % half of the excitation width (minus 1)
loc = 0.33333;              % location of the raised cosine
width = halfWidth*2 + 1;    % full width
startLoc = floor(NS * loc) - halfWidth + 1; % start location of the exciation

%% Initialise states (simply supported boundary conditions)
% string
u = zeros(NS-1, 1);
u(startLoc:startLoc+width-1) = u(startLoc:startLoc+width-1) + amp * hann(width);
uPrev = u; % only set initial displacement, not initial velocity

% plate
w = zeros(Nw, 1); 
wPrev = w;

%% Initialise (linear) interpolation and spreading operators

% string
xcu = floor(connLocS * NS);
alphaU = connLocS * NS - xcu;
xcu = xcu + 1;          % +1 because of 1-based matlab

Iu = zeros(1, NS-1);    % interpolation operator
Iu(xcu) = 1 - alphaU;
Iu(xcu+1) = alphaU;

Ju = 1/hS * Iu';        % spreading operator

% plate
xcwX = floor(connLocPX * Nxw);
alphaWX = connLocPX * Nxw - xcwX;
xcwX = xcwX + 1;        % +1 because of 1-based matlab
xcwY = floor(connLocPY * Nyw);
alphaWY = connLocPY * Nyw - xcwY;
xcwY = xcwY + 1;        % +1 because of 1-based matlab

Imat = zeros(Nyw, Nxw); % interpolation operator
Imat(xcwY, xcwX) = (1 - alphaWX) * (1 - alphaWY);
Imat(xcwY, xcwX+1) = alphaWX * (1 - alphaWY);
Imat(xcwY+1, xcwX) = (1 - alphaWX) * alphaWY;
Imat(xcwY+1, xcwX+1) = alphaWX * alphaWY;

Iw = reshape(Imat, 1, Nw);
Jw = 1/hP^2 * Iw';      % spreading operator

% Initialise output vector
out = zeros(lengthSound, 1);

%% Initialisation of vectors storing energy values (for quicker computation)
% string 
kinEnergyU = zeros(lengthSound, 1);
potEnergyU = zeros(lengthSound, 1);
totEnergyU = zeros(lengthSound, 1);

q0U = zeros(lengthSound, 1);
q1U = zeros(lengthSound, 1);
qTotU = 0;

% plate 
kinEnergyW = zeros(lengthSound, 1);
potEnergyW = zeros(lengthSound, 1);
totEnergyW = zeros(lengthSound, 1);

q0W = zeros(lengthSound, 1);
q1W = zeros(lengthSound, 1);
qTotW = 0;

% connection 
connEnergy = zeros(lengthSound, 1);
totEnergyC = zeros(lengthSound, 1);

qC = zeros(lengthSound, 1);
qTotC = 0;

totEnergy = zeros(lengthSound, 1);

%% Main Loop
for n = 1:lengthSound
    
    %% Calculate [Int]ermediate system states without force term 
    uInt = BSoverA * u + CSoverA * uPrev;
    wInt = BPoverA * w + CPoverA * wPrev;
    
    %% Retrieve connection locations without force terms
    uStar = Iu * uInt;
    wStar = Iw * wInt;

    % Relative displacement at n and n-1
    eta = Iu * u - Iw * w;
    etaPrev = Iu * uPrev - Iw * wPrev;

    %% Calculate force
    rPlus = K1 / 4 + K3 * eta^2 / 2 + R / (2*k);
    rMinus = K1 / 4 + K3 * eta^2 / 2 - R / (2*k);
    
    f = (uStar - wStar + K1 / (2 * rPlus) * eta + rMinus / rPlus * etaPrev) ...
            / (1/rPlus + Iu * Ju * k^2 / (rhoS * A * (1+sig0S * k)) + Iw * Jw * k^2 / (rhoP * HP * (1+sig0P * k)));
    
    %% Add forces to scheme
    uNext = uInt - Ju * k^2 * f / (rhoS * A * (1 + sig0S * k));
    wNext = wInt + Jw * k^2 * f / (rhoP * HP * (1 + sig0P * k));
       
    % Retrieve output
    out(n) = Iu * uNext;
      
    %% Calculate energy
    
    % index variable used for the damping terms
    idx = n - (1 * (n~=1));
    
    % eta^{n+1} is used to calculate the (energetic) damping term of the connection 
    etaNext = Iu * uNext - Iw * wNext;

    % String
    kinEnergyU(n) = rhoS * A * hS / 2 * sum((1/k * (u - uPrev)).^2);
    potEnergyU(n) = T / 2 * hS * sum((Dxp * [0; u]) .* (Dxp * [0; uPrev])) ...
        + ES * IS * hS / 2 * sum((DxxS * u) .* (DxxS * uPrev));

    % String damping 
    q0U(n) = 2 * sig0S * rhoS * A * hS * sum((1/(2*k) * (uNext - uPrev)).^2);
    q1U(n) = - 2 * sig1S * rhoS * A * hS * sum(1/(2*k) * (uNext - uPrev)...
        .* (1/k * (DxxS * u - DxxS * uPrev)));
    qTotU = qTotU + k * (q0U(idx) + q1U(idx));

    totEnergyU(n) = kinEnergyU(n) + potEnergyU(n) + qTotU;

    %% Plate
    kinEnergyW(n) = rhoP * HP / 2 * hP^2 * sum((1/k * (w-wPrev)).^2);
    potEnergyW(n) = Dvar * hP^2 / 2 * sum((D * w) .* (D * wPrev));
    
    % Plate damping 
    q0W(n) = 2 * sig0P * rhoP * HP * hP^2 * sum((1/(2*k) * (wNext - wPrev)).^2);
    q1W(n) = -2 * sig1P * rhoP * HP * hP^2 * sum(1/(2*k) * (wNext - wPrev)...
        .* (1/k * (D * w - D * wPrev)));
    idx = n - (1 * (n~=1));
    qTotW = qTotW + k * (q0W(idx) + q1W(idx));

    totEnergyW(n) = kinEnergyW(n) + potEnergyW(n) + qTotW;

    %% Connection
    connEnergy(n) = K1 / 8 * (eta + etaPrev)^2 + K3 / 4 * (eta^2 * etaPrev^2);
    
    % Connection damping
    qC(n) = R * (1/(2*k) * (etaNext - etaPrev))^2;
    qTotC = qTotC + k * qC(idx);
    
    totEnergyC(n) = connEnergy(n) + qTotC;
    
    %% Total Energy
    totEnergy(n) = totEnergyU(n) + totEnergyW(n) + totEnergyC(n);
    
    %% Plot system in 3D
    if drawThings && mod(n, drawSpeed) == 0
        
        % system states
        hold off;
        
        %string 
        u_s = plot3([0:NS] + xcwX +1 - xcu, xcwY * ones(NS+1), [0;u;0] + offset, 'r', 'Linewidth', 2);
        hold on;
        
        % plate
        mesh(reshape(w, Nyw, Nxw))
        colormap gray;
        caxis([-0.05, 0.05]);
        
        % connection
        plot3([xcwX, xcwX]+1, [xcwY, xcwY], [Iu * u + offset, Iw * w], 'color', [0, 0.85, 0], 'Linewidth', 2)

        % plot settings
        zlim([-0.2, amp + offset])
        xticks([])
        yticks([])
        zticks([])
        title("$n = " + n + "$", 'interpreter', 'latex', 'Fontsize', 16)
        view (20, 35);
        box on
        grid on
        
        % draw
        drawnow;
    end

    % Update States
    uPrev = u;
    u = uNext;
        
    wPrev = w;
    w = wNext;
    
end

%% Plot normalised energy (should be within machine precision)
plot(totEnergy(1:n) / totEnergy(1) - 1)