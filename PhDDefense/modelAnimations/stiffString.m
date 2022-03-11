%{
    Implementation of the damped stiff string including an energy analysis. 
    This implementation accompanies the PhD Thesis: "The Emulated Ensemble" 
    by Silvin Willemsen (more specifically Chapter 4).
    CC 3.0 Silvin Willemsen 2021.
%}

close all;
clear all;

%% Plotting options
drawThings = false;
drawEnergy = false; % Either plot the normalised energy or the string state
drawSpeed = 1;
if drawThings
    recordVid = true;
else
    recordVid = false;
end

%% Initialise variables
fs = 44100;         % Sample rate [Hz]
k = 1 / fs;         % Time step [s]
if recordVid
    lengthSound = 400 * fs / 44100;   % Length of the simulation (1 second) [samples]             
else
    lengthSound = fs;
end

if drawThings
    figure('Position', [489 644 560 213])
    if recordVid
        slowdown = 1;
        loops = lengthSound * slowdown;
        M(loops) = struct('cdata',[],'colormap',[]);
        frame = 1;
    end
end

% Material properties and geometry
L = 1;              % Length [m]
r = 5e-4;           % Radius [m]
A = pi * r^2;       % Cross-sectional area [m^2] (circular cross-section)
rho = 7850;         % Material density [kg / m^3]
E = 2e11;           % Young's modulus [Pa]
I = pi * r^4 / 4;   % Area moment of inertia [m^4]
T = 300;           % Tension [N]

% Damping coefficients
sig0 = 2;         % Frequency-independent damping [s^{-1}]
sig1 = 0.005;       % Frequency-dependent damping [m^2/s]

% Scheme coefficients
c = sqrt(T / (rho * A));            % Wave speed [m/s]
c = 220.5;
kappa = sqrt(E * I / (rho * A));    % Stiffness coefficient [m^2/s]

% Grid spacing and number of intervals
h = sqrt(1/2 * (c^2*k^2 + 4*sig1*k ...
    + sqrt((c^2*k^2 + 4*sig1*k)^2 + 16*kappa^2*k^2)));
N = floor(L/h);     % Number of intervals between grid points
h = L / N;          % Recalculation of grid spacing based on integer N

% Update coefficients
lambdaSq = c^2 * k^2 / h^2;
muSq = kappa^2 * k^2 / h^4;

% Boundary conditions ([c]lamped, [s]imply supported or [f]ree)
bc = "s";            

% Change N to the usable range
Norig = N;          % save original N
if bc == "c"        % Clamped: range reduces by 2 at each boundary 
    N = N - 4;
elseif bc == "s"    % Simply supported: range reduces by 1 at each boundary
    N = N - 2;
elseif bc == "f"    % Free: full range is used
    N = N;
end 
    
%% Initialise state vectors
uNext = zeros(N+1, 1); 
u = zeros(N+1, 1);

%% Initialise scheme matrices
Id  = eye(N+1);         % identity matrix

% Dxx and Dxxxx
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

Acoeff = (1 + sig0 * k);
B = 2 * Id + c^2 * k^2 * Dxx - kappa^2 * k^2 * Dxxxx + 2 * sig1 * k * Dxx;
C = -(1 - sig0 * k) * Id - 2 * sig1 * k * Dxx;

% calculate divisions of matrices before main loop for faster implementation
BoverA = B / Acoeff;
CoverA = C / Acoeff;

%% Initial conditions (raised cosine)
temporalExcitation = false;

ratio = 0.2;
loc = floor(ratio * (Norig+1));       % Center location
halfWidth = floor((Norig+1) / 20);    % Half-width of raised cosine
width = 2 * halfWidth;      % Full width
rcX = 0:width;              % x-locations for raised cosine
rc = 0.5 - 0.5 * cos(2 * pi * rcX / width); % raised cosine
    f = zeros(lengthSound, 1);

if temporalExcitation
    pluck = true;
    excitationLength = 30; %samples
    if pluck
        q = 1;
        excitationLength = excitationLength / 2;
    else
        q = 2;
    end
    f(1:excitationLength) = 0.01 - 0.01 * cos(q * pi * (1:excitationLength) / excitationLength);

end


excitation = zeros(length(u), 1);
% excitation(loc) = 1;
excitation(loc-halfWidth : loc+halfWidth) = rc;
if ~temporalExcitation
    u = excitation;
end
% Set initial velocity to zero
uPrev = u;

% Output location
outLoc = round(0.85 * (N+1));

%% Initialise matrices for energy
if bc == 'c'
    Dxp = sparse(1:Norig-2, 1:Norig-2, -ones(1, Norig-2), Norig-2, Norig-2) + ...
            sparse(1:Norig-3, 2:Norig-2, ones(1, Norig-3), Norig-2, Norig-2);
    DxxE = sparse(2:N+3, 1:N+2, ones(1, N+2), N+3, N+3) + ...
        sparse(1:N+3, 1:N+3, -2 * ones(1, N+3), N+3, N+3) + ... 
        sparse(1:N+2, 2:N+3, ones(1, N+2), N+3, N+3);
    DxxE = DxxE / h^2;

elseif bc == 's'
    Dxp = sparse(1:Norig, 1:Norig, -ones(1, Norig), Norig, Norig) + ...
            sparse(1:Norig-1, 2:Norig, ones(1, Norig-1), Norig, Norig);
elseif bc == 'f'
    Dxp = sparse(1:Norig+1, 1:Norig+1, -ones(1, Norig+1), Norig+1, Norig+1) + ...
    sparse(1:Norig, 2:Norig+1, ones(1, Norig), Norig+1, Norig+1);
    Dxp(end, :) = [];
    
    DxxE = sparse(2:N+1, 1:N, ones(1, N), N+1, N+1) + ...
        sparse(1:N+1, 1:N+1, -2 * ones(1, N+1), N+1, N+1) + ... 
        sparse(1:N, 2:N+1, ones(1, N), N+1, N+1);
    DxxE(1, :) = [];
    DxxE(end, :) = [];

    DxxE = DxxE / h^2;
end
Dxp = Dxp / h;

% Initialise energy vectors
kinEnergy = zeros(lengthSound, 1); % kinetic energy
potEnergy = zeros(lengthSound, 1); % potential energy
totEnergy = zeros(lengthSound, 1); % hamiltonian (total energy)
q0 = zeros(lengthSound, 1);        % frequency independent damping
q1 = zeros(lengthSound, 1);        % frequency dependent damping
qTot = 0;                          % total damping (summed form)

% initalise output vector
out = zeros(lengthSound, 1);

%% Simulation loop
for n = 1:lengthSound
        
    % Update equation 
    uNext = BoverA * u + CoverA * uPrev + temporalExcitation * f(n) * excitation;
    
    % Retrieve output
    out(n) = u(outLoc);
    
    %% Energy
    if bc == "c" % clamped
        
        % energy in the system
        kinEnergy(n) = rho * A * h / 2 * sum((1/k * (u - uPrev)).^2);
        potEnergy(n) = T / 2 * h * sum((Dxp * [0; u]) .* (Dxp * [0; uPrev]))...
                        + E * I * h / 2 * sum((DxxE * [0; u; 0]) .* (DxxE * [0; uPrev; 0]));

        % damping
        q0(n) = 2 * sig0 * rho * A * h * sum((1/(2*k) * (uNext - uPrev)).^2);
        q1(n) = - 2 * sig1 * rho * A * h * sum(1/(2*k) * ([0; uNext; 0] - [0; uPrev; 0])...
            .* (1/k * (DxxE * [0; u; 0] - DxxE * [0; uPrev; 0])));
        idx = n - (1 * (n~=1));
        qTot = qTot + k * (q0(idx) + q1(idx));

        % total energy 
        totEnergy(n) = kinEnergy(n) + potEnergy(n) + qTot;

    elseif bc == "s" % simply supported
        
        % energy in the system
        kinEnergy(n) = rho * A * h / 2 * sum((1/k * (u - uPrev)).^2);
        potEnergy(n) = T / 2 * h * sum((Dxp * [0; u]) .* (Dxp * [0; uPrev])) ...
            + E * I * h / 2 * sum((Dxx * u) .* (Dxx * uPrev));
        
        % damping
        q0(n) = 2 * sig0 * rho * A * h * sum((1/(2*k) * (uNext - uPrev)).^2);
        q1(n) = - 2 * sig1 * rho * A * h * sum(1/(2*k) * (uNext - uPrev)...
            .* (1/k * (Dxx * u - Dxx * uPrev)));
        idx = n - (1 * (n~=1));
        qTot = qTot + k * (q0(idx) + q1(idx));

        % total energy 
        totEnergy(n) = kinEnergy(n) + potEnergy(n) + qTot;

    elseif bc == "f" % free
        
        % energy in the system
        scaling = ones(N+1, 1);
        scaling(1) = 0.5;
        scaling(end) = 0.5;

        kinEnergy(n) = rho * A * h / 2 * sum(scaling .* (1/k * (u - uPrev)).^2);
        potEnergy(n) = T / 2 * h * sum((Dxp * u) .* (Dxp * uPrev)) ...
            + E * I * h / 2 * sum((DxxE * u) .* (DxxE * uPrev));

        % damping
        q0(n) = 2 * sig0 * rho * A * h * sum(scaling .* (1/(2*k) * (uNext - uPrev)).^2);
        q1(n) =  - 2 * sig1 * rho * A * h * sum(1/(2*k) * scaling .* (uNext - uPrev)...
            .* (1/k * (Dxx * u - Dxx * uPrev)));
        idx = n - (1 * (n~=1));
        qTot = qTot + k * (q0(idx) + q1(idx));
        
        % total energy 
        totEnergy(n) = kinEnergy(n) + potEnergy(n) + qTot;

    end
    
    %% Plot string state
    if drawThings && mod(n, drawSpeed) == 0
        
        if ~drawEnergy % plot string state
            
            % add boundary points for plot
            if bc == "c"
                uPlot = [0; 0; u; 0; 0];
            elseif bc == "s"
                uPlot = [0; u; 0];
            elseif bc == "f"
                uPlot = u;
            end

            % plot string
            hold off;
            plot((0:Norig), uPlot, 'k', 'Linewidth', 1.5)
%             hold on;
%             scatter(outLoc-1, uPlot(outLoc), 80, 'r', 'Marker', 'o', 'Linewidth', 2)
            % plot settings
            xlim([0, Norig])
            ylim([-1, 1])

            xticks([0 Norig])
            xticklabels({'$0$','$N$'})

            xlabel('$l$', 'interpreter', 'latex')
            ylabel('$u_l^n$', 'interpreter', 'latex')
%             title("String state")

            set(gca, 'Fontsize', 16, 'tickLabelInterpreter', 'latex', 'Linewidth', 2)

        else % plot normalised energy (should be within machine precision)
            
            plot((totEnergy(1:n) - totEnergy(1)) / totEnergy(1), 'k')
            title("Normalised energy")
            set(gca, 'Fontsize', 16, 'tickLabelInterpreter', 'latex', 'Linewidth', 2)

        end
        set(gca, 'Position', [0.0125 0.0188 0.9768 0.9577])
        set(gcf, 'color', 'w')
        axis off
        drawnow;
        for j = 1:slowdown
            M(frame) = getframe (gcf);
            frame = frame + 1;
        end
    end
    
    %% Update system states
    uPrev = u;
    u = uNext;
    
end
%%

if recordVid
    while isempty(M(end).cdata)
        M(end) = [];
    end
    v = VideoWriter('stiffStringHighStiffness.mp4', 'MPEG-4');
    open(v)
    writeVideo(v, M);
    close(v)
else
    figure('Position', [173 578 507 220])

    t = (1:lengthSound)/fs;
    plot(t, out, 'k', 'Linewidth', 1.5)


    ylim([-0.6, 0.6]);
    yLim = ylim;
    xLab = xlabel("$t$", 'interpreter', 'latex');
    yLab = ylabel("$u_{\fontsize{7}{7}\selectfont\textrm{out}}^n$", 'interpreter', 'latex');

    xLab.Position(2) = yLim(1) - (yLim(2) - yLim(1)) * 0.15;
    xLim = xlim;
    yLab.Position(1)  = xLim(1) -0.05 * (xLim(2) - xLim(1));
    % grid on;

    set(gca, 'Linewidth', 1.5, 'Fontsize', 16, ...
        'Position', [0.0871 0.2000 0.8875 0.7250], ...
        'TickLabelInterpreter', 'latex')
end