close all;
clear all;

%% Initialise variables
fs = 44100;             % sample rate [Hz]
k = 1 / fs;             % time step [s]
lengthSound = fs;       % length of the simulation (1 second) [samples]

f0 = 0/(k*pi);               % fundamental frequency [Hz]
omega0 = 2 * pi * f0;   % angular (fundamental) frequency [Hz]
M = 1;                  % mass [kg]
K = omega0^2 * M;       % spring constant [N/m]

%% initial conditions (u0 = 1, d/dt u0 = 0)
u = 1;                  
uPrev = 1;

% initialise output vector
out = zeros(lengthSound, 1);

%% Simulation loop
for n = 1:lengthSound
    
    % Update equation 
    uNext = (2 - K * k^2 / M) * u - uPrev; 
    
    out(n) = u;
    if mod(n, 100) == 0
        plot(out(1:n))
        drawnow;
    end
    % Update system states
    uPrev = u;
    u = uNext;
end


%% Plotting
figure('Position', [173 578 827 220])

t = 1:lengthSound;
% plot(t, zeros(length(t), 1), '--', 'Linewidth', 2, 'color', [0.5, 0.5, 0.5])
% hold on;
subp1 = subplot(1, 2, 1)
plot(t, out, 'k', 'Linewidth', 2)

xlim([0, 500])
ylim([-1.1, 1.1])
xLab = xlabel("$n$", 'interpreter', 'latex');
yLab = ylabel("$u^n$", 'interpreter', 'latex');

xLab.Position(2) = -1.35;
yLab.Position(1) = -30;
% grid on;

set(gca, 'Linewidth', 1.5, 'Fontsize', 16, ...
    'Position', [0.0532 0.2000 0.4115 0.7591], ...
    'TickLabelInterpreter', 'latex')
subplot(1, 2, 2)
dbFFT = 20 * log10(abs(fft(out)));
plot([440, 440], [0, dbFFT(440)], '--', 'color', [0.5, 0.5, 0.5], 'Linewidth', 2)
hold on;
plot(0:lengthSound-1, dbFFT, 'k', 'Linewidth', 2);
xlim([0, 1000])
% grid on;
xticks([0, 300, 440, 600, 900])
set(gca, 'Linewidth', 1.5, 'Fontsize', 16,...
    'Position', [0.5532 0.2000 0.4115 0.7591], ...
    'TickLabelInterpreter', 'latex')

xLab2 = xlabel("$f$ [Hz]", 'interpreter', 'latex');
yLab2 = ylabel("Magnitude [dB]", 'interpreter', 'latex');

xLab2.Position(1) = 750;
xLab2.Position(2) = -8;
yLab2.Position(1) = -70;



set(gcf, 'color', 'w')