close all; 
clear all;
%% Pulse train generator
fs = 44100;             % Sample rate [Hz]
lengthSound = fs;       % Length of the sound [samples]
f = 440;                % Pulse train frequency
dutyC = 1;           % Duty cycle [0-1]
amp = 1;                % Amplitude

%% Create input signal
for n = 0:lengthSound
    if mod(n, fs / f) <= fs / f * dutyC
        vIn(n+1) = amp *sin(f * pi / dutyC * mod(n, fs / f) / fs);
    else
        vIn(n+1) = 0;
    end
end
figure('Position', [173 578 786 220])
plot(0:lengthSound, vIn, 'k', 'Linewidth', 2)

xlim([0, 900])
yLim = ylim;
ylim([yLim * 1.1])

xLab2 = xlabel("$n$", 'interpreter', 'latex');
yLab2 = ylabel("$v_\textrm{\fontsize{7}{7}\selectfont in}$", 'interpreter', 'latex', 'Fontsize', 20);

yLim = ylim;
xLab2.Position(2) = yLim(1) -0.07 * (yLim(2) - yLim(1));

yticks([])
xticks(0:100:900)
set(gca, 'Linewidth', 1.5, 'Fontsize', 16,...
    'Position', [0.0356 0.1500 0.9338 0.8091], ...
    'TickLabelInterpreter', 'latex')
set(gcf, 'color', 'w')
