clear all
close all;

% samplerate
fs = 44100;
k = 1/fs;

% user-defined variables
f0 = 150;
L = 1;
c = f0 * 2 * L;
sig0 = 0;
sig1 = 0.00;

% calculate grid spacing
h = sqrt(c^2 * k^2 + 4 * sig1 * k);
N = floor (L / h); % number of sections
h = L / N;

% define lambda squared
lambdaSq = (c * k / h)^2

%% initialise state vectors
uNext = zeros(N+1, 1);
u = zeros(N+1, 1);
u(3:9) = hann(7);
uPrev = u;


lengthSound = fs;

range = 2:N; % l = [1, ..., N-1]
outLoc = floor(N/pi);
% for loop
for n = 1:lengthSound
    % update equation
    uNext(range) = (2 * u(range) - uPrev(range) ...
        + lambdaSq * (u(range+1) - 2 * u(range) + u(range-1)) + sig0 * k * uPrev(range) ...
        + 2 * sig1 * k / h^2 * (u(range+1) - 2 * u(range) + u(range-1) - uPrev(range+1) + 2 * uPrev(range) - uPrev(range-1))) / (1+sig0 * k);
%     uNext(range) = 2 * u(range) - uPrev(range) ...
%         + lambdaSq * (u(range+1) - 2 * u(range) + u(range-1))- 2* sig0 * k * (u(range) - uPrev(range));

    out(n) = u(outLoc);
%     plot(u);
%     drawnow;
    
    uPrev = u;
    u = uNext;
end
plot ([0:lengthSound-1]/fs, out)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel("$u^n_{" + num2str(outLoc) + "}$", 'interpreter', 'latex')

set(gca, 'Linewidth', 2, 'Fontsize', 16)
set(gcf, 'color', 'w')

figure;
spectrogram(out,512,64,512, 44100, 'yaxis');
set(gca, 'Fontsize', 16, 'FontName', 'times')

