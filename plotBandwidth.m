%{
    File calculating bandwidth manually. Proves that this follows the same
    curve as \f_{max}= f_s/\pi * sin^{-1}(\lambda).
%}

clear all;
close all;

quality = 1;
[fftOut, ~] = calcBandwidth(quality, 2, 0);


figure('Position', [100, 400, 400, 300]);
plot(10 * (0:length(fftOut)-1), 20 * log10(fftOut), 'k', 'Linewidth', 1);
% plot(10 * (0:length(fftOut)-1), fftOut, 'k', 'Linewidth', 1);

xlabel("$f$ (in Hz)", 'interpreter', 'latex');
ylabel("Magnitude Response (in dB)", 'interpreter', 'latex')
ylim([-60,60])
grid on;
title("$\lambda = " + num2str(quality) + "$", 'interpreter', 'latex');
set(gca, 'FontSize', 16, 'Linewidth', 2);
set(gcf,'color','w');