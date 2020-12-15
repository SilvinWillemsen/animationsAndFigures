%{
    Plot dispersion. Combinations for figure:
        quality = 1;        nIdx = 112;
        quality = 0.9;      nIdx = 128;
        quality = 1.001;    nIdx = 112;
%}

clear all;
close all;

quality = 1.001;
nIdx = 112;
figure('Position', [100, 400, 400, 300]);
calcBandwidth(quality, 5, nIdx);

xlim([0, 50])
xlabel("$l$", 'interpreter', 'latex');
ylabel("$u_l^n$", 'interpreter', 'latex')
grid on;
% title("$\lambda = " + num2str(quality) + "\quad n = " + num2str(nIdx) + "$", 'interpreter', 'latex', 'Fontsize', 20);
title("$\lambda > 1$", 'interpreter', 'latex', 'Fontsize', 20);
set(gca, 'FontSize', 16, 'Linewidth', 2);
set(gcf,'color','w');