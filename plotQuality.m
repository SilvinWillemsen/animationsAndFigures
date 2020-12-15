%{
    File calculating bandwidth manually. Proves that this follows the same
    curve as \f_{max}= f_s/\pi * sin^{-1}(\lambda).
%}

clear all;
close all;
j = 1; 
stepSize = 0.005;
res = zeros(1/stepSize, 1);
range = stepSize:stepSize:1;
% for i = range
%     [~, res(j)] = calcBandwidth(i, 2, 1);
%     j = j+1;
% end

figure('Position', [100, 400, 400, 235]);
% plot(range, res, 'k', 'Linewidth', 2);
% hold on;
plot(range, 2/pi * asin(range), 'k', 'Linewidth', 2)
xlabel("$\lambda$", 'interpreter', 'latex');
ylabel("$f_\textrm{\fontsize{7}{0}\selectfont max} (\times0.5f_\textrm{\fontsize{7}{0}\selectfont s}$)", 'interpreter', 'latex')
grid on;
set(gca, 'FontSize', 16, 'Linewidth', 2, 'xtick', linspace(0, 1, 6));
set(gcf,'color','w');