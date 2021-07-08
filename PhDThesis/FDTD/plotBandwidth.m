%{
    File calculating bandwidth manually. Proves that this follows the same
    curve as \f_{max}= f_s/\pi * sin^{-1}(\lambda).
%}

clear all;
close all;

quality = 1;
combinefigs = true;
qualityVec = ([1, 0.9, 0.5]);

figure('Position', [100, 400, 400, 300]);

if combinefigs
    range = 1:length(qualityVec);
else
    range = 1
end

scaling = 1.8;
offset = 0.35;
textVals = {"$\lambda = 1$", "$\lambda = 0.9$", "$\lambda = 0.5$"};
letters = {"(a)", "(b)", "(c)"};
xlabelSave = [];

topSpace = 0.005;
bottomSpace = 0.16;
rightSpace = 0.005;

leftSpace = 0.15;

for i = range
    subplot(3, 1, i)
    [fftOut, ~] = calcBandwidth(qualityVec(i), 2, 0);
    plot(0.01 * (0:length(fftOut)-1), 20 * log10(fftOut), 'k', 'Linewidth', 1);
    xlim([0, 24])
    ylim([-50, 80])

    grid on
    
    height = (1-bottomSpace) / 3;
    
    set(gca, 'FontSize', 16, 'Linewidth', 1.5, ...
        'Position', [leftSpace, 1 - (height * i)- topSpace, 1-leftSpace - rightSpace, height]);
    if i == 2
        ylabel("Magnitude Response (in dB)", 'interpreter', 'latex')
    end
    set(gca, 'XTick', (0:5) * 5);
    set(gca, 'YTick', [-40, 0, 40]);

    if i ~= 3
        set(gca, 'xticklabel', ['', '', '', '', '', '']);
    end
    test = xlim;
    text(test(2) - 1, 60, letters{i}, ...
        'interpreter', 'latex', 'Fontsize', 16, ...
        'horizontalAlignment', 'center', ...
        'verticalAlignment', 'middle')
    text(test(2) * 0.5, 60, textVals{i}, ...
        'interpreter', 'latex', 'Fontsize', 16, ...
        'horizontalAlignment', 'center', ...
        'verticalAlignment', 'middle')

end
% plot(10 * (0:length(fftOut)-1), fftOut, 'k', 'Linewidth', 1);

xlabel("$f$ (in kHz)", 'interpreter', 'latex');

grid on
% title("$\lambda = " + num2str(quality) + "\quad n = " + num2str(nIdx) + "$", 'interpreter', 'latex', 'Fontsize', 20);
% title("$\lambda > " + num2str(quality) + "$", 'interpreter', 'latex', 'Fontsize', 20);
set(gcf,'color','w');