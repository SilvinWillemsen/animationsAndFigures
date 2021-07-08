%{
    Plot dispersion. Combinations for figure:
        quality = 1;        nIdx = 112;
        quality = 0.9;      nIdx = 128;
        quality = 1.001;    nIdx = 112;
%}

clear all;
close all;

combinefigs = true;
qualityVec = fliplr([1, 0.8, 1.001]);
if combinefigs
    range = 1:length(qualityVec);
else
    range = 1
end
figure('Position', [100, 400, 400, 300]);
scaling = 1.8;
offset = 0.35;
textVals = {"$\lambda > 1$", "$\lambda < 1$", "$\lambda = 1$"};
letters = {"(c)", "(b)", "(a)"};

for i = range
    quality = qualityVec(i);
    if quality == 0.8
        nIdx = 110;
    else
        nIdx = 112;
    end
    calcBandwidth(quality, 5, nIdx, i * scaling);
    text(47.8, (i + offset) * scaling, letters{i}, ...
        'interpreter', 'latex', 'Fontsize', 16, ...
        'horizontalAlignment', 'center', ...
        'verticalAlignment', 'middle')
    text(25, (i + offset) * scaling, textVals{i}, ...
        'interpreter', 'latex', 'Fontsize', 16, ...
        'horizontalAlignment', 'center', ...
        'verticalAlignment', 'middle')

    if combinefigs
        hold on;
    end
end
xlim([0, 50])
xlabel("$l$", 'interpreter', 'latex');
ylabel("$u_l^n$", 'interpreter', 'latex')
grid on;

ylim([0.5, (length(qualityVec) + 0.5)] * scaling)
ytickLocs = [(1+0.5) * scaling, (2+0.5) * scaling];
ylabelSave = ['', ''];
xlabelSave = [];
xtickLocs = [];
% title("$\lambda = " + num2str(quality) + "\quad n = " + num2str(nIdx) + "$", 'interpreter', 'latex', 'Fontsize', 20);
% title("$\lambda > " + num2str(quality) + "$", 'interpreter', 'latex', 'Fontsize', 20);
set(gca, 'FontSize', 16, 'Linewidth', 2, ...
     'yticklabel', ylabelSave, 'YTick', ytickLocs, ...
     'xticklabel', xlabelSave, 'XTick', xtickLocs, ...
     'Position', [0.0625 0.08 0.92 0.90]);
set(gcf,'color','w');