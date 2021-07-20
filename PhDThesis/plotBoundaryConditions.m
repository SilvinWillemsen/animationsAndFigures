close all;
clear all;

figure('Position', [403 604 305 194])
gcaPos = [0.1180 0.1237 0.8566 0.8351];
width = 60;
range = (0:width-1)/width;
lW = 1.5;
bc = "ss";
if bc == "c"
    plot(range, hann(width), 'k', 'Linewidth', lW)
    hold on;
    plot(range, -hann(width), '--', 'color', [0.5, 0.5, 0.5], 'Linewidth', lW)
    
elseif bc == "ss"
    hannSS = 2*hann(width*2) - 1;
    hannSS = hannSS(hannSS >= 0);
    plot(range, hannSS, 'k', 'Linewidth', lW)
    hold on;
    plot(range, -hannSS, '--', 'color', [0.5, 0.5, 0.5], 'Linewidth', lW)
    
elseif bc == "f"
    hannSS = 2*hann(width*2) - 1;
    hannSS = hannSS(hannSS >= 0);

    plot(range, hannSS * 2 - 1 , 'k', 'Linewidth', lW)
    hold on;
    plot(range, -hannSS * 2 + 1, '--', 'color', [0.5, 0.5, 0.5], 'Linewidth', lW)
    
end

set(gca, 'Linewidth', lW, 'Position', gcaPos, 'Fontsize', 16)
set(gcf, 'color', 'w')
