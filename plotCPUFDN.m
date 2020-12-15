close all;
clear all;
data = [11.0, 18.4, 3.11;...
        12.0, 19.8, 5.41;...
        15.2, 22.7 -- 7.9;...
        22.2, 28.6 -- 13.3;...
        40.2, 46.1 -- 30.4;...
        100.1, 110.5, 92.5];
    
data1 = data(:,1);
data2 = data(:,2);
data3 = data(:,3);
range1 = 1:6;
range2 = 2.^(1:6);
semilogx(2.^(1:6), data(:,1), 'r', 'Linewidth', 2);
hold on;
semilogx(2.^(1:6), data(:,2), 'b', 'Linewidth', 2);
semilogx(2.^(1:6), data(:,3), 'k', 'Linewidth', 2);
legend('Unfixed (EQ)', 'Unfixed (IR)', 'Fixed', 'Location', 'northwest');
xticks([2.^(1:6)])
xlim([2, 64])
set(gca, 'Linewidth', 2, 'Fontsize', 16)
xlabel('FDN order');
ylabel('CPU (%)');

grid off;
grid on;
