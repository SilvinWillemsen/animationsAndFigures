clear all;
close all;
clc;
width = 100;
N = 500;
excitationLength = 10;
drawLength = 14;
excitationFunction = zeros (N, drawLength);
one = ones(N,1);
loc = N/2; 

figure('Renderer', 'painters', 'Position', [100 100 754 300])
AxesH = axes;
drawnow;
InSet = get(gca, 'TightInset');
InSet(1) = InSet(1) + 0.05;
offset = 2;
q = 1;
set(gca, 'Position', [InSet(1:2), 1-InSet(1)+0.075, 1-InSet(2)+InSet(4)]);
for i = -offset:drawLength - 1
    if i >= 0 && i <= excitationLength
        excitationForce = (1 - cos(q * pi * i / excitationLength)) * 0.5;
    else
        excitationForce = 0;
    end
    
    excitationFunction(loc - width / 2:loc + width / 2, i + 1 + offset) = (1 - cos(2 * pi * [0:width] / width)) * 0.5 .* excitationForce;
    plot3(0:N-1, one * i, excitationFunction(:, i + 1 + offset), 'k', 'LineWidth', 2);
    hold on;
end
view(60,40)
xticks([]);
yticks([]);
zticks([]);
ylim([-offset-1, drawLength])
xlabel("$x$", 'interpreter', 'latex', 'Position',[N/2 -offset-1 0])
ylabel("$t$", 'interpreter', 'latex', 'Position', [N + 125 (drawLength - offset - 2)/2 0]);
zlabel("$e(x)f(t)$", 'interpreter', 'latex')
%, 'XColor', 'white', 'YColor', 'white', 'ZColor', 'white'
scatter3(loc - 1, 5, 0, 400, [0, 0.5, 0], '.');
text(loc-1.2, 3.75, "$x_0$", 'Fontsize', 20, 'interpreter', 'latex', 'Color', [0, 0.5, 0]);
scatter3(N, 0, 0, 400, 'b', '.');
text(N*1.1,-1.5, "$t_0$", 'Fontsize', 20, 'interpreter', 'latex', 'Color', 'b');
plot3([loc - width / 2 - 1, loc + width / 2 - 1], [excitationLength - 3 excitationLength - 3], [0 0], 'r', 'LineWidth', 2);
text(loc - 1.2, excitationLength - 4.25, "$x_\textrm{\fontsize{7}{7}\selectfont w}$", 'Fontsize', 20, 'interpreter', 'latex', 'Color', 'r');
plot3([N N], [0 excitationLength], [0 0], 'b', 'LineWidth', 2);
text(N + 50, excitationLength/2, "$t_\textrm{\fontsize{7}{7}\selectfont d}$", 'Fontsize', 20, 'interpreter', 'latex', 'Color', 'b');


xlim([0 N+125])
set(gca, 'Fontsize', 20, 'Position', [0.0660 0.0250 0.9233 0.9658])
grid on;
box on;