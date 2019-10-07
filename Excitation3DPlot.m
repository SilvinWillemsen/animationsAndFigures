clear all;
close all;
clc;
width = 100;
N = 500;
excitationLength = 10;
drawLength = 15;
excitationFunction = zeros (N, drawLength);
one = ones(N,1);
loc = N/2; 

figure('Renderer', 'painters', 'Position', [100 100 700 300])
AxesH = axes;
drawnow;
InSet = get(gca, 'TightInset');
InSet(1) = InSet(1) + 0.05;
set(gca, 'Position', [InSet(1:2), 1-InSet(1)+0.075, 1-InSet(2)+InSet(4)]);
for i = 0:drawLength - 1
    if i <= excitationLength
        excitationForce = (1 - cos(1 * pi * i / excitationLength)) * 0.5;
    else
        excitationForce = 0;
    end
    
    excitationFunction(loc - width / 2:loc + width / 2, i + 1) = (1 - cos(2 * pi * [0:width] / width)) * 0.5 .* excitationForce;
    plot3(0:N-1, one * i, excitationFunction(:, i + 1), 'k', 'LineWidth', 2);
    hold on;
end
view(60,40)
xticks([]);
yticks([]);
zticks([]);
xlabel("$x$", 'interpreter', 'latex', 'Position',[N/2 0 0])
ylabel("$t$", 'interpreter', 'latex', 'Position', [N + 125 drawLength/2 0]);
zlabel("$E_eF_e$", 'interpreter', 'latex')
%, 'XColor', 'white', 'YColor', 'white', 'ZColor', 'white'
scatter3(loc - 1, 5, 0, 400, [0, 0.5, 0], '.');
text(loc-1.2, 3.75, "$x_e$", 'Fontsize', 20, 'interpreter', 'latex', 'Color', [0, 0.5, 0]);
plot3([loc - width / 2 - 1, loc + width / 2 - 1], [excitationLength - 3 excitationLength - 3], [0 0], 'r', 'LineWidth', 2);
text(loc - 1.2, excitationLength - 4.25, "$w_e$", 'Fontsize', 20, 'interpreter', 'latex', 'Color', 'r');
plot3([N N], [0 excitationLength], [0 0], 'b', 'LineWidth', 2);
text(N + 50, excitationLength/2, "$d_e$", 'Fontsize', 20, 'interpreter', 'latex', 'Color', 'b');


xlim([0 N+125])
set(gca, 'Fontsize', 20)
grid on;