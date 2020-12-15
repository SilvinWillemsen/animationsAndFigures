close all;
clear all;
N = 100;
alpha = zeros(N,1);
x = (1:N) / N;
for j = 0:0.01:0.99
    epsilon = [linspace(1, 1, N/4), ((cos(2 * pi * [1:N/2]/(N/2))+1) * 0.5) * j + (1-j), linspace(1, 1, N/4)];
%     epsilon = [linspace(1, 1, N/3), linspace(0.5, 0.5, N/3), linspace(1, 1, N/3)];

    subplot(2,1,1)
    plot(x, epsilon, 'Linewidth', 2)
    ylabel('$\epsilon(x)$','interpreter', 'latex');
    xlabel('$x$','interpreter', 'latex');
    ylim([-0.1 1.1])
    set(gca, 'Fontsize', 16);
    
    
    for j = 1:N
        alpha(j)= sum(epsilon(1:j))/sum(epsilon);
    end

    subplot(2,1,2)
    hold off;
    scatter (x, ones(N,1));
    hold on;
    scatter(alpha, zeros(N,1));
    leg1 = legend('$x$', '$\alpha$');
    set(leg1,'Interpreter','latex');
    ylim([-1, 2]);
    xlabel('Scaled location')
    set(gca, 'Fontsize', 16);
    drawnow;
end