close all;
x = 1;          % starting point
eps = 1e-4;     % threshold
% if the threshold has not been crossed before this amount of iterations, do not iterate more
maxIterations = 100;
xRange = -0.5:0.01:1.5;
figure('Position', [440 450 295 348])

% loop until a maximum number of iterations
for i = 1:maxIterations
    
    % calculate next iteration (Eq. (8.11)
    xNext = x - (exp(x) - 1) / (exp(x));
    hold off;
    plot([min(xRange), max(xRange)], [0, 0], 'color', [0.75, 0.75, 0.75], 'Linewidth', 1.5)
    hold on;
    plot([0, 0], [-4, exp(xRange(end))], 'color', [0.75, 0.75, 0.75], 'Linewidth', 1.5)
    plot(xRange, exp(xRange) - 1, 'k', 'Linewidth', 2)
    b = exp(x) - exp(x) * x - 1;
    plot(xRange, exp(x) * xRange + b, 'r', 'Linewidth', 1.5)
    scatter(x, exp(x) - 1, 'b', 'filled');
    scatter(xNext, 0, 'g', 'filled');
    if i == 1
        text(x - 0.15, exp(x) - 1+0.15, '$x_i$', 'Fontsize', 18, 'interpreter', 'latex')
        text(xNext + 0.05, -0.15, '$x_{i+1}$', 'Fontsize', 18, 'interpreter', 'latex')
    end
    xlim([-0.5000    1.30])
    ylim([-1.5    3.5])
    xlabel('$x$', 'interpreter', 'latex')
    ylabel('$f(x)$', 'interpreter', 'latex')
    title("$i =\ $" + i, 'interpreter', 'latex')
    set(gca, 'Position', [0.1541 0.1379 0.7956 0.7871], 'Linewidth', 2, ...
        'Fontsize', 16)
    grid on
    drawnow;
    
    % threshold check (Eq. (8.12)
    if abs(xNext - x) < eps
        break; % break out of the for loop
    end
    % update the value of x
    x = xNext;
end
disp("The root of f(x) is at x = " + xNext)