if firstPlot
    figure('Position', [173 578 827 220])
    color = 'b'
else
    color = 'r'
end
% create Q matrix (one-step form)
Q = [Amat \ B, Amat \ C;
     eye(N+1), zeros(N+1)];
 
% obtain complex frequencies
s = 1/k * log(eig(Q));

% obtain positive frequencies and sort them
s = s(imag(s) >= 0);
[~, order] = sort(imag(s));
s = s(order);
subplot(121)
if firstPlot
    plot(0, 0)
end
hold on
scatter(1:length(s), imag(s)/(2*pi), color, 'Linewidth', 2)
grid on
xlabel("Mode number $p$", 'interpreter', 'latex')
ylabel("$f_p$ [Hz]", 'interpreter', 'latex')
xlim([1, length(s)])
title("Modal frequency")
set(gca, 'Fontsize', 16, 'Linewidth', 2, 'FontName', 'times', ...
    'Position', [0.0632 0.2091 0.4115 0.6909])

subplot(122)
if firstPlot
    plot(0, 0)
end
hold on
scatter(1:length(s), real(s), color, 'Linewidth', 2)
grid on
xlim([1, length(s)])
title("Damping per mode")
xlabel("Mode number $p$", 'interpreter', 'latex')
ylabel("$\sigma_p$ [s$^{-1}$]", 'interpreter', 'latex')
set(gca, 'Fontsize', 16, 'Linewidth', 2, 'FontName', 'times', ...
    'Position', [0.5732 0.2091 0.4115 0.6909])
