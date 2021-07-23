close all
figure('Position', [173 562 731 236])

t = 0:1/44100:0.022;
x = sawtooth(2*pi*220*t+0.85 * pi,7/8);
plot(t, x, 'k', 'Linewidth', 2)
ylim([-1.5, 1.5])
% xticks([])
% yticks([])
xlim([t(1), t(end)])
xlabel("$t$", 'interpreter', 'latex')
ylabel("$u(x_\textrm{\fontsize{7}{7}\selectfont out},t)$", 'interpreter', 'latex')
grid on
set(gca, 'Fontsize', 16, 'Linewidth', 2, 'Fontname', 'times', ...
    'Position', [0.0438 0.1780 0.9425 0.7754])