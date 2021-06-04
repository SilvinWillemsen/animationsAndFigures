figure('Position', [440 563 560 235])

vrel = -1:0.001:1;
a = 100;

phi = sqrt(2 * a) .* vrel .* exp(-a * vrel.^2 + 1/2);
plot(vrel, phi, 'k', 'Linewidth', 2)
ylim([-1.1, 1.1])
grid on;
xlabel("$v_\textrm{\fontsize{7}{0}\selectfont rel}$", 'interpreter', 'latex')
ylabel("$\Phi(v_\textrm{\fontsize{7}{0}\selectfont rel})$", 'interpreter', 'latex')

set(gca, 'Linewidth', 2, 'Fontsize', 16);