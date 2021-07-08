figure('Position', [508 569 581 229])
s0 = 10^4;
v = -0.5:0.001:0.5;
muC = 0.3;
muS = 0.8;
fN = 5;
fC = muC * fN;
fS = muS * fN;

vS = 0.1;

zss = sign(v) / s0 .* (fC  + (fS - fC) * exp(-(v/vS).^2));
plot(v, zss, 'k', 'Linewidth', 2)
ylim([-5e-4, 5e-4])
xticks([-0.5:0.1:0.5])
yticks([-5e-4:2.5e-4:5e-4])
xlabel('$v$', 'interpreter', 'latex', 'Fontsize', 16)
ylab = ylabel('$z_\textrm{\fontsize{7}{7}\selectfont ss}(v)$', 'interpreter', 'latex', 'Fontsize', 16)
ylab.Position(1) = -0.55
grid on
set(gca, 'Fontsize', 16, 'Linewidth', 2, 'Position', [0.0935 0.2031 0.8844 0.6913])