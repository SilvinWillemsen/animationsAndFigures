stringState = zeros(100, 1);
stringState(40:60) = hann(21);
plot(stringState, 'k', 'Linewidth', 2)
plot([zeros(50,1); stringState], 'k', 'Linewidth', 2)
hold on;
plot(zeros(50,1), 'r', 'Linewidth', 2)

ylim([-1.5,1.5])
axis off
set(gcf, 'color', 'w')