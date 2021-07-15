
N = 1000;
x=0:1/(N-1):1;
startVal = (ceil(0.7 * length(x)) - 1) / ceil(0.7 * length(x));
y = [linspace(0, 1, floor(0.3 * N)), ...
     linspace(startVal, 0, ceil(0.7 * N))]

scatter(x, y)
ylim([-1.1, 1.1])

width = N / 10;
loc = N / 5;
startLoc = floor(loc - (width/2));
y = zeros(N, 1);
y(startLoc:startLoc+width) = hann(width+1);

width = N / 5;
loc = 7 * N / 10;
startLoc = floor(loc - (width/2));
y2 = zeros(N, 1);

y2(startLoc:startLoc+width) = -0.5 * hann(width+1);

plot(y)
hold on;
plot(y2)
ylim([-1.1, 1.1])
