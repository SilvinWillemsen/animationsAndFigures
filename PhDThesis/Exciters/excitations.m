clear all;
% close all

%% raised cosine
%continuous
L = 2;
x = 0:0.01:L;
xw = 0.2;
x0 = 0.3;
xs = x0 - xw/2;
xe = x0 + xw/2;

idx = 1;
e = zeros(size(x))';
for xx = x
    if xx >= xs && xx <= xe
        e(idx) = 0.5 - 0.5 * cos(2*pi * ((xx)-xs) / xw);
    end
    idx = idx + 1;
end

plot(x, e)

%% discrete
w = 11;
N = 100;
l = 0:N;
l0 = 5;
ls = l0 - floor(w/2);
le = l0 + floor(w/2);
E = zeros(size(l))';
idx = 1;
for ll = l
    if ll >= ls && ll <= le
        E(idx) = 0.5 - 0.5 * cos(2*pi * (ll-ls) / w);
    end
    idx = idx + 1;
end
plot(l, E)

%% pluck triangle
%continuous
L = 2;
x = 0:0.01:L;
idx = 1;
eamp = 1;
x0 = 1
for xx = x
    if xx <= x0
        eTri(idx) = eamp / x0 * xx;
    else
        eTri(idx) = eamp / (x0 - L) * (xx - L);

    end
    idx = idx + 1;
end
plot(eTri)

%% discrete
N = 100;
l0 = 40;
l = 0:N;
idx = 1;
for ll = l
    if ll <= l0
        Etri(idx) = eamp / l0 * ll;
    else
        Etri(idx) = eamp / (l0 - N) * (ll - N);

    end
    idx = idx + 1;
end
plot(Etri)


%% 2DRaised cos
Lx = 2;
Ly = 1;
x = 0:0.01:Lx;
y = 0:0.01:Ly;
rw = 0.1;
x0 = 0.3;
y0 = 0.25;

xs = x0 - rw/2;
xe = x0 + rw/2;
ys = y0 - rw/2;
ye = y0 + rw/2;

idxX = 1;
idxY = 1;
e = zeros(length(y), length(x));
for yy = y
    for xx = x
        if xx >= xs && xx <= xe && yy >= ys && yy <= ye
            e(idxY, idxX) = 0.25 * (1-cos(2 * pi * (xx - xs) / rw)) * (1-cos(2 * pi * (yy - ys) / rw));   
        else
            e(idxY, idxX) = 0;
        end
        idxX = idxX + 1;
    end
    idxY = idxY + 1;
    idxX = 1;
end
imagesc(e)

%% discrete
Nx = 100;
Ny = 100;
l = 0:Nx;
m = 0:Ny;
l0 = floor(0.3 * Nx);
m0 = floor(0.75 * Ny);
rwDisc = 10;
ls = l0 - floor(rwDisc/2);
le = l0 + floor(rwDisc/2);
ms = m0 - floor(rwDisc/2);
me = m0 + floor(rwDisc/2);


e = zeros(Ny-1, Nx-1);
idxX = 1;
idxY = 1;

for mm = m
    for ll = l
        if ll >= ls && ll <= le && mm >= ms && mm <= me
            e(idxY, idxX) = 0.25 * (1-cos(2 * pi * (ll - ls) / rwDisc)) * (1-cos(2 * pi * (mm - ms) / rwDisc));   
        else
            e(idxY, idxX) = 0;
        end
        idxX = idxX + 1;
    end
    idxY = idxY + 1;
    idxX = 1;
end
subplot(211)
imagesc(e)

%% using hann windows
% e = zeros(Ny-1, Nx-1);
% % discrepancy is due to the 1-based nature of matlab
% e(ms+1:me+1, ls+1:le+1) = hann(rwDisc+1) * hann(rwDisc+1)';
% subplot(212)
figure('Position', [489 604 578 253])
view(20,35);
surf(hann(50) * hann(50)')
colormap gray;
caxis([-0.5, 1]);
box on
xticks([])
yticks([])
zticks([])
set(gca, 'Fontsize', 16, 'tickLabelInterpreter', 'latex', 'Linewidth', 1.5, ...
    'Position', [0.0700 0.1100 0.9250 0.8150])

yLab = ylabel('$l$', 'interpreter', 'latex');
% yLab.Position = [24.9123143849772,-3.42032538420517,-0.049501856347039];

xLab = xlabel('$m$', 'interpreter', 'latex');
% xLab.Position = [51.9344402355818,30.6832563101389,-0.01185608391491];
zlabel('$E_{(l,m), \textrm{\fontsize{7}{7}\selectfont rc}}$', 'interpreter', 'latex')
