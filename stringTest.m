clear all;
close all;
fs = 8000;
t = 1:fs + 543;
y = sin(2 * pi * t / fs) + sin(3 * pi * t / fs) + sin(6.5 * pi * t / fs);

N = length(t);
uNext = zeros(N, 1);
u = y';
uPrev = y';
lambdaSq = 1;
N = N - 4;
B = sparse(2:N, 1:N-1, lambdaSq * ones(1, N-1), N, N) + ...
    sparse(1:N, 1:N, (2 - 2 * lambdaSq) * ones(1, N), N, N) + ...
    sparse(1:N-1, 2:N, lambdaSq * ones(1, N-1), N, N);
C = sparse(1:N, 1:N, -1 * ones(1, N), N, N);
N = N + 4;
figure('Renderer', 'painters', 'Position', [200 200 600 200])
loops = 3001;
M(100) = struct('cdata',[],'colormap',[]);
set(gca, 'Color', 'white')
frame = 1;
continuous = false;
bc = 'ss';
for n = 1:loops
    uNext(3:N-2) = B * u(3:N-2) + C * uPrev(3:N-2);
%     uNext(172:N-172) = B(172:N-172, 172:N-172) * u(172:N-172) + C(172:N-172, 172:N-172) * uPrev(172:N-172);
    uPrev = u;
    u = uNext;
    if mod(n, 1) == 0
        if continuous
            plot(uNext, 'k', 'LineWidth', 5)
            s = floor((n*1000 - 1) / fs) / 1000;
            if mod(s * 1000, 1) == 1
                string = ['t = ', num2str(s), '0 s'];
            else
                string = ['t = ', num2str(s), ' s'];
            end
        else
            clf;
            scatter(1:171:N + 171, [uNext(1:171:N); 0], 400, 'k', '.')
            if strcmp(bc, 'ss')
                hold on;
                scatter([1-171, N+171], -[uNext(172) uNext(8380)], 40, 'white');
            elseif strcmp(bc, 'clamped')
                hold on;
                scatter([1-171, N+171], [0 0], 400, 'k', '.');
%                 scatter([1-2*171, N+2*171], [0 0], 40, 'k', 'MarkerEdgeColor',[0 0 0]);
            end
                string = ['n = ', num2str(n)];
        end
%         title (string, 'horizontalAlignment', 'center');
        axis off;
        set(gca, 'FontSize', 15, 'Color', 'white')
        set(gcf, 'Color', 'w')
        xlim([0, 9000])
        ylim([-3, 3])
        drawnow;
        M(frame) = getframe(gcf);
        frame = frame + 1;
    end
end
% fig = figure;
% movie(fig,M,2)
v = VideoWriter('stringTestSS.mp4', 'MPEG-4');
open(v)
writeVideo(v, M);
close(v)