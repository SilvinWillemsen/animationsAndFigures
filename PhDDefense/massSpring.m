close all

plotMassSprings = true;

fs = 441000;
lengthSound = 5000;
k = 1/fs;
c = 300;
h = c * k;
L = 1.0211;
N = floor(L/h);
h = L/N;
loc = 0.2;
locN = floor(loc*N+1);
offsetAmp1 = 1/locN;
offsetAmp2 = 1/(N+1-locN);
u = [linspace(offsetAmp1, 1, locN), linspace(1 - offsetAmp2, 0, N+1-locN)]';
uNext = zeros(N+1, 1);
uPrev = u;
recordVid = true;
slowdown = 0.25;
if recordVid
    loops = lengthSound * slowdown;
    M(loops) = struct('cdata',[],'colormap',[]);
    frame = 1;
end
stepsize = 100;
for n = 1:lengthSound 
    uNext(2:end-1) = 2 * u(2:end-1) - uPrev(2:end-1) ...
        + c^2 * k^2 / h^2 * (u(3:end) - 2 * u(2:end-1) + u(1:end-2));
    if recordVid
        if mod(n, 100) == 1
            hold off
            xlim([0, 1])
            ylim([-1, 1])
            if plotMassSprings
            spr = Spring(0.02, 5);
            for i = 1:stepsize:length(u)-stepsize
                point1 = [i/length(u), u(i)];
                point2 = [(i+stepsize)/length(u), u(i+stepsize)];
                [x, y] = spr.getSpr(point1, point2);
                plot (x, y, 'k')
                hold on

            end
        %     plot([(1:50:length(u)-50); (1:50:length(u)-50)+50], [u(1:50:end-50)';u((1:50:end-50) + 50)'], 'k')
            scatter((0:stepsize:length(u)-1)/length(u), u(1:stepsize:end), 40, 'w', 'filled', 'MarkerEdgeColor', 'k', 'Linewidth', 2)

            else
                plot((0:stepsize:length(u)-1)/length(u), u(1:stepsize:end), 'k', 'Linewidth', 1, 'Marker', '.', 'MarkerSize', 20)
            end
            xlim([0, 1])
            ylim([-1.2, 1.2])
            axis off
            set(gcf, 'color', 'w')
            drawnow;
            for j = 1:24
                M(frame) = getframe (gcf);
                frame = frame + 1;
            end
        end
    end
    uPrev = u;
    u = uNext;
end
if recordVid
    while isempty(M(end).cdata)
        M(end) = [];
    end
    if plotMassSprings
        v = VideoWriter('massSpring.mp4', 'MPEG-4');
    else
        v = VideoWriter('FDTD.mp4', 'MPEG-4');
    end
    open(v)
    writeVideo(v, M);
    close(v)
end