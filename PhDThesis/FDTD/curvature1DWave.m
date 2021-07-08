clear all;
close all;

figure('Position', [0, 300, 600, 200])
fs = 44100;
k = 1/fs;

plotDiscrete = false;
lengthSound = 300;

c = 1470/10;
h = c*k;
L = 1;
N = floor(L/h);
% h = L/N;

lambdaSq = (c*k/h)^2;

uNext = zeros(N+1, 1);
u = zeros(N+1, 1);

widthU = 101;

halfWidthU = floor(widthU/2);

u(floor(N/2)-halfWidthU : floor(N/2) + halfWidthU) = hann(widthU); 
uPrev = u;

recordVid = false;

if recordVid
    loops = 500;
    M(lengthSound) = struct('cdata',[],'colormap',[]);
    frame = 1;
end

range = 2:N;
for n = 1:lengthSound
    uNext(range) = 2*u(range) - uPrev(range) + lambdaSq * (u(range+1) - 2*u(range) + u(range-1));

    %% drawthings
    if mod(n,1) == 0
        if plotDiscrete
            plot(0:N, u, 'k', 'Marker', '.', 'Markersize', 20);
            xlabel('$l$', 'interpreter', 'latex')
            ylabel("$u_l^n$", 'interpreter', 'latex')
            title("$n = " + num2str(n-1)+"$", 'interpreter', 'latex');
        else
            hold off;
            plot((0:N) / N, u, 'k', 'Linewidth', 2);
            xlabel('$x$ (m)', 'interpreter', 'latex')
            ylabel("$u(x,t)$", 'interpreter', 'latex')

        end
        ylim([-1, 1])

    %     grid on;
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
        set(gcf, 'color', 'w')
    
        j = 1;
        for i = 0:5:N
            if i ~= 0 && i ~= N
                curv(j) = (u(i+1) - 2 * u(i) + u(i-1)) * 200;
                if curv(j) ~= 0
                    color = [curv(j) + 0.5, 1-(abs(curv(j)) + 0.5), 1-(curv(j)+0.5)];
                    myArrow([i / N, i / N], [0, curv(j)*1] + u(i), 10 * abs(curv(j)), 20 * abs(curv(j)), 20 * abs(curv(j)), color);
                end
                j = j+1;
            end
        end

        drawnow;
        if recordVid
            M(frame) = getframe(gcf);
            frame = frame + 1;
        end
%         pause(0.1);
    end
    uPrev = u;
    u = uNext;

end
if recordVid
    M = M(1:n-1);
    v = VideoWriter('curvature1DWave.mp4', 'MPEG-4');
    open(v)
    writeVideo(v, M);
    close(v)
end