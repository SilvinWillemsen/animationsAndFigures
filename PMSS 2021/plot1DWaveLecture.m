clear all;
close all;

figure('Position', [0, 300, 600, 200])
fs = 44100;
k = 1/fs;

drawSpeed = 1;

plotDiscrete = false;
plotBoth = false;
if plotDiscrete
    lengthSound = 60;
    c = 1470;
    widthU = 7;
else
    lengthSound = 600;
    c = 1470/10;
    widthU = 101;
end
h = c*k;
L = 1;
N = floor(L/h);
% h = L/N;

lambdaSq = (1*c*k/h)^2;

uNext = zeros(N+1, 1);
u = zeros(N+1, 1);
wNext = zeros(N+1, 1);
w = zeros(N+1, 1);

widthW = 11;

halfWidthU = floor(widthU/2);
halfWidthW = floor(widthW/2);

u(floor(N/2)-halfWidthU+1 : floor(N/2) + halfWidthU+1) = hann(widthU); 
uPrev = u;

w(floor(N/2)-halfWidthW : floor(N/2) + halfWidthW) = hann(widthW); 
wPrev = w;

range = 2:N;

recordVid = true;

if recordVid
    loops = 300;
    M(loops) = struct('cdata',[],'colormap',[]);
    frame = 1;
end

for n = 1:lengthSound
    uNext(range) = 2*u(range) - uPrev(range) + lambdaSq * (u(range+1) - 2*u(range) + u(range-1));
    wNext(range) = 2*w(range) - wPrev(range) + lambdaSq * (w(range+1) - 2*w(range) + w(range-1));

%     uNext(end)= 2*u(end) - uPrev(end) + lambdaSq * (2 * u(end-1) - 2*u(end));
%     uNext(1)= 2*u(1) - uPrev(1) + lambdaSq * (2 * u(2) - 2*u(1));

    output(n) = u(end-3);
    %% drawthings
    if mod(n,drawSpeed) == 0
        if plotDiscrete
            plot(0:N, u, 'k', 'Marker', '.', 'Markersize', 20);
            xLab = xlabel('$l$', 'interpreter', 'latex');
            xLab.Position(2) = -1.1;
            ylabel("$q_l^n$", 'interpreter', 'latex')
            title("$n = " + num2str(n-1)+"$", 'interpreter', 'latex');
        else
            if plotBoth
                hold off;
                plot((0:N) / N, u, 'b', 'Linewidth', 2);
                hold on;
                plot((0:N) / N, w, 'r', 'Linewidth', 1);
            else
                plot((0:N) / N, u, 'k', 'Linewidth', 2);
            end
            xLab = xlabel('$x$', 'interpreter', 'latex')
            xLab.Position(2) = -1.1;
%             legend('$u$', '$w$', 'interpreter', 'latex')
            ylabel("$q$", 'interpreter', 'latex')
            title("$t = " + num2str(floor(10000*(n-1)/fs)/10, 2)+"$ ms", 'interpreter', 'latex');
            xticklabels({'$0$', '$L$'})
            xticks([0, 1])
            set(gca, 'ticklabelInterpreter', 'latex')
            
            %             title("$n = " + num2str(n-1)+"$", 'interpreter', 'latex');

        end
        ylim([-1, 1])
        xticklabels({'$0$', '$N$'})
        xticks([0, N])
        set(gca, 'ticklabelInterpreter', 'latex')

    %     grid on;
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
        set(gcf, 'color', 'w')
        axis off
        drawnow;
%         pause(0.1)
        if recordVid
            M(frame) = getframe(gcf);
            frame = frame + 1;
        end
 
    end
    if plotDiscrete && recordVid
        for tt = 1:4
            M(frame) = getframe(gcf);
            frame = frame + 1;
        end
    end
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;

end
if recordVid
    M = M(1:loops-1);
    if plotDiscrete
        v = VideoWriter('1DWaveLambda1005.mp4', 'MPEG-4');
    else
        v = VideoWriter('1dwaveDAFxCont.mp4', 'MPEG-4');
    end
    open(v)
    writeVideo(v, M);
    close(v)
end