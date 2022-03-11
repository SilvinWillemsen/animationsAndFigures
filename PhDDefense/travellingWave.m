close all
fs = 44100;
lengthSound = 500;
k = 1/fs;
c = 300;
h = c * k;
L = 1;

N = floor (L/h);
h = L/N;
loc = 0.2;
locN = floor(loc*N+1);
offsetAmp1 = 1/locN;
offsetAmp2 = 1/(N+1-locN);
lWave = 0.5 * [linspace(offsetAmp1, 1, locN), linspace(1 - offsetAmp2, 0, N+1-locN)];
rWave = 0.5 * [linspace(0, 1, locN), linspace(1 - offsetAmp2, offsetAmp2, N+1-locN)];
% rWave = fliplr(rWave)

recordVid = true;
slowdown = 2;
if recordVid
    loops = lengthSound * slowdown;
    M(loops) = struct('cdata',[],'colormap',[]);
    frame = 1;
end

for n = 1:lengthSound   
%     hold off
%     plot(lWave-1.5, '.', 'color', [0.75, 0.75, 0.75], 'Linewidth', 0.25);
%     hold on;
%     plot(rWave-1.5, '.', 'color', [0.75, 0.75, 0.75], 'Linewidth', 0.25);
    plot(lWave+rWave, 'k', 'Linewidth', 2)
    ylim([-2.5, 1])
    axis off;
    set(gcf, 'color', 'w')
    drawnow;
    if recordVid
        for i = 1:slowdown
            M(frame) = getframe (gcf);
            frame = frame + 1;
        end
    end
    
    ends = [lWave(1), rWave(end)];
    lWave = [lWave(2:N+1), -ends(2)];
    rWave = [-ends(1), rWave(1:N)];

end

if recordVid
    M = M(1:loops-1);
    v = VideoWriter('travellingWavesOrig.mp4', 'MPEG-4');
    open(v)
    writeVideo(v, M);
    close(v)
end