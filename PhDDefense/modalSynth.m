close all
clear all;
fs = 44100;
lengthSound = 5000;
k = 1/fs;
c = 300;
h = c * k;
L = 1;
f0 = c / (2 * L);
P = 0.2;
amp = 1;
numModes = 10000;
m = 1:numModes;
l = (0:0.01:L)';
C = 2 * amp ./ (m.^2 * pi^2 * P * (1-P)) .* sin (pi * m * P);

plotState = true;
recordVid = true;
slowdown = 0.5;
if recordVid
    loops = lengthSound * slowdown;
    M(loops) = struct('cdata',[],'colormap',[]);
    frame = 1;
end

if plotState
    for n = 1:2:lengthSound
        u = sum(C .* cos (f0 * m * n / fs) .* sin(m .* pi .* l / L), 2);
        hold off
        plot(u, 'k', 'Linewidth', 2)
        hold on
         for mm = 1:10
            uSin = C(mm) * cos (f0 * mm * n / fs) .* sin(mm .* pi .* l / L);
            plot1 = plot(uSin, 'k', 'Linewidth', 1.5);
            plot1.Color(4) = 0.2;
            hold on
            axis off
        end
%         ylim([-1, 1])
        ylim([-2.5, 1])
        axis off
        set(gcf, 'color', 'w')
        drawnow;
%         if recordVid
%             for i = 1:slowdown
                M(frame) = getframe (gcf);
                frame = frame + 1;
%             end
%         end
    end
else
    figure('Position', [489 195 338 662])
    for n = 1:2:lengthSound
        hold off
        for mm = 1:10
            uSin = 1/(mm) * cos (f0 * mm * n / fs) .* sin(mm .* pi .* l / L);
            plot(uSin-1.5*mm, 'color', 'k', 'Linewidth', 1.5)
            hold on
            axis off
        end
        xlim([0, length(l)])
        ylim([-16, 1])
        set(gca, 'Position', [0.0296 0.0017 0.9556 0.9948])
        set(gcf, 'color', 'w')
        drawnow
        if recordVid
%             for i = 1:slowdown
                M(frame) = getframe (gcf);
                frame = frame + 1;
%             end
        end
    end
end
% 0.75 * [1-abs(C(mm)), 1-abs(C(mm)), 1-abs(C(mm))]
if recordVid
    M = M(1:loops-1);
    v = VideoWriter('modalOrig.mp4', 'MPEG-4');
    open(v)
    writeVideo(v, M);
    close(v)
end