clear all;
close all;
fs = 8000;
t = (1:fs)' + 543;
y = sin(2 * pi * t / fs) + sin(3 * pi * t / fs) + sin(6.5 * pi * t / fs);

N = length(t);

y = zeros(N, 1);
startPoint = floor(N/3);
endPoint = startPoint + 1000;
y(startPoint : endPoint-1) = 5 * hann(endPoint-startPoint);
uNext = zeros(N, 1);

u = y;
uPrev = y;

lambdaSq = 1;
N = N - 4;
B = sparse(2:N, 1:N-1, lambdaSq * ones(1, N-1), N, N) + ...
    sparse(1:N, 1:N, (2 - 2 * lambdaSq) * ones(1, N), N, N) + ...
    sparse(1:N-1, 2:N, lambdaSq * ones(1, N-1), N, N);
C = sparse(1:N, 1:N, -1 * ones(1, N), N, N);
N = N + 4;
figure('Renderer', 'painters', 'Position', [200 200 1000 300])
loops = 3001;
set(gcf, 'Color', 'white')
continuous = false;
bc = 'ss';
j = 0;
        
subplot(1, 2, 1)

set(gca, 'FontSize', 15, 'Color', 'white', 'Position', [0.02 0.10 0.45 0.8])
dist = 150;
axisOff = true;

for n = 1:loops
    uNext(3:N-2) = B * u(3:N-2) + C * uPrev(3:N-2);
%     uNext(172:N-172) = B(172:N-172, 172:N-172) * u(172:N-172) + C(172:N-172, 172:N-172) * uPrev(172:N-172);
    uPrev = u;
    u = uNext;
%     if mod(n,1196) == 1 && n <= 1196*2
%         subplot(1, 2, 1)
%         c = (0.5 * cos(2 * pi * n/1196) + 0.5)^5 * 0.9 + 0.1;
% %         c = abs((n-1096)/dist) / (1296/dist);
%         plot(uNext - (n-1196)/dist, 'color', [0, 0, 0, c], 'LineWidth', 2)
%         hold on;
%         s = floor((n*1000 - 1) / fs) / 1000;
%         if mod(s * 1000, 1) == 1
%             string = ['t = ', num2str(s), '0 s'];
%         else
%             string = ['t = ', num2str(s), ' s'];
%         end
%         if axisOff
%             axis off;
%         end
% %             set(gcf, 'Color', 'w')
%         xlim([0, 9000])
%         ylim([-8, 14])
%         drawnow;                        
%         continuous = false;
%     end

    if mod(n,1196) == 1
        j = j + 1;
        subplot(3, 2, j)

        plot(uNext, 'k', 'LineWidth', 2);
        if axisOff
            axis off;
        end
        set(gca, 'FontSize', 15, 'Color', 'white', 'Position', [0.05, 1-0.15*(j+1), 0.45, 0.2])
        xlim([0, 9000])
        ylim([-1, 5])
        drawnow;

    end
    if mod(n, 1196) == 1
        j = j + 1;
        subplot(3, 2, j)
        scatter(1:171:N + 171, [uNext(1:171:N); 0], 160, 'k', '.')
%         title (string, 'horizontalAlignment', 'center');
        if axisOff
            axis off;
        end
        set(gca, 'FontSize', 15, 'Color', 'white', 'Position', [0.55, 1-0.15*j, 0.45, 0.2])
        xlim([0, 9000])
        ylim([-1, 5])
        drawnow;
    end
end
