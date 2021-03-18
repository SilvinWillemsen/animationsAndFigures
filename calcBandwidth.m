function [fftOut, lastPeak] = calcBandwidth (quality, halfWidth, plotthings, height)
    fs = 44100;
    k = 1/fs;
    lengthSound = fs*0.1;

    L = 1;
    N = 50;
    h = L / N;
    c = h/k;

    lambdaSq = (quality * c * k / h)^2;

    uNext = zeros(N+1, 1);
    u = zeros(N+1, 1);
    u(round(N/2 - halfWidth) : round(N/2 + halfWidth)) = hann(2*halfWidth + 1);
%     u(floor(N/2)) = 1;
    uPrev = u;

    range = 2:N;
    out = zeros(lengthSound, 1);
    for n = 1:lengthSound
        uNext(range) = 2 * u(range) - uPrev(range) + lambdaSq * (u(range+1) - 2 * u(range) + u(range-1));
        uPrev = u;
        u = uNext;
        out(n) = uNext(floor(N/3));
        if (plotthings > 1 && n == plotthings) || n == -1
            plot(0:N, uNext + height, 'k', 'Linewidth', 1.5);
            ylim([-0.8, 0.8])
            drawnow;
        end
    end
    fftOut = abs(fft(out));
    fftOut = fftOut(1:lengthSound/2);
    if plotthings == 1
        hold off
        plot(fftOut)
        hold on;
        plot(lengthSound / (2*pi) * asin(lambda) * lengthSound / 2, [0, max(fftOut(1:lengthSound/2))], 'k')
        drawnow;
    end
    test = find(fftOut > 10);
    % if isempty(test)
    %     disp("> 50")
    %     test = find(fftOut > 50);
    % end
    if isempty(test)
        lastPeak = [];
        disp("No peak found");
        return;
    end
    lastPeak = test(end) / (lengthSound / 2);
%     soundsc(out, fs)
    if mod(floor(quality * 1000) / 10, 10) == 0
        quality
    end
end