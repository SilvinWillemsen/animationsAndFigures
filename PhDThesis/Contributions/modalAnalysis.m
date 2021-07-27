% close all;
firstIteration = true;

if firstIteration
%     clear all;
    firstIteration = true;
else
    figure('Position', [200, 200, 1200, 450])
    set(gcf, 'color', 'w');
end
figure('Position', [556 485 493 372]);
loopingN = false;
loopNStart = 15; % also use for Ninit an Nend
loopNend = 20;
plotModeShapesBool = ~firstIteration;

lowPassConnection = false;
lpExponent = 30;

modeToPlot = 8; % if -1 plot all modes
limSubplots = 3;
if loopingN
    range = loopNStart:0.01:loopNend;
else
    range = 1;
end
if plotModeShapesBool
    loopAmount = 100;
elseif loopingN
    loopAmount = 1000;
else
    loopAmount = 1000;
end
alfExpStart = 6;
alfExpEnd = 10;
alfExp = alfExpStart:(alfExpEnd-alfExpStart)/loopAmount:alfExpEnd;
alfExp = ones(loopAmount,1) * 100;
if firstIteration
    loopAmountRange = 1:loopAmount;
else
    loopAmountRange = 1:(loopAmount+1);
end

% choose interpolation
interpolation = "quadratic";

plotMulti = false;

%{
    Number from the right boundary (quite important, switches between
    different techniques)
    -1: Adding to the center alternating between left and right string.
    0: Interpolated boundary
    1: Right string has a single moving point. Using simply supported boundary condition
    2: Right string has two moving points. When trying to solve the cubic
    interpolation, w_2 is always 0 (that's why this is a bit different)
    >3: (Expected behaviour) Selects where to add points (to left string).
%}

numFromBound = -1;

% sinc settings
even = false; % for sinc interpolation only
shifted = false; % for odd sinc interpolation only
overrideSincWidth = 2;

%{
Decides whether to use the full range or not 
    0: fullSinc is false 
    1: include all moving points
    2: include boundaries as well
    3: include virtual grid points
%}
fullSinc = 0;
fSCentered = true; % "true" only works if numFromBound == -1

for Nloop = range
    fs = 44100;                 % sample rate
    k = 1/fs;                       % time step
    
    if loopingN
        Ninit = Nloop;          % number of intervals (!)
        Nend = Nloop;
    else
        Ninit = loopNStart;
        Nend = loopNend;
    end
    
    h = 1/Ninit;                % grid spacing
    c = h/k;                    % wave speed
    N = floor(1/h);             % recalculate (should be exactly equal to Ninit)
    NinitSave = Ninit;
    
    cInit = h/k;
    cEnd = 1/(Nend * k);
    cVec = linspace (cInit, cEnd, loopAmount + 1);
    
    if interpolation == "none"
        h = 1/N;
    end
    
    lambdaSq = (cInit * k / h)^2;
    
    
    % Create B-matrix
    BFull = zeros(N, N); % matrix for scheme without boundaries, plus one overlap
    
    
    if numFromBound == -1 % add alternatively to both in the center
        M = ceil(N/2);
        Mw = floor(N/2);
    else
        M = N-numFromBound;
        Mw = numFromBound;
    end
    
    if interpolation == "none"
        Dxxu = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));
        BFull(1:N, 1:N) = 2 * eye(N) + lambdaSq * Dxxu;
    else
        Dxxu = (sparse(2:M, 1:M-1, ones(1, M-1), M, M) + ...
            sparse(1:M, 1:M, -2 * ones(1, M), M, M) + ...
            sparse(1:M-1, 2:M, ones(1, M-1), M, M));
        BFull(1:M, 1:M) = 2 * eye(M) + lambdaSq * Dxxu;

        Dxxw = (sparse(2:Mw, 1:Mw-1, ones(1, Mw-1), Mw, Mw) + ...
            sparse(1:Mw, 1:Mw, -2 * ones(1, Mw), Mw, Mw) + ...
            sparse(1:Mw-1, 2:Mw, ones(1, Mw-1), Mw, Mw));

        BFull((M+1):end, (M+1):end) = 2 * eye(Mw) + lambdaSq * Dxxw;

        BFullInit = BFull;
    end
    maxNumberOfPoints = max(Ninit, Nend);
    if firstIteration
        modesSave = zeros(loopAmount, floor(maxNumberOfPoints));
    end
    NPrev = N;
    j = 1;
    if firstIteration
        loopStart = zeros(ceil(abs(Nend-Ninit)), 1);
        loopStart(j) = 1;
    end
    addLastPoint = true;
    
    for i = loopAmountRange
        h = cVec(i) * k;
        Ninit = 1/h;
        N = floor(Ninit);
        alf  = Ninit - N;
        alfSave(i) = alf;
        if interpolation == "none"
            h = 1/N;
        end
        
        lambdaSq = (cVec(i) * k / h)^2;
        
        if N ~= NPrev
%             modeToPlot = modeToPlot + 1;
            BFull = zeros(N, N);
            if i~= 1 && firstIteration
                j = j + 1;
                loopStart(j) = i;
            end
            if i == loopAmount
                addLastPoint = false;
            end
            
            if numFromBound == -1 % add alternatively to both in the center
                M = ceil(N/2);
                Mw = floor(N/2);
            else
                M = N-numFromBound;
                Mw = numFromBound;
            end
            
            Dxxu = (sparse(2:M, 1:M-1, ones(1, M-1), M, M) + ...
                sparse(1:M, 1:M, -2 * ones(1, M), M, M) + ...
                sparse(1:M-1, 2:M, ones(1, M-1), M, M));
            BFull(1:M, 1:M) = 2 * eye(M) + lambdaSq * Dxxu;
            
            Dxxw = (sparse(2:Mw, 1:Mw-1, ones(1, Mw-1), Mw, Mw) + ...
                sparse(1:Mw, 1:Mw, -2 * ones(1, Mw), Mw, Mw) + ...
                sparse(1:Mw-1, 2:Mw, ones(1, Mw-1), Mw, Mw));
            
            BFull((M+1):end, (M+1):end) = 2 * eye(Mw) + lambdaSq * Dxxw;
            
            BFullInit = BFull;
        end
        NPrev = N;
        
        if interpolation == "none"
            
            Dxxu = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));
            BFull(1:N, 1:N) = 2 * eye(N) + lambdaSq * Dxxu;
            
        elseif interpolation == "linear"
            if numFromBound == 1
                BFull(M, M+1) = alf;
                BFull(M + 1, M-1 : M) = [(1-alf), alf];
            else
                BFull(M, (M+1):(M + 2)) = [alf, (1-alf)];
                BFull(M + 1, (M-1) : M) = [(1-alf), alf];
            end
        elseif interpolation == "quadratic"
            ip = [(alf - 1)/(alf + 1), ...
                  1, ...
                  -(alf - 1)/(alf + 1)];
            if numFromBound == 1
                BFull(M, (M):(M + 1)) = BFullInit(M, (M):(M + 1)) + ip(1:2);
                BFull(M+1, (M-1):(M+1)) = BFullInit(M+1, (M-1):(M+1)) + fliplr(ip);
            else
                BFull(M, (M):(M + 2)) = BFullInit(M, (M):(M + 2)) + ip;
                BFull(M+1, (M-1):(M+1)) = BFullInit(M+1, (M-1):(M+1)) + fliplr(ip);
            end
        elseif interpolation == "cubic"
            ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
            Ainv = inv([1, -ip(4); -ip(4), 1]);
            if numFromBound == 2
                BFull(M, (M + 1):(M + 2)) = ip(3:-1:2) * Ainv(1, 1);
                BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + ip(1:3) * Ainv(1, 2);
                BFull(M+1, (M + 1):(M + 2)) = BFullInit(M+1, (M + 1):(M + 2)) + ip(3:-1:2) * Ainv(2, 1);
                BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
            elseif numFromBound == 1
                BFull(M, M+1) = (ip(3) - ip(1)) * Ainv(1, 1);
                BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + ip(1:3) * Ainv(1, 2);
                BFull(M+1, M+1) = BFullInit(M+1, M+1) + (ip(3) - ip(1)) * Ainv(2, 1);
                BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
            elseif numFromBound == 0
                hLocs = 1:-h:0;
                if hLocs(end) >= h/2
                    alf = (h-hLocs(end))/ hLocs(end);
                    %                     uVirtual = -alf * u(1);
                    BFull(N, N) = BFullInit(N, N) - alf;
                    
                else
                    alf = (2*hLocs(end)) / h;
                    %                     uVirtual = -(alf * u(1) + (1-alf) * u(2));
                    BFull(N, N-1:N) = BFullInit(N, N-1:N) - [(1-alf), alf];
                end
                %                 BFull(N, N) = -1;
            else
                BFull(M, (M + 1):(M + 3)) = ip(3:-1:1) * Ainv(1, 1);
                BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + ip(1:3) * Ainv(1, 2);
                BFull(M+1, (M + 1):(M + 3)) = BFullInit(M+1, (M + 1):(M + 3)) + ip(3:-1:1) * Ainv(2, 1);
                BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
            end
        elseif interpolation == "altCubic"
            ip = [(alf - 1) / (alf + 2), ...
                 (alf + 1) / 2, ...
                1-alf, ...
                alf * (alf - 1) / (2 * (alf + 2))];
            if numFromBound == 2
                BFull(M, M:(M + 2)) = BFullInit(M, M:(M + 2)) + ip(1:3);
                BFull(M+1, (M-2):(M+1)) = BFullInit(M+1, (M-2):(M+1)) + fliplr(ip);
            elseif numFromBound == 1
                BFull(M, M:(M + 1)) = BFullInit(M, M:(M + 1)) + ip(1:2) - [0, ip(4)];
                BFull(M+1, (M-2):(M+1)) = BFullInit(M+1, (M-2):(M+1)) + fliplr(ip);
            elseif numFromBound == 0
                disp ("undefined");
            else
                BFull(M, M:(M + 3)) = BFullInit(M, M:(M + 3)) + ip;
                BFull(M+1, (M-2):(M+1)) = BFullInit(M+1, (M-2):(M+1)) + fliplr(ip);
            end 
%             BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
        elseif interpolation == "shiftedCubic"
            ip  = [-(alf*(alf - 1))/((alf + 1)*(alf + 2)), ...
                   (2*(alf - 1))/(alf + 1), ...
                   2/(alf + 1), ...
                   -(2*(alf - 1))/((alf + 1)*(alf + 2))];
            if numFromBound == 1
                BFull(M, M-1:(M + 1)) = BFullInit(M, M-1:(M + 1)) + ip(1:3);
                BFull(M+1, (M-1):(M+1)) = BFullInit(M+1, (M-1):(M+1)) + fliplr(ip(2:4));
            else
                BFull(M, M-1:(M + 2)) = BFullInit(M, M-1:(M + 2)) + ip;
                BFull(M+1, (M-1):(M+2)) = BFullInit(M+1, (M-1):(M+2)) + fliplr(ip);
            end
        elseif interpolation == "quartic"
            ip = [-(alf*(alf - 1))/((alf + 2)*(alf + 3)), ...
                (2*(alf - 1))/(alf + 2), ...
                1, ...
                -(2*(alf - 1))/(alf + 2), ...
                (alf*(alf - 1))/((alf + 2)*(alf + 3))];
            if numFromBound == 2
                BFull(M, (M-1):(M + 2)) = BFullInit(M, (M-1):(M + 2)) + ip(1:4);
                BFull(M+1, (M-2):(M+2)) = BFullInit(M+1, (M-2):(M+2)) + fliplr(ip);
            elseif numFromBound == 1
                BFull(M, (M-1):(M + 1)) = BFullInit(M, (M-1):(M + 1)) + ip(1:3) - [0, 0, ip(5)];
                BFull(M+1, (M-2):(M+1)) = BFullInit(M+1, (M-2):(M+1)) + fliplr(ip(2:5));
            else
                BFull(M, (M-1):(M + 3)) = BFullInit(M, (M-1):(M + 3)) + ip;
                BFull(M+1, (M-2):(M+2)) = BFullInit(M+1, (M-2):(M+2)) + fliplr(ip);
            end
        elseif interpolation == "six"
            ip = [(alf*(alf - 1)*(alf + 1))/((alf + 3)*(alf + 4)*(alf + 5)), ...
                  -(3*alf*(alf - 1))/((alf + 3)*(alf + 4)), ...
                  (3*(alf - 1))/(alf + 3), ...
                  1, ...
                  -(3*(alf - 1))/(alf + 3), ...
                  (3*alf*(alf - 1))/((alf + 3)*(alf + 4)), ...
                  -(alf*(alf - 1)*(alf + 1))/((alf + 3)*(alf + 4)*(alf + 5))];
             BFull(M, (M-2):(M + 4)) = BFullInit(M, (M-2):(M + 4)) + ip;
             BFull(M+1, (M-3):(M+3)) = BFullInit(M+1, (M-3):(M+3)) + fliplr(ip);

        elseif interpolation == "highEven"
             ip =  [-(alf*(alf - 1)*(alf + 1)*(alf + 2)*(alf + 3)*(alf + 4))/((alf + 6)*(alf + 7)*(alf + 8)*(alf + 9)*(alf + 10)*(alf + 11)), ...
                     (6*alf*(alf - 1)*(alf + 1)*(alf + 2)*(alf + 3))/((alf + 6)*(alf + 7)*(alf + 8)*(alf + 9)*(alf + 10)), ...
                                        -(15*alf*(alf - 1)*(alf + 1)*(alf + 2))/((alf + 6)*(alf + 7)*(alf + 8)*(alf + 9)), ...
                                                             (20*alf*(alf - 1)*(alf + 1))/((alf + 6)*(alf + 7)*(alf + 8)), ...
                                                                                -(15*alf*(alf - 1))/((alf + 6)*(alf + 7)), ...
                                                                                                  (6*(alf - 1))/(alf + 6), ...
                                                                                                                        1, ...
                                                                                                 -(6*(alf - 1))/(alf + 6), ...
                                                                                 (15*alf*(alf - 1))/((alf + 6)*(alf + 7)), ...
                                                            -(20*alf*(alf - 1)*(alf + 1))/((alf + 6)*(alf + 7)*(alf + 8)), ...
                                         (15*alf*(alf - 1)*(alf + 1)*(alf + 2))/((alf + 6)*(alf + 7)*(alf + 8)*(alf + 9)), ...
                    -(6*alf*(alf - 1)*(alf + 1)*(alf + 2)*(alf + 3))/((alf + 6)*(alf + 7)*(alf + 8)*(alf + 9)*(alf + 10)), ...
  (alf*(alf - 1)*(alf + 1)*(alf + 2)*(alf + 3)*(alf + 4))/((alf + 6)*(alf + 7)*(alf + 8)*(alf + 9)*(alf + 10)*(alf + 11))];
            ipRange = M + (-(length(ip)-1)/2:(length(ip)-1)/2);
            BFull(M, ipRange + 1) = BFullInit(M, ipRange + 1) + ip;
            BFull(M+1, ipRange) = BFullInit(M+1, ipRange) + fliplr(ip);
        elseif interpolation == "sinc"
            includeUMp1AndWm1 = true;
            
            if alf < 1e-6
                alf = alf + 1e-6;
            end
            
            alphaBand = 0.8; % relative bandwidth range
            bmax = alphaBand*pi;
            
            if numFromBound == -1
                sincWidth = floor(N / 2) - 1;
                sincWidth = 2;
            else
%                 sincWidth = numFromBound+1;
            end
            if overrideSincWidth ~= 0
                sincWidth = overrideSincWidth; % override sincwidth
            end
            if fullSinc == 0
                if includeUMp1AndWm1
                    if even
                        xUMp1 = [-sincWidth:-1, -1:sincWidth-1]';
                        xUMp1 = xUMp1 + [zeros(sincWidth, 1); alf * ones(sincWidth+1, 1)];
                    else
                        if shifted
                            xUMp1 = [-sincWidth:-1, -1:sincWidth-2]';
                            xUMp1 = xUMp1 + [zeros(sincWidth, 1); alf * ones(sincWidth, 1)];
                        else
                            xUMp1 = [-sincWidth+1:-1, -1:sincWidth-1]';
                            xUMp1 = xUMp1 + [zeros(sincWidth-1, 1); alf * ones(sincWidth+1, 1)];
                        end
                    end
                else
                    xUMp1 = (-sincWidth:sincWidth-1)';
                    xUMp1 = xUMp1 + [zeros(sincWidth, 1); alf * ones(sincWidth, 1)];
                end
            elseif fullSinc == 1
                xUMp1 = (1:N)' - M - 1;
                xUMp1 = xUMp1 + [zeros(M, 1); alf * ones(N-M, 1) - 1];
            elseif fullSinc == 2
                xUMp1 = (0:N+1)' - M - 1;
                xUMp1 = xUMp1 + [zeros(M+1, 1); alf * ones(N-M+1, 1) - 1];
            elseif fullSinc == 3
                xUMp1 = (-1:N+2)' - M - 1;
                xUMp1 = xUMp1 + [zeros(M+2, 1); alf * ones(N-M+2, 1) - 1];
%                 xUMp1 = (-1:N+2)' - M - 2;
%                 xUMp1 = xUMp1 + [zeros(M+3, 1); alf * ones(N-M+1, 1) - 1];
            end
            if numFromBound == -1 && fSCentered && fullSinc ~= 0
                xUMp1 = xUMp1(2:end);
            end
            
            iLen = length (xUMp1); % length of interpolation (N in stefans implementation)
            bU = (sin(bmax*xUMp1)./xUMp1);
            if sum(isnan(bU))
                idxIsNan = find(isnan(bU));
%                 bU(idxIsNan-2) = -bmax;
%                 bU(idxIsNan-1) = -bmax;
                bU(idxIsNan) = bmax;

            end 
            distU = xUMp1*ones(1,iLen)-ones(iLen,1)*xUMp1';    % distance matrix between points
            AU = sin(bmax*distU)./distU;
            AU(1+(iLen+1)*[0:iLen-1]') = bmax;         % collection of sinc functions with centers at grid point locations
%             AU(isnan(AU)) = bmax;
            aU = AU\bU; %optimal coefficients
%             if alf == 0
%                 idxAu = find(round(aU) == 1);
%                 aU(idxAu-2) = -1;
%                 aU(idxAu-1) = 1;
%                 aU(idxAu) = 1;
%             end
%             aU = aU / sum(aU);
%             sum(aU)
            if fullSinc == 0
                aW = flipud(aU);

%                 if includeUMp1AndWm1
%                     if even
%                         xWm1 = [-sincWidth+1:1, 1:sincWidth]';
%                         xWm1 = xWm1 - [alf * ones(sincWidth+1, 1); zeros(sincWidth, 1)];
%                     else
%                         xWm1 = [-sincWidth+1:1, 1:sincWidth-1]';
%                         xWm1 = xWm1 - [alf * ones(sincWidth, 1); zeros(sincWidth, 1)];
%                     end
%                 else
%                     xWm1 = (-sincWidth+1:sincWidth)';
%                     xWm1 = xWm1 - [alf * ones(sincWidth, 1); zeros(sincWidth, 1)];
%                 end
            elseif fullSinc == 1
                xWm1 = (1:N)' - M;
                xWm1 = xWm1 - [alf * ones(M, 1) - 1; zeros(N-M, 1)];
            elseif fullSinc == 2
                xWm1 = (0:N+1)' - M;
                xWm1 = xWm1 - [alf * ones(M+1, 1) - 1; zeros(N-M+1, 1)];

            elseif fullSinc == 3
                xWm1 = (-1:N+2)' - M;
                xWm1 = xWm1 - [alf * ones(M+2, 1) - 1; zeros(N-M+2, 1)];
%                 xWm1 = (-1:N+2)' - M;
%                 xWm1 = xWm1 - [alf * ones(M+3, 1); ones(N-M+1, 1)];
            end
            if numFromBound == -1 && fSCentered && fullSinc ~= 0
                xWm1 = xWm1(1:end-1);
            end
            if fullSinc ~= 0
                bW = (sin(bmax*xWm1)./xWm1);
                if sum(isnan(bW))
                    bW(isnan(bW)) = bmax;
                end 
                distW = xWm1*ones(1,iLen)-ones(iLen,1)*xWm1';    % distance matrix between points
                AW = sin(bmax*distW)./distW;
                AW(1+(iLen+1)*[0:iLen-1]') = bmax;         % collection of sinc functions with centers at grid point locations
                AW(isnan(AW)) = bmax;
                aW = AW\bW; %optimal coefficients
            end
            
            
%             if mod(i, 10) == 0
%                 hold off;
%                 plot(xUMp1 + M + 1, bU)
%                 hold on;
%                 uRange = min(xUMp1):0.001:max(xUMp1);
%                 plot (uRange + M + 1, sin(bmax*uRange)./(uRange));
%                 
%                 plot(xWm1 + M - 1, bW)
%                 wRange = min(xWm1):0.001:max(xWm1);
%                 plot (wRange + M - 1, sin(bmax*wRange)./(wRange));
%                 drawnow;
%             end
            
            if numFromBound == -1 && fSCentered
                inputRange1 = 2:N;
                inputRange2 = 1:N-1;

            else
                inputRange1 = 1:N;
                inputRange2 = 1:N;

            end
            if fullSinc == 1
                BFull(M, inputRange1) = BFullInit(M, inputRange1) + aU';
                BFull(M+1, inputRange2) = BFullInit(M+1, inputRange2) + aW';
            elseif fullSinc == 2
                BFull(M, inputRange1) = BFullInit(M, inputRange1) + aU(2:end-1)';
                BFull(M+1, inputRange2) = BFullInit(M+1, inputRange2) + aW(2:end-1)';
    
            elseif fullSinc == 3
                BFull(M, inputRange1) = BFullInit(M, inputRange1) + aU(3:end-2)' - [aU(1), zeros(1, length(aU)-6), aU(end)];
                BFull(M+1, inputRange2) = BFullInit(M+1, inputRange2) + aW(3:end-2)' - [aW(1), zeros(1, length(aW)-6), aW(end)];
            elseif M + sincWidth == N+1 % use boundary condition
                if includeUMp1AndWm1
                    if even
                        BFull(M, M-sincWidth+1:end) = BFullInit(M, M-sincWidth+1:end) + aU(1:end-2)' - [zeros(1, length(aU)-3), aU(end)];
                        BFull(M+1, (M-sincWidth : end)) = BFullInit(M+1, (M-sincWidth : end)) + aW(1:end-1)';
                    else
                        BFull(M, M-sincWidth+2:end) = BFullInit(M, M-sincWidth+2:end) + aU(1:end-2)' - [zeros(1, length(aU)-3), aU(end)];
                        BFull(M+1, (M-sincWidth : end)) = BFullInit(M+1, (M-sincWidth : end)) + aW';
                    end
                else
                    BFull(M, M-sincWidth+2:end) = BFullInit(M, M-sincWidth+2:end) + aU(1:end-2)' - [zeros(1, length(aU)-3), aU(end)];
                    BFull(M+1, (M-sincWidth: end-1)) = BFullInit(M+1, (M-sincWidth: end-1)) + aW(1:end-1)';
                end
            elseif M + sincWidth == N % use boundary condition
                if includeUMp1AndWm1
                    if even
                        BFull(M, M-sincWidth+1 : end) =  BFullInit(M, M-sincWidth+1 : end) + aU(1:end-1)';
                        BFull(M+1, M-sincWidth : end) = BFullInit(M+1, M-sincWidth : end) + aW';
                    else
                        BFull(M, M-sincWidth+2 : end) =  BFullInit(M, M-sincWidth+2 : end) + aU(1:end-1)';
                        BFull(M+1, M-sincWidth : end-1) = BFullInit(M+1, M-sincWidth : end-1) + aW';
                    end     
%                     BFull(M, M-sincWidth+1:end) = BFullInit(M, M-sincWidth+1:end) + aU(1:end-1)';
%                     BFull(M+1, (M-sincWidth : end)) = BFullInit(M+1, (M-sincWidth : end)) + aW';
                else
                    BFull(M, M-sincWidth+2:end) = BFullInit(M, M-sincWidth+2:end) + aU(1:end-1)';
                    BFull(M+1, (M-sincWidth: end-1)) = BFullInit(M+1, (M-sincWidth: end-1)) + aW';
                end
            else
%                 BFull(M, (M + 1):(M + 3)) = aU(1:3);
                if includeUMp1AndWm1
                    if even
                        sincRange = (M-sincWidth : M+sincWidth);
                        BFull(M, sincRange+1) =  BFullInit(M, sincRange+1) + aU';
                        BFull(M+1, sincRange) = BFullInit(M+1, sincRange) + aW';
                    else
                        sincRange = (M-sincWidth+1 : M+sincWidth);

                        if shifted
                            BFull(M, sincRange) =  BFullInit(M, sincRange) + aU';
                            BFull(M+1, sincRange) = BFullInit(M+1, sincRange) + aW';
                        else
                            BFull(M, sincRange+1) =  BFullInit(M, sincRange+1) + aU';
                            BFull(M+1, sincRange-1) = BFullInit(M+1, sincRange-1) + aW';
                        end
                    end                    
                else
                    sincRange =[M-sincWidth:M-1, M+1:M+sincWidth];
                    BFull(M, sincRange+1) =  BFullInit(M, sincRange+1) + aU';
                    BFull(M+1, sincRange) = BFullInit(M+1, sincRange) + aW';
                end
%                 if  mod(i, 100) == 0 && i>0 
%                     hold off   
%                     plot(aU)
%                     hold on;
%                     plot(aW)
% %                     pause(0.5)
%                     title(alf)
%                     ylim([-1,1])
%                     drawnow;
%                 end
%                 BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
            end
%             if mod(i, 1) == 0
%                 subplot(2,1,1)
%                 plot (aU)
%                 subplot(2,1,2)
%                 plot (aW)
%                 drawnow;
%             end
        elseif interpolation == "altSinc"
            alphaBand = 1; % relative bandwidth range
            bmax = alphaBand*pi / h;
            
            xIpLocs = (-alf + (-1:2)) * h;
            sincIp = sin(bmax * xIpLocs) ./ (bmax * xIpLocs);
            sincIp(isnan(sincIp)) = 1;
%             sincIp = fliplr(sincIp);
            Ainv = inv([1, -sincIp(4); -sincIp(4), 1]);
            
            if numFromBound == 1
                BFull(M, M+1) = (sincIp(3) - sincIp(1)) * Ainv(1, 1);
                BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + sincIp(1:3) * Ainv(1, 2);
                BFull(M+1, M+1) = BFullInit(M+1, M+1) + (sincIp(3) - sincIp(1)) * Ainv(2, 1);
                BFull(M+1, (M-2 : M)) =  sincIp(1:3) * Ainv(2, 2);
                
            else
                
                BFull(M, (M + 1):(M + 3)) = sincIp(3:-1:1) * Ainv(1, 1);
                BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + sincIp(1:3) * Ainv(1, 2);
                BFull(M+1, (M + 1):(M + 3)) = BFullInit(M+1, (M + 1):(M + 3)) + sincIp(3:-1:1) * Ainv(2, 1);
                BFull(M+1, (M-2 : M)) =  sincIp(1:3) * Ainv(2, 2);
            end            
            if mod(i, 100) == 0
                hold off;
                plot (xIpLocs, sincIp)
                hold on;
                plotRange = (min(xIpLocs) : 0.001 : max(xIpLocs));
                plot (plotRange, sin(bmax * plotRange) ./ (bmax * plotRange))
                title(alf)
                drawnow;
            end
        end
        if lowPassConnection
%                 lpVec = 0.5 * [-(1-alf)^(lpExponent), (1-alf)^(lpExponent)];
            lpVec = 0.5 * [-cos((alf) * pi/2)^lpExponent, cos((alf) * pi/2)^lpExponent];
            BFull(M, M:M+1) = BFull(M, M:M+1) + lpVec;
            BFull(M+1, M:M+1) = BFull(M+1, M:M+1) - lpVec;
        end
        % imagesc(BFull)
        % drawnow;          
        [~, D, W] = eig(BFull, 'vector');
        if plotModeShapesBool
    %         order = Ninit:-1:1;
            plotModeShapes;
        end
        if firstIteration
            modesSave(i, 1:N) = sort(1/(2 * pi * k) * acos (1/2 * D));
        end
    end
    if addLastPoint && firstIteration
        loopStart(j+1) = loopAmount;
    end
    
    if firstIteration
        hold off;
        plotModesSave;
    end

    drawnow;
end
% 
% for i = 1:10:loopAmount
%     plot((diff(real(modesSave(i, ~isnan(modesSave(i, :)))))));
%     ylim([0.75 * real(modesSave(i, 1)) , real(modesSave(i, 1))]);
%     title(alfExp(i));
%     drawnow;
% end

% xData = 1:loopAmount;
% goodness = zeros(floor(loopNStart), 1);
% rms = zeros(floor(loopNStart), 1);
% sse = zeros(floor(loopNStart), 1);
% 
% figure;
% for i = 1:loopNStart
%     [~, got] = fit (xData(~isnan(modesSave(:,i)))', modesSave(~isnan(modesSave(:,i)),i), 'poly1');
%     goodness(i) = got.adjrsquare;
%     rms(i) = got.rmse;
%     sse(i) = got.sse;
%     
%     %     got
%     %     eval(['mode', num2str(i), ' = modesSave(:,', num2str(i), ');']);
% end
% subplot(211)
% plot(goodness)
% subplot(212)
% plot(sse)


