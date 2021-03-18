close all;
clear all;

plotSubs = true; % if true it makes subplots, otherwise it animates
projIsPerspective = false; %perspective or not
vertical = true;

% if plotSubs
    if vertical
        width = 400;
        figure1 = figure('Position', [0, 0, width - 20, ((width) / 4 + 20) * 1.35/0.18 ]);
    else
        height = 400;
        figure1 = figure('Position', [0, 0, height / 4 * 1.35/0.18, height - 20]);
    end
% else
%     figure1 = figure('Position', [0, 700, 1400, 300]);
% end
trombaString = load ('/Users/SilvinW/repositories/TrombaMarina/Builds/MacOSX/build/Debug/stringState.csv');
bridge = load ('/Users/SilvinW/repositories/TrombaMarina/Builds/MacOSX/build/Debug/bridgeState.csv');
body = load ('/Users/SilvinW/repositories/TrombaMarina/Builds/MacOSX/build/Debug/bodyState.csv');
sampleNumber = load ('/Users/SilvinW/repositories/TrombaMarina/Builds/MacOSX/build/Debug/sampleNumber.csv');

% load matfiles
% load trombaString.mat
% load bridge.mat
% load body.mat
% load sampleNumber.mat



N = length(trombaString(1, :));
Ny = size(body, 1) / length(sampleNumber);
Nx = size(body, 2);

connLocString = 196 - 2;
bridgeLocPlateX = 107 + 1;
bridgeLocPlateY = floor(0.135 / 0.18 * Ny)  + 1;


% mean value of string state
stringScaling = 30;
testMean = mean(mean(((ones(size(trombaString))) * 0.5 * 0.18 + (-trombaString + 5e-6) * stringScaling)));

scaling = 10; % scaling for platestate

% use different range if plotSubs == true
if plotSubs
    j = 1;
    range = 18:20:90;
else
    range = 1:size(trombaString, 1)-1;
end

%% loop
for i = range
    %subplots
    if plotSubs
        if vertical
            test = subplot(1, 4, j);
        else
            test = subplot(4, 1, j);
        end
        if projIsPerspective
            test.Position = [-0.05, 1-j*0.24, 1.09, 0.18];
        else
            if vertical
                test.Position = [(j-1)*0.25 + 0.03, 0.02,  0.18, 0.94]
            else
                test.Position = [0, 1-j*0.24, 1, 0.18];
            end
        end

        j = j + 1;
    end
    
    % plate x plot indices
    idxs = (1:Nx) / Nx * 1.35 + (1.90 - 1.35);

    [bX, bY, bZ] = sphere;
    rad = 8e-3;
    bXBridge = bX * rad + connLocString / N * 1.90 + 0.3/Nx * 1.90;
    bYBridge = bY * rad + 0.5 * 0.18 + (-bridge(i) + 5e-6) * stringScaling;
    bZBridge = bZ * rad * 0.02 + 5e-5;
    
    if projIsPerspective
        bXLeft = bX * rad + 2/N * 1.90;
        bXRight = bX * rad + 1.90 - 2/N;
    else
        bXLeft = bX * rad + 1/N * 1.90;
        bXRight = bX * rad + 1.90 - 1/N;
    end
    
    bYLeft = bY * rad + 0.5 * 0.18;
    bYRight = bYLeft;
    
    bZLeft = bZBridge * 0.02 + rad;
    bZRight = bZBridge * 0.02 + rad;
    
%     bXPlate = bX * rad + connLocString / N * 1.90;
%     bYPlate = bY * rad + bridgeLocPlateY / Ny * 0.18;
%     bZPlate = bZ * rad * 0.02 + 5e-6;

    hold off;
    % surf string using data from tubeplot (needs tubeplot.m)
    % string is 1.90 m long
    [x, y, z] = tubeplot([(1:N) / (N) * 1.90; (ones(1, N) * 0.5 * 0.18 + (-trombaString(i, :) + 5e-6) * stringScaling) ; ones(1, N) * 5e-6], 5e-3);
    surf(x, y, 0.02*z, (y * 0 - 0.8e-7) * stringScaling);%(y - testMean) * 2e-5 * stringScaling + 2e-7 * stringScaling);
    shading interp
    %colors
    colormap(gray);
%     shading interp
    
    hold on;
    
    % bridge
%     surf(bXBridge, bYBridge, bZBridge, (bYBridge - testMean) * 0.00005 * scaling + 9e-7 * scaling)
    surf(bXBridge, bYBridge, bZBridge, bYBridge*0 + 8e-3)

    %boundaries
    surf(bXLeft, bYLeft, bZLeft, bYLeft * 0 - 10)
    surf(bXRight, bYRight, bZRight, bYLeft * 0 - 10)

    
    % surf plate (plate is 0.18m wide and 1.35m long and lies under the right part of the string (from x = 0.55 until x = 1.90) )
    surf(idxs, (1:Ny) / Ny * 0.18, body(1+i*Ny : (i+1) * Ny, :)*scaling, 'Facealpha', 0.7);
    
    % surf bridgeloc on plate
%     surf(bXPlate, bYPlate, bZPlate);%, (bZPlate - testMean) * 0.00005 * scaling + 8e-7 * scaling)

    xScaling = 8e-4;
    [x1, y1, z1] = tubeplot([-10:10; 10:-1:-10 ; zeros(1, 21)] * xScaling, 2e-3);
    [x2, y2, z2] = tubeplot([-10:10; -10:10 ; zeros(1, 21)] * xScaling, 2e-3);
%     figure;
    x1 = x1 + connLocString / N * 1.90 + 0.3/Nx * 1.90;
    x2 = x2 + connLocString / N * 1.90 + 0.3/Nx * 1.90;
    
    y1 = y1 + bridgeLocPlateY / Ny * 0.18;
    y2 = y2 + bridgeLocPlateY / Ny * 0.18;
    
    z1 = z1 + 5e-4;
    z2 = z2 + 5e-4;
    surf(x1, y1, z1);
    surf(x2, y2, z2);

    % title
    title("t = " + num2str(round(sampleNumber (i) / 44.100)) + " [ms]", 'Fontsize', 15);

    % plot limits 
    xlim([0, 1.91])
    ylim([0, 0.18])
    zlim([-5e-3, 0.02 * 5])
    
    % color
    colormap(gray);
    caxis([-1e-6 1e-6] * scaling)
    shading interp

    % view
%     view([6.43279145378544 40.6120183749305])
    view([90, 90])
    % extra settings
    axis off;
    if projIsPerspective
        set(gca, 'Projection', 'perspective');%, 'Position', [-0.1, -0.1, 1, 1])
    end
    set(gcf, 'color', [0.75, 0.75, 0.75])

    if j == 2
        mArrow3([1/3 * 1.90,0.1*0.18,0.01], [1/3 * 1.90,0.55 * 0.18,0.01], ...
            'stemWidth', 2e-3, ...
            'tipWidth', 8e-3);
    end
    drawnow;
    pause(0.01);
end