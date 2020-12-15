close all;
clear all;

plotSubs = true; % if true it makes subplots, otherwise it animates

if plotSubs
    height = 400;
    figure1 = figure('Position', [0, 0, height / 4 * 1.35/0.18, height - 20]);
else
    figure1 = figure('Position', [0, 700, 1400, 300]);
end
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
connLocString = 196 - 1;
bridgeLocPlateX = 107 + 1;
bridgeLocPlateY = 12 + 1;

Ny = size(body, 1) / length(sampleNumber);

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
        test = subplot(4, 1, j);
        test.Position = [-0.05, 1-j*0.24, 1.09, 0.18];
        j = j + 1;
    end
    
    % plate x plot indices
    idxs = (1:134) / 134 * 1.35 + (1.90 - 1.35);

    [bX, bY, bZ] = sphere;
    rad = 8e-3;
    bXBridge = bX * rad + connLocString / N * 1.90;
    bYBridge = bY * rad + 0.5 * 0.18 + (-bridge(i) + 5e-6) * stringScaling;
    bZBridge = bZ * rad + 5e-6;
    
    bXLeft = bX * rad + 2/N * 1.90;
    bXRight = bX * rad + 1.90 - 2/N;
    
    bYLeft = bY * rad + 0.5 * 0.18;
    bYRight = bYLeft;
    
    bZLeft = bZBridge * 0.02 + rad;
    bZRight = bZBridge * 0.02 + rad;
    
    bXPlate = bX * rad + connLocString / N * 1.90;
    bYPlate = bY * rad + bridgeLocPlateY / Ny * 0.18;
    bZPlate = bZ * rad * 0.02 + 5e-6;

    hold off;
    % surf string using data from tubeplot (needs tubeplot.m)
    % string is 1.90 m long
    [x, y, z] = tubeplot([(1:N) / (N) * 1.90; (ones(1, N) * 0.5 * 0.18 + (-trombaString(i, :) + 5e-6) * stringScaling) ; ones(1, N) * 5e-6], 5e-3);
    surf(x, y, 0.02*z, (y - testMean) * 2e-5 * stringScaling + 2e-7 * stringScaling);
   
    
    %colors
    colormap(gray);
%     shading interp
    
    hold on;
    
    % bridge
%     surf(bXBridge, bYBridge, bZBridge, (bYBridge - testMean) * 0.00005 * scaling + 9e-7 * scaling)
    surf(bXBridge, bYBridge, bZBridge, bYBridge*0 + 8e-6)

    %boundaries
    surf(bXLeft, bYLeft, bZLeft, bYLeft * 0 - 10)
    surf(bXRight, bYRight, bZRight, bYLeft * 0 - 10)

    
    % surf plate (plate is 0.18m wide and 1.35m long and lies under the right part of the string (from x = 0.55 until x = 1.90) )
    surf(idxs, (1:17) / 17 * 0.18, body(1+i*Ny : (i+1) * Ny, :)*scaling, 'Facealpha', 0.7);
    
    % surf bridgeloc on plate
    surf(bXPlate, bYPlate, bZPlate);%, (bZPlate - testMean) * 0.00005 * scaling + 8e-7 * scaling)

    % title
    title("t = " + num2str(round(sampleNumber (i) / 44.100)) + " [ms]", 'Fontsize', 16);

    % plot limits 
    xlim([0, 1.90])
    ylim([0, 0.18])
    zlim([-5e-3, 0.02 * 5])
    
    % color
    colormap(gray);
    caxis([-1e-6 1e-6] * scaling)
    shading interp

    % view
%     view([6.43279145378544 40.6120183749305])
    view([0, 90])
    % extra settings
    axis off;
    set(gca, 'Projection', 'perspective');%, 'Position', [-0.1, -0.1, 1, 1])
    set(gcf, 'color', [0.75, 0.75, 0.75])

    drawnow;
    pause(0.01);
end