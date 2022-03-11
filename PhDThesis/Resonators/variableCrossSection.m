close all;
clear all;
figure('Position', [440 574 695 224])

discrete = false;
x = 0:1/1000:1;
freq = 2;
offset = 1;
amp = 0.25;

tube = amp * cos(2 * freq * pi * x) + offset

plot(x, tube, 'k', 'Linewidth', 2)
hold on;
plot(x, -tube, 'k', 'Linewidth', 2)
xlim([-0.05, 1.05])

stepsize = 1;
idxJump = length(x)-1 * stepsize
l = 0;
lLoc = 1 / (2*stepsize)
for xloc = 0:stepsize:1
    plot([xloc, xloc], [tube(round(xloc*(length(x)-1))+1), -tube(round(xloc*(length(x)-1))+1)], 'k', 'Linewidth', 1)
    if discrete
        switch (l)
            case 0
                txt = "$\bar S_0$";
            case lLoc - 1
                txt = "$\bar S_{l-1}$";
            case lLoc
                txt = "$\bar S_l$";
            case lLoc + 1
                txt = "$\bar S_{l+1}$";
            case 1 / stepsize 
                txt = "$\bar S_N$";
            otherwise
                txt = "";

        end
    else
        switch (l)
            case 0
                txt = "$S(0)$";
            case 1 / stepsize 
                txt = "$S(L)$";
            otherwise
                txt = "";

        end
    end
    l = l + 1;
    text(xloc + 0.005, 0, txt, 'interpreter', 'latex', 'Fontsize', 16)
    
end

if discrete
    l = 0;
    for xloc = stepsize/2:stepsize:1
        plot([xloc, xloc], [tube(round(xloc*(length(x)-1))+1), -tube(round(xloc*(length(x)-1))+1)], '--k', 'Linewidth', 1)

        switch (l)
            case 0
                txt = "$S_{1/2}$";
            case lLoc - 1
                txt = "$S_{l-1/2}$";
            case lLoc
                txt = "$S_{l+1/2}$";
            case 1/stepsize - 1
                txt = "$S_{N-1/2}$";

            otherwise
                txt = "";
        end

        l = l + 1;

        text(xloc + 0.005, 0, txt, 'interpreter', 'latex', 'Fontsize', 16)
    end
    offset = 0.075
    plot([2*stepsize, 3*stepsize], [1.5, 1.5], 'k', 'Linewidth', 1.5)
    plot([2*stepsize, 2*stepsize], [1.5+offset, 1.5-offset], 'k', 'Linewidth', 1.5)
    plot([3*stepsize, 3*stepsize], [1.5+offset, 1.5-offset], 'k', 'Linewidth', 1.5)
    text(2.5*stepsize, 1.7, '$h$', 'interpreter', 'latex', 'horizontalAlignment', 'center', 'Fontsize', 18)
end

axis off
set(gca, 'Position', [0 0 0.97 1])
set(gcf, 'color', 'w')