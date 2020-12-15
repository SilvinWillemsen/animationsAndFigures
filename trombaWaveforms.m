clear all;
close all;

yTextOffset = 400;
[dataStringBody, fs] = audioread("trombaPureStatesString.m4a");
[dataBridgeBody, fs] = audioread("bridgeLeftBodyRight.wav");
[dataPureStates, fs] = audioread("trombaPureStates.m4a");
% plot(dataStringBody)
% hold on;
figure('Renderer', 'painters', 'Position', [100 100 600 300])
plot(dataPureStates)
range = (150300:151500) + 35;
um = plot(range / fs, dataPureStates(range,1) / 1000, '-.', 'color', [0.20,0.20,0.20], 'Linewidth', 1.2);
hold on;
up = plot(range / fs, dataPureStates(range,2) / 1000, 'color', 'k', 'Linewidth', 1.2);
xlim([range(1), range(end)] / fs);
yticks([(-1:2) .* 1e-5]);
ylim([-1.1e-5, 2.1e-5])
grid on
uoff = plot([range(1),range(end)] / fs, [5e-6, 5e-6], '--', 'color', [0.5,0.5,0.5], 'Linewidth', 1);

% text(3.4325, 6.5e-6, '$u_{off}$', 'interpreter', 'latex', 'Fontsize', 20, 'color', [0.5, 0.5, 0.5] )% hold on;
legend([um, up, uoff],["$w(t)$","$z(x_\textrm{\fontsize{7}{0}\selectfont sm},y_\textrm{\fontsize{7}{0}\selectfont sm},t)$", "$w_\textrm{\fontsize{7}{0}\selectfont off}$"], 'interpreter', 'latex', 'Fontsize', 14, 'Position', [0.7720 0.74 0.2097 0.2304])
text((range(1) - (range(floor(0.065 * length(range))) - range(1))) / fs, 0.5e-5, "Displacement (m)", 'Rotation', 90, 'Fontsize', 20, 'Horizontalalignment', 'center', 'fontname','times')
text(range(floor(length(range)*0.5)) / fs, -1.7e-5, "Time (s)", 'Fontsize', 20, 'Horizontalalignment', 'center', 'fontname','times')
title("Bridge and Body Displacement", 'Position', [3.422 2.1216e-05 1.4211e-14])
set(gca, 'Linewidth', 2, 'Fontsize', 20, 'Position', [0.075, 0.18, 0.92, 0.73], 'fontname','times');  % Set it to times
% scatter(1:length(dataPureStates), dataPureStates(:,2))



% initRange = 156900:161000;
% plot(dataStringBody(initRange, 1))
% hold on;
% % plot(dataStringBody(initRange, 2))
% plot(dataBridgeBody(initRange+147, 1))
% plot(dataBridgeBody(initRange+147, 2))
% range = 1:6000;
% subplot(3,1,1)
% plot(dataStringBody(initRange,1), 'k', 'Linewidth', 2);
% title("String")
% set(gca, 'Fontsize', 16, 'Linewidth', 2, 'Position', [0.12, 0.73, 0.85, 0.22]);
% xlim([0,4100]);
% xticks([])
% grid on
% text(-yTextOffset, 0, "$u_s$", 'Rotation', 90, 'interpreter', 'latex', 'Fontsize', 26, 'Horizontalalignment', 'center')
% 
% subplot(3,1,2)
% plot(dataBridgeBody(initRange+147,1), 'k', 'Linewidth', 2);
% title("Bridge")
% set(gca, 'Fontsize', 16, 'Linewidth', 2, 'Position', [0.12, 0.43, 0.85, 0.22]);
% xlim([0,4100]);
% xticks([])
% grid on
% text(-yTextOffset, 0, "$u_m$", 'Rotation', 90, 'interpreter', 'latex', 'Fontsize', 26, 'Horizontalalignment', 'center')
% 
% 
% subplot(3,1,3)
% newData = dataBridgeBody(initRange+147,2);
% plot(initRange / fs, dataBridgeBody(initRange+147,2), 'k', 'Linewidth', 2)
% title("Body")
% set(gca, 'Fontsize', 16, 'Linewidth', 2, 'Position', [0.12, 0.13, 0.85, 0.22]);
% xlim([initRange(1) / fs, initRange(end) / fs]);
% grid on;
% 
% 
% text(3.6, -0.095, "Time (s)", 'Fontsize', 16)
% text((initRange(1)-yTextOffset) / fs, -0.04, "$u_p$", 'Rotation', 90, 'interpreter', 'latex', 'Fontsize', 26)
