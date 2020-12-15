close all;
clear all;

load bowedStringData.mat

x_b = 0.2;
x_mp = 0.8;
y_mp = 0.75;

x = 0:20;
y = 0:10;
z = 1;

x_b_coord = floor(length(u1) * x_b);
x_mp_coord = floor(length(x) * x_mp);
y_mp_coord = floor(length(y) * y_mp);

figure('Renderer', 'painters', 'Position', [100 100 600 300])
hold on;
% grid on;
view(25,14);

grid on

% Bowed string
box on
u_s = plot3([0:length(u1)-1] / (length(u1) - 1) * (length(x)-1), ones(length(u1)) * y_mp_coord, (u1-1e-5) * 8000 + 2, 'b', 'Linewidth', 2);

%Bridge
u_m = scatter3(x_mp_coord, y_mp_coord, 2, 100, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[145/255, 255/255 , 145/255]);
% text(x_mp_coord + 0.4, y_mp_coord, 2.1, '$u_{m}$', 'interpreter', 'latex', 'Fontsize', 18);

% Bowing position
scatter3((x_b_coord * (length(x) - 1.8)) / (length(u1) - 1), y_mp_coord, 2 ... %(u1(x_b_coord)-1e-5) * 8000 + 2
    , 20, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k');
plot3([x_b_coord, x_b_coord] * (length(x) - 1.8) / (length(u1) - 1), [y_mp_coord, y_mp_coord], [2, (u1(x_b_coord)-1e-5) * 8000 + 2], '--', 'color', 'k', 'Linewidth', 1.5);
text(x_b_coord * (length(x) - 1) / (length(u1) - 1) + 0.2, y_mp_coord, 2.1, '$\chi_\textrm{\fontsize{7}{0}\selectfont b}$', 'interpreter', 'latex', 'Fontsize', 18);

%Finger coordinate on string
scatter3((length(u1) / 2 * length(x)) / (length(u1) - 1), y_mp_coord, 2 ... %(u1(x_b_coord)-1e-5) * 8000 + 2
    , 20, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k');
text((length(u1) / 2 * length(x)) / (length(u1) - 1) + 0.2, y_mp_coord, 2.1, '$\chi_\textrm{\fontsize{7}{0}\selectfont f}$', 'interpreter', 'latex', 'Fontsize', 18);

%\eta_{mp}
plot3([x_mp_coord, x_mp_coord], [y_mp_coord, y_mp_coord], [2, 1], '-.', 'color', 'k', 'Linewidth', 1.5);
text(x_mp_coord + 0.2, y_mp_coord, 1.5, '$\eta_\textrm{\fontsize{7}{0}\selectfont mp}$', 'interpreter', 'latex', 'Fontsize', 18);

%Bridge coordinate on plate
scatter3(x_mp_coord, y_mp_coord, 1, 20, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k');
text(x_mp_coord + 0.4, y_mp_coord, 1.15, '$(x_\textrm{\fontsize{7}{0}\selectfont mp}, y_\textrm{\fontsize{7}{0}\selectfont mp})$', 'interpreter', 'latex', 'Fontsize', 18)

% bridge coordinate on string
scatter3(x_mp_coord, y_mp_coord-0.001, 2, 20, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k');
text(x_mp_coord + 0.4, y_mp_coord, 2.1, '$\chi_\textrm{\fontsize{7}{0}\selectfont sm}$', 'interpreter', 'latex', 'Fontsize', 18)

% plate
Xnames = {'$0$'; '$N_x$'};
Ynames = {'$0$'; '$N_y$'};
Znames = {'$0$'; '$w_\textrm{\fontsize{7}{0}\selectfont off}$'};
for zz = 1:length(z)
    for l = 1:length(y)
        plot3([x(1) x(end)],[y(l) y(l)], [z(zz) z(zz)], 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
    end
    for m = 1:length(x)
        u_p = plot3([x(m) x(m)],[y(1) y(end)], [z(zz) z(zz)], 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
    end
end

xlim([0 length(x)-1])
ylim([0 length(y)-1])
zlim([1 2.5])

set(gca, 'Position', [0.07 0.08 0.86 0.94], 'Projection', 'perspective')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'xtick',[0,length(x)-1],'xticklabel', Xnames, 'FontSize', 18)

% text([1,length(x)-1], repmat(-0.5, length(Xnames), 1), repmat(0.8, length(Xnames), 1), Xnames, 'horizontalAlignment', 'center', 'Fontsize', 18, 'interpreter', 'latex')
set(gca,'ytick',[0, length(y)-1],'yticklabel', Ynames, 'FontSize', 18) 
set(gca,'ztick',[0, 2],'zticklabel', Znames, 'FontSize', 18) 

legend([u_s(1), u_m, u_p], ["$u$", "$w$", "$z$"], 'interpreter', 'latex', 'Location', 'northeast');