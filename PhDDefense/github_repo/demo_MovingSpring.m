%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Project:     2D Springs Plot
%
% Filename:    demo_MovingSpring.m
%
% Author:      Xue-She Wang
%
% E-mail:      xueshe.wang@duke.edu 
%              xueshe.wang@gmail.com
%
% Description: This demo file makes an animation of three moving points
%              connectted by two different springs. 
%
% Revision:    2018/11/20
%              - first release.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% Create instances of two springs
spr1 = Spring(0.8, 25); % 1st input: radius, 2nd input: number of coils
spr2 = Spring(0.5, 10);

% Initialize figure
pt_fig = plot(0, '.', 'MarkerSize', 20);
hold on
axis equal
spr1_fig = plot(0, 'k', 'LineWidth', 1);
spr2_fig = plot(0, 'k', 'LineWidth', 0.5);

% Set time of animation
t = linspace(0, 10, 300);

for i = 1 : length(t)
    
    % Update points
    pt1 = [cos(t(i)), 2 * sin(t(i))];
    pt2 = [10 * cos(2 * t(i)), (7 + 3 * sin(6 * t(i))) * sin(2 * t(i))];
    pt3 = [2 * sin(3 * t(i)), -1];
    
    % Update springs
    [spr1_x, spr1_y] = spr1.getSpr(pt1, pt2);
    [spr2_x, spr2_y] = spr2.getSpr(pt1, pt3);
    
    % Update points plotting
    set(pt_fig, 'XData', [pt1(1), pt2(1), pt3(1)], ...
                'YData', [pt1(2), pt2(2), pt3(2)]);
    
    % Update springs plotting
    set(spr1_fig, 'XData', spr1_x, 'YData', spr1_y);
    set(spr2_fig, 'XData', spr2_x, 'YData', spr2_y);
            
    axis([-10, 10, -10, 10])
    pause(0.03)
    
end