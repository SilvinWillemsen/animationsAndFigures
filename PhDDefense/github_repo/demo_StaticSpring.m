%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Project:     2D Springs Plot
%
% Filename:    demo_StaticSpring.m
%
% Author:      Xue-She Wang
%
% E-mail:      xueshe.wang@duke.edu 
%              xueshe.wang@gmail.com
%
% Description: This demo file plots a simple static spring connecting two
%              arbitrary points.
%
% Revision:    2018/11/20
%              - first release.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% Set radius of spring
R = 2;

% Set spring's coil number
num_coil = 30;

% Create an intance of spring
spr = Spring(R, num_coil);

% Create two points
point1 = [3, -8];
point2 = [-2, 6];

% Get coordinates of spring
[x, y] = spr.getSpr(point1, point2);

% Plot spring
plot(x, y, 'LineWidth', 1)
axis([-10, 10, -10, 10])
axis equal

