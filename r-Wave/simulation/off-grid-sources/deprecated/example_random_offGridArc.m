% Random Off-grid Arc Source Example
%
% This example generates a randomly oriented off-grid arc source. It plots
% the source points and beam axis, along with a marker indicating the
% arc's centre.
%
% ABOUT:
%     author      - Elliott Wise
%     date        - 15th May 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2013 Elliott Wise

% create a kgrid
Nx = 32;
dx = 1/Nx;
kgrid = makeGrid(Nx, dx, Nx, dx);

% randomly place the disc source within this grid
arc_pos = rand(1, 2) - 0.5;
radius = 0.3*rand(1);
diameter = (rand(1) + 1) * radius;
focus_pos = rand(1, 2)-0.5;

% compute and plot the disc source points
plot_arc = true;
arc = offGridArc(kgrid, arc_pos, radius, diameter, focus_pos, 'Plot', plot_arc);

% add indicators for the disc centre, and mask extent
hold on
% mask
v = linspace(0, max(abs(arc(:))), 10);
contour(kgrid.y_vec, kgrid.x_vec, abs(arc), v)
% parameters
plot(arc_pos(2), arc_pos(1), 'ro', 'MarkerFaceColor', 'red')
plot(focus_pos(2), focus_pos(1), 'go', 'MarkerFaceColor', 'green')
axis equal
xlim(kgrid.y_size/2*[-1 1])
ylim(kgrid.x_size/2*[-1 1])
