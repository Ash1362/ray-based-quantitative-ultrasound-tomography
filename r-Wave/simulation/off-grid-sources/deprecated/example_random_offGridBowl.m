% Random Off-grid Bowl Source Example
%
% This example generates a randomly oriented off-grid bowl source. It plots
% the source points and beam axis, along with markers indicating the input
% parameters (bowl_pos, focus_pos, radius, etc.
%
% ABOUT:
%     author      - Elliott Wise
%     date        - 28th March 2017
%     last update - 29th March 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2013 Elliott Wise

% create a kgrid
N = 64;
dx = 1/N;
dy = 1.1/N;
dz = 1.2/N;
kgrid = makeGrid(N, dx, N, dy, N, dz);

% randomly place the bowl source within this grid
bowl_pos = rand(1, 3) - 0.5;
radius = rand(1);
diameter = (rand(1) + 1) * radius;
focus_pos = rand(1, 3)-0.5;

% compute and plot the bowl source points and beam axis
plot_bowl = true;
bowl = offGridBowl(kgrid, bowl_pos, radius, diameter, focus_pos, plot_bowl);

% add markers indicating the bowl position parameters
hold on
plot3(bowl_pos(2), bowl_pos(1), bowl_pos(3), 'go', 'MarkerFaceColor', 'green')
plot3(focus_pos(2), focus_pos(1), focus_pos(3), 'go', 'MarkerFaceColor', 'green')

% print the bowl size parameters
fprintf('Bowl radius:   %.2f\n', radius)
fprintf('Bowl diameter: %.2f\n', diameter)

% overlay two isosurfaces of the mask to make sure it evaluated in the right place
% partial mask
p = patch(isosurface(kgrid.y_vec,kgrid.x_vec,kgrid.z_vec,abs(bowl),0.1));
set(p, 'FaceAlpha', 0.5, 'FaceColor', 'red', 'EdgeColor', 'none');
% full mask
p = patch(isosurface(kgrid.y_vec,kgrid.x_vec,kgrid.z_vec,abs(bowl),0));
set(p, 'FaceAlpha', 0.2, 'FaceColor', 'blue', 'EdgeColor', 'none');
camlight, lighting flat
