% Random Off-grid Line Source Example
%
% This example generates a random off-grid line source. It plots the source
% points, along with markers indicating the input parameters.
%
% ABOUT:
%     author - Elliott Wise
%     date   - 29th March 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2013 Elliott Wise

% number of simulation dimensions
dim = 3;

% create a kgrid
Nx = 64;
dx = 1/Nx;
switch dim
    case 1, kgrid = makeGrid(Nx, dx);
    case 2, kgrid = makeGrid(Nx, dx, Nx, dx);
    case 3, kgrid = makeGrid(Nx, dx, Nx, dx, Nx, dx);
end

% randomly place the line source within this grid
startpoint = rand(1, dim) - 0.5;
endpoint = rand(1, dim) - 0.5;

% compute and plot the line source points
plot_line = true;
line = offGridLine(kgrid, startpoint, endpoint, plot_line);

% add indicators for the line position parameters, and mask extent
hold on
switch dim
    case 1
        plot(kgrid.x_vec, line, 'b:') % mask
        plot(startpoint, 1, 'ro')     % parameters
        plot(endpoint, 1, 'ro')
        ylim([-0.5 max(abs(line(:)))])
    case 2
        % mask
        v = linspace(0, max(abs(line(:))), 10);
        contour(kgrid.y_vec, kgrid.x_vec, abs(line), v)
        % parameters
        plot(startpoint(2), startpoint(1), 'ro', 'MarkerFaceColor', 'red')
        plot(endpoint(2), endpoint(1), 'ro', 'MarkerFaceColor', 'red')
    case 3
        % parameters
        plot3(startpoint(2), startpoint(1), startpoint(3), 'go', 'MarkerFaceColor', 'green')
        plot3(endpoint(2), endpoint(1), endpoint(3), 'go', 'MarkerFaceColor', 'green')
        % partial mask
        p = patch(isosurface(kgrid.y_vec,kgrid.x_vec,kgrid.z_vec,abs(line),0.1));
        set(p, 'FaceAlpha', 0.5, 'FaceColor', 'red', 'EdgeColor', 'none');
        % full mask
        p = patch(isosurface(kgrid.y_vec,kgrid.x_vec,kgrid.z_vec,abs(line),0));
        set(p, 'FaceAlpha', 0.2, 'FaceColor', 'blue', 'EdgeColor', 'none');
        camlight, lighting flat
end

