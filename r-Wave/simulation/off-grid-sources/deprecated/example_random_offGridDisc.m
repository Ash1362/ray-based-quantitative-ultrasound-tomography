% Random Off-grid Disc Source Example
%
% This example generates a randomly oriented off-grid disc source. It plots
% the source points and beam axis, along with a marker indicating the
% disc's centre.
%
% ABOUT:
%     author - Elliott Wise
%     date   - 30th March 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2013 Elliott Wise

% number of simulation dimensions
dim = 3;

% create a kgrid
Nx = 64;
dx = 1/Nx;
switch dim
    case 2, kgrid = makeGrid(Nx, dx, Nx, dx);
    case 3, kgrid = makeGrid(Nx, dx, Nx, dx, Nx, dx);
end

% randomly place the disc source within this grid
centre = rand(1, dim) - 0.5;
radius = 0.3*rand(1);
if dim > 2
    orientation = rand(1, 3);
else
    orientation = [];
end

% compute and plot the disc source points
plot_disc = true;
disc = offGridDisc(kgrid, centre, radius, orientation, plot_disc);

% add indicators for the disc centre, and mask extent
hold on
switch dim
    case 2
        % mask
        v = linspace(0, max(abs(disc(:))), 10);
        contour(kgrid.y_vec, kgrid.x_vec, abs(disc), v)
        % parameters
        plot(centre(2), centre(1), 'ro', 'MarkerFaceColor', 'red')
    case 3
        % parameters
        plot3(centre(2), centre(1), centre(3), 'go', 'MarkerFaceColor', 'green')
        % partial mask
        p = patch(isosurface(kgrid.y_vec,kgrid.x_vec,kgrid.z_vec,abs(disc),0.1));
        set(p, 'FaceAlpha', 0.5, 'FaceColor', 'red', 'EdgeColor', 'none');
        % full mask
        p = patch(isosurface(kgrid.y_vec,kgrid.x_vec,kgrid.z_vec,abs(disc),0));
        set(p, 'FaceAlpha', 0.2, 'FaceColor', 'blue', 'EdgeColor', 'none');
        camlight, lighting flat
end
