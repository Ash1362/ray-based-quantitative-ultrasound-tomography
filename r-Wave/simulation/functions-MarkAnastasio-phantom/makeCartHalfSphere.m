function [half_sphere, normal_vectors] = makeCartHalfSphere(radius, num_points, ...
    center_pos, plot_half_sphere, z_bn)
%MAKECARTHALFSPHERE     Create a 3D cartesian half sphere.
%
% DESCRIPTION:
%       makeCartHalfSphere creates a 3 x num_points array of the Cartesian
%       coordinates of points evenly distributed over a half sphere using the
%       Golden Section Spiral method. 
%
% USAGE:
%       sphere = makeCartHalfSphere(radius, num_points)
%       sphere = makeCartHalfSphere(radius, num_points, center_pos)
%       [sphere,normal_vectors] = makeCartHalfSphere(radius, num_points, center_pos, plot_sphere)
%
% INPUTS:
%       radius          - sphere radius [m]
%       num_points      - number of points on the half sphere
%
% OPTIONAL INPUTS
%       center_pos       - [x, y, z] position of the circle center [m] 
%                          (default = [0, 0, 0])
%       plot_half_sphere - Boolean controlling whether the Cartesian points
%                          are plotted (default = false)
%       z_bn             - bounds on z coordinate to exclude too low and high 
%                          points. Denoted in percentages [lb ub], such that the
%                          the first point will be at z=-radius*lb and the last
%                          point will be at z=-radius*(1-ub) (default = [0.5 / num_points, 0])
%
% OUTPUTS:
%       half_sphere    - 3 x num_points array of Cartesian coordinates
%       normal_vectors - 3 x num_points array of normal vectors to the half
%                        sphere surface
%
% ABOUT:
%       modified by     - Felix Lucka
%       date            - 17.02.2017
%       last update     - 24.05.2018
%       
% See also cart2grid, makeCartCircle, makeCartSphere, makeSphere


% check for off_z input
if nargin < 5
    z_bn = [0.5 / num_points, 0];
end

% check for plot_half_sphere input
if nargin < 4
    plot_half_sphere = false;
end


% check for center_pos input
if nargin < 3
    cx = 0;
    cy = 0;
    cz = 0;
else
    cx = center_pos(1);
    cy = center_pos(2);
    cz = center_pos(3);    
end

% generate angle functions using the Golden Section Spiral method
inc = pi * (3-sqrt(5));
k = 0:num_points-1;
z = k * (1 - sum(z_bn))/(num_points-1) + z_bn(1) - 1 ;
r = sqrt(1 - (z.^2));
phi = k * inc;

% create the unit sphere
half_sphere = [cos(phi).*r; sin(phi).*r; z];

% normal vectors 
if nargout > 1
    normal_vectors = half_sphere;
end

% scale up to radius
half_sphere = radius.*half_sphere;

% offset if needed
half_sphere(1, :) = half_sphere(1, :) + cx;
half_sphere(2, :) = half_sphere(2, :) + cy;
half_sphere(3, :) = half_sphere(3, :) + cz;

% plot results1
if plot_half_sphere
    
    % select suitable axis scaling factor
    [x_sc, scale, prefix] = scaleSI(max(half_sphere(:)));  %#ok<ASGLU>
    
    % create the figure
    figure;
    plot3(half_sphere(1,:)*scale, half_sphere(2,:)*scale, half_sphere(3,:)*scale, '.');
    xlabel(['x [' prefix 'm]']);
    ylabel(['y [' prefix 'm]']);
    zlabel(['z [' prefix 'm]']);
    axis equal;
    grid on;
    box on;
    
    if nargout > 1
        hold on;
        quiver3(half_sphere(1,:)*scale, half_sphere(2,:)*scale, half_sphere(3,:)*scale,...
            normal_vectors(1,:), normal_vectors(2,:), normal_vectors(3,:))
    end
end
