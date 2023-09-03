function [attenuation_rays, nan_binary, receiver_indices, caustic_rays,...
    ray_direction_x, ray_direction_y, ray_direction_z] = calcGeomAttenuationRays(...
    ray_position_x, ray_position_y, ray_position_z, ray_position_x_left,...
    ray_position_y_left, ray_position_z_left, ray_position_x_right,...
    ray_position_y_right, ray_position_z_right, ray_time, ray_spacing,...
    do_get_direction)
%CALCGEOMATTENUATIONRAYS computes the geometrical attenuation including
%refraction effects along the rays
%
%   % DESCRIPTION:
%       calcGeomAttenuationRays computes the geometrical attenuation
%       including the refraction effects on the rays' sampled points
%       with respect to a reference point, which is chosen the second point
%       along each ray.
% USAGE:
%
%
% INPUTS:
%       ray_position_x        - the x position of the points on the main
%                               (linked) rays
%       ray_position_y        - the y position of the points on the main
%                               (linked) rays
%       ray_position_z        - the z position of the points on the main
%                               (linked) rays. (Empty, for 2D case)
%       ray_position_x_left   - the x position of the points on the left
%                               auxiliary ray. (Empty, if the method for
%                               computing rays' Jacobian is 'raylinked')
%       ray_position_y_left   - the y position of the points on the left
%                               auxiliary ray. (Empty, if the method for
%                               computing rays' Jacobian is 'raylinked')
%       ray_position_z_left   - the z position of the points on the left
%                               auxiliary ray. (Empty, for 2D case or if the
%                               method for computing rays' Jacobian is 'raylinked')
%       ray_position_x_right  - the x position of the points on the right
%                               auxiliary ray. (Empty, if the method for
%                               computing rays' Jacobian is 'raylinked')
%       ray_position_y_right  - the y position of the points on the right
%                               auxiliary ray. (Empty, if the method for
%                               computing rays' Jacobian is 'raylinked')
%       ray_position_z_right  - the z position of the points on the right
%                               auxiliary ray. (Empty, for 2D case or if the
%                               method for computing rays' Jacobian is 'raylinked')
%       ray_time              - the accumulated time delays [s] along the
%                               main (linked) rays, ie. acoustic length
%                               divided by the reference sound speed
%       ray_spacing           - ray spacing [m]
%       distance_integer_rays - the number of rays between the main linked ray
%                               and the two other linked rays used for
%                               computing the Jacobian of the main ray,
%                               instead of the auxiliary rays. Not applied, if the
%                               method for computing the ray's Jacobian is
%                               set 'auxiliary'
%       greens_method         - the mathod for approximating the Green's
%                               function
%       do_get_direction      - Boolean controllong whether the cartesian
%                               direction of the ray is computed or not.
%
% OUTPUTS:
%       attenuation_rays     - the geomterical attenuation on the main rays' points
%                              with respect to the chosen reference point,
%                              ie. the second point on each ray
%       nan_binary           - the binary mask, which is true for rays'
%                              points outside the detection surface(ring).
%                              the rays' parameters for those points will
%                              be set nan (not computed).
%       receiver_indices     - a vector containing the index of columns
%                              correponding to the reception points in
%                              matrices containing the rays' parameters,
%                              eg. ray_poition_x....
%       caustic_rays         - the cumulative times the sign of the rays' Jacobian
%                              along the rays has been changed
%       ray_direction_x      - the x Cartesian direction of the rays
%       ray_direction_y      - the y Cartesian direction of the rays
%       ray_direction_z      - the y Cartesian direction of the rays
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 30.03.2020
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian



% extrapolate the rays beyond the receivers linearly for ensuring that the
% arc length on the auxiliary rays or neighboring rays is sufficiently large
% for computing the Jacobian on last points on the main rays
[ray_position_x, ray_position_y, ray_position_z, ray_position_x_left, ray_position_y_left,...
    ray_position_z_left, ray_position_x_right, ray_position_y_right, ray_position_z_right,...
    nan_binary, receiver_indices] =...
    extrapolateRays(ray_position_x, ray_position_y, ray_position_z,...
    ray_position_x_left, ray_position_y_left, ray_position_z_left,...
    ray_position_x_right, ray_position_y_right, ray_position_z_right,...
    ray_spacing);

% compute the rays' Jacobian
[jacobian_rays, caustic_rays, ray_direction_x, ray_direction_y, ray_direction_z] =...
    calcGeomAttenuationJacobian(ray_position_x, ray_position_y, ray_position_z,...
    ray_position_x_left, ray_position_y_left, ray_position_z_left,...
    ray_position_x_right, ray_position_y_right, ray_position_z_right,...
    nan_binary, do_get_direction);

if isempty(ray_position_z)
    
% compute the geometrical portion of the acoustic attenuation
attenuation_rays = sqrt( (jacobian_rays(:, 2)./jacobian_rays) .*...
    (ray_time./ray_time(:, 2)) );

else
    
% compute the geometrical portion of the acoustic attenuation
attenuation_rays = sqrt(jacobian_rays(:, 2)./jacobian_rays) .*...
    (1./ray_spacing);    
    
end

end