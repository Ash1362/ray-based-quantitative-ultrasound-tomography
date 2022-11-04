function [rays_position_x, rays_position_y, nan_binary, receiver_indices] =...
    correctRayLastPositions(rays_position_x, rays_position_y, ray_spacing,...
    auxiliary_ray)
%CORRECTRAYLASTPOSITIONS makes the last ray spacing just before the reception point
%the smae as the ray spacing along the ray
%
%   % DESCRIPTION:
%       correctRayLastPositions moves the last point along the ray, ie. the
%       reception point, so that the spacing between the reception point
%       and the ray's point just before that equals the ray spacing along
%       the ray

%
% USAGE:
%
%
% INPUTS:
%       ray_position_x        - the x position of the points on the rays
%       ray_position_y        - the y position of the points on the rays
%       ray_spacing           - the ray spacing
%       auxiliary_method      - the method used for including the auxiliary
%                               rays in computing the geometrical
%                               attenuation. This can be set 'paraxial', or
%                               'angle_perturbation'. (This is set nan, if
%                               the auxiliary rays are empty variable.)
% OUTPUTS:
%       ray_position_x        - the corrected x position of the points
%       ray_position_y        - the corrected y position of the points
%       nan_binary            - the binary mask, which is true for the
%                               ray's point outside the detection
%                               surface(ring), and will be given nan values.
%       receiver_indices      - a vector containing the index of columns
%                               correponding to the reception points in
%                               matrices containing the rays' parameters,
%                               eg. ray_poition_x....
%      auxiliary_ray          - Boolean indicating whrther the ray is
%                               auxiliary or not

%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 30.03.2020
%
% This script is part of the r-Wave Tool-box (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian

% get the number of transducers
num_transducer = size(rays_position_x, 1);

% get the maximum number of points along the rays
num_point = size(rays_position_x, 2);

% find the nans in the ray positions
nan_binary = isnan(rays_position_x);

% get an ascending order for columns
[~, column_ascending_orders] = ndgrid(1:num_transducer, 1:num_point);

% get the index of columns for all ray points by removing nans
column_ascending_orders = ~ nan_binary .* column_ascending_orders;

% get the index of columns for the rays' points matching the receivers
[~, receiver_indices] = max(column_ascending_orders, [], 2);

% correct the last position along the ray, if the method for using
% auxiliary rays is not 'paraxial', or i
if ~auxiliary_ray % || strcmp(auxiliary_method, 'paraxial'))
    
    % get the indices of transducers
    rows = (1:num_transducer)';
    
    % get the binary vector of receivers with a minimal distance from the emitter
    % greater than the ray spacing
    binaries = receiver_indices > 2;
    
    % get the linear indices for the point matching the receivers on the
    % position matrix for all rays' points, ie. rays_poition_x or
    % rays_position_y
    receiver_linear_indices = rows(binaries) + num_transducer * (receiver_indices(binaries) - 1);
    
    % make the distance between the last point before the reception point
    % and the reception point equal to the ray spacing. This is needed for accurately
    % computing the geomterical attenuation along the rays
    rays_position_x(receiver_linear_indices) =  2 * rays_position_x(receiver_linear_indices - num_transducer)...
        - rays_position_x(receiver_linear_indices - 2* num_transducer);
    rays_position_y(receiver_linear_indices) =  2* rays_position_y(receiver_linear_indices - num_transducer)...
        - rays_position_y(receiver_linear_indices - 2* num_transducer);
    
    % find the indices for the rays smaller than ray spacing
    % This occurs if the distance of emitter-receiver pair is
    % very close and is less than ray spacing
    indices_short_rays = find(~binaries);
    if ~isempty(indices_short_rays)
        
        for ind_transducer = indices_short_rays
            
            % get the distance of the ray
            cart_lastpoints_x = rays_position_x(ind_transducer, 2)....
                - rays_position_x(ind_transducer, 1);
            cart_lastpoints_y = rays_position_y(ind_transducer, 2)....
                - rays_position_y(ind_transducer, 1);
            distance_lastpoints = norm([cart_lastpoints_x, cart_lastpoints_y]);
            
            % linearly extrapolate the second ray to the ray spacing along each
            % Cartesian coordinate
            rays_position_x(ind_transducer, 2) = rays_position_x(ind_transducer, 1)...
                + ray_spacing/distance_lastpoints * cart_lastpoints_x;
            rays_position_y(ind_transducer, 2) = rays_position_y(ind_transducer, 1)...
                + ray_spacing/distance_lastpoints * cart_lastpoints_y;
        end
        
    end
    
end




end