function [ray_position_x, ray_position_y, ray_position_z, nan_binary, receiver_indices] =...
    correctRayLastPosition(ray_position_x, ray_position_y, ray_position_z, ray_spacing,...
    auxiliary_ray)
%CORRECTRAYLASTPOSITION makes the last ray spacing just before the reception point
%the same as the ray spacing along the ray
%
%   % DESCRIPTION:
%       correctRayLastPosition moves the last point on the ray forward, i.e.,
%       the reception point forward so that the spacing between the reception
%       point and the ray's point just before that be equal to the ray spacing
%       along the ray

%
% USAGE:
%
%
% INPUTS:
%       ray_position_x        - the x cartesian position of the points on the rays
%       ray_position_y        - the y cartesian position of the points on the rays
%       ray_position-z        - the z cartesian position of the points on the rays
%       ray_spacing           - the ray spacing [m]
%       auxiliary_ray         - Boolean determining the ray is auxiliary or not
% OUTPUTS:
%       ray_position_x        - the corrected x cartesian position of the points on the rays
%       ray_position_y        - the corrected y cartesian position of the points on the rays
%       ray_position-z        - the corrected z cartesian position of the points on the rays
%       nan_binary            - the binary mask, which is true for the
%                               ray's point outside a chosen binary mask.
%       receiver_indices      - a vector containing the index of columns
%                               correponding to the reception points in
%                               matrices containing the rays' parameters,
%                               eg. ray_poition_x
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 30.03.2020
%
% This script is part of the r-Wave toolbox.
% Copyright (c) 2022 Ashkan Javaherian

% get the number of dimensions
if isempty(ray_position_z)
    dim = 2;
else
    dim = 3;
end



% get the number of transducers
num_transducer = size(ray_position_x, 1);

% get the maximum number of points along the rays
num_point = size(ray_position_x, 2);

% find the nans in the matries defining rays' parameters, e.g.,
% rays' x position 
nan_binary = isnan(ray_position_x);

% get an ascending order for columns
[~, column_ascending_orders] = ndgrid(1:num_transducer, 1:num_point);

% get the index of columns for all ray points by removing nans
column_ascending_orders = ~nan_binary .* column_ascending_orders;

% get the index of point matching the receiver for each linked ray, i.e.,
% the index of column associated with the reception point for each row 
[~, receiver_indices] = max(column_ascending_orders, [], 2);


%%=========================================================================
% CORRECT THE LAST RAY SPACING ALONG THE RAY
%==========================================================================
% correct the last ray spacing along the ray, if the ray is not auxiliary
if ~auxiliary_ray 
    
    % get the indices of transducers
    rows = (1:num_transducer)';
    
    % get the binary vector of receivers with a minimal distance from the emitter
    % greater than the ray spacing
    receiver_binaries = receiver_indices > 2;
    
    % get the linear indices for the points matching the receivers on the
    % position matrix
    receiver_linear_indices = rows(receiver_binaries) +...
        num_transducer * (receiver_indices(receiver_binaries) - 1);
    
    % make the distance between the reception point and the last point before 
    % the reception point equal to the ray spacing. This is needed for accurately
    % computing the geometrical attenuation along the rays
    % x coordinate
    ray_position_x(receiver_linear_indices) =...
        2 * ray_position_x(receiver_linear_indices - num_transducer)...
        - ray_position_x(receiver_linear_indices - 2 * num_transducer);
    % y coordinate
    ray_position_y(receiver_linear_indices) =...
        2 * ray_position_y(receiver_linear_indices - num_transducer)...
        - ray_position_y(receiver_linear_indices - 2 * num_transducer);
    % z coordinate
    if dim == 3
        ray_position_z(receiver_linear_indices) =...
            2 * ray_position_z(receiver_linear_indices - num_transducer)...
            - ray_position_z(receiver_linear_indices - 2 * num_transducer);
    end
    
    % find the linear indices for the reception point of rays smaller than
    % ray spacing
    % This occurs if the distance of emitter-receiver pair is
    % very close and is less than ray spacing
    indices_short_rays = find(~receiver_binaries);
    
    if ~isempty(indices_short_rays)
        
        for ind_transducer = indices_short_rays
            
            % make the distance between the two rays' points equal to the
            % ray spacing
            cart_lastpoint_x = ray_position_x(ind_transducer, 2)....
                - ray_position_x(ind_transducer, 1);
            cart_lastpoint_y = ray_position_y(ind_transducer, 2)....
                - ray_position_y(ind_transducer, 1);
            if dim == 3
                cart_lastpoint_z = ray_position_z(ind_transducer, 2)....
                    - ray_position_z(ind_transducer, 1);
            end
            
            % get the current distance [m] of the two rays' points
            switch dim
                case 2
            distance_lastpoints = norm([cart_lastpoint_x, cart_lastpoint_y]);
                case 3
            distance_lastpoints = norm([cart_lastpoint_x, cart_lastpoint_y,...
                cart_lastpoint_z]);
            end
                    
            % linearly extrapolate the second ray to corrcet the ray spacing along each
            % Cartesian coordinate
            ray_position_x(ind_transducer, 2) = ray_position_x(ind_transducer, 1)...
                + ray_spacing/distance_lastpoints * cart_lastpoint_x;
            ray_position_y(ind_transducer, 2) = ray_position_y(ind_transducer, 1)...
                + ray_spacing/distance_lastpoints * cart_lastpoint_y;
            
            if dim == 3
                ray_position_z(ind_transducer, 2) = ray_position_z(ind_transducer, 1)...
                    + ray_spacing/distance_lastpoints * cart_lastpoint_z;
            end
                
        end
        
    end
    
end
    

end