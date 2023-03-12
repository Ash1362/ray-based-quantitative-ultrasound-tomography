function [ray_position_x, ray_position_y, ray_position_x_left, ray_position_y_left, ...
    ray_position_x_right, ray_position_y_right, nan_binary, receiver_indices] = extrapolateRays(...
    ray_position_x, ray_position_y, ray_position_x_left, ray_position_y_left,...
    ray_position_x_right, ray_position_y_right, ray_spacing, distance_integer_rays)
%EXTRAPOLATERAYS extrapolates the rays for computing the derivative of
%rays' points with respect to the initial angles
%
%   % DESCRIPTION:
%       extrapolateRays add points to the rays beyond the reception points
%       using a linear extrapolation from the previous rays' points. This
%       is required for an accurate computation of the derivate of the
%       rays' points matching the reception points with respect to the
%       initial angles.
%
% USAGE:
%
%
% INPUTS:
%       ray_position_x        - the x position of the points on the main
%                               (linked) rays
%       ray_position_y        - the y position of the points on the main
%                               (linked) rays
%       ray_position_x_left   - the x position of the points on the left
%                               auxiliary ray. Empty variable, if the
%                               method for computing the ray's Jacobian is
%                               set 'raylinked'
%       ray_position_y_left   - the y position of the points on the left
%                               auxiliary ray. Empty variable, if the
%                               method for computing the ray's Jacobian is
%                               set 'raylinked'
%       ray_position_x_right  - the x position of the points on the right
%                               auxiliary ray. Empty variable, if the
%                               method for computing the ray's Jacobian is
%                               set 'raylinked'
%       ray_position_y_right  - the y position of the points on the right
%                               auxiliary ray. Empty variable, if the
%                               method for computing the ray's Jacobian is
%                               set 'raylinked'
%       ray_spacing           - ray spacing [m]
%       distance_integer_rays - the number of rays between the main linked ray
%                               and the two other linked rays used for
%                               computing the Jacobian of the main ray,
%                               instead of the auxiliary rays. Not applied, if the
%                               method for computing the ray's Jacobian is
%                               set 'auxiliary'
%                               

%
% OUTPUTS:
%       ray_position_x        - the corrected x position of the points on the main
%                               (linked) rays
%       ray_position_y        - the corrected y position of the points on the main
%                               (linked) rays
%       ray_position_x_left   - the corrected x position of the points on the left
%                               auxiliary ray. Empty variable, if the
%                               method for computing the ray's Jacobian is
%                               set 'raylinked'
%       ray_position_y_left   - the corrected y position of the points on the left
%                               auxiliary ray. Empty variable, if the
%                               method for computing the ray's Jacobian is
%                               set 'raylinked'
%       ray_position_x_right  - the corrected x position of the points on the right
%                               auxiliary ray. Empty variable, if the
%                               method for computing the ray's Jacobian is
%                               set 'raylinked'
%       ray_position_y_right  - the corrected y position of the points on the right
%                               auxiliary ray. Empty variable, if the
%                               method for computing the ray's Jacobian is
%                               set 'raylinked'
%       nan_binary            - the binary mask, which is true for the
%                               ray's point outside the detection
%                               surface(ring),and will be given nan values.
%       receiver_indices      - a vector containing the index of columns
%                               correponding to the reception points in
%                               matrices containing the rays' parameters,
%                               eg. ray_poition_x....

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 30.03.2020
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian


% get the approach taken for computing the geomterical attenuation
if isempty(ray_position_x_left)
    
    % if the left auxiliary ray is empty, the used auxiliary method is
    % 'raylinked'. Using this approach, the geomterical attenuation will be
    % computed using the adjacent rays, instead of the auxiliary rays.
    attenuation_geom_method = 'raylinked';
    
    % set the auxiliary method nan (not used)
    auxiliary_method = nan;
else
    
    if isempty(ray_position_x_right)
        
        % If the right auxiliary ray is empty, the used auxiliary method
        % is 'paraxial'
        auxiliary_method = 'paraxial';
    else
        
        auxiliary_method = 'angle_perturbation';
    end
    
    % if the left auxiliary ray is nonempty, the used auxiliary method is
    % 'auxiliary'
    attenuation_geom_method = 'auxiliary';
    
end


% get the number of transducers
num_transducer = size(ray_position_x, 1);

% Make the last spacing along the rays before interception by the receiver
% equal to ray spacing, and also get the index of point along each ray matching
% the receiver. Also, compute a binary mask indicating to all the rays' points
% inside the detection ring
[ray_position_x, ray_position_y, nan_binary, receiver_indices] = correctRayLastPositions(...
    ray_position_x, ray_position_y, ray_spacing, false);


% do the same for the auxiliary rays
if ~isempty(ray_position_x_left)
    [ray_position_x_left, ray_position_y_left, ~, receiver_indices_left]...
        = correctRayLastPositions(ray_position_x_left, ray_position_y_left, ray_spacing,...
        true);
end

if ~isempty(ray_position_x_right)
    [ray_position_x_right, ray_position_y_right, ~, receiver_indices_right]...
        = correctRayLastPositions(ray_position_x_right, ray_position_y_right, ray_spacing,...
        true);
end


switch attenuation_geom_method
    
    case 'raylinked'
        
        % The arc length of the neighboring rays are not necessarily equal. Therefore,
        % for computing the geometrical attenuation on the last points along
        % the rays with larger arc lengths, the arc length of the neigboring
        % rays must be artificially enlarged using extrapolation.
        
        % compute the number of the required extrapolation points for each
        % linked ray
        receiver_indices_extended = [receiver_indices(end-distance_integer_rays:end);...
            receiver_indices; receiver_indices(1:distance_integer_rays+1)];
        num_extrapolation_points = movmax(receiver_indices_extended, 2*distance_integer_rays+3)...
            - receiver_indices_extended;
        num_extrapolation_points = num_extrapolation_points(distance_integer_rays+2:end-distance_integer_rays-1) + 1;
        
        
        % extend the x coordinate of the rays
        [ray_position_x] = addPointsRays(ray_position_x, receiver_indices,...
            num_extrapolation_points, num_transducer);
        
        % extend the y coordinate of the rays
        [ray_position_y] = addPointsRays(ray_position_y, receiver_indices,...
            num_extrapolation_points, num_transducer);
        
        % allocate empty variable for the outputs that are not used
        ray_position_x_left = [];
        ray_position_y_left = [];
        
        ray_position_x_right = [];
        ray_position_y_right = [];
        
        
    case 'auxiliary'
        
        
        if strcmp(auxiliary_method, 'angle_perturbation')
            
            % get the number of required extrapolation points for the main (linked)
            % ray, and the left and right auxiliary rays
            receiver_indices_max = max(max(receiver_indices_left, receiver_indices),...
                receiver_indices_right)+1;
            
            num_extrapolation_points_left = receiver_indices_max - receiver_indices_left;
            num_extrapolation_points_right = receiver_indices_max - receiver_indices_right;
            num_extrapolation_points = receiver_indices_max - receiver_indices;
            
            % extend the x coordinate of the right auxiliary rays
            ray_position_x_right = addPointsRays(ray_position_x_right,...
                receiver_indices_right, num_extrapolation_points_right, num_transducer);
            
            
            % extend the y coordinate of the right auxiliary rays
            ray_position_y_right = addPointsRays(ray_position_y_right,...
                receiver_indices_right, num_extrapolation_points_right, num_transducer);
            
            
        else
            
            % set the number of extrapolation points 1
            num_extrapolation_points_left = ones(num_transducer, 1);
            num_extrapolation_points = ones(num_transducer, 1);
            
            % allocate empty variables for the ray position for the right auxiliary ray
            ray_position_x_right = [];
            ray_position_y_right = [];
            
        end
        
        
        % extend the x coordinate of the main (linked) rays
        ray_position_x = addPointsRays(ray_position_x, receiver_indices,...
            num_extrapolation_points, num_transducer);
        
        % extend the y coordinate of the main (linked) rays
        ray_position_y = addPointsRays(ray_position_y, receiver_indices,...
            num_extrapolation_points, num_transducer);
        
        % extend the x coordinate of the left auxiliary rays
        ray_position_x_left = addPointsRays(ray_position_x_left,...
            receiver_indices_left, num_extrapolation_points_left, num_transducer);
        
        % extend the y coordinate of the left auxiliary rays
        ray_position_y_left = addPointsRays(ray_position_y_left,...
            receiver_indices_left, num_extrapolation_points_left, num_transducer);
        
end


end