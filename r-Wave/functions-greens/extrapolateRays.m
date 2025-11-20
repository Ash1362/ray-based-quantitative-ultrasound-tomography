function [ray_position_x, ray_position_y, ray_position_z, ray_position_x_left,...
    ray_position_y_left, ray_position_z_left, ray_position_x_right,...
    ray_position_y_right, ray_position_z_right, nan_binary, receiver_indices] =...
    extrapolateRays(ray_position_x, ray_position_y, ray_position_z,...
    ray_position_x_left, ray_position_y_left, ray_position_z_left,...
    ray_position_x_right, ray_position_y_right, ray_position_z_right, ray_spacing)
%EXTRAPOLATERAYS extrapolates the rays for later computing the derivative of
%rays with respect to the initial angles
%
%   % DESCRIPTION:
%       extrapolateRays add points to the rays beyond the reception points
%       using a linear extrapolation. This is required for computing derivate
%       of rays' points with respect to the initial angles.
%
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
%       ray_spacing           - ray spacing [m]
%
%
%
% OUTPUTS:
%       ray_position_x        - the extrapolated x position of the points on
%                               the main (linked) rays
%       ray_position_y        - the extrapolated y position of the points on
%                               the main (linked) rays
%       ray_position_z        - the extrapolated z position of the points on
%                               the main (linked) rays. (Empty, for 2D case)
%       ray_position_x_left   - the extrapolated x position of the points on
%                               the left auxiliary ray. (Empty, if the method
%                               for computing rays' Jacobian is 'raylinked')
%       ray_position_y_left   - the extrapolated y position of the points on
%                               the left auxiliary ray. (Empty, if the method for
%                               computing rays' Jacobian is 'raylinked')
%       ray_position_z_left   - the extrapolated z position of the points on
%                               the left auxiliary ray. (Empty, for 2D case or
%                               if the method for computing rays' Jacobian is
%                               'raylinked')
%       ray_position_x_right  - the extrapolated x position of the points on
%                               the right auxiliary ray. (Empty, if the method
%                               for computing rays' Jacobian is 'raylinked')
%       ray_position_y_right  - the extrapolated y position of the points on
%                               the right auxiliary ray. (Empty, if the method
%                               for computing rays' Jacobian is 'raylinked')
%       ray_position_z_right  - the extrapolated z position of the points on
%                               the right auxiliary ray. (Empty, for 2D case or
%                               if the method for computing rays' Jacobian is
%                              'raylinked')
%       nan_binary            - the binary mask, which is true for the
%                               ray's point outside a binary mask, and will
%                               be given nan values for rays' parameters.
%       receiver_indices      - a vector containing the index of rays'
%                               points matching the reception points, i.e.,
%                               the index of columns in matrices for the rays'
%                               parameters. e.g. ray_poition_x....
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 30.03.2020
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian

% get the number of dimensions
if isempty(ray_position_z)
    dim = 2;
else
    dim = 3;
end

% get the approach for computing the geomtetrical attenuation
if isempty(ray_position_x_left)
    
    % if the left auxiliary ray is empty, the used method for tracing and
    % using auxiliary rays for computing rays' Jacobian  is 'raylinked'.
    % Using this approach, rays' jacobian will be computed using the
    % neighboring rays by finite differences, instead of tracing auxiliary rays.
    attenuation_geom_method = 'raylinked';
    
    % set the auxiliary method nan (not used)
    auxiliary_method = nan;
else
    
    if isempty(ray_position_x_right)
        
        % If the right auxiliary ray is empty, the used aproach for tracing
        % auxiliary rays and computing the rays' Jacobian  is 'paraxial'.
        auxiliary_method = 'paraxial';
        
    else
        
        % If the right auxiliary ray is not empty (given), the used approach
        % for tracing auxiliary rays and computing rays' Jacobian is
        % based on tracing auxiliary rays, independent of the main (linked)
        % ray.
        auxiliary_method = 'angle_perturbation';
    end
    
    % if the left auxiliary ray is not empty, the used auxiliary method is
    % 'auxiliary'
    attenuation_geom_method = 'auxiliary';
    
end

% give an error if the inputs are given such that the method for computing
% rays' trajectory is not set paraxial
if dim == 3 && ~strcmp(auxiliary_method, 'paraxial')
    error(['For 3D case, the method for tracing rays and computing the'...
          'Jacobian of rays must be set paraxial.'])
end

% get the number of transducers
num_transducer = size(ray_position_x, 1);

% Make the last spacing along the rays before interception by the receiver
% equal to ray spacing, and also get the index of the ray's points which match
% the receiver (reception point). Also, compute a binary mask indicating to all
% rays' points inside the detection ring
[ray_position_x, ray_position_y, ray_position_z, nan_binary, receiver_indices] =...
    correctRayLastPosition(ray_position_x, ray_position_y, ray_position_z,...
    ray_spacing, false);


% do the same for the auxiliary rays
if ~isempty(ray_position_x_left)
    [ray_position_x_left, ray_position_y_left, ray_position_z_left, ~,...
        receiver_indices_left] = correctRayLastPosition(...
        ray_position_x_left, ray_position_y_left,...
        ray_position_z_left, ray_spacing, true);
end

if ~isempty(ray_position_x_right)
    [ray_position_x_right, ray_position_y_right, ray_position_z_right,  ~,...
        receiver_indices_right] = correctRayLastPosition(...
        ray_position_x_right, ray_position_y_right,...
        ray_position_z_right, ray_spacing, true);
end


switch attenuation_geom_method
    
    case 'raylinked'
        
        % The arc length of the neighboring rays may not be equal. Therefore,
        % for computing rays' Jacobian on the last points along the rays using
        % neighboring rays, the arc length of the neigboring rays must be
        % increased via extrapolation.
        
        % determine the required number of extrapolation points for each
        % linked ray
        receiver_indices_extended = [receiver_indices(end);...
            receiver_indices; receiver_indices(1)];
        
        % get the required number of extrapolation points
        num_extrapolation_points = movmax(receiver_indices_extended, 3)...
            - receiver_indices_extended;
        
        % increase the required number of extrapolation points by 1,
        % because one more point will be used for computing the derivative
        % of the rays' position with respect to the arc length using a 
        % centered difference approach
        num_extrapolation_points = num_extrapolation_points(2:end-1) + 1;
        
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
            
            % get the required number of extrapolation points for the main (linked)
            % ray, and the left and right auxiliary rays
            % add 
            receiver_indices_max = max(max(receiver_indices_left, receiver_indices),...
                receiver_indices_right) + 1;
            
            % get the required number of extrapolation points
            
            % left auxiliary ray
            num_extrapolation_points_left = receiver_indices_max - receiver_indices_left;
            
            % right auxiliary ray
            num_extrapolation_points_right = receiver_indices_max - receiver_indices_right;
            
            % linked (reference) ray
            num_extrapolation_points = receiver_indices_max - receiver_indices;
            
            % extrapolate the x coordinate of the right auxiliary rays
            ray_position_x_right = addPointsRays(ray_position_x_right,...
                receiver_indices_right, num_extrapolation_points_right, num_transducer);
            
            % extrapolate the y coordinate of the right auxiliary rays
            ray_position_y_right = addPointsRays(ray_position_y_right,...
                receiver_indices_right, num_extrapolation_points_right, num_transducer);
            
            
        else
            
            % If the auxiliary rays are paraxial, set the number of extrapolation points 1
            num_extrapolation_points_left = ones(num_transducer, 1);
            num_extrapolation_points = ones(num_transducer, 1);
            
            % allocate empty variables for the ray position for the right auxiliary ray
            ray_position_x_right = [];
            ray_position_y_right = [];
            
        end
        
        
        % extend the x coordinate of the linked (reference) rays
        ray_position_x = addPointsRays(ray_position_x, receiver_indices,...
            num_extrapolation_points, num_transducer);
        
        % extend the y coordinate of the linked (reference) rays
        ray_position_y = addPointsRays(ray_position_y, receiver_indices,...
            num_extrapolation_points, num_transducer);
        
        % extend the x coordinate of the left auxiliary rays
        ray_position_x_left = addPointsRays(ray_position_x_left,...
            receiver_indices_left, num_extrapolation_points_left, num_transducer);
        
        % extend the y coordinate of the left auxiliary rays
        ray_position_y_left = addPointsRays(ray_position_y_left,...
            receiver_indices_left, num_extrapolation_points_left, num_transducer);
        
        if dim == 3
            
        % extend the z coordinate of the linked (reference) rays
        ray_position_z = addPointsRays(ray_position_z, receiver_indices,...
            num_extrapolation_points, num_transducer);
        
        % extend the z coordinate of the left auxiliary rays
        ray_position_z_left = addPointsRays(ray_position_z_left,...
            receiver_indices_left, num_extrapolation_points_left, num_transducer);
        end  
        
end


end