function [parameters_grid, parameters_receiver, nan_grid_binary, caustic_number,...
    receiver_order, caustic_receiver, directions_grid, interp_receiver] = calcParametersGreens(...
    ray_position, ray_time, ray_absorption, ray_position_left, ray_position_right,...
    sound_speed_relative, grid_x, grid_y, ray_spacing, detec_radius, mask_radius_relative,...
    field_mode, greens_method, do_get_direction, interp_receiver)
%CALCPARAMETERSGREENS computes the parameters of the Green's function on the
% grid points, as well as receievers
%
%   % DESCRIPTION:
%       calcParametersGreens computes the cumulative time delays, the geometrical
%       attenuation including refraction effects, cumulative acoustic absorption,
%       and stores them in matrix for all grid points inside a given binary
%       mask, as well as all reception points. For the grid points, relative
%       sound speed (sound speed over the sound speed in water) is also added
%       as the last column to this matrix.
%
%
% USAGE:
%
%
% INPUTS:
%       ray_position         - the position of the points on the main
%                               (linked) rays in all Cartesian coordinates
%       ray_time             - the accumulated time delays [s] along the
%                               main (linked) rays, ie. acoustic length
%                               divided by the sound speed in water
%                               (reference sound speed)
%       ray_absorption       - the accumulated acoustic absorption long the
%                              main (linked) rays
%       ray_position_left    - the position on the left auxiliary rays in all
%                              Cartesian coordinates
%       ray_position_right   - the position on the right auxiliary rays in all
%                              Cartesian coordinates
%       sound_speed_relative - the sound speed on the grid points over the
%                              sound speed in water
%       grid_x               - the x coordinate of the grid points
%       grid_y               - the y coordinate of the grid points
%       ray_spacing          - the scalar value representing the ray
%                              spacing [m]
%       detec_radius         - the radius [m] of the detection ring
%                              (surface)
%       mask_radius_relative - the radius of the mask over the radius of
%                              the detection ring (surface)
%       field_mode           - the field which can be 'forward' or
%                              'adjoint'
%       greens_method        - the approach for computing the Green's
%                              function. The only supported approach is
%                              'frequency', and other approaches are depreacted.
%       directions_grid       - the ray directions on the computational
%                               grid
%       do_get_direction      - Boolean controllong whether the cartesian
%                               direction of the ray is computed or not.
%       interp_receiver      - a struct containing the interpolation
%                              operator for lineraly interpolating the
%                              last point of the rays on the detection ring (surface)
%                              to the the reception points. This is used if
%                              the ray linking is not done, and the rays
%                              are traced with equidistant initial angles.
%                              Avoiding ray linking is not recommended
%                              strongly, but included here as a benchmark
%                              for showing the effectiveness of the ray
%                              linking on the computed pressure field on
%                              the receivers. This contains the fields:
%       'angles'            - the angular position of the receivers
%       'indices'           - the ascending order of the angular poition of
%                            the receivers

% OUTPUT
%       parameters_grid      - a matrix with rows the grid points inside
%                              the binary mask, and columns: time delays,
%                              geoemetrical attenutaion, acoustic
%                              absorption (if not zero), and relative sound speed
%       parameters_receivers- a matrix with rows the receivers, and columns:
%                                time delays, geoemetrical attenutaion, acoustic
%                                absorption (if not zero), and relative sound speed
%       nan_grid_binary        - a binary mask indicating to the grid
%                                points outside the given binary mask
%       caustic_number        - the cumulative integer number of the caustics
%                                from the rays emanated from the emitter
%       receiver_order        - a vector indicating the order of receivers in
%                               the variable 'parameters_receivers'.
%       caustic_receiver      - the cumulative integer number of caustics before the ray is
%                               intercepted by the reception point
%       field_mode            - whether the parameters are computed for the
%                            forward field or the adjoint field. This can be set
%                            'forward' or 'adjoint'
%       greens_method         - the method for approximating the Green's function
%                               (only 'frequency' is supported, and other approaches
%                                are deprecated.)
%       polar direction       - the polar perturbed direction of the rays from
%                               a perturbation to the initial position
%       interp_receiver    - A struct containing the information required for interpolation operator
%                            for linearly interpolating the last point of the
%                            rays on the detection ring (surface) to the
%                            reception points. An empty variable, if not applied.
%                            This contains the fields:
%       'angles'           - the angular position of the receivers
%       'indices'          - the ascending order of the angular poition of
%                            the receivers
%       'lastpoint_angles' - the angular position of the last point of the
%                            rays
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 30.03.2020
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian


if isscalar(ray_absorption)
    
    % a scalar absorption coefficient as
    % the input ignores the acoustic absorption,
    % and make the accumulated acoustic absorption along
    % the rays zero
    do_absorption = false;
    ray_absorption = 0;
else
    do_absorption = true;
end

% choose the attenuation_geom_method for computing the geometrical attenuation along the rays
if isempty(ray_position_left)
    
    % compute the ray's Jacobian with respect to the initial angle
    % only using the linked (main) rays (not recommended)
    attenuation_geom_method = 'raylinked';
    
    % set the auxiliary method nan (not used)
    auxiliary_method = nan;
    
else
    
    % compute the ray's Jacobian with repect to the initial angle
    % using the auxiliary rays (recommended)
    attenuation_geom_method = 'auxiliary';
    
    if isempty(ray_position_right)
        
        % set the auxiliary method 'paraxial'
        auxiliary_method = 'paraxial';
    else
        
        % set the auxiliary method 'angle perturbation'
        auxiliary_method = 'angle_perturbation';
    end
    
    
end




% get the position of the points on the rays along each Cartesian coordinate
ray_position_x = ray_position(1:2:size(ray_position, 1), :);
ray_position_y = ray_position(2:2:size(ray_position, 1), :);
clear ray_position

switch attenuation_geom_method
    
    
    case 'auxiliary'
        
        % the distance between rays for computing the geometrical attenuation
        distance_integer_rays = [];
        
        % get the position of the points along the left reference rays along
        % each Cartesian coordinate
        ray_position_x_left = ray_position_left(1:2:size(ray_position_left, 1), :);
        ray_position_y_left = ray_position_left(2:2:size(ray_position_left, 1), :);
        clear ray_position_left
        
        
        switch auxiliary_method
            
            case 'angle_perturbation'
                
                % get the position of the points along the right reference rays along
                % each Cartesian coordinate
                ray_position_x_right = ray_position_right(1:2:size(ray_position_right, 1), :);
                ray_position_y_right = ray_position_right(2:2:size(ray_position_right, 1), :);
                clear ray_position_right
                
            case 'paraxial'
                
                % allocate empty variables to the coordinates of the right
                % auxiliary ray, if the auxiliary method is paraxial
                ray_position_x_right = [];
                ray_position_y_right = [];
                
              
        end
        
    case 'raylinked'
        
        % the distance of rays used for approximating the Greens
        switch field_mode
            case'forward'
                distance_integer_rays = 0;
            case 'adjoint'
                distance_integer_rays = 0;
        end
        
        % allocate empty variables to the position of points for the reference
        % rays, if the attenuation_geom_method is 'raylinked'
        % position of the points along the left reference rays for each Cartesian coordinate
        ray_position_x_left = [];
        ray_position_y_left = [];
        
        % position of the points along the right reference rays for each Cartesian coordinate
        ray_position_x_right = [];
        ray_position_y_right = [];
end


switch greens_method
    case 'analytic'
        % calculate the attenuation along the rays
        [ray_attenuation, nan_binary, receiver_indices, caustic_rays,...
            ray_direction_x, ray_direction_y] = calcGeomAttenuationRays(...
            ray_position_x, ray_position_y, ray_position_x_left, ray_position_y_left,...
            ray_position_x_right, ray_position_y_right, ray_time, ray_spacing,...
            distance_integer_rays, greens_method, do_get_direction);
    otherwise
        notImpErr
end


% discard the emission point from the rays' parameters
ray_position_x(:, 1) = [];
ray_position_y(:, 1) = [];
ray_time(:, 1) = [];
ray_attenuation(:, 1) = [];
nan_binary(:,1) = [];
if do_absorption
    ray_absorption(:, 1) = [];
end
if ~isempty(ray_direction_x)
    % remove the initial point from the x direction along the ray
    ray_direction_x(:, 1) = [];
    % remove the initial point from the y direction along the ray
    ray_direction_y(:, 1) = [];
   
end



switch field_mode
    case 'forward'
        
        % calculate the linear indices of receivers, considering that the first
        % column (emission point) has been removed from the rays
        receiver_indices = sub2ind(size(ray_position_x), 1:size(ray_position_x, 1),...
            receiver_indices' - 1);
        [receiver_binaries, receiver_order]= ismember(find(~nan_binary), receiver_indices);
        
        receiver_order = receiver_order(receiver_binaries);
        
    case 'adjoint'
        
        
        if strcmp(attenuation_geom_method, 'raylinked')  &&   any(strcmp(greens_method, {'time', 'frequency2'}))
            
            % because the adjoint field is calculated only for a cone with axis
            % the receiver-emitter pair, the geometrical attenuation for the
            % surface of the cone cannot be comuted, and must be discarded.
            % Note that for computing the adjoint field, distance_integer_rays is
            % always set zero
            ray_position_x([1, end], :) = [];
            ray_position_y([1, end], :) = [];
            ray_time([1, end], :) = [];
            ray_attenuation([1, end], :) = [];
            nan_binary([1, end], :) = [];
            if do_absorption
                ray_absorption([1, end], :) = [];
            end
            
        end
        
        receiver_order = [];
        
        
end



% remove the nans from the rays' parameters
ray_position_x = ray_position_x(~nan_binary);
ray_position_y = ray_position_y(~nan_binary);
ray_time = ray_time(~nan_binary);
ray_attenuation = ray_attenuation(~nan_binary);
if do_absorption
    ray_absorption = ray_absorption(~nan_binary);
end
if ~isempty(ray_direction_x)
    ray_direction_x = ray_direction_x(~nan_binary);
    ray_direction_y = ray_direction_y(~nan_binary);
end


% get a binary mask for confining the rays points to the ROI mask
ray_mask = ray_position_x.^2 + ray_position_y.^2 <= (mask_radius_relative * detec_radius)^2;


% get the number of ray points inside the defined mask
num_raypoints = nnz(ray_mask);

% construct the ray-to-grid interpolation matrix
[gridpoint_indices, tri, bc, nan_grid_binary] = gridDataFast2DModified(ray_position_x(ray_mask),...
    ray_position_y(ray_mask), grid_x, grid_y);

interpolation_matrix = sparse(gridpoint_indices(:), tri(:), bc(:), nnz(nan_grid_binary), num_raypoints);

% interpolate the parameters on the rays' points to the grid points
if do_absorption
    
    % include acoustic absorption
    parameters_grid = interpolation_matrix * [ray_time(ray_mask),...
        ray_attenuation(ray_mask), ray_absorption(ray_mask)];
    
else
    
    % exclude the acoustic absorption
    parameters_grid = interpolation_matrix * [ray_time(ray_mask), ray_attenuation(ray_mask)] ;
    
end



% add the relative sound speed to the parameters
parameters_grid = [parameters_grid, sound_speed_relative(nan_grid_binary)];





if do_get_direction
    
    % ensure the ray directions be unit
    ray_direction = 1./sqrt(ray_direction_x(ray_mask).^ 2 + ray_direction_y(ray_mask).^2) .* [ray_direction_x(ray_mask),...
        ray_direction_y(ray_mask)];
    
    % get the Cartesian directions along the rays
    directions_grid = interpolation_matrix * ray_direction;
    
    % ensure the ray directions on the grid points be unit
    directions_grid(:, 1:2) = 1./sqrt(directions_grid(:, 1).^2 + directions_grid(:, 2).^2) .*...
        directions_grid(:, 1:2);
    
else
    directions_grid = [];
end



%if strcmp(attenuation_geom_method, 'auxiliary')

% remove the number of caustics on the emission points
caustic_rays(:, 1) = [];

if strcmp(field_mode, 'adjoint') &&  strcmp(attenuation_geom_method, 'raylinked')...
        &&   any( strcmp(greens_method, {'time', 'frequency2'}) )
    caustic_rays([1, end], :) = [];
end

caustic_rays = caustic_rays(~nan_binary);

% interpolate the integer number
% of changes in the rays' Jacobian before reaching that point
% from rays to the grid.
% (We used linear interpolation, but isn't neighboring interpolation
% more suitable than linear interpolation for this case?
caustic_number = round(interpolation_matrix *...
    caustic_rays(ray_mask));


if ~any(caustic_number)
    caustic_number = 0;
end


switch field_mode
    case 'forward'
        
        % get the integer number of caustics on the receivers
        caustic_receiver = caustic_rays(receiver_binaries);
        if ~any(caustic_receiver)
            caustic_receiver = 0;
        end
        
    case 'adjoint'
        caustic_receiver = [];
end



%else

%  caustic_numbers = [];
%  caustic_number_transducers = [];

%end



switch field_mode
    
    case 'forward'
        
        % get the parameters on the receivers
        if do_absorption
            
            % include the acoustic absorption
            parameters_receiver = [ray_time(receiver_binaries),...
                ray_attenuation(receiver_binaries), ray_absorption(receiver_binaries)];
            
        else
            
            % exculde the acoustic absorption
            parameters_receiver = [ray_time(receiver_binaries), ray_attenuation(receiver_binaries)];
            
        end
        
        
    case 'adjoint'
        
        % allocate empty variables
        parameters_receiver = [];
        % directions_receiver = [];
        
end

if ~isempty(interp_receiver)
    
    % The last point of the rays in the the polar coordinate
    % This will be used if the rays are not traced by ray linking, and the
    % rays are initialised by an equal angular spacing, and they are then
    % interpolated to the receivers using a 1D linear interpolation on the ring
    interp_receiver.lastpoint_angles = cart2pol(ray_position_x(receiver_binaries),...
        ray_position_y(receiver_binaries));
    
end


end