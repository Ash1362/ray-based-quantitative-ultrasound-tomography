function [parameters_grid, parameters_receiver, nan_grid_binary, caustic_number,...
    receiver_order, caustic_receiver, directions_grid, interp_receiver] = calcParametersGreens(...
    ray_position, ray_time, ray_absorption, ray_position_left, ray_position_right,...
    sound_speed_relative, grid_pos, ray_spacing, detec_radius, mask_radius_relative,...
    field_mode, do_get_direction, interp_receiver)
%CALCPARAMETERSGREENS computes the parameters of the Green's function on the
% grid points and receievers
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
%                              (linked) rays in all Cartesian coordinates
%       ray_time             - the accumulated time delays [s] along the
%                              main (linked) rays, i.e., the acoustic length
%                              divided by the reference sound speed
%       ray_absorption       - the accumulated acoustic absorption along the
%                              main (linked) rays
%       ray_position_left    - the position of the left auxiliary rays in all
%                              Cartesian coordinates
%       ray_position_right   - the position of the right auxiliary rays in all
%                              Cartesian coordinates
%       sound_speed_relative - the sound speed on the grid points over the
%                              sound speed in water
%       grid_pos             - the position [m] of the grid points
%       ray_spacing          - the scalar value representing the ray
%                              spacing [m]
%       detec_radius         - the radius [m] of the detection ring
%                              (surface)
%       mask_radius_relative - the radius of the mask over the radius of
%                              the detection ring (surface)
%       field_mode           - the mode for approximating the pressure. This
%                              can be 'forward' or 'adjoint'
%       do_get_direction     - Boolean controllong whether the cartesian
%                               direction of the ray is computed or not.
%       interp_receiver      - a struct containing the interpolation
%                              operator for lineraly interpolating the
%                              last point of the rays, i.e., the interception
%                              point of the ray by the detection ring (or
%                              surface), onto the the reception points.
%                              This is used if the ray linking is not done,
%                              and the rays are traced with equidistant initial
%                              angles. Avoiding ray linking is not recommended
%                              strongly, but included here as a benchmark
%                              for showing the effectiveness of the ray
%                              linking on improving accuracy of the approximated
%                              pressure field on the receivers.
%                              This contains the fields:
%       'angles'             - the angular position of the receivers
%       'indices'            - the ascending order of the angular poition of
%                              the receivers

% OUTPUT
%       parameters_grid      - a matrix with rows the grid points inside
%                              the binary mask, and columns: time delays,
%                              geoemetrical attenutaion, acoustic
%                              absorption (if not zero), and relative sound speed
%       parameters_receivers - a matrix with rows the receivers, and columns:
%                              time delays, geoemetrical attenutaion, acoustic
%                              absorption (if not zero), and relative sound speed
%       nan_grid_binary      - a binary mask indicating to the grid points
%                              outside the given binary mask
%       caustic_number       - the cumulative integer number of the caustics
%                              along the rays
%       receiver_order       - a vector indicating the order of receivers in
%                              the variable 'parameters_receivers'.
%       caustic_receiver     - the cumulative integer number of caustics before
%                              the ray is intercepted by the reception point
%       direction_grid       - the Cartesian direction of the rays
%       interp_receiver      - A struct containing the parameters for linearly
%                              interpolating the last point of the rays on the
%                              detection ring (surface) onto the reception points.
%                              An empty variable, if not applied.
%                              This contains the fields:
%       'angles'             - the angular position of the receivers
%       'indices'            - the ascending order of the angular poition of
%                              the receivers
%       'lastpoint_angles'   - the angular position of the last point of the
%                              rays
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 30.03.2020
%
% This script is part of the r-Wave toolbox
% Copyright (c) 2022 Ashkan Javaherian


if isscalar(ray_absorption)

    % if the absorption coefficient on the ray is given as a scalar, the acoustic
    % absorption and dispersion will be neglected.
    do_absorption = false;
    ray_absorption = 0;
else

    do_absorption = true;
end

% determine the method for computing the rays' Jacobian, which is used
% for computing the geometrical attenuation
if isempty(ray_position_left)

    % compute the ray's Jacobian with respect to the initial angle
    % using only the linked (main) rays and without the need for auxiliary
    % rays using a centred finite difference scheme
    attenuation_geom_method = 'raylinked';

    % set the auxiliary method nan
    auxiliary_method = nan;

    % the rays' direction perturbations are not required.
    do_get_direction_perturbation = false;

else

    % compute the ray's Jacobian with repect to the initial angle
    % using the auxiliary rays (recommended)
    attenuation_geom_method = 'auxiliary';

    if isempty(ray_position_right)

        % set the auxiliary method 'paraxial'
        auxiliary_method = 'paraxial';

        % the rays' direction perturbations are not required.
        do_get_direction_perturbation = false;

    else

        if (size(ray_position_right, 1) == size(ray_time, 1))

            % set the auxiliary method 'paraxial'
            auxiliary_method = 'paraxial';

            % using the paraxial ray, the input variable ray_position_right is used
            % for storing the perturbed angles of rays' direction because of a
            % perturbation to the initial position. This matrix should be
            % num_receiver x num_raypoints, the same as size of the input variable
            % ray_time
            do_get_direction_perturbation = true;

        else

            if (size(ray_position_right, 1) == size(ray_position_left, 1))

                % set the auxiliary method 'angle perturbation'
                auxiliary_method = 'angle_perturbation';

                % the rays' direction perturbations are not required.
                do_get_direction_perturbation = false;

            else

                % give an error becuse the size of the input variable
                % ray_position_right is not correct.
                error('The size of ray_position_right is not correct.')
            end

        end

    end

end


if ~do_get_direction

    % if the rays' directions are not required, the rays' directions
    % perturbations due to the perturbation to the initial positions
    % are not required as well
    do_get_direction_perturbation = false;

end




% get the number of dimensions using the matrix for cartesian position of
% the grid points
dim = size(grid_pos, 2);

% get the number of dimensions from the cartesian ray poistions
if rem(size(ray_position, 1), dim) > 0
    error(['The number of dimensions for the Cartesian position of the grid points'...
        'and that for the position of the sampled points along the rays are not consistent.'])
end

% get the position of the points on the rays along each Cartesian coordinate
ray_position_x = ray_position(1:dim:end, :);
ray_position_y = ray_position(2:dim:end, :);
if dim == 3
    ray_position_z = ray_position(3:dim:end, :);
end

clear ray_position

switch attenuation_geom_method

    case 'auxiliary'

        % For 3D , get the position of the points along the left auxiliary rays
        % along each Cartesian coordinate
        % For 2D, get the corrss product of the two auxiliary rays along each
        % cartesian coordinate
        ray_position_x_left = ray_position_left(1:dim:end, :);
        ray_position_y_left = ray_position_left(2:dim:end, :);
        if dim == 3
            ray_position_z_left = ray_position_left(3:dim:end, :);
        end

        clear ray_position_left


        switch auxiliary_method

            case 'angle_perturbation'

                % get the position of the points along the right auxiliary rays along
                % each Cartesian coordinate (for 2d)
                ray_position_x_right = ray_position_right(1:dim:end, :);
                ray_position_y_right = ray_position_right(2:dim:end, :);
                if dim == 3
                    ray_position_z_right = ray_position_right(3:dim:end, :);
                end

                clear ray_position_right

            case 'paraxial'

                % allocate empty variables to the position of the right
                % auxiliary ray, if the auxiliary method is paraxial
                ray_position_x_right = [];
                ray_position_y_right = [];
                if dim == 3
                    ray_position_z_right = [];
                end

        end

    case 'raylinked'

        % allocate empty variables to the position of points on the
        % auxiliary rays, if the method for computing the ray's Jaacobian is
        % raylinked.
        % left auxiliary rays
        ray_position_x_left = [];
        ray_position_y_left = [];

        % right_auxiliary rays
        ray_position_x_right = [];
        ray_position_y_right = [];

        if dim == 3

            % the auxiliary rays along the z coordinate
            ray_position_z_left = [];
            ray_position_z_right = [];
        end
end


if dim == 2
    ray_position_z = [];
    ray_position_z_left = [];
    ray_position_z_right = [];
end


%%=========================================================================
% COMPUTE THE GEOMETRICAL ATTENUATION
%==========================================================================
% compute the geometrical attenuation along the rays
[ray_attenuation, nan_binary, receiver_indices, caustic_rays,...
    ray_direction_x, ray_direction_y, ray_direction_z] = calcGeomAttenuationRays(...
    ray_position_x, ray_position_y, ray_position_z, ray_position_x_left,...
    ray_position_y_left, ray_position_z_left, ray_position_x_right,...
    ray_position_y_right, ray_position_z_right, ray_time, ray_spacing,...
    do_get_direction);


% discard the emission point from the rays' parameters
% x
ray_position_x(:, 1) = [];
% y
ray_position_y(:, 1) = [];
% z
if ~isempty(ray_position_z)
    ray_position_z(:, 1) = [];
end

% time delays
ray_time(:, 1) = [];

% geometrical attenuation
ray_attenuation(:, 1) = [];

% nans
nan_binary(:,1) = [];

% acoustic absorption
if do_absorption
    ray_absorption(:, 1) = [];
end

% caustics
caustic_rays(:, 1) = [];

% rays' directions
if ~isempty(ray_direction_x)
    ray_direction_x(:, 1) = [];
end

if ~isempty(ray_direction_y)
    ray_direction_y(:, 1) = [];
end

if ~isempty(ray_direction_z)
    ray_direction_z(:, 1) = [];
end

if do_get_direction_perturbation
    ray_position_right(:, 1) = [];
end

switch field_mode
    case 'forward'

        % compute the linear indices of receivers, considering that the
        % emission point (the first column) has been discarded from the rays
        %receiver_indices = sub2ind(size(ray_position_x), 1:size(ray_position_x, 1),...
        %    receiver_indices' - 1);
        receiver_indices = sub2ind(size(ray_position_x), 1:size(ray_position_x, 1),...
            receiver_indices') - size(ray_position_x, 1);

        [receiver_binaries, receiver_order]= ismember(find(~nan_binary), receiver_indices);

        % get the order of receivers
        receiver_order = receiver_order(receiver_binaries);

    case 'adjoint'


        %if strcmp(attenuation_geom_method, 'raylinked')

        % because the adjoint field is calculated only for a cone with axis
        % the receiver-emitter pair, the geometrical attenuation for the
        % surface of the cone cannot be comuted, and must be discarded.
        % Note that for computing the adjoint field, distance_integer_rays is
        % always set zero
        %ray_position_x([1, end], :) = [];
        %ray_position_y([1, end], :) = [];
        %ray_time([1, end], :) = [];
        %ray_attenuation([1, end], :) = [];
        %nan_binary([1, end], :) = [];
        %if do_absorption
        %    ray_absorption([1, end], :) = [];
        %end

        %end

        receiver_order = [];


end

% discard the nans in the matrices for the rays' parameters
% x
ray_position_x = ray_position_x(~nan_binary);
% y
ray_position_y = ray_position_y(~nan_binary);
% z
if ~isempty(ray_position_z)
    ray_position_z = ray_position_z(~nan_binary);
end

% time delays
ray_time = ray_time(~nan_binary);

% geometrical attenuation
ray_attenuation = ray_attenuation(~nan_binary);

% accumulated acoustic absorption
if do_absorption
    ray_absorption = ray_absorption(~nan_binary);
end

% caustics
caustic_rays = caustic_rays(~nan_binary);

% rays' directions
if ~isempty(ray_direction_x)
    ray_direction_x = ray_direction_x(~nan_binary);
end

if ~isempty(ray_direction_y)
    ray_direction_y = ray_direction_y(~nan_binary);
end

if ~isempty(ray_direction_z)
    ray_direction_z = ray_direction_z(~nan_binary);
end

if do_get_direction_perturbation
    ray_position_right = ray_position_right(~nan_binary);
end


% get the binary mask for the sampled points on the rays
switch dim
    case 2
        ray_mask = ray_position_x.^2 + ray_position_y.^2 < (mask_radius_relative * detec_radius)^2;
    case 3
        ray_mask = ray_position_x.^2 + ray_position_y.^2 + ray_position_z.^2 < ...
            (mask_radius_relative * detec_radius)^2;
end

% get the number of ray points inside the defined mask
num_raypoints = nnz(ray_mask);

%%=========================================================================
% CONSTRUCT THE RAY-TO-GRID INTERPOLATION MATRIX
%==========================================================================
% create the sparse interpolation matrix for interpolating the parameters
% of Green's function on rays' points onto the grid points
switch dim
    case 2

        % get the indices and coefficients for constructing ray-to-grid interpolation matrix
        [gridpoint_indices, tri, bc, nan_grid_binary] = gridDataFast2DModified(...
            ray_position_x(ray_mask), ray_position_y(ray_mask),...
            grid_pos(:, 1), grid_pos(:, 2));

    case 3

        % get the indices and coefficients for constructing ray-to-grid interpolation matrix
        [gridpoint_indices, tri, bc, nan_grid_binary] = gridDataFast3DModified(...
            ray_position_x(ray_mask), ray_position_y(ray_mask),...
            ray_position_z(ray_mask), grid_pos(:, 1), grid_pos(:, 2),...
            grid_pos(:, 3));
end

% construct the interpolation matrix
interpolation_matrix = sparse(gridpoint_indices(:), tri(:), bc(:),...
    nnz(nan_grid_binary), num_raypoints);

% interpolate the parameters on the rays' points onto the grid points
if do_absorption

    % get the parameters of Green's function on the grid points
    % (include the accumulated acoustic absorption)
    parameters_grid = interpolation_matrix * [ray_time(ray_mask),...
        ray_attenuation(ray_mask), ray_absorption(ray_mask)];

else

    % get the parameters of Green's function on the grid points
    % (ignore the accumulated acoustic absorption)
    parameters_grid = interpolation_matrix * [ray_time(ray_mask), ray_attenuation(ray_mask)] ;

end


% add a column containing the relative sound speeds to the matrix for
% parameters of the Green's function on the grid points
parameters_grid = [parameters_grid, sound_speed_relative(nan_grid_binary)];


if do_get_direction

    switch dim
        case 2

            % get the matrix for rays' directions
            ray_direction = [ray_direction_x(ray_mask), ray_direction_y(ray_mask)];

        case 3

            % get the matrix for rays' directions
            ray_direction = [ray_direction_x(ray_mask), ray_direction_y(ray_mask),...
                ray_direction_z(ray_mask)];

    end

    % make the rays' directions unit vector
    ray_direction = 1./vecnorm(ray_direction, 2, 2) .* ray_direction;

    % interpolate the rays' directions from rays to the grid points
    directions_grid = interpolation_matrix * ray_direction;

    % make the interpolated rays' directions on the grid points unit
    % vectors
    directions_grid = 1./ vecnorm(directions_grid, 2 , 2) .*...
        directions_grid;

    if do_get_direction_perturbation

        % add acolumn for ray's direction perturbations in the polar coordinate
        % due to perturbations to initial positions of the rays to the
        % ray's directions in the Cartesian coordinates
        directions_grid = [directions_grid, interpolation_matrix * ray_position_right(ray_mask)];

    end

else

    directions_grid = [];

end



%if strcmp(attenuation_geom_method, 'auxiliary')

% remove the number of caustics on the emission points


%if strcmp(field_mode, 'adjoint') &&  strcmp(attenuation_geom_method, 'raylinked')...
%    &&   any( strcmp(greens_method, {'time', 'frequency2'}) )
% caustic_rays([1, end], :) = [];
% end

% caustic_rays = caustic_rays(~nan_binary);

% interpolate the integer acccumulated number of changes in the sign of rays' Jacobian
% from rays to the grid. (This interpolation makes sense, because the trajectory of
% rays are computed on the smoothed refractive index.)
caustic_number = round(interpolation_matrix *...
    caustic_rays(ray_mask));

%if ~any(caustic_number)
%    caustic_number = 0;
%end


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

        % get the parameters of Green's function on the receivers
        if do_absorption

            % consruct the matrix (include the acoustic absorption)
            parameters_receiver = [ray_time(receiver_binaries),...
                ray_attenuation(receiver_binaries), ray_absorption(receiver_binaries)];

        else

            % consruct the matrix (ignore the acoustic absorption)
            parameters_receiver = [ray_time(receiver_binaries), ray_attenuation(receiver_binaries)];

        end


    case 'adjoint'

        % allocate empty variables
        parameters_receiver = [];


end

if ~isempty(interp_receiver)

    switch dim
        case 2

            % The last point of the rays in the the polar coordinate
            % This will be used if the rays are not traced by ray linking, and the
            % rays are initialised by an equal angular spacing, and they are then
            % interpolated to the receivers using a 1D linear interpolation on the ring
            interp_receiver.lastpoint_angles = cart2pol(ray_position_x(receiver_binaries),...
                ray_position_y(receiver_binaries));
        case 3
            error('Not implemented yet!')
    end

end


end


