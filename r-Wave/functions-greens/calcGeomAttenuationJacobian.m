function [jacobian, caustic_rays, direction_x, direction_y, direction_z] =...
    calcGeomAttenuationJacobian(position_x, position_y, position_z, position_x_left,...
    position_y_left, position_z_left, position_x_right, position_y_right, position_z_right,...
    nan_binary, do_get_direction)
%CALCGEOMATTENUATIONJACOBIAN computes the Jacobian of rays
%
%   % DESCRIPTION:
%       calcGeomAttenuationJacobian computes the rays' Jacobian and the
%       cumulative number of caustics along the rays. The rays' Jacobian
%       will be used for computing the geometrical portion of the
%       attenuation (amplidue decay) which includes refraction effects.
%       The cumulative number of caustics along the rays determines the number
%       of pi/2 shifts must be applied to the phase on the rays' points.
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
%       nan_binary            - the binary mask, which is true for the
%                               rays's points outside the detection surface (ring)
%       do_get_direction      - Boolean controllong whether the cartesian
%                               direction of the rays is given as an output or not.

%       
% OUTPUTS:
%       jacobian              - the absolute Jacobian of the rays on the
%                               sampled points 
%       caustic_rays          - the cumulative number of caustics along
%                               the rays
%       direction_x           - the x Cartesian direction of the rays
%       direction_y           - the y Cartesian direction of the rays
%       direction_z           - the z Cartesian direction of the rays

%       
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 07.04.2021
%
% This script is part of the r-Wave toolbox 
% Copyright (c) 2022 Ashkan Javaherian


% get the number of dimensions
if isempty(position_z)
    dim = 2;
else
    dim = 3;
end

if isempty(position_x_left)
    
    % If the position of the left auxiliary ray is empty, the approach for
    % computing the rays' jacobian uses only the linked rays, and auxiliary
    % rays are not used. In addition, computing rays' Jacobian using only
    % the linked rays is not applicable to 3D case.
    if dim == 3
        error(['The geometrical attenuation for 3D case must be computed using'...
            'auxiliary rays approximated using paraxial rays tracing.']);
    end
    
    % Compute the derivative of x position of the points on each linked ray
    % with respect to the initial angle of the ray using the x position of
    % the equivalent points on the neighboring left and right linked rays
    perturbed_x_pos_angle  = 1/2 * ([position_x(2:end, :); position_x(1,:)]...
        - [position_x(end ,:); position_x(1:end-1,:)]);
    perturbed_x_pos_angle(nan_binary) = nan;
    
    % Compute the derivative of y position of the points on each linked ray
    % with respect to the initial angle of the ray using the y position of
    % the equivalent points on the neighboring left and right linked rays
    perturbed_y_pos_angle  = 1/2 * ([position_y(2:end, :); position_y(1, :)]...
        - [position_y(end ,:); position_y(1:end-1,:)]);
    perturbed_y_pos_angle(nan_binary) = nan;
    
else
    
    % For paraxial rays, the below will give
    % perturbed_x_angle =  position_x_left
    if isempty(position_x_right)
        position_x_right = -position_x_left;
    end
     if isempty(position_y_right)
        position_y_right = -position_y_left;
    end
    
    % compute the perturbed x position of the points on each linked ray
    % because of a perturbation to the initial angle of the ray
    % (for 3D case, the cross product of the left and right auxiliary rays)
    perturbed_x_pos_angle =  1/2 * (position_x_left - position_x_right);
    
    % compute the perturbed y position of the points on each linked ray
    % because of a perturbation to the initial angle of the ray
    % (for 3D case, the cross product of the left and right auxiliary rays)
    perturbed_y_pos_angle =  1/2 * (position_y_left - position_y_right);
    
    if dim == 3
        
        if isempty(position_z_right)
            position_z_right = -position_z_left;
        end
        
        % compute the perturbed z position of the points on each linked ray
        % because of a perturbation to the initial angle of the ray
        % (for 3D case, the cross product of the left and right auxiliary rays)
        perturbed_z_pos_angle =  1/2 * (position_z_left - position_z_right);
        
    end
  
    

end

% compute the derivative of the x position of the points on each linked ray
% with respect to the arc length of the ray
%direction_x = [(position_x(:, 2) - position_x(:, 1)),...
%    1/2 * (position_x(:, 3:end) - position_x(:, 1:end-2)),...
%    position_x(:, end) - position_x(:, end-1)];

% compute the derivative of the y position of the points on each linked ray
% with respect to the arc length of the ray
%direction_y = [(position_y(:, 2) - position_y(:, 1)),...
%    1/2 * (position_y(:, 3:end) - position_y(:, 1:end-2)),...
%    position_y(:, end) - position_y(:, end-1)];


direction_x = [position_x(:, 2:end) - position_x(:, 1:end-1), position_x(:, end) - position_x(:, end-1)];
direction_y = [position_y(:, 2:end) - position_y(:, 1:end-1), position_y(:, end) - position_y(:, end-1)];


switch dim
    case 2
        
        % compute the the rays' Jacobian
        jacobian = perturbed_x_pos_angle .* direction_y ...
            - perturbed_y_pos_angle .* direction_x;
    case 3
        
        % compute the derivative of the z position of the points on each linked ray
        % with respect to the arc length of the ray
        direction_z =  [(position_z(:, 2) - position_z(:, 1)),...
            1/2 * (position_z(:, 3:end) - position_z(:, 1:end-2)),...
            position_z(:, end) - position_z(:, end-1)];
        
        % compute the rays' Jacobian as the inner product of the rays'
        % direction and the cross product of the perturbation vectors
        % because of perturbation to the initial angles of the ray
        jacobian = perturbed_x_pos_angle .* direction_x ...
            + perturbed_y_pos_angle .* direction_y...
            + perturbed_z_pos_angle .* direction_z;
        
        
end


if ~ do_get_direction
    
    % if perturbed direction of the rays from perturbation to the initial positions
    % is not needed to compute, make the directions empty varibale.
    direction_x = [];
    direction_y = [];
    
end

if dim == 2
    direction_z = [];
end
    
%%=========================================================================
% COMPUTE THE ACCUMULATED NUMBER OF CAUSTICS
%==========================================================================
% It should be emphasized that because the refractive index (or wavenumber) is
% smoothed, the determination of caustics should be very approximate.
% compute the accumulated number of changes in the sign of the Jacobian
% along the rays

% allocate a matrix
jacobian_changes = [zeros(size(jacobian, 1), 1),...
    jacobian(:, 2:end).* jacobian(:, 1:end-1) < 0];

% get the rays, for which the Jacobian changes at least once
caustic_indices = any(jacobian_changes, 2);

% get the number of rays with caustic
num_caustic = nnz(caustic_indices);

% allocate a zero matrix for the accumulated number of caustics on the
% rays (rows) and points on the rays (columns)
caustic_rays = zeros(size(jacobian));

if  num_caustic > 0
    
    % display for testing
    disp(['The number of the rays for which the sign of the Jacobian is changed are:'...
      num2str(num_caustic)])
    
    % Get the Jacobian changes for the rays with at least one caustic
    jacobian_caustic = jacobian_changes(caustic_indices, :);
    
    % allocate a matrix for the accumulated number of caustics on the points on the
    % rays with at least one caustic
    caustic_number_matrix = zeros(size(jacobian_caustic));
    
    % get the index of rays and the index of points on the rays inclusing
    % caustics
    [caustic_rows, caustic_columns] = find(jacobian_caustic);
    
    % compute the accumulated number of changes in the Jacobian of the rays
    % the acumulated number of caustics along the rays
    for i = 1: length(caustic_rows)
         caustic_number_matrix(caustic_rows(i), caustic_columns(i):end) =...
             caustic_number_matrix(caustic_rows(i), caustic_columns(i):end) + 1;
    end
    
    % get the required number of pi/2 shifts on rays' sampled points based
    % on accumulated number of changes in sign of rays' Jacobian
    caustic_rays(caustic_indices, :) = caustic_number_matrix;
    
end

% get the absolute Jacobian 
jacobian = abs(jacobian);

% make the values ouside the binary mask nan
jacobian(nan_binary) = nan;

end