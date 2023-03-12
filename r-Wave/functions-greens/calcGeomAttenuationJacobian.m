function [jacobian_absolute, caustic_rays, direction_x, direction_y]...
    = calcGeomAttenuationJacobian(position_x, position_y, position_x_left,...
    position_y_left, position_x_right, position_y_right, nan_binary,...
    distance_integer_rays, do_get_direction)
%CALCGEOMATTENUATIONJACOBIAN computes the Jacobian of rays
%
%   % DESCRIPTION:
%       calcGeomAttenuationJacobian computes the absolute Jacobian of rays
%       and the cumulative number of Caustics on the rays' points. The
%       absolute Jacobian of the rays will be used for calculating the
%       geomterical portion of the attenuation (amplidue decay) which
%       includes the refraction effects. The cumulative number of caustics 
%       along the rays determines the number of pi/2 shifts must be applied 
%       to the computed pressure field on the rays' points.
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
%       nan_binary            - the binary mask, which is true for the
%                               rays's points outside the detection surface (ring)
%       distance_integer_rays - the number of rays between the main linked ray
%                               and the two other linked rays used for
%                               computing the Jacobian of the main ray,
%                               instead of the auxiliary rays. Not applied, if the
%                               method for computing the ray's Jacobian is
%                               set 'auxiliary'
%       do_get_direction      - Boolean controllong whether the cartesian
%                               direction of the ray is given as an output or not.

%       
% OUTPUTS:
%       jacobian_absolute     - the absolute value of the Jacobian of the
%                               rays on the smapled points along the rays
%       caustic_rays          - the cumulative number of caustics alosng
%                               the rays
%       direction_x        -  the x Cartesian direction of the rays
%       direction_y        -  the y Cartesian direction of the rays

%       
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 30.03.2020
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian


if isempty(position_x_left)
    
    % Compute the derivative of x position of the points on each linked ray
    % with respect to the initial angle of the corresponding linked ray using the x position of
    % the equivalent points on the neighboring left and right linked rays
    jacobian_x_angle = 1/2 * ([position_x(distance_integer_rays+2:end, :); position_x(1:distance_integer_rays+1,:)]...
        - [position_x(end-distance_integer_rays:end ,:); position_x(1:end-distance_integer_rays-1,:)]);
    jacobian_x_angle(nan_binary) = nan;
    
    % Compute the derivative of y position of the points on each linked ray
    % with respect to initial angle of the corresponding linked ray using the y position of
    % the equivalent points on the neighboring left and right linked rays
    jacobian_y_angle = 1/2 * ([position_y(distance_integer_rays+2:end, :); position_y(1:distance_integer_rays+1,:)]...
        - [position_y(end-distance_integer_rays:end ,:); position_y(1:end-distance_integer_rays-1,:)]);
    jacobian_y_angle(nan_binary) = nan;
    
else
    
    % make the right auxiliary ray zero, if not used.
    if isempty(position_x_right)
        position_x_right =0;
    end
     if isempty(position_y_right)
        position_y_right =0;
    end
    
    
    
  %  if ~strcmp(auxiliary_method, 'paraxial')
    % Compute the derivative of x position of the points on each linked ray
    % with respect to the initial angle of the corresponding linked ray using the x position of
    % the equivalent points on the left and right auxiliary rays
    jacobian_x_angle =  1/2 * (position_x_right - position_x_left);
    
    % Compute the derivative of y position of the points on each linked ray
    % with respect to the initial angle of the corresponding linked ray using the y position of
    % the equivalent points on the left and right auxiliary rays
    jacobian_y_angle =  1/2 * (position_y_right - position_y_left);
     
  %  end
    
    
end

% Compute the derivative of x position of the points on each linked ray
% with respect to the arc length along the ray
direction_x =  [(position_x(:, 2) - position_x(:, 1)),...
    1/2 * (position_x(:, 3:end) - position_x(:, 1:end-2)),...
    position_x(:, end) - position_x(:, end-1)];

% Compute the derivative of y position of the points on each linked ray
% with respect to the arc length along the ray
direction_y =  [(position_y(:, 2) - position_y(:, 1)),...
    1/2 * (position_y(:, 3:end) - position_y(:, 1:end-2)),...
    position_y(:, end) - position_y(:, end-1)];

% direction_x = [position_x(:, 2:end) - position_x(:, 1:end-1), position_x(:, end) - position_x(:, end-1)];
% direction_y = [position_y(:, 2:end) - position_y(:, 1:end-1), position_y(:, end) - position_y(:, end-1)];

% compute the Jacobian of the ray points as the determinant of the Jacobian matrix
%if any(strcmp(auxiliary_method, {'angle_perturbation', nan}))

jacobian_determinant = jacobian_x_angle .* direction_y ...
    - jacobian_y_angle .* direction_x;
% else
    
    % compute the cross product of the auxiliary rays
    
 %   position_x_left .* position_y_right - position_y_left .* position_x_right
    
    
    % compute the inner product of the derivative with respect to the arc
    % length and the cross product of the auxiliary rays
  %  jacobian_determinant = jacobian_x_s .* auxilary_cross_product(:,:,1) ...
  %      + jacobian_y_s .* auxilary_cross_product(:,:,2);

        
% end

if ~ do_get_direction
    
    % if perturbed direction of the rays from perturbation to the initial positions
    % is not required, make the directions empty varibale.
    direction_x = [];
    direction_y = [];
    
end
    

% Find the accumulated number of changes in the sign of the Jacobian
% along the rays

% Allocate a matrix
jacobian_determinant_changes = [zeros(size(jacobian_determinant, 1), 1),...
    jacobian_determinant(:, 2:end).* jacobian_determinant(:, 1:end-1) < 0];

% Find the rays, for which the Jacobian changes at least once
caustic_indices = any(jacobian_determinant_changes, 2);

% Get the number of rays with caustic
num_caustic = nnz(caustic_indices);

% allocate a zer matrix for the accumulated number of caustics on the
% rays (rows) and points on the rays (columns)
caustic_rays = zeros(size(jacobian_determinant));

if  num_caustic > 0
    
    % display for testing
    % disp(['The number of the rays for which the sign of the Jacobian is changed are:'...
    %  num2str(num_caustic)])
    
    % Get the Jacobian changes for the rays with at least one cautic
    jacobian_determinant_caustic = jacobian_determinant_changes(caustic_indices, :);
    
    % Allocate a matrix for the accumulated number of caustics on the points on the
    % rays with at least one caustic
    caustic_number_matrix = zeros(size(jacobian_determinant_caustic));
    
    % Get the index of rays and the index of point on each ray on which
    % caustic occurs
    [caustic_rows, caustic_columns] = find(jacobian_determinant_caustic);
    
    % Compute the accumulated number of changes in the Jacobian of the rays
    for i = 1: length(caustic_rows)
        caustic_number_matrix(caustic_rows(i), caustic_columns(i):end) =...
            caustic_number_matrix(caustic_rows(i), caustic_columns(i):end) + 1;
    end
    
    % Get the number of pi/2 shift for each ray point based on the number
    % of the changes in the Jacobian of the rays on the points along the ray
    caustic_rays(caustic_indices, :) = caustic_number_matrix;
    
end

% the absolute of the Jacobian 
jacobian_absolute = abs(jacobian_determinant);

% remove values from nans (the points outside the detection surface)
jacobian_absolute(nan_binary) = nan;

end