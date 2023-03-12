function [pos, interp_coeff_vec] = calcRayDualUpdate(refractive, ray_fields_params,...
    cartesian_position_emitter, direction, xvec, yvec, zvec, pos_grid_first,...
    pos_grid_end, dx, ds, grid_size, dim, detec_radius, mask, calc_coeffs)
% CALCRAYDUAL UPDATE traces a ray using a Dual update algorithm, and stores
% the coefficients for the ray-to-grid interpolation.
%
% DESCRIPTION:
% calRayDualUpdate computes the position of the sampled points along the ray
% and their coefficientsfor an interpolation onto the grid. The rays are traced
% given an initial position and an initial direction for the ray.
% The interpolation coefficients will be used for an integration of the acoustic
% length along the ray given a refractive index distribution on the grid.
% The ray tracing is done using the Dual Update algorithm, and the ray-to-grid
% interpolation (and vice versa) is done using a 'Bilinear' or 'Bspline' approach.
% Please see Section 3.3.1., Algorithm 1, A Javaherian, F Lucka F and B T Cox,
% 2020, Refraction-corrected ray-based inversion for three-dimensional ultrasound
% tomography of the breast, Inverse Problems, Vol. 36 125010.
% USAGE:
%
%
% INPUTS:
%       refractive        - the refractive index on the grid points
%       ray_fields_params - a struct with fields the direction gradients
%                           for an interpolation using a'Bilinear'
%                           approach, or parameters for choosing indices of
%                           the grid points and their associated coefficients for an
%                           interpolation using a 'Bspline' approach.
%                           The fields are:
%       'refractive_gradient_x' - discretised refractive index gradient along x
%       'refractive_gradient_y' - discretised refractive index gradient along y
%       'refractive_gradient_z' - discretised refractive index gradient along z
%                              using a'Bspline' interpolation, this
%                              includes the matrices for intertpolation:
%       'raytogrid_indices_x'   - x indices for B-spline interpolation
%       'raytogrid_indices_y'   - y indices for B-spline interpolation
%       'raytogrid_indices_z'   - z indices for B-spline interpolation
%     'raytogrid_coeff_matrix'  - matrix for calculating B-spline
%                                 interpolation coefficients of the field
%     'raytogrid_coeff_derivative_matrix' - matrix for calculating B-spline
%                                 interpolation coefficients of the directional
%                                 gradientsof the field
%       cartesian_position_emitter - a dim x 1 Cartesian position of the emitter
%       direction        - a dim x 1 vector of the cartesian initial direction along
%                         the ray
%       xvec             - the x vector of the grid points
%       yvec             - the y vector of the grid points
%       zvec             - the z vector of the grid points
%       pos_grid_first   - a dim x 1 Cartesian position of the first index of
%                         the grid
%       pos_grid_end     - a dim x 1 Cartesian position of the end index of
%                         the grid
%       dx               - a scalar representing the grid spacing, the same
%                          along all the Cartesian coordinates [m]
%       ds               - a scalar representing the ray spacing [m]
%       grid_size        - the size of the grid
%       dim              - the number of dimensions of the grid
%       detec_radius     - if it is scalar value, the radius of the detection
%                          ring (2D), or detection surface(3D), if it is
%                          string '1', or '2', it will show the scenario
%                          for testing the ray tracing algorithm using a
%                          Maxwell fish-eye lens phantom. The scenario '1'
%                          and '2' is for measuring the accuracy in computing the
%                          acoustic length, and deviation of the rays'
%                          trajectory from an expected circular trajectory.
%       mask             - a binary mask used for calculating the
%                          coefficients. The sound speed ouside the binary
%                          mask is assumed homogeneous water.
%       calc_coeffs      - a Boolean controlling whether the
%                          interpolation coefficients are stored or
%                          not. This Boolean is set true only for a linked
%                          ray after solving the ray linking problem.

% OPTIONAL INPUTS:
%
% OUTPUTS:
%
%      pos              - a dim x 1 vector of the cartesian position of the
%                         end point along the ray
%      interp_coeff_vec - the vector of interpolation coefficients for the
%                         rays' points along the ray within the binary mask.
%                         Only if the purpose is testing the ray tracing algorithm
%                         and using scenario for measuring the deviation from a
%                         circular trajectory using a Maxwell’s ‘fish-eye lens’
%                         phantom, this variable is used for storing the position of
%                         ray points, not interpolation coefficients.

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.06.2022
%
% This script is part of the r-Wave Toolbox
% Copyright (c) 2022 Ashkan Javaherian


% get the Boolean controlling whether the task is testing ray tracing, or
% it is construction of the system matrix for image reconstruction
switch class(detec_radius)
    case 'double'
        test_raytracing = false;
    case 'char'
        test_raytracing = true;
        if ~strcmp(detec_radius, {'1', '2'})
            error('The scenario must be a string as 1 or 2.')
        end
end


if isfield(ray_fields_params, 'refractive_gradient_x')
    
    % get the grid-to-ray interpolation method
    interp_method = 'Bilinear';
    
    % get the gradient of the refractive index field along the x coordinate
    refractive_gradient_x = ray_fields_params.refractive_gradient_x;
    % get the gradient of the refractive index field along the x coordinate
    refractive_gradient_y = ray_fields_params.refractive_gradient_y;
    % get the gradient of the refractive index field along the x coordinate
    refractive_gradient_z = ray_fields_params.refractive_gradient_z;
    
else
    
    % get the grid-to-ray interpolation method
    interp_method = 'Bspline';
    
    % get parameters for the Bspline interpolation
    raytogrid_indices_x = ray_fields_params.raytogrid_indices_x;
    raytogrid_indices_y = ray_fields_params.raytogrid_indices_y;
    raytogrid_indices_z = ray_fields_params.raytogrid_indices_z;
    raytogrid_coeff_matrix = ray_fields_params.raytogrid_coeff_matrix;
    raytogrid_coeff_derivative_matrix = ray_fields_params.raytogrid_coeff_derivative_matrix;
    
end



if calc_coeffs
    
    if ~test_raytracing || strcmp(detec_radius, '1')
        
        % allocate a sparse vector for storing the interpolation coefficients
        interp_coeff_vec = sparse(prod(grid_size), 1);
    else
        
        interp_coeff_vec = [];
        
    end
    
    
end

% remove the last index from the size of the grid (not used)
grid_size = grid_size(1: dim-1);

% normalise direction - the direction along the ray must be a unit vector
direction =  1/norm(direction) * direction;

% initialise the position of the ray with the emission point
pos = cartesian_position_emitter;

if test_raytracing
    switch detec_radius
        case '1'
            if dim == 2
                cartesian_position_emitter = -cartesian_position_emitter;
            end
        case '2'
            interp_coeff_vec = pos;
    end
end


switch interp_method
    case 'Bilinear'
        [indices, coeff, n, dn] = interpLocal(pos, xvec, yvec, zvec, ...
            pos_grid_first, dx, grid_size, dim, mask, refractive,...
            refractive_gradient_x, refractive_gradient_y, refractive_gradient_z);
    case 'Bspline'
        [indices, coeff, n, dn] = interpLocalBspline(pos, xvec, yvec, zvec, ...
            pos_grid_first, dx, grid_size, dim, mask, refractive,...
            raytogrid_indices_x, raytogrid_indices_y,...
            raytogrid_indices_z, raytogrid_coeff_matrix, ...
            raytogrid_coeff_derivative_matrix, true);
end

if isempty(coeff)
    h = 0;
else
    
    % compute the second-order term
    h =  ds/n * (dn - (dn' * direction) * direction);
    
    % integrate the refractive index along the ray to calculate the
    % acoustic length, if requested
    if  calc_coeffs && ~strcmp(detec_radius, '2')
        interp_coeff_vec(indices) = interp_coeff_vec(indices) + 1/2 * coeff';
    end
    
end


%  update direction
g = direction + h/2;

% normalise the update direction
g = 1/norm(g) * g;

% store the current position
pos_previous = pos;

% update the current position of the ray
pos = pos_previous + ds * g;



if strcmp(detec_radius, '2')
    
    % get the position (interp_coeff_vec is used as position for this
    % case.)
    interp_coeff_vec = [interp_coeff_vec, pos];
end


% if the purpose is not testing ray tracing, the ray must be terminated inside
% the detection area (volume). if the purpose is testing ray tracing, the ray is
% terminated when the distance of the current point of the ray with the expected
% last point is smaller than ray spacing.
% For the latter termination criterion, the B-spline interpolation, for which the number of
% required neighboring grid points are more than Bilinear is considered.
while (test_raytracing && norm(pos - cartesian_position_emitter) > ds - (1e-10)  &&...
        all(pos > pos_grid_first + dx & pos < pos_grid_end - 2 * dx)) ||...
        (~test_raytracing &&  norm(pos) - detec_radius <= 1e-10    &&...
        all(pos > pos_grid_first + dx & pos < pos_grid_end - 2 * dx))
    
    % compute the interpolation coefficients, their corresponding indices
    % on the grid. Then use them for approximating the refractive index and the
    % vector of its  directional gradient on the ray point
    switch interp_method
        case 'Bilinear'
            [indices, coeff, n, dn] = interpLocal(pos, xvec, yvec, zvec, ...
                pos_grid_first, dx, grid_size, dim, mask, refractive, ...
                refractive_gradient_x, refractive_gradient_y, refractive_gradient_z);
        case 'Bspline'
            [indices, coeff, n, dn] = interpLocalBspline(pos, xvec, yvec, zvec, ...
                pos_grid_first, dx, grid_size, dim, mask, refractive,...
                raytogrid_indices_x, raytogrid_indices_y,...
                raytogrid_indices_z, raytogrid_coeff_matrix, ...
                raytogrid_coeff_derivative_matrix, true);
    end
    
    
    
    % update the direction
    direction = direction + h;
    
    % normalise the direction along the ray
    direction = 1/norm(direction) * direction;
    
    
    if isempty(coeff)
        h = 0;
    else
        % calculate the secon-order term
        h =  ds/n * (dn - (dn' * direction) * direction);
        % integrate the refractive index along the ray to calculate the
        % acoustic length, if requested
        if calc_coeffs && ~strcmp(detec_radius,  '2')
            interp_coeff_vec(indices) = interp_coeff_vec(indices) + coeff';
        end
    end
    
    %  update direction
    g = direction + h/2;
    
    % normalise the update direction
    g = 1/norm(g) * g;
    
    % store the current position
    pos_previous = pos;
    
    % update the current position of the ray
    pos = pos_previous + ds * g;
    
    if strcmp(detec_radius, '2')
        
        % get the position
        interp_coeff_vec = [interp_coeff_vec, pos];
    end
    
end

if test_raytracing
    
    % set the last position as the initial position
    pos = cartesian_position_emitter;
    % correct the last ray spacing
    ds_final = norm(pos - pos_previous);
    
    if strcmp(detec_radius, '2')
        
        % get the position
        interp_coeff_vec(:, end) = pos;
    end
    
else
    
    %  the last point must be on the detection surface
    [pos , ds_final] = calcIntersectBall(pos_previous, g, detec_radius);
    
end


if calc_coeffs
    
    if ~strcmp(detec_radius, '2')
        
        % multiply by the ray spacing to get the integral of the refractive index
        % (acoustic length) along the ray
        interp_coeff_vec = ds * interp_coeff_vec;
        
        % correct the integration for the last point
        if ~isempty(coeff)
            interp_coeff_vec (indices) = interp_coeff_vec (indices) + ...
                1/2 * (ds_final - ds) * coeff';
        end
        
        if  all (pos > pos_grid_first + dx  &  pos < pos_grid_end - 2 * dx)
            % compute the interpolation coefficients, their corresponding indices
            % on the grid.
            switch interp_method
                case 'Bilinear'
                    
                    [indices, coeff, ~, ~] = interpLocal(pos, xvec, yvec, zvec, ...
                        pos_grid_first, dx, grid_size, dim, mask, refractive);
                    
                case 'Bspline'
                    
                    [indices, coeff, ~, ~] = interpLocalBspline(pos, xvec, yvec, zvec,...
                        pos_grid_first, dx, grid_size, dim, mask, refractive,...
                        raytogrid_indices_x, raytogrid_indices_y,...
                        raytogrid_indices_z, raytogrid_coeff_matrix, ...
                        raytogrid_coeff_derivative_matrix, false);
                    
            end
            
            
            if ~isempty(coeff)
                interp_coeff_vec(indices) = interp_coeff_vec(indices) +...
                    1/2 * ds_final * coeff';
            end
            
            
        end
        
    end
    
else
    
    % give an empty variable, if interpolation coefficients are not required
    interp_coeff_vec = [];
    
end



end