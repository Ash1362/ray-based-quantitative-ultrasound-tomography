function [pos, interp_coeff_vec] = calcRayRungeKutta2ndBspline(....
    refractive, cartesian_position_emitter, direction,...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds,...
    grid_size, dim, detec_radius, mask, raytogrid_indices_x,...
    raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix, ...
    raytogrid_coeff_derivative_matrix, calc_coeffs)
%CALCRAYRUNGEKUTTA2NDBSPLINE traces a ray using a secon-order Runge-Kutta algorithm,
% and stores the coefficients for the ray-to-grid interpolation. The interpolation 
% is done using a Bspline approach.
%
% DESCRIPTION:
% calRayRungeKutta2ndBspline calculates positions along the ray and their interplation
% coefficients, given an initial position and an initial direction along the
% ray. The interpolation coefficients will be used for an integration of
% the acoustic length along the ray given a discretised refractive index
% distribution. The ray tracing is done using a second-order variant of the
% Runge-Kutta algorithm, also Known as Heun's method, and the interpolations
% between the grid and rays are done using a Bspline approach.
%
% USAGE:
%
%
% INPUTS:
%       refractive       - disretised refractive index                             
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
%                         along all the Cartesian coordinates [m]
%       ds               - a saclar representing the ray spacing [m]
%       grid_size        - the size of the grid
%       dim              - the number of dimensions of the grid
%       mask             - a binary mask used for calculating the
%                          coefficients. The sound speed ouside the binary
%                          mask is assumed homogeneous water.                
%      raytogrid_indices_x -  x indices for B-spline interpolation
%      raytogrid_indices_y -  y indices for B-spline interpolation                                  
%      raytogrid_indices_z -  z indices for B-spline interpolation
%      raytogrid_coeff_matrix     - matrix for calculating B-spline
%                                   interpolation coefficients of the field
%      raytogrid_coeff_derivative_matrix - matrix for calculating B-spline
%                                         interpolation coefficients of the
%                                         directional gradients of the field 
%       calc_coeffs           - a Boolean controlling whether the
%                               interpolation coefficients are stored or
%                               not. This Boolean is set true only for a linked
%                               ray after solving the ray linking problem.
% OPTIONAL INPUTS:
%
% OUTPUTS:
%
%      pos                - a dim x 1 vector of the cartesian poistion of the 
%                           end point along the ray
%      interp_coeff_vec - the vector of interpolation coefficients for the ray's points
%                         along the ray within the binary mask

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.06.2022
%
% This script is part of the r-Wave Toolbox 
% Copyright (c) 2020 Ashkan Javaherian

if calc_coeffs
    
    % allocate a sparse vector for storing the interpolation coefficients
    interp_coeff_vec = sparse(prod(grid_size), 1);
end
% remove the last index from the size of the grid (not used)
grid_size = grid_size(1: dim-1);




% initialise the position of the ray with the emission point
pos = cartesian_position_emitter;

% calculate the interpolation coefficients and their corresponding indices
% on the grid
[indices, coeff, n, dn] = interpLocalBspline(pos, xvec, yvec, zvec, ...
    pos_grid_first, dx, grid_size, dim, mask, refractive, raytogrid_indices_x,...
    raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix, ...
    raytogrid_coeff_derivative_matrix, true);

if isempty(coeff)
    dn = 0; n = 1;
else
    
    % integrate the refractive index along the ray to calculate the
    % acoustic length, if requested
    if calc_coeffs
        interp_coeff_vec(indices) = interp_coeff_vec(indices) + 1/2 * coeff';
    end
end


% normalise direction - the direction along the ray must be a unit vector
% multiplied by refractive index
direction =  n/norm(direction) * direction;

% calculate the position component of k1
k1_x = 1/n * direction;



% calculate the position component of k2
direction_k2 = direction + ds * dn;

pos_k2 = pos + ds * k1_x;


% calculate the interpolation coefficients and their corresponding indices
% on the grid
[~, coeff_k2, n_k2, dn_k2] = interpLocalBspline(pos_k2, xvec, yvec, zvec, ...
    pos_grid_first, dx, grid_size, dim, mask, refractive, raytogrid_indices_x,...
    raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix, ...
    raytogrid_coeff_derivative_matrix, true);


if isempty(coeff_k2)
     dn_k2 = 0; n_k2 = 1;
end


% normalise direction - the direction along the ray must be a unit vector
% multiplied by refractive index
direction_k2 =  n_k2/norm(direction_k2) * direction_k2;


k2_x = 1/n_k2 * direction_k2;



% store the current position
pos_previous = pos;


% update the current position of the ray
pos = pos_previous + ds/2 * (k1_x + k2_x);


% the ray must be terminated inside the detection, and also the
% pointsincluded in the interpolation must not exceed the compuational grid
% The second stopping criterion ensures that  all the grid points included in the
% Bspline interpolation for the current point pos and the next point pos_k2
% do not exceed the edges of the computational grid.
while   norm(pos) - detec_radius <= 1e-10    &&   all(pos > pos_grid_first + 2 * dx & pos < pos_grid_end - 3 * dx)
   
    
    % update the direction
    direction = direction + ds/2 * (dn + dn_k2);
    
    
    % calculate the interpolation coefficients and their corresponding indices
    % on the grid
    [indices, coeff, n, dn] = interpLocalBspline(pos, xvec, yvec, zvec, ...
        pos_grid_first, dx, grid_size, dim, mask, refractive, raytogrid_indices_x,...
        raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix, ...
        raytogrid_coeff_derivative_matrix, true);
    
    
    if isempty(coeff)
        dn = 0; n = 1;
    else
        
        % integrate the refractive index along the ray to calculate the
        % acoustic length, if requested
        if calc_coeffs  
            interp_coeff_vec(indices) = interp_coeff_vec(indices) + coeff';
        end
    end
 
    
    
    % normalise direction - the direction along the ray must be a unit vector
    % multiplied by refractive index
    direction =  n/norm(direction) * direction;
    
    % calculate the position component of k1
    k1_x = 1/n * direction;
    
    
    % calculate the position component of k2
    direction_k2 = direction + ds * dn;
    
    pos_k2 = pos + ds * k1_x;
   
    % calculate the interpolation coefficients and their corresponding indices
    % on the grid
    [~, coeff_k2, n_k2, dn_k2] = interpLocalBspline(pos_k2, xvec, yvec, zvec, ...
        pos_grid_first, dx, grid_size, dim, mask, refractive, raytogrid_indices_x,...
        raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix, ...
        raytogrid_coeff_derivative_matrix, true);
    
    
    if isempty(coeff_k2)
         dn_k2 = 0; n_k2 = 1;
    end
    
    
    % normalise direction - the direction along the ray must be a unit vector
    % multiplied by refractive index
    direction_k2 =  n_k2/norm(direction_k2) * direction_k2;
    
    
    k2_x = 1/n_k2 * direction_k2;
    
    
    
    % store the current position
    pos_previous = pos;
    
    
    % update the current position of the ray
    pos = pos_previous + ds/2 * (k1_x + k2_x);
    
    
end


%  the last point must be on the detection surface
[pos , ds_final] = calcIntersectBall(pos_previous, 1/2 * (k1_x+k2_x), detec_radius);


if calc_coeffs
    
    % multiply by the ray spacing to get the integral of the refractive index
    % (acoustic length) along the ray
    interp_coeff_vec = ds * interp_coeff_vec;
    
    % correct the integration for the last point
    if ~isempty(coeff)
        interp_coeff_vec(indices) = interp_coeff_vec(indices) + ...
            1/2 * (ds_final - ds) * coeff';
    end
    
    if  all (pos > pos_grid_first +  dx  &  pos < pos_grid_end - 2 * dx )
             
        [indices, coeff, ~, ~] = interpLocalBspline(pos, xvec, yvec, zvec,...
            pos_grid_first, dx, grid_size, dim, mask, refractive, raytogrid_indices_x,...
            raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix, ...
            raytogrid_coeff_derivative_matrix, false);
        
        if ~isempty(coeff)
            interp_coeff_vec(indices) = interp_coeff_vec(indices) +...
                1/2 * ds_final * coeff';
        end
        
    end
else
    
    % give an empty variable, if interpolation coefficients are not required
    interp_coeff_vec = [];
    
end





end