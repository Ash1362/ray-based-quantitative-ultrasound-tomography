function [pos, interp_coeff_vec] = calcRayCharacteristicsBspline(...
    refractive, cartesian_position_emitter, direction,...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds,...
    grid_size, dim, detec_radius, mask, raytogrid_indices_x,...
    raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix, ...
    raytogrid_coeff_derivative_matrix, calc_coeffs)
%CALCRAY traces a ray and store its associated interpolation coefficients
%
% DESCRIPTION:
% calRay calculates positions along the ray and their interplation
% coefficients, given an initial position and an initial direction along the
% ray. the interpolation coefficients will be used for an integration of
% acoustic length along the ray gioven a discretised refractive index distribution
% USAGE:
%
%
% INPUTS:
%       refractive            - disretised refractive index
%       refractive_gradient_x - discretised refractive index gradient along x
%       refractive_gradient_y - discretised refractive index gradient along y
%       refractive_gradient_z - discretised refractive index gradient along z
%       cartesian_position_emitter - a dim x 1 Cartesian position of the emitter
%       direction        - a dim x 1 vector of the cartesian initial direction along
%                         the ray
%       xvec             - the x vector of grid points
%       yvec             - the y vector of grid points
%       zvec             - the z vector of grid points
%       pos_grid_first   - a dim x 1 Cartesian position of the first index of
%                         the grid
%       pos_grid_end     - a dim x 1 Cartesian position of the end index of
%                         the grid
%       dx               - a scalar representing the grid spacing, the same
%                         along all the Cartesian coordinates [m]
%       ds               - a saclar representing the ray spacing [m]
%       grid_size        - the size of the grid
%       dim              - the dimension of the medium
%       mask             - a binary mask used for calculating the
%                          coefficients.
%                          the refractive index in water can be assumed the
%                          same, and thus the associated line integarals are
%                          cancelled for a difference imaging.)
%       calc_coeffs           - a boolean controlling whether the
%                               interpolation coefficients are stored or
%                               not (This is often set true only for a linked
%                               ray after solving the ray linking problem.)
% OPTIONAL INPUTS:

%
% OUTPUTS:
%
%      pos                - a dim x 1 vector of the cartesian poistion of the
%                           end point along the ray
%      interp_coeff_vec - the vector of interpolation coefficients for the points
%                         along the ray within the mask

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian and Ben Cox

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
    pos_grid_first, dx, grid_size, dim, mask, refractive,...
    raytogrid_indices_x, raytogrid_indices_y,...
    raytogrid_indices_z, raytogrid_coeff_matrix, ...
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




% update the direction
direction = direction + ds * dn;




% store the current position
pos_previous = pos;
% update the current position of the ray
pos = pos_previous + ds/n * direction;


% the ray must be terminated inside the detection, and also the
% pointsincluded in the interpolation must not exceed the compuational grid
while   norm(pos) - detec_radius <= 1e-10    &&   all (pos > pos_grid_first + dx  &  pos < pos_grid_end - 2 * dx )
    
    % calculate the interpolation coefficients and their corresponding indices
    % on the grid
    [indices, coeff, n, dn] = interpLocalBspline(pos, xvec, yvec, zvec, ...
        pos_grid_first, dx, grid_size, dim, mask, refractive,...
        raytogrid_indices_x, raytogrid_indices_y,...
        raytogrid_indices_z, raytogrid_coeff_matrix, ...
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
    
    
    
    
    % update the direction
    direction = direction + ds * dn;
    
    % store the current position
    pos_previous = pos;
    % update the current position of the ray
    pos = pos_previous + ds/n * direction;
    
    
end


%  the last point must be on the detection surface
[pos , ds_final] = calcIntersectBall(pos_previous, 1/n * direction, detec_radius);


if calc_coeffs
    
    % multiply by the ray spacing to get the integral of the refractive index
    % (acoustic length) along the ray
    interp_coeff_vec = ds * interp_coeff_vec;
    
    
    
    % correct the integration for the last point
    if ~isempty(coeff)
        interp_coeff_vec (indices) = interp_coeff_vec (indices) + ...
            1/2 * (ds_final - ds) * coeff';
    end
    
   if all (pos > pos_grid_first + dx  &  pos < pos_grid_end - 2 * dx )
        
            % calculate the interpolation coefficients and their corresponding indices
            % on the grid
            [indices, coeff, ~, ~] = interpLocalBspline(pos, xvec, yvec, zvec,...
                pos_grid_first, dx, grid_size, dim, mask, refractive,...
                raytogrid_indices_x, raytogrid_indices_y,...
                raytogrid_indices_z, raytogrid_coeff_matrix, ...
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