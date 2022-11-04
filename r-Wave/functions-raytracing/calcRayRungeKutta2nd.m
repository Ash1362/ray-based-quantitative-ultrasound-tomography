function [pos, interp_coeff_vec] = calcRayRungeKutta2nd(refractive, refrcative_gradient_x,...
    refractive_gradient_y,refractive_gradient_z, cartesian_position_emitter,...
    direction, xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds, grid_size,...
    dim, detec_radius, mask, calc_coeffs)
%CALCRAYRUNGEKUTTA2ND traces a ray using a secon-order Runge-Kutta algorithm,
% and stores the coefficients for the ray-to-grid interpolation. The interpolation 
% is done using a Biliear approach.
%
% DESCRIPTION:
% calRayRungeKutta2nd calculates positions along the ray and their interplation
% coefficients, given an initial position and an initial direction along the
% ray. The interpolation coefficients will be used for an integration of
% the acoustic length along the ray given a discretised refractive index
% distribution. The ray tracing is done using a second-order variant of the
% Runge-Kutta algorithm, also Known as Heun's method, and the interpolations
% between the grid and rays are done using a Bilinear approach.
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
%       detec_radius     - the radius of the detection ring (2D), or
%                          detection surface(3D)
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
%      pos               - a dim x 1 vector of the cartesian poistion of the 
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
[indices, coeff, n, dn] = interpLocal(pos, xvec, yvec, zvec, ...
    pos_grid_first, dx, grid_size, dim, mask, refractive,...
    refrcative_gradient_x, refractive_gradient_y, refractive_gradient_z);


if isempty(coeff)  
        n = 1; dn = 0;
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
[~, coeff_k2, n_k2, dn_k2] = interpLocal(pos_k2, xvec, yvec, zvec, ...
    pos_grid_first, dx, grid_size, dim, mask, refractive,...
    refrcative_gradient_x, refractive_gradient_y, refractive_gradient_z);


if isempty(coeff_k2)  
        n_k2 = 1; dn_k2 = 0;
end


% normalise direction - the direction along the ray must be a unit vector
% multiplied by refractive index
direction_k2 = n_k2/norm(direction_k2) * direction_k2;


k2_x = 1/n_k2 * direction_k2;


% store the current position
pos_previous = pos;

% update the current position of the ray
pos = pos_previous + ds/2 * (k1_x + k2_x);

% the ray must be terminated inside the detection
while   norm(pos) - detec_radius <= 1e-10    &&   all(pos >= pos_grid_first + dx & pos <= pos_grid_end - dx)
   
    
    
   % update the direction
    direction = direction + ds/2 * (dn + dn_k2);

% calculate the interpolation coefficients and their corresponding indices
% on the grid
[indices, coeff, n, dn] = interpLocal(pos, xvec, yvec, zvec, ...
    pos_grid_first, dx, grid_size, dim, mask, refractive,...
    refrcative_gradient_x, refractive_gradient_y, refractive_gradient_z);


if isempty(coeff)  
        n = 1; dn = 0;
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
[~, coeff_k2, n_k2, dn_k2] = interpLocal(pos_k2, xvec, yvec, zvec, ...
    pos_grid_first, dx, grid_size, dim, mask, refractive,...
    refrcative_gradient_x, refractive_gradient_y, refractive_gradient_z);


if isempty(coeff_k2)  
        n_k2 = 1; dn_k2 = 0;
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
[pos , ds_final] = calcIntersectBall(pos_previous, 1/2 * (k1_x + k2_x), detec_radius);


if calc_coeffs
    
    % multiply by the ray spacing to get the integral of the refractive index
    % (acoustic length) along the ray
    interp_coeff_vec = ds * interp_coeff_vec;
      
    % correct the integration for the last point
    if ~isempty(coeff)
        interp_coeff_vec (indices) = interp_coeff_vec (indices) + ...
            1/2 * (ds_final - ds) * coeff';
    end
    
    if all (pos >= pos_grid_first  &  pos <= pos_grid_end )
        
       
            % calculate the interpolation coefficients and their corresponding indices
            % on the grid
            [indices, coeff, ~, ~] = interpLocal(pos, xvec, yvec, zvec, ...
                pos_grid_first, dx, grid_size, dim, mask, refractive);
            
            
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