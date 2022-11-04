function [pos, interp_coeff_vec, acoustic_length_ray, num_points] = testSecondOrderRay(refractive,...
    params, cartesian_position_emitter, direction, ray_grid, pos_grid_first, ds, varargin)
%TESTSECONDORDERRAY measures the accuracy of ray tracing using the second
%order approach
%
% DESCRIPTION:
%
% USAGE:
% testSecondOrderRay measures the accuracy of ray tracing using the second
% order method on a 'Maxwell Fish-eye' phantom. the measures are 'ray-path'
% (the position of the points  along the ray and their) and 'acoustic length'
% (acoustic length). The second-order method includes 'Dual-update' (proposed approach)
% and 'Mixed-step'(prototype approach), and the interpolation can be done using
% 'Bilinear' or 'Bspline' methods
%
% INPUTS:
%       refractive            - disretised refractive index
%       params                - a struct containing the fields:
%                               if interpolation method is 'Bilinear',
%       'x' - discretised refractive index gradient along x
%       'y' - discretised refractive index gradient along y
%       'z' - discretised refractive index gradient along z
%                               if interpolation method is 'Bspline-Catmull-Rom'
%                               or 'Bspline-Denis',
%       'x' -  x indices for B-spline interpolation
%       'y' -  y indices for B-spline interpolation                                  
%       'z' -  z indices for B-spline interpolation
%       'coeff_matrix'            - matrix for calculating B-spline
%                                  interpolation coefficients of the field
%       'coeff_derivative_matrix' - matrix for calculating B-spline
%                                  interpolation coefficients of the
%                                  directional gradients of the field 
%          
%       cartesian_position_emitter - a dim x 1 Cartesian position of the emitter
%       direction        - a dim x 1 vector of the cartesian initial direction along
%                         the ray
%       ray_grid         - a struct for defining the computational grid
%       pos_grid_first   - a dim x 1 Cartesian position of the first index of
%                         the grid
%       ds               - a saclar representing the spacing along the ray [m]
%       calc_coeffs     - a boolean controlling whether the
%                         interpolation coefficients are stored or
%                         not (This is often set true only for a linked
%                         ray after solving the ray linking problem.)
% OPTIONAL INPUTS:

%
% OUTPUTS:
%
%      pos              - a dim x 1 vector of the cartesian poistion of the
%                         end point along the ray
%      interp_coeff_vec - the vector of interpolation coefficients for the points
%                         along the ray within the mask
%      acoustic_length_ray - the acoustic length along the ray for a full
%                            rotation
%      num_points          - the number of points along the ray

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2020 Ashkan Javaherian 



para = [];
para.metric = 'acoustic-length';
para.ray_method = 'mixed-step';
para.raytogrid_interp = 'Bspline-Denis';

if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


% get the parameters of the computational grid
dim = ray_grid.dim;
xvec = ray_grid.x_vec;
yvec = ray_grid.y_vec;
zvec = ray_grid.z_vec;
dx = ray_grid.dx;


switch dim 
    case 2   
        grid_size = [ray_grid.Nx, ray_grid.Ny];
    case 3    
        grid_size = [ray_grid.Nx, ray_grid.Ny, ray_grid.Nz];
end

% define the mask true for the whole computational grid
mask = true(prod(grid_size), 1);


% get the parameters
switch para.raytogrid_interp
    case 'Bilinear'
        refractive_gradient_x = params.x;
        refractive_gradient_y = params.y;
        refractive_gradient_z = params.z;
    case {'Bspline-Catmull-Rom', 'Bspline-Denis'}
        raytogrid_indices_x = params.x;
        raytogrid_indices_y = params.y;
        raytogrid_indices_z = params.z;
        raytogrid_coeff_matrix = params.raytogrid_coeff_matrix;
        raytogrid_coeff_derivative_matrix = params.raytogrid_coeff_derivative_matrix;
end


if strcmp(para.metric, 'acoustic-length')
    
    % allocate a sparse vector for storing the interpolation coefficients
    interp_coeff_vec = sparse(prod(grid_size), 1);
end
% remove the last index from the size of the grid (not used)
grid_size = grid_size(1: dim-1);



% normalise direction - the direction along the ray must be a unit vector
direction = 1/norm(direction) * direction;

% initialise the counter for the number of ray points
num_points = 1;
% initialise the position of the ray with the emission point
pos(:, num_points) = cartesian_position_emitter;

% calculate the interpolation coefficients and their corresponding indices
% on the grid
switch para.raytogrid_interp
    case 'Bilinear'
        [indices, coeff, n, dn] = interpLocal(pos(:, num_points), xvec, yvec, zvec, ...
            pos_grid_first, dx, grid_size, dim, mask, refractive,...
            refractive_gradient_x, refractive_gradient_y, refractive_gradient_z);
    case {'Bspline-Catmull-Rom', 'Bspline-Denis'}
        [indices, coeff, n, dn] = interpLocalBspline(pos(:, num_points), xvec, yvec, zvec, ...
            pos_grid_first, dx, grid_size, dim, mask, refractive,...
            raytogrid_indices_x, raytogrid_indices_y, raytogrid_indices_z,...
            raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix, true);
end





% calculate the secon-order term
h =  ds/n * (dn - (dn' * direction) * direction);

% integrate the refractive index along the ray to calculate the
% acoustic length, if requested
if strcmp(para.metric, 'acoustic-length')
    interp_coeff_vec(indices) = interp_coeff_vec(indices) + (1/2) * coeff';
end


switch para.ray_method
    case 'dual-update'
        
        % update the direction
        direction_ray = direction + h/2;
        % normalise the direction along the ray
        direction_ray = 1/norm(direction_ray) * direction_ray;
        % update the current position of the ray
        pos(:, num_points+1) = pos(:, num_points)+ ds * direction_ray;
        % update the direction for updating the secon-order derivative
        direction = direction + h;
        % normalise the direction for updating the second-order derivative
        direction = direction /norm(direction);
        
        
    case 'mixed-step'
        
        % update the direction
        direction = direction + h/2;
        % normalise the direction along the ray, which is the same as the direction for
        % updating the second-order derivative
        direction = 1/norm(direction) * direction;
        % update the current position of the ray
        pos(:, num_points+1) = pos(:, num_points)+ ds * direction;
        
end

% increase the counter for the number of ray points
num_points = num_points + 1;

% terminte the ray before interception with the emission point after a full
% rotation
while   norm(pos(:, num_points) - cartesian_position_emitter)  > ds    ||  num_points < 3
  
    % calculate the interpolation coefficients and their corresponding indices
    % on the grid
    switch para.raytogrid_interp
        case 'Bilinear'
            [indices, coeff, n, dn] = interpLocal(pos(:, num_points), xvec, yvec, zvec, ...
                pos_grid_first, dx, grid_size, dim, mask, refractive,...
                refractive_gradient_x, refractive_gradient_y, refractive_gradient_z);
        case {'Bspline-Catmull-Rom', 'Bspline-Denis'}
            [indices, coeff, n, dn] = interpLocalBspline(pos(:, num_points), xvec, yvec, zvec, ...
                pos_grid_first, dx, grid_size, dim, mask, refractive,...
                raytogrid_indices_x, raytogrid_indices_y, raytogrid_indices_z,...
                raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix, true);
    end
    
    
    
    
    
    
    % calculate the secon-order term
    h =  ds/n * (dn - (dn' * direction) * direction);
    % integrate the refractive index along the ray to calculate the
    % acoustic length, if requested
    if  strcmp(para.metric, 'acoustic-length')
        interp_coeff_vec(indices) = interp_coeff_vec(indices) + coeff';
    end
    
    switch para.ray_method
        case 'dual-update'
            
            % update the direction
            direction_ray = direction + h/2;
            % normalise the direction along the ray
            direction_ray = 1/norm(direction_ray) * direction_ray;
            % update the current position of the ray
            pos(:, num_points+1) = pos(:, num_points)+ ds * direction_ray;
            % update the direction for updating the secon-order derivative
            direction = direction + h;
            % normalise the direction for updating the second-order derivative
            direction = direction /norm(direction);
            
            
        case 'mixed-step'
            
            % update the direction
            direction = direction + h;
            % normalise the direction along the ray, which is the same as the direction for
            % updating the second-order derivative
            direction = 1/norm(direction) * direction;
            % update the current position of the ray
            pos(:, num_points+1) = pos(:, num_points)+ ds * direction;
            
    end
    
    % increase the counter for the number of ray point
    num_points = num_points + 1;
    
    
end

switch para.metric
    
    case 'acoustic-length'
        
        %  match the position of the last point  point must be on the detection surface
        ds_final = norm(cartesian_position_emitter - pos(:, end));
        
        
        
        % multiply by the ray spacing to get the integral of the refractive index
        % (acoustic length) along the ray
        interp_coeff_vec = ds * interp_coeff_vec;
        
        
        
        % correct the integration for the last point
        interp_coeff_vec (indices) = interp_coeff_vec (indices) + ...
            1/2 * (ds_final - ds) * coeff';
        
        
        
        % calculate the interpolation coefficients and their corresponding indices
        % on the grid for the last point, which is the position of the
        % emission point
        switch para.raytogrid_interp
            case 'Bilinear'
                [indices, coeff, ~, ~] = interpLocal(cartesian_position_emitter, xvec, yvec, zvec, ...
                    pos_grid_first, dx, grid_size, dim, mask, refractive);
            case {'Bspline-Catmull-Rom', 'Bspline-Denis'}
                [indices, coeff, ~, ~] = interpLocalBspline(cartesian_position_emitter, xvec, yvec, zvec, ...
                    pos_grid_first, dx, grid_size, dim, mask, refractive,...
                    raytogrid_indices_x, raytogrid_indices_y, raytogrid_indices_z,...
                    raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix, false);
        end
        
        interp_coeff_vec(indices) = interp_coeff_vec(indices) +...
            (1/2) * ds_final * coeff';
        
        
        acoustic_length_ray = interp_coeff_vec' * refractive(:);
        
    case 'ray-path'
        
        % give an empty variable, if the interpolation coefficients
        % are not required
        interp_coeff_vec = [];
        acoustic_length_ray = nan;
        
end

        
end