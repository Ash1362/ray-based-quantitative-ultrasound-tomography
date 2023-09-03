function [indices, coeff, refractive_target, refractive_gradient_target] = interpLocalBspline(pos_point, ...
    xvec, yvec, zvec, pos_grid_first, dx, grid_size, dim, mask, refractive,...
    raytogrid_indices_x, raytogrid_indices_y, raytogrid_indices_z,...
    raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix, calc_gradient_coeffs)
%INTERPLOCAL gives an approximation of a discretised field at a target point using
% bilinear interpolation.
%
% DESCRIPTION:
% interpLocal provides an approximation of a field and its directional gradients
% at a given point.
%
% USAGE:
%
%
% INPUTS:
%       pos_point       - a dim x 1 vector of the Cartesian position of the
%                         target point
%       xvec            - the x vector of grid points
%       yvec            - the y vector of grid points
%       zvec            - the z vector of grid points

%       pos_grid_first  - a dim x 1 Cartesian position of the first index of
%                         the grid
%       dx              - a scalar representing the grid spacing, the same
%                         along all the Cartesian coordinates
%       grid_size       - the size of the grid along x (2D case) and x and
%                         y (3D case)
%       dim             - the dimension of the medium
%       mask            - a binary mask used for calculating the
%                         coefficients for difference imaging between
%                         phantom and water
%      refractive            - disretised refractive index
%      raytogrid_indices_x -  x indices for B-spline interpolation
%      raytogrid_indices_y -  y indices for B-spline interpolation
%      raytogrid_indices_z -  z indices for B-spline interpolation
%      raytogrid_coeff_matrix     - matrix for calculating B-spline
%                                   interpolation coefficients of the field
%      raytogrid_coeff_derivative_matrix - matrix for calculating B-spline
%                                         interpolation coefficients of the
%                                         directional gradients of the field
%      calc_gradient_coeffs   - Boolean controlling wheter the interpolation
%                               coefficients for directional gradients are
%                               calculated.

% OPTIONAL INPUTS:

%
% OUTPUTS:
%
%       indices         - the indices of the voxel at which the target point resides
%       coeff           - the interpolation coefficient of the indices of the voxel
%                         for the target point
%       refractive_targe  - the interpolated refractive index on the
%                         target point
%       refractive_gradient_target - the dim * 1 vector of the directional
%                                    gradients of the refractive index
%                                    at the target point

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian 



[indices, coeff, coeff_derivative_x, coeff_derivative_y, coeff_derivative_z] =...
    calcCoeffsIndicesBspline(pos_point, xvec, yvec, zvec,...
    pos_grid_first, grid_size, dx, dim, mask, raytogrid_indices_x,...
    raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix, ...
    raytogrid_coeff_derivative_matrix);

   


if isempty(coeff)
    
    refractive_target = [];
    refractive_gradient_target = [];
    
else
    
    % interpolate from the discretised grid to the target point
    refractive_target = coeff * refractive(indices);
    
    
    if calc_gradient_coeffs
        
        switch dim
            
            case 2
                
                
                refractive_gradient_target = [coeff_derivative_x * refractive(indices);...
                    coeff_derivative_y * refractive(indices)];
                
                
            case 3
                
                
                refractive_gradient_target = [coeff_derivative_x * refractive(indices);
                    coeff_derivative_y * refractive(indices);
                    coeff_derivative_z * refractive(indices)];
                
                
        end
        
        
        
    else
        
        refractive_gradient_target = [];
        
    end
end

end

