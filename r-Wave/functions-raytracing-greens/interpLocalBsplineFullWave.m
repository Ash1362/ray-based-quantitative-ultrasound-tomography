function [indices, coeff, refractive_target, refractive_gradient_target,...
    refractive_second_gradient_target, refractive_nonsmoothed_target, absorption_coeff_target] = ...
    interpLocalBsplineFullWave(pos_point, xvec, yvec, zvec, pos_grid_first,...
    dx, grid_size, dim, mask, refractive, refractive_nonsmoothed, absorption_coeff,...
    raytogrid_indices_x, raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix,...
    raytogrid_coeff_derivative_matrix, raytogrid_coeff_second_derivative_matrix,...
    calc_gradient_coeffs, paraxial)
%INTERPLOCALBSPLINEFULLWAVE gives an approximation of a discretised field at
% a target point using B-spline interpolation.
%
% DESCRIPTION:
% interpLocalBsplineFullWave provides an approximation of a field and its
% directional gradients at a given off-grid point.
%
% USAGE:
%
%
% INPUTS:
%       pos_point       - a dim x 1 vector of the Cartesian position [m] of the
%                         target point
%       xvec            - the x vector of grid points
%       yvec            - the y vector of grid points
%       zvec            - the z vector of grid points
%       pos_grid_first  - a dim x 1 Cartesian position of the first index of
%                         the grid
%       dx              - a scalar representing the grid spacing, the same
%                         along all the Cartesian coordinates
%       grid_size       - the 1 x dim grid size along all dimensiosn
%       dim             - the dimensions of the grid
%       mask            - a binary mask used for ray tracing
%       refractive      - disretised smoothed refractive index for computing the
%                        trajectory of ray
%       refractive_nonsmoothed - the disretised nonsmoothed refractive index,
%                                which is integrated along the rays
%       absorption_coeff   - the discretised absorption coefficient, which
%                            is integrated along the rays.
%      raytogrid_indices_x -  x indices for B-spline interpolation
%      raytogrid_indices_y -  y indices for B-spline interpolation
%      raytogrid_indices_z -  z indices for B-spline interpolation
%      raytogrid_coeff_matrix - matrix for calculating B-spline interpolation
%                               coefficients of the field
%      raytogrid_coeff_derivative_matrix - matrix for calculating B-spline
%                                         interpolation coefficients of the
%                                         directional gradients of the field
%      calc_gradient_coeffs   - Boolean controlling wheter the interpolation
%                               coefficients for directional gradients are
%                               calculated.
%      paraxial               - Boolean indicating whether the traced ray is
%                               paraxial or not.

% OPTIONAL INPUTS:

%
% OUTPUTS:
%
%       indices         - the indices of the voxel at which the target point resides
%       coeff           - the interpolation coefficient of the indices of the voxel
%                         for the target point
%       refractive_target - a scalar representing the interpolated value of the
%                         smoothed refractive index from the grid points to 
%                         the current point on the ray
%       refractive_gradient_target - the dim x 1 vector of the interpolated
%                                    directional gradients of the refractive index
%                                    on the target point
%       refractive_second_gradient_target - the dim x dim matrix of the interpolated
%                                           directional second-order derivatives of the
%                                           refractive index on the target point
%       refractive_nonsmoothed_target - a scalar representing the interpolated value of the
%                         nonsmoothed refractive index from the grid points to 
%                         the current point on the ray 
%       absorption_coeff_target - a scalar representing the interpolated value of the
%                         nonsmoothed absorption coefficient from the grid points to 
%                         the current point on the ray 
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian



[indices, coeff, coeff_derivative_x, coeff_derivative_y, coeff_derivative_z,...
    coeff_second_derivative_xx, coeff_second_derivative_yy,...
    coeff_second_derivative_xy, coeff_second_derivative_xz, coeff_second_derivative_yz,...
    coeff_second_derivative_zz] =...
    calcCoeffsIndicesBspline1st2ndOrder(pos_point, xvec, yvec, zvec,...
    pos_grid_first, grid_size, dx, dim, mask, raytogrid_indices_x,...
    raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix, ...
    raytogrid_coeff_derivative_matrix, raytogrid_coeff_second_derivative_matrix,...
    paraxial);

if ~isempty(coeff)
    % interpolate from the discretised grid to the target point
    refractive_target = coeff * refractive(indices);
    % interpolate the nonsmoothed refractive index distribution from the
    % discretised grid to the target point
    refractive_nonsmoothed_target = coeff * refractive_nonsmoothed(indices);
    
    if isscalar(absorption_coeff)
        absorption_coeff_target = 0;
    else
        absorption_coeff_target = coeff * absorption_coeff(indices);
    end
    
    if calc_gradient_coeffs
        
        % get the refractive index on the grid points included in the
        % interpolation
        refractive_indices = refractive(indices);
        
        switch dim
            
            case 2
                
                % get the the vector for the first-order derivative of the
                % refractive index
                refractive_gradient_target = [coeff_derivative_x * refractive_indices;...
                    coeff_derivative_y * refractive_indices];
                
                if paraxial
                    
                % get the the matrix for the second-order derivative of the
                % refractive index
                refractive_second_gradient_xy =  coeff_second_derivative_xy * refractive_indices;
                refractive_second_gradient_target = [coeff_second_derivative_xx * refractive_indices,...
                    refractive_second_gradient_xy;...
                    refractive_second_gradient_xy, coeff_second_derivative_yy * refractive_indices];
                end
                
                
            case 3
                
                
                % get the the vector for the first-order derivative of the
                % refractive index
                refractive_gradient_target = [coeff_derivative_x * refractive_indices;
                    coeff_derivative_y * refractive_indices;
                    coeff_derivative_z * refractive_indices];
                
                
                % get the the matrix for the second-order derivative of the
                % refractive index
                if paraxial
                % get the off-diagonal components of the second-order
                % derivative matrix
                refractive_second_gradient_xy =  coeff_second_derivative_xy * refractive_indices;
                refractive_second_gradient_xz =  coeff_second_derivative_xz * refractive_indices;
                refractive_second_gradient_yz =  coeff_second_derivative_yz * refractive_indices;
                
                % get the off-diagonal components of the second-order
                % derivative matrix
                refractive_second_gradient_target = [coeff_second_derivative_xx * refractive_indices,...
                    refractive_second_gradient_xy, refractive_second_gradient_xz;...
                    refractive_second_gradient_xy, coeff_second_derivative_yy * refractive_indices,...
                    refractive_second_gradient_yz; refractive_second_gradient_xz,...
                    refractive_second_gradient_yz, coeff_second_derivative_zz * refractive_indices];
                end
                
        end
        
        
        
    else
        
        refractive_gradient_target = [];
        refractive_second_gradient_target = [];
    end
else
    
    refractive_target = [];
    refractive_gradient_target = [];
    refractive_second_gradient_target = [];
    refractive_nonsmoothed_target = [];
    absorption_coeff_target = [];
    
end

% if the traced ray is not paraxial ray, the second-order gradient is not
% required
if ~paraxial
    refractive_second_gradient_target = [];
end

end

