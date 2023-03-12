function [indices, coeff, refractive_target, refractive_gradient_target] = interpLocal(pos_point, ...
    xvec, yvec, zvec, pos_grid_first, dx, grid_size, dim, mask, refractive,...
    refractive_gradient_x, refractive_gradient_y, refractive_gradient_z)
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
%       refractive            - disretised refractive index
%       refractive_gradient_x - the discretised refractive index gradient along x
%       refractive_gradient_y - the discretised refractive index gradient along y
%       refractive_gradient_z - the discretised refractive index gradient along z
% OPTIONAL INPUTS:

%
% OUTPUTS:
%
%       indices         - the indices of the voxel at which the target point resides
%       coeff           - the interpolation coefficient of the indices of the voxel
%                         for the target point
%       refractive      - the interpolated refractive index on the
%                         target point
%       refractive_gradient_indices - the dim * 1 vector of the interpolated
%                                     directional gradients of the refractive index
%                                     on the target point

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian



[indices, coeff] = calcCoeffsIndices(pos_point, xvec, yvec,...
        zvec, pos_grid_first, grid_size, dx, dim, mask);




if isempty(coeff)
    
    % avoid interpolation, if the coeeficients are not collected
    refractive_target = [];
    refractive_gradient_target = [];
    
else
    
    % interpolate from the discretised grid to the target point
    refractive_target = coeff * refractive(indices);
    
    if nargin > 10    %  && any(mask(indices))
        
        if dim == 2
            refractive_gradient_target = [coeff * refractive_gradient_x(indices);...
                coeff * refractive_gradient_y(indices)];
        else
            refractive_gradient_target = [coeff * refractive_gradient_x(indices);
                coeff * refractive_gradient_y(indices);
                coeff*  refractive_gradient_z(indices) ];
        end
        
    else
        
        % avoid calculating the directional gradient on the ray point
        refractive_gradient_target = [];
        
    end
    
    
    
end




end