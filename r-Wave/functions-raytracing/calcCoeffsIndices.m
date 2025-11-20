function [indices, coeff] = calcCoeffsIndices(pos_point, xvec, yvec, zvec, pos_grid_first, grid_size, dx,  dim, mask)
%CALCULATECOEFFSINDICES calculates the indices and their corresponding
%interpolation coefficients
%
%
% calcCoeffsIndices calculates the the indices of the grid points and their
% corresponding interpolation coefficients for the voxel encompassing a target
% point. A local ''bilinear'' interpolation is used.
%
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
%       grid_spacing    - a scalar representing the grid spacing, the same 
%                         along all the Cartesian coordinates
%       grid_size       - the size of the grid along x (2d case) and x and
%                         y (3D case)
%       dim             - the dimension of the medium
%       mask            - a binary mask used for calculating the
%                         coefficients for difference imaging between
%                         phantom and water
% OPTIONAL INPUTS:

%        
% OUTPUTS:
%
%       indices         - the indices of the voxel at which the target point resides
%       coeff           - the interpolation coefficient of the indices of the voxel
%                         for the target point
%                                 
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 


[indices, index_x, index_y, index_z] = calcInterpIndices(pos_point, pos_grid_first, dx, grid_size, dim);

% Calculate the coefficients if at least one grid point from the grid cell 
% encompassing the ray point are inside the binary mask
if all(indices - numel(mask)<0)
if any(mask(indices))
    
    
    
    switch dim
        
        case 2
            
            
            xd   = (pos_point(1) - xvec(index_x))/dx;
            yd   = (pos_point(2) - yvec(index_y))/dx;
            xdyd =  xd * yd;
            
            
            coeff = [1 - xd - yd + xdyd, xd - xdyd, yd - xdyd, xdyd];
            
            
        case 3
            
            
            xd = (pos_point(1) - xvec(index_x))/dx;
            yd = (pos_point(2) - yvec(index_y))/dx;
            zd = (pos_point(3) - zvec(index_z))/dx;
            xdyd   = xd*yd;
            ydzd   = yd*zd;
            xdzd   = xd*zd;
            xdydzd = xdyd*zd;
            
            
            
            coeff = [1 - xd - yd - zd + xdyd + ydzd + xdzd - xdydzd, ...
                + zd - ydzd - xdzd + xdydzd, ...
                + yd - xdyd - ydzd + xdydzd, ...
                + ydzd - xdydzd, ...
                + xd - xdyd - xdzd + xdydzd, ...
                + xdzd - xdydzd, ...
                + xdyd - xdydzd, ...
                + xdydzd];
            
            
    end
    
    
else
    
    coeff = [];
    
end
else
   coeff = [];
end



end
