function [indices, index_x, index_y, index_z] = calcInterpIndicesBspline(...
    pos_point, pos_grid_first, dx, grid_size, dim, raytogrid_indices_x,...
    raytogrid_indices_y, raytogrid_indices_z)
%CALCULATEINTERPINDICES calculates the indices of the voxel at which
%a target point is placed


% calculateInterpIndices calculates the indices of the voxel at which a
% target point exist

% USAGE:
%
%
% INPUTS:
%       pos_point       - a dim x 1 vector of the Cartesian position of the
%                         target point
%       pos_grid_first  - a dim x 1 Cartesian position of the first index of
%                         the grid
%       grid_spacing    - a scalar representing the grid spacing, the same
%                         along all the Cartesian coordinates
%       grid_size       - the size of the grid along x (2d case) and x and
%                         y (3D case)
%       dim             - the dimension of the medium
%      raytogrid_indices_x -  x indices for B-spline interpolation
%      raytogrid_indices_y -  y indices for B-spline interpolation
%      raytogrid_indices_z -  z indices for B-spline interpolation
% OPTIONAL INPUTS:

%
% OUTPUTS:
%
%       indices         - the indices of the grid points forming the voxel
%                         in a Matlab column-wise format
%       index_x         - the x index of the grid point that is the first index
%                         of the voxel
%       index_y         - the y index of the grid point that is the first index
%                         of the voxel
%       index_z         - the z index of the grid point that is the first index
%                         of the voxel

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian 

index = round(1 + floor((pos_point - pos_grid_first)/dx));
index_x = index(1);
index_y = index(2);
indices_x = index_x + raytogrid_indices_x;
indices_y = index_y + raytogrid_indices_y;    


switch dim
    case 2
        index_z = [];
       
        % A fast implementation of sub2ind
        indices  = indices_x + (indices_y-1) * grid_size(1);
        
        
        
    case 3
        index_z = index(3);
        
        indices_z = index_z + raytogrid_indices_z;
        % A fast implementation of sub2ind
        indices = indices_x + (indices_y-1) * grid_size(1) + (indices_z-1) * prod(grid_size);
        
end


end