function [binary_mask, interp_matrix] = interpNeighborUST(comp_grid,...
    transducer_position)
%INTERPNEIGHBORUST CONSTRUCTS A SPARSE MATRIX FOR NEIGHBORING INTERPOLATION
%
% DESCRIPTION:
%       interpNeighborUST constructs a sparse matrix for interpolation
%       between the transducers and the grid points
%
% USAGE:
%      

%
% INPUTS:
%       comp_grid   - the computational grid
%       transducer_position - the dim x N_t Cartesian position of the
%                             transducers, which are assumed points

%
% OPTIONAL INPUTS:    
%      
%      
%
% OUTPUTS:
%      binary_mask    - binary mask of size of the grid, and is true for all
%                       the uniion of all grid points contributed into the
%                       neigboring interpolation
%      interp_matrix -  a sparse matrix for neigboring interpolation between
%                       the grid points and the transducers, which are
%                       assumed points
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 15.12.2019
%       last update     - 15.12.2019
%       
%
% 
% This function is part of the r-Wave Toolbox.
% Copyright (C) 2022 Ashkan Javaherian 
%%


% get the number of dimensions
dim = comp_grid.dim;

% get the number of grid points along te Cartesian coordinates
Nx  = comp_grid.Nx;
Ny  = comp_grid.Ny;

switch dim
    case 2
N    =  Nx*Ny;
Nxyz = [Nx,Ny];          
    case 3
Nz = comp_grid.Nz;
N    = Nx*Ny*Nz;
Nxyz = [Nx,Ny,Nz];
end

% get the number of transducers
num_transducer = size(transducer_position, 2);


% get the linear index of grid points
grid_point_indices = reshape(1:N, Nxyz);
switch dim
    case 2
        % allocate a false binary matrix for the binary mask for
        % interpolation coefficients for all transducers
        binary_mask = false(Nx, Ny);
        
        % get the d-dimensional Cartesian coordinates of the grid points 
        [x, y] = ndgrid(kgrid.x_vec, kgrid.y_vec);
        
        % get the handle function for interpolation to each grid point
        interp_map = griddedInterpolant(x , y , grid_point_indices, 'nearest');
    case 3
        % allocate a false binary matrix for the binary mask for
        % interpolation coefficients for all transducers
        binary_mask = false( Nx, Ny, Nz);
        
        % get the d-dimensional Cartesian coordinates of the grid points 
        [x, y, z] = ndgrid(kgrid.x_vec, kgrid.y_vec, kgrid.z_vec);
        
        % get the handle function for interpolation to each grid point
        interp_map = griddedInterpolant(x, y, z, grid_point_indices, 'nearest');
end


% fill the binary mask
indices = interp_map(transducer_position.');
binary_mask(indices) = true;

% get the true indices of the binary mask
binary_mask_indices  = find(binary_mask(:));
                
% get the location of each index in the true binary indices
[~, loc] = ismember(indices, binary_mask_indices);

% get the sparse matrix for neighboring interpolation
interp_matrix = sparse(1:num_transducer, loc, ones(1, num_transducer), num_transducer, ...
                            nnz(binary_mask));
end           