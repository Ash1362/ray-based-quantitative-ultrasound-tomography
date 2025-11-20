function [recon_grid] = makeReconstructionGrid(grid_spacing, maximum_position, varargin)
%MAKERECONSTRUCTIONGRID makes a grid for image reconstruction

% DESCRIPTION:
%      makeReconstructionGrid makes a grid for image reconstruction 
%      in ultrasound computed tomography, when ray-based approaches
%      are used for image reconstruction
% USAGE:
%      
%
% INPUTS:
%       grid_spacing      - the grid spacing along coordinates x, y and z
%       maximum_position  - a scalar representing an esimated maximum position of 
%                           the transducers from the origin of Cartesian
%                           coordinate, which is the centre of the
%                           detection volume [m]. For a circular or
%                           hemispherical detection geometry, this will be
%                           radius of the detection circle (surface).
%
% OPTIONAL INPUTS:
%       'grid_expansion'     - a 1 x dim vector for extension of the grid
%                              beyond the detection surface
%       'reconstruction_geometry' - the detection surface that can be 
%                              'real' (like Pammoth) or
%                              'origin-centred' (like kWave)
%       'z_pos_height'       - the maximum z position of the grid points, if
%                              the reconstruction geometry is set 'real'
%       
%
% OUTPUTS:
%       recon_grid        - reconstruction grid in the form of a struct array
%                           with fields:
%       'x_vec'           - the vector of position of the grid points along
%                           x direction
%       'y_vec'           - the vector of position of the grid points along
%                           y direction
%       'z_vec'           - the vector of position of the grid points along
%                           z direction
%       'Nx'              - the number of grid points along x direction
%       'Ny'              - the number of grid points along y direction
%       'Nz'              - the number of grid points along z direction
%       'x'               - dim-D position of the grid points along x
%       'y'               - dim-D position of the grid points along y
%       'z'               - dim-D position of the grid points along z
%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 


para.grid_expansion = [1.02, 1.02, 1.02];
para.reconstruction_geometry = 'real';
para.z_pos_height = 5e-3;


% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


% get the number of the dimensions 
recon_grid.dim = size(grid_spacing, 2);

% grid spacing
recon_grid.dx = grid_spacing(1);
recon_grid.dy = grid_spacing(2);
if recon_grid.dim == 3
    recon_grid.dz = grid_spacing(3);
end
% a vector of position of the grid points along Cartesian coordinate x
recon_grid.x_vec = (- para.grid_expansion(1)*maximum_position: grid_spacing(1) : + para.grid_expansion(1)*maximum_position).';
% the number of the grid points along coordinate x
recon_grid.Nx = length(recon_grid.x_vec);
if recon_grid.dim > 1
    % a vector of position of the grid points along Cartesian coordinate y
    recon_grid.y_vec = (- para.grid_expansion(2)*maximum_position: grid_spacing(2) : + para.grid_expansion(2)*maximum_position).';
    % the number of the grid points along coordinate x
    recon_grid.Ny = length(recon_grid.y_vec);
    
    if recon_grid.dim > 2
        switch para.reconstruction_geometry
            case 'real'
                
                % a vector of position of the grid points along coordinate z
                recon_grid.z_vec = (- para.grid_expansion(3)*maximum_position:...
                    grid_spacing(3): para.z_pos_height).';
                % the number of the grid points along coordinate x
                recon_grid.Nz = length(recon_grid.z_vec);
                
            case 'origin-centred'
                % a vector of position of the grid points along coordinate z
                recon_grid.z_vec = (- para.grid_expansion(3)*maxium_position:...
                    grid_spacing(3): + para.grid_expansion(3)*maximum_position).';
                % the number of the grid points along coordinate x
                recon_grid.Nz = length(recon_grid.z_vec);
            otherwise
        end
    end
end



% dim-D grid coordinates based on the coordinates contained in the
% positions vectors
switch recon_grid.dim
    case 1
        recon_grid.x = x_vec;
    case 2
        [recon_grid.x, recon_grid.y] = ndgrid(recon_grid.x_vec, recon_grid.y_vec);
    case 3
        [recon_grid.x, recon_grid.y, recon_grid.z] = ndgrid(recon_grid.x_vec, recon_grid.y_vec, recon_grid.z_vec);
end


% the size of the grid
switch recon_grid.dim
    case 2
        recon_grid.size = [recon_grid.Nx, recon_grid.Ny];
    case 3
        recon_grid.size = [recon_grid.Nx, recon_grid.Ny, recon_grid.Nz];
end
        


end