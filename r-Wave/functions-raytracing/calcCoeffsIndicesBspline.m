function [indices, coeff, coeff_derivative_x, coeff_derivative_y, coeff_derivative_z] ...
    = calcCoeffsIndicesBspline(pos_point, xvec, yvec, zvec,...
    pos_grid_first, grid_size, dx,  dim, mask, raytogrid_indices_x,...
    raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix, ...
    raytogrid_coeff_derivative_matrix)
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
%      raytogrid_indices_x -  x indices for B-spline interpolation
%      raytogrid_indices_y -  y indices for B-spline interpolation
%      raytogrid_indices_z -  z indices for B-spline interpolation
%      raytogrid_coeff_matrix - matrix for calculating B-spline
%                                   interpolation coefficients of the field
%      raytogrid_coeff_derivative_matrix - matrix for calculating B-spline
%                                         interpolation coefficients of the
%                                         directional gradients of the field
% OPTIONAL INPUTS:

%
% OUTPUTS:
%
%       indices         - the indices of the voxel at which the target point resides
%       coeff           - the interpolation coefficient of the indices used
%                         for the interpolationof the refractive index
%       coeff_derivative_x - the interpolation coefficient of the indices used
%                            for the interpolation of the directional derivative
%                            of the refractive index along x
%       coeff_derivative_y - the interpolation coefficient of the indices used
%                            for the interpolation of the directional derivative
%                            of the refractive index along y
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian 


[indices, index_x, index_y, index_z] = calcInterpIndicesBspline(pos_point,...
    pos_grid_first, dx, grid_size, dim, raytogrid_indices_x,...
    raytogrid_indices_y, raytogrid_indices_z);



% Calculate the coefficients if at least one grid point from the grid cell 
% encompassing the ray point are inside the binary mask
if all(indices - numel(mask)<0)
if any(mask(indices))
    
    % get the distance of the off-grid point and the first grid point
    % divided by the grid spacing along 
    % x coordinate
    xd = (pos_point(1) - xvec(index_x))/dx;
    % y coordinate
    yd = (pos_point(2) - yvec(index_y))/dx;
    
    % get the polynomial of the interpolation of the map along
    % x coordinate
    polynomial_x = raytogrid_coeff_matrix * [xd^3; xd^2; xd; 1];
    % y coordinate
    polynomial_y = raytogrid_coeff_matrix * [yd^3; yd^2; yd; 1];
    
    % get the polynomial of the interpolation of the map along
    % x coordinate
    polynomial_derivative_x = raytogrid_coeff_derivative_matrix *[xd^2; xd; 1];
    % y coordinate
    polynomial_derivative_y = raytogrid_coeff_derivative_matrix *[yd^2; yd; 1];
    
    
    
    switch dim
        
        case 2
            
            coeff = (vectorise(polynomial_x * polynomial_y'))';
            coeff_derivative_x = (vectorise(polynomial_derivative_x * polynomial_y'))';
            coeff_derivative_y = (vectorise(polynomial_x * polynomial_derivative_y'))';
            coeff_derivative_z = [];
            
        case 3
            
            zd = (pos_point(3) - zvec(index_z))/dx;
            polynomial_z = raytogrid_coeff_matrix * [zd^3; zd^2; zd; 1];
            polynomial_derivative_z = raytogrid_coeff_derivative_matrix *[zd^2; zd; 1];
            
            
            
            % get the interpolation coefficients for interpolation of the
            % map
            coeff = vectorise(polynomial_x * polynomial_y'.*...
                reshape(polynomial_z, 1, 1, [])).';
            
            % get the interpolation coefficients for interpolation of the
            % the first derivative along
            % x coordinate
            coeff_derivative_x = vectorise(polynomial_derivative_x * polynomial_y'.*...
                reshape(polynomial_z, 1, 1, [])).';
          
            % y coordinate
            coeff_derivative_y = vectorise(polynomial_x * polynomial_derivative_y'.*...
                reshape(polynomial_z, 1, 1, [])).';
            
            % z coordinate
            coeff_derivative_z = vectorise(polynomial_x * polynomial_y'.*...
                reshape(polynomial_derivative_z, 1, 1, [])).';
            
            
            %p_xy = polynomial_x * polynomial_y';
            %p_dxy = polynomial_derivative_x * polynomial_y';
            %p_xdy = polynomial_x * polynomial_derivative_y';
            
            %coeff = zeros(4,4,4);
            %coeff_derivative_x = zeros(4,4,4);
            %coeff_derivative_y = zeros(4,4,4);
            %coeff_derivative_z = zeros(4,4,4);
            
            %for i = 1:4
            %    coeff(:,:,i) = polynomial_z(i) * p_xy;
            %    coeff_derivative_x(:,:,i) = polynomial_z(i) * p_dxy;
            %   coeff_derivative_y(:,:,i) = polynomial_z(i) * p_xdy;
            %    coeff_derivative_z(:,:,i) = polynomial_derivative_z(i) * p_xy;
            %end
            
            %coeff = (coeff(:))';
            %coeff_derivative_x = (coeff_derivative_x(:))';
            %coeff_derivative_y = (coeff_derivative_y(:))';
            %coeff_derivative_z = (coeff_derivative_z(:))';
    end
    
    
else
    
    coeff = [];
    coeff_derivative_x = [];
    coeff_derivative_y = [];
    coeff_derivative_z = [];
    
end

else
    
    coeff = [];
    coeff_derivative_x = [];
    coeff_derivative_y = [];
    coeff_derivative_z = []; 
    
end



end
