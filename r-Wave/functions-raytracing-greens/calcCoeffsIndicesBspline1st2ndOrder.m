function [indices, coeff, coeff_derivative_x, coeff_derivative_y, coeff_derivative_z,...
    coeff_second_derivative_xx, coeff_second_derivative_yy,...
    coeff_second_derivative_xy, coeff_second_derivative_xz, coeff_second_derivative_yz,...
    coeff_second_derivative_zz] = calcCoeffsIndicesBspline1st2ndOrder(...
    pos_point, xvec, yvec, zvec,...
    pos_grid_first, grid_size, dx,  dim, mask, raytogrid_indices_x,...
    raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix, ...
    raytogrid_coeff_derivative_matrix, raytogrid_coeff_second_derivative_matrix,...
    paraxial)
%CALCULATECOEFFSINDICESBSPLINE1ST2NDORDER calculates the indices and their corresponding
%interpolation coefficients
%
%
% calcCoeffsIndicesBSpline1st2ndOrder calculates the the indices of the grid points and their
% corresponding interpolation coefficients for the voxel encompassing a target
% point. A local B-spline interpolation is used.
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
%                               interpolation coefficients of the field
%      raytogrid_coeff_derivative_matrix - matrix for calculating B-spline
%                                         interpolation coefficients of the
%                                         directional gradients of the field
%      paraxial        - Boolean indicating whether the traced ray is paraxial or not
% OPTIONAL INPUTS:

%
% OUTPUTS:
%
%       indices         - the indices of the voxel at which the target point resides
%       coeff           - the interpolation coefficient of the indices used
%                         for the interpolationof the refractive index
%       coeff_derivative_x - the first derivative along the x coordinate of
%                            the interpolation coefficients
%       coeff_derivative_x - the first derivative along the y coordinate of
%                            the interpolation coefficients
%       coeff_derivative_z - the first derivative along the z coordinate of
%                            the interpolation coefficients
%       coeff_second_derivative_xx - the second derivative along the xx cooordinates
%                                    of the interpolation coefficients
%       coeff_second_derivative_yy - the second derivative along the yy cooordinates
%                                    of the interpolation coefficients
%       coeff_second_derivative_xy - the second derivative along the xy cooordinates
%                                    of the interpolation coefficients
%       coeff_second_derivative_xz - the second derivative along the xz cooordinates
%                                    of the interpolation coefficients
%       coeff_second_derivative_yz - the second derivative along the yz cooordinates
%                                    of the interpolation coefficients
%       coeff_second_derivative_zz - the second derivative along the zz cooordinates
%                                    of the interpolation coefficients
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian



% get the indices of the first grid point included in the interpolation along x,
% y, and z dimensions
[indices, index_x, index_y, index_z] = calcInterpIndicesBspline(pos_point,...
    pos_grid_first, dx, grid_size, dim, raytogrid_indices_x,...
    raytogrid_indices_y, raytogrid_indices_z);

% compute the interpolation coefficients, if at least one grid point from the grid cell
% encompassing the ray point are located inside the binary mask for ray
% tracing
if all(indices - numel(mask)<0)
    if any(mask(indices))
        
        % get the distance of the off-grid point and the first grid point
        % divided by the grid spacing along the x coordinate
        xd = (pos_point(1) - xvec(index_x))/dx;
        
        % get the distance of the off-grid point and the first grid point
        % divided by the grid spacing along the y coordinate
        yd = (pos_point(2) - yvec(index_y))/dx;
        
        % get the interpolation polynomials along the x coordinate
        polynomial_x = raytogrid_coeff_matrix * [xd^3; xd^2; xd; 1];
        
        % get the interpolation polynomials along the y coordinate
        polynomial_y = raytogrid_coeff_matrix * [yd^3; yd^2; yd; 1];
        
        % get the first derivative along the x coordinate of the interpolation polynomials
        polynomial_derivative_x = raytogrid_coeff_derivative_matrix * [xd^2; xd; 1];
        
        % get the first derivative along the y coordinate of the interpolation polynomials
        polynomial_derivative_y = raytogrid_coeff_derivative_matrix * [yd^2; yd; 1];
        
        if paraxial
            
            % get the second derivative along the x coordinate of the interpolation polynomials
            polynomial_second_derivative_x = raytogrid_coeff_second_derivative_matrix * [xd; 1];
            
            % get the second derivative along the y coordinate of the interpolation polynomials
            polynomial_second_derivative_y = raytogrid_coeff_second_derivative_matrix * [yd; 1];
        end
        
        switch dim
            
            case 2
                
                % get the interpolation coefficients
                coeff = (vectorise(polynomial_x * polynomial_y'))';
                
                % get the first derivative along the x coordinate of the interpolation coefficients
                coeff_derivative_x = (vectorise(polynomial_derivative_x * polynomial_y'))';
                
                % get the first derivative along the y coordinate of the interpolation coefficients
                coeff_derivative_y = (vectorise(polynomial_x * polynomial_derivative_y'))';
                
                % get the first derivative along the z coordinate of the interpolation coefficients
                coeff_derivative_z = [];
                
                
                if paraxial
                    
                    % get the second derivative along the xx coordinates of the interpolation coefficients
                    coeff_second_derivative_xx = (vectorise(polynomial_second_derivative_x * polynomial_y')).';
                    % get the second derivative along the yy coordinates of the interpolation coefficients
                    coeff_second_derivative_yy = (vectorise(polynomial_x * polynomial_second_derivative_y')).';
                    % get the second derivative along the xy coordinates of the interpolation coefficients
                    coeff_second_derivative_xy = (vectorise(polynomial_derivative_x * polynomial_derivative_y')).';
                    
                else
                    
                    % allocate empty variable
                    coeff_second_derivative_xx = [];
                    coeff_second_derivative_yy = [];
                    coeff_second_derivative_xy = [];
                    
                end
                
                % get the second derivative along the xz coordinates of the interpolation coefficients
                coeff_second_derivative_xz = [];
                % get the second derivative along the yz coordinates of the interpolation coefficients
                coeff_second_derivative_yz = [];
                % get the second derivative along the zz coordinates of the interpolation coefficients
                coeff_second_derivative_zz = [];
                
                
                
            case 3
                
                % get the distance of the target point and the first grid point
                % divided by the grid spacing along the z coordinate
                zd = (pos_point(3) - zvec(index_z))/dx;
                
                % get the interpolation polynomials along the z coordinate
                polynomial_z = raytogrid_coeff_matrix * [zd^3; zd^2; zd; 1];
                
                % get the first derivative along the z coordinate of the interpolation polynomials
                polynomial_derivative_z = raytogrid_coeff_derivative_matrix * [zd^2; zd; 1];
                
                % get the second derivative along the z coordinate of the interpolation polynomials
                polynomial_second_derivative_z = raytogrid_coeff_second_derivative_matrix * [zd; 1];
                
                % get the interpolation coefficients for interpolation of the
                % map
                coeff = vectorise(polynomial_x * polynomial_y' .*...
                    reshape(polynomial_z, 1, 1, []) ).';
                
                % get the interpolation coefficients for interpolation of the
                % the first derivative along
                % x coordinate
                coeff_derivative_x = vectorise(polynomial_derivative_x * polynomial_y' .*...
                    reshape(polynomial_z, 1, 1, []) ).';
                
                % y coordinate
                coeff_derivative_y = vectorise(polynomial_x * polynomial_derivative_y' .*...
                    reshape(polynomial_z, 1, 1, []) ).';
                
                % z coordinate
                coeff_derivative_z = vectorise(polynomial_x * polynomial_y' .*...
                    reshape(polynomial_derivative_z, 1, 1, []) ).';
                
                
                
                % get the interpolation polynomials along the x and x=y
                % coordinates
                % p_xy = polynomial_x * polynomial_y';
                % p_dxy = polynomial_derivative_x * polynomial_y';
                % p_xdy = polynomial_x * polynomial_derivative_y';
                
                % coeff = zeros(4,4,4);
                % coeff_derivative_x = zeros(4,4,4);
                % coeff_derivative_y = zeros(4,4,4);
                % coeff_derivative_z = zeros(4,4,4);
                
                %for i = 1:4
                %    coeff(:,:,i) = polynomial_z(i) * p_xy;
                %    coeff_derivative_x(:,:,i) = polynomial_z(i) * p_dxy;
                %    coeff_derivative_y(:,:,i) = polynomial_z(i) * p_xdy;
                %    coeff_derivative_z(:,:,i) = polynomial_derivative_z(i) * p_xy;
                %end
                
                %coeff = (coeff(:))';
                %coeff_derivative_x = (coeff_derivative_x(:))';
                %coeff_derivative_y = (coeff_derivative_y(:))';
                %coeff_derivative_z = (coeff_derivative_z(:))';
                
                
                
                if paraxial
                    
                    % get the second derivative of the map along
                    % xx coordinates
                    coeff_second_derivative_xx = vectorise(polynomial_second_derivative_x * polynomial_y'.*...
                        reshape(polynomial_z, 1, 1, [])).';
                    % xy coordinates
                    coeff_second_derivative_xy = vectorise(polynomial_derivative_x * polynomial_derivative_y' .*...
                        reshape(polynomial_z, 1, 1, [])).';
                    % xz coordinates
                    coeff_second_derivative_xz = vectorise(polynomial_derivative_x * polynomial_y'.*...
                        reshape(polynomial_derivative_z, 1, 1, [])).';
                    
                    % yy coordinates
                    coeff_second_derivative_yy = vectorise(polynomial_x * polynomial_second_derivative_y'.*...
                        reshape(polynomial_z, 1, 1, [])).';
                    % yz coordinates
                    coeff_second_derivative_yz = vectorise(polynomial_x * polynomial_derivative_y'.*...
                        reshape(polynomial_derivative_z, 1, 1, [])).';
                    
                    % zz coordinates
                    coeff_second_derivative_zz = vectorise(polynomial_x * polynomial_y'.*...
                        reshape(polynomial_second_derivative_z, 1, 1, [])).';
                
                    
                    % get the second derivative along the yy coordinates of the interpolation coefficients
                    %     coeff_second_derivative_yy(:,:,i) = polynomial_z(i) * p_xdy2;
                    % get the second derivative along the xy coordinates of the interpolation coefficients
                    %     coeff_second_derivative_xy(:,:,i) = polynomial_z(i) * p_dxdy;
                    % get the second derivative along the xz coordinates of the interpolation coefficients
                    %     coeff_second_derivative_xz(:,:,i) = polynomial_derivative_z(i) * p_dxy;
                    % get the second derivative along the yz coordinates of the interpolation coefficients
                    %     coeff_second_derivative_yz(:,:,i) = polynomial_derivative_z(i) * p_xdy;
                    % get the second derivative along the zz coordinates of the interpolation coefficients
                    %     coeff_second_derivative_zz(:,:,i) = polynomial_second_derivative_z * p_xy;
                    
                    
                    
                    
                    
                    
                    
                    % p_dx2y = polynomial_second_derivative_x * polynomial_y';
                    % p_xdy2 = polynomial_x * polynomial_second_derivative_y';
                    % p_dxdy = polynomial_derivative_x * polynomial_derivative_y';
                    
                    % coeff_second_derivative_xx = zeros(4,4,4);
                    % coeff_second_derivative_yy = zeros(4,4,4);
                    % coeff_second_derivative_xy = zeros(4,4,4);
                    % coeff_second_derivative_xz = zeros(4,4,4);
                    % coeff_second_derivative_yz = zeros(4,4,4);
                    % coeff_second_derivative_zz = zeros(4,4,4);
                    
                    
                    % for i = 1:4
                    
                    % get the second derivative along the xx coordinates of the interpolation coefficients
                    %     coeff_second_derivative_xx(:,:,i) = polynomial_z(i) * p_dx2y;
                    % get the second derivative along the yy coordinates of the interpolation coefficients
                    %     coeff_second_derivative_yy(:,:,i) = polynomial_z(i) * p_xdy2;
                    % get the second derivative along the xy coordinates of the interpolation coefficients
                    %     coeff_second_derivative_xy(:,:,i) = polynomial_z(i) * p_dxdy;
                    % get the second derivative along the xz coordinates of the interpolation coefficients
                    %     coeff_second_derivative_xz(:,:,i) = polynomial_derivative_z(i) * p_dxy;
                    % get the second derivative along the yz coordinates of the interpolation coefficients
                    %     coeff_second_derivative_yz(:,:,i) = polynomial_derivative_z(i) * p_xdy;
                    % get the second derivative along the zz coordinates of the interpolation coefficients
                    %     coeff_second_derivative_zz(:,:,i) = polynomial_second_derivative_z * p_xy;
                    % end
                    
                    
                else
                    
                    % allocate empty variable
                    coeff_second_derivative_xx = [];
                    coeff_second_derivative_yy = [];
                    coeff_second_derivative_xy = [];
                    coeff_second_derivative_xz = [];
                    coeff_second_derivative_yz = [];
                    coeff_second_derivative_zz = [];
                    
                    
                end
                
                
        end
        
        
    else
        
        coeff = [];
        coeff_derivative_x = [];
        coeff_derivative_y = [];
        coeff_derivative_z = [];
        coeff_second_derivative_xx = [];
        coeff_second_derivative_yy = [];
        coeff_second_derivative_xy = [];
        coeff_second_derivative_xz = [];
        coeff_second_derivative_yz = [];
        coeff_second_derivative_zz = [];
        
    end
    
else
    
    
    coeff = [];
    coeff_derivative_x = [];
    coeff_derivative_y = [];
    coeff_derivative_z = [];
    coeff_second_derivative_xx = [];
    coeff_second_derivative_yy = [];
    coeff_second_derivative_xy = [];
    coeff_second_derivative_xz = [];
    coeff_second_derivative_yz = [];
    coeff_second_derivative_zz = [];
end




end