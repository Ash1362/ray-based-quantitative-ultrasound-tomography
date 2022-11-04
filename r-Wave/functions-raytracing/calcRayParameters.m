function [polar_direction_residual, interp_coeff_vec, cartesian_pos_endpoint,...
    polar_direction_endpoint] = calcRayParameters(refractive,refractive_gradient_x,...
    refractive_gradient_y, refractive_gradient_z, cartesian_position_emitter,...
    polar_direction_receiver, polar_initial_direction, xvec, yvec, zvec,...
    pos_grid_first, pos_grid_end, dx, ds, grid_size, dim, detec_radius,...
    mask, raytogrid_indices_x, raytogrid_indices_y, raytogrid_indices_z,...
    raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix, calc_coeffs,...
    raylinking_method, interp_method)
%CALCRAYPARAMETERS traces a ray and calculate the information about the end
%point of the ray
%
% DESCRIPTION:
% calRayParameters traces a ray and calculate the information about the end
% point of the ray. This isnformation will be used for solving the ray
% linking inverse problem
% USAGE:
%
%
% INPUTS:
%       refractive            - disretised refractive index
%       refractive_gradient_x - discretised refractive index gradient along x
%       refractive_gradient_y - discretised refractive index gradient along y
%       refractive_gradient_z - discretised refractive index gradient along z
%       cartesian_position_emitter - a dim x 1 cartesian position of the emitter
%       polar_direction_receiver   - a (dim-1) x 1 vector of the polar direction from emitter to
%                                    receiver
%       polar_initial_direction    - a (dim-1) x 1 vector of the polar initial
%                                    direction of the ray
%       xvec                 - the x vector of grid points
%       yvec                 - the y vector of grid points
%       zvec                 - the z vector of grid points
%       pos_grid_first       - a dim x 1 Cartesian position of the first index of
%                               the grid
%       pos_grid_end         - a dim x 1 Cartesian position of the end index of
%                              the grid
%       dx                   - a scalar representing the grid spacing, the same
%                              along all the Cartesian coordinates [m]
%       ds                   - a saclar representing the ray spacing [m]
%       grid_size            - the size of the grid
%       dim                  - the dimension of the medium
%       mask                 - a binary mask used for calculating the
%                              interpolation coefficients
%      raytogrid_indices_x   -  x indices for B-spline interpolation
%      raytogrid_indices_y   -  y indices for B-spline interpolation                                  
%      raytogrid_indices_z   -  z indices for B-spline interpolation
%      raytogrid_coeff_matrix - matrix for calculating B-spline
%                               interpolation coefficients of the field
%      raytogrid_coeff_derivative_matrix - matrix for calculating B-spline
%                                         interpolation coefficients of the
%                                         directional gradients of the field                              
%       calc_coeffs          - a boolean controlling whether the
%                              interpolation coefficients are stored or
%                              not (This is often set true only for a linked
%                              ray after solving the ray linking problem.)
%       interp_method        - method for interpolation
% OPTIONAL INPUTS:
%
% OUTPUTS:
%       polar_direction_residual - dim-1 x 1 matrix of polar residual, the discrepancy
%                                  between the polar unit vector from emitter
%                                  to the last point of the ray and a polar
%                                  unit vector from emitter to the receiver
%       interp_coeff_vec         - interpolation coefficients for the grid
%                                  points inside the mask
%       cartesian_pos_endpoint   - dim x 1 matrix of cartesian position of the end point of the ray
%       polar_direction_endpoint - dim-1 x 1 matrix of polar directions of 
%                                  geometrical unit vectors from the emitter
%                                  to the end point of the rays

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2020 Ashkan Javaherian 


% allocate a dim x 1 vector for
cartesian_initial_direction = zeros(dim, 1);

% transform the initial direction from polar to casrtesian coordinate
switch dim
    case 2
        
        switch raylinking_method
            case 'Secant'
                [cartesian_initial_direction(1), cartesian_initial_direction(2)] = ...
                    pol2cart(polar_initial_direction, 1);
            case 'Regula-Falsi'
                % Using 'Regula falsi' method, the origin of polar coordinate 
                % is set along a line segment from emitter to the
                % origin, so the the polar direction of a geometrical vector
                % starting from the emitter will be the angle of the
                % corresponding vector with a line segment connecting the 
                % emitter to the origin
                 cartesian_initial_direction = rotateRodriguez(- cartesian_position_emitter,...
                    polar_initial_direction);
        end
        
        
        
    case 3
        
        [cartesian_initial_direction(1), cartesian_initial_direction(2),...
            cartesian_initial_direction(3)] = sph2cart(polar_initial_direction(1),...
            polar_initial_direction(2), 1);
        
        
end

% trace the ray
% calculate the cartesian points along the ray and their associated
% interpolation coefficients
switch interp_method
    case 'Bilinear'
[cartesian_pos_endpoint, interp_coeff_vec] = calcRay(...
    refractive, refractive_gradient_x, refractive_gradient_y,...
    refractive_gradient_z, cartesian_position_emitter, cartesian_initial_direction,...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds,...
    grid_size, dim, detec_radius, mask, calc_coeffs);
    case 'Bspline'
[cartesian_pos_endpoint, interp_coeff_vec] = calcRayRungeKutta2ndBspline(... % calcRayBspline(... %   
    refractive, cartesian_position_emitter, cartesian_initial_direction, ...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds,...
    grid_size, dim, detec_radius, mask,...
    raytogrid_indices_x, raytogrid_indices_y,...
     raytogrid_indices_z, raytogrid_coeff_matrix, ...
     raytogrid_coeff_derivative_matrix, calc_coeffs);
end

% the cartesian direction of the geometrical vector linking the emitter
% to the last point of the ray
cartesian_direction_endpoint = cartesian_pos_endpoint - cartesian_position_emitter;


% the polar direction of the geometrical vector linking the emitter
% to the last point of the ray
switch dim
    case 2
        
        switch raylinking_method
            
            case 'Secant'
                [polar_direction_endpoint, ~] = cart2pol(cartesian_direction_endpoint(1),...
                    cartesian_direction_endpoint(2));
                
            case 'Regula-Falsi'
                polar_direction_endpoint = calcDirectionalAngle(-cartesian_position_emitter,...
                    cartesian_direction_endpoint);
                
        end
        
        
    case 3
        polar_direction_endpoint = zeros(dim-1, 1);
        [polar_direction_endpoint(1), polar_direction_endpoint(2)] =...
            cart2sph(cartesian_direction_endpoint(1),...
            cartesian_direction_endpoint(2),...
            cartesian_direction_endpoint(3));
end


% calculate the polar residual, the discrepancy between the polar
% unit vector from emitter to the receiver and polar unit vector from
% emitter to the last point of the ray, if requested
if ~isempty(polar_direction_receiver)
    polar_direction_residual = polar_direction_endpoint - polar_direction_receiver;
    
    % wrap the first component into the interval [-pi,+pi] radian
    polar_direction_residual(1) = wrapToPi(polar_direction_residual(1));
else
    polar_direction_residual = [];
end





end