function [polar_direction_residual, interp_coeff_vec, cartesian_pos_endpoint,...
    polar_direction_endpoint] = calcRayLinkForward(ray_interp_coeffs, refractive,...
    cartesian_position_emitter, polar_direction_receiver, polar_initial_direction,...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds, grid_size, dim, detec_geom,...
    mask, calc_coeffs, raylinking_method, raytracing_method)
% CALCRAYLINKFORWARD traces a ray and gives the information about the end
% point of the ray
%
% DESCRIPTION:
% calRayLinkForward traces a ray and calculate the information about the end
% point of the ray. In other words, this function solves the forward
% operators for iteratvely solving the inverse problem of ray
% linking. For the optimal ray after ray linking, this function gives the 
% interpolation coefficients along the ray.

% USAGE:
%
%
% INPUTS:
%       ray_interp_coeffs   - a struct with fields the direction gradients
%                           for an interpolation using a'Bilinear'
%                           approach, or parameters for choosing indices of
%                           the grid points and their associated coefficients for an
%                           interpolation using a 'Bspline' approach.
%                           The fields are:
%       'refractive_gradient_x' - discretised refractive index gradient along x
%       'refractive_gradient_y' - discretised refractive index gradient along y
%       'refractive_gradient_z' - discretised refractive index gradient along z
%                              using a'Bspline' interpolation, this
%                              includes the matrices for intertpolation:
%       'raytogrid_indices_x'   - x indices for B-spline interpolation
%       'raytogrid_indices_y'   - y indices for B-spline interpolation                                  
%       'raytogrid_indices_z'   - z indices for B-spline interpolation
%     'raytogrid_coeff_matrix'  - matrix for calculating B-spline
%                                 interpolation coefficients of the field
%     'raytogrid_coeff_derivative_matrix' - matrix for calculating B-spline
%                                 interpolation coefficients of the directional 
%                                 gradientsof the field
%       refractive          - the smoothed refractive index
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
%       detec_geom           - the geometry of the detection surface
%                              with fields:
%      'radius_circle'       - the radius [m] of the circular (2D) or 
%                              hemi-spherical (3D) detection surface, or
%      'radius_cylinder'     - the radius [m] of the cylinder in x-y plane,
%                               or
%      'line_coeff'          - the coefficients [a,b,c] for equation ax+by=c of
%                              line or an intersection of a plane with x-y
%                              plane
%       mask                 - a binary mask used for calculating the
%                              interpolation coefficients                             
%       calc_coeffs          - a boolean controlling whether the
%                              interpolation coefficients are stored or
%                              not (This is often set true only for a linked
%                              ray after solving the ray linking problem.)
%       raylink_method     -  method for ray linking. For 2D case, this can be
%                              'Regula-Falsi' or 'Secant', and for 3D case, this
%                              can be 'Quasi-Newton'. 'Regula-Falsi' converges
%                              well with initial guess far from true, but it
%                              converges solwly. 'Secant' and 'Quasi-Newton'
%                              are fast, but converge badly for initial guesses
%                              far from true, and are therefore used
%                              through iteratively reconstruction of the
%                              sound speed, where the linked ray for each
%                              iteration is used as initial guess for ray
%                              linking for the next iteration.
%       raytracing_method   - the method for ray tracing, which can be
%                             'Mixed-step', 'Dual-update', 'Characteristics',
%                             or 'Runge-kutta-2nd'.


%        
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
%       date            - 30.11.2019
%       last update     - 07.09.2021
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 


% allocate a dim x 1 vector for
cartesian_initial_direction = zeros(dim, 1);

% transform the initial direction from polar to Cartesian coordinate
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


% compute the ray and its end point if the ray is computed as part of ray linking,
% or the ray-to-grid interpolation coefficients if the ray is already
% linked between an emitter-receiver pair.
switch raytracing_method
    case 'Mixed-step'
        
[cartesian_pos_endpoint, interp_coeff_vec] = calcRayMixedStep(refractive,...
    ray_interp_coeffs, cartesian_position_emitter, cartesian_initial_direction,...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds, grid_size, dim,...
    detec_geom, mask, calc_coeffs);

    case 'Dual-update'
        
[cartesian_pos_endpoint, interp_coeff_vec] = calcRayDualUpdate(refractive,...
    ray_interp_coeffs, cartesian_position_emitter, cartesian_initial_direction,...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds, grid_size, dim,...
    detec_geom, mask, calc_coeffs);

    case 'Characteristics'
        
[cartesian_pos_endpoint, interp_coeff_vec] = calcRayCharacteristics(refractive,...
    ray_interp_coeffs, cartesian_position_emitter, cartesian_initial_direction,...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds, grid_size, dim,...
    detec_geom, mask, calc_coeffs);

    case 'Runge-kutta-2nd'
        
[cartesian_pos_endpoint, interp_coeff_vec] = calcRayRungeKutta2nd(refractive,...
    ray_interp_coeffs, cartesian_position_emitter, cartesian_initial_direction,...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds, grid_size, dim,...
    detec_geom, mask, calc_coeffs); 
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