function [polar_direction_residual, cartesian_pos_endpoint,...
    polar_direction_endpoint, ray_positions, ray_acoustic_length,...
    ray_absorption, rayspacing_receiver, ray_angles_perturbation] = ...
    calcRayLinkForwardFullWave(ray_interp_coeffs, refractive, refractive_nonsmoothed,...
    absorption_coeff, cartesian_position_emitter, polar_direction_receiver,...
    polar_initial_direction, xvec, yvec, zvec, pos_grid_first,...
    pos_grid_end, dx, ds, grid_size, dim, detec_geom, mask, calc_coeffs,...
    raylinking_method, interp_method, auxiliary_ray, auxiliary_method, do_direction_angle)
%CALCRAYLINKFORWARDFULLWAVE traces a ray and calculates the information
% about the ray.
%
% DESCRIPTION:
% calRayLinkForwardFullWave traces a ray and calculates the information
% about the end point of the ray. In other words, this function solves the
% forward operators for iteratvely solving the inverse problem of ray
% linking. Only for the linked (optimal ray after ray linking),
% this function computes the accumulated parameters along the rays, as well
% as the information required for computing the geometerical attenuation.

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
%                                 interpolation coefficients for the first-order 
%                                 gradient of the field
%     'raytogrid_coeff_second_derivative_matrix' - matrix for calculating B-spline
%                                 interpolation coefficients for the second_order 
%                                 gradient of the field
%       refractive          - the smoothed refractive index sued for
%                             computing rays' trajectories
%       refractive_nonsmoothed - the nonsmoothed refractive index matrix used
%                                for integration along the linked rays
%       absorption_coeff       - the nonsmoothed absorption coefficient matrix
%                                used for integration along the linked rays
%       cartesian_position_emitter - a dim x 1 cartesian position of the emitter
%       polar_direction_receiver   - a (dim-1) x 1 vector of the polar direction from emitter to
%                                    receiver
%       polar_initial_direction    - a (dim-1) x 1 vector of the polar initial
%                                    direction of the ray
%       rotation_matrix       - the rotation matrix for applying perturbation to the
%                               initial position and direction of the ray
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
%                              accumulated parameters along the ray are
%                              computed and stored or not.
%       calc_coeffs          - a boolean controlling whether the
%                              parameters of the Green's function are
%                              computed or not. This is set true for the
%                              linked (optimal) ray.
%       interp_method        - method for interpolation
%       auxiliary_ray        - Boolean indicating whether the ray to be
%                              traced is auxiliary ray or not
%       auxiliary_method     - the method for tracing the auxiliary rays


% OPTIONAL INPUTS:
%
% OUTPUTS:
%       polar_direction_residual - dim-1 x 1 matrix of polar residual, the discrepancy
%                                  between the polar unit vector from emitter
%                                  to the last point of the ray and a polar
%                                  unit vector from emitter to the receiver
%       cartesian_pos_endpoint   - dim x 1 matrix of cartesian position of the end point of the ray
%       polar_direction_endpoint - dim-1 x 1 matrix of polar directions of
%                                  geometrical unit vectors from the emitter
%                                  to the end point of the rays
%       ray_positions - the dim x num_raypoints Cartesian position of the
%                       rays' points. This parameter is computed only for
%                       the linked ray (optimal ray after ray linking), or
%                       the auxiliary rays.
%       ray_acoustic_length - the acoustic length along the ray. This
%                       parameter is only computed for the linked ray.
%       ray_absorption - the accumulated acoustic absorption along the
%                        ray (only computed for the linked ray)
%       rayspacing_receiver - the ray spacing between the points n and n-1
%                        along the ray, where n is the index of the end point,
%                        i.e., the point on the receiver. (only computed
%                        for the linked ray.)
%       ray_angles_perturbation - the perturbation to the ray's directions 
%                         in the polar coordinates because of a perturbation
%                         to the initial position. 
%                         This output variable is computed if auxiliary_ray
%                         is true, and auxiliary_method is set 'paraxial'.
%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2020
%       last update     - 30.12.2020
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian 

if ~(auxiliary_ray && strcmp(auxiliary_method, 'paraxial'))

    % if the traced ray is not auxiliary, or the approach for tracing auxiliary 
    % ray is not 'paraxial', the rays' direction perturbations are not required for computing.
    do_direction_angle = false;

end

if dim == 3 &&  auxiliary_ray && ~any(strcmp(auxiliary_method, {'paraxial1', 'paraxial2'}))
    error(['If the number of dimensions are three, the auxiliary method must be set'...
       'paraxial1 of paraxial2']) 
end




if size(cartesian_position_emitter, 2) < 2

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
        
         % get the Cartesian initial direction
        [cartesian_initial_direction(1), cartesian_initial_direction(2),...
            cartesian_initial_direction(3)] = sph2cart(polar_initial_direction(1),...
            polar_initial_direction(2), 1);
        
        
end

else
     
    % get the Cartesian initial direction
    cartesian_initial_direction = 1/ds * (cartesian_position_emitter(:, 2)...
        - cartesian_position_emitter(:, 1));
    
end




 % get the initial perturbation to the position of the ray 
 % That must be normal to the initial direction of the ray 
if auxiliary_ray
    switch auxiliary_method
        case 'paraxial'
                    cartesian_initial_direction_perturbation = [-cartesian_initial_direction(2);...
                        cartesian_initial_direction(1)];
                    
        case 'paraxial1'
            
                    cartesian_initial_direction_perturbation = [0; cartesian_initial_direction(3);...
                        -cartesian_initial_direction(2)];
                    
                    auxiliary_method = 'paraxial';
                    
        case 'paraxial2'
            
                  cartesian_initial_direction_perturbation = [cartesian_initial_direction(3); 0;...
                        -cartesian_initial_direction(1)];
                   auxiliary_method = 'paraxial'; 
                   
        otherwise
            
           % the perturbation to the initial cartesian direction of the ray
           cartesian_initial_direction_perturbation = []; 

    end
else
    
    % get an empty variable for the perturbation vector for the
    % initial direction
    cartesian_initial_direction_perturbation = [];
    
    % set the auxiliary method nan
    auxiliary_method = nan;
    
end



switch interp_method
    case 'Bilinear'
        error(['Bilinear interpolation is not used for the Greens function approach'...
            'The user must use the Bspline-Denis approach.'])
        
    case {'Bspline'}
       
        % compute the ray and its end point if the ray is computed for solving
        % the forward operators for iteratively solving the inverse problem of
        % ray linking. Only for the linked ray (optimal ray), get the parameters
        % of the Green's function along the ray
        [cartesian_pos_endpoint, ray_positions, ray_acoustic_length,...
            ray_absorption, ray_angles_perturbation, rayspacing_receiver] = calcRayRungeKutta2ndBsplineFullWave(...
            ray_interp_coeffs, refractive, refractive_nonsmoothed, absorption_coeff,...
            cartesian_position_emitter, cartesian_initial_direction,...
            cartesian_initial_direction_perturbation, xvec, yvec, zvec, pos_grid_first, ...
            pos_grid_end, dx, ds, grid_size, dim, detec_geom, mask,...
            calc_coeffs, auxiliary_ray, auxiliary_method, do_direction_angle);
       
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