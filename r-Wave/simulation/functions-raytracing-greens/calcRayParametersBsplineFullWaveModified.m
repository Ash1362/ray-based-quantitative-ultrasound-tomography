function [polar_direction_residual, cartesian_pos_endpoint,...
    polar_direction_endpoint, ray_positions, ray_acoustic_length,...
    ray_absorption, rayspacing_receiver] = ...
    calcRayParametersBsplineFullWaveModified(refractive, cartesian_position_emitter,...
    polar_direction_receiver, polar_initial_direction, rotation_matrix,...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds, grid_size, dim,...
    detec_radius, mask, raytogrid_indices_x, raytogrid_indices_y,...
    raytogrid_indices_z, raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix, ...
    raytogrid_coeff_second_derivative_matrix, refractive_nonsmoothed,...
    absorption_coeff, calc_coeffs, raylinking_method, interp_method,...
    auxiliary_ray, auxiliary_method)
%CALCRAYPARAMETERS traces a ray and calculate the information about the end
%point of the ray
%
% DESCRIPTION:
% calRayParameters traces a ray and calculate the information about the end
% point of the ray. This information will be used for solving the ray
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
%                              the refractive index in water can be assumed the
%                              same, and thus the associated line integarals are
%                              cancelled for a difference imaging.)
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
%       auxiliary_ray        - Boolean indicating whether the ray to be
%                              traced is auxiliary ray or not
%       auxiliary_method     - the method for tracing the auxiliary rays
%       do_perturb_initial_psoition - Boolean controlling whether the
%                              perturbed direction of the ray from perturbation
%                              to the initial position is computed or not
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

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian and Ben Cox

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

if auxiliary_ray
    switch auxiliary_method
        case 'paraxial'
            
            % get the perturbation vector for the initial direction of the
            % ray, if the method for tracing auxiliary rays is paraxial
            cartesian_initial_direction_perturbation = rotation_matrix...
                * cartesian_initial_direction;
            
            % adjust the amplitude of the perturbation to the direction
            cartesian_initial_direction_perturbation = ...
                - cartesian_initial_direction_perturbation/...
                norm(cartesian_initial_direction_perturbation);
            
            
        case 'angle_perturbation'
            
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



% trace the ray
% calculate the cartesian points along the ray and their associated
% interpolation coefficients
switch interp_method
    case 'Bilinear'
        error(['Bilinear interpolation is not used for the Greens function approach'...
            'The user must use the Bspline-Denis approach.'])
        
    case {'Bspline'}
       % [cartesian_pos_endpoint, interp_coeff_vec, ray_positions, ray_acoustic_length,...
       %     ray_absorption, ray_directions, rayspacing_receiver] = calcRayRungeKutta2ndBsplineFullWaveModified5(...
       %    refractive, cartesian_position_emitter, cartesian_initial_direction, ...
       %     cartesian_initial_direction_perturbation, initial_derivative_position,...
       %     xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds, grid_size, dim, detec_radius,...
       %     mask, raytogrid_indices_x, raytogrid_indices_y, raytogrid_indices_z,...
       %     raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix,...
       %     raytogrid_coeff_second_derivative_matrix, refractive_nonsmoothed,...
       %     absorption_coeff, calc_coeffs, auxiliary_ray, auxiliary_method);
       
       
       [cartesian_pos_endpoint, ray_positions, ray_acoustic_length,...
           ray_absorption, rayspacing_receiver] = calcRayRungeKutta2ndBsplineFullWaveModifiedd(...
           refractive, cartesian_position_emitter, cartesian_initial_direction, ...
           cartesian_initial_direction_perturbation, xvec, yvec, zvec, pos_grid_first, ...
           pos_grid_end, dx, ds, grid_size, dim, detec_radius, mask,...
           raytogrid_indices_x, raytogrid_indices_y, raytogrid_indices_z,...
           raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix,...
           raytogrid_coeff_second_derivative_matrix, refractive_nonsmoothed,...
           absorption_coeff, calc_coeffs, auxiliary_ray, auxiliary_method);
       
 
       
       
       
    %    [cartesian_pos_endpoint, interp_coeff_vec, ray_positions, ray_acoustic_length,...
    %        ray_absorption, rayspacing_receiver] = calcRayRungeKutta2ndBsplineFullWaveModified(...
    %        refractive, cartesian_position_emitter, cartesian_initial_direction, ...
    %        cartesian_initial_direction_perturbation,...
    %        xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds, grid_size, dim, detec_radius,...
    %       mask, raytogrid_indices_x, raytogrid_indices_y, raytogrid_indices_z,...
    %        raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix,...
    %        raytogrid_coeff_second_derivative_matrix, refractive_nonsmoothed,...
    %        absorption_coeff, calc_coeffs, auxiliary_ray, auxiliary_method);
        
        %[cartesian_pos_endpoint, interp_coeff_vec, ray_positions, ray_acoustic_length,...
        %   ray_absorption, rayspacing_receiver] = calcRayRungeKutta2ndBsplineFullWave(...
        %   refractive, cartesian_position_emitter, cartesian_initial_direction, ...
        %   xvec, yvec, zvec, pos_grid_first, pos_grid_end, dx, ds,...
        %   grid_size, dim, detec_radius, mask, raytogrid_indices_x, raytogrid_indices_y,...
        %   raytogrid_indices_z, raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix,...
        %   raytogrid_coeff_second_derivative_matrix, refractive_nonsmoothed,...
        %   absorption_coeff, calc_coeffs, auxiliary_ray);
end

% the cartesian direction of rhe geometrical vector linking the emitter
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