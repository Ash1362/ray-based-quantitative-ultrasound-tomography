function [pos, ray_positions, ray_acoustic_length, ray_absorption, ds_final] =...
    calcRayRungeKutta2ndBsplineFullWaveModifiedd(refractive, cartesian_position_emitter,...
    direction, direction_perturbation, xvec, yvec, zvec,...
    pos_grid_first, pos_grid_end, dx, ds, grid_size, dim, detec_radius, mask,...
    raytogrid_indices_x, raytogrid_indices_y, raytogrid_indices_z, raytogrid_coeff_matrix,...
    raytogrid_coeff_derivative_matrix, raytogrid_coeff_second_derivative_matrix,...
    refractive_nonsmoothed, absorption_coeff, calc_coeffs, auxiliary_ray, auxiliary_method)
%calcRayBsplineFullWave traces a ray, stores the ray-to-grid interpolation
%coefficients, and integrate the field's parameters along the rays
%
% DESCRIPTION:
% calcRayBsplineFullWave computes the trajectory of rays, and the ray-to-grid
% interplation coefficients for the computed points, given an initial position and an
% initial direction for the ray. The interpolation coefficients will be stored
% in a sparse matrix, and will be used for image reconstruction based on time-of-flights.
% The Cartesinan position of ray's points, together with the accumulated information along the ray,
% will be stored for image reconstruction of the sund speed using Green's
% inversion approach, which will incorporate the scattering effects into the
% image reconstruction.
%
% INPUTS:
%       refractive       - disretised refractive index
%       cartesian_position_emitter - a dim x 1 Cartesian position of the emitter
%       direction        - a dim x 1 vector of the cartesian initial direction along
%                         the ray
%       xvec             - the x vector of grid points
%       yvec             - the y vector of grid points
%       zvec             - the z vector of grid points
%       pos_grid_first   - a dim x 1 Cartesian position of the first index of
%                         the grid
%       pos_grid_end     - a dim x 1 Cartesian position of the end index of
%                         the grid
%       dx               - a scalar representing the grid spacing [m]
%       ds               - a saclar representing the ray spacing [m]
%       grid_size        - the grid size
%       dim              - the dimension of the medium
%       mask             - a binary mask used for ray tracing
%      raytogrid_indices_x -  x indices for B-spline interpolation
%      raytogrid_indices_y -  y indices for B-spline interpolation
%      raytogrid_indices_z -  z indices for B-spline interpolation
%      raytogrid_coeff_matrix     - matrix for calculating B-spline
%                                   interpolation coefficients of the field
%      raytogrid_coeff_derivative_matrix - matrix for calculating B-spline
%                                         interpolation coefficients of the
%                                         directional gradients of the field
%      refractive_nonsmoothed - the original (nonsmoothed) refractive index
%                                for integrating the field's parameters
%                                along the rays
%      absorption_coeff       - the matrix of absorption coefficient. If it
%                              is set a scaler, the acoustic absorption is
%                              fully neglected.
%      calc_coeffs           - a boolean controlling whether the
%                               interpolation coefficients are stored or
%                               not (This is often set true only for the optimal
%                               linked ray after solving the ray linking problem.)
%       auxiliary_ray        - Boolean indicating whether the ray to be
%                              traced is auxiliary ray or not
%       auxiliary_method     - the method for tracing the auxiliary rays
% OPTIONAL INPUTS:
%
% OUTPUTS:
%
%      pos              - a dim x 1 vector of the Cartesian position of the
%                           end point along the ray
%      interp_coeff_vec - the vector of ray-to-grid interpolation coefficients for the points
%                         along the ray within the mask
%      ray_positions    - the Cartesian position of the ray's points
%      ray_acoustic_length - the accumulated acoustic length on the ray's points
%      ray_absorption   - the accumulated acoustic absorption on the rays'
%                         points
%      ds_final         - the final ray spacing between the receiver and
%                         the ray's point just next to that

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian


paraxial = auxiliary_ray & strcmp(auxiliary_method, 'paraxial');

if ~isscalar(absorption_coeff)
    do_absorption = true;
else
    do_absorption = false;
    ray_absorption = [];
end


if size(cartesian_position_emitter, 2) > 1
    
    % set the Boolean for computing the auxiliary adjoint ray true,
    % if the auxiliary adjoint ray is paraxial.  For this case, the
    % paraxial ray is traced along a given reference ray
    auxiliary_adjoint_paraxial = true;
    
    % initialise the point index
    point_index = 1;
    
else
    
   % the ray is not an auxiliary paraxial adjoint ray 
   auxiliary_adjoint_paraxial = false; 
    
end

if ~paraxial && auxiliary_adjoint_paraxial
    error(['If the position of the main ray is given, the ray must be auxiliary,'...
        'and the approach for tracing the ray must be set paraxial.'])
end



% remove the last index from the size of the grid (not used)
grid_size = grid_size(1: dim-1);


if auxiliary_adjoint_paraxial
    
    % initialise the position of the ray with the emission point
    ray_positions = cartesian_position_emitter;
    
    % get the current position
    pos = ray_positions(:, point_index);
    
    % increase the index of the ray point by 1
    point_index = point_index + 1;
    
    % get the index of the last point along the ray (the index of receiver)
    lastpoint_index = find(~isnan(ray_positions(1, :)), 1, 'last');
    
else
    
    % get the current position
    pos = cartesian_position_emitter;
    
    % set the point index 0 such that the stopping criterion of the while loop
    % is always satisfied.
    point_index = 0;
    
    % get the receiver index
    lastpoint_index = inf;
    
end

clear cartesian_position_emitter

% calculate the interpolation coefficients, interpolated parameters on the ray's point,
% and their corresponding indices on the grid
[~, coeff, n, dn, d2n, ns, as] = interpLocalBsplineFullWaveModified(pos, xvec, yvec, zvec, ...
    pos_grid_first, dx, grid_size, dim, mask, refractive, refractive_nonsmoothed,...
    absorption_coeff, raytogrid_indices_x, raytogrid_indices_y,...
    raytogrid_indices_z, raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix,...
    raytogrid_coeff_second_derivative_matrix, true, paraxial);

if isempty(coeff)
    n = 1; dn = 0; d2n = 0; ns = 1;
    if do_absorption
        as = 0;
    end
end


if calc_coeffs || auxiliary_ray
    
    if paraxial
        
        % get the first column of a dim * num_ray_pos vector for storing the
        % perturbed positions because of perturbation to the initial
        % direction
        ray_positions_perturbation = [0; 0];
        
    else
        
        % get the first column of a dim * num_ray_pos vector for storing the
        % positions
        ray_positions = pos;
        
    end
    
    if calc_coeffs
        
        % initialise the acoustic length along the ray
        ray_acoustic_length = 0;
        
        if do_absorption
            
            % initialise the acoustic absorption along the ray
            ray_absorption = 0;
        end
    end
    
end

% normalise direction - the direction along the ray must be a unit vector
% multiplied by the refractive index
direction = n/norm(direction) * direction;

% calculate the position component of q1
q1_x = 1/n * direction;

if ~auxiliary_adjoint_paraxial
    
    % calculate the second term for updating the direction
    direction_2 = direction + ds * dn;
end

% calculate the second term for updating the position
pos_2 = pos + ds * q1_x;

% calculate the interpolation coefficients for the second point
[~, coeff_2, n_2, dn_2, d2n_2, ~, ~] = interpLocalBsplineFullWaveModified(pos_2, xvec, yvec, zvec, ...
    pos_grid_first, dx, grid_size, dim, mask, refractive, refractive_nonsmoothed,...
    absorption_coeff, raytogrid_indices_x, raytogrid_indices_y,...
    raytogrid_indices_z, raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix,...
    raytogrid_coeff_second_derivative_matrix, true, paraxial);

if isempty(coeff_2)
    n_2 = 1; dn_2 = 0; d2n_2 = 0;
end


if ~auxiliary_adjoint_paraxial
    
    % normalise direction - the direction along the ray must be a unit vector
    % multiplied by the refractive index
    direction_2 = n_2/norm(direction_2) * direction_2;
    
    % calculate the position component of q2
    q2_x = 1/n_2 * direction_2;
    
    % store the current position
    pos_previous = pos;
    
    % update the current position of the ray based on the trapezoid rule
    pos = pos_previous + ds/2 * (q1_x + q2_x);
    
else
    
    % get the current position
    pos = ray_positions(:, point_index);
    
    % increase the index of ray point by 1
    point_index = point_index + 1;
    
end


if paraxial
    
    % calculate the perturbation position component of q1
    q1_perturb_x = 1/n * direction_perturbation;
    
    % initialise the direction perturbation component of q1
    % consider that the initial perturbation vector is [0,0]^T.
    q1_perturb_k = [0; 0];
    
    % calculate the second term for updating the perturbation direction vector.
    % consider that the initial perturbation vector is [0,0]^T.
    direction_perturbation_2 = direction_perturbation;
    
    % calculate the perturbation position component of q2
    q2_perturb_x = 1/n_2 * direction_perturbation_2;
    
    % get the current position perturbation
    pos_perturb_previous = [0; 0];
    
    % calculate the second term for updating the perturbation position vector
    % consider that the initial perturbation vector is [0,0]^T.
    pos_perturb_2 = ds * q1_perturb_x;
    
    
    % update the position perturbation of the ray based on the trapezoid rule
    pos_perturb =  ds/2 * (q1_perturb_x + q2_perturb_x);
   
    
end

% the ray must be terminated inside the detection surface
% if the ray is paraxial adjoint ray, the ray is stopped on the index associated with the last point on the ray 
while (norm(pos) - detec_radius <= -1e-10    && ...
            all(pos > pos_grid_first + dx & pos < pos_grid_end - 2 * dx))  &&  (point_index < lastpoint_index)
          
    
        if ~auxiliary_adjoint_paraxial
            
            % update the direction
            direction = direction + ds/2 * (dn + dn_2);
            
        else
            
            % get the direction
            direction = ray_positions(:, point_index) - pos;
            
        end
    
    
    
    if paraxial
        
        % update the the direction perturbation component of q2
        q2_perturb_k = (1/n_2^2 * (dn_2 * dn_2') + 1/n_2 * d2n_2) * pos_perturb_2;
        
        % update the direction perturbation based on the Trapezoid rule
        direction_perturbation = direction_perturbation...
            + ds/2 * (q1_perturb_k + q2_perturb_k);
        
    end
    
    % calculate the interpolation coefficients, interpolated parameters on the ray's point,
    % and their corresponding indices on the grid
    [~, coeff, n, dn, d2n, ns, as] = interpLocalBsplineFullWaveModified(pos, xvec, yvec, zvec, ...
        pos_grid_first, dx, grid_size, dim, mask, refractive, refractive_nonsmoothed,...
        absorption_coeff, raytogrid_indices_x, raytogrid_indices_y,...
        raytogrid_indices_z, raytogrid_coeff_matrix, ...
        raytogrid_coeff_derivative_matrix, raytogrid_coeff_second_derivative_matrix,...
        true, paraxial);
    
    
    if isempty(coeff)
        n = 1; dn = 0; d2n = 0; ns = 1;
        if do_absorption
            as = 0;
        end
       
    end
    
    if calc_coeffs || auxiliary_ray
        
        if paraxial
            
            % add the current perturbed position vector to the associated matrix
            ray_positions_perturbation = [ray_positions_perturbation, pos_perturb];
            
        else
            
            % add the current position vector to the associated matrix, if
            % the ray is not auxiliary, or is an auxiliary ray which is not
            % paraxial
            ray_positions = [ray_positions, pos];
            
        end
        
        if calc_coeffs
            
            % add the local acoustic length to the associated matrix
            ray_acoustic_length = [ray_acoustic_length, ds * ns];
            
            if do_absorption
                
                % add the local absorption to the associated matrix
                ray_absorption = [ray_absorption, ds * as];
            end
        end
        
    end
    
    
    % normalise direction - the direction along the ray must be a unit vector
    % multiplied by refractive index
    direction = n/norm(direction) * direction;
    
    % calculate the position component of q1
    q1_x = 1/n * direction;
    
    if ~auxiliary_adjoint_paraxial
        
        % calculate the second term for updating the direction
        direction_2 = direction + ds * dn;
    end
    
    % calculate the second term for updating the position
    pos_2 = pos + ds * q1_x;
    
    
    % calculate the interpolation coefficients for the second point
    [~, coeff_2, n_2, dn_2, d2n_2, ~, ~] = interpLocalBsplineFullWaveModified(pos_2, xvec, yvec, zvec, ...
        pos_grid_first, dx, grid_size, dim, mask, refractive, refractive_nonsmoothed,...
        absorption_coeff, raytogrid_indices_x, raytogrid_indices_y,...
        raytogrid_indices_z, raytogrid_coeff_matrix, ...
        raytogrid_coeff_derivative_matrix, raytogrid_coeff_second_derivative_matrix,...
        true, paraxial);
    
    if isempty(coeff_2)
        n_2 = 1; dn_2 = 0; d2n_2 = 0;
    end
    
    if ~auxiliary_adjoint_paraxial
        
        % normalise direction - the direction along the ray must be a unit vector
        % multiplied by the refractive index
        direction_2 = n_2/norm(direction_2) * direction_2;
        
        % calculate the position component of q2
        q2_x = 1/n_2 * direction_2;
        
        % store the current position
        pos_previous = pos;
        
        % update the current position of the ray based on the trapezoid rule
        pos = pos_previous + ds/2 * (q1_x + q2_x);
        
    else
        
        % get the current position
        pos = ray_positions(:, point_index);
        
        % increase the index of ray point by 1
        point_index = point_index + 1;
        
    end

    if paraxial
        
        % calculate the perturbation position component of q1
        q1_perturb_x = 1/n * direction_perturbation;
        
        % update the direction perturbation component of q1
        q1_perturb_k = (1/n^2 * (dn * dn') + 1/n * d2n) * pos_perturb;
        
        % calculate the second term for updating the perturbation direction vector.
        direction_perturbation_2 = direction_perturbation + ds * q1_perturb_k;
        
        % calculate the perturbation position component of q2
        q2_perturb_x = 1/n_2 * direction_perturbation_2;
        
        % calculate the second term for updating the perturbation position vector
        pos_perturb_2 =  pos_perturb + ds * q1_perturb_x;
        
        % get the current position perturbation
        pos_perturb_previous = pos_perturb;
        
        % update the position perturbation of the ray based on the trapezoid rule
        pos_perturb =  pos_perturb_previous + ds/2 * (q1_perturb_x + q2_perturb_x);
        
             
    end
    
  
end


%  the last point must be on the detection surface except for auxiliary
%  rays
if ~auxiliary_ray || paraxial
    
    if ~auxiliary_adjoint_paraxial
    % if the ray is not auxiliary, or the ray is an auxiliary ray which is
    % paraxial, the last point of the ray position must be on the detection
    % ring (surface)
    [pos , ds_final] = calcIntersectBall(pos_previous, 1/2 * (q1_x + q2_x),...
        detec_radius);
    else
        
        % compute the final ray spacing from the last two points along the
        % given main ray
        ds_final = norm(ray_positions(:, point_index)...
            - ray_positions(:, point_index-1));
    end
    
    if paraxial
        
        % update the last position perturbation of the ray based on the trapezoid rule
        pos_perturb =  pos_perturb_previous + ds_final/2 * (q1_perturb_x + q2_perturb_x);
        
    end
    
else
    
    % keep the last position and ray spacing
    ds_final = ds;
    
end


if all (pos > pos_grid_first + dx  &  pos < pos_grid_end - 2 * dx )
    
    
    % calculate the interpolation coefficients, interpolated parameters on the ray's point,
    % and their corresponding indices on the grid
    [~, coeff, ~, ~, ~, ns, as] = interpLocalBsplineFullWaveModified(pos, xvec, yvec, zvec, ...
        pos_grid_first, dx, grid_size, dim, mask, refractive, refractive_nonsmoothed,...
        absorption_coeff, raytogrid_indices_x, raytogrid_indices_y,...
        raytogrid_indices_z, raytogrid_coeff_matrix, ...
        raytogrid_coeff_derivative_matrix, raytogrid_coeff_second_derivative_matrix,...
        false, paraxial);
    
    if isempty(coeff)
        ns = 1;
        if do_absorption
            as  = 0;
        end
    end
    
    
    if calc_coeffs || auxiliary_ray
        
        if paraxial
            
            % add the last perturbed position vector to the associated matrix
            ray_positions_perturbation = [ray_positions_perturbation, pos_perturb];
            
        else
            % add the last position vector to the associated matrix
            ray_positions = [ray_positions, pos];
            
        end
        
        
        if calc_coeffs
            
            % add the local acoustic length of the last point to the associated
            % vector
            ray_acoustic_length = [ray_acoustic_length, ds_final * ns];
            
            if do_absorption
                
                % add the local absorption of the last point to the associated vector
                ray_absorption = [ray_absorption, ds_final * as];
                
            end
            
        end
        
        
    end
    
end


if calc_coeffs
    
    % accumulate the acoustic length
    ray_acoustic_length = cumsum(ray_acoustic_length);
    
    if do_absorption
        
        % accumulate the absorption
        ray_absorption = cumsum(ray_absorption);
        
    end
    
else
    
    % give an empty variable
    ray_acoustic_length = [];
    ray_absorption = [];
    
    if ~auxiliary_ray
        
        % if the ray is not a linked ray, and therefore, calc_coeffs is set
        % false, and also the ray is not auxiliary, give an empty variable for
        % the ray position
        ray_positions = [];
    end
    
end

if paraxial
    
    % give the perturbed ray poitions to the ray positions
    ray_positions = ray_positions_perturbation;
    
end


end