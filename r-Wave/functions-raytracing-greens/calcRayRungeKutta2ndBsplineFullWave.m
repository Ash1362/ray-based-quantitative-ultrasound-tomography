function [pos, ray_positions, ray_acoustic_length, ray_absorption,...
    ray_angles_perturbation, ds_final] =...
    calcRayRungeKutta2ndBsplineFullWave(ray_interp_coeffs, refractive,...
    refractive_nonsmoothed, absorption_coeff, cartesian_position_emitter, direction,...
    direction_perturbation, xvec, yvec, zvec, pos_grid_first, pos_grid_end,...
    dx, ds, grid_size, dim, detec_geom, mask, calc_coeffs, auxiliary_ray, ...
    auxiliary_method, do_direction_angle)
%calcRayRungeKutta2ndBsplineFullWave traces a ray, stores the ray-to-grid interpolation
%coefficients, and integrate the field's parameters along the rays
%
% DESCRIPTION:
% calcRayBsplineFullWave computes the trajectory of rays, and the ray-to-grid
% interplation coefficients for the computed points, given an initial position and an
% initial direction for the ray. The interpolation coefficients will be stored
% in a sparse matrix, and will be used for image reconstruction based on time-of-flights.
% The Cartesian position of ray's points, together with the accumulated information along the ray,
% will be stored for image reconstruction of the sund speed using Green's
% inversion approach, which will incorporate the scattering effects into the
% image reconstruction.
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
%       calc_coeffs      - a boolean controlling whether the
%                          parameters of the Green's function are
%                          computed or not. This is set true for the
%                          linked (optimal) ray.
%       auxiliary_ray    - Boolean indicating whether the ray to be
%                          traced is auxiliary ray or not
%       auxiliary_method - the method for tracing the auxiliary rays
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
%      ray_angles_perturbation - the perturbation to the ray's directions 
%                         in the polar coordinates because of a perturbation
%                         to the initial position
%      ds_final         - the final ray spacing between the receiver and
%                         the ray's point just next to that

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 05.04.2023
%
% This script is part of the r-Wave toolbox.
% Copyright (c) 2022 Ashkan Javaherian



if isfield(detec_geom, 'line_coeff')

    % set the sign for the terminating criterion for ray tracing
    if  detec_geom.line_coeff(1:2) * cartesian_position_emitter(1:2) - ...
            detec_geom.line_coeff(3) < 0
        sgn = 1;
    else
        sgn = -1;
    end

end


paraxial = auxiliary_ray & strcmp(auxiliary_method, 'paraxial');

if ~paraxial

    % the perturbation to ray's directions because of the perturbation to
    % initial position are computed.
    do_direction_angle = false;
    
end

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

    % initialise the position of the ray with the emission point
    ray_positions = cartesian_position_emitter;

    clear cartesian_position_emitter

    % get the current position
    pos = ray_positions(:, point_index);

    % increase the index of the ray point by 1
    point_index = point_index + 1;

    % get the index of the last point along the ray (the index of receiver)
    lastpoint_index = find(~isnan(ray_positions(1, :)), 1, 'last');

else

    % the ray is not an auxiliary paraxial adjoint ray
    auxiliary_adjoint_paraxial = false;

    % get the current position
    pos = cartesian_position_emitter;

    % set the point index 0 such that the stopping criterion of the while loop
    % is always satisfied for the starting point.
    point_index = 2;

    % get the receiver index
    lastpoint_index = inf;

end

if ~paraxial && auxiliary_adjoint_paraxial
    error(['If the position of the main ray is given, the ray must be auxiliary,'...
        'and the approach for tracing the ray must be set paraxial.'])
end


% get the interpolation parameters
if ~isfield(ray_interp_coeffs, 'refractive_gradient_x')

    % the indices of the grid points for interpolation along the Cartesian
    % coordinates
    raytogrid_indices_x = ray_interp_coeffs.raytogrid_indices_x;
    raytogrid_indices_y = ray_interp_coeffs.raytogrid_indices_y;
    raytogrid_indices_z = ray_interp_coeffs.raytogrid_indices_z;

    % the refractive index
    raytogrid_coeff_matrix = ray_interp_coeffs.raytogrid_coeff_matrix;

    % the first-order gradient of the refractive index
    raytogrid_coeff_derivative_matrix = ray_interp_coeffs.raytogrid_coeff_derivative_matrix;

    % the second-order gradient of the refractive index
    raytogrid_coeff_second_derivative_matrix = ray_interp_coeffs.raytogrid_coeff_second_derivative_matrix;

else

    % give an error
    error(['For the greens approach, the grid-to-ray interpolation must be'...
        'done using the Bspline approach.'])

end

% remove the last index from the size of the grid (not used)
grid_size = grid_size(1: dim-1);

% calculate the interpolation coefficients, interpolated parameters on the ray's point,
% and their corresponding indices on the grid
[~, coeff, n, dn, d2n, ns, as] = interpLocalBsplineFullWave(pos, xvec, yvec, zvec, ...
    pos_grid_first, dx, grid_size, dim, mask, refractive, refractive_nonsmoothed,...
    absorption_coeff, raytogrid_indices_x, raytogrid_indices_y,...
    raytogrid_indices_z, raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix,...
    raytogrid_coeff_second_derivative_matrix, true, paraxial);

%if isempty(coeff)
    n = 1; dn = zeros(dim, 1); d2n = 0; ns = 1;
    if do_absorption
        as = 0;
    end
%end



% normalise direction - the direction along the ray must be a unit vector
% multiplied by the refractive index
direction = n/norm(direction) * direction;

% calculate the position component of q1
q1_x = 1/n * direction;


if calc_coeffs || auxiliary_ray

    if paraxial

        % get the first column of a dim * num_ray_pos matrix for storing the
        % perturbed positions because of perturbation to the initial
        % direction
        ray_positions_perturbation = zeros(dim, 1);

        if do_direction_angle

            if dim > 2
                Error('Computing the ray perturbation is not supported for 3D case.')
            end

            % get the perturbation to initial position. This is done for the
            % second paraxial system of equations which is solved for
            % computing the perturbation to the angle of the ray's
            % directions because of a perturbation to the initial position
            % of the ray
            pos_perturb2 = [-pos(2); pos(1)];

            % set the perturbation to initial ray's direction zero
            direction_perturbation2 = zeros(2,1); 
            
            % get the first element of a 1 * num_ray_pos vector for storing the
            % angle of the perturbed ray's directions because of perturbation to
            % the initial position
            ray_angles_perturbation = 1/(n*(point_index - 1)) * [-q1_x(2), q1_x(1)] *...
                direction_perturbation2;

        end

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

% calculate the second term for updating the direction
direction_2 = direction + ds * dn;


% calculate the second term for updating the position
pos_2 = pos + ds * q1_x;

% calculate the interpolation coefficients for the second point
[~, coeff_2, n_2, dn_2, d2n_2, ~, ~] = interpLocalBsplineFullWave(pos_2, xvec, yvec, zvec, ...
    pos_grid_first, dx, grid_size, dim, mask, refractive, refractive_nonsmoothed,...
    absorption_coeff, raytogrid_indices_x, raytogrid_indices_y,...
    raytogrid_indices_z, raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix,...
    raytogrid_coeff_second_derivative_matrix, true, paraxial);

%if isempty(coeff_2)
    n_2 = 1; dn_2 = zeros(dim, 1); d2n_2 = 0;
%end


if ~auxiliary_adjoint_paraxial

    % normalise direction - the direction along the ray must be a unit vector
    % multiplied by the refractive index
    direction_2 = n_2/norm(direction_2) * direction_2;

    % calculate the position component of q2
    q2_x = 1/n_2 * direction_2;

    % store the current position
    pos_previous = pos;

    % update the current position of the ray based on the trapezoid rule
    % pos = pos_previous + ds/2 * (q1_x + q2_x);

    % get the update direction for the position
    direction_corrected = (q1_x+q2_x)/norm(q1_x+q2_x);

    % update the current position of the ray using a trapezoid rule
    pos = pos_previous + ds * direction_corrected;


else

    % get the current position
    pos = ray_positions(:, point_index);

end


if paraxial

    % calculate the first common matrix term, which is added for updating perturbations
    % to both space and direction vectors
    term_matrix = 1/n^2 * direction * dn';

    % calculate the perturbation position component of q1
    q1_perturb_x = 1/n * direction_perturbation;

    % initialise the direction perturbation component of q1
    % consider that the initial perturbation vector is [0,0]^T.
    q1_perturb_k = term_matrix' * direction_perturbation;

    % calculate the second term for updating the perturbation direction vector.
    % consider that the initial perturbation vector is [0,0]^T.
    direction_perturbation_2 = direction_perturbation + ds * q1_perturb_k;

    % calculate the second common term for updating the perturbation position vector
    % consider that the initial perturbation vector is [0,0]^T.
    pos_perturb_2 = ds * q1_perturb_x;

    % calculate the second common matrix term, which is added for updating perturbations
    % to both space and slowness vectors
    term_matrix_2 = 1/n_2^2 * direction_2  * dn_2';

    % calculate the perturbation position component of q2
    q2_perturb_x = - term_matrix_2 * pos_perturb_2 + 1/n_2 * direction_perturbation_2;

    % get the current position perturbation
    pos_perturb_previous = zeros(dim, 1);

    % update the position perturbation of the ray based on the trapezoid rule
    pos_perturb = ds/2 * (q1_perturb_x + q2_perturb_x);

    % update the second common matrix for q2
    term2_matrix_2 = -1/n_2 * (dn_2 * dn_2') + d2n_2;


    if do_direction_angle

        % compute the second paraxial system of equations which are solved for
        % computing the perturbation to the angle of the ray's directions
        % because of a perturbation to the initial position of the ray

        % calculate the second common matrix term, which is added for updating perturbations
        % to both space and direction vectors
        term2_matrix = -1/n * (dn * dn') + d2n;

        % calculate the perturbation position component of q1
        q1_perturb_x2 = - term_matrix * pos_perturb2 + ...
            1/n * direction_perturbation2;

        % update the direction perturbation component of q1
        q1_perturb_k2 = term2_matrix * pos_perturb2 + term_matrix' *...
            direction_perturbation2; 

        % calculate the second term for updating the perturbation direction
        % vector, considering that the initial perturbation to the
        % direction vector is zero
        direction_perturbation_22 = direction_perturbation2 + ds * q1_perturb_k2;

        % calculate the second term for updating the perturbation position vector
        pos_perturb_22 = pos_perturb2 + ds * q1_perturb_x2;

        % update the the direction perturbation component of q2
        q2_perturb_k2 =  term2_matrix_2 * pos_perturb_22 +...
            term_matrix_2' * direction_perturbation_22;

        % store the last perturbation to the ray's direction
        direction_perturbation_previous2 = direction_perturbation2;

        % update the direction perturbation based on the Trapezoid rule
        direction_perturbation2 = direction_perturbation_previous2 ...
            + ds/2 * (q1_perturb_k2 + q2_perturb_k2);


    end

end

% the ray must be terminated inside the detection surface
% if the ray is paraxial adjoint ray, the ray is stopped on the index associated with the last point on the ray
while   ((isfield(detec_geom, 'radius_circle') && norm(pos) - detec_geom.radius_circle <= 1e-10) ||...
        (isfield(detec_geom, 'radius_cylinder') && norm(pos(1:2)) - detec_geom.radius_cylinder <= 1e-10) ||...
        (isfield(detec_geom, 'line_coeff') && sgn * (detec_geom.line_coeff(1:2) * pos(1:2) - ...
        detec_geom.line_coeff(3)) <= 1e-10))...
        && all(pos > pos_grid_first + 2 * dx & pos < pos_grid_end - 3 * dx)  &&  (point_index < lastpoint_index)

    % increase the index of ray point by 1
    point_index = point_index + 1;


    if ~auxiliary_adjoint_paraxial

        % update the direction
        direction = direction + ds/2 * (dn + dn_2);

    else

        % get the direction
        direction = ray_positions(:, point_index) - pos;

    end



    if paraxial


        % update the the direction perturbation component of q2
        q2_perturb_k = term2_matrix_2 * pos_perturb_2 +...
            term_matrix_2' * direction_perturbation_2;

        % update the direction perturbation based on the Trapezoid rule
        direction_perturbation = direction_perturbation...
            + ds/2 * (q1_perturb_k + q2_perturb_k);

        if do_direction_angle

            % compute the second paraxial system of equations which are solved for
            % computing the perturbation to the angle of the ray's directions
            % because of a perturbation to the initial position of the ray

            % calculate the perturbation position component of q2
            q2_perturb_x2 = - term_matrix_2 * pos_perturb_22 + 1/n_2 * direction_perturbation_22;

            % update the position perturbation of the ray based on the trapezoid rule
            pos_perturb2 =  pos_perturb2 + ds/2 * (q1_perturb_x2 + q2_perturb_x2);

        end


    end

    % calculate the interpolation coefficients, interpolated parameters on the ray's point,
    % and their corresponding indices on the grid
    [~, coeff, n, dn, d2n, ns, as] = interpLocalBsplineFullWave(pos, xvec, yvec, zvec, ...
        pos_grid_first, dx, grid_size, dim, mask, refractive, refractive_nonsmoothed,...
        absorption_coeff, raytogrid_indices_x, raytogrid_indices_y,...
        raytogrid_indices_z, raytogrid_coeff_matrix, ...
        raytogrid_coeff_derivative_matrix, raytogrid_coeff_second_derivative_matrix,...
        true, paraxial);


    if isempty(coeff)
        n = 1; dn = zeros(dim, 1); d2n = 0; ns = 1;
        if do_absorption
            as = 0;
        end

    end

    % normalise direction - the direction along the ray must be a unit vector
    % multiplied by refractive index
    direction = n/norm(direction) * direction;

    % calculate the position component of q1
    q1_x = 1/n * direction;


    if calc_coeffs || auxiliary_ray

        if paraxial

            % add the current perturbed position vector to the associated matrix
            ray_positions_perturbation = [ray_positions_perturbation, pos_perturb];

            if do_direction_angle 

                % get the angle of the last (current) perturbation to the
                % ray's direction vector
                ray_angles_perturbation = [ray_angles_perturbation,...
                    1/(n*(point_index - 1)) * [-q1_x(2), q1_x(1)] * direction_perturbation2];

            end

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


    % calculate the second term for updating the direction
    direction_2 = direction + ds * dn;
    

    % calculate the second term for updating the position
    pos_2 = pos + ds * q1_x;


    % calculate the interpolation coefficients for the second point
    [~, coeff_2, n_2, dn_2, d2n_2, ~, ~] = interpLocalBsplineFullWave(pos_2, xvec, yvec, zvec, ...
        pos_grid_first, dx, grid_size, dim, mask, refractive, refractive_nonsmoothed,...
        absorption_coeff, raytogrid_indices_x, raytogrid_indices_y,...
        raytogrid_indices_z, raytogrid_coeff_matrix, ...
        raytogrid_coeff_derivative_matrix, raytogrid_coeff_second_derivative_matrix,...
        true, paraxial);

    if isempty(coeff_2)
        n_2 = 1; dn_2 = zeros(dim, 1); d2n_2 = 0;
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
        % pos = pos_previous + ds/2 * (q1_x + q2_x);

        % get the update direction for the position
        direction_corrected = (q1_x+q2_x)/norm(q1_x+q2_x);

        % update the current position of the ray using a trapezoid rule
        pos = pos_previous + ds * direction_corrected;

    else

        % get the current position
        pos = ray_positions(:, point_index);

    end

    if paraxial

        % calculate the first common matrix term, which is added for updating perturbations
        % to both space and slowness vectors
        term_matrix = 1/n^2 * direction * dn';

        % calculate the first common matrix term, which is added for updating perturbations
        % to both space and slowness vectors
        term2_matrix = -1/n * (dn * dn') + d2n;

        % calculate the perturbation position component of q1
        q1_perturb_x = - term_matrix * pos_perturb + 1/n * direction_perturbation;

        % update the direction perturbation component of q1
        q1_perturb_k = term2_matrix * pos_perturb +...
            term_matrix' * direction_perturbation;

        % calculate the second term for updating the perturbation direction vector.
        direction_perturbation_2 = direction_perturbation + ds * q1_perturb_k;

        % calculate the second term for updating the perturbation position vector
        pos_perturb_2 = pos_perturb + ds * q1_perturb_x;

        % calculate the second common scalar term, which is added for updating perturbations
        % to both space and slowness vectors
        term_matrix_2 = 1/n_2^2 * direction_2 * dn_2';

        % calculate the perturbation position component of q2
        q2_perturb_x = - term_matrix_2 * pos_perturb_2 + 1/n_2 * direction_perturbation_2;

        % get the current position perturbation
        pos_perturb_previous = pos_perturb;

        % update the position perturbation of the ray based on the trapezoid rule
        pos_perturb =  pos_perturb_previous + ds/2 * (q1_perturb_x + q2_perturb_x);

        % update the second common matrix for q2
        term2_matrix_2 = -1/n_2 * (dn_2 * dn_2') + d2n_2;

        if do_direction_angle

            % calculate the perturbation position component of q1
            q1_perturb_x2 = - term_matrix * pos_perturb2 + 1/n * direction_perturbation2;

            % update the direction perturbation component of q1
            q1_perturb_k2 = term2_matrix * pos_perturb2 +...
                term_matrix' * direction_perturbation2;

            % calculate the second term for updating the perturbation direction vector.
            direction_perturbation_22 = direction_perturbation2 + ds * q1_perturb_k2;

            % calculate the second term for updating the perturbation position vector
            pos_perturb_22 = pos_perturb2 + ds * q1_perturb_x2;

            % update the the direction perturbation component of q2
            q2_perturb_k2 = term2_matrix_2 * pos_perturb_22 +...
                term_matrix_2' * direction_perturbation_22;

            % store the last perturbation to the ray's direction
            direction_perturbation_previous2 = direction_perturbation2;

            % update the direction perturbation based on the Trapezoid rule
            direction_perturbation2 = direction_perturbation_previous2 ...
                + ds/2 * (q1_perturb_k2 + q2_perturb_k2);

        end

    end


end


%  the last point must be on the detection surface except for auxiliary
%  rays
if ~auxiliary_ray || paraxial

    if ~auxiliary_adjoint_paraxial

        % if the ray is not auxiliary, or the ray is an auxiliary ray which is
        % paraxial, the last point of the ray position must be on the detection
        % ring (surface)
        [pos, ds_final] = calcLineIntersect(pos_previous, direction_corrected,...
            detec_geom);

    else

        % compute the final ray spacing from the last two points along the
        % given main ray
        ds_final = norm(ray_positions(:, point_index)...
            - ray_positions(:, point_index-1));
    end

    if paraxial

        % update the last position perturbation of the ray based on the trapezoid rule
        pos_perturb = pos_perturb_previous + ds_final/2 * (q1_perturb_x + q2_perturb_x);

        if do_direction_angle

            % update the last angle of the perturbation to the direction of the ray
            direction_perturbation2 = direction_perturbation_previous2 + ...
                ds_final/2 * (q1_perturb_k2 + q2_perturb_k2);
        end

    end

else

    % keep the last position and ray spacing
    ds_final = ds;

end


if all (pos > pos_grid_first + 2 * dx  &  pos < pos_grid_end - 3 * dx)


    % calculate the interpolation coefficients, interpolated parameters on the ray's point,
    % and their corresponding indices on the grid
    [~, coeff, n, ~, ~, ns, as] = interpLocalBsplineFullWave(pos, xvec, yvec, zvec, ...
        pos_grid_first, dx, grid_size, dim, mask, refractive, refractive_nonsmoothed,...
        absorption_coeff, raytogrid_indices_x, raytogrid_indices_y,...
        raytogrid_indices_z, raytogrid_coeff_matrix, ...
        raytogrid_coeff_derivative_matrix, raytogrid_coeff_second_derivative_matrix,...
        false, paraxial);

    if isempty(coeff)
        ns = 1;

        if do_direction_angle
            n = 1;
        end
        
        if do_absorption
            as  = 0;
        end
    end


    if calc_coeffs || auxiliary_ray

        if paraxial

            % add the last perturbed position vector to the associated matrix
            ray_positions_perturbation = [ray_positions_perturbation, pos_perturb];

            if do_direction_angle
        
                % get the angle of the last (current) perturbation to the
                % ray's direction vector
                ray_angles_perturbation = [ray_angles_perturbation,...
                     1/(n*(point_index - 1)) * [-q1_x(2), q1_x(1)] * direction_perturbation2];

            end

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


if ~do_direction_angle

    % give an empty variable
    ray_angles_perturbation = [];

end


end