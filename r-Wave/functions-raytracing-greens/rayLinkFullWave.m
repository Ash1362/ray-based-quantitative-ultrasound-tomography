function [optimal_polar_initial_direction_allreceivers,...
    cartesian_position_endpoint_allreceivers, num_rays_allreceivers, ray_positions_allreceivers,...
    acoustic_length_allreceivers, absorption_allreceivers, rayspacing_allreceivers,...
    ray_positions_auxilary_left_allreceivers, ray_positions_auxilary_right_allreceivers,...
    adjoint_ray_positions_auxilary_left_allreceivers, adjoint_ray_positions_auxilary_right_allreceivers] =...
    rayLinkFullWave(ray_interp_coeffs, refractive, refractive_nonsmoothed, absorption_coeff,...
    cartesian_position_emitter, cartesian_position_receivers, polar_direction_allreceivers,...
    polar_initial_direction_allreceivers, xvec, yvec, zvec, pos_grid_first,...
    pos_grid_end, grid_spacing, ray_spacing, grid_size, dim, detec_geom,...
    mask, num_points, num_emitter, para)
%RSAYLINKFULLWAVE links the rays between a single emitter and all reeivers
%
% DESCRIPTION:
%           rayLinkFullWave links the rays emenated from a single emitter to
%           all receivers, and stores all the information corresponding to the
%           the linked rays. The ray linking is done iteratively often using a
%           given initial guess, or alternatively using straight rays

% USAGE:
%
%
% INPUTS:
%       ray_interp_coeffs    - a struct containing the specified variables for ray tracing
%                              using a 'Bilinear' intrerpolation, this includes the directional
%                              gradients of the refrective index
%                              distribution
%       'refractive_gradient_x' - discretised refractive index gradient along x
%       'refractive_gradient_y' - discretised refractive index gradient along y
%       'refractive_gradient_z' - discretised refractive index gradient along z
%                              using a'Bspline' interpolation, this
%                              includes the matrices for intertpolation:
%       'raytogrid_indices_x'   - x indices for B-spline interpolation
%       'raytogrid_indices_y'   - y indices for B-spline interpolation
%       'raytogrid_indices_z'   - z indices for B-spline interpolation
%       'raytogrid_coeff_matrix'  - matrix for calculating B-spline
%                                 interpolation coefficients of the field
%       'raytogrid_coeff_derivative_matrix' - matrix for calculating B-spline
%                                 interpolation coefficients of the directional
%                                 gradientsof the field
%       refractive             - the smoothed dispersive refractive index used for
%                                computing the trajectory of rays
%       refractive_nonsmoothed - the nonsmoothed refractive index matrix used
%                                for integration along the linked rays
%       absorption_coeff       - the nonsmoothed absorption coefficient matrix
%                                used for integration along the linked rays
%       cartesian_position_emitter - a dim x 1 cartesian position of the emitter
%       cartesian_position_receivers - a dim x num_receiver cartesian
%                                       position of all receivers
%       polar_direction_allreceivers - a (dim-1) x num_receiver vector of the
%                                      polar direction of straight rays
%                                      from emitter to all the receivers
%       polar_initial_direction_allreceivers - a (dim-1) x num_receiver matrix
%                                      of the polar initial direction of the
%                                      rays for all receivers
%       xvec                 - the x vector of grid points
%       yvec                 - the y vector of grid points
%       zvec                 - the z vector of grid points
%       pos_grid_first       - a dim x 1 Cartesian position of the first index of
%                              the grid
%       pos_grid_end         - a dim x 1 Cartesian position of the end index of
%                              the grid
%       grid_spcing          - a scalar representing the grid spacing [m]
%       ray_spacing          - a scalar representing the ray spacing [m]
%       grid_size            - the size of the grid
%       dim                  - the dimension of the medium
%       detec_radius         - the radius [m] of the excitaion/dection ring (or
%                              hemisphere)
%       mask                 - a binary mask for ray tracing
%       num_points           - a fixed maximum number for the maximum number of
%                              points on the auxiliary rays
%                              this should be larger than the maximum points between
%                              the rays linking emitter-receiver pairs.
%       num_emitter          - the number of emitters
%       para                 - a struct containing the fields:
%       'varepsilon'         - the stopping criterion for ray linking inverse problem
%       'max_iter'           - the maximum number of iterations for ray
%                              linking inverse problem
%       'raylink_method'     - method for ray linking
%       'interp_method'      - method for interpolation
% OPTIONAL INPUTS:
%
% OUTPUTS:
%       system_matrx_single_emitter - num_receiver x num_gridpoints  sparse
%                                    metrix containing the raytogrid interpolation
%                                    coefficients for all grid points
%       optimal_polar_initial_direction_allreceivers - dim-1 x num_receiver matrix
%                                               of the polar initial direction
%                                               for the linked rays. this can be
%                                               used as an initial guess for the
%                                               next UST iteration
%       cartesian_position_endpoint_allreceivers - dim x num_receiver matrix
%                                                  of Cartesian position of
%                                                  the end point of the ray
%                                                  for all receivers
%       num_rays_allreceivers               - num_receiver x 1 vector of the number of ray tracing
%                                              for ray linking between the
%                                              emitter and each of all receivers
%       ray_positions_allreceivers          - the (dim*num_receiver) x
%                                             num_points matrix containing
%                                             the Cartesian position of the
%                                             points on the linked rays
%       acoustic_length_allreceivers        - the num_receiver x num_points
%                                             matrix containing the
%                                             accumulated acoustic length
%                                             including dispersion effects
%                                             on the linked rays
%      absorption_allreceivers              - the num_receiver x num_points
%                                             matrix containing the
%                                             accumulated acoustic
%                                             absorption on the linked rays
%      rayspacing_allreceivers              - a num_receiver x 1 vector of
%                                             the spacing [m] between the
%                                             receivers and
%                                             the point on the correponding
%                                             linked ray just before the
%                                             receiver
%                                             (Those are different from
%                                             the fixed spacing for rays' points.)
%   ray_positions_auxilary_left_allreceivers - the Cartesian position of points on
%                                              left auxiliry ray for the
%                                              forward field
%   ray_positions_auxilary_right_allreceivers- the Cartesian position of points on
%                                              right auxiliry ray for the
%                                              forward field
%   adjoint_ray_positions_auxilary_left_allreceivers- the Cartesian position
%                                              of points on left auxiliry ray for the
%                                              adjoint field
%   adjoint_ray_positions_auxilary_right_allreceivers- the Cartesian position of points on
%                                              right auxiliry ray for the adjoint field
%   ray_directions_auxilary_allreceivers    - the perturbed direction of
%                                             the forward auxiliary rays after
%                                             perturbation to the initial
%                                             position
%   adjoint_ray_directions_auxilary_allreceivers - the perturbed direction of
%                                             the adjoint auxiliary rays after
%                                             perturbation to the initial
%                                             position
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 14.01.2023
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian.

% calculate the accumulated acosutic absorption, if the absorption
% coefficient map is nonzero
do_absorption = ~isscalar(absorption_coeff);

% number of receivers
num_receiver = size(polar_direction_allreceivers, 2);

% a binary vector, which is true when the distance of the emitter to
% receivers are sufficiently large
binary_distances = vecnorm(cartesian_position_receivers - cartesian_position_emitter) > ray_spacing;

%%=========================================================================
% ALLOCATE VECTORS/MATRICES FOR STORING THE INFORMATION FOR ALL RECEIVERS
%==========================================================================

% the polar initial direction of the linked rays
optimal_polar_initial_direction_allreceivers = nan(dim-1, num_receiver);

% the cartesian position of the end point of the linked rays
cartesian_position_endpoint_allreceivers = nan(dim, num_receiver);

% the number of traced rays for ray linking the emitter to each of the receivers
num_rays_allreceivers = nan(num_receiver, 1);

% the last spacing along the ray (the end point is the position of the
% receiver.)
rayspacing_allreceivers = nan(num_receiver, 1);

% position along the ray
ray_positions_allreceivers = nan(dim * num_receiver, num_points);

% the acoustic length along the ray
acoustic_length_allreceivers = nan(num_receiver, num_points);

% the accumulated acoustic absorption along the ray
if do_absorption
    absorption_allreceivers = nan(num_receiver, num_points);
else
    absorption_allreceivers = [];
end


if para.auxiliary_ray
    
    % matrices for the position of the auxiliary rays
    
    % forward left auxiliary ray
    ray_positions_auxilary_left_allreceivers = nan(dim * num_receiver, num_points);
    
    % adjoint left auxiliary ray
    adjoint_ray_positions_auxilary_left_allreceivers = nan(dim * num_receiver, num_points);
    
    if ~strcmp(para.auxiliary_method, 'paraxial')
        
        % forward right auxiliary ray
        ray_positions_auxilary_right_allreceivers = nan(dim * num_receiver, num_points);
        
        % adjoint right auxiliary ray
        adjoint_ray_positions_auxilary_right_allreceivers = nan(dim * num_receiver, num_points);

        
    else
 
        if para.do_direction_angle

        % the perturbation to the forward rays' directions in the polar coordinates because
        % of a perturbation to the initial positions
        ray_positions_auxilary_right_allreceivers  = nan(num_receiver, num_points);

        % the perturbation to the backward rays' directions in the polar coordinates because
        % of a perturbation to the initial positions
        adjoint_ray_positions_auxilary_right_allreceivers = nan(num_receiver, num_points);

        else

        % give empty variables    
        ray_positions_auxilary_right_allreceivers= [];
        adjoint_ray_positions_auxilary_right_allreceivers = [];
        
        end

        
    end
    
else
    
    % give empty variables for the positions of the auxiliary rays, if
    % auxiliary rays are not required
    ray_positions_auxilary_left_allreceivers = [];
    ray_positions_auxilary_right_allreceivers = [];
    adjoint_ray_positions_auxilary_left_allreceivers = [];
    adjoint_ray_positions_auxilary_right_allreceivers = [];
    

end

% define a handle function for tracing the ray, and obtain the information
% this function handle solves the forward problem of ray linking inverse
% problem
solve_ray = @(polar_initial_direction, polar_direction_receiver, calc_coeffs, auxiliary_ray)...
    calcRayLinkForwardFullWave(ray_interp_coeffs, refractive, refractive_nonsmoothed,...
    absorption_coeff, cartesian_position_emitter, polar_direction_receiver,...
    polar_initial_direction, xvec, yvec, zvec, pos_grid_first, pos_grid_end,...
    grid_spacing, ray_spacing, grid_size, dim, detec_geom, mask, calc_coeffs,...
    para.raylinking_method, para.interp_method, auxiliary_ray, [], []);

%% ========================================================================
% INITIAL GUESS
%==========================================================================

% choose initial guess for the left and right polar
% initial direction.

if strcmp(para.raylinking_method, 'Regula-Falsi')
    
    % By a ray initialisation using a 'Local' approach, the left and right
    % polar intial directions (initial interval) are chosen the same for
    % all receivers, i.e. approximately -pi/2, +pi/2
    % with respect to a geomterical vector from the emitter to the
    % centre of the detection circle. Both initial angles are chosen a bit
    % smaller in order to ensure the initial directions of the rays send the rayd
    % inside the detection circle, approximately tangent to the periphery of the circle
    polar_initial_direction_left  =  - pi/2 * (1 - (1/num_receiver) + para.varepsilon);
    polar_initial_direction_right =  + pi/2 * (1 - (1/num_receiver) + para.varepsilon);
    
    % calculate the polar direction of unit geometrical
    % vectors from emitter to the last point of the rays
    % solved by the left and right polar initial directions
    % because these two rays rae tangent to the circle,
    % they are intercepted by the the detection surface
    % very soon, and thus they are very short, and give
    % a maximal range for the left and right initial
    % directions
    [~, ~, polar_direction_endpoint_left] = feval(solve_ray,...
        polar_initial_direction_left, [], false, false);
    [~, ~, polar_direction_endpoint_right] = feval(solve_ray,...
        polar_initial_direction_right, [], false, false);
    
end


if strcmp(para.raylinking_method, 'Secant')

    % get the large bound as tangent directions to the detection ring
    % in the polar coordinate
    % bound_large = wrapToPi(cart2pol(-cartesian_position_emitter(1),...
    %    -cartesian_position_emitter(2)) + [-pi/2, pi/2]);

    % get the emitter position in the polar coordinate
    [emitter_angle, emitter_radius] = cart2pol(cartesian_position_emitter(1),...
        cartesian_position_emitter(2));

    % get the angle spacing in terms of the ray spacing along the periphery
    % of the detection ring
    ray_spacing_angle = ray_spacing/(1.2 *emitter_radius);
    % ray_spacing_angle = pi / num_receiver;
    
    % get the cartesian position of the end point of a vector representing
    % the lower bound for the initial direction of the ray
    [pos_low_x, pos_low_y] = pol2cart(emitter_angle + ray_spacing_angle,...
        emitter_radius);

    % get the cartesian position of the end point of a vector representing
    % the upper bound for the initial direction of the ray
    [pos_up_x, pos_up_y] = pol2cart(emitter_angle - ray_spacing_angle,...
        emitter_radius);

    % get the bounds for the initial angle of the rays initiaised on the 
    % emitter position
    bound_large = wrapToPi([cart2pol(pos_low_x - cartesian_position_emitter(1),...
        pos_low_y - cartesian_position_emitter(2)),...
        cart2pol(pos_up_x - cartesian_position_emitter(1),...
        pos_up_y - cartesian_position_emitter(2))]);

end

% solve the ray linking inverse problem for each individual emitter-receiver pair
for ind_receiver = 1 : num_receiver
  
    %disp( ['Number of receiver:'  num2str(ind_receiver)] )

    if binary_distances(ind_receiver)
        
        
        switch para.raylinking_method
            
            case 'Regula-Falsi'
                
                polar_direction_initial_guess(1) = polar_initial_direction_left;
                polar_direction_initial_guess(2) = polar_initial_direction_right;
                res_initial_guess(1) = polar_direction_endpoint_left...
                    - polar_direction_allreceivers(ind_receiver);
                res_initial_guess(2) = polar_direction_endpoint_right...
                    - polar_direction_allreceivers(ind_receiver);
                
                % calculate the path of rays using the Regula-falsi approach
                [polar_initial_direction, cartesian_position_endpoint, num_rays,...
                    ray_positions, ray_acoustic_length, ray_absorption, rayspacing_receiver] =...
                    solveRegulaFalsiFullWave(solve_ray, polar_direction_allreceivers(:, ind_receiver),...
                    polar_direction_initial_guess, res_initial_guess, para.varepsilon, para.max_iter);
                
                
            case 'Secant'
                [polar_initial_direction, cartesian_position_endpoint, num_rays,...
                    ray_positions, ray_acoustic_length, ray_absorption, rayspacing_receiver] =...
                    solveSecantFullWave(solve_ray, polar_direction_allreceivers(ind_receiver),...
                    polar_initial_direction_allreceivers(ind_receiver), bound_large,...
                    para.varepsilon, para.max_iter);
                
            case 'Newton'
                
                % give an error, because taking Newton's approach for ray
                % linking is unnecessarily costly
                error(['The Newton approach for ray linking is deprecated.'])
                
                % ray linking using Newton's method
                [polar_initial_direction, cartesian_position_endpoint, num_rays]=...
                    solveNewton(solve_ray, polar_direction_allreceivers(:, ind_receiver),...
                    polar_initial_direction_allreceivers(:, ind_receiver), para.varepsilon, para.max_iter);
                
            case 'Quasi-Newton'
                
                % give an error
                error(['For 3D case, only the time-of-flight-based approach'...
                    'is included for this code version.'])
                
        end


    %else
        
        % the traced ray is a straight line travelling through
        % the water and directly intercepted by the receiver.
        % Therefore, the ray linking inverse problem is not solved
        % for this emitter-receiver pair, and because the interpolation
        % coefficients are the same as those in the system matrix for
        % the homogeneous water, the coefficients will be ignored in the system
        % matrix, and are not stored.
        %if ~strcmp(para.raylinking_method, 'Regula-Falsi')
        %    polar_initial_direction = polar_initial_direction_allreceivers(:, ind_receiver);
        %end
        
        % give nan to the parameters of ray and associated Green's function
     %   cartesian_position_endpoint = nan;
     %   num_rays = 0;
     %   ray_positions = nan;
     %   ray_acoustic_length = nan;
     %   ray_absorption = nan;
     %   polar_initial_direction = nan;
        
    %end
    
    
    % fill the matrix for updated polar initial direction for each receiver
    % using 'Regula-Falsi' method, the optimal polar
    % initial direction will not be used as the initial
    % guess for the next UST iteration
    if ~strcmp(para.raylinking_method, 'Regula-Falsi')
        optimal_polar_initial_direction_allreceivers(:, ind_receiver) =...
            polar_initial_direction;
    end
    
    % get the ray positions for all the receivers
    ray_positions_allreceivers(ind_receiver * dim - (dim - 1): ind_receiver * dim,...
        1:size(ray_positions, 2) ) = ray_positions;
    
    % get the accumulated acoustic length for all the receivers
    acoustic_length_allreceivers(ind_receiver, 1:size(ray_acoustic_length, 2) )...
        = ray_acoustic_length;
    
    if do_absorption
        
        %  get the accumulated acoustic absorption for all the receivers
        absorption_allreceivers(ind_receiver, 1:size(ray_absorption, 2) )...
            = ray_absorption;
    end
    
    if para.auxiliary_ray
        
     
        %%=================================================================
        % CALCULATE THE FORWARD AUXILIARY RAYS
        %==================================================================
        
        % get the handle function for the forward auxiliary rays
        solve_auxiliary_ray = @(polar_initial_direction, cartesian_position_start,...
            calc_coeffs, auxiliary_ray, auxiliary_method)...
            calcRayLinkForwardFullWave(ray_interp_coeffs, refractive, refractive_nonsmoothed,...
            absorption_coeff, cartesian_position_start, [], polar_initial_direction,...
            xvec, yvec, zvec, pos_grid_first, pos_grid_end, grid_spacing,...
            ray_spacing, grid_size, dim, detec_geom, mask, calc_coeffs,...
            para.raylinking_method, para.interp_method, auxiliary_ray,...
            auxiliary_method, para.do_direction_angle);
        
        
        
        switch para.auxiliary_method
            case 'angle_perturbation'
                
                % get the perturbed initial angle for the left
                % and right auxiliary rays
                angle_left = polar_initial_direction - para.reference_angle;
                angle_right = polar_initial_direction + para.reference_angle;
                
                % calculate the position of the right auxiliary ray
                [~, ~, ~, ray_positions_auxilary_left, ~, ~, ~] = feval(solve_auxiliary_ray,...
                    angle_left, cartesian_position_emitter, false, true,...
                    para.auxiliary_method);
                
                % calculate the position of the right auxiliary ray
                [~, ~, ~, ray_positions_auxilary_right, ~, ~, ~] = feval(solve_auxiliary_ray,...
                    angle_right, cartesian_position_emitter, false, true,...
                    para.auxiliary_method);
                
                
            case 'paraxial'
                
                switch dim
                    
                    case 2
                        
                        % calculate the perturbed position for the first auxiliary
                        % ray
                        [~, ~, ~, ray_positions_auxilary_left, ~, ~, ~, ray_angles] = feval(solve_auxiliary_ray,...
                            polar_initial_direction, cartesian_position_emitter, false, true, 'paraxial');
                        
                    case 3
                        

                        error('Not implemented yet.')
                        % calculate the perturbed position for the first auxiliary
                        % ray
                        [~, ~, ~, ray_positions_auxilary_left, ~, ~, ~] = feval(solve_auxiliary_ray,...
                            polar_initial_direction, cartesian_position_emitter, false, true, 'paraxial1');
                        
                        % calculate the perturbed position for the second auxiliary
                        % ray
                        [~, ~, ~, ray_positions_auxilary_right, ~, ~, ~] = feval(solve_auxiliary_ray,...
                            polar_initial_direction, cartesian_position_emitter, false, true, 'paraxial2');
                        
                        % get the cross product of rays' pertrurbed positions
                        ray_positions_auxilary_left = cross(ray_positions_auxilary_left,...
                            ray_positions_auxilary_right, 1);
                        
                end
                
        end
        
        %%=====================================================================
        % CALCULATE THE ADJOINT AUXILIARY RAYS
        %======================================================================
        % get the initial direction of the adjoint ray in the polar
        % coordinate
        switch para.raylinking_method
            case 'Regula-Falsi'
                
                adjoint_polar_initial_direction = calcDirectionalAngle(...
                    [-cartesian_position_receivers(:, ind_receiver);  0],...
                    [ray_positions(:, end-1) - ray_positions(:, end); 0]);
            case 'Secant'
                
                [adjoint_polar_initial_direction, ~] = cart2pol(...
                    ray_positions(1, end-1) - ray_positions(1, end),...
                    ray_positions(2, end-1) - ray_positions(2, end));
            case 'Quasi-Newton' 
                    
                % give an error
                error(['For 3D case, only the time-of-flight-based approach'...
                    'is included for this code version.'])
                
        end
        
        
        % get the initial angles for the adjoint auxiliary rays
        switch para.auxiliary_method
            case 'angle_perturbation'
                
                % polar direction for the left auxiliary ray
                angle_left = adjoint_polar_initial_direction - para.reference_angle;
                % polar direction for the right auxiliary ray
                angle_right = adjoint_polar_initial_direction + para.reference_angle;
                
                % calculate the position of the left auxiliary adjoint ray.
                [~, ~, ~, adjoint_ray_positions_auxilary_left, ~, ~, ~] = feval(solve_auxiliary_ray,...
                    angle_left, cartesian_position_receivers(:, ind_receiver), false, true,...
                    para.auxiliary_method);
                
                % calculate the position of the right auxiliary adjoint ray.
                [~, ~, ~, adjoint_ray_positions_auxilary_right, ~, ~, ~] = feval(solve_auxiliary_ray,...
                    angle_right, cartesian_position_receivers(:, ind_receiver), false, true,...
                    para.auxiliary_method);
                
            case 'paraxial'
                
                % polar direction for the left auxiliary ray
                angle_left = adjoint_polar_initial_direction;
                
                % get the amplitude for the perturbation to the
                % initial position
                [ray_positions_adjoint] = calcRayAdjointParaxial(ray_positions,...
                    ray_spacing, rayspacing_receiver);
                
                switch dim
                    case 2
                        
                        % calculate the perturbed position for the first auxiliary
                        % adjoint ray
                        [~, ~, ~, adjoint_ray_positions_auxilary_left, ~, ~, ~, adjoint_ray_angles] = ...
                            feval(solve_auxiliary_ray, angle_left, ray_positions_adjoint,...
                            false, true, 'paraxial');
                        
                    case 3
                        

                        error('Not implemented yet.')

                        % calculate the perturbed position for the first auxiliary
                        % adjoint ray
                        [~, ~, ~, adjoint_ray_positions_auxilary_left, ~, ~, ~] = ...
                            feval(solve_auxiliary_ray, angle_left, ray_positions_adjoint,...
                            false, true, 'paraxial1');
                        
                        % calculate the perturbed position for the second auxiliary
                        % adjoint ray
                        [~, ~, ~, adjoint_ray_positions_auxilary_right, ~, ~, ~] = ...
                            feval(solve_auxiliary_ray, angle_left, ray_positions_adjoint,...
                            false, true, 'paraxial2');
                        
                        % get the cross product of adjoint rays' pertrurbed
                        % positions.
                        adjoint_ray_positions_auxilary_left = cross(adjoint_ray_positions_auxilary_left,...
                            adjoint_ray_positions_auxilary_right, 1);
                end
                
        end
        
        
        %%=====================================================================
        % FILL THE MATRICES FOR THE AUXILIARY RAYS FOR THE CURRENT RECEIVER
        %======================================================================
        % get the starting index for filling the matrix for the position of the
        % auxiliary rays
        receiver_index_position = ind_receiver * dim - (dim-1);
        
        % get the ray positions for the forward auxiliary rays
        % left auxiliary ray
        ray_positions_auxilary_left_allreceivers(receiver_index_position:ind_receiver * dim,...
            1:size(ray_positions_auxilary_left, 2) ) = ray_positions_auxilary_left;
        
        % get the ray positions for the adjoint auxiliary rays
        % left auxiliary ray
        adjoint_ray_positions_auxilary_left_allreceivers(receiver_index_position:ind_receiver * dim,...
            1:size(adjoint_ray_positions_auxilary_left, 2) ) = adjoint_ray_positions_auxilary_left;
        
        
        if strcmp(para.auxiliary_method, 'paraxial') 

            if para.do_direction_angle

            % get the perturbation to the forward ray because of a perturbation to the
            % initial position for the current receiver
            ray_positions_auxilary_right_allreceivers(ind_receiver,...
                1:size(ray_angles, 2)) = 2 * pi / num_emitter * ray_angles;

            % get the perturbation to the backward ray because of a perturbation to the
            % initial position for the current receiver
            adjoint_ray_positions_auxilary_right_allreceivers(ind_receiver,...
                1:size(adjoint_ray_angles, 2)) = 2 * pi/ num_receiver * adjoint_ray_angles;

            end

        else
            
            % right auxiliary ray
            ray_positions_auxilary_right_allreceivers(receiver_index_position:ind_receiver * dim,...
                1:size(ray_positions_auxilary_right, 2) ) = ray_positions_auxilary_right;
            
            % right auxiliary ray
            adjoint_ray_positions_auxilary_right_allreceivers(receiver_index_position:ind_receiver * dim,...
                1:size(adjoint_ray_positions_auxilary_right, 2) ) = adjoint_ray_positions_auxilary_right;
            
        end
        
    end
    
    % fill the matrix's row for cartesian position of the end point of the linked ray for each receiver
    cartesian_position_endpoint_allreceivers(:, ind_receiver) = cartesian_position_endpoint;
    
    % fill the matrix's row for the number of solved rays for ray
    % linking for each receiver
    num_rays_allreceivers(ind_receiver) = num_rays;
    
    % get the last ray spacings for all receivers
    rayspacing_allreceivers(ind_receiver) = rayspacing_receiver;
    
    end
end

% using 'Regula-Falsi' method, the optimal polar
% initial direction will not be used as the initial
% guess for the next UST iteration
if strcmp(para.raylinking_method, 'Regula-Falsi')
    optimal_polar_initial_direction_allreceivers = [];
end



end