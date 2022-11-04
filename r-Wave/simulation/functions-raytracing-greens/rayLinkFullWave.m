function [optimal_polar_initial_direction_allreceivers,...
    cartesian_position_endpoint_allreceivers, num_rays_allreceivers, ray_positions_allreceivers,...
    acoustic_length_allreceivers, absorption_allreceivers, rayspacing_allreceivers, ray_positions_auxilary_left_allreceivers,...
    ray_positions_auxilary_right_allreceivers, adjoint_ray_positions_auxilary_left_allreceivers,...
    adjoint_ray_positions_auxilary_right_allreceivers, ray_directions_auxilary_allreceivers,...
    adjoint_ray_directions_auxilary_allreceivers] = rayLinkFullWave(ray_fields_params,...
    refractive, refractive_nonsmoothed, absorption_coeff, cartesian_position_emitter,...
    cartesian_position_receivers, polar_direction_allreceivers,...
    polar_initial_direction_allreceivers, xvec, yvec, zvec, pos_grid_first,...
    pos_grid_end, grid_spacing, ray_spacing, grid_size, dim, detec_radius,...
    mask, num_points, para)
%RSAYLINKFULLWAVE links the rays between a single emitter and all reeivers
%
% DESCRIPTION:
%           rayLinkFullWave links the rays emenated from a single emitter to
%           all receivers, and stores all the information correspinding to the
%           the linked rays. The ray linking is done iteratively often using a
%           given initial guess, or alternatively using straight rays

% USAGE:
%
%
% INPUTS:
%       ray_fields_params      - a struct containing the specified variables for ray tracing
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
%       refractive_nonsmoothed - the nonsmoothed refractive index used
%                                for integration of the parameters along the
%                                linked rays
%       absorption_coeff        - the absorption coefficient matrix
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
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian and Ben Cox

% calculate the accumulated acosutic absorption, if the absorption
% coefficient map is nonzero
do_absorption = ~isscalar(absorption_coeff);


% number of receivers
num_receiver = size(polar_direction_allreceivers, 2);


% a binary vector, which is true when the distances of the emeitter to
% receivers are sufficiently large
binary_distances = vecnorm(cartesian_position_receivers - cartesian_position_emitter) >  0; % para.minimum_distance;



switch para.interp_method
    case 'Bilinear'
        
        % get the directional gradient fields
        refractive_gradient_x = ray_fields_params.refractive_gradient_x;
        refractive_gradient_y = ray_fields_params.refractive_gradient_y;
        refractive_gradient_z = ray_fields_params.refractive_gradient_z;
        
        % allocate empty variables for parameters for spline interpolation
        raytogrid_indices_x = [];
        raytogrid_indices_y = [];
        raytogrid_indices_z = [];
        raytogrid_coeff_matrix = [];
        raytogrid_coeff_derivative_matrix = [];
        
    case {'Bspline'}
        
        % allocate empty variables for the directional gradient fields
        % refractive_gradient_x = [];
        % refractive_gradient_y = [];
        % refractive_gradient_z = [];
        
        % get parameters for spline interpolation
        raytogrid_indices_x = ray_fields_params.raytogrid_indices_x;
        raytogrid_indices_y = ray_fields_params.raytogrid_indices_y;
        raytogrid_indices_z = ray_fields_params.raytogrid_indices_z;
        raytogrid_coeff_matrix = ray_fields_params.raytogrid_coeff_matrix;
        raytogrid_coeff_derivative_matrix = ray_fields_params.raytogrid_coeff_derivative_matrix;
        raytogrid_coeff_second_derivative_matrix = ray_fields_params.raytogrid_coeff_second_derivative_matrix ;
end


% the indices of rows of nonzero elements of a sparse system matrix, and
% represents the indices of receivers
% ind_receivers = [];

% the indices of columns of nonzeros elements of a sparse system matrix,
% and represents the indices of grid points in a Matlab column-wise order
% ind_gridpoints = [];

% the values of nonzeros elements of the sparse system matrix, and includes the
% coefficient for interpolation between the rays
% linking the emitter to receivers (row indices)and grid points (column indices)
% val_coeffs = [];

% the polar initial direction of the linked rays
optimal_polar_initial_direction_allreceivers = zeros(dim-1, num_receiver);

% allocate matrices for storing the information from ray linking
% the cartesian position of the end point of the linked rays
cartesian_position_endpoint_allreceivers = zeros(dim, num_receiver);

% allocate a vector for the number of traced rays for ray linking
% between the vector and each of the receivers
num_rays_allreceivers = zeros(num_receiver, 1);

% % allocate a vector for the last spacing along the ray with the end point
% position of the receiver
rayspacing_allreceivers = zeros(num_receiver, 1);

% allocate vector for poistion of time delays of ray points
ray_positions_allreceivers = nan(dim*num_receiver, num_points);
acoustic_length_allreceivers = nan(num_receiver, num_points);
if do_absorption
    absorption_allreceivers = nan(num_receiver, num_points);
else
    absorption_allreceivers = [];
end



if para.auxiliary_ray
    
    % allocate matrices for the position of the auxiliary rays
    % forward left auxiliary ray
    ray_positions_auxilary_left_allreceivers = nan(dim * num_receiver, num_points);
    
    % adjoint left auxiliary ray
    adjoint_ray_positions_auxilary_left_allreceivers = nan(dim * num_receiver, num_points);
    
    
    
    if ~ strcmp(para.auxiliary_method, 'paraxial')
        
        % forward right auxiliary ray
        ray_positions_auxilary_right_allreceivers = nan(dim * num_receiver, num_points);
        % adjoint right auxiliary ray
        adjoint_ray_positions_auxilary_right_allreceivers = nan(dim * num_receiver, num_points);
        
        % give empty variables for the positions of the auxiliary rays, if
        % d direction/d x(s_0)$ is not required
        ray_directions_auxilary_allreceivers = [];
        adjoint_ray_directions_auxilary_allreceivers = [];
        
        
    else
        
        % allocate empty variables to the right auxiliary ray
        ray_positions_auxilary_right_allreceivers = [];
        adjoint_ray_positions_auxilary_right_allreceivers = [];
        
        if  para.do_perturb_initial_position
            
            % allocate matrices for the polar directions of the auxiliary rays
            % for computing % d direction/d x(s_0)$
            % forward auxiliary ray
            ray_directions_auxilary_allreceivers = nan((dim)*num_receiver, num_points);
            
            % adjoint auxiliary ray
            adjoint_ray_directions_auxilary_allreceivers = nan((dim)*num_receiver, num_points);
            
            
        else
            
            % give empty variables for the positions of the auxiliary rays, if
            % d direction/d x(s_0)$ is not required
            ray_directions_auxilary_allreceivers = [];
            adjoint_ray_directions_auxilary_allreceivers = [];
            
        end
        
        
    end
    
else
    
    % give empty variables for the positions of the auxiliary rays, if
    % auxiliary rays are not required
    ray_positions_auxilary_left_allreceivers = [];
    ray_positions_auxilary_right_allreceivers = [];
    adjoint_ray_positions_auxilary_left_allreceivers = [];
    adjoint_ray_positions_auxilary_right_allreceivers = [];
    
    % give empty variables for the directions of the auxiliary rays, if
    % auxiliary rays are not required
    ray_directions_auxilary_allreceivers = [];
    adjoint_ray_directions_auxilary_allreceivers = [];
    
    
end

% define a handle function for tracing the ray, and obtain the information
% this function handle solves the forward problem of ray linking inverse
% problem
solve_ray = @(polar_initial_direction, polar_direction_receiver, calc_coeffs, auxiliary_ray)...
    calcRayParametersBsplineFullWaveModified(...
    refractive, cartesian_position_emitter, polar_direction_receiver, polar_initial_direction, [], ...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, grid_spacing, ray_spacing, ...
    grid_size, dim, detec_radius, mask, raytogrid_indices_x, raytogrid_indices_y,...
    raytogrid_indices_z, raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix, ...
    raytogrid_coeff_second_derivative_matrix, refractive_nonsmoothed, absorption_coeff, calc_coeffs,...
    para.raylinking_method, para.interp_method, auxiliary_ray, []);




%% ========================================================================
% INITIAL GUESS
%==========================================================================

% choose initial guess for the left and right polar
% initial direction.
switch para.raylinking_initialisation
    case 'Local'
        
        switch para.raylinking_method
            
            case 'Regula-Falsi'
                
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
                [~, ~, polar_direction_endpoint_left] = feval(solve_ray, polar_initial_direction_left, [], false, false);
                [~, ~, polar_direction_endpoint_right] = feval(solve_ray, polar_initial_direction_right, [], false, false);
                
            otherwise
                
                % By a ray initialisation using a 'Local' approach,
                % the
                polar_direction_initial_guess_allreceivers = [];
                res_initial_guess_allreceivers = [];
                
        end
        
    case 'Global'
        
        
        Error('Not implemented yet!')
        
        % Using 'Regula Falsi' approach for ray linking, together with
        % a ray initialisation using a 'Global' approach,
        % the left and right polar intial directions (initial interval)
        % are chosen individually for each receiver.
        % A number of rays with even angles are simultaneously solved. For each receiver,
        % a pair of rays with endpoints closest to the receiver,
        % and encompass the receiver are chosen as the
        % initial guess for the the left and right polar
        % initial directions. This pair will give opposite
        % signs for the residual corresponding to the polar direction of
        % the end point of the rays, and forms a very tight (optimal) initial interval
        % for applying the 'Regula Falsi' method.
        
        % Using other methods, the initial guess for the direction
        % of the ray will be the direction of the ray for the optimal ray
        % calculated from the last UST iteration.
        % Only for the first iteration, the initial guess for the initial direction
        % of the ray is chosen the direction of a straight line segment linking emitter to
        % the receiver
        
        [polar_direction_initial_guess_allreceivers, res_initial_guess_allreceivers] =...
            initialiseAnglesFullWave(solve_ray, polar_direction_allreceivers,...
            dim, 1, para.varepsilon, para.raylinking_method, false);
end






% solve the ray linking inverse problem for each individual emitter-receiver pair
for ind_receiver = 1 : num_receiver
    
    
    % avoid ray linking for a sufficiently close emitter-receiver pair,
    % because the rays only travel through the water for
    % this pair
    if binary_distances(ind_receiver)
        
        
        switch para.raylinking_method
            
            case 'Regula-Falsi'
                
                % the initial left and right initial directions
                switch para.raylinking_initialisation
                    case 'Local'
                        
                        polar_direction_initial_guess(1) = polar_initial_direction_left;
                        polar_direction_initial_guess(2) = polar_initial_direction_right;
                        res_initial_guess(1) = polar_direction_endpoint_left...
                            - polar_direction_allreceivers(ind_receiver);
                        res_initial_guess(2) = polar_direction_endpoint_right...
                            - polar_direction_allreceivers(ind_receiver);
                    case 'Global'
                        
                        polar_direction_initial_guess(1) = polar_direction_initial_guess_allreceivers(1, ind_receiver);
                        polar_direction_initial_guess(2) = polar_direction_initial_guess_allreceivers(2, ind_receiver);
                        res_initial_guess(1) = res_initial_guess_allreceivers(1, ind_receiver);
                        res_initial_guess(2) = res_initial_guess_allreceivers(2, ind_receiver);
                end
                
                % calculate the path of rays using the Regula-falsi approach
                [polar_initial_direction, cartesian_position_endpoint, num_rays,...
                    ray_positions, ray_acoustic_length, ray_absorption, rayspacing_receiver] = solveRegulaFalsiFullWave(solve_ray,...
                    polar_direction_allreceivers(:, ind_receiver), polar_direction_initial_guess,...
                    res_initial_guess, para.varepsilon, para.max_iter);
                
                
            case 'Secant'
                [polar_initial_direction, cartesian_position_endpoint, num_rays,...
                    ray_positions, ray_acoustic_length, ray_absorption, rayspacing_receiver]=...
                    solveSecantFullWave(solve_ray, polar_direction_allreceivers(ind_receiver),...
                    polar_initial_direction_allreceivers(ind_receiver), para.varepsilon, para.max_iter);
                
            case 'Newton'
                
                % give an error, because taking Newton's approach for ray
                % linking is unnecessarily costly
                Error(['Newtons approach is not used for ray linking for the Greens approch'...
                    'because of its high computational cost.'])
                
                % ray linking using Newton's method
                [polar_initial_direction, cartesian_position_endpoint, num_rays]=...
                    solveNewton(solve_ray, polar_direction_allreceivers(:, ind_receiver),...
                    polar_initial_direction_allreceivers(:, ind_receiver), para.varepsilon, para.max_iter);
                
            case 'Quasi-Newton'
                
                % ray linking using Quasi-Newton method
                link_args = {'Method', 'Good-Broyden', 'initial_derivative', 'finite-difference',...
                    'smooth', true};
                
                [polar_initial_direction, cartesian_position_endpoint, num_rays,...
                    ray_positions, ray_acoustic_length, ray_absorption, rayspacing_receiver]=...
                    solveQuasiNewtonFullwave(solve_ray, polar_direction_allreceivers(:, ind_receiver),...
                    polar_initial_direction_allreceivers(:, ind_receiver), para.varepsilon, para.max_iter, [], link_args{:});
        end
        
    else
        
        % the traced ray is a straight line travelling through
        % the water and directly intercepted by the receiver.
        % Therefore, the ray linking inverse problem is not solved
        % for this emitter-receiver pair, and because the interpolation
        % coefficients are the same as those in the system matrix for
        % the homogeneous water, the coefficients will be ignored in the system
        % matrix, and are not stored.
        if ~strcmp(para.raylinking_method, 'Regula-Falsi')
            polar_initial_direction = polar_initial_direction_allreceivers(:, ind_receiver);
        end
        
        cartesian_position_endpoint = nan;
        num_rays = 0;
        ray_positions = nan;
        ray_acoustic_length = nan;
        ray_absorption = nan;
        
        
    end
    
    
    
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
            calc_coeffs, auxiliary_ray, rotation_matrix)...
            calcRayParametersBsplineFullWaveModified(...
            refractive, cartesian_position_start, [], polar_initial_direction,...
            rotation_matrix, xvec, yvec, zvec, pos_grid_first,...
            pos_grid_end, grid_spacing, ray_spacing, grid_size, dim, detec_radius,...
            mask, raytogrid_indices_x, raytogrid_indices_y, raytogrid_indices_z,...
            raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix,...
            raytogrid_coeff_second_derivative_matrix, refractive_nonsmoothed,...
            absorption_coeff, calc_coeffs, para.raylinking_method, para.interp_method,...
            auxiliary_ray, para.auxiliary_method);
        
        
        
        switch para.auxiliary_method
            case 'angle_perturbation'
                
                % allocate empty variable for the rotation matrices
                rotation_matrix_left = [];
                rotation_matrix_right = [];
                
                
                % get the perturbed initial angle for the left
                % and right auxiliary rays
                angle_left = polar_initial_direction - para.reference_angle;
                angle_right = polar_initial_direction + para.reference_angle;
                
                
                % calculate the position of the right auxiliary ray
                % for paraxial rays, perturbed positions, (and perturbed directions, if
                % requested)
                [~, ~, ~, ray_positions_auxilary_right, ~, ~, ~] = feval(solve_auxiliary_ray,...
                    angle_right, cartesian_position_emitter, false, true,...
                    rotation_matrix_right);
                
                
                
            case 'paraxial'
                
                switch dim
                    
                    case 2
                        
                        % get the two initial orthogonal perturbation vectors
                        % by rotating the initial direction of the main ray by
                        % -pi/4 and pi/4 Radians
                        
                        % get rotation angle in Radians
                        angle_rotation = pi/2;
                        
                        % get the rotation matrix for the left auxiliary ray
                        % This matrix, which enforces pi/2 radians rotation,
                        % is used for ensuring that the initial direction
                        % perturbation is normal to the initial direction.
                        rotation_matrix_left = [cos(-angle_rotation), -sin(-angle_rotation);...
                            sin(-angle_rotation), cos(-angle_rotation)];
                        
                    case 3
                        
                        % get the rotation matrix for the left auxiliary ray
                        rotation_matrix_left = [0, 0, 0;...
                            0, 0, 1;...
                            0, -1, 0];
                        
                        % get the rotation matrix for the right auxiliary ray
                        rotation_matrix_right = [0, 0, 1;...
                            0, 0, 0,;...
                            -1, 0, 0];
                        
                        
                end
                
                % get the initial polar direction of the ray before
                % rotation
                angle_left = polar_initial_direction;
                % angle_right = polar_initial_direction;
                
        end
        
        
        % calculate the position of the left auxiliary ray.
        % for paraxial rays, perturbed positions, (and perturbed directions, if
        % requested)
        [~, ~, ~, ray_positions_auxilary_left, ~, ~, ~] = feval(solve_auxiliary_ray,...
            angle_left, cartesian_position_emitter, false, true,...
            rotation_matrix_left);
        
        
        
        
        %%=====================================================================
        % CALCULATE THE ADJOINT AUXILIARY RAYS
        %======================================================================
        % get the initial direction of the adjoint ray in the polar
        % coordinate
        switch para.raylinking_method
            case 'Regula-Falsi'
                adjoint_polar_initial_direction = calcDirectionalAngle([-cartesian_position_receivers(:, ind_receiver);  0],...
                    [ray_positions(:, end-1) - ray_positions(:, end); 0]);
            case 'Secant'
                [adjoint_polar_initial_direction, ~] = cart2pol(ray_positions(1, end-1) - ray_positions(1, end),...
                    ray_positions(2, end-1) - ray_positions(2, end));
        end
        
        
        % get the initial angles for the adjoint auxiliary rays
        switch para.auxiliary_method
            case 'angle_perturbation'
                
                % polar direction for the left auxiliary ray
                angle_left = adjoint_polar_initial_direction - para.reference_angle;
                % polar direction for the right auxiliary ray
                angle_right = adjoint_polar_initial_direction + para.reference_angle;
                
                
                % calculate the position of the right auxiliary adjoint ray.
                [~, ~, ~, adjoint_ray_positions_auxilary_right, ~, ~, ~] = feval(solve_auxiliary_ray, angle_right,...
                    cartesian_position_receivers(:, ind_receiver), false, true,...
                    rotation_matrix_right);
                
                
                % calculate the position of the left auxiliary adjoint ray.
                [~, ~, ~, adjoint_ray_positions_auxilary_left, ~, ~, ~] = feval(solve_auxiliary_ray,...
                    angle_left, cartesian_position_receivers(:, ind_receiver),...
                    false, true, rotation_matrix_left);
                
                
            case 'paraxial'
                
                % polar direction for the left auxiliary ray
                angle_left = adjoint_polar_initial_direction;
                
                % get the amplitude for the perturbation to the
                % initial position
                % rotation_matrix_left.amplitude = para.perturbation_position_amplitude;
                [ray_positions_adjoint] = calcRayAdjointParaxial(ray_positions, ray_spacing, rayspacing_receiver);
                
                % calculate the position of the left auxiliary adjoint ray.
                [~, ~, ~, adjoint_ray_positions_auxilary_left, ~, ~, ~] = feval(solve_auxiliary_ray,...
                    angle_left, ray_positions_adjoint,...
                    false, true, rotation_matrix_left);
                
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
        
        
        if ~strcmp(para.auxiliary_method, 'paraxial')
            
            % right auxiliary ray
            ray_positions_auxilary_right_allreceivers(receiver_index_position:ind_receiver * dim,...
                1:size(ray_positions_auxilary_right, 2) ) = ray_positions_auxilary_right;
            
            
            % right auxiliary ray
            adjoint_ray_positions_auxilary_right_allreceivers(receiver_index_position:ind_receiver * dim,...
                1:size(adjoint_ray_positions_auxilary_right, 2) ) = adjoint_ray_positions_auxilary_right;
            
        else
            if  para.do_perturb_initial_position
                
                
                % get the starting index for filling the matrix for the polar direction of the
                % auxiliary rays
              %  receiver_index_direction = ind_receiver * (dim-1) - (dim-2);
                receiver_index_direction = ind_receiver * dim - (dim-1);
                % if the method for computing the auxiliary rays is paraxial, and
                % if $d direction/d x(s_0)$ is requested, fill the rows associated
                % with the current receiver
                
                % get the ray directions for the forward auxiliary rays
                % left auxiliary ray
                ray_directions_auxilary_allreceivers(receiver_index_direction:ind_receiver * (dim),...
                    1:size(ray_directions_auxilary, 2) ) = ray_directions_auxilary;
                
                
                % get the ray directions for the adjoint auxiliary rays
                % left auxiliary ray
                adjoint_ray_directions_auxilary_allreceivers(receiver_index_direction:ind_receiver * (dim),...
                    1:size(adjoint_ray_directions_auxilary, 2) ) = adjoint_ray_directions_auxilary;
                
                
            end
            
            
        end
        
        
    end
    
    
    % fill the matrix for catesian position of the end point of the linked ray for each receiver
    cartesian_position_endpoint_allreceivers(:, ind_receiver) = cartesian_position_endpoint;
    
    % fill the matrix for the number of solved rays for ray
    % linking for each receiver
    num_rays_allreceivers(ind_receiver) = num_rays;
    
    rayspacing_allreceivers(ind_receiver) = rayspacing_receiver;
   
    % add the indices of columns for the current receiver to nonzero column indices for all receivers
   % ind_gridpoints = [ind_gridpoints; columns];
    % add the indices of rows (the index of the current receiver) to nonzero row indices for all receivers
   % ind_receivers = [ind_receivers; ind_receiver * ones(size(columns))];
    % add the values (interpolation coefficient) of the
    % nonzeros elements of the sparse matrix with row and
    % columns indicated by 'ind_receivers' and 'ind_gridpoints',
    % respectively.
   % val_coeffs = [val_coeffs; vals];
    
end

% using 'Regula-Falsi' method, the optimal polar
% initial direction will not be used as the initial
% guess for the next UST iteration
if strcmp(para.raylinking_method, 'Regula-Falsi')
    optimal_polar_initial_direction_allreceivers = [];
end




end