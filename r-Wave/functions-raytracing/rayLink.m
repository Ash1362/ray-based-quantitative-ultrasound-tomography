function [system_matrix_single_emitter, optimal_polar_initial_direction_allreceivers,...
    cartesian_position_endpoint_allreceivers, num_rays_allreceivers] = rayLink(ray_fields_params,...
    refractive, cartesian_position_emitter, cartesian_position_receivers, polar_direction_allreceivers,...
    polar_initial_direction_allreceivers, xvec, yvec, zvec, pos_grid_first,...
    pos_grid_end, grid_spacing, ray_spacing, grid_size, dim, detec_radius, mask, para)
%RSAYLINK constructs the system matrix for a single emitter
%
% DESCRIPTION:
%           rayLink solves the ray linking inverse problem for all receivers
%           corresponding to an emitter, and store the information about the end
%           point of the linked rays. Subsequently, a num_receiver x num_gridpoints
%           sparse metrix, which contains the interpolation coefficients for the linked rays
%           for all receivers, are cponstructed
% USAGE:
%
%
% INPUTS:
%       ray_fields_params      - a struct containing the specified variables for ray tracing
%                              using a 'Bilinear' intrerpolation, this includes the directional
%                              gradients of the refrective index
%                              distribution:
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
%       cartesian_position_emitter - a dim x 1 cartesian position of the emitter
%       cartesian_position_receivers - a dim x num_receiver cartesian
%                                       position of all receivers
%       polar_direction_allreceivers - a (dim-1) x num_receiver vector of the polar direction of 
%                                    unit geometrical vectors from emitter
%                                    to all the receivers
%       polar_initial_direction_allreceivers - a (dim-1) x num_receiver matrix of the polar initial
%                                    direction of the rays for all receivers
%       xvec                 - the x vector of grid points
%       yvec                 - the y vector of grid points
%       zvec                 - the z vector of grid points
%       pos_grid_first       - a dim x 1 Cartesian position of the first index of
%                               the grid
%       pos_grid_end         - a dim x 1 Cartesian position of the end index of
%                              the grid
%       dx                   - a scalar representing the grid spacing, the same
%                              along all the Cartesian coordinates [m]
%       ray_spacing          - a saclar representing the ray spacing [m]
%       grid_size            - the size of the grid
%       dim                  - the dimension of the medium
%       mask                 - a binary mask with the same size as the grid 
%                              for an indication to the grid points
%                              associated with inhomogeneities (object)
%       ray_spacing          - the spacing along the ray [m]
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
%       cartesian_position_endpoint_allreceivers - dim x num_receiver matrix of cartesian position 
%                                              of the end point of the ray
%                                              for all receivers
%       num_rays_allreceivers               - num_receiver x 1 vector of the number of ray tracing
%                                              for ray linking between the
%                                              emitter and each of all receivers
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2020 Ashkan Javaherian 


% number of receivers
num_receiver = size(polar_direction_allreceivers, 2);


% a binary vector, which is true when the distances of the emeitter to
% receivers are sufficiently large
switch para.binaries_emitter_receiver
    case 'distances'
        binaries_emitter_receiver = vecnorm(cartesian_position_receivers - cartesian_position_emitter) > para.minimum_distance;
    case 'open_angle'
        binaries_emitter_receiver = (calcAngleEmitterReceiver(cartesian_position_emitter,...
            cartesian_position_receivers, []) < para.open_angle)' & ...
            vecnorm(cartesian_position_receivers - cartesian_position_emitter) > 0;
end



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
        
    case 'Bspline' 
        
       % allocate empty variables for the directional gradient fields 
        refractive_gradient_x = [];
        refractive_gradient_y = [];
        refractive_gradient_z = [];
        
        % get parameters for spline interpolation
        raytogrid_indices_x = ray_fields_params.raytogrid_indices_x;
        raytogrid_indices_y = ray_fields_params.raytogrid_indices_y;
        raytogrid_indices_z = ray_fields_params.raytogrid_indices_z;
        raytogrid_coeff_matrix = ray_fields_params.raytogrid_coeff_matrix;
        raytogrid_coeff_derivative_matrix = ray_fields_params.raytogrid_coeff_derivative_matrix;
end


% the indices of rows of nonzero elements of a sparse system matrix, and
% represents the indices of receivers
ind_receivers = [];

% the indices of columns of nonzeros elements of a sparse system matrix,
% and represents the indices of grid points in a Matlab column-wise order
ind_gridpoints = [];

% the values of nonzeros elements of the sparse system matrix, and includes the
% coefficient for interpolation between the rays
% linking the emitter to receivers (row indices)and grid points (column indices)
val_coeffs = [];

% the polar initial direction of the linked rays
optimal_polar_initial_direction_allreceivers = zeros(dim-1, num_receiver);

% allocate matrices for storing the information from ray linking
% the cartesian position of the end point of the linked rays
cartesian_position_endpoint_allreceivers = zeros(dim, num_receiver);

% a vector containg the number of number of traced rays for ray linking
% between the vector and each of the receivers
num_rays_allreceivers = zeros(num_receiver, 1);

% define a handle function for tracing the ray, and obtain the information
% this function handle solves the forward problem of ray linking inverse
% problem
solve_ray = @(polar_initial_direction, polar_direction_receiver, calc_coeffs)calcRayParameters(...
    refractive, refractive_gradient_x, refractive_gradient_y, refractive_gradient_z,...
    cartesian_position_emitter, polar_direction_receiver, polar_initial_direction,...
    xvec, yvec, zvec, pos_grid_first, pos_grid_end, grid_spacing, ray_spacing, ...
    grid_size, dim, detec_radius, mask, raytogrid_indices_x, raytogrid_indices_y,...
    raytogrid_indices_z, raytogrid_coeff_matrix, raytogrid_coeff_derivative_matrix, ...
    calc_coeffs, para.raylinking_method, para.interp_method);

%% ========================================================================
% INITIAL GUESS
%==========================================================================

switch para.raylinking_method
    
    
    case 'Regula-Falsi'
        
        % choose initial guess for the left and right polar
        % initial direction.
        switch para.raylinking_initialisation
            case 'Local'
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
                [ ~, ~, ~, polar_direction_endpoint_left] = feval(solve_ray, polar_initial_direction_left, [], false);
                [ ~, ~, ~, polar_direction_endpoint_right] = feval(solve_ray, polar_initial_direction_right, [], false);
                
            case 'Global'
                % By a ray initialisation using a 'Global' approach,
                % the left and right polar intial directions (initial interval)
                % are chosen individually for each receiver. A number of rays
                % with even angles are simultaneously solved. For each receiver,
                % a pair of rays with endpoints closest to the receiver,
                % and encompass the receiver are chosen as the
                % initial guess for the the left and right polar
                % initial directions. This pair will give opposite
                % signs for the residual corresponding to the polar direction of
                % the end point of the rays, and forms a very tight (optimal) initial interval
                % for applying the 'Regula Falsi' method.
                [polar_direction_initial_guess_allreceivers, res_initial_guess_allreceivers] =...
                    initialiseAngles(solve_ray, polar_direction_allreceivers,...
                    dim, 1, para.varepsilon, para.raylinking_method);
        end
        
    otherwise
        
        
        % choose initial guess for the left and right polar
        % initial direction.
        switch para.raylinking_initialisation
            
            
            case 'Local'
                
                % By a ray initialisation using a 'Local' approach,
                % the optimal ray from solving the ray linking inverse problem
                % for the last previous UST iteration will be used
                % for the initial guess for solving the ray linking
                % inverse problem for the current UST iteration.
                polar_direction_initial_guess_allreceivers = [];
                res_initial_guess_allreceivers = [];
            case 'Global'
                
                
                % By a ray initialisation using a 'Global'
                % approach, a number of rays with even initial
                % directions are first solved. For each receiver, the ray
                % with the closest distance to the receiver may be
                % used as the initial guess. But for deciding on
                % that, it is first compared with the ray with
                % an initial guess from the optimal solution of the
                % ray linking for the last previous UST iteration. The
                % closer ray will be used as initial guess for
                % the polar initial direction of the ray for solving the
                % ray linking inverse problem in the current UST iteration.
                [polar_direction_initial_guess_allreceivers, res_initial_guess_allreceivers] =...
                    initialiseAngles(solve_ray, polar_direction_allreceivers,...
                    dim, 1, para.varepsilon, para.raylinking_method);
        end
        
        
end


% solve the ray linking inverse problem for each individual emitter-receiver pair
for ind_receiver = 1 : num_receiver
    
    % disp(ind_receiver)
    % avoid ray linking for a sufficiently close emitter-receiver pair,
    % becuse the rays only travel through the water for
    % this pair
    if binaries_emitter_receiver(ind_receiver)
        
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
                
                [interp_coeff_vec, ~ , cartesian_position_endpoint, num_rays] = solveRegulaFalsi(solve_ray,...
                    polar_direction_allreceivers(:, ind_receiver), polar_direction_initial_guess,...
                    res_initial_guess, para.varepsilon, para.max_iter);
                
                
                
            case 'Secant'
                [interp_coeff_vec, polar_initial_direction, cartesian_position_endpoint, num_rays]=...
                    solveSecant(solve_ray, polar_direction_allreceivers(ind_receiver),...
                    polar_initial_direction_allreceivers(ind_receiver), para.varepsilon, para.max_iter);
                
            case 'Newton'
                % ray linking using Newton's method
                [interp_coeff_vec, polar_initial_direction, cartesian_position_endpoint, num_rays]=...
                    solveNewton(solve_ray, polar_direction_allreceivers(:, ind_receiver),...
                    polar_initial_direction_allreceivers(:, ind_receiver), para.varepsilon, para.max_iter);
                
            case 'Quasi-Newton'
                
                
                % ray linking using Quasi-Newton method
                link_args = {'Method', 'Good-Broyden', 'initial_derivative', 'finite-difference',...
                    'smooth', true};
                
                [interp_coeff_vec, polar_initial_direction, cartesian_position_endpoint, num_rays]=...
                    solveQuasiNewton(solve_ray, polar_direction_allreceivers(:, ind_receiver),...
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
        interp_coeff_vec = sparse(prod(grid_size), 1);
        cartesian_position_endpoint = nan;
        num_rays = 0;
        
        
        
    end
    
    
    
    % fill the matrix for updated polar initial direction for each receiver
    % using 'Regula-Falsi' method, the optimal polar
    % initial direction will not be used as the initial
    % guess for the next UST iteration
    if ~strcmp(para.raylinking_method, 'Regula-Falsi')
        optimal_polar_initial_direction_allreceivers(:, ind_receiver) = polar_initial_direction;
    end
    
    % fill the matrix for catesian position of the end point of the linked ray for each receiver
    cartesian_position_endpoint_allreceivers(:, ind_receiver) = cartesian_position_endpoint;
    
    % fill the matrix for the number of solved rays for ray
    % linking for each receiver
    num_rays_allreceivers(ind_receiver) = num_rays;
    
    
    % find the indices of columns and values of the nonzero elements of the
    % system matrix from the vector of interpolation coefficient for an individual receiver
    % for all the grid points
    [columns, ~, vals]  = find(interp_coeff_vec);
    % add the indices of columns for the current receiver to nonzero column indices for all receivers
    ind_gridpoints = [ind_gridpoints; columns];
    % add the indices of rows (the index of the current receiver) to nonzero row indices for all receivers
    ind_receivers = [ind_receivers; ind_receiver * ones(size(columns))];
    % add the values (interpolation coefficient) of the
    % nonzeros elements of the sparse matrix with row and
    % columns indicated by 'ind_receivers' and 'ind_gridpoints',
    % respectively.
    val_coeffs = [val_coeffs; vals];
    
    
end



% using 'Regula-Falsi' method, the optimal polar
% initial direction will not be used as the initial
% guess for the next UST iteration
if strcmp(para.raylinking_method, 'Regula-Falsi')
    optimal_polar_initial_direction_allreceivers = [];
end


% construct the sparse system matrix for all
system_matrix_single_emitter = sparse(ind_receivers, ind_gridpoints,...
    val_coeffs, num_receiver, prod(grid_size));

            

end
