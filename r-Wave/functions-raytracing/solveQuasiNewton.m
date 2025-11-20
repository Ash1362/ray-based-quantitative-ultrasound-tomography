function [interp_coeff_vec, polar_initial_direction, cartesian_pos_endpoint,...
    num_rays] = solveQuasiNewton(solve_ray, polar_direction_receiver, polar_direction_initial_guess,...
    varepsilon, max_iter, res_initial_guess, varargin)
%SOLVEQUSAINEWTON solves a ray linking inverse problem using Quasi-Newton
%methods
%
% DESCRIPTION:
% solveQuasiNewton solves the ray linking inverse problem using Quasi-Newton
% methods. This approach, which is used only for ray linking in 3D case,
% solves the inverse problem of ray linking in the form of a nonlinear system
% of equations using derivative-free method.
% USAGE:
%
%
% INPUTS:
%      solve_ray                - the function handle for solving the ray and the
%                                 parameters associated with the end point of the ray
%      polar_direction_receiver - the dim-1 x 1 polar direction of the straight line
%                                 connecting emitter to receiver
%      polar_direction_initial_guess - the dim-1 x 1 initial guess for the polar initial direction of
%                                       the ray. 
%      varepsilon                - stopping tolerance for ray linking                              
%      max_iter                  - the maximum number of iterations for ray linking
%      res_initial_guess         - the dim-1 x 1 residual obtained for the end point of the ray initialised
%                                 by the initial guess for the polar direction

% OPTIONAL INPUTS:
%      Method                   - the Quasi-Newton approach used for ray
%                                 linking (default : 'Good-Broyden')
%      initial_derivative       - first derivative of the end point of the
%                                 ray with respect to the initial polar
%                                 direction of the ray. This can be an
%                                 identity matrix ('identity'), or computed 
%                                 using finite diferences using a large
%                                 perturbation for ensuring the stability
%                                 ('finite-difference')
%                                 (default : 'finite-difference')

%                                 
% OUTPUTS:
%       interp_coeff_vec         - interpolation coefficients for the grid
%                                  points inside the mask
%       polar_initial_direction -  dim-1 x 1 optimal polar initial direction
%                                  unit vector from emitter to the receiver
%       cartesian_pos_endpoint   - dim x 1 vector of the Cartesian position of the
%                                  end point of the ray
%       num_rays                 - the number of rays for solving the ray
%                                  linking inverse problems

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 



para = [];
para.Method = 'Good-Broyden';
para.initial_derivative = 'finite-difference';
para.smooth = true;

if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% trace the first ray, calculate the end point and its residual, and
% associated interpolation coefficients

if isempty(res_initial_guess)
    % calculate and store the residual for the initial guess
    [polar_direction_residual, interp_coeff_vec, cartesian_pos_endpoint, polar_direction_endpoint] = ...
        feval (solve_ray, polar_direction_initial_guess, polar_direction_receiver, true);
    res_initial_guess = polar_direction_residual;
    
else
    [~, interp_coeff_vec, cartesian_pos_endpoint, polar_direction_endpoint] = ...
        feval(solve_ray, polar_direction_initial_guess, [], true);
    polar_direction_residual = res_initial_guess;
end


% calculate the norm of residual for the initial guess
norm_res_initial_guess = norm(res_initial_guess);


% initialise the counter for calculating the number
% of rays for solving the ray linking inverse
% problem
num_rays = 1;
% update the polar initial direction of the ray
polar_initial_direction = polar_direction_initial_guess;
% update the residual

if norm_res_initial_guess > varepsilon
  
    % choose the bounds ( must be 0.2)
    bounds = [polar_initial_direction - 0.5,...
        polar_initial_direction + 0.5];
    
    % a scalar controlling the projection of an update of the initial angle
    % outside the bounds to the interior of the feasible. this parameter
    % controls the distance of the projected initial angle from the associated bound (fixed)
    zeta = 0.5;
    
    % the maximum permissible reduction for the amplitude of search direction
    % by a projection to the feasible set (safe-guard). An update very close to a bound may
    % need to be significanly reduced in amplitude in order to remain a feasible set.
    % However, very small-size step size may cause instability for the next iteration.
    % so we enforce a bound on the maximum reduction of the search direction
    kappa = 1e-6;
    
    
    if para.smooth
        
        %  a scalar controlling the update of the derivative matrix (initial value)
        tau = 1;
        % a parameter for controlling the change of tau in the adapative smoothing
        % loop
        switch para.Method
            case 'Good-Broyden'
                % an increment factor for moving tau from 1
                eta = 0.01;
            case 'BFGS'
                % a reduction weight for tau inside the loop
                eta = 0.5;
        end
        % a scalar upper bound for the maximum to minimim singular value
        % of the derivative matrix (fixed)
        sigma_max = 1e4;
        % a scalar lower bound for the minimim singular value of the derivative matrix
        % (fixed)
        sigma_min = 1e-4;
    else
        
        % allocate empty variables for smoothing parameters,
        % if smoothing of the derivative matrix is not requested
        tau = [];
        eta = [];
        sigma_max = [];
        sigma_min = [];
        
    end
    
    
    switch para.initial_derivative
        case 'finite-difference'
            % calculate a derivative matrix roughly using a very large perturbation in
            % order to obtain a well-posed initial guess matrix
            derivative_matrix = calcRayJacobian(solve_ray, polar_initial_direction,...
                polar_direction_endpoint, 1e-6);
        case 'identity'
            derivative_matrix = eye(size(polar_direction_endpoint, 1));
    end
    
    % increase the number of rays with 2 because of calculating the Jacobian matrix
    % in the above line
    num_rays = num_rays + 2;
    
    % update the residual norm (start by the initial guess)
    norm_res = norm_res_initial_guess;
    
   
    varepsilon_search_direction = 1e-2 * varepsilon;
    
    while norm_res > varepsilon    &&  num_rays < max_iter
        
        % update the polar initial direction, polar direction of the end
        % point of the ray, and residual
        [polar_initial_direction, polar_direction_residual, derivative_matrix,...
            norm_res, norm_search_direction] = updateQuasiNewton(solve_ray,...
            polar_direction_receiver, polar_initial_direction, ...
            polar_direction_residual, derivative_matrix, bounds, zeta,...
            kappa, tau, eta, sigma_max, sigma_min, para);
        
       % increase the number of rays with 1 because of tracing a ray
       % for solving the forward operator of ray linking
       num_rays = num_rays + 1; 
       
       if norm_search_direction < varepsilon_search_direction
           num_rays = max_iter;
           break;
       end
        
        
    end
    
    % if the last residul is greater than the residual for
    % the initial guess, replace the last initial
    % direction and its residual by those for the
    % initial guess
    if norm_res > norm_res_initial_guess
        
        polar_initial_direction = polar_direction_initial_guess;
   
        % allocate empty variable for interpolation coefficients
        interp_coeff_vec = [];
    else   
        
    % If ray linking is successful, compute the interpolation coefficients for the linked ray
    [ ~ , interp_coeff_vec, cartesian_pos_endpoint] = feval(solve_ray, polar_initial_direction,...
        polar_direction_receiver, true);
    
    end
    
end


end