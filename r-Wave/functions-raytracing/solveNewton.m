function [interp_coeff_vec, polar_initial_direction, cartesian_pos_endpoint,...
    num_rays] = solveNewton(solve_ray, polar_direction_receiver, polar_direction_initial_guess,...
    varepsilon, max_iter, res_initial_guess)
%SOLVENEWTON solves a ray linking inverse problem using the Newton's method
%
% DESCRIPTION:
% solveNewton solves the ray linking inverse problem using the Newton's
% method. This approach is used only for ray linking in 3D case.
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
%
% OUTPUTS:
%       interp_coeff_vec         - interpolation coefficients for the grid
%                                  points inside the mask
%       polar_initial_direction -  dim-1 x 1 optimal polar initial direction
%                                  unit vector from emitter to the receiver
%       cartesian_pos_endpoint   - dim x 1 cartesian position of the
%                                  end point of the ray
%       num_rays                 - the number of rays traced for solving the ray
%                                  linking inverse problem

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2020 Ashkan Javaherian





% trace the first ray, calculate the end point and its residual, and
% associated interpolation coefficients

if nargin < 6
    % calculate and store the residual for the initial guess
    [polar_direction_residual, interp_coeff_vec, cartesian_pos_endpoint, polar_direction_endpoint] = ...
        feval (solve_ray, polar_direction_initial_guess, polar_direction_receiver, true);
    res_initial_guess = polar_direction_residual;
    
else
    [~, interp_coeff_vec, cartesian_pos_endpoint, polar_direction_endpoint] = ...
        feval(solve_ray, polar_direction_initial_guess, [], true);
    polar_direction_residual = res_initial_guess;
end
% store the polar direction of th ened point for the initial guess
polar_direction_endpoint_initial_guess = polar_direction_endpoint;

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
   
    
%  a scalar controlling perturbations applied to the polar initial direction
%  for calculating the Jacobian matrix (fixed)
tau = 1e-15;

% a scalar upper bound for the maximum to minimim singular value 
% of the Hessian matrix (fixed)
sigma_max = 1e4;

% a scalar lower bound for the minimim singular value of the Hessian matrix
% (fixed)
sigma_min = 1e-6;
  
% set the bollean indicating an ill-conditioned Hessian matrix false
bad_hessian = false;


% update the residual norm (start by the initial guess)
norm_res = norm_res_initial_guess;
    
    while norm_res > varepsilon  &&   num_rays < max_iter
        
        % the Hessian matrix is ill-conditioned
        if bad_hessian 
            
        % increase tau by 10
        tau = 10 * tau;
        % restart the iterations using the initial guess for the initial polar direction
        % and its associated polar direction for the end point of the ray
        polar_initial_direction = polar_direction_initial_guess;
        polar_direction_endpoint = polar_direction_endpoint_initial_guess;
        end
        
        % update the polar initial direction, polar direction of the end
        % point of the ray, and residual
        [polar_initial_direction, polar_direction_endpoint,...
        polar_direction_residual, norm_res, bad_hessian, num_rays] = updateNewton(solve_ray,...
        polar_direction_receiver, polar_initial_direction, polar_direction_endpoint,...
        polar_direction_residual, norm_res, num_rays, tau, sigma_max, sigma_min);
    
       
       
        
    end
    
    % if the last residul is greater than the residual for
    % the initial guess, replace the last initial
    % direction and its residual by those for the
    % initial guess
    if norm_res > norm_res_initial_guess
        polar_initial_direction = polar_direction_initial_guess;
    end
    
    % calculate the interpolation coefficients for the linked ray
    [ ~ , interp_coeff_vec, cartesian_pos_endpoint] = feval(solve_ray, polar_initial_direction,...
        polar_direction_receiver, true);
    
end



end