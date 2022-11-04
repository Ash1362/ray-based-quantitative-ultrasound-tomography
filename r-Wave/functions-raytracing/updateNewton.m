function [updated_polar_initial_direction, updated_polar_direction_endpoint,...
    updated_polar_direction_residual, updated_norm_res, bad_hessian, num_rays] = updateNewton(solve_ray,...
    polar_direction_receiver, polar_initial_direction, polar_direction_endpoint,...
    polar_direction_residual, norm_res, num_rays, tau, sigma_max, sigma_min)
%UPDATENEWTON solves an iteration of ray linking inverse problem using the Newton's method
%
% DESCRIPTION:
% updateNewtont solves an iteration of ray linking inverse problem using
% the Newton's method. This approach is used only for ray linking in 3D case.
% USAGE:
%
%
% INPUTS:
%      solve_ray                - the function handle for solving the ray and the
%                                 parameters associated with the end point of the ray
%      polar_direction_receiver - the polar direction of a geometrical
%                                 vector from emitter to receiver 
%      polar_initial_direction  - a 2 x 1 vector of the polar initial
%                                 direction of the last update of the ray
%      polar_direction_endpoint - a 2 x 1 vector of the polar direction of a geometrical
%                                 unit vector from emitter to the last point of
%                                 the updated ray  
%      polar_direction_residual - a 2 x 1 vector of the polar residual, the discrepancy
%                                 between the polar unit vector from emitter to 
%                                 the last point of the ray and a polar unit vector from emitter to the receiver
%      norm_res                 - the norm of residual
%      num_rays                 - the total number of traced rays from the
%                                 first iteration of ray linking
%      tau                      - a scalar controlling perturbations applied to the 
%                                 two components of the current polar initial direction
%      sigma_max               - a scalar bound for the maximum to minimim singular value 
%                                of the Hessian matrix
%      sigma_min               - a scalar bound for the minimim singular value of the Hessian matrix
% OPTIONAL INPUTS:
%
% OUTPUTS:
%      updated_polar_initial_direction  - a 2 x 1 vector of the update of polar initial
%                                 direction of the last update of the ray
%      updated_polar_direction_endpoint - a 2 x 1 vector of the update of polar 
%                                 direction of a geometrical a unit vector from 
%                                 emitter to the last point of the updated ray  
%      updated_polar_direction_residual - a 2 x 1 vector of the update of polar residual, 
%                                 the discrepancy between the polar unit vector from emitter to 
%                                 the last point of the ray and a polar unit vector from emitter
%                                 to the receiver
%      updated_norm_res                 - the updated norm of residual 
%      bad_hessian                      - a binary variable indiating whether the singular values
%                                        of the updated Hessian matrix is
%                                        ill-conditioned
%      num_rays                 - the total number of traced rays from the first iteration 
%                                 of ray linkinhg to after the current iteration


% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2020 Ashkan Javaherian 



% calculate the Jacobian matrix
[jacobian_matrix] = calcRayJacobian(solve_ray, polar_initial_direction,...
    polar_direction_endpoint, tau);

% increase the number of rays with 2 because of calculating the Jacobian matrix
% in the above line
num_rays = num_rays + 2;

% calculate the Hessian matrix
hessian_matrix = jacobian_matrix' * jacobian_matrix;

% check if all elements of the Hessian matrix are finite
if all(isfinite(hessian_matrix(:)))
    
    % compute the gradient of the objective function, the right-hand-side of 
    % the normal equation
    grad = - jacobian_matrix' * polar_direction_residual;
    % update the bound for the minimal singular value by including the gradient of the
    % objective function
    lambda = min(norm(grad), sigma_min);
    
    % calculate the singular values of the Hessian matrix
    singular_values = eig(hessian_matrix);
    
    
    hessian_smooth_iteration = 0;
    % a loop for ensuring the singular values of the Hessian matrix are
    % within our chosen range. this loop is terminated when the
    % following conditions are not satisfied.The maximum number of
    % iterations of the lopp prevents the perturantion tau becomes so large
    % (A small perturbation may cause instability in the presence of singularities,
    % and also large perturnation does not give a precise Jacobian matrix)
    while hessian_smooth_iteration < 5  && ( any(singular_values < lambda) ||  max(singular_values)/min(singular_values)>sigma_max )
        
        % increase the parameter tau by 10
        tau = 10 * tau;
        % update the Jacobian matrix using the updated tau
        [jacobian_matrix] = calcRayJacobian(solve_ray, polar_initial_direction,...
            polar_direction_endpoint, tau);
        
        % increase the number of rays with 2 because of calculating the Jacobian matrix
        % in the above line
        num_rays = num_rays + 2;

        
        % update the Hessian matrix
        hessian_matrix = jacobian_matrix' * jacobian_matrix;
        
        % terminate the loop if any components of the Hessian matrix
        % becomes infinite (a safe-guard to avoid instability)
        
        if   all(isfinite(hessian_matrix(:)))
            
            % update the singulr values of the Hessian matrix
            singular_values = eig(hessian_matrix);
            hessian_smooth_iteration = hessian_smooth_iteration + 1;
        else
            break;
        end
        
    end
    
    
end



if    all(isfinite(hessian_matrix(:)))
    
    
    % check if the smoothing loop was terminated before the maximum number
    % of iteration, and give a desired Hessian matrix with singular values
    % satisfying the conditions for well-posedness
    if ~any(singular_values < lambda)  &&  max(singular_values)/min(singular_values) < sigma_max
        
        % a boolean determining whether the singular values of the Hessian matrix
        % do not satisfy our conditions, and will thus not give a descent
        % direction
        bad_hessian = false;
        % upodate the right-hand-side of the normal equation using the updated
        % Jacobian matrix
        grad = - jacobian_matrix' * polar_direction_residual;
        % calculate the search direction 
        search_direction = hessian_matrix \ grad;
        
        
        % applying the line search technqu for satisfying the Armijo condition
        % initialse the step length (fixed)
        step_length = 1;
        % choose a factor for reduction of the step length (fixed)
        reduction_factor = 0.5;
        % choose a parameter for controlling the Armijo condition (fixed)
        mu = 1e-4;
        % apply the iterative line-search technique.
        % the maximum number of iterations for an iterative reduction of the step length by
        % a factor 0.5 is set 5, and is always fixed.
        
        for ind_linesearch = 1 : 5
            
            % update the polar initial direction of the ray using the
            % updated step length
            updated_polar_initial_direction = polar_initial_direction...
                + step_length * search_direction;
            
            % calculate the ray using the updated polar initial direction,
            % and update the end point of the ray
            [updated_polar_direction_residual, ~, ~, updated_polar_direction_endpoint]...
                = feval(solve_ray, updated_polar_initial_direction, polar_direction_receiver, false);
            
            % increase the number of rays with 1 because of tracing a ray
            % for solving the forward operator of ray linking
            num_rays = num_rays + 1;

            % update the norm of the residual
            updated_norm_res = norm(updated_polar_direction_residual);
            
            % terminate the line-search loop, if the Armijo condition is
            % satisfied, otherwise continue the loop with a step length
            % reduced by 0.5
            if updated_norm_res <= norm_res + mu * step_length * grad' * search_direction
                break;
            end
            % reduce the step length by 0.5
            step_length = reduction_factor * step_length;
            
            
        end
        
    else
        
        % a boolean determining whether the singular values of the Hessian matrix
        % do not satisfy our conditions, and will thus not give a descent
        % direction
        bad_hessian = true;
        % do not accept the updated initial direction, and replace it by the
        % last update, if the Hessian matrix is not sufficiently
        % well-conditioned
        updated_polar_initial_direction = polar_initial_direction;
        updated_polar_direction_endpoint = polar_direction_endpoint;
        updated_polar_direction_residual = polar_direction_residual;
        updated_norm_res = norm_res;
        
    end
    
    
else
    
    
    % avoid any updates if at least one of of the elements of the Hessian matrix is infinite
    bad_hessian = true;
    % do not accept the updated initial direction, and replace it by the
    % last update, if the Hessian matrix is not sufficiently
    % well-conditioned
    updated_polar_initial_direction = polar_initial_direction;
    updated_polar_direction_endpoint = polar_direction_endpoint;
    updated_polar_direction_residual = polar_direction_residual;
    updated_norm_res = norm_res;
    
end


end