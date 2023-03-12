function [polar_initial_direction, cartesian_pos_endpoint,...
    num_rays, ray_positions, ray_acoustic_length, ray_absorption, rayspacing_receiver] =...
    solveSecantFullWave(solve_ray, polar_direction_receiver, polar_direction_initial_guess,...
    varepsilon, max_iter, res_initial_guess)
%SOLVESECANT solves a ray linking inverse problem using the Secant method
%
% DESCRIPTION:
% solveSecant solves the ray linking inverse problem using the 'Secant'
% method. This approach is used only for ray linking in 2D case.
% USAGE:
%
%
% INPUTS:
%      solve_ray                - the function handle for solving the ray and the
%                                 parameters associated with the end point of the ray
%      polar_direction_receiver - the polar direction of the straight line
%                                 connecting emitter to receiver
%      polar_direction_initial_guess - the initial guess for the polar initial direction of
%                                       the ray. 
%      varepsilon                - stopping tolerance for ray linking                              
%      max_iter                  - the maximum number of iterations for ray linking
%      res_initial_guess         - the residual obtained for the end point of the ray initialised
%                                 by the initial guess for the polar direction

% OPTIONAL INPUTS:
%
% OUTPUTS:
%       polar_initial_direction -  optimal polar initial direction
%                                  unit vector from emitter to the receiver
%       cartesian_pos_endpoint   - the cartesian position of the end point of the ray
%       num_rays                 - the number of rays traced for solving the ray
%                                  linking inverse problems
%       ray_positions            - the Cartesian position [m] of the points on
%                                  the linked ray
%       ray_acoustic_length      - the accumulated acoustic length on the
%                                  points for the linked ray
%       ray_absorption           - the accumulated acoustic absorption on
%                                  the points for the linked ray
%       rayspacing_receiver      - the ray spacing [m] for the two last points 
%                                  on the linked ray


% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian 


% trace the first ray, calculate the end point and its residual, and
% associated interpolation coefficients

% define the paprameters for applying constraints
kappa = 1e-6;
zeta = 0.5;

if nargin < 6
    % calculate and store the residual for the initial guess
    [polar_direction_residual, cartesian_pos_endpoint, ~, ...
        ray_positions, ray_acoustic_length, ray_absorption, rayspacing_receiver] = ...
        feval (solve_ray, polar_direction_initial_guess, polar_direction_receiver, true, false);
    res_initial_guess = polar_direction_residual;
    
else
    [~, cartesian_pos_endpoint, ~, ...
        ray_positions, ray_acoustic_length, ray_absorption, rayspacing_receiver] = ...
        feval(solve_ray, polar_direction_initial_guess, [], true, false);
    polar_direction_residual = res_initial_guess;
end
% initialise the counter for calculating the number
% of rays for solving the ray linking inverse
% problem
num_rays = 1;
% update the polar initial direction of the ray
polar_initial_direction = polar_direction_initial_guess;


if abs(polar_direction_residual) > varepsilon
    
    % choose the bounds
    bounds = [polar_initial_direction - 0.1,...
        polar_initial_direction + 0.1];
    
    
    % initialise the search direction with a value in range of the stopping
    % tolerance
    search_direction = varepsilon;
    
    
    while abs(polar_direction_residual) > varepsilon  &&   num_rays < max_iter
       
        
        % check to see if an update of the initial direction using the 
        % current search direction does not exceed the lower bound. if so,
        % reduce the size of the search direction so that the update
        % will be inside the feasible set
        if polar_initial_direction + search_direction < bounds(1)
            % define a vector psi
            psi = zeta * (bounds(1) - polar_initial_direction)...
                ./ search_direction;
            % modify the search direction for a projection to the interior of the feasible set
            search_direction = sign(psi) .* max(abs(psi), kappa).*...
                search_direction;
        end
        % check to see if an update of the initial direction using the 
        % current search direction does not exceed the upper bound. if so,
        % reduce the size of the search direction so that the update
        % will be inside the feasible set
        if polar_initial_direction + search_direction > bounds(2)
            % define a vector psi
            psi = zeta * (bounds(2) - polar_initial_direction)...
                ./ search_direction;
            % modify the search direction for a projection to the interior of the feasible set
            search_direction = sign(psi) .* max(abs(psi), kappa).*...
                search_direction;
        end
        
        
        
        % update the initial direction of the ray using the updated
        % constarined search direction
        polar_initial_direction = polar_initial_direction + ...
            search_direction;
        
        % store the updated value of the residual (will be used for
        % updating the search direction
        previous_polar_direction_residual = ...
            polar_direction_residual;
        
        
        % trace the ray and calculate the residual from
        % the end point of the ray
        polar_direction_residual = feval(solve_ray, polar_initial_direction,...
            polar_direction_receiver, false, false);
        % increase the counter
        num_rays = num_rays + 1;
       
        % update the search direction using 'Secant'
        % method
        search_direction = search_direction / (previous_polar_direction_residual -...
            polar_direction_residual) * polar_direction_residual;
        
       
    end
    
    % if the last residul is greater than the residual for
    % the initial guess, replace the last initial
    % direction and its residual by those for the
    % initial guess
    if abs(polar_direction_residual) > abs(res_initial_guess)
        polar_initial_direction = polar_direction_initial_guess;
    end
    
    % calculate the interpolation coefficients for the linked ray
    [ ~, cartesian_pos_endpoint, ~, ray_positions, ray_acoustic_length,...
        ray_absorption, rayspacing_receiver] = feval(solve_ray, polar_initial_direction,...
        polar_direction_receiver, true, false);
    
end



end