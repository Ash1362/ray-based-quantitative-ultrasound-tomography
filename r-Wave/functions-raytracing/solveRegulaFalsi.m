function [interp_coeff_vec, polar_initial_direction, cartesian_pos_endpoint,...
    num_rays] = solveRegulaFalsi(solve_ray, polar_direction_receiver, ...
    polar_direction_initial_guess, res_initial_guess, varepsilon, max_iter)
%SOLVEREGULAFALSI solves a ray linking inverse problem using the Regula Falsi method
%
% DESCRIPTION:
% solveRegulaFalsi method solves the ray linking inverse problem using the 'Regula Falsi'
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
%                                       the ray. This pair provides the initial guess for the
%                                       intervals used for applying Regula Falsi
%                                       approach
%      res_initial_guess         - the residual obtained for the end point of the ray initialised
%                                 by the initial guess for the polar direction
%      varepsilon                - stopping tolerance for ray linking                              
%      max_iter                  - the maximum number of iterations for ray linking

% OPTIONAL INPUTS:
%
% OUTPUTS:
%       interp_coeff_vec         - interpolation coefficients for the grid
%                                  points inside the mask
%       polar_initial_direction -  optimal polar initial direction
%                                  unit vector from emitter to the receiver
%       cartesian_pos_endpoint   - the cartesian position of the
%                                  end point of the ray
%       num_rays                 - the number of rays traced for solving the ray
%                                  linking inverse problems

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 

% trace the first ray, calculate the end point and its residual, and
% associated interpolation coefficients

if abs(res_initial_guess(1) - res_initial_guess(2))< varepsilon
    
    % trace a single ray for calculation of the interpolation coefficient.
    % avoid doing ray linking
    [~, interp_coeff_vec, cartesian_pos_endpoint] = ...
        feval(solve_ray, polar_direction_initial_guess(1), [], true);
    
    polar_initial_direction = polar_direction_initial_guess(1);
    
    % count only a single ray
    num_rays = 1;
    
else
    
    
    
    % get the polar initial direction of the rays for the intervals
    polar_initial_direction_left = polar_direction_initial_guess(1);
    polar_initial_direction_right = polar_direction_initial_guess(2);
    
    % get the corresponding residulas for the intervals
    polar_residual_left = res_initial_guess(1);
    polar_residual_right = res_initial_guess(2);
    
    
    % update the polar initial direction
    polar_initial_direction  = (polar_initial_direction_left * polar_residual_right...
        - polar_initial_direction_right * polar_residual_left)/...
        (polar_residual_right - polar_residual_left);
    % update the residual
    [polar_residual, interp_coeff_vec, cartesian_pos_endpoint] = ...
        feval(solve_ray, polar_initial_direction, polar_direction_receiver, true);
    
    % initialise the counter for calculating the number
    % of rays for solving the ray linking inverse
    % problem
    num_rays = 1;
    
    
    if abs(polar_residual) > varepsilon
        % allocate vectors for the initial direction that gives minimum
        % residual during the ray linking
        polar_initial_direction_minimum = polar_initial_direction;
        polar_residual_minimum = polar_residual;
        
        
        
        while abs(polar_residual) > varepsilon  &&   num_rays < max_iter
            
            % update the intervals
            if polar_residual_left * polar_residual < 0
                polar_initial_direction_right = polar_initial_direction;
                polar_residual_right = polar_residual;
            elseif  polar_residual_left * polar_residual > 0
                polar_initial_direction_left = polar_initial_direction;
                polar_residual_left = polar_residual;
            else    
            end
            
            
            % update the polar initial direction
            polar_initial_direction  = (polar_initial_direction_left * polar_residual_right...
                - polar_initial_direction_right * polar_residual_left)/...
                (polar_residual_right - polar_residual_left);
            % update the residual
            [polar_residual, ~, ~] = ...
                feval(solve_ray, polar_initial_direction, polar_direction_receiver, false);
            
            % increase the counter
            num_rays = num_rays + 1;
            
            % replace the minimum residual with the current residual, if the current
            % residual is smaller than the previously calculated minimal
            % residual during the iterative process
            if abs(polar_residual) < abs(polar_residual_minimum)
                polar_initial_direction_minimum = polar_initial_direction;
                polar_residual_minimum = polar_residual;
            end
            
            
        end
        
        
        
        % if the residual does not reach the tolerance after the maximum number of iterations,
        % replace it by the initial direction that has given the minimal residual,
        % and use that for calculating the interpolation coefficients
        if polar_residual > varepsilon
            polar_initial_direction = polar_initial_direction_minimum;
        end
        
        
        % calculate the interpolation coefficients for the linked ray
        [~, interp_coeff_vec, cartesian_pos_endpoint,~] = ...
            feval(solve_ray, polar_initial_direction, [], true);
    end
    
    
    
    
    
    
end


end