function [jacobian_matrix] = calcRayJacobian(solve_ray, polar_initial_direction,...
    polar_direction_endpoint, tau)
%CALCRAYJACOBIAN calculates the Jacobian matrix for a ray linking inverse
% problem based on Newton's method 
%
% DESCRIPTION:
% calcRayJacobian calculates the Jacobian matrix for solving a ray linking
% inverse problem for a 3D case based on Newton's method. The Jacobian 
% applies perturbation to the two components of the polar initial direction
% of the ray in order to calculate the perturbations of the two components of the
% polar direction of a geometerical vector from emitter to the last point
% of the ray on the detection surface. The perurbation to polar initial
% directions are applied using forward finite difference method.

% USAGE:
%
%
% INPUTS:
%      solve_ray                - the function handle for solving the ray and the
%                                 parameters associated with the end point of the ray
%      polar_initial_direction  - a 2 x 1 vector of the polar initial direction of the
%                                reference ray for which the Jacobian matrix is calculated
%      polar_direction_endpoint - a 2 x 1 vector of the polar direction of a geometrical vector
%                                 from the emitter to last point of the ray
%      tau                      - a scalar controlling perturbations applied to the two 
%                                 two components of the polar initial direction
% OPTIONAL INPUTS:
%
% OUTPUTS:
%       Jacobian                - a 2 x 2 Jacobian matrix representing the perturbations
%                                 applied to the polar direction of a geomterical vector
%                                 from emitter to the last point of the ray from 
%                                 the perturnations applied to the polar initial direction

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 



 n = length(polar_initial_direction);
 jacobian_matrix = zeros(n);
 
 
 for ind = 1:n
     
     % adjust the perturbation
     if polar_initial_direction(ind) == 0
         h = tau;
     else
         h = tau*sign(polar_initial_direction(ind))*max(abs(polar_initial_direction(ind)),...
             norm(polar_initial_direction,1)/n);
     end
     
     % apply perturbation to the polar initial direction
     perturbed_polar_initial_direction = polar_initial_direction;
     perturbed_polar_initial_direction(ind) = polar_initial_direction(ind) + h;

     
     % [~, ~, ~, perturbed_polar_ditrection_endpoint] = feval(solve_ray,...
     %    perturbed_polar_initial_direction, [], false);
     % perturbation =  perturbed_polar_ditrection_endpoint - polar_direction_endpoint;
     
     % the perturbation must be in range [-pi,pi] in order to ovoid
     % discontinities in transitions at -pi and +pi
     %  perturbation(1) = wrapToPi(perturbation(1));
     
     % get the perturbation to the end point of the ray
     perturbation = feval(solve_ray,...
         perturbed_polar_initial_direction, polar_direction_endpoint, false);
     
     % get a column of the Jacobian of the end point of the ray to the
     % initial direction in the polar coordinates
     jacobian_matrix(:, ind)= perturbation/h;
 end

 
end

