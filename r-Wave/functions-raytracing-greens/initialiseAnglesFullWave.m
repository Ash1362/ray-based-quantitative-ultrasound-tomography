function [polar_directions_initial_guess, res_initial_guess] = ...
    initialiseAnglesFullWave(solve_ray, polar_direction_receivers, dim, interval_factor,...
    varepsilon, raylinking_method, reference_ray)
%INITIALISEANGLESFullWave chooses initial directions closest to the receiver
%
% DESCRIPTION:
% initialiseAnglesFullWave initialises the ray linking problem simulataneously for all
% receivers using a 'Global' approach. It chooses a set of initial directons
% giving interception points (end points for rays) closest to each receiver
% in terms of the polar direction of a geomtrical vector from emitter to
% those points (and receiver)
% Alternatively using 'Regula Falsi' method, which is only used for 2D case,
% two initial directions giving interception points close to and containing
% the receiver are chosen.
% Using 'Regula Falsi' method, the polar directions of the interception points
% (and the receuiver) are angles between geomterical vectors from the emitter to
% those interception points and another reference vector from the emitter
% to the centre of the detection surface, which is origin of the cartesian coordinate
% by convention
% USAGE:
%
%
% INPUTS:
%       solve_ray                 - a function for solving a ray continaing the
%                                   the required parameters
%       polar_direction_receivers - a 1 x num_receiver vector containing
%                                  the polar direction of geometrical
%                                  vectors from the emitter to the receivers
%      dim                       - the dimension of the medium
%      interval_factor           - a coefficient controlling the distance
%                                   between the initial polar directions,
%                                   from which the closest indices for each
%                                   receiver are chosen
%      varepsilon               - a small scalar used as a tolerance
%      raylinking_method         - method for ray linking
% OPTIONAL INPUTS:

%
% OUTPUTS:
%      polar_directions_initial_guess - a dim-1 x num_receiver vector of the initial guess 
%                                      for the polar initial directions of the rays  
%      res_initial_guess         - If given, the dim-1 x num_receiver residual obtained 
%                                 for the end point of the rays initialised by the initial guess

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian

% choose a set of initial polar directions for the rays
angles = polar_direction_receivers(:, 1:interval_factor:end);
% For 2D case, using the 'Regula Falsi' method, sort angles from -pi/2 to
% +pi/2 in order to avoid the discontinuity (the jump from -pi/2 to + pi/2)
% this set of anlges will be used as initial guess for solving the ray linking problem
if strcmp(raylinking_method, 'Regula-Falsi')
    switch dim
        case 2
            [angles_ascend, ~ ] = sort(angles);
        case 3
            error('Regula Falsi method is not used for 2D case.');
    end
else
    angles_ascend = angles;
end


% number of angles
num_angles = size(angles, 2);


% trace a set of rays using the choosen initial directions, and calculate the
% polar directions of the ened point of the rays (interception points).
% Using 'Regula Falsi' method, the polar ditrections will be angles between
% the geomterical vectors from the emitter to the interception points and a
% reference vector from the emitter to the centre of the detection surface
% (the origin of the cartesian coordinate)


% allocate a dim-1 x num_angle matrix for storing the polar direction of the
% interception points of the rays solved by the chosen initial directions
polar_direction_endpoints = zeros(dim-1, num_angles);


for i = 1 : num_angles
    [~, ~, ~, polar_direction_endpoints(:, i)] = feval(solve_ray,...
        angles_ascend(:, i), [], false, reference_ray);
end

% For each receiver, choose the polar initial directions that give the
% closest endpoints to the receivers in terms of the polar directions from the emitter
if strcmp(raylinking_method, 'Regula-Falsi')
    
    % Using 'Regula Falsi' method, for each receiver, choose a pair of initial
    % angles that give the closest end points to the receiver with different
    % sign in terms of the anglar distance
    [indices_initial_guess, res_initial_guess] = findClosestAnglesRegulaFalsi(...
        polar_direction_endpoints, polar_direction_receivers, dim, varepsilon);
    
else
    [indices_initial_guess, res_initial_guess] = findClosestAngles(polar_direction_endpoints,...
        polar_direction_receivers, dim);
end

% choose the initial direction for initialisation of the ray linking
% inverse problem using a 'Global' approach
polar_directions_initial_guess = angles_ascend(indices_initial_guess);



end