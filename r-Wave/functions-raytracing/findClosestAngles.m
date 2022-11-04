function [indices, res_indices] = findClosestAngles(polar_direction_endpoints,...
    polar_direction_receivers, dim)
%FINDCLOSESTANGLESREGULAFALSI finds the index of the initial angles closest to the angular
%position of the receiver
%
% DESCRIPTION:
% findClosestAnglesRegulaFalsi choose the index of initial direction for each receiver
% that give the closest end points to the receiver in terms of the anglar distance
% of a geometerical vector from emitter to the receiver
% this is used for an initialisation of the ray linking invers problem
% smilutaneously for all receivers using a 'Global' approach
% USAGE:
%
%
% INPUTS:
%      solve_ray                  - a function for solving a ray continaing the
%                                  the required parameters
%       polar_direction_endpoints - a dim-1 x num_receiver matrix containg the
%                                   polar direction of geometrical vectors
%                                   from the emitter to the last points of
%                                   the rays
%       polar_direction_receivers - a dim-1 x num_receiver matrix containing
%                                  the polar direction of geometrical
%                                  vectors from the emitter to the receivers
%      dim                        - the dimesion of the medium

% OPTIONAL INPUTS:

%
% OUTPUTS:
%
%      indices    - a 1 x num_receiver vector containing the chosen indices
%                   from a set of chosen initial directions for initialisation
%                   of the ray linking problem for each receiver
%      res_indices - a dim-1 x num_receiver matrix containing the distance of the
%                    end points of the rays corresponding to the chosen indices
%                    to the receivers in term of the polar directions from
%                    the emitter


% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2020 Ashkan Javaherian 



% the number of receivers
num_receiver  = size(polar_direction_receivers, 2);
% allocate matrices for storing the closest indices and their corresponding
% residual, the polar distances between the end point of the rays and receivers
% in terms of the polar direction from the emitter to those points
indices      = zeros(1, num_receiver );
res_indices  = zeros(dim-1, num_receiver );

for i_r = 1:num_receiver
    
    % calculate the distance of the interception point to the receiver in terms
    % of the unit polar direction from the emitter to those points
    polar_residual = polar_direction_endpoints - polar_direction_receivers(:, i_r);
    
    % find the closest index and its corresponding residual polar distance between the
    % end point of the rays and the receivers). this is used as an initial guess for ray linking
    [~,  ix] = min(vecnorm(polar_residual));
    res = polar_residual(:, ix);
    indices(i_r) = ix;
    res_indices(:, i_r) = res;
    
end


end