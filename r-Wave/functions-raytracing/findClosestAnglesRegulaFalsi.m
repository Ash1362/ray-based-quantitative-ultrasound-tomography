function [indices, res_indices] = findClosestAnglesRegulaFalsi(polar_direction_endpoints,...
    polar_direction_receivers, dim, varepsilon)
%FINDCLOSESTANGLESREGULFALSI finds the index of the initial angles closest to the angular
%position of the receiver
%
% DESCRIPTION:
% findClosestAnglesRegulaFalsi chooses two indices of initial angles for each receiver
% that give the closest end points to the receiver with different signs
% in terms of the anglar distance from emitter to the receiver so that
% the angular interval defined by the end points of the ray for each receiver
% contains the angular direction of a geomterical vector from emitter to that receiver
% USAGE:
%
%
% INPUTS:
%       polar_direction_endpoints - a 1 x num_receiver vector containg the
%                                   polar direction of geometrical vectors
%                                   from the emitter to the last points of
%                                   the rays
%       polar_direction_receivers - a 1 x num_receiver vector containing
%                                   the polar direction of geometrical
%                                   vectors from the emitter to the receivers
%                                  to the centre of the detection surface (origin).
%      varepsilon                 - a small scalar used as a tolerance
%      dim                        - the dimesion of the medium

% OPTIONAL INPUTS:

%
% OUTPUTS:
%
%      indices    - a 2 x num_receiver matrix containing the chosen two indices
%                   from the chosen initial directions for initialisation
%                   of the ray linking problem for each receiver
%      res_indices - a 2 x num_receiver matrix containing the discrepancy of the
%                    end points of the rays corresponding to the chosen indices
%                    and the receivers in terms of the polar direction of
%                    unit geomterical vectors from the emitter to those
%                    end points and receivers


% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2020 Ashkan Javaherian 


% give an error if this function is used for 3D case
if dim == 3
    error(['FindClosestAngles.m is only used for solving the ray linking problem'...
        'using Regula Falsi approah in 2D case.']);
end

% the number of receivers
num_receiver  = size(polar_direction_receivers, 2);
% allocate matrices for storing the closest indices and their corresponding
% distances to the receivers in terms of the polar direction from the emitter
indices = zeros(dim, num_receiver);
res_indices  = zeros(dim, num_receiver);

for i_r = 1:num_receiver
    
    % calculate the distance of the interception point to the receiver in terms
    % of the unit polar direction from the emitter
    res = polar_direction_endpoints - polar_direction_receivers(i_r);
    
    % find the closest index and its corresponding polar distance of the
    % interception point to the receiver
    [~,ix1] = min(abs(res));
    res1 = res(ix1);
    
    
    if abs(res1) < varepsilon
        
        % if the interception point is sufficiently close to the receiver
        % within a tolerance, choose the second index the same as the first
        % index.
        ix2 = ix1;
    else
        
        % find the closest index to ix1 for which the angular distance of
        % the end point of the ray to the receiver has different sign
        signs_before_ix1 = flip(sign(res1) .* sign(res(1: ix1-1) ));
        signs_after_ix1  =      sign(res1) .* sign(res(ix1+1: end));
        index_before  = find(signs_before_ix1 < 0, 1, 'first');
        index_after   = find(signs_after_ix1 < 0, 1, 'first');
        
        
        if isempty(index_before) && isempty(index_after)
          %  ix2 = max(1, ix1-1);
          %  ix1 = min(size(polar_direction_endpoints, 2), ix1+1);
            ix1 = 1;
            ix2 = size(polar_direction_endpoints, 2);
   
        elseif isempty(index_after)
            ix2 = ix1 - index_before;
        elseif isempty(index_before)
            ix2 = ix1 + index_after;    
            
        else
            if index_before <= index_after
                ix2 = ix1 - index_before;
            else
                ix2 = ix1 + index_after;
            end
        end
        
    end
    
  
   
    
    % fill the first and second rows of allocated matrices for the indices and
    % thier angular distances with the chosen index with negative and positive
    % angular distances, respectively
    if sign(res1) <= 0
        indices(:, i_r) = [ix1; ix2];
        res_indices(:, i_r) = [res(ix1); res(ix2)];
    else
        indices(:, i_r) = [ix2; ix1];
        res_indices(:, i_r) = [res(ix2) ; res(ix1)];
    end
    
end



end

