function [ray_position_extrapolated] = addPointsRays(ray_position, receiver_indices,...
    num_extrapolation_points, num_ray)
%ADDPOINTSRAYS adds a number of points to the end of rays using a linear
%extrapolation
%
%     DESCRIPTION:
%       addPointsRays adds a given number of points to the end point of rays using a
%       linear extrapolation. The spacing of the added points equals the ray spacing.
%       This is done for computing the derivative of the last rays' points with respect
%       to the initial angle or with respect to the arc length, because the auxilary rays
%       may have fewer points than the linked (reference) ray, and therefore, the rays'
%       Jacobian on the last points of the linked (reference) ray cannot be computed
%       without increase of points on the auxiliary rays.

%
% USAGE:
%       
%
% INPUTS:
%       ray_position             - the position of the points on the rays
%                                  in a single Cartesian coordinate
%       receiver_indices         - a vector containing the index of columns
%                                  correponding to the reception points in 
%                                  matrices containing the rays' parameters
%       num_extrapolation_points - a vector containing the number of points
%                                  must be added to the last point of each ray
%       num_ray                  - the number of rays
% OUTPUTS:
%       ray_position_extrapolated - the position of the points on the corrected rays
%                                  in a single Cartesian coordinate
%       
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 26.04.2023
%
% This script is part of the r-Wave toolbox 
% Copyright (c) 2022 Ashkan Javaherian

ray_position_extrapolated = ray_position;

for ind_ray = 1:num_ray

    if all(~isnan(ray_position(ind_ray, 1)))
        
    % get the index of the last point on the ray, i.e. the reception point
    % which has been shifted forward for correcting the ray spacing
    receiver_index = receiver_indices(ind_ray);
    
    % add a given number of points beyond the last point using a linear
    % extrapolation
    ray_position_extrapolated(ind_ray,...
        receiver_index:receiver_index + num_extrapolation_points(ind_ray)) =...
        cumsum([ray_position(ind_ray, receiver_index),...
        (ray_position(ind_ray, receiver_index)...
        -ray_position(ind_ray, receiver_index-1))...
        * ones(1, num_extrapolation_points(ind_ray))]);

    end

end



end