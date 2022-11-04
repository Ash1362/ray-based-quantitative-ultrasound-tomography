function [ray_position_extrapolated] = addPointsRays(ray_position, receiver_indices,...
    num_extrapolation_points, num_ray)
%ADDPOINTSRAYS adds a number of points to the rays using a linear
%extrapolation
%
%     DESCRIPTION:
%       addPointsRays add a given number of points to the rays using a
%       linear extrapolation. The spacing of the added points equals the 
%       ray spacing. This is required for computing the derivative of 
%       the last rays' points with respect to the initial angle, because 
%       the auxilary rays may not be large enough in arc length (may be shorter
%       than the main rays), and must be increased in arc length. 

%
% USAGE:
%       
%
% INPUTS:
%       ray_position             - the position of the points on the rays
%                                  in a single Cartesian coordinate
%       receiver_indices         - a vector containing the index of columns
%                                 correponding to the reception points in 
%                                 matrices containing the rays' parameters,
%                                 eg. ray_poition_x.... 
%       num_extrapolation_points - a vector containing the number of points
%                                  to be added to each ray
%       num_ray                  - the number of rays
% OUTPUTS:
%       ray_position_extrapolated - the position of the points on the corrected rays
%                                  in a single Cartesian coordinate
%       
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 30.03.2020
%
% This script is part of the r-Wave Tool-box (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian

ray_position_extrapolated = ray_position;

for ind_ray = 1:num_ray
    
    % get the index of the last point on the ray (the point associated with
    % the reception point that has been shifted for correcting the ray spacing)
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