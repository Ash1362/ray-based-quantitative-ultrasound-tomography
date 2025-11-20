function [angle] = calcDirectionalAngle(cartesian_direction_first, ...
    cartesian_direction_last)
%CALCDIRECTIONALANGLE calculates the directional angle between two sets of 
%geometrical vectors at their intersection point
%
% DESCRIPTION:
% calDirectionalAngle calculates the cartesian directional anglea between two
% sets of geoemtrical vectors at their intersection points. The angle is calculated in 
% radians between [-pi,+pi] in a counterclockwise direction from the first
% to the last vector. If two (mathematical) vectors are given, the angles
% are computed componentwise.

% INPUTS:
%       cartesian_direction_first - dim x n Cartesian direction of the first
%                                   geometric vectors
%       cartesian_direction_last  - dim x n Cartesian direction of the last
%                                   geometric vectors 
% OPTIONAL INPUTS:

%
% OUTPUTS:    
%       angle           - the angle between geometrical vectors

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 

% get the size of the first array
sz = size(cartesian_direction_first);

% compute the angles between the two-dimensional vectors included in each array
% and reshape the resulting vector of angles to the size of the original array
angle = atan2(cartesian_direction_first(1,:) .* cartesian_direction_last(2,:)... 
            - cartesian_direction_first(2,:) .* cartesian_direction_last(1,:),...
              cartesian_direction_first(1,:) .* cartesian_direction_last(1,:)... 
            + cartesian_direction_first(2,:) .* cartesian_direction_last(2,:));

if length(sz) > 2

% make the fisrt dimension one
sz(1) = 1;

% make the dimension of the nagle the same az size of the the originla
% array, but make the fisrt dimension one.
angle = reshape(angle, sz);

end

end