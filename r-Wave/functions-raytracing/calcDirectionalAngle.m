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
% Copyright (c) 2020 Ashkan Javaherian 


angle = atan2(cartesian_direction_first(1,:) .* cartesian_direction_last(2,:)... 
            - cartesian_direction_first(2,:) .* cartesian_direction_last(1,:),...
              cartesian_direction_first(1,:) .* cartesian_direction_last(1,:)... 
            + cartesian_direction_first(2,:) .* cartesian_direction_last(2,:));
                   

end