function [pos_intersection, distance] = calcIntersectBall(pos_start,...
    direction, radius)
% CALCINTERSECTBALL calculates the intersection of a line with
% a Ball.

% DESCRIPTION:
% calcIntersectBall calculates the intersection a line with a circular
% (spherial) detection surface.      
%
% USAGE:
%      [rs] = solveIntersectionLineBall(r,ds,R)
%
% INPUTS:
%   
%       pos_start - dim x 1 cartesian position of the starting point of the line 
%       direction - dim x 1direction of the line 
%       radius    - radius of the ball
%
% OPTIONAL INPUTS:
%        
%
% OUTPUTS:
%       pos_intersection  - dim x 1 cartesian position of the intersection 
%                           point of the line with the ball
%                           (detection circle or hemi-sphere)
%       d                 - distance between the starting and end point


% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 

% make the direction a unit vector
direction = direction / norm(direction);

% calculate the coefficients for solving a second-order
% equation for the the intersection
b    =   direction' * pos_start;
c    =   pos_start' * pos_start- radius^2;
crit =   b^2 - c;



% calculate the distance of the intersection point from the starting point
% of the line
if crit < 0
    
    error (['The line does not intersect with the Ball.'... 
     'The initial angle of the ray is incorrect.']);
    
else
    
    if crit == 0
        distance = - b;
    else       
d1 = -b + sqrt(crit);
d2 = -b - sqrt(crit);
[~,ix] = min(abs([d1, d2]));

if ix == 1
    distance = d1;
else
    distance = d2;   
end

    end
    
end

% calculate the intersection point 
pos_intersection = pos_start + distance * direction;



end