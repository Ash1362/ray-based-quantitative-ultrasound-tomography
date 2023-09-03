function [pos_intersection, distance] = calcLineIntersect(pos_start,...
    direction, detec_geom)
% CALCLINEINTERSECT calculates the intersection point of a straight line
% with a line or surface

% DESCRIPTION:
% calcLineIntersect calculates the intersection point of a straight line
% with a line or surface
%
% USAGE:
%
%
% INPUTS:
%
%       pos_start - dim x 1 cartesian position of the starting point of the line
%       direction - dim x 1 direction of the line
%       detec_geom       - a struct defining the geometry of the detection
%                          surface with fields:
%      'radius_circle'   - the radius [m] of the circular (2D) or
%                          hemi-spherical (3D) detection surface, or
%      'radius_cylinder' - the radius [m] of the cylinder in x-y plane
%      'slope'           - the slope of the linear (2D) or the planar
%                          detection surface (3D) in x-y plane
%      'y_intercept'     - the y interception of the linear detection
%                          surface or the planar detection surface in x-y
%                          plane, or
%                          a char wich can be either '1', or '2', and
%                          determines the scenario for testing the ray tracing
%                          algorithm using a Maxwell fish-eye lens phantom. The
%                          scenarios '1' and '2' is for measuring the accuracy
%                          in computing the acoustic length, and deviation of the rays'
%                          trajectory from an expected circular trajectory,
%                          respectively.
%
% OPTIONAL INPUTS:
%
%
% OUTPUTS:
%       pos_intersection  - dim x 1 cartesian position of the intersection
%                           point of the straight line with the line or
%                           surface
%       distance          - distance between the starting and end point


% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 14.07.2019
%       last update     - 14.02.2023
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian

% get the number of dimensions
dim = length(pos_start);

% make the direction a unit vector
direction = direction / norm(direction);

if isfield(detec_geom, 'radius_circle')   ||  isfield(detec_geom, 'radius_cylinder')
    
    if isfield(detec_geom, 'radius_circle')
        dim1 = dim;
    elseif isfield(detec_geom, 'radius_cylinder')
        dim1 = dim-1;
    end
    
    % compute the coefficients associated with a second-orderequation for the intersection
    b = direction(1:dim1)' * pos_start(1:dim1);
    c = pos_start(1:dim1)' * pos_start(1:dim1)- detec_geom.radius_circle^2;
    crit = b^2 - c;
    
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
    
else
    
    if  isfield(detec_geom, 'line_coeff')
        
        % calculate the intersection point of a line with another line or plane
        % using the Cramer's rule
        switch dim
            case 2
                
                % get the coefficients [a,b, c] for equation ax+by=c
                % of the intersecting straight line for forming the Cramer
                % matrix
                % see https://byjus.com/maths/equation-line/
                cramer_matrix = [1/direction(1),-1/direction(2), [1/direction(1),-1/direction(2)] * pos_start;...
                    detec_geom.line_coeff];
                
                % get the intersection point
                % see https://www.cuemath.com/algebra/cramers-rule/
                
                % allocate a zero vector for the intersection point
                pos_intersection = zeros(dim, 1);
                
                % allocate a dim x 1 vector for the intersection point
                pos_intersection(1) = det(cramer_matrix(:, [3,2]))/det(cramer_matrix(:, 1:2));
                
                % get the y position of the intersection point
                pos_intersection(2) = det(cramer_matrix(:, [1,3]))/det(cramer_matrix(:, 1:2));
                
                
            case 3
                error('Not implemented yet.')
        end
        
        % compute the distance of the intersection point and the
        % starting position of the intersecting line
        distance = norm(pos_intersection - pos_start);
        
    else
        
        
    end
end


end