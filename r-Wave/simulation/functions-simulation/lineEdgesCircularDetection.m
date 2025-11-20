function [start_edge_pos,end_edge_pos, slope] = lineEdgesCircularDetection(tr_pos, radius )
%LINEEDGESCIRCULARDETECTION calculates the starting and ending points of line transducers 
%
% DESCRIPTION:
%          lineEdgesCircularDetection calculates the starting and ending points of a set of line
%          segments that are tangent to a circle. These lines define the geometry of
%          a set of line transducers along a circular detection surface
%          (line)
%
% USAGE:
%         
%
% INPUTS:
%       tr_pos         - 2/3 x N array of cartesian points containing the centers of
%                        the transducers
%       radius         - the half-length of the line transducers 

% OPTIONAL INPUTS:
%
%
% OUTPUTS:
%       start_edge_pos -  the position of the starting point of the line segments
%                (transducers) in Cartesian coordinates [m]               
%       end_edge_pos   -  the position of the ending point of the line segments
%                (transducers) in Cartesian coordinates [m].
%       slope          -  the slope of geometrical vectors starting from the centre of the
%                (circular) detection line to the centre of the line
%                segments (transducers)

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 15.12.2019
%       last update     - 15.12.2019
%
% This function is part of the r-Wave Toolbox 
% Copyright (C) 2020 Ashkan Javaherian 

if size(tr_pos, 1)> 2
    error('A Line tranducer cannot be used for 3D geometry.')
end
slope   = - tr_pos(1,:) ./ tr_pos(2,:);
dlx = radius ./ sqrt(slope .^ 2 + 1) ;
dly = slope .* dlx;

% if dlx (resp. dly) is parallel to the Cartesian Coordinates,
% set dlx = radius;
dlx(isnan(dlx))= radius;
dly(isnan(dly))= radius;


start_edge_pos = cat(1, tr_pos(1, :) - dlx, tr_pos(2, :) - dly);    
end_edge_pos   = cat(1, tr_pos(1, :) + dlx, tr_pos(2, :) + dly); 

% make the order of the calculated edges of the lines transducers 
% anti-clock wise with respect to the origin for all transducers
% the order of the two edges of each line is calculated using the
% the third component of the cross product of the associated geometrical 
% vectors strating from the origin that are extended to three dimensions
start_edge_padded = padarray(start_edge_pos, [1 0], 'post');
end_edge_padded   = padarray(end_edge_pos  , [1 0], 'post');
% for each line, make the direction from the starting point to the end point 
% counter-clockwise
cross_products = cross(start_edge_padded, end_edge_padded);
% the lines for which the start and end points must be reversed
reverse_indices = cross_products(3,:) < 0;
clear start_edge_padded end_edge_padded cross_products
% swap the starting and ending points, if needed
[end_edge_pos(:,reverse_indices), start_edge_pos(:,reverse_indices)] = deal(start_edge_pos(:,reverse_indices), end_edge_pos(:,reverse_indices));

    
end