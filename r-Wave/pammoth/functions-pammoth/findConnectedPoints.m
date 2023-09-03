function [connected_points, num_connected_points] = findConnectedPoints(positions)
%FINDCONNECTEDpointS finds the neighboring points in a set of points
%
% DESCRIPTION:
%     findConnectedpoints finds the neighboring points in a set of
%     points on a surface using triangulation
%
% USAGE:
%     
%
% INPUTS:
%     positions     - the position of points
%
% OUTPUTS:
%     connected_points  - a cell array containing the index of the connected points to each points
%     num_connected_points  - a vector containg the index of the connected points to each points
%
% ABOUT:
%     author        - Ashkan Javaherian
%     date          - 05.08.2020
%     last update   - 12.12.2020
%
%
% This function is part of the r-Wave Toolbox.
% Copyright (C) 2021 Ashkan Javaherian
%

% the number of points 
num_point = size(positions, 2);
% apply a triangulation to the first two Cartesian coordinates of the positions
triangles = delaunayTriangulation(positions(1:2, :)');
% display the triangulation
figure;trimesh(triangles.ConnectivityList, triangles.Points(:,1), triangles.Points(:,2));

v = vertexAttachments(triangles);
connected_points= cell(1, num_point);
num_connected_points = zeros(1, num_point);
for ind_point = 1:num_point
  connected_points{ind_point} = unique(triangles.ConnectivityList(v{ind_point}, :));
   num_connected_points(ind_point) = size(connected_points{ind_point}, 1);
end

end

