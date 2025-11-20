function [ray_position_adjoint] = calcRayAdjointParaxial(ray_position, ray_spacing,...
    rayspacing_receiver)
%CALCRAYSADJOINTPARAXIAL computes the auxiliary adjoint ray using the
%paraxial approach
%
% DESCRIPTION:
%       calcRaysAdjointParaxial computes the auxiliary adjoint ray using
%       the paraxial approach. The main adjoint ray is not traced. Instead,
%       the psoition of the ray points along the main forward ray is
%       reversed, and is used for computing the adjoint auxiliary ray.
%
% USAGE:
%      
%
% INPUTS:
%       ray_position     - a dim x num_raymaxpoints Cartesian position of the points along
%                          the forward main (linked) rays. Here, dim is the dimension 
%                          of the medium, and num_maxraypoints is a scalar value representing
%                          the maximum permissible number of ray points for linking the 
%                          emitter-receiver pairs.
%       ray_spacing         - a scalar value representing the spacing [m] along the rays
%       rayspacing_receiver - a vector containing the last sapcing along the linked rays,
%                             i.e. the spacing between the reception point and the point
%                             just before that
%
% Optional INPUTS:
%              
% 
% OUTPUTS:               - a dim x num_raymaxpoints Cartesian position of the points along
%                          the reversed ray.
%       

%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 27.11.2019
%       last update     - 14.01.2023
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian

% get the number of dimensions
dim = size(ray_position, 1);

% get the nans
nan_binary = isnan(ray_position(1,:));

% get the index of point on the receiver
receiver_index = find(~isnan(ray_position(1,:)), 1, 'last');

% reverse the x cartesian positions of the sampled points on the linked ray
ray_position_x_adjoint = reverseRayParameter(ray_position(1, :),...
   rayspacing_receiver/ray_spacing, receiver_index, nan_binary, false);

% reverse the y cartesian position of the sampled points on the linked ray
ray_position_y_adjoint = reverseRayParameter(ray_position(2, :),...
   rayspacing_receiver/ray_spacing, receiver_index, nan_binary, false);

% make a matrix containing the reversed cartesian position of the rays
% along x and y coordinates
ray_position_adjoint = [ray_position_x_adjoint;...
    ray_position_y_adjoint];


if dim == 3
           
% reverse the z cartesian positions 
ray_position_z_adjoint = reverseRayParameter(ray_position(3, :),...
   rayspacing_receiver/ray_spacing, receiver_index, nan_binary, false);
          
% add the reversed z cartesian position to the rays' positions
ray_position_adjoint = [ray_position_adjoint;...
    ray_position_z_adjoint];
       
end



end

