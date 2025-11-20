function [field_smoothed] = smoothField(field, num_points, value_background)
%SMOOTHFIELD applies smoothing on a field
%
% DESCRIPTION:
% smoothField applies smoothing on a field
%
%
% USAGE:
%
%
% INPUTS:
%       field      - a field
%       num_points - the number of points along each Cartesian coordinate used for averaging 
%
% OPTIONAL INPUTS:
%
% OUTPUTS:
%       field_amoothed  - the smoothed field
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 

if nargin < 3
    value_background = 1;
end




if rem(num_points, 2) == 0
    error(' The smoothing window size must be an odd natural scalar.');
end
% dispaly the chosen smoothing window size, an odd natural scalar
disp(['The smoothing window size is ' num2str(num_points)]);

% the half size of the averaging window
num_halfpoints = (1/2) * (num_points-1);

field_smoothed = field;

 switch length(size(field))
     case 2
         field_smoothed = 1/(num_points^2) * convn(field, true(num_points, num_points), 'same');
         field_smoothed([1:num_halfpoints, end-num_halfpoints+1:end], : ) = value_background;
         field_smoothed( :, [1:num_halfpoints, end-num_halfpoints+1:end]) = value_background;
     case 3
         field_smoothed = 1/(num_points^3) * convn(field, true(num_points, num_points, num_points), 'same');
         field_smoothed([1:num_halfpoints, end-num_halfpoints+1:end], :, : ) = value_background;
         field_smoothed(:, [1:num_halfpoints, end-num_halfpoints+1:end], :) = value_background;
         field_smoothed(:,:, [1:num_halfpoints, end-num_halfpoints+1:end]) = value_background;
 end
  

