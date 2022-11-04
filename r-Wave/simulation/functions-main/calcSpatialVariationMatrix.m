function [spatial_variation] = calcSpatialVariationMatrix(field_direction,...
    field_background, mask, grid_size, smoothing_parameter, gamma, grid_spacing)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% get the matrix of the sound speed used for computing the diffusion
% coefficient
field_background_grid = reshape(field_background, grid_size);

% get the matrix of the update direction field on which the directional
% gradients act
field_direction_grid = zeros(grid_size);
field_direction_grid(mask) = field_direction; 

% compute the directional gradients and diffusio coefficient
[gradient_x, gradient_y, gradient_z, diffusion_coeff] = ...
    calcDirectionalGradient(field_direction_grid, field_background_grid,...
    smoothing_parameter, grid_spacing);

% enforce the divergence operator 
[spatial_variation] = - calcDivergence(diffusion_coeff .* gradient_x,...
    diffusion_coeff .* gradient_y, diffusion_coeff .* gradient_z, mask, grid_spacing) +...
    gamma * field_direction;
end