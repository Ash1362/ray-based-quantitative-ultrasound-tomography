function [gradient_total_variation] = calcDivergence(field_x, field_y, field_z,...
     mask, grid_spacing)
%   CALCDIVERGENCE computes the action of a divergence operator on
%   a field
%   
%
% DESCRIPTION:
%      calcDivergence computes the action of a divergece operator on a field.
%      That will be used for computing the gradient of L1 norm of total variation
%      of a field as a part of total variation regularization using a Paeron Malik method. 
%      The total variation of a function is in the form J(x)= ||DX||_1.
%      The gradient will be in the form \nabla J(X) = M[X_b]X, where 
%      M[X]= D^T c[X_b] D X with c[X_b] is the diffusion coefficient.
%      This fucntion computes the action of D^T. 

%
% USAGE:
%     
%
% INPUTS:
%       field_x          - the directional gradient along the Cartesian
%                          coordinate x 
%       field_y          - the directional gradient along the Cartesian
%                          coordinate y
%       field_z          - the directional gradient along the Cartesian
%                          coordinate z
%       mask             - the binary mask for image reconstruction
%       grid_spacing     - the grid spacing [m]


 
% OPTIONAL INPUTS:
%
 
% OUTPUTS:
%       gradient_total_variation - the action of the divergence operator on the firld

% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 18.03.2020
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian


switch ndims(field_x)
    
    case 2
               
%spatial_variation_x = cat(1, field_x(1, :), diff(field_x(1:end-1, :), 1, 1), - field_x(end, :));
%spatial_variation_y = cat(2, field_y(:, 1), diff(field_y(:, 1:end-1), 1, 2), - field_y(:, end));

spatial_variation_x = cat(1, field_x(1, :), diff(field_x, 1, 1) );
spatial_variation_y = cat(2, field_y(:, 1), diff(field_y, 1, 2) );

gradient_total_variation = - (spatial_variation_x(mask) + spatial_variation_y(mask))/grid_spacing;


    case 3
        
%spatial_variation_x = cat(1, field_x(1, :, :), diff(field_x(1:end-1, :, :), 1, 1), - field_x(end, :, :));
%spatial_variation_y = cat(2, field_y(:, 1, :), diff(field_y(:, 1:end-1, :), 1, 2), - field_y(:, end, :));
%spatial_variation_z = cat(3, field_z(:, :, 1), diff(field_z(:, :, 1:end-1), 1, 3), - field_z(:, :, end));

spatial_variation_x = cat(1, field_x(1, :, :), diff(field_x, 1, 1));
spatial_variation_y = cat(2, field_y(:, 1, :), diff(field_y, 1, 2));
spatial_variation_z = cat(3, field_z(:, :, 1), diff(field_z, 1, 3));

gradient_total_variation = - (spatial_variation_x(mask) + spatial_variation_y(mask) + spatial_variation_z(mask))/grid_spacing;

end



end