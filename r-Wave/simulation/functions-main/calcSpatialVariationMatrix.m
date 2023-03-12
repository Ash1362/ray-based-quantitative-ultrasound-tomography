function [gradient_total_variation] = calcSpatialVariationMatrix(field,...
    field_background, mask, grid_size, smoothing_parameter, gamma, grid_spacing)
%   CALCSPATIALVARIATIONMATRIX computes the gradient of an L1 norm of the total
%   variation function.
%   
%
% DESCRIPTION:
%      calcSpatialVariationMatrix computes the gradient of an L1 norm of the total
%      variation function, when a total varaiation reegularisation using a
%      Paeron Malik method is used. 
%      The total variation of a function is in the form J(x)= ||DX||_1.
%      The gradient will be in the form \nabla J(X) = M[X]X, where 
%      M[X]= D^T c[X_b] D X with c[X_b] the diffusion coefficient.
%      Here D gives the directional gradients, and its adjoint D^T is
%      a divergnece operator.

%
% USAGE:
%     
%
% INPUTS:
%       field            - the field X on which the directional gradient
%                          operator D acts
%       field_background - the fixed field x_b used for computing the diffusion
%                          coefficient C
%       mask             - the binary mask for image reconstruction
%       grid_size        - the number of grid points along the Cartesian coordinates
%       smoothing_parameter - the smoothing parameter in the dominator of
%                           the diffusion coefficient digonal matrix
%       gamma             - parameter for ahifting the eigenvalues of the matrix M
%                           for making the matrix M positive definite
%       grid_spacing      - the grid spacing [m]

 
% OPTIONAL INPUTS:
%
 
% OUTPUTS:
%       gradient_total_variation - the action of M[X_b] on X

% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 18.03.2020
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian


% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 18.03.2020
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian


error(['The preconditioning approach is under development. Please set this'...
   'option false for now!'])

% get the matrix of the sound speed used for computing the diffusion
% coefficient
field_background_grid = reshape(field_background, grid_size);

% get the matrix of the update direction field on which the directional
% gradients act
field_direction_grid = zeros(grid_size);
field_direction_grid(mask) = field; 

% compute the directional gradients and diffusio coefficient
[gradient_x, gradient_y, gradient_z, diffusion_coeff] = ...
    calcDirectionalGradient(field_direction_grid, field_background_grid,...
    smoothing_parameter, grid_spacing);

% enforce the divergence operator 
[gradient_total_variation] = - calcDivergence(diffusion_coeff .* gradient_x,...
    diffusion_coeff .* gradient_y, diffusion_coeff .* gradient_z, mask, grid_spacing) +...
    gamma * field;
end