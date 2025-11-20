function [gradient_x, gradient_y, gradient_z, diffusion_coefficient] = ...
    calcDirectionalGradient(field, field_background,...
    smoothing_parameter, grid_spacing)
%   CALCDIRECTIONALGRADIENT computes the directional gradients of field,
%   and the diffusion coefficient
%   
%
% DESCRIPTION:
%      calcDirectionalGradient computes the directional gradients of a
%      field and the diffusion coefficient. They will be used for computing
%      the gradient of L1 norm of total variation of a field as a part of
%      total variation regularization using a Paeron Malik method. 
%      The total variation of a function is in the form J(x)= ||DX||_1.
%      The gradient will be iteratively computed in the form \nabla J(X) = M[X_b]X,
%      where M[X]= D^T c[X_b] D with c[X_b] is the diffusion coefficient.
%      This fucntion computes the action of c[X_b] D X on 
%
% USAGE:
%     
%
% INPUTS:
%       field             - the field X on which the operator M acts (X)
%       field_background  - the field X_b used for computing the diffusion
%                           coefficient C [X_b]
%       smoothing_parameter - the smoothing paparmeter in the dominator of
%                           the diffusion coefficient
%       grid_spacing      - the grid spacing [m] 

 
% OPTIONAL INPUTS:
%
 
% OUTPUTS:
%       gradient_x       - the directional gradient along the Cartesian
%                          coordinate x 
%       gradient_y       - the directional gradient along the Cartesian
%                          coordinate y
%       gradient_z       - the directional gradient along the Cartesian
%                          coordinate z
%       diffusion coefficient - the first diagonal of a diffusion
%                          coefficient matrix
%       above
% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 18.03.2020
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian


switch ndims(field_background)
    case 2
        
        
        % compute the diffusion coefficient using the background field
        gradient_x = padarray(diff(field_background, 1, 1), [1, 0], 0, 'post')/ grid_spacing;
        gradient_y = padarray(diff(field_background, 1, 2), [0, 1], 0, 'post')/ grid_spacing;
        
        diffusion_coefficient = 1./ (1 + (gradient_x.^2 +...
            + gradient_y.^2)./(smoothing_parameter^2)  );
        
        % compute and overwite the directional gradients using the update direction field
        gradient_x = padarray(diff(field, 1, 1), [1, 0], 0, 'post')/ grid_spacing;
        gradient_y = padarray(diff(field, 1, 2), [0, 1], 0, 'post')/ grid_spacing;
        
        
        gradient_z = 0;
        
    case 3
        
        % compute the diffusion coefficient using the background field
        gradient_x = padarray(diff(field_background, 1, 1), [1, 0, 0] , 0, 'post')/ grid_spacing;
        gradient_y = padarray(diff(field_background, 1, 2), [0, 1, 0] , 0, 'post')/ grid_spacing;
        gradient_z = padarray(diff(field_background, 1, 3), [0, 0, 1] , 0, 'post')/ grid_spacing;
        
        
        diffusion_coefficient = 1./ ( 1 + (gradient_x.^2 +...
            + gradient_y.^2  + gradient_z.^2)./(smoothing_parameter^2) );
        
        % compute and overwite the directional gradients using the update direction field
        gradient_x = padarray(diff(field, 1, 1), [1, 0, 0] , 0, 'post')/ grid_spacing;
        gradient_y = padarray(diff(field, 1, 2), [0, 1, 0] , 0, 'post')/ grid_spacing;
        gradient_z = padarray(diff(field, 1, 3), [0, 0, 1] , 0, 'post')/ grid_spacing;
        
        
end

end