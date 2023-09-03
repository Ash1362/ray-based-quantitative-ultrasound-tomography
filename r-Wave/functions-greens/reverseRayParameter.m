function [ray_parameter_adjoint] = reverseRayParameter(ray_parameter, interp_coeff,...
    receiver_indices, nan_binary, parameter_binary)
%REVERSERAYPARAMETER reverses the position and accumulated parameters
%along the linked forward rays
%
% DESCRIPTION:
%       reverseRayParameter reverses the Cartesian positions and the accumulated parameres
%       along the forwards rays initilised from an emitter for computing the adjoint rays
%       starting from the receivers and intecepted by the amitter after
%       travelling through the medium inside the detection surface (ring)

%
% USAGE:
%      
%
% INPUTS:
%       ray_parameter     - a num_receiver x num_maxraypoints positions
%                           (in a single Cartesian coordinate), or parmeters
%                           on the sampled points along the forward rays linking an
%                           emitter to all receivers
%       interp_coeff      - a num_receiver x 1 vector of the interpolation 
%                           coefficients for all the linked rays
%       receiver_indices  - a num_receiver x 1 vector of index of the
%                           reception (last) points on the linked rays 
%       nan_binary        - the binary matrix indicating nans in the matrix of 
%                           positions (parmeters) of the points on the
%                           linked rays
%       parameter_binary  - a boolean determimh whether the input 'ray_parameter'
%                           contains the Cartesian positions, or
%                           accumulated parmeters
%         
% Optional INPUTS:
%              
% 
% OUTPUTS:
%       ray_parameter_adjoint - a num_receiver x num_maxraypoints positions
%                           (in a single Cartesian coordinate), or parmeters
%                           on the sampled points along the adjoint rays linking 
%                           all receivers to the emitter

%       
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian 

     
% allocate a zero matrix for the adjoint rays with the same size as the
% matrix for the forward rays
ray_parameter_adjoint = zeros(size(ray_parameter));

% get the number of receivers
num_receiver = length(receiver_indices);

% fil-in the matrix with the nans the same as for the matrix for the
% forward rays
ray_parameter_adjoint(nan_binary) = nan;

% interpolate the ray parameters from the forward rays to the adjoint
% rays using the given interpolation cioefficients
ray_parameter_adjoint(:, 2:end) = interp_coeff.* ray_parameter(:, 2:end)...
    + (1 - interp_coeff) .* ray_parameter(:, 1:end-1);

% make the ray parameter for the first points for the adjoint matrix the same as
% those for the forward matrix
ray_parameter_adjoint(:, 1) = ray_parameter(:, 1);


% get the linear indices for the receivers
receiver_indices_linear = sub2ind(size(ray_parameter), (1:num_receiver)', receiver_indices);


% make the parameters for the reception points in the adjoint matrix the same as
% those in the forward matrix
ray_parameter_adjoint(receiver_indices_linear) = ray_parameter(receiver_indices_linear);


% check to see if all the interpolated parameters are finite
if any(~isfinite(ray_parameter_adjoint)-nan_binary)
    error(['The interpolated parameters onto the adjoints rays contain infinite values.']);
end

% reverse the Cartesian positions (parameters) on the points along each linked ray
% such that they are initiliased from receivers and ended on the emitter
    for ind_receiver = 1 : num_receiver
        ray_parameter_adjoint(ind_receiver, 1:receiver_indices(ind_receiver)) =...
            flip(ray_parameter_adjoint(ind_receiver, 1:receiver_indices(ind_receiver)));
    end
   
    % If an accumulated parameter is being reversed
    if parameter_binary
        ray_parameter_adjoint = ray_parameter_adjoint(:,1)...
            - ray_parameter_adjoint;
    end

  
end