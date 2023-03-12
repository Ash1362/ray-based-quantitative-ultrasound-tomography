function [directions] = calcDerivativeAngleToInitialPosition(directions)
%CALCDERIVATIVEANGLETOINITIALPOSITION computes the derivative of the
%direction of the ray to the initial position
%
% DESCRIPTION:
%       calcDerivativeAngleToInitialPosition computes the derivative of
%       the direction of the ray in the polar coordinates, i.e., perturbation to 
%       angle of the ray, with respect to the initial position of the ray. 
%       For each forward (resp. adjoint) ray, this derivative is computed using
%       finite differences and using the the rays initialised from the adjacent 
%       emitters (resp. receivers).

%
% USAGE:
%      
%
% INPUTS:
%       directions        - a 1 x num_emitter (resp. 1 x num_receiver) cell 
%                           each containing the the Cartesian direction of the
%                           forward (resp. adjoint) rays initialised from an emitter
%                           (resp. receiver)
%
% Optional INPUTS:
%              
% 
% OUTPUTS:
%       directions     - a 1 x num_emitter (resp. 1 x num_receiver) cell 
%                           each containing the the Cartesian direction of the
%                           forward (resp. adjoint) rays initialised from an emitter
%                           (resp. receiver), and the derivative of the
%                           angle of the rays with respect to the initial
%                           position of the ray
%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 23.08.2020
%       last update     - 14.01.2023
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian 



directions = cat(3, directions{:});
directions = cat(3, directions(:,:, end), directions, directions(:, :, 1));
directions = cat(2, directions(:,:, 2:end-1), 1/2 * acos(dot(directions(:,:, 3:end), directions(:,:, 1:end-2), 2)));

directions = permute(squeeze(mat2cell(directions, size(directions, 1), size(directions, 2),...
    ones(1, size(directions, 3)))), [2, 1]);
end

