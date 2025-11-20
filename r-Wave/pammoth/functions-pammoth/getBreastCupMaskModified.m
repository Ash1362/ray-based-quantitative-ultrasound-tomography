function mask = getBreastCupMaskModified(model, cup_index, z_offset)
%FUNCTIONTEMPLATE is a template for a function describtion
%
% DETAILS: 
%   Compared to the original function 'getBreastCupMask.m', the coefficients
%   are modified such that the cup does not match the real cup, but will be
%   a mask larger than the cup for image reconstruction.
%   functionTemplate.m can be used as a template 
%
% USAGE:
%   x = functionTemplate(y)
%
% INPUTS:
%   y - bla bla
%
% OPTIONAL INPUTS:
%   z    - bla bla
%   para - a struct containing further optional parameters:
%       'a' - parameter a
%
% OUTPUTS:
%   x - bla bla
%
% ABOUT:
%       author          - Felix Lucka
%       date            - 19.12.2019
%       last update     - 19.12.2019
%
% See also
%  % Ashkan Javaherian added a cup size 9 as the binary mask for
%  image reconstruction for all cup sizes.
if(nargin < 3)
    z_offset = 0;
end

f = getBreastCupFunctionModified(cup_index);

[X,Y,Z] = ndgrid(model.x_vec, model.y_vec, model.z_vec - z_offset);
R       = sqrt(X.^2 + Y.^2);
clear X Y
mask    = (Z > -f(R));

end