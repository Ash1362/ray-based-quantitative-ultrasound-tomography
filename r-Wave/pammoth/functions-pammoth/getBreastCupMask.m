function mask = getBreastCupMask(model, cup_index, z_offset)
%FUNCTIONTEMPLATE is a template for a function describtion
%
% DETAILS: 
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

if(nargin < 3)
    z_offset = 0;
end

f = getBreastCupFunction(cup_index);

[X,Y,Z] = ndgrid(model.x_vec, model.y_vec, model.z_vec - z_offset);
R       = sqrt(X.^2 + Y.^2);
clear X Y
mask    = (Z > -f(R));

end