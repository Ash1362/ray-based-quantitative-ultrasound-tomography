function [f, d_f, dd_f, e2s] = getBreastCupFunctionModified(cup_index)
%GETBREASTCUPFUNCTIONMODIFIED returns a modified version of the 1D function mapping axial displacement to
%breast cup height. 
%
% DETAILS: 
%   Compared to the original function 'getBreastCupFunction.m', the coefficients are modified 
%   such that the cup does not match the real cup, but will be a mask
%   larger than the cup for image reconstruction.
% The cup shape is found by rotation of a curve f(x) defined for positive 
% x around the z-axis. The curve is composed of an ellipse, a cubic 
% polynomial and a horizontal line:
%   - the ellipse is determined by 'radii' (a, b) and defines the curve f
%     between 0 and e2s, the point where the cubic polynomial starts. The
%     point e2s is determined by f(x) = h
%   - the cubic polynomial defines the curve f between e2s and w and is 
%     determined by the requirements:
%        - the polynomial and its first derivative are continous at e2s
%        - the polynomial is 0 in x = w
%        - the polynomia's first derivative in x = w is -s
%
% USAGE:
%   cup_fun = getBreastCupFunction(1)
%
% INPUTS:
%   cup_index - number of the breast cup size (from 1 to 8)
%
% OUTPUTS:
%   cup_fun    - function handle mapping distance from the z-axis (x=y=0) to 
%                breast cup height
%   d_cup_fun  - first derivative of the above function
%   dd_cup_fun - first derivative of the above function
%         e2s  - x value for which cubic and ellipse join
%
% ABOUT:
%       author          - Felix Lucka, Ashkan Javaherian
%                        
%
%       date            - 25.11.2019
%       last update     - 25.11.2019
% See also

% table with the cup parameters
%    a         b         h         w         s
% 0.070000  0.035000  0.017500  0.105000  0.100000;
cup_para = ...
[0.074600  0.047600  0.031600  0.155000  0.315200;
 0.076400  0.052100  0.033000  0.155000  0.315200;
 0.078200  0.056600  0.034500  0.155000  0.315200;
 0.080000  0.061100  0.036000  0.155000  0.315200;
 0.082000  0.066100  0.038000  0.155000  0.315200;
 0.084200  0.071100  0.040000  0.155000  0.315200;
 0.086400  0.078100  0.042000  0.155000  0.315200;
 0.088600  0.085100  0.044500  0.155000  0.315200];

% select those that are needed 
a = cup_para(cup_index, 1);
b = cup_para(cup_index, 2);
h = cup_para(cup_index, 3);
w = cup_para(cup_index, 4);
s = cup_para(cup_index, 5);

% compute point where the ellipse and cubic polynomial join
e2s          = sqrt(abs((h^2 - b^2) / (b^2/a^2)));
% set up the ellipse and its derivatives
ellipse      = @(x) sqrt(b^2 - b^2/a^2 * x.^2) .* (abs(x) <= e2s);
d_ellipse    = @(x) - ((b/a)^2 * x) ./ sqrt(b^2 - b^2/a^2 * x.^2) .* (abs(x) <= e2s);
dd_ellipse   = @(x) (((b/a)^4 * x.^2) ./ (b^2 - b^2/a^2 * x.^2)^(3/2) - ...
                     ((b/a)^2) ./ sqrt(b^2 - b^2/a^2 * x.^2)) .* (abs(x) <= e2s);
                 
% set up linear system to find the coefficients of the third cubic polynomial
A_cubic     = [e2s^3   e2s^2 e2s 1; % the polynomial is continous at e2s
               3*e2s^2 2*e2s 1   0; % the polynomial's first derivative is continous at e2s
               w^3     w^2   w   1; % the polynomial is 0 at w
               3*w^2   2*w   1   0];% the polynomial's first derivative is -s at w
b_cubic     = [h; d_ellipse(e2s); 0; -s];

% solve to get coefficients
cubic_coef  = A_cubic \ b_cubic;
cubic       = @(x) (cubic_coef(1) * x.^3 + cubic_coef(2) * x.^2 + ...
                    cubic_coef(3) * x + cubic_coef(4)) .* (x > e2s & x < w);
d_cubic     = @(x) (3*cubic_coef(1) * x.^2 + 2*cubic_coef(2) * x + ...
                      cubic_coef(3)) .* (x > e2s & x < w);
dd_cubic    = @(x) (6*cubic_coef(1) * x + 2*cubic_coef(2)) .* (x > e2s & x < w);

% set up function and derivatives
f    = @(x) ellipse(abs(x))    + cubic(abs(x));
d_f  = @(x) sign(x) .* (d_ellipse(abs(x))  + d_cubic(abs(x)));
dd_f = @(x) dd_ellipse(abs(x)) + dd_cubic(abs(x));

end