function x = vectorise(x)
%VECTORISE vetorises a matrix
%
% DESCRIPTION:
%       vectorise(x) gives x(:)
%
% USAGE:
%       
%
% INPUTS:
%       x - a matrix (multidimensional array)
%                 
%
% OUTPUTS:
%       x  - a stacked vector of the entries of x 
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.19.2020
%       last update     - 30.19.2020
%
% See also reshape

x = x(:);
