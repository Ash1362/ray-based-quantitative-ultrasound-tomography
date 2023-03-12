function a = intManipulation(a, mode)
%INTMANIPULATION performs different manipulations on integers like
%replacing them with the smallest even number larger then them
%
% DESCRIPTION:
%       intManipulation can be used to modify an input integer to a
%       close-by integer fullfilling a certain condition
%
% USAGE:
%       a = intManipulation(a,mode)
%
% INPUTS:
%       a    - integer
%       mode - mode, choose one of the following manipulations:
%           'evenUp'   - return smallest even number larger than a 
%           'oddUp'    - return smallest odd number larger than a 
%           'evenDown' - return largest even number smaller than a 
%           'oddDown'  - return largest odd number smaller than a 
%
% OUTPUTS:
%       a - modified integer
%
% ABOUT:
%       author          - Felix Lucka
%       date            - 3rd March 2017
%       last update     - 3rd March 2017

switch mode
    case 'evenUp'
        a = a + mod(a, 2);
    case 'oddUp'
        a = a + mod(a + 1, 2);
    case 'evenDown'
        a = a - mod(a, 2);
    case 'oddDown'
        a = a - mod(a + 1, 2);
    otherwise
        error('invalid mode')
end