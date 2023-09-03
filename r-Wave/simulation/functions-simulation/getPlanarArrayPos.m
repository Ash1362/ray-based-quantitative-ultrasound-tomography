function [transducer_position] = getPlanarArrayPos(...
    array_transducer_number, array_size, array_center_polar, rot_radius, mode)
% GETPLANARARRAYPOS computes the position of the transduces along a rotating line
%(or planar) array
%
% DESCRIPTION:
%           getPlanarArrayPos computes the Cartesian position of the transduces
%           along a rotating line (or planar) array along the x-y planes and around
%           the z axis
%          
%          
%
% USAGE:
%        
%
% INPUTS:
%       array_transducer_number  - 1 x dim-1 number of transducers along the x-z coordinates
%       array_size               - 1 x dim-1 size [m] of the array along the x-z coordinates
%       array_center_polar       - 1 x num_pos polar positions of the centre of arrays
%       rot_radius               - the radius [m] of a circle along which the centre of array rotates
%       mode                     - the mode, which can be 'emitter ' or 'receiver'

% OPTIONAL INPUTS:
%
%
% OUTPUTS:
%       cart_pos_transducer     - if mode is set 'receiver' a 1 x num_pos cell array,
%                                 each containing the Cartesian position of the transducers
%                                 for each angular position of the centre
%                                 of array. if mode is set 'emitter', the a
%                                 dim x num_emitter Cartesian position of emitters 

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 31.01.2023
%       last update     - 31.01.2023
%
% 

% get the number of dimensions
dim = length(array_transducer_number) + 1;

% get the number of position of the centre of the array along a circular
% rotation axis
num_pos = length(array_center_polar);

% get the Cartesian position of the centre of arrays
[array_center_cart_x, array_center_cart_y] = pol2cart(array_center_polar, rot_radius);


switch dim
    case 2
        
% get the cartesian position of the centre of the starting and ending
% transduces for each position of the array
[transducer_pos_start, transducer_pos_end, ~] = lineEdgesCircularDetection([array_center_cart_x;...
    array_center_cart_y], 1/2 * array_size);

% allocate a cell array for the Cartesian position of the transducers
transducer_position = cell(1, num_pos);

for ind_pos = 1 : num_pos
   
% get the Cartesian position of the transducers for each position of the centre of the array
transducer_cartesian_position = zeros(2, array_transducer_number);

% get the x positions 
transducer_cartesian_position(1, :) = linspace(transducer_pos_start(1, ind_pos),...
    transducer_pos_end(1, ind_pos), array_transducer_number);

% get the y positions
transducer_cartesian_position(2, :) = linspace(transducer_pos_start(2, ind_pos),...
    transducer_pos_end(2, ind_pos), array_transducer_number);

transducer_position{ind_pos} = transducer_cartesian_position;

end

% make the positions a matrix if the transducers are emitter
if strcmp(mode, 'emitter')
    transducer_position = cell2mat(transducer_position);
end
    

end