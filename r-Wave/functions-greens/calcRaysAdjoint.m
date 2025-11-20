function [ray_position_adjoint, ray_time_adjoint, ray_absorption_adjoint,...
    ray_position_left_adjoint, ray_position_right_adjoint] = calcRaysAdjoint(...
    ray_position, ray_time, ray_absorption, ray_position_left, ray_position_right,...
    ray_spacing, rayspacing_receiver)
%CALCRAYSADJOINT computes the rays on which the adjoint pressure fields
%are computed
%
% DESCRIPTION:
%       calcRaysAdjoint computes the rays on which the adjoint pressure fields are
%       computed. The rays for the adjoint fields are not traced, but are 
%       computed by reversing the accumulated parameters along the rays used
%       for computing the forward field.

%
% USAGE:
%      
%
% INPUTS:
%       ray_position     - a 1 x num_emitter cell containing (dim * num_receiver) x
%                          num_maxraypoints Cartesian position of the points along
%                          the forward main (linked) rays. Here, dim is the dimension 
%                          of the medium, num_receiver is the number of receivers, 
%                          num_maxraypoints is a scalar representing the maximum permissible
%                          number of ray points for linking the emitter-receiver pairs.
%       ray_time         - a 1 x num_emitter cell containing num_receiver x num_raypoints 
%                          accumulated time delays on the points along the forward rays.
%       ray_absorption   - a 1 x num_emitter cell containing num_receiver x num_raypoints 
%                          accumulated acoustic absorption on the points along
%                          the forward rays.
%       ray_position_left - a 1 x num_emitter cell containing (dim * num_receiver) x
%                           num_maxraypoints Cartesian position of the points along
%                           the left auxiliary ray for comuting the adjoint field. 
%                           Here, dim is the dimension of the medium, num_receiver is the number of receivers, 
%                           num_maxraypoints is a scalar representing the maximum permissible
%                           number of ray points for linking the emitter-receiver pairs. 
%                           Note that he auxiliary rays for computing the adjoint rays
%                           are traced (not reversed), because the geomterical attenuation
%                           along the rays is not reversible.                      
%       ray_position_right - a 1 x num_emitter cell containing (dim * num_receiver) x
%                           num_maxraypoints Cartesian position of the points along
%                           the right auxiliary ray. Here, dim is the dimension 
%                           of the medium, num_receiver is the number of receivers, 
%                           num_maxraypoints is a scalar representing the maximum permissible
%                           number of ray points for linking the emitter-receiver pairs.
%                           Note that he auxiliary rays for computing the adjoint rays
%                           are traced (not reversed), because the geomterical attenuation
%                           along the rays is not reversible.
%       ray_spacing        - a scalar representing the spacing [m] along the rays
%       rayspacing_receiver - a vector containing the last sapcing along the
%                             linked rays, ie. the spacing between the
%                             reception point and the point just before that
%
% Optional INPUTS:
%              
% 
% OUTPUTS:
%       ray_position_adjoint - a 1 x num_receiver cell array containing (dim * num_emitter) x
%                          num_maxraypoint Cartesian position of the points along
%                          the forward main (linked) rays. Here, dim is the dimension 
%                          of the medium, num_receiver is the number of receivers, 
%                          num_maxraypoint is a scalar representing the maximum permissible
%                          number of ray points for linking the emitter-receiver pairs.
%       ray_time_adjoint   - a 1 x num_receiver cell array containing num_emitter x num_raypoints 
%                          accumulated time delays on the points along the forward rays.
%       ray_absorption_adjoint - a 1 x num_receiver cell containing num_emitter x num_raypoints 
%                          accumulated acoustic absorption on the points along the forward rays.
%       ray_position_left_adjoint - a 1 x num_receiver cell containing (dim * num_emitter) x
%                           num_maxraypoint Cartesian position of the points along
%                           the left auxiliary ray for comuting the adjoint field. 
%                           Here, dim is the dimension of the medium, num_receiver is the number of receivers, 
%                           num_maxraypoint is a scalar representing the maximum permissible
%                           number of ray points for linking the emitter-receiver pairs. 
%                           Note that the auxiliary rays for computing the adjoint rays
%                           are traced (not reversed), because the geomterical attenuation
%                           along the rays is not reversible.                      
%       ray_position_right_adjoint - a 1 x num_receiver cell containing (dim * num_emitter) x
%                           num_maxraypoint Cartesian position of the points along
%                           the right auxiliary ray. Here, dim is the dimension 
%                           of the medium, num_receiver is the number of receivers, 
%                           num_maxraypoint is a scalar representing the maximum permissible
%                           number of ray points for linking the emitter-receiver pairs.
%                           Note that the auxiliary rays for computing the adjoint rays
%                           are traced (not reversed), because the geomterical attenuation
%                           along the rays is not reversible.
%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 23.08.2020
%       last update     - 14.04.2023
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian 



% determine the method for computing the rays' Jacobian, which is used
% for computing the geometrical attenuation
if isempty(ray_position_left{1})

    % compute the ray's Jacobian with respect to the initial angle
    % using only the linked (main) rays and without the need for auxiliary
    % rays using a centred finite difference scheme
    attenuation_geom_method = 'raylinked';

    % set the auxiliary method nan
    auxiliary_method = nan;

     % the rays' direction perturbations are not required.
     do_get_direction_perturbation = false;

else

    % compute the ray's Jacobian with repect to the initial angle
    % using the auxiliary rays (recommended)
    attenuation_geom_method = 'auxiliary';

    if isempty(ray_position_right{1})

        % set the auxiliary method 'paraxial'
        auxiliary_method = 'paraxial';

        % the rays' direction perturbations are not required.
        do_get_direction_perturbation = false;

    else

        if (size(ray_position_right{1}, 1) == size(ray_time{1}, 1))

            % set the auxiliary method 'paraxial'
            auxiliary_method = 'paraxial';

            % using the paraxial ray, the input variable ray_position_right is used
            % for storing the perturbed angles of rays' direction because of a
            % perturbation to the initial position. This matrix should be
            % num_receiver x num_raypoints, the same as size of the input variable
            % ray_time
            do_get_direction_perturbation = true;

        else

            if (size(ray_position_right{1}, 1) == size(ray_position_left{1}, 1))

            % set the auxiliary method 'angle perturbation'
            auxiliary_method = 'angle_perturbation';

            % the rays' direction perturbations are not required.
            do_get_direction_perturbation = false;

            else

                % give an error becuse the size of the input variable
                % ray_position_right is not correct.
                error('The size of ray_position_right is not correct.')
            end

        end

    end

end



% get the number of emitters
num_emitter = length(ray_position);

% get the number of receivers
num_receiver = size(ray_time{1}, 1);

% get the maximum number of points along the rays
num_maxraypoint = size(ray_position{1}, 2);

% get the number of dimensions
dim = rem(size(ray_position{1}, 1), 2) + 2;

if dim == 3 && strcmp(auxiliary_method, 'angle_perturbation')
    error(['For the 3D case, the auxiliary rays cannot be traced with perturbing the'...
        'initial angles, but must be traced using the paraxial approach.'])
end

if dim == 3 && strcmp(attenuation_geom_method, 'raylinked')
    error(['For the 3D case, the geomterical attenuation cannot be computed'...
        'using only linked rays, and auxiliary rays must be traced using the paraxial approach.'])
end
    
% allocate cell arrays for storing the Cartesian position and time delays 
% along the reversed rays
% x position
ray_position_x_adjoint = cell(1, num_emitter);
% y position
ray_position_y_adjoint = cell(1, num_emitter);

if dim == 3
    
    % z position
    ray_position_z_adjoint = cell(1, num_emitter);
end

% accumulated time delays
ray_time_adjoint = cell(1, num_emitter);

% accumulated acoustic absorption
ray_absorption_adjoint = cell(1, num_emitter);

% determine acoustic absorption (and dispersion) are included or not 
do_absorption = ~isscalar(ray_absorption{1});

if ~do_absorption
ray_absorption_adjoint = cell(1, num_receiver);
for ind_receiver = 1: num_receiver
    ray_absorption_adjoint{ind_receiver} = 0;
end
end


% allocate cell arrays for storing the Cartesian position along the reversed auxiliary rays
if strcmp(attenuation_geom_method, 'auxiliary')
    
    % x position
    ray_position_x_left_adjoint = cell(1, num_emitter);
    
    % y position
    ray_position_y_left_adjoint = cell(1, num_emitter);
    
    if dim == 3
        
        % z position
        ray_position_z_left_adjoint = cell(1, num_emitter);
    end
    
    

if strcmp(auxiliary_method, 'angle_perturbation')
    
    % if the method used for tracing auxiliary rays is 'paraxial', the right
    % auxiliary rays are not given, but if the used method for tracing
    % auxiliary rays is 'angle_perturbation', the right auxiliary ray must
    % also be given,
    ray_position_x_right_adjoint = cell(1, num_emitter);
    ray_position_y_right_adjoint = cell(1, num_emitter);
    
       
end
    
end




%parfor (ind_emitter  = 1 : num_emitter , para.num_worker)
for ind_emitter = 1:num_emitter   

% get the nans (check to make sure the command line is true!)
nan_binary = isnan(ray_time{ind_emitter});

% get the index of the reception point on each main (linked) ray

% get ascending order for columns
[~, column_ascending_orders] = ndgrid(1:num_receiver, 1:num_maxraypoint);

% set the order for nans zeros
column_ascending_orders = ~nan_binary .* column_ascending_orders;

% get the indices of the columns associated with the reception points 
[~, receiver_indices] = max(column_ascending_orders, [], 2);

% The samplied points on the rays for computing the adjoint fields must be
% equidistant and are initilised from the reception point. However, the
% spacing between the reception point and the point just before that on each
% forward linked ray is different from the scalar ray spacing. Therefore, 
% the positions and the correponding accumulated parameters on the adjoint rays are
% corrected by 1D interpolation along the rays after reversion.

% get the interpolation cofficient for 1D interpolation for all linked rays.
% This is the distance between the reception point and the point just befor that
% over the ray spacing. 
interp_coeff = rayspacing_receiver{ind_emitter}./ray_spacing;

% reverse the x position 
ray_position_x_adjoint{ind_emitter} = reverseRayParameter(ray_position{ind_emitter}(1:dim:end, :),...
    interp_coeff, receiver_indices, nan_binary, false);

% reverse the y position 
ray_position_y_adjoint{ind_emitter} = reverseRayParameter(ray_position{ind_emitter}(2:dim:end, :),...
    interp_coeff, receiver_indices, nan_binary, false);

if dim == 3
    
    % reverse the z position
    ray_position_z_adjoint{ind_emitter} = reverseRayParameter(ray_position{ind_emitter}(3:dim:end, :),...
        interp_coeff, receiver_indices, nan_binary, false);
end

% reverse the accumulated time delays
ray_time_adjoint{ind_emitter} = reverseRayParameter(ray_time{ind_emitter},...
    interp_coeff, receiver_indices, nan_binary, true);

% reverse the accumulated acoustic absorption
if do_absorption
    ray_absorption_adjoint{ind_emitter} = reverseRayParameter(ray_absorption{ind_emitter},...
        interp_coeff, receiver_indices, nan_binary, true);
end

end

% convert the cells for the parameters the rays to a three dimensional matrix
ray_position_x_adjoint = cell2mat(permute(ray_position_x_adjoint, [1, 3, 2]));
ray_position_y_adjoint = cell2mat(permute(ray_position_y_adjoint, [1, 3, 2]));

if dim == 3
    ray_position_z_adjoint = cell2mat(permute(ray_position_z_adjoint, [1, 3, 2]));
end

ray_time_adjoint = cell2mat(permute(ray_time_adjoint, [1, 3, 2]));

ray_position_adjoint = zeros(dim*num_emitter, num_maxraypoint, num_receiver);

% change the order of the dimensions of matrices to num_emitter x num_raypoints x
% num_receiver, because for calculation of the adjoint fields, the emitters
% and receivers are interchenged
ray_position_adjoint(1:dim:end, :,:) = permute(ray_position_x_adjoint, [3, 2, 1]);
ray_position_adjoint(2:dim:end, :,:) = permute(ray_position_y_adjoint, [3, 2, 1]);
if dim == 3
    ray_position_adjoint(3:dim:end, :,:) = permute(ray_position_z_adjoint, [3, 2, 1]);
end
ray_time_adjoint = permute(ray_time_adjoint, [3, 2, 1]);

% convert the three dimensional matrices for the position of points on the main rays to cell arrays
% with num_receiver cells corresponding to the implicit excitations for the adjoint field
ray_position_adjoint = squeeze(mat2cell(ray_position_adjoint, dim * num_emitter,...
    num_maxraypoint, ones(num_receiver, 1)) );  
ray_time_adjoint = squeeze(mat2cell(ray_time_adjoint, num_emitter,...
    num_maxraypoint, ones(num_receiver, 1)) );

if do_absorption
    ray_absorption_adjoint = cell2mat(permute(ray_absorption_adjoint, [1, 3, 2]));
    ray_absorption_adjoint = permute(ray_absorption_adjoint, [3, 2, 1]);
    ray_absorption_adjoint = squeeze(mat2cell(ray_absorption_adjoint, num_emitter,...
        num_maxraypoint, ones(num_receiver, 1)) );
end

%%=========================================================================
% REVERSE THE AUXILIARY RAYS
%==========================================================================
switch attenuation_geom_method
    
    case 'auxiliary'
        
%parfor (ind_emitter  = 1 : num_emitter , para.num_worker)         
for ind_emitter = 1 : num_emitter
    
% rearrange the left auxiliary ray for the adjoint field

% get the x position of the points on the left auxiliary rays
ray_position_x_left_adjoint{ind_emitter} = ray_position_left{ind_emitter}(1:dim:end, :);
% get the y position of the points on the left auxiliary rays
ray_position_y_left_adjoint{ind_emitter} = ray_position_left{ind_emitter}(2:dim:end, :);

if dim == 3
    
    % get the z position of the points on the left auxiliary rays
    ray_position_z_left_adjoint{ind_emitter} = ray_position_left{ind_emitter}(3:dim:end, :);
end


if strcmp(auxiliary_method, 'angle_perturbation')

% get the position of the points on the right auxiliary rays in the Cartesian
% coordinates
ray_position_x_right_adjoint{ind_emitter} = ray_position_right{ind_emitter}(1:dim:end, :);
ray_position_y_right_adjoint{ind_emitter} = ray_position_right{ind_emitter}(2:dim:end, :);
end


end
      
% convert the cells for the parameters of the adjoint rays to a three dimensional matrix
% x cartesian position
ray_position_x_left_adjoint = cell2mat(permute(ray_position_x_left_adjoint, [1, 3, 2]));

% y cartesian position
ray_position_y_left_adjoint = cell2mat(permute(ray_position_y_left_adjoint, [1, 3, 2]));

if dim == 3
    
   % z cartesian position
   ray_position_z_left_adjoint = cell2mat(permute(ray_position_z_left_adjoint, [1, 3, 2])); 
end

% allocate zero matrices for the position on the rays
ray_position_left_adjoint = zeros(dim * num_emitter, num_maxraypoint, num_receiver);


% change the order of the dimensions of matrices to num_emitter x num_raypoints x
% num_receiver, because for computing the adjoint fields, the emitters and receivers 
% are swapped
ray_position_left_adjoint(1:dim:end, :,:) = squeeze(permute(ray_position_x_left_adjoint, [3, 2, 1]));
ray_position_left_adjoint(2:dim:end, :,:) = squeeze(permute(ray_position_y_left_adjoint, [3, 2, 1]));
if dim == 3
    ray_position_left_adjoint(3:dim:end, :,:) = squeeze(permute(ray_position_z_left_adjoint, [3, 2, 1]));
end
    
% convert the three dimensional matrices for the position of points on the left auxiliary rays to cell arrays
% with num_receiver cells corresponding to the implicit excitations for the adjoint field
ray_position_left_adjoint = squeeze(mat2cell(ray_position_left_adjoint, dim * num_emitter,...
    num_maxraypoint, ones(num_receiver, 1)) ); 


switch auxiliary_method
    
    case 'angle_perturbation'

        
        % repeat the procedure for the right auxiliary rays the same as for the left auxiliary
        % rays

        % make the cell array a matrix of size num_receiver x num_raypoints
        % x num_emitter for each cartesian coordinate
        ray_position_x_right_adjoint = cell2mat(permute(ray_position_x_right_adjoint, [1, 3, 2]));
        ray_position_y_right_adjoint = cell2mat(permute(ray_position_y_right_adjoint, [1, 3, 2]));
        
        ray_position_right_adjoint = zeros(dim * num_emitter, num_maxraypoint, num_receiver);
        
        ray_position_right_adjoint(1:dim:end, :,:) = squeeze(permute(ray_position_x_right_adjoint, [3, 2, 1]));
        ray_position_right_adjoint(2:dim:end, :,:) = squeeze(permute(ray_position_y_right_adjoint, [3, 2, 1]));
        
        
        % convert the three dimensional matrices for the position of points on the right auxiliary rays to cell arrays
        % with num_receiver cells corresponding to the implicit excitations for the adjoint field
        ray_position_right_adjoint = squeeze(mat2cell(ray_position_right_adjoint, dim *num_emitter,...
            num_maxraypoint, ones(num_receiver, 1)) );

    case 'paraxial'

        if do_get_direction_perturbation

        % make the cell array a matrix of size num_receiver x num_raypoints
        % x num_emitter
        ray_position_right_adjoint = cell2mat(permute(ray_position_right, [1, 3, 2]));

        % To test
        ray_position_right_adjoint = permute(ray_position_right_adjoint, [3, 2, 1]);


        % convert the three-dimensional matrices for the perturbed direction vector of the rays in
        % in the polar coordinate to a cell array with num_receiver cells corresponding to the
        % implicit excitations for computing the adjoint field. Each cell
        % contains  a matrix of size num_emitter x num_raypoints containg
        % the perturbed ray's direction vectors in the polar coordinate
        % because of perturbation to the initial position of the backward
        % (adjoint) ray. The initial positions of the backward rays are
        % receivers
        ray_position_right_adjoint = squeeze(mat2cell(ray_position_right_adjoint, num_emitter,...
            num_maxraypoint, ones(num_receiver, 1)) );


        else
       
        % get a cell array of empty variables
        ray_position_right_adjoint = cell(1, num_receiver);

        end

end

    case 'raylinked'
      
ray_position_left_adjoint = cell(1, num_receiver);
ray_position_right_adjoint = cell(1, num_receiver);
 
end

end