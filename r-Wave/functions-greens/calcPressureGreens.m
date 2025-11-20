function [pressure_grid, pressure_receiver] = calcPressureGreens(...
    approximate_pressure, parameters_grid, caustic_number, nan_grid_binary,...
    pressure_source, source_mode, parameters_receiver, caustic_receiver,...
    receiver_order, num_frequency, num_receiver)
%CALCPRESSUREGREENS computes the pressure field on the grid points and
% receivers using the Green's function
%
%   % DESCRIPTION:
%       calcPressureGreens computes the pressure field on the grid points
%       and receivers using the Green's function for the heterogeneous and
%       absorbing media. The parameters of the Green's function on the grid
%       points and receivers are given.
%
%
% USAGE:
%
%
% INPUTS:
%       approximate_pressure - the handle function for approximating the
%                             pressure field given a pressure source
%                             and the parameters of the Green's function
%       parameters_grid     - a matrix with rows the grid points inside
%                             the binary mask, and columns: time delays,
%                             geoemetrical attenuation including the refraction
%                             effects, acoustic absorption (if not zero),
%                             and relative sound speed (the sound speed over
%                             the sound speed in water)
%       nan_grid_binary     - a binary mask indicating to the grid
%                             points outside the given binary mask
%       pressure_source     - pressure source
%       source_mode         - the source mode for computing the pressure
%                             field on the grids. This can be 'normal' or
%                             'time-reversal'. The method for computing
%                             pressure on the receivers is always 'normal'.
%       caustic_number      - the cumulative integer number of the caustics
%                             along the rays 
%       parameters_receiver-  a matrix with rows the receivers, and columns:
%                             time delays, geometrical attenuation including
%                             refraction effects, acoustic absorption (if not zero),
%                             and relative sound speed. If empty, the mode will be
%                             automatically set 'adjoint'
%       caustic_receiver    - the cumulative times the sign of the rays'
%                             Jacobian has been changed before the ray is
%                             intercepted by the receiver (reception point).
%       receiver_order      - a vector indicating the order of the receivers in
%                             rows of the matrix 'parameters_receiver'.
%       num_frequency       - the number of frequencies
%
%
%
% OUTPUTS:
%       pressure_grid     - the pressure field on the grid points
%       pressure_receiver - the pressure field on the receivers (if the
%                            mode is 'forward')

%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.03.2020
%       last update     - 30.03.2020
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian

if isempty(parameters_receiver)
    
    % If the parameters on the receivers are not given
    if ~isempty(caustic_receiver) || ~isempty(receiver_order)
        error('The inputs for the receivers must be all empty, or nonempty.')
    end
    
    % If the inputs associated with the receivers are not given,
    % set the mode 'adjoint'
    field_mode = 'adjoint';
    
else
    
    % If the inputs associated with the receivers are given,
    % set the mode 'adjoint'
    field_mode = 'forward';
    
end

if size(parameters_grid, 2) < 4
    
    % if the number of columns is 3 or less, ignore the acoustic absorption
    % and dispersion
    do_absorption = false;
    ray_absorption = 0;
    
else
    
    % include the acoustic absoprion and dispersion
    do_absorption = true;
    
end

if ~isempty(parameters_receiver)
    
    % the number of columns in the 'parameters_receiver' must be always one
    % fewer than 'parameters_grid', because the columns for 'sound_speed_relative'
    % is always 1, and is not stored.
    if size(parameters_receiver, 2) ~= size(parameters_grid, 2)-1
        error(['The size of the two matrices for prameters on the grid and receivers'...
            'are not consistent.']);
    end
end


%if strcmp(field_mode, 'forward')
    
    % get the number of receivers
    %num_receiver = 128;
    % num_receiver = size(parameters_receiver, 1);
%end

% get the number of grid points
num_gridpoints = length(nan_grid_binary);

% allocate a matrix for the pressure on the grid points
pressure_grid = zeros(num_gridpoints, num_frequency);

% reverse the values of the last column, if the approach for computing the
% Green's function is time-reversal, because this coefficient is used for 
% computing the amplitude, which must be reversed.
if strcmp(source_mode, 'time-reversal')
     parameters_grid(:, end) = 1./parameters_grid(:, end);
end
    
% compute the pressure on the grid points given the parameters of the
% Green's function on the grid
if do_absorption
    
    % include acoustic absorption and dispersion
    pressure_grid(nan_grid_binary, :) = sqrt(parameters_grid(:, 4))  .*... 
        approximate_pressure(pressure_source, parameters_grid(:, 1), parameters_grid(:, 2),...
        parameters_grid(:, 3), caustic_number, source_mode);
    
else
    
    % ignore acoustic absorption and dispersion
    pressure_grid(nan_grid_binary, :) = sqrt(parameters_grid(:, 3))  .*... 
        approximate_pressure(pressure_source, parameters_grid(:, 1), parameters_grid(:, 2),...
        ray_absorption, caustic_number, source_mode);
    
end



switch field_mode
    case 'forward'
        
        % approximate the unsorted pressure on the receivers using the Green's
        % function approaches
        if do_absorption
            
            % include the acoustic absorption and dispersion
            pressure_receivers_unsorted = approximate_pressure(pressure_source, parameters_receiver(:, 1),...
                parameters_receiver(:, 2), parameters_receiver(:, 3), caustic_receiver, 'normal');
        else
            
            % ignore the acoustic absorption and dispersion
            pressure_receivers_unsorted = approximate_pressure(pressure_source, parameters_receiver(:, 1),...
                parameters_receiver(:, 2), 0, caustic_receiver, 'normal');
        end
        
        % allocate a matrix for the pressure field on the receivers
        pressure_receiver = zeros(num_receiver, size(pressure_grid, 2));
        
        % get the pressure time series on the receivers
        pressure_receiver(receiver_order ,:) = ...
            pressure_receivers_unsorted;
          
    case 'adjoint'
        
        % allocate an empty variable for the pressure field on the receivers
        pressure_receiver = [];
               
end



end