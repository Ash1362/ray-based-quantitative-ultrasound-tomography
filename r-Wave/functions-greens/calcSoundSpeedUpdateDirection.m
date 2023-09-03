function [sound_speed_update_direction, objective_function, computation_time,...
    forward_signal, jacobian_signal] = calcSoundSpeedUpdateDirection(...
    pressure_data, pressure_source_time, time_array, omega, parameters_grid, parameters_grid_adjoint,...
    parameters_receiver, nan_grid_binary, nan_grid_binary_adjoint, caustic_number,...
    caustic_number_adjoint, receiver_order, caustic_receiver, num_rays, perturbation, directions_grid,...
    directions_grid_adjoint, sound_speed, absorption, mask, random_residual,...
    omega_spacing, step_length, varargin)
%calcSoundSpeedUpdateDirection computes the update direction of the sound speed
% using the Green's approach
%
% DESCRIPTION:
%       calcSoundSpeedUpdateDirection computes the update direction of the
%       sound speed. The update can be the gradient of the residual norm,
%       a perturbation update reprsenting th action of the inverse of the Hessian
%       matrix on the gradient. The third approach is the backprojection,
%       in which the Hessian matrix is diagonalized, and is directly
%       inversed.
%
%
% USAGE:
%
%
% INPUTS:
%       pressure_data   - a num_receiver x num_time x num_emitter matrix
%                         containing the object-in-water pressure data
%       pressure_source_time - the pressure source in time
%       time_array      - the time array [s]
%       omega           - the vector of angular frequencies [rad/s]
%       parameters_grid - a num_emitter x 1 cell array for the approximated
%                         parameters of the Green's function on the grid
%                         points after being produced by emitters. The
%                         rows of the matrix corresponding to each emitter
%                         represent the grid points inside a binary mask,
%                         and columns represent time delays, geometrical attenutaion,
%                         acoustic absorption (if not zero), and relative sound speed
%                         for the forward pressure field produced by the emitter
%       parameters_grid_adjoint - a num_receiver x 1 cell array for the approximated
%                         parameters of the Green's function on the grid
%                         points after being produced by receivers. The
%                         rows of the matrix corresponding to each receiver
%                         represent the grid points inside a binary mask,
%                         and columns represent time delays, geometrical attenutaion,
%                         acoustic absorption (if not zero), and relative sound speed
%                         for the adjoint pressure field produced by the
%                         receiver
%       parameters_receivers - a num_emitter x 1 cell array each containing a matrix
%                              with rows the receivers, and columns: time delays,
%                              geoemetrical attenutaion, acoustic absorption (if not zero),
%                              and relative sound speed for the forward pressure field
%                              produced by the emitter on all the receivers
%       nan_grid_binary      - a num_emitter x 1 cell array each containing
%                              the grid points outside a binary mask, which are
%                              given nan values for the parameters for the forward
%                              pressure field
%       nan_grid_binary_adjoint - a num_receiver x 1 cell array each containing
%                                 the grid points outside a binary mask, which are
%                                 given nan values for the parameters for the adjoint
%                                 pressure field
%       caustic_number       - a num_emitter x 1 cell array each containing
%                              the cumulative integer number of the caustics
%                              on the rays initialised from the emitter. The numbers are
%                              interpolated on the grid points
%       caustic_number_adjoint  - a num_receiver x 1 cell array each containing
%                                 the cumulative integer number of the caustics
%                                 from the rays emanated from the receiver on
%                                 each grid point inside the mask
%       receiver_order        - a vector indicating the order of receivers in
%                               the variable 'parameters_receivers'.
%       caustic_receiver     - a num_emitter x 1 cell array each containing
%                              the cumulative integer number of caustics
%                              along the rays initialised from emitters and
%                              linking to all the receivers
%       num_rays             - a num_emitter x 1 cell array containing num_receiver x 1
%                             number of rays used for linking each
%                             emitter-receiver pair
%       directions_grid       - a num_emitter x 1 cell array each containing
%                              the direction of the forward rays
%       directions_grid_adjoint - a num_receiver x 1 cell array each
%                              containing the direction of the adjoint rays
%       perturbation       - a perturbation in the sound speed on the grid
%                           points
%       sound_speed        - the sound speed distribution on the grid points
%       absorption         - a struct containing:
%       'coeff'             the absorption coefficient [dB MHz^{-y} cm^{-1}] on the grid points
%       'power'           - the scalar representing the exponent power
%       sound_speed_background - the background (only-water) sound speed
%       mask               - the binary mask, which is true for the grid
%                          points included in the image reconstruction
%       random_residual    - a random residual on which the action of the
%                           adjoint of the Frechet derivative is tested
%                           through an inner product test. For image reconstruction.
%                            this variable is set an empty variable.
%       omega_spacing     - the spacing [rad s^{-1}] of the angular frequency
%       step_length       - the step length for computing the update direction
%
%
% Optional INPUTS:
%       num_workers                - the number of workers for parallel programming
%       task                       - the task done by the m-file function.
%                                    This can be 'forward', 'gradient', 'hessian', or
%                                    'backprojection'
%       max_rays                   - the maximum number of rays used for solving the
%                                    ray linking inverse problems
%       greens_domain              - the domain for approximating the
%                                    pressure field using Green's function,
%                                    this can be 'time', or 'frequency'
%       greens_method              - the method for calculating the pressure
%                                    field, This can be 'analytic',...
%       record_output_signals      - two binaries the first controlling whether
%                                    the forward pressure field recorded on the receivers
%                                    is stored and is given as an output of the function,
%                                    and the second controlling whether the action of the
%                                    Frechet derivative on the sound speed perturbation
%                                    $ J \delta c$ s recorded on the receivers is given as
%                                    an output of the function. The latter will be used for
%                                    checking the accuracy of the Frechet derivative and its adjoint
%                                    through an inner product test.
%      deconvolve_source           - Bolean controlling whether the source is deconvolved or not
%      deconvolution_parameter     - the regularisation parameter for
%                                    Tikhonov-based deconvolution
%
% OUTPUTS:
%           sound_speed_upadte_direction - the sound speed update direction
%                                           on the grid points
%           objective_function   - the half squared L2 norm of the residual
%           computation_time     - the time [s] for running the m-file function
%           forward_signal       - the forward pressure field recorded on the
%                                  receivers in the frequency domain. Empty variable,
%                                  if the first component of the optional input
%                                  'record_output_signals' is set false.
%           jacobian_signal      - the action of the Frechet derivative (Jacobian
%                                  matrix on the sound speed perturbation.
%                                  This is recorded on the receivers in the frequency domain.
%                                  Empty variable, if the second component of the optional input
%                                  'record_output_signals' is set false.
%
%
%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.04.2020
%       last update     - 21.12.2021
%
% This script is part of the r-Wave toolbox.
% Copyright (c) 2022 Ashkan Javaherian

para = [];
para.num_workers = 16;
para.task = 'gradient';
para.max_rays = 1000;
para.greens_domain = 'frequency';
para.greens_method = 'analytic';
para.record_output_signals = [false, false];
para.emitter_index = 1;
para.deconvolve_source = true;
para.deconvolution_parameter = 0;


if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


% if the vector of angular frequencies is not given, set the optional
% inputs 'greens_domain' and 'greens_method' time
if isempty(omega)
    error('The time-domain approach is not implemented yet.')
    para.greens_domain = 'time';
    para.greens_method = 'time';
end


if strcmp(para.task, {'backprojection'})
    if  isempty(directions_grid{1})  ||  isempty(directions_grid_adjoint{1})
        Error(['Using the backprojection approach,'...
            'the information for the direction of the rays must be given as an input.'])
    end
end


tstart = tic;

% get the size of the time arrays
num_time = length(time_array);

% get the number of emitters
num_emitter = length(receiver_order);

% get the number of receivers
num_receiver = length(receiver_order{1});


% Convert the simulated or measured data to a cell array, if the data is in the
% form of a three-dimesnional matrix
if any(strcmp(para.task, {'forward', 'gradient','backprojection'}))
    
    if isempty(pressure_data)
        
        % If the task is set forward or gradient, the pressure data must be
        % nonempty.
        error(['The pressure data must be nonempty, if the task is set forward,'...
            'gradient or backprojection.'])
        pressure_data = cell(num_emitter, 1);
    else
        if ~iscell(pressure_data)
            
            % check the size of the given forward pressure data
            if any(size(pressure_data) ~= [num_receiver, num_time, num_emitter])
                error('The size of the given forward pressure data is not consistent.')
            end
            pressure_data = mat2cell(pressure_data, num_receiver, size(pressure_data, 2),...
                ones(num_emitter, 1));
        end
    end
end



% some input checks
switch para.task
    
    case 'hessian'
        
        if  isempty(perturbation)
            error('If the task is set hessian, the perturbation vector must not be empty.')
        end
        
        if isempty(random_residual)
            
            
            % all emitters are used for image reconstruction
            emitter_indices = 1: num_emitter;
            
            % If the random residual is empty, the purpose of implementing the
            % function is image reconstruction, and therefore, the action of the
            % Jacobian on the perturbation does not need to be stored
            para.record_output_signals(2) = false;
            
            
        else
            
            % If the purpose is an inner product test, only a single emitter will be used
            emitter_indices = para.emitter_index;
            
            % deactivate the parallel programming by setting the number of
            % workers zero
            para.num_workers = 0;
            
            % If the random residual is not cell array, check the size and
            % and convert it to a cell array
            if ~iscell(random_residual)
                
                % check the size of the given random residual
                if any(size(random_residual) ~= [num_receiver, num_time, num_emitter])
                    disp('The size of the given random residual is not consistent.')
                end
                
                random_residual = mat2cell(random_residual, num_receiver,...
                    size(random_residual, 2), ones(num_emitter, 1));
            end
            
            % If the random residual is nonempty, the purpose of implementing the
            % function is an inner product test, and therefore, the action of the
            % Jacobian on the perturbation does not need to be stored.
            para.record_output_signals(2) = true;
            
        end
        
    case {'forward', 'gradient', 'backprojection'}
        
        if   ~all([isempty(perturbation), isempty(random_residual)])
            
            % If the task is set 'forward' or 'gradient', both the perturbation
            % vector and random residual must be empty.
            error (['If the task is forward, gradient or backprojection, both perturbation vector and'...
                'the random residual must be set as empty inputs.'])
            
        else
            
            % all emitters are used for image reconstruction
            emitter_indices = 1: num_emitter;
            
        end
end

if any(strcmp(para.task, {'gradient', 'backprojection', 'hessian'}))
    
    % Compute the adjoint field if the task is set either 'gradient', 'backprojection' or
    % 'hessian'
    calc_adjointfield = true;
else
    
    % Do not compute the adjoint field if the task is set 'forward'
    calc_adjointfield = false;
end

% Get the sound speed on the grid points inside the mask for image
% reconstruction
sound_speed = sound_speed(mask);

% Get the power law exponent
y = absorption.power;

if size(parameters_grid{1}, 2) < 4
    
    % If parameters_grid should has less than four columns, absorption map is
    % on the grid points is not given, and is assumed zero.
    absorption = zeros(nnz(mask), 1);
    
else
    
    % Confine the absorption to the grid points inside the mask for image
    % reconstruction, and convert the unit from [dB MHz^{-y} cm^{-1}] to
    % [napers (rad/s)^{-y} m^{-1}]
    absorption = db2neper(absorption.coeff(mask), y);
    
end



switch para.greens_method
    
    case  'analytic'
        
        % square angular frequencies
        omega_squared = omega.^2;
        
        % fractional power plus 1  (y+1)
        omega_frac_plus = omega.^(y+1);
        
        % fractional power minus 1 (y-1_
        omega_frac_minus = omega.^(y-1);
        
    case {'time'}
        
        Error('Not implemented yet!')
        
end

% get the excitation pulse
switch para.greens_domain
    
    case 'time'
        
    case 'frequency'
        
        switch para.greens_method
            
            
            case  'analytic'
                
                if para.deconvolve_source
                    
                    % deconvolve the pressure time series from the
                    % excitation pulse
                    transform_frequency = @(data) deconvolve(pressure_source_time,...
                        data, time_array, omega, para.deconvolution_parameter);
                    
                    % make the pressure source 1, if the data is deconvolved
                    % from the pressure source
                    pressure_source = 1;
                    
                    
                else
                    
                    % Define a handle function for computing the Discrete
                    % Fourier Transform (DFT) of a set of time series
                    transform_frequency = @(data) discreteFourierTransform(...
                        data, omega, time_array);
                    
                    % Get the DFT of the pressure source
                    pressure_source = transform_frequency([pressure_source_time,...
                        zeros(1, num_time-length(pressure_source_time))]);
                    
                    
                end
                
        end
        
end


% Define a handle function for approximating the pressure field using the
% Green's function
approx_pressure = @(source_pulse, time_delays, geom_attenuation, absorption,...
    caustic_number, source_mode) approxPressureGreens(source_pulse, time_delays,...
    geom_attenuation, absorption, caustic_number, omega, y, [], source_mode);


% Define a handle function for approximating the pressure field on the grid points and
% the receivers using the Green's function
calc_pressure = @(approx_pressure, parameters_grid, caustic_number, nan_grid_binary,...
    pressure_source, parameters_receiver, caustic_receiver, receiver_order,...
    source_mode) calcPressureGreens(approx_pressure, parameters_grid, caustic_number,...
    nan_grid_binary, pressure_source, source_mode, parameters_receiver,...
    caustic_receiver, receiver_order, length(omega));

% Initialise the update direction of the sound speed by zero
sound_speed_update = 0;

% Initialise the wave number update of the sound speed by zero
slowness_squared_update = 0;

% Initialise the objcetive function, i.e., the half square of the L2 norm of
% the residual, by zero
objective_function = 0;

if para.record_output_signals(1)
    
    % Allocate a cell array for the forward signals on the receivers, if
    % required
    forward_signal = cell(num_emitter, 1);
else
    
    % allocate an empty variable if the forward pressures on the receivers are not stored
    forward_signal = [];
    
end


if para.record_output_signals(2)
    
    % Allocate a matrix for the action of the Frechet derivative on the perturbation
    % of the sound speed, $J \delta c$. This will be used for testing the Frechet derivative
    % and its adjoint using an inner product test. However, for solving the linearised subproblem,
    % the actions $J^T J \delta c_l$  are computed by matrix-vector products without the need for
    % storage of $J \delta c$, so this variable eill not be used.
    jacobian_signal = cell(num_emitter, 1);
    
else
    
    % allocate an empty variable if the action of the Frechet derivative
    % on the sound speed perturbation is not stored
    jacobian_signal = [];
    
end

parfor (ind_emitter = emitter_indices, para.num_workers)
% for ind_emitter = 1: num_emitter
    
    %% ============================================================
    % FORWARD PRESSURE FIELD
    %==============================================================
    % display the current emitter for which the update direction of
    % the sound speed is computed.
    disp(['Computing the update direction of the sound speed for emitter:'...
        num2str(ind_emitter)]);
    
    % approximate the forward pressure field on the grid points (using the
    % requested source model) and the receivers (using the normal source mode),
    % given the rays' parameters on the grid points and the receivers.
    switch para.task
        case {'forward','gradient'}
            
            [pressure_grid, pressure_receiver] = calc_pressure(approx_pressure,...
                parameters_grid{ind_emitter}, caustic_number{ind_emitter},...
                nan_grid_binary{ind_emitter}, pressure_source,...
                parameters_receiver{ind_emitter}, caustic_receiver{ind_emitter},...
                receiver_order{ind_emitter}, 'normal');
            
        case 'backprojection'
            
            [pressure_grid, pressure_receiver] = calc_pressure(approx_pressure,...
                parameters_grid{ind_emitter}, caustic_number{ind_emitter},...
                nan_grid_binary{ind_emitter}, pressure_source,...
                parameters_receiver{ind_emitter}, caustic_receiver{ind_emitter},...
                receiver_order{ind_emitter}, 'time-reversal');
            
        case 'hessian'
            
            [pressure_grid, ~] = calc_pressure(approx_pressure,...
                parameters_grid{ind_emitter}, caustic_number{ind_emitter},...
                nan_grid_binary{ind_emitter}, pressure_source, [],[],[],...
                'normal');
            
    end
    
    
    if any(strcmp(para.task, {'forward', 'gradient', 'backprojection'}))
        
        % calculate the residual, i.e., the discrepancy of the the approximated
        % pressure using Green's function on the receivers and the measured (simulated) data on the receivers
        residual = pressure_receiver - transform_frequency(pressure_data{ind_emitter});
        
        % exclude from the residual the rays that either travel through only water,
        % or fail to reach the reception point after ray linking
        residual(num_rays{ind_emitter}== 1 | num_rays{ind_emitter} == para.max_rays, :) = 0;
        % residual(num_rays{ind_emitter} == para.max_rays, :) = 0;
        
        % calculate the L2 norm of the residual for the current emitter, and add it to
        % the objective function
        objective_function = objective_function + norm(residual, 'fro')^2;
        
    end
    
    
    if para.record_output_signals(1)
        
        % get the forward pressure field recorded on the receivers
        forward_signal{ind_emitter} = pressure_receiver;
        
    end
    
    
    if calc_adjointfield
        
        
        % Initialise the adjoint pressure field for the current emitter
        pressure_grid_adjoint = 0;
        
        % Allocate a matrix for the action of the Frechet derivative
        % (Jacobian matrix) on the sound speed perturbation for the
        % current emitter
        %  jacobian_signal_receiver = zeros(num_receiver, num_frequency);
        
        for ind_receiver = 1 : num_receiver
            
            
            switch para.task
                
                case 'hessian'
                    
                    %% ============================================================
                    % ACTION OF THE JACOBIAN ON THE PERTURBATION (J \delta c_l)
                    %==============================================================
                    % Compute the action of the Frechet derivative on
                    % the sound speed perturbation vector, i.e. the product
                    % of the Jacobian matrix on the perturbation vector.
                    
                    % Multiply the forward pressure at every grid point by the scattering
                    % potential
                    if num_rays{ind_emitter}(ind_receiver) < para.max_rays
                        
                        % compute p(x, x_e) x \Upsilon x
                        % perturbation
                        res_signal = - 2 * perturbation(nan_grid_binary_adjoint{ind_receiver}) .*...
                            pressure_grid(nan_grid_binary_adjoint{ind_receiver}, :) .* ...
                            (( 1 ./ sound_speed(nan_grid_binary_adjoint{ind_receiver}).^3 .* omega_squared +...
                            (tan(pi * y/2) + 1i)  ./ sound_speed(nan_grid_binary_adjoint{ind_receiver}).^2 .*...
                            absorption(nan_grid_binary_adjoint{ind_receiver}) .* omega_frac_plus));
                        
                        % Approximate the frequency-domain forward pressure emanated from the grid points
                        % including nans on the receiver
                        [res_signal, ~] = calc_pressure(approx_pressure,...
                            parameters_grid_adjoint{ind_receiver},...
                            caustic_number_adjoint{ind_receiver},....
                            nan_grid_binary_adjoint{ind_receiver}, res_signal,...
                            [], [], [], 'normal');
                        
                        % add a line for test
                        %  res_signal(num_rays{ind_emitter}== 1 | num_rays{ind_emitter} == para.max_rays, :) = 0;
                        
                        % Add the frequency-domain signals for the non-nan grid points
                        % for getting the action of the Frechet derivative
                        % (Jacobian matrix) on the perturbation vector for
                        % the current emitter-receiver pair
                        res_signal = sum(res_signal(nan_grid_binary_adjoint{ind_receiver}, :));
                        
                    else
                        
                        res_signal = 0;
                        
                    end
                    
                    
                    
                    if  para.record_output_signals(2)
                        
                        % Store the action of the Frechet derivative (Jacobian matrix)
                        % on the sound speed perturbation for the current receiver, if the purpose is an
                        % an inner product test on the Frechet derivative and its adjoint
                        jacobian_signal{ind_emitter}(ind_receiver, :) = res_signal;
                    end
                    
                    
                    if ~isempty(random_residual)
                        
                        % If the random residual is nonempty, and
                        % therefore, the purpose is an inner
                        % product test, transform the random residual
                        % to the frequency domain
                        res_signal = transform_frequency(random_residual{ind_emitter}(ind_receiver, :));
                        
                    end
                    
                    % compute the complex conjugate of the signal
                    if strcmp(para.task, 'hessian')
                        if isempty(random_residual)
                            
                            % If the random residual is nonempty, and
                            % therefore, the purpose is computing the sound speed update,
                            % transform the random residual
                            % to the frequency domain
                            res_signal = conj(res_signal);
                        end
                    end
                    
                case 'gradient'
                    
                    % get the complex conjugate of the residual for the current receiver
                    res_signal = conj(residual(ind_receiver, :));
                    
                    
                case 'backprojection'
                    
                    % get the residual for the current receiver
                    res_signal = residual(ind_receiver, :);
                    
            end
            
            if any(res_signal)
                
                
                %% ============================================================
                % ACTION OF THE ADJOINT OF THE JACOBIAN ON THE PERTURBED PRESSURE
                % OR RESIDUAL (J^T x P_{res}), or (J^T x \delta p_{res})
                %==============================================================
                % Apply backprojection on the signal on the receiver for getting
                % the adjoint pressure field on the grid points
                switch para.task
                    case 'backprojection'
                        [pressure_grid_adjoint_receiver, ~] = calc_pressure(approx_pressure,...
                            parameters_grid_adjoint{ind_receiver}, caustic_number_adjoint{ind_receiver},....
                            nan_grid_binary_adjoint{ind_receiver}, res_signal, [], [], [],...
                            'time-reversal');
                        
                    otherwise
                        [pressure_grid_adjoint_receiver, ~] = calc_pressure(approx_pressure,...
                            parameters_grid_adjoint{ind_receiver}, caustic_number_adjoint{ind_receiver},....
                            nan_grid_binary_adjoint{ind_receiver}, res_signal, [], [], [],...
                            'normal');
                end
                
                
                switch para.task
                    case {'gradient', 'hessian'}
                        
                        coeff_adjoint = 1;
                        
                    case{'backprojection'}
                        
                        % compute the squared amplitude of the slowness vector
                        % for the current emitter-receiver pair times
                        % the perturbed direction of the adjoint ray
                        % emanated from the current receiver from
                        % perturbation to the position of the
                        % receiver times the absoute of the
                        % mean amplitude  
                        angle = calcDirectionalAngle(directions_grid{ind_emitter}(:, 1:2).',...
                            - directions_grid_adjoint{ind_receiver}(:, 1:2).').';
                        
                        % get the coefficients associated with the weigting
                        % function diagonalising the Hessian matrix
                        % coeff_adjoint = 4 ./ sound_speed.^2 .* cos(angle/2).^2; 
                        coeff_adjoint = 4./sound_speed(nan_grid_binary_adjoint{ind_receiver}) .* cos(angle/2).^2;                       
                        coeff_adjoint = coeff_adjoint .* directions_grid{ind_emitter}(:, 3) .*...
                            directions_grid_adjoint{ind_receiver}(:, 3);
                        
                end
                
                % add the backprojected pressure field from the current receiver to
                % the backprojected pressure for all the receivers
                pressure_grid_adjoint = pressure_grid_adjoint...
                    + coeff_adjoint .* pressure_grid_adjoint_receiver;
            end
            
        end
        
        
        % Mutiply the forward pressure field by the adjoint
        % pressure field on each grid point
        pressure_forward_adjoint =  pressure_grid .* pressure_grid_adjoint;
        
        switch para.task
            case {'gradient', 'hessian'}
                
                % Multiply the product of the forward and adjoint pressure
                % field by $\Upsilon(c)$ on each grid point, and add the
                % resulting update direction on the grid points for the current
                % emitter to update direction for all the emitters
                sound_speed_update = sound_speed_update...
                    + (-2 ./ sound_speed.^3) .* ((pressure_forward_adjoint * omega_squared')...
                    + (tan(pi*y/2) + 1i) * sound_speed .* absorption .*...
                    (pressure_forward_adjoint * omega_frac_plus'));
                
                
            case 'backprojection'
                
                % Multiply the product of the forward and adjoint pressure
                % field by $\omega \Upsilon(c)$ on each grid point, and add the
                % resulting update direction on the grid points for the current
                % emitter to update direction for all the emitters
                %wave_number_update = wave_number_update + 1/2 * ...
                %    sum(pressure_forward_adjoint .* ...
                %    sound_speed .*...
                %    (omega_squared./sound_speed + (tan(pi*y/2) - 1i) * absorption .* omega_frac) ./...
                %    (omega_squared./sound_speed + 2 * tan(pi*y/2) * absorption .* omega_frac), 2);
                
                % get the matrix of the absoroption map obeying the frequency power law with
                % exponent power -1, i.e., \alpha_0 \omega^(y-1)
                alpha1 = absorption .* omega_frac_minus;
                
                % update the slowness square for the current emitter
                slowness_squared_update =  slowness_squared_update + ...
                    sum(pressure_forward_adjoint .* ...
                    (1./sound_speed + tan(pi*y/2) * alpha1) .* ...
                    (1./sound_speed + y * tan(pi*y/2) * alpha1) ./ ...
                    ((1./sound_speed + (tan(pi*y/2) + 1i) * alpha1) .* omega), 2);
        end
        
        
    end
    
end


if calc_adjointfield
    
    
    switch para.task
        
        case 'hessian'
            
            % get the sound speed update direction
            sound_speed_update_direction = sound_speed_update;
            
        case 'gradient'
            
            % Multiply the computed sound speed update direction
            % by the inverse of the grid area for accounting for spatial sampling
            % and the step length
            sound_speed_update_direction = step_length *...
                sound_speed_update;
            
        case  'backprojection'
            
            % Multiply the computed sound speed update direction
            % by the spacing of the angular frequencies for accounting for
            % integration over angular frequencies, by the inverse of the grid area
            % for accounting for spatial sampling and $1/(2pi)^3$ times the step length
            sound_speed_update_direction = step_length * omega_spacing/...
                ((2*pi)^3) *  slowness_squared_update;
    end
    
else
    
    sound_speed_update_direction = [];
    
end


% multiply the objective function by 1/2 by convention
objective_function = 1/2 * objective_function;


% compute the time for running the m-file function
computation_time = toc(tstart);


end