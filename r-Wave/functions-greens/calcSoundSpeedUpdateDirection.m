function [sound_speed_update_direction, objective_function, computation_time,...
    forward_signal, jacobian_signal] = calcSoundSpeedUpdateDirection(...
    data_freq, source_freq, time_array, omega, parameters_grid, parameters_grid_adjoint,...
    parameters_receiver, nan_grid_binary, nan_grid_binary_adjoint, caustic_number,...
    caustic_number_adjoint, receiver_order, caustic_receiver, binary_data_level,...
    perturbation, directions_grid, directions_grid_adjoint, sound_speed, absorption,...
    mask, random_residual, omega_spacing, step_length, varargin)
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
%       data_freq      - a num_emitter x num_freq x num_receiver matrix
%                        representing the frequency-domain ultrasound data
%                        for the chosen frequencies
%       source_freq     - the source in the frequency domain, if the
%                         ultarsound time traces have been deconvolved from
%                         the source, the source should be set 1 for all
%                         frequencies.
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
%       binary_data_level    - a num_emitter x 1 cell array, each containing
%                              a num_receiver x 1 vector representing the
%                              signals associated with the accepted
%                              emitter-receiver pairs for the current
%                              frequency. The accepted emitter-receivers
%                              pairs for each freqnency level are chosen
%                              based on a minimum distance between the
%                              emitter-receiver pairs and removing
%                              outliers.
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
%                           this variable is set an empty variable.
%       omega_spacing      - the spacing [rad s^{-1}] of the angular frequency
%       detec_radius       - the radius [m] of the detection ring
%       step_length        - the step length for computing the update direction
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
%       do_scaling                 - Boolean controlling whether scaling is
%                                    enforced on data and solution spaces or not.
%                                    For simulation studies for which the amplitude
%                                    of the source is known, scaling is not required.
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
% para.max_rays = 1000;
para.greens_domain = 'frequency';
para.greens_method = 'analytic';
para.do_scaling = false;
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
            'the information for the direction of the rays must be given.'])
    end
end


tstart = tic;

% get the size of the time arrays
num_time = length(time_array);

% get the number of emitters
num_emitter = length(parameters_grid);

% get the number of receivers
num_receiver = length(parameters_grid_adjoint);

% get the number of frequencies
num_freq = length(omega);

% check size of the observed data
if any(strcmp(para.task, {'forward', 'gradient','backprojection'}))

    if isempty(data_freq)

        % If the task is set forward or gradient, the frequency-domain
        % pressure data must be nonempty.
        error(['The frequency-domain pressure data must be nonempty,'...
            'if the task is set forward, gradient or backprojection.'])

    else

        % check the size of the given forward pressure data
        if any(size(data_freq) ~= [num_receiver, num_freq, num_emitter])
            error(['The pressure data size is not consistent with the cell arrays'...
                'containing the ray parameters.'])
        end

    end
end

if ~isempty(binary_data_level)

    % convert the binary data from a cell array to a three-dimensional matrix
    binary_data_level = cat(3, binary_data_level{:});

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

                % make the given random residual a cell array
                random_residual = squeeze(mat2cell(random_residual, num_receiver,...
                    size(random_residual, 2), ones(num_emitter, 1)));

                % If the random residual is nonempty, the purpose of implementing the
                % function is an inner product test, and therefore, the action of the
                % Jacobian on the perturbation does not need to be stored.
                para.record_output_signals(2) = true;

            end

            % if random residual is nonempty, the purpose is an inner product test,
            % and therefore, the action of the Jacobian matrix on
            % the perturbation must be stored.
            para.record_output_signals(2) = true;

        end


    case {'forward', 'gradient', 'backprojection'}

        if   ~(isempty(perturbation) && isempty(random_residual))

            % If the task is set 'forward' or 'gradient', both the perturbation
            % vector and random residual must be empty.
            error (['If the task is forward, gradient or backprojection,'...
                'both perturbation vector and the random residual'...
                'must be set empty inputs.'])

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

    % If the task is set 'forward', do not compute the adjoint field
    calc_adjointfield = false;

end

% Get the sound speed on the grid points inside the mask for image
% reconstruction
sound_speed = sound_speed(mask);

% get the power law exponent
y = absorption.power;

if size(parameters_grid{1}, 2) < 4

    % If parameters_grid has less than four columns, absorption map is
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

        % squared angular frequencies
        omega_squared = omega.^2;

        % fractional power plus 1  (y+1) of angular frequencies
        omega_frac_plus = omega.^(y+1);

        % fractional power minus 1 (y-1) of angular frequencies
        omega_frac_minus = omega.^(y-1);

    case {'time'}

        Error('Not implemented yet!')

end

%%=========================================================================
% FOR THE BACKWARD RAYS, CHANGE THE CELL ARRAYS TO SINGLE-FORMAT MATRICES
%==========================================================================
% create a 3D logic-format matrix for the Nan grid points
nan_grid_binary_adjoint = logical(cat(3, nan_grid_binary_adjoint{:}));

% create a 3D single-format matrix for the rays' parameters on the grid points
parameters_grid_adjoint = single(cat(3, parameters_grid_adjoint{:}));

% create a 3D single-format matrix for the accumulated rays' caustics on the grid points
caustic_number_adjoint = single(cat(3, caustic_number_adjoint{:}));

% create a 3D single-format matrix for the accumulated rays' caustics on the grid points
directions_grid_adjoint = single(cat(3, directions_grid_adjoint{:}));




% Define a handle function for approximating the pressure field using the
% Green's formulae
approx_pressure = @(source_pulse, time_delays, geom_attenuation, absorption,...
    caustic_number, source_mode) approxPressureGreens(source_pulse, time_delays,...
    geom_attenuation, absorption, caustic_number, omega, y, [], source_mode);


% Define a handle function for approximating the pressure field on the grid points and
% the receivers using the Green's function
calc_pressure = @(approx_pressure, parameters_grid, caustic_number, nan_grid_binary,...
    pressure_source, parameters_receiver, caustic_receiver, receiver_order,...
    source_mode) calcPressureGreens(approx_pressure, parameters_grid, caustic_number,...
    nan_grid_binary, pressure_source, source_mode, parameters_receiver,...
    caustic_receiver, receiver_order, num_freq, num_receiver);

% Initialise the update direction of the sound speed by zero
sound_speed_update = 0;

% Initialise the squared slowness update of the sound speed by zero
slowness_squared_update = 0;

% Allocate a cell array for the forward signals on the grid points
pressure_grid = cell(num_emitter, 1);

% Allocate a cell-array for the forward signals on the receivers
pressure_receiver = cell(num_emitter, 1);

if para.record_output_signals(2)

    % Allocate a matrix for the action of the Frechet derivative (Jacobian matrix) on the perturbation
    % of the sound speed, i.e., $J \delta c$. This will be used for testing the Frechet derivative
    % and its adjoint using an inner product test. However, for solving the linearised subproblem,
    % the actions $J^T J \delta c$  are computed by matrix-vector products without the need for
    % storage of $J \delta c$, so this variable will not be stored.
    jacobian_signal = cell(num_emitter, 1);

else

    % allocate an empty variable if the action of the Frechet derivative
    % on the sound speed perturbation is not stored.
    jacobian_signal = [];

end


parfor (ind_emitter = emitter_indices, para.num_workers)
    %   for ind_emitter = emitter_indices

    %% ============================================================
    % COMPUTE FORWARD PRESSURE FIELD
    %==============================================================
    % display the current emitter for which the update direction of
    % the sound speed is computed.
    disp(['Computing the forward portion of update direction of sound'...
        'speed for emitter:' num2str(ind_emitter)]);

    % approximate the forward Green's function (or pressure field) on the
    % grid points and on the receivers, given the rays' parameters on the
    % grid points and the receivers.
    switch para.task
        case {'forward','gradient'}

            % Compute the forward Green's function (or pressure field) on the
            % grid points and on the receivers
            [pressure_grid{ind_emitter}, pressure_receiver{ind_emitter}] = calc_pressure(approx_pressure,...
                parameters_grid{ind_emitter}, caustic_number{ind_emitter},...
                nan_grid_binary{ind_emitter}, source_freq,...
                parameters_receiver{ind_emitter}, caustic_receiver{ind_emitter},...
                receiver_order{ind_emitter}, 'normal');


        case 'backprojection'

            % Compute the forward Green's function on receivers and the
            % inversed forward Green's function on the grid points
            [pressure_grid{ind_emitter}, pressure_receiver{ind_emitter}] = calc_pressure(approx_pressure,...
                parameters_grid{ind_emitter}, caustic_number{ind_emitter},...
                nan_grid_binary{ind_emitter}, source_freq,...
                parameters_receiver{ind_emitter}, caustic_receiver{ind_emitter},...
                receiver_order{ind_emitter}, 'time-reversal');

        case 'hessian'

            % Compute the forward Green's function only on the receivers
            [pressure_grid{ind_emitter}, ~] = calc_pressure(approx_pressure,...
                parameters_grid{ind_emitter}, caustic_number{ind_emitter},...
                nan_grid_binary{ind_emitter}, source_freq, [],[],[],...
                'normal');
    end

end

%%=========================================================================
% CLEAR THE PARAMETERS NOT REQUIRED
%==========================================================================
% clear the parameters that are not required blow
clear parameters_grid caustic_number parameters_receiver

if any(strcmp(para.task, {'forward', 'gradient', 'backprojection'}))

    if isempty(random_residual)

        % convert the binary data from a cell array to a three-dimensional matrix
        pressure_receiver = cat(3, pressure_receiver{:});

        % apply the binary data on the computed pressures on the receivers
        pressure_receiver = binary_data_level .* pressure_receiver;

        if para.do_scaling

            % Compute the scaling factor for amplitudes
            % The scaling factor is a matrix of size num_receiver x num_freq
            % For excluded receivers, it will have value 0/0=nan.
            scaling_factor = dot(pressure_receiver, data_freq, 3)./...
                dot(pressure_receiver, pressure_receiver, 3);

        else

            % set the cslaing factor a ones matrix
            scaling_factor = ones(num_receiver, num_freq);

        end

        % compute the frequency-domain residual, the discrepancy of the scaled
        % computed frequency-domain pressure on the receivers and the frequency-domain
        % measured data
        residual = (scaling_factor .* pressure_receiver)...
            - (binary_data_level .* data_freq);

        % calculate the objective function, the squared L2 norm of the residual
        % Nans correponding to the excluded receivers are not accounted for.
        objective_function = norm(residual(isfinite(residual)))^2;

        % make the three-dimensional matrix for the residual a cell array
        residual = squeeze(mat2cell(residual, num_receiver, num_freq,...
            ones(num_emitter, 1)));

        if strcmp(para.task, 'backprojection')
            scaling_factor = 1./scaling_factor;
        end

    end


else


    if strcmp(para.task, 'hessian')

        % allocate empty cell array
        residual = cell(num_emitter, 1);

    end

    % allocate ones for scaling factors, if the purpose
    % is an inner product test
    scaling_factor = ones(num_receiver, num_freq);

end

if para.record_output_signals(1)

    % get the forward pressure field recorded on the receivers
    forward_signal = pressure_receiver;
else

    % allocate an empty variable if the forward pressures on the receivers are not stored
    forward_signal = [];

end

%%=========================================================================
% CLEAR THE PARAMETERS NOT REQUIRED
%==========================================================================
% clear the parameters that are not required blow
clear parameters_grid caustic_number pressure receiver





if calc_adjointfield

   parfor (ind_emitter = emitter_indices, para.num_workers)
         %  for ind_emitter = emitter_indices


        switch para.task
            
            case'hessian'

        % compute p(x, x_e) \Upsilon x
        % perturbation
        res_signal_emitter = - 2 * perturbation.* pressure_grid{ind_emitter} .* ...
            (( 1 ./ sound_speed.^3 .* omega_squared +...
            (tan(pi * y/2) + 1i)  ./ sound_speed.^2 .*...
            absorption.* omega_frac_plus));


        % make the nans zero
        % nan values may occur for the pressure field approximated on the grid
        % points that lie outside the mask for ray tracing. So the interpolation
        % of the rays' parameters onto these grid points will be nan.
        res_signal_emitter(~isfinite(res_signal_emitter))= 0;

            case 'backprojection'
                
         % get the rays' directions for the current emitter    
        directions_grid_emitter = directions_grid{ind_emitter};
        
        end

        disp(['Computing the backprojected portion of update direction of the sound'...
            'speed for emitter:' num2str(ind_emitter)]);


        % Initialise the adjoint pressure field for the current emitter
        pressure_grid_adjoint = 0;

        % Allocate a matrix for the action of the Frechet derivative
        % (Jacobian matrix) on the sound speed perturbation for the
        % current emitter
        %  jacobian_signal_receiver = zeros(num_receiver, num_frequency);

        for ind_receiver = 1 : num_receiver


            % get the nan binaries
            nan_binaries = nan_grid_binary_adjoint(:,:, ind_receiver);

            switch para.task

                case 'hessian'

                    %% ============================================================
                    % ACTION OF THE JACOBIAN ON THE PERTURBATION (J \delta c_l)
                    %==============================================================
                    % Compute the action of the Frechet derivative on
                    % the sound speed perturbation vector, i.e. the product
                    % of the Jacobian matrix on the perturbation vector.

                    % Multiply the forward pressure at every grid point by the scattering
                    % potential, if the binary data is 1 for the current
                    % emitter-receiver pair


                    if all(binary_data_level(ind_receiver,:, ind_emitter))

                        
                        % Approximate the frequency-domain forward pressure sccatterd
                        % on the grid points and propagating to the receiver
                        [res_signal, ~] = calc_pressure(approx_pressure,...
                            parameters_grid_adjoint(:,:,ind_receiver),...
                            caustic_number_adjoint(:,:,ind_receiver),....
                            nan_binaries, res_signal_emitter,...
                            [], [], [], 'normal');


                        % compute the action of the Frechet derivative
                        % (Jacobian matrix) on the perturbation vector for
                        % the current emitter-receiver pair
                        res_signal = sum(res_signal);

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
                        % product test, get the residual for the current
                        % emitter-receiver pair
                        res_signal = random_residual{ind_emitter}(ind_receiver, :);

                    end

                    % compute the complex conjugate of the signal
                    if strcmp(para.task, 'hessian')

                        if isempty(random_residual)

                            % If the random residual is nonempty, and
                            % therefore, the purpose is computing the sound speed update,
                            % transform the random residual to the frequency domain
                            res_signal = scaling_factor(ind_receiver, :) .* conj(res_signal);

                        end

                    end

                case 'gradient'

                    % get the complex conjugate of the residual for the current receiver
                    res_signal = conj(residual{ind_emitter}(ind_receiver, :));


                case 'backprojection'

                    % get the residual for the current receiver
                    res_signal = residual{ind_emitter}(ind_receiver, :);

            end

            if any(res_signal)


                %% ============================================================
                % ACTION OF THE ADJOINT OF THE JACOBIAN ON THE PERTURBED PRESSURE
                % OR RESIDUAL (J^T x P_{res}), or (J^T x \delta p_{res})
                %==============================================================
                % compute the backprojected Green's function (or pressure
                % field) emitted from the receiver for the current
                % emitter-receiver pair on the grid points
                switch para.task
                    case 'backprojection'

                        % compute the inversed backprojected Greens' function
                        % or pressure on the grid points
                        [pressure_grid_adjoint_receiver, ~] = calc_pressure(approx_pressure,...
                            parameters_grid_adjoint(:,:,ind_receiver), caustic_number_adjoint(:,:,ind_receiver),....
                            nan_binaries, res_signal, [], [], [],...
                            'time-reversal');

                    otherwise

                        % compute the backprojected Greens' function or pressure
                        % on the grid points
                        [pressure_grid_adjoint_receiver, ~] = calc_pressure(approx_pressure,...
                            parameters_grid_adjoint(:,:,ind_receiver), caustic_number_adjoint(:,:,ind_receiver),....
                            nan_binaries, res_signal, [], [], [],...
                            'normal');
                end


                switch para.task

                    case {'gradient', 'hessian'}

                        % get 1 for the coefficient enforced on the
                        % backprojected field
                        coeff_adjoint = 1;

                    case{'backprojection'}

                        % compute the squared amplitude of the slowness vector
                        % for the current emitter-receiver pair times
                        % the perturbed direction of the backprojected ray
                        % initialized on the current emitter and receiver after a
                        % perturbation to the position of the emitter and
                        % receiver
                        angle = calcDirectionalAngle(directions_grid_emitter(:, 1:2).',...
                            - directions_grid_adjoint(:,1:2,ind_receiver).').';


                        % get the coefficients associated with the weighting
                        % function diagonalising the Hessian matrix
                        coeff_adjoint = 4./sound_speed .* cos(angle/2).^2 .*...
                             abs(directions_grid_emitter(:, end)) .*...
                            abs(directions_grid_adjoint(:, end, ind_receiver));
                        
                end


                % add the backprojected pressure field from the current receiver to
                % the backprojected pressure for all the receivers
                pressure_grid_adjoint = pressure_grid_adjoint +...
                    coeff_adjoint .* ...
                    (scaling_factor(ind_receiver, :) .*...
                    pressure_grid_adjoint_receiver);
            end

        end


        % Mutiply the forward pressure field by the adjoint
        % pressure field on each grid point
        pressure_forward_adjoint =  pressure_grid{ind_emitter} .*...
            pressure_grid_adjoint;

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
                ((2*pi)^3) * slowness_squared_update;

    end

else

    sound_speed_update_direction = [];

end

if strcmp(para.task, 'hessian')

    % give an empty variable if the taks is computing the Hessian matrix
    objective_function = [];
else
    % multiply the objective function by 1/2 by convention
    objective_function = 1/2 * objective_function;
end

% compute the time for running the m-file function
computation_time = toc(tstart);


end