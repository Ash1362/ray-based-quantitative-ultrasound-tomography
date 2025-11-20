function [img, out, para] = reconstructGreensImage(data, time_array, recon_grid,...
    emitter, receiver, sound_speed_water, img, ray_initial_angles, simulation_prop,...
    varargin)
%RECONSTRUCTGREENSIMAGE reconstructs the image of the sound speed using the
%Green's approach
%
%
% DESCRIPTION:
%           reconstructGreensImage reconstructs the image of the sound
%           speed using the Green's approach
%
%
% USAGE:
%
%
% INPUTS:
%       data              - the ultrasound time series measured from the
%                           object inside water. This is a matrix of size
%                           num_receiver x num_time x num_emitter with
%                           num_emitter the number of emitters, num_time
%                           the number of time instants for measurement
%                           and num_receiver the number of receivers.
%       time_array        - a time array [s] of size 1 x num_time on which the
%                           the UST time series are measured.
%       recon_grid        - the computational grid for image reconstruction.
%       emitter           - a struct which defines the properties of the
%                           excitation as follows: This includes the
%                           fields 'positions', 'pulse', 'pulse_duration',
%                           and 'shot_time'
%       emitter.positions - 2/3 x num_emitter array of Cartesian position
%                           of the center of the emitter objects
%       emitter.pulse     - a vector of size 1 x num_te with num_te <= num_time.
%       emitter.shot_time - the first arrival of the excitation pulse, which
%                           can be negative [sec]
%       receiver          - a struct which includes the fields
%       receiver.positions - 2/3 x N array of Cartesian position of the
%                           centre of the receiver objects. For rotating
%                           ultrasound systems, the position of the receivers
%                           is changed for different excitations. For the rotating case,
%                           this is a cell array with length num_angle (the number of rotation
%                           angles) with each cell a 2/3 x num_transducer array of
%                           Cartesian position of transducers for each
%                           angle. (num_receiver = num_angle x num_transducer)
%       sound_speed_water   - the sound speed [m/s] in only water
%       img                 - a matrix representing the initial guess of the sound speed
%                           image for iterative image reconstruction
%       ray_initial_angles - the num_emitter x 1 cell array each corresponding
%                           to an emitter, and containing the initial angle
%                           of the rays initialised from the emitter and linked
%                           to all recevers
%      simulation_prop      - a struct array containing the information and properties of
%                             the k-wave simulation. This struct includs the fields:
%      'PML'                - a  1 x dim vector of the number of PML layers
%                             along each dimension
%      'f_max'              - the maximum frequency [Hz] supported by the
%                             computational grid
%      'detec_radius'       - the radius of the detection ring (surface)
%      'z_offset'           - the z offset [m] between the grid for k-Wave simulation
%                             and the grid for image reconstruction. The
%                             image reconstruction is done on a detection surface
%                             (ring) with centre at the origin. (empty for
%                             2D case)
%      'data_path'          - the data path for storing the simulated UST
%                             data and the computed TOFs
%      'x'                  - the x position of the grid for the k-Wave
%                             simulation
%      'y'                  - the y position of the grid for the k-Wave
%                             simulation
%      'z'                  - the z position of the grid for the k-Wave
%                             simulation
%      'sound_speed'        - the ground truth sound speed distribution [m/s]
%                             (used only for evaluation purposes, not for
%                              image reconstruction)
%      'sound_speed_ref'    - the sound speed [m/s] in only water
%      'alpha_coeff'        - the absorption coefficient [dBMHz^{-y} cm^{-1}]
%      'alpha_power'        - the exponent power of the acoustic absorption used for
%                             the k-Wave simulation
%      't_array'            - the time array [s] containing time instants on
%                             which the synthetic data for object in water
%                             and only water are computed.
%
%
%
% OPTIONAL INPUTS:
%      'num_worker_pool'   - The number of workers for parallel programming
%      'reconstruct_image' - Boolean contrlling whether the image is reconstructed
%                            or loaded
%                            (default: true)
%      'deconvolve_source' - Boolean controlling whether the measured signals
%                            are deconvolved from the pressure source or
%                            not. (default:true)
%      'deconvolution_parameter' - get the regularisation parameter for
%                            deconvolution (default : 0)
%       'do_scaling'       - Boolean controlling whether scaling is
%                            enforced on data and solution spaces or not.
%                            For simulation studies for which the amplitude
%                            of the source is known, scaling is not required.
%      'absorption_map'    - the known absorption map, which can be 'true',
%                            'homogeneous' or 'none' (deafult:true)
%      'noise_level'         - an estimated level of noise [dB] for the data
%                              (default : 40)
%      'f_range'               - a 1 by 2 matrix for the frequnecy [Hz] range
%                              This is chosen based on the frequency spectrum
%                              of the excitation pulse or measured signals.
% the supported by the grid for data simulation,
%                              or the transducers
%                              (default: [0.2e6,1.72e6])
%      'optimisation_approach' - The optimisation approach, which can be
%                              either 'hessian' or 'backprojection'.
%                              (default: 'backprojection')
%      'emittersegments_to_nworkers' - The number of segments for dividing emitters to the number of
%                            cpu workers.  For avoiding memory leakage during the
%                            parallel programming, the emitters are first
%                            divided into a number of segments, and then
%                            the segment are consecutively run over the threads.
%      'z_pos_height'      - the maximum z position [m] of the grid points, if
%                            the reconstruction geometry is set 'real'
%                            (Default: 1e-2)
%      'mask_coeff'        - a coefficient which defines the normalised radius of
%                            the circular binary mask with radius = mask_coef * detec_radius
%                            for image reconstruction. Here, detec_radius
%                            is the radius of the detection surface, which
%                            must be chosen greater 1. (Default : 1.03)
%      'minimum_distance_coeff' - a scalar coefficient indicating the minimum distance
%                            of the transducers for including in the image
%                            reconstruction. The emitter-receiver pairs with distances
%                            less than minimum_distance = minimum_distance_coeff * detec_radius
%                            The measured signals associated with close emitter-receiver
%                            pairs may be deteriorated by the effects of directivity
%                            of the transducers. In addition, rays linking the
%                            close emitter-receiver pairs only travel inside water.
%      'sound_speed_direction_constraint' - A maximum range in which the update direction
%                            of the sound speed can be.
%      'smoothing_window_size_initial' - the initial size of the smoothing window applied
%                            on the updates of the refractive index.
%      'raylinking_method' - the method for ray linking. This can be either
%                            'Regula-Falsi' or 'Secant' for 2D case, or
%                            'Quasi-Newton' for 3D case.
%      'max_raylinking_iter' - the maximum permissible number of iterations for ray
%                             linking
%      'raylinking_threshold' - the threshold distance [m] for stopping the
%                               ray linking iterations. If the distance
%                               between the interception point of the ray
%                               by the detection surface (ring) and the centre
%                               of the receiver becomes less than
%                               raylinking_threshod, the ray linking
%                               algorithm will be stopped.
%      'auxiliary_angle'      - the perturbation enforced on the initial angle of
%                               the linked rays for tracing auxiliary rays
%      'raytogrid_spacing'    - the spacing of the sampled points along the ray
%                               to the grid spacing (Default: 1/2 (2D case),
%                               1 (3D case))
%      'gridtoray_interp'     - the approach for interpolating from grid
%                               points to ray
%      'attenuation_geom_method' - the approach for computing the
%                               geometrical portion of the attenuation
%      'absorption_coeff_homogeneous' - a scalar value representing a homogeneous
%                               absorption coefficient for the whole breast
%      'do_priorconditioning' - Boolean controlling whether the priorconditining
%                            is included or not. (default :false)
%      'filter_data'          - Boolean controlling whether the data is
%                               filtered in the frequency domain or not.
%                               (default : false)
%      'cut_off'              - A scalar value representing the cut off factor
%                               of a filter used for removing ringing artefact
%                               from the sound speed update direction computed
%                               after soving each linearized subproblem. (default:1)
%      'order'                - A scalar value representing the order of a filter
%                               used for removing ringing artefact from the sound speed
%                               update direction computed after soving each linearized
%                               subproblem. (default:inf)
%      'outliers_fraction'    - the fraction of the emitter-receiver pairs
%                               excluded from image reconstruction for a specific
%                               frequency, because their magnitudes are outliers.
%                               (default: 0.01)
%      'multiply_data_window' - Boolean controlling whether the time series
%                               are multiplied by an expoentially-varying
%                               time window or not (default : false)
%      'do_hom_greens'       - A boolean controlling whether the Green's
%                               functions are approximated using
%                               homogeneous Green's functions or not
%
%
%
% OUTPUTS:
%      img                    - the reconstructed image of the sound speed
%      recon_grid             - the computational grid for image reconstruction
%      tof_discrepancy        - the num_reciver x num_emitter matrix of the
%                               the discrepancy of the time-of-flights
%                               for the object-in-water data and only-water
%                               data
%      ray_initial_angles     - the initial angle of the linked (optimal) rays for all
%                               emitter-receiver pairs after ray linking applied to the
%                               last update of the sound speed
%      out                    - a struct containing the stored data of
%                               the iterative reconstruction. This
%                               struct includes the fields:
%      'norm_res'             - the L2 norm of res updates for the cg algorithm
%                               and the L2 norm of update directions of the sound speed
%                               for the steepest descent algorithm at each inner iteration
%      'residual_norm'        - the first component of norm_res
%      'num_rays'             - the number of rays for ray linking between
%                               the emitter-receiver pairs
%      'system_matrix_time'   - the cpu time for construction of the system
%                               matrix
%      'update_time'          - the cpu time for solving each linear subproblem
%       'imgs'                - all the reconstructed images of the sound speed
%                               after solving the linearised subproblems
%                               If the sound_speed_ground_truth is given,
%                               this struct can also includes the fields:
%       'relative_error'      - the L2 norm of the discrepancy of the
%                               updates of the sound speed and the ground
%                               truth over the L2 norm of the disrepancy of
%                               water and the ground truth time 100.
%        'rmse'               - the rmse of the reconstructed images with
%                               respect to the ground truth.
%        para                 - the optional inputs used for the image
%                               reconstruction

%
%
% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 05.10.2022
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian

para.num_worker_pool = 16;
para.reconstruct_image = true;
para.deconvolve_source = false;
para.deconvolution_parameter = 0;
para.do_scaling = false;
para.absorption_map = 'true';
para.noise_level = 40;
para.f_range = [0.2e6, 1.72e6]; % [0.2e6, 1.65e6];
para.optimisation_approach = 'backprojection';
para.raylinking_method = 'Secant';
para.max_raylinking_iter = 1000;
para.raylinking_threshold = 1e-6;
para.auxiliary_angle = pi/(2 *180);
para.raytogrid_spacing = 1;
para.shot_angles = 'raylinking';
para.attenuation_geom_method = 'auxiliary';
para.auxiliary_method = 'paraxial';
para.absorption_coeff_homogeneous = 0.5;
% para.stopping_tol = 1e-3;
para.minimum_distance_coeff = 0.7368;
para.sound_speed_direction_constraint = inf;
para.do_priorconditioning = false;
para.filter_data = false;
para.cut_off = 1;
para.order = inf;
para.outliers_fraction = 0.01;
para.multiply_data_window = false;
para.do_hom_greens = false;
para.do_direction_angle = true;

% get the grid spacing [m], which is assumed the same along all the
% Cartesian coordinates
grid_spacing = recon_grid.x_vec(2) - recon_grid.x_vec(1);


% choose the smoothing window size for the supported  grid spacing [m]
% Give an error if the grid spacing is not set 2mm, 1mm or 0.5 mm
if  grid_spacing < (5e-4)-(1e-10)
    error('The grid spacing smaller than 0.5mm is not necessary for a ray-based method.')
elseif  abs(grid_spacing - (5e-4))< 1e-10
    para.smoothing_window_size_initial = 13;
elseif  abs(grid_spacing - (7.5e-4))< 1e-10
    para.smoothing_window_size_initial = 9;
elseif abs(grid_spacing - (1e-3))< 1e-10
    para.smoothing_window_size_initial = 7;
elseif abs(grid_spacing - (2e-3))< 1e-10
    para.smoothing_window_size_initial = 5;
elseif grid_spacing > (2e-3)+(1e-10)
    error('The grid spacing must not be larger than 2mm.')
else
    error('The supported grid spacings are 2mm, 1mm, 0.75mm or 0.5mm.')
end

% get the number of the dimensions
dim = size(emitter.positions, 1);

% A scalar for choosing the ratio of the radius of a circular (spherical) mask
% representing the Region-Of-Interest (ROI) for image reconstruction to the
% radius of the detection ring (surface)
switch dim
    case 2
        para.mask_coeff = 0.95;
    case 3
        para.mask_coeff = 0.83;
end

% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


switch para.optimisation_approach

    case 'hessian'

        % the number of discretised frequencies within the frequency range
        para.num_frequency = 180;          % (Default = 180)

        % the number of discretised frequencies for each frequency level (linearisation)
        para.num_frequency_per_level = 2;  % (Default = 2)

        % The number of frequency levels. Each frequency level contains
        % a fixed number of evenly spaced sampled frquencies.
        if para.do_hom_greens
            para.num_frequency_level = 25;     % (Default: 25)
        else
            para.num_frequency_level = 60;     % (Default: 50)
        end


        if para.noise_level > 32

            % the number of inner iterations for solving each linearised
            % subproblem when the estimated SNR for UST data is more than
            % 32 dB.
            para.num_iter = 15;

        elseif para.noise_level > 27

            % the number of inner iterations for solving each linearised
            % subproblem when the estimated SNR for UST data is between
            % 27-32 dB.
            para.num_iter = 12;

        else

            % the number of inner iterations for solving each linearised
            % subproblem when the estimated SNR for UST data is smaller
            % than 27 dB.

            % reducing the number of iterations for solving the
            % linearised subproblems is an implicit way for
            % enforcing regularisation.
            para.num_iter = 10;    % (Default: 10)
        end

        if para.do_hom_greens

            % If the sound speed is assumed homogeneous (only-water),
            % the number of inner iterations should be multiplied by 2/3,
            % leading to an early stopping criterion for solving each linearized
            % subproblem, and implicitly acts as a regularizer.
            para.num_iter = round(2/3 * para.num_iter);

        end


    case {'gradient','backprojection'}

        % the number of discretised frequencies within the frequency range
        para.num_frequency = 180;         % (Default = 190)

        % the number of discretised frequencies for each frequency level (linearisation)
        para.num_frequency_per_level = 2; % (Default = 2)

        % The number of frequency levels. Each frequency level contains
        % a fixed number of evenly spaced sampled frquencies.
        % get the number of frequency levels used in the iterative algorithms
        % A single update direction will be computed for each frequency level.
        para.num_frequency_level = 90;    % (Default = 95)

end



switch para.optimisation_approach
    case 'gradient'
        error('The gradient approach for optimisation is deprecated.')
        para.step_length = 1e-6;
    case 'hessian'
        para.step_length = (2.75e1)/grid_spacing; % default : 2.8e1/grid_spacing;
    case 'backprojection'
        para.step_length = 1.2e-1; % default : 1.2e-1;
end


if para.do_priorconditioning

    % get the number of iterations for solving the linearised problem for
    % inverting the preconditioning matrix
    para.preconditioning_iter = 100;

    % get the stooping tolerance for solving the linearised problem for
    % inverting the preconditioning matrix
    para.preconditioning_tol = 1e-3;

    % get the smoothing parameter for total-variation preconditioning
    para.smoothing_parameter = 1e-3;

    % the scalar parameter added to the preconditioning matrix for ensuring
    % the positive definiteness
    para.gamma = 1e-2;

end

% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


% get the frequency indices used for image reconstruction
para.frequency_level = 1 : para.num_frequency_level;


% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% Using the optimisation approach based on backprojection, the direction of
% the rays must be computed.
do_get_direction = strcmp(para.optimisation_approach, 'backprojection');

if ~strcmp(para.optimisation_approach, 'backprojection')

    % if the optimisation approach is not 'backprojection', the perturbed
    % rays' direction becuase of perturbation to initial positions are not
    % requitred.
    para.do_direction_angle = false;
end
      
% For checking the progress during running the m-file function, the reconstructed
% images after each 10 frequency levels are visualized and stored in a
% temporary folder.
% Add the reconstructed image after each 10 frequency levels in a temporary
% folder.
image_directory = 'results/results_temp/';

% make the path for saving the images
makeDirectory(image_directory);

%%=========================================================================
% SET THE PARAMETERS FOR THE SIMULATED DATA
% =========================================================================
% get the number of emitters
num_emitter = size(emitter.positions, 2);

% get the number of receivers
if iscell(receiver.positions)

    % get the number of receivers for the rotational setting
    num_receiver = size(receiver.positions{1}, 2);
else

    % get the number of receivers for the fixed setting
    num_receiver = size(receiver.positions, 2);
end


% get the initial guess
if length(unique(img(:))) == 1
    initial_guess =  'homogeneous';
else

    initial_guess = 'heterogeneous';

    if isempty(ray_initial_angles)
        error(['Using a heterogeneous image as the initial guess, the initial angles of the bent rays must be given.'...
            'This should be available through image reconstruction using time-of-flights.'])
    end

end


% check the smoothing window size
if rem(para.smoothing_window_size_initial, 2) == 0 && para.smoothing_window_size_initial > 3
    error('The smoothing window size must be an odd natural number, if it is larger than 3.')
end


% get the Boolean controlling whether the auxiliary rays are traced, or not
switch para.attenuation_geom_method
    case 'auxiliary'
        trace_auxiliary_ray = true;
    case 'raylinked'
        trace_auxiliary_ray = false;
end


% choose an angular frequency range [rad/s]
frequency_range = 2 * pi * para.f_range;

% the vector of angular frequencies [rad/s]
omega = linspace(frequency_range(1),...
    frequency_range(2), para.num_frequency);

% get the spacing of the angular frequency
omega_spacing = mean(diff(omega));

% choose the indices of angular frequencies for the first frequency level
frequency_indices = 1:para.num_frequency_per_level;

% get the radius of the detection surface (ring)
detec_radius = norm(emitter.positions(:, 1));

% get the minimum distance [m] for accepting the emitter-receiver pairs
% included in the image reconstruction
minimum_distance = para.minimum_distance_coeff * detec_radius;

% get the grid-to-ray interpolation approach
gridtoray_interp = 'Bspline';


%% ========================================================================
% GET THE DISTANCE BETWEEN EMITTER-RECEIVER PAIRS
%==========================================================================
% give an empty variable to the field rotation_indices, if it does not
% exist, i.e., the position of receivers are fixed with changes in
% excitations.
if ~isfield(emitter, 'rotation_indices')

    % add an empty field 'rotation_indices' to struct emitter.
    emitter.rotation_indices = [];

end

% Calculate the distance between emitter-receiver pairs
%distance_emitters_receivers = calculateDistanceEmitterReceiver(emitter.positions,...
%    receiver.positions, emitter.rotation_indices);

% get the spacing [m] of the sampled points on the rays
ray_spacing = grid_spacing * para.raytogrid_spacing;


% display the chosen ray-to-grid spacing
disp(['The ray to grid spacing is:' num2str(para.raytogrid_spacing)])


%%=========================================================================
% PREPROCESS THE ULTRASOUND DATA AND TRANSFORM IT TO THE FREUENCY DOMAIN
%==========================================================================

% get the optional inputs for preprocessing the ultrasound data set and
% transform it to the frequency domain
preprocess_data_args = {'num_workers', para.num_worker_pool, 'deconvolve_source',...
    para.deconvolve_source, 'filter_data', para.filter_data, 'outliers_fraction',...
    para.outliers_fraction, 'multiply_data_window', para.multiply_data_window};

% preprocess the ultrasound data set and transform it to the frequency
% domain
[data, source_freq, binary_data] = preprocessGreensData(data,...
    time_array, omega, emitter, receiver, sound_speed_water, minimum_distance,...
    preprocess_data_args{:});



%% ========================================================================
% INTERPOLATE THE MEDIUM'S PROPERTIES ONTO THE GRID FOR RAY TRACING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the sound speed, and absorption coefficient maps, and interpolate them
% from the grid for data simulation onto the grid for ray tracing

if ~isempty(simulation_prop)

    switch recon_grid.dim
        case 2

            if isfield(simulation_prop, 'sound_speed')
                % interpolate the sound speed from the grid for data simulation
                % onto the grid for ray tracing and computing the Green's function
                sound_speed_phantom = interpn(simulation_prop.x, simulation_prop.y,...
                    simulation_prop.sound_speed, recon_grid.x, recon_grid.y);
                sound_speed_phantom(~isfinite(sound_speed_phantom)) = sound_speed_water;
            else
                sound_speed_phantom = [];
            end

            if isfield(simulation_prop, 'alpha_coeff')

                % interpolate the absorption coefficient from the grid for data simulation
                % onto the grid for ray tracing and computing the Green's function
                absorption_phantom = interpn(simulation_prop.x, simulation_prop.y,...
                    simulation_prop.alpha_coeff, recon_grid.x, recon_grid.y);
                absorption_phantom(~isfinite(absorption_phantom)) = 0;

            elseif isfield(simulation_prop, 'alpha_coeff_real')

                % get an estimated homogeneous absorption coefficient
                % [dB MHz^{-y} cm^{-1}] within the breast
                absorption_phantom = simulation_prop.alpha_coeff_real;
            else

                absorption_phantom = 0;

            end


        case 3

            if isfield(simulation_prop, 'sound_speed')
                % interpolate the sound speed from the grid for data simulation
                % onto the grid for ray tracing and computing the Green's function
                sound_speed_phantom = interpn(simulation_prop.x, simulation_prop.y,...
                    simulation_prop.z - simulation_prop.z_offset, simulation_prop.sound_speed,...
                    recon_grid.x, recon_grid.y, recon_grid.z, 'spline');
                sound_speed_phantom(:, :, end) = sound_speed_water;
            else
                sound_speed_phantom = [];
            end

            % interpolate the absorption coefficient
            if isfield(simulation_prop, 'alpha_coeff')

                % interpolate the absorption coefficient from the grid for data simulation
                % onto the grid for ray tracing and computing the Green's function
                absorption_phantom = interpn(simulation_prop.x, simulation_prop.y,...
                    simulation_prop.z - simulation_prop.z_offset, simulation_prop.alpha_coeff,...
                    recon_grid.x, recon_grid.y, recon_grid.z, 'spline');
                absorption_phantom(:, :, end) = 0;

            elseif isfield(simulation_prop, 'alpha_coeff_real')

                % get an estimated homogeneous absorption coefficient
                % [dB MHz^{-y} cm^{-1}] within the breast
                absorption_phantom = simulation_prop.alpha_coeff_real;
            else

                absorption_phantom = 0;

            end

    end


    %% ========================================================================
    % GET THE ABSORPTION COEFFICIENT (ONLY FOR SIMULATION STUDIES)
    %==========================================================================
    % For experimental studies, the field 'alpha_coeff_real' is used, and
    % therefore, abosrption_map is not affected by para.absorption_map
    if isfield(simulation_prop, 'alpha_coeff')

        switch para.absorption_map
            case 'homogeneous'

                % set the absorption coefficient inside the breast homogeneous and
                % 0.5 dB MHz^{-y} cm^{-1}
                absorption_phantom(absorption_phantom > 0) = para.absorption_coeff_homogeneous;

            case 'none'

                % ignore the acoustic absorption in the image reconstruction
                absorption_phantom = 0;
        end

    end


    % make a struct array containing the absorption coefficient map and the
    % true given exponent factor
    absorption = struct('coeff', absorption_phantom, 'power', simulation_prop.alpha_power);


else

    % give an empty variable for the sound speed of phantom
    sound_speed_phantom = [];

    % ignore abosrption, if simulation_prop is not given.
    absorption = struct('coeff', 0, 'power', 1.4);

end


% get the absorption in terms of neper
absorption.coeff_neper = db2neper(absorption.coeff, absorption.power);

%%=========================================================================
% GET THE INITIAL GUESS
%==========================================================================
switch initial_guess

    case 'heterogeneous'

        if para.frequency_level(1) == 1

            % get the initial sound speed image as an image reconstructed by
            % the inversion approach using Time-of-flights
            % apply the smoothing window on the initial guess
            img = smoothField(img, para.smoothing_window_size_initial - 2, sound_speed_water);

        end

    case 'homogeneous'

        % get a homogeneous sound speed for all grid points as the initial
        % guess
        img = mean(img(:)) * ones(recon_grid.size);

        % if the initial guess for the sound speed is set homogeneous and water,
        % make the initial angle of the rays an empty variable, so the
        % rays will be straight lines which are directly intercepted by the receivers
        ray_initial_angles = [];

end


%%=========================================================================
% GET THE BINARY MASK FOR IMAGE RECONSTRUCTION
%==========================================================================
switch dim
    case 2

        % get a binary mask for image reconstruction. The sound speed for the grid
        % points outside this mask will be assumed the background sound
        % speed(water). This does not need to be the same as the mask used for image
        % reconstruction using Time-of-flights
        mask_recon = recon_grid.x.^2 + recon_grid.y.^2 < ((para.mask_coeff-0.05) * detec_radius)^2; % (Default = para.mask_coeff-0.05)

        % get another mask, a bit larger than the mask for image reconstruction, for ray tracing
        mask_raytracing = recon_grid.x.^2 + recon_grid.y.^2 < (para.mask_coeff * detec_radius)^2;

    case 3
        NotImpErr
end

% set the sound speed outside the binary mask the background sound speed
% (water), remind that img is the sound speed image reconstructed by
% Time-of-Flights
img(~mask_recon) = sound_speed_water;



if strcmp(para.optimisation_approach, 'backprojection' )

    % Using the backprojection (Hessian-free) approach for reconstructing
    % the sound speed, the Green's functions are never computed by
    % assumption of homogeneous medium.
    para.do_hom_greens = false;
end

% initialize the sound speed image used for approximating the Green's function
% parameters. If the optional input 'do_hom_greens' is set true, this sound speed
% map is always used.
img_update = sound_speed_water * ones(recon_grid.size);

% get the Cartesian position of the grid points inside ROI mask
grid_x = recon_grid.x(mask_recon);
grid_y = recon_grid.y(mask_recon);


% Define a handle function for approximating the parameters of the Green's
% function on the grid points and receivers
calc_parameters = @(ray_position_emitter, ray_time_emitter, ray_absorption_emitter,...
    ray_position_emitter_left, ray_position_emitter_right, sound_speed_relative,...
    field_mode) calcParametersGreens(....
    ray_position_emitter, ray_time_emitter, ray_absorption_emitter, ray_position_emitter_left,...
    ray_position_emitter_right, sound_speed_relative, [grid_x, grid_y], ray_spacing,...
    detec_radius, para.mask_coeff, field_mode, do_get_direction, []);


% get the number of frequencies for each linearisation
num_frequency_per_level = para.num_frequency_per_level;

% get the number of linearisations
num_frequency_level = para.num_frequency_level;

if strcmp(para.optimisation_approach, 'hessian')

    % get the number of iterations
    num_iter = para.num_iter;
end

%% ========================================================================
% DEFINE THE OPTIONAL INPUTS
%==========================================================================
% define the optional inputs for ray tracing
ray_args = {'nworker_pool', para.num_worker_pool, 'interp_method', gridtoray_interp,...
    'raylinking_method', para.raylinking_method, 'raylinking_initialisation',...
    'Local', 'max_iter', para.max_raylinking_iter, 'raytogrid_spacing', para.raytogrid_spacing,...
    'auxiliary_method', para.auxiliary_method, 'auxiliary_ray', trace_auxiliary_ray,...
    'reference_angle', para.auxiliary_angle, 'max_num_points_factor', 1.3,...
    'do_direction_angle', para.do_direction_angle,...
    'varepsilon', para.raylinking_threshold, 'smoothing_window_size',...
    para.smoothing_window_size_initial, 'angular_frequency_centre', mean(omega(frequency_indices))};

% define the optional inputs for computing the update of the sound speed.
% update_args{4} will be replaced during the loop, so make sure the order
% of elaments in update_args must not be changed.
update_args = {'num_workers', para.num_worker_pool, 'task', para.optimisation_approach,...
    'deconvolve_source', para.deconvolve_source, 'deconvolution_parameter',...
    para.deconvolution_parameter, 'do_scaling', para.do_scaling,...
    'record_output_signals', [false, false]};

%%=========================================================================
% ALLOCATE VARIBALES FOR ITERATIVE IMAGE RECONSTRUCTION
%==========================================================================
out = [];
% allocate a matrix for storing the gradient of the objective function inside the mask for image reconstruction
% at the different frequency levels
out.gradient_obj = zeros(nnz(mask_recon), num_frequency_level);

% allocate a vector for the CPU time spent for calculation of the trajetory of rays
% and the accumulated parameters at the frequency levels
out.ray_cpu_time = zeros(num_frequency_level, 1);

% allocate a vector for the norm of the gradient at the frequency levels
out.gradient_norm = zeros(num_frequency_level,1);

% allocate a vector for the norm of the update direction at the frequency levels
out.update_norm = zeros(num_frequency_level,1);

% allocate a vector for the CPU time for calculating the L2 norm of the gradient
% at the frequency levels
out.gradient_cpu_time = zeros(num_frequency_level, 1);

% allocate a vector for the objective function at the frequency levels
out.objective_function = zeros(num_frequency_level, 1);


if strcmp(para.optimisation_approach, 'hessian')

    % allocate a vector for the norm of the residuals associated with the
    % the linearised objective function
    out.norm_res = zeros(num_frequency_level, num_iter);

    % allocate a vector for the CPU time for calculation of the action of
    % Hessian $H$ on the update perturbation of the sound speed  $\delta c_l$, ie. $H \delta c_l$
    out.hessian_cpu_time = zeros(num_frequency_level, num_iter);

end


if ~isempty(sound_speed_phantom)

    % allocate a matrix for the relative error of the reconstructed images at
    % the frequency levels and the iteration at each frequency level
    switch para.optimisation_approach
        case {'gradient', 'hessian'}

            out.relative_error = zeros(num_frequency_level, num_iter);

        case {'backprojection'}

            out.relative_error = zeros(num_frequency_level, 1);
    end
end


% allocate cell arrays
parameters_grid = cell(1, num_emitter);
parameters_receiver = cell(1, num_emitter);
nan_grid_binary = cell(1, num_emitter);
caustic_number = cell(1, num_emitter);
caustic_receiver = cell(1, num_emitter);
receiver_order = cell(1, num_emitter);
directions_grid = cell(1, num_emitter);
parameters_grid_adjoint = cell(1, num_receiver);
nan_grid_binary_adjoint = cell(1, num_receiver);
caustic_number_adjoint = cell(1, num_receiver);
directions_grid_adjoint = cell(1, num_receiver);

% display the step length
disp(['The step length, which is fixed for all frequency levels, is:'...
    num2str(para.step_length)]);

for frequency_level = para.frequency_level


    % store the reconstructed sound speed image at the last level.
    % This will be used for postprocessing the computed update direction
    % at this level.
    img_previous = img;


    if para.do_hom_greens

        % the medium on which the ray tracing is performed is assumed
        % homogeneous only water, so the traced ray are always straight
        % rays initialized on emitters towards the receivers, and ray
        % linking is not required.
        ray_initial_angles = [];

    else

        % get the sound speed image used for approximating the Green's function
        % parameters. If the optional input 'do_hom_greens' is set false,
        % this sound speed map is set the current sound speed map.
        img_update = img;

    end

    % get the frequency indices
    frequency_indices = ((frequency_level-1) * num_frequency_per_level + 1) :...
        (frequency_level * num_frequency_per_level);

    % if the smoothing window size is smaller than 3, the smoothing is
    % gradually reduced by an increase in the frequency
    if para.smoothing_window_size_initial < 3.1

        ray_args{end-2} = ray_args{end-2} - 0.1;

    else

        % update the frequency for ray tracing. The frequency affects the trajectory of rays
        % through acoustic dispersion effects.
        ray_args{end} = mean(omega(frequency_indices));

       
        % set the smoothing window for ray tracing. Start with larger
        % smoothing window and reduce it by progressing from low
        % frequencies to high frequencies
        % Default:
        if ray_args{end}/(2*pi)< 4.5e5
            ray_args{end-2} = 13;
        elseif ray_args{end}/(2*pi)< 7.0e5
            ray_args{end-2} = 11;
        elseif ray_args{end}/(2*pi)< 9.5e5
            ray_args{end-2} = 9;
        else
            ray_args{end-2} = 7;
        end

        % If grid spacing is set 0.75mm, increase the size of smoothing
        % wwindoe by 4 for all frequency levels.
        if  abs(grid_spacing - (7.5e-4))< 1e-10
            ray_args{end-2} = ray_args{end-2} + 4;
        end


    end


    % display the current angular frequency interval
    disp(['The current frequency interval is:'...
        num2str( 1/(2*pi) * [omega(frequency_indices(1)),...
        omega(frequency_indices(num_frequency_per_level))], '%1.5e')]);

    % Compute the trajectory of the rays for the chosen frequency, as well as
    % the required parameters on the rays, and if the ray linking method is set
    % 'Secant', update the initial angles of the rays such that the optimal initial angles
    % after ray linking at the frequency level 'n-1' is used as the initial guess for ray
    % linking at the frequency level 'n'
    [~, num_rays, ray_initial_angles, ray_position, ray_time, ray_absorption, rayspacing_receiver,...
        ray_position_left, ray_position_right, adjoint_ray_position_left,...
        adjoint_ray_position_right, out.ray_cpu_time(frequency_level)] =...
        computeRaysParameters(recon_grid, sound_speed_water./img_update,...
        absorption.coeff, emitter.positions, receiver.positions,...
        ray_initial_angles, mask_raytracing, emitter.rotation_indices, [], ray_args{:});



    % Compute the parameters on the grid points and the receivers for the
    % forward field
    parfor (ind_emitter = 1 : num_emitter, para.num_worker_pool)
        %for ind_emitter = 1 : num_emitter


        [parameters_grid{ind_emitter}, parameters_receiver{ind_emitter},...
            nan_grid_binary{ind_emitter}, caustic_number{ind_emitter},...
            receiver_order{ind_emitter}, caustic_receiver{ind_emitter},...
            directions_grid{ind_emitter}, ~] = calc_parameters(...
            ray_position{ind_emitter}, ray_time{ind_emitter}, ray_absorption{ind_emitter},...
            ray_position_left{ind_emitter}, ray_position_right{ind_emitter},...
            img_update(mask_recon)/sound_speed_water, 'forward');

        ray_position_left{ind_emitter} = [];
        ray_position_right{ind_emitter} = [];


    end




    if strcmp(para.optimisation_approach, 'backprojection')

        % get the perturbation to the angle of the forward rays because of
        % a perturbation to the initial position (position of emitter), by
        % assuming traight ray, i.e., no refraction
        % This derivative is computed using the chain rule.
        [directions_grid] = calcDerivativeAngleToInitialPosition(...
            directions_grid, dim, [grid_x, grid_y], emitter.positions.',...
            'adjacent');

    end


    % Reverse the forward rays (trajectories and the parameters on the rays) for computing the adjoint rays
    [ray_position, ray_time, ray_absorption, ray_position_left, ray_position_right] =...
        calcRaysAdjoint(ray_position, ray_time, ray_absorption,...
        adjoint_ray_position_left, adjoint_ray_position_right, ray_spacing,...
        rayspacing_receiver);



    % Compute the parameters on the grid points and the receivers for the
    % backward field
    parfor (ind_receiver = 1 : num_receiver, para.num_worker_pool)
        %for ind_receiver = 1 : num_receiver

        [parameters_grid_adjoint{ind_receiver}, ~, nan_grid_binary_adjoint{ind_receiver},...
            caustic_number_adjoint{ind_receiver}, ~, ~, directions_grid_adjoint{ind_receiver}, ~] =...
            calc_parameters(ray_position{ind_receiver}, ray_time{ind_receiver},...
            ray_absorption{ind_receiver}, ray_position_left{ind_receiver},...
            ray_position_right{ind_receiver}, img_update(mask_recon)/sound_speed_water,...
            'adjoint');


        ray_position{ind_receiver} = [];
        ray_time{ind_receiver} = [];
        ray_absorption{ind_receiver} = [];
        ray_position_left{ind_receiver} = [];
        ray_position_right{ind_receiver} = [];


    end




    if strcmp(para.optimisation_approach, 'backprojection')

        % get the perturbation to the angle of the backward rays because of
        % a perturbation to the initial position (position of emitter), by
        % assuming traight ray, i.e., no refraction
        % This derivative is computed using the chain rule.
        [directions_grid_adjoint] = calcDerivativeAngleToInitialPosition(...
            directions_grid_adjoint, dim, [grid_x, grid_y], receiver.positions.',...
            'adjacent');

    end


    if strcmp(para.optimisation_approach, 'backprojection')

        % get the squared slowness
        slowness_squared = 1./(img.^2);

    end


    % get the binary data for the current frequency level

    % allocate a cell array of length num_emitter for the binary data map
    % for the current level
    binary_data_level = cell(num_emitter, 1);

    for ind_emitter = 1: num_emitter

        % get a binary vector whose component is zero for bad linked rays
        binary_linked_rays = num_rays{ind_emitter} < para.max_raylinking_iter-1;

        % get the binary data vector for the current frequency level and
        % emitter after enforcing the zeros corresponding to the
        % bad-linked rays
        binary_data_level{ind_emitter} = binary_data{ind_emitter}(:, frequency_indices) .*...
            binary_linked_rays;

    end



    switch para.optimisation_approach
        case {'gradient','hessian'}
            update_args{4} = 'gradient';
        case 'backprojection'
            update_args{4} = 'backprojection';
    end

    % compute the gradient of the objective function, ie. the action of the
    % adjoint of Frechet derivative on the residual
    [out.gradient_obj(:, frequency_level), out.objective_function(frequency_level),...
        out.gradient_cpu_time(frequency_level), ~ , ~] = ...
        calcSoundSpeedUpdateDirection(squeeze(data(:, frequency_indices, :)), source_freq, time_array,...
        omega(frequency_indices), parameters_grid, parameters_grid_adjoint, parameters_receiver,...
        nan_grid_binary, nan_grid_binary_adjoint, caustic_number, caustic_number_adjoint, receiver_order,...
        caustic_receiver, binary_data_level, [], directions_grid, directions_grid_adjoint,...
        img, absorption, mask_recon, [], omega_spacing, para.step_length, update_args{:});

    % get the real part of the gradient
    out.gradient_obj(:, frequency_level)= real(out.gradient_obj(:, frequency_level));

    % compute the L2 norm of the gradient
    out.gradient_norm(frequency_level) = norm(out.gradient_obj(:, frequency_level));

    % display the results for the current iteration
    disp(['The frequency level (linearised subproblem) is:' num2str(frequency_level)])
    disp(['The objective function is:' num2str(out.objective_function(frequency_level), '%1.2f')])
    disp(['The norm of the gradient is:' num2str(out.gradient_norm(frequency_level), '%1.5e')])
    disp(['The CPU time for computing the rays was:' num2str(out.ray_cpu_time(frequency_level)/60, '%1.2f') 'min'])
    disp(['The CPU time for calculating gradient of the sound speed was:' num2str(out.gradient_cpu_time(frequency_level)/60, '%1.2f') 'min'])


    if any(strcmp(para.optimisation_approach, 'hessian'))

        % switch the task for the sound speed update to 'hessian'
        update_args{4} = para.optimisation_approach;

        % initial the cg residual by the minus gradient
        residual_cg = - out.gradient_obj(:, frequency_level);


        % get the preconditioned res
        if para.do_priorconditioning

            % define the handle function for computing the gradient of the total variation
            preconditioner = @(field) calcSpatialVariationMatrix(field,...
                img, mask_recon, recon_grid.size, para.smoothing_parameter,...
                para.gamma, 1);

            % enforce the preconditioner on the residual through a cg algorithm
            [preconditioned_residual_cg, ~, ~] = pcg(preconditioner, residual_cg, ...
                para.preconditioning_tol, para.preconditioning_iter);

            % apply the parameter reducing the ill-condition of the
            % Hessian matrix
            % preconditioned_residual_cg =  preconditioned_residual_cg - para.lambda * img(mask_recon);

        else
            preconditioned_residual_cg = residual_cg;
        end


        % initialise the update direction by the residual
        update_direction = preconditioned_residual_cg;


        for iter = 1 : num_iter

            if ~isempty(sound_speed_phantom)

                % compute the relative error
                out.relative_error(frequency_level, iter) = 100 * norm(img(mask_recon)...
                    - sound_speed_phantom(mask_recon))/...
                    norm(sound_speed_phantom(mask_recon) - sound_speed_water);
            end

            % get the norm of the residual
            out.norm_res(frequency_level, iter) = norm(residual_cg);

            % calculate the action of the Hessian matrix on the update direction
            [hessian_by_update_direction, ~, out.hessian_cpu_time(frequency_level, iter), ~ , ~] =...
                calcSoundSpeedUpdateDirection(squeeze(data(:, frequency_indices, :)), source_freq, time_array,...
                omega(frequency_indices), parameters_grid, parameters_grid_adjoint, parameters_receiver,...
                nan_grid_binary, nan_grid_binary_adjoint, caustic_number, caustic_number_adjoint, receiver_order,...
                caustic_receiver, binary_data_level, update_direction, directions_grid, directions_grid_adjoint,...
                img, absorption, mask_recon, [], omega_spacing, para.step_length, update_args{:});


            % get the real part of the action of the Hessian matrix on the sound speed
            % update
            hessian_by_update_direction = real(hessian_by_update_direction);

            % update the cg step length 1
            step_length_cg = (residual_cg' * preconditioned_residual_cg)...
                / (update_direction' * hessian_by_update_direction);

            % update the reconstructed image of the sound speed inside the mask
            img(mask_recon) = img(mask_recon) +  step_length_cg * update_direction;

            % store the current cg residual
            residual_cg_previous = residual_cg;

            % store the current preconditioned cg residual
            preconditioned_residual_cg_previous = preconditioned_residual_cg;

            % update the cg residual
            residual_cg = residual_cg_previous - step_length_cg * hessian_by_update_direction;


            % get the preconditioned res
            if para.do_priorconditioning

                % enforce the preconditioner on the update direction through a cg algorithm
                [preconditioned_residual_cg, ~, ~] = pcg(preconditioner, residual_cg,...
                    para.preconditioning_tol, para.preconditioning_iter);
            else
                preconditioned_residual_cg = residual_cg;
            end


            % update the cg step length 2
            beta = (residual_cg' * preconditioned_residual_cg)....
                / (residual_cg_previous' * preconditioned_residual_cg_previous);

            % update the update direction
            update_direction = preconditioned_residual_cg + beta * update_direction;

            % display the results for the current inner iteration
            disp(['The level is:' num2str(frequency_level) '. The inner iteration is:' num2str(iter)])
            disp(['The relative error is:' num2str(out.relative_error(frequency_level, iter), '%1.2f') '%'])
            disp(['The norm of the residual is:' num2str(out.norm_res(frequency_level, iter), '%1.5e')])
            disp(['The CPU time for computing the action of the Hessian on the upadte direction is:'...
                num2str(out.hessian_cpu_time(frequency_level, iter)/60, '%1.2f') 'min'])

            % check th stopping criterion for nonmonotone convergence
            % check if the norm of residual for all the three last
            % solutions of the current linearised subproblem is greater
            % than the last previous update before those three updates
            if iter > 3
                if   all(out.norm_res(frequency_level, iter-3)< out.norm_res(frequency_level, iter-2:iter))
                    break
                end
            end

        end

    else

        % get the update direction of the wavenumber as the minus gradient
        update_direction = - out.gradient_obj(:, frequency_level);

        % update the reconstructed image of the wave number inside the mask
        slowness_squared(mask_recon) = slowness_squared(mask_recon) + update_direction;

        % update the reconstructed image of the sound speed inside the mask
        % from the updated wavenumber
        img = slowness_squared.^(-1/2);

        if ~isempty(sound_speed_phantom)

            % compute the relative error
            out.relative_error(frequency_level) = 100 * norm(img(mask_recon) - sound_speed_phantom(mask_recon))/...
                norm(sound_speed_phantom(mask_recon) - sound_speed_water);

            disp(['The relative error is:' num2str(out.relative_error(frequency_level), '%1.2f') '%'])
        end

    end


    % remove rininging artefact from the updated direction
    img_direction = ringingRemovalFilt(recon_grid.x_vec.', recon_grid.y_vec.',...
        img-img_previous, sound_speed_water, 1/(2*pi) * mean(omega(frequency_indices)),...
        para.cut_off, para.order);

    % update the sound speed image
    img = img_previous + img_direction;

    % update the sound speed image
    img = img_previous + max(min(img-img_previous, + para.sound_speed_direction_constraint),...
        - para.sound_speed_direction_constraint);

    % compute the norm of the update direction at the frequency levels
    out.update_norm(frequency_level) = norm(img(mask_recon)-img_previous(mask_recon));

  
    if rem(frequency_level, 10) < 0.5
        figure;imshow(img, [min(img(:)), max(img(:))]);
        axis image;colorbar;

        % For checking the progress during running, the reconstructed images
        % after each 10 frequency levels are visualized and stored in a
        % temporary folder.
        saveas(gca, [image_directory 'fig' '_' 'level_'...
            num2str(frequency_level)], 'fig');
    end



end

out.ray_initial_angles = ray_initial_angles;


end