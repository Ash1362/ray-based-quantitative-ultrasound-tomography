function [ray_position, ray_time, ray_absorption, ray_position_left, ray_position_right,...   
adjoint_ray_position_left, adjoint_ray_position_right, ray_spacing,...
rayspacing_receiver, relative_discrepancy_absorbing, relative_discrepancy_nonabsorbing] =...
    validateGreensHeterogeneous(data_object, data_water, time_array, emitter, receiver,...
     simulation_prop, tof_discrepancy, plot_directory, varargin)
%VALIDATEGREENSHETEROGENEOUS compares the pressure field approximated using Green's function
% with the pressure field simulated by the k-Wave for a heterogeneous medium.
%   
%
% DESCRIPTION:
%             validateGreensHeterogeneous compares the approximated pressure field
%             by ray approximation to heterogeneous Green's function with the 
%             pressure field simulated by k-Wave on the receivers after being 
%             produced by emitters and travelling across a digital breast
%             phantom
%      
%
% USAGE:
%     
%
% INPUTS:
%       data_object       - the ultrasound time series measured from the 
%                           object inside water. This is a matrix of size
%                           num_receiver x num_time x num_emitter with
%                           num_emitter the number of emitters, num_time
%                           the number of time instants for measurement
%                           and num_receiver the number of receivers.
%       data_water        - the ultrasound time series measured from water.
%                           This is a matrix of size num_receiver x num_time
%                           x num_emitter with num_emitter the number of emitters,
%                           num_time the number of time instants for measurement
%                           and num_receiver the number of receivers.
%       time_array        - a time array of size 1 x num_time 
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
%      tof_discrepancy     - the num_receiver x num_emitter matrix of the
%                            discrepancy of time-of-flights beween the
%                            object-in-water data and only-water data.
%                            This will be used, if the plan is to remove 
%                            the scatted waves from the measured (k-Wave
%                            simulated data).
%      plot_directory      - the directory for saving the plots
%  
%      
 
% OPTIONAL INPUTS:
%      'num_worker_pool'   - the number of workers for parallel programming
%      'save_plots'        - Boolean controlling whether the plots are saved
%                            or not. (default : false)
%      'deconvolve_source' - Boolean controlling whether the measured signals
%                            are deconvolved from the pressure source or
%                            not. (default:false)
%      'deconvolution_parameter' - get the regularisation parameter for
%                            deconvolution (default : 1e-6)
%      'choose_receiver'   - the approach for choosing the receiver for
%                           plotting. The receiver can be chosen by index 
%                           ('index'), or alternatively by distances to
%                           the chosen emitter ('distance')
%      'emitter_index'    - the index of a chosen emitter for analysis
%      'receiver_index'   - the index of a chosen receiver for analysis
%      'num_discretised_frequencies' - the number of discretised frequencies
%                           within the interval (default :50)
%      'grid_spacing'     - the grid spacing [m]
%      'z_pos_height'      - the maximum z position [m] of the grid points, if
%                            the reconstruction geometry is set 'real'
%                            (Default: 1e-2)
%      'mask_coeff'        - a coefficient which defines the normalised radius of
%                            the circular binary mask with radius = mask_coef * detec_radius
%                            for image reconstruction. Here, detec_radius
%                            is the radius of the detection surface, which
%                            must be chosen greater 1. (Default : 1.03)
%      'shot_angles'      - the method for choosing the initial angle,
%                           The initial angles can be determined by ray linking
%                           ('raylinking'), or can be set aligning the straight lines
%                           from emitter to all receivers ('even_spaced_angle')
%                           (default = 'raylinking') The second approach is
%                           only used as a benchmark.
%      'attenuation_geom_method' - the method for calculating the refraction-induced
%                           attenuation. This can be set 'auxiliary' or 'raylinked'
%                           (default = 'auxiliary')
%                           The latter is only used as benchmark for demonstrating
%                           the importance of ray linking. For the latter approach,
%                           the pressure is then linearly inpterpolated from the end point
%                           of the rays on the detection ring (surface) onto the reception
%                           points.
%      'auxiliary_method' - the approach for tracing auxiliary rays, which will
%                           be used for computing the geomterical attenuation.
%                           This can be set 'paraxial' or 'angle_perturbation'
%                           (default = 'paraxial').
%      'angle_perturb_auxiliary' - the perturbation to the initial angle of
%                           the ray for computing the amiplitude. It will be used,
%                           only if 'auxiliary_method' is set 'angle_perturbation'.
%                           (default : pi/(2*180) [rad])
%      'gridtoray_interp' - the approach for interpolation from grid to the sampled
%                           points on the ray (Default: 'Bspline')
%      'max_raylinking_iter' - the maximum number of iterations for ray
%                             linking
%      'raytogrid_spacing ' - the spacing of the sampled points along the
%                           rays over the grid spacing (default: 1) 
%      'remove_scattered_waves' - Boolean controlling whether the sacttered waves are
%                           removed from the k-Wave simulated signals or
%                           not. (default: false)
%      'filter_signals'    - Boolean controlling whether the time-series are filtered in the frequency
%                           domain, or not. (default : false)
%      'amplitude_plot_approach' - the approach for plotting amplitudes.
%                           This can be set 'maximum', or 'single-frequency'.
%      'frequency_single' - the interested single frequency for analysis

%       
%
% OUTPUTS:
%      ray_position                             - the num_emitter x 1 cell array
%                                                 for the position of the sampled
%                                                 points along the rays
%      ray_time                                 - the num_emitter x 1 cell array
%                                                 for the accumulated time delays
%                                                 on the sampled points along
%                                                 the ray
%      ray_absorption                           - the num_emitter x 1 cell array
%                                                 for the accumulated acoustic 
%                                                 absorption on the sampled
%                                                 points along the rays
%      ray_position_left                        - the position of the
%                                                 sampled points along the left
%                                                 auxiliary ray
%      ray_position_right                       - the position of the
%                                                 sampled points along the right auxiliary
%                                                 ray
%      adjoint_ray_position_left                - the position of the
%                                                 sampled points along the left auxiliary
%                                                 adjoint ray
%      adjoint_ray_position_right                - the position of the
%                                                 sampled points along the right auxiliary
%                                                 adjoint ray
%      ray_spacing                              - a scalar value representing the 
%                                                 spacing [m] of the
%                                                 sampled points along the
%                                                 rays
%      ray_spacing_receiver                     - the spacing [m] of the
%                                                 last two points along the linked rays
%                                                 for all emitter-receiver pairs.                                                
%      relative_discrepancy_absorbing           - the percentage relative
%                                                 discrepancy between the
%                                                 signals simulated by
%                                                 k-Wave and the signals
%                                                 approximated by the
%                                                 Green's function for the 
%                                                 object in water on all
%                                                 receivers after being
%                                                 produced by all emitters,
%                                                 when the acoustic absorption
%                                                 and dispersion are included.
%      relative_discrepancy_nonabsorbing        - the percentage relative
%                                                 discrepancy between the
%                                                 signals simulated by
%                                                 k-Wave and the signals
%                                                 approximated by the
%                                                 Green's function for the 
%                                                 object in water on all
%                                                 receivers after being
%                                                 produced by all emitters,
%                                                 when the acoustic absorption
%                                                 and dispersion are ignored.
   
%    
%
% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 10.08.2022
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian 

para.num_worker_pool = 16;
para.save_plots = false;
para.deconvolve_source = false;
para.deconvolution_parameter = 1e-6;
para.num_discretised_frequencies = 50;
para.grid_spacing = 1e-3;
para.z_pos_height = 1e-2;
para.mask_coeff = 1.05; 
para.shot_angles = 'raylinking';
para.attenuation_geom_method = 'auxiliary';
para.auxiliary_method = 'paraxial';
para.angle_perturb_auxiliary = pi/(2*180);
para.gridtoray_interp = 'Bspline';
para.max_raylinking_iter = 1000;
para.raytogrid_spacing = 1;
para.remove_scattered_waves = false;
para.filter_signals = false;
para.amplitude_plot_approach = 'single-frequency';
para.frequency_single = 1e6;
para.choose_receiver = 'index';
para.emitter_index = 1;

% get the number of the receivers
num_receiver = size(receiver.positions, 2);

% get the index of the receiver
para.receiver_index = round(num_receiver/2);

% get the number of the dimensions 
dim = size(emitter.positions, 1);

% get the preferred distance of the chosen emitter-receiver pairs in terms
% of the radius of the detection ring
switch dim
    case 2
        para.distance_range = 2 * [0.93, 0.95];
    case 3
        para.distance_range = 2 * [0.85, 0.87];
end
   
% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end



% choose the smoothing window size based on the grid spacing [m]
if para.grid_spacing < 1e-3
    error('The grid spacing smaller than 1mm is not necessary for a ray-based method.')
elseif para.grid_spacing < 1.5e-3
    para.smoothing_window_size = 7;
elseif para.grid_spacing < 2.01e-3
    para.smoothing_window_size = 5;
else
    error('The grid spacing must not be larger than 2mm.')
end

if rem(para.smoothing_window_size, 2) == 0 && para.smoothing_window_size > 3
    error('The smoothing window size must not be an even value, if it is larger than 3.')
end
 
if para.mask_coeff < 1
    error(['Using the absolute approach, the mask_coeff, which is the normalised radius'...
        'of the binary mask for image reconstruction, must be larger than the radius of'...
        'the detection surface(ring).'])
end

% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% get the index of the chosen emitter for plotting
emitter_index = para.emitter_index;

% get the index of the chosen receiver for plotting
receiver_index = para.receiver_index;


% display the number of emitters and receivers
disp(['The number of the chosen emitter is:' num2str(emitter_index)]);
disp(['The number of the chosen receiver is:' num2str(receiver_index)]);

% the method for initilisation of ray linking
raylinking_initialisation = 'Local';
if strcmp(raylinking_initialisation, 'Global')
    error('The Global approach for initialisation of the ray linking is not supported yet.')
end

% get the number of emitters
num_emitter = size(emitter.positions, 2);

% check the position of the transducers used for the k-Wave simulation 
% and image reconstruction is consistent.
if ~isempty(simulation_prop)
    switch simulation_prop.detection_geom
        case 'sphere'
            
            % get the radius of the detection surface (ring)
            detec_radius = norm(emitter.positions(:, 1));
            
            if abs(detec_radius - simulation_prop.detec_radius)> 1e-10
                error(['The position of the transducers for data simulation and'...
                    'image reconstruction is not consistent.'])
            end
        case 'cylinder'
            
            % get the radius of the detection surface (ring)
            detec_radius = norm(emitter.positions(1:2, 1));

            if abs(detec_radius - simulation_prop.detec_radius)> 1e-10
                error(['The position of the transducers for data simulation and'...
                    'image reconstruction is not consistent.'])
            end
            
    end
end



% choose a range for the angular frequency [rad/s]
frequency_range = 2 * pi * [0, simulation_prop.f_max];

% get the length of the time array
num_time = length(time_array);
        
% get the time spacing
dt = time_array(2) - time_array(1);
     
% get the angular frequency array supported by the grid for data
% simulation
omega = linspace(0, 2*pi * simulation_prop.f_max, para.num_discretised_frequencies);
        
% get the angular frequencies in the chosen frequency range
omega = omega(omega > frequency_range(1) &...
            omega <= frequency_range(2));
        
% get the pressure source in the frequency domain
pressure_source = discreteFourierTransform(emitter.pulse,...
            omega, time_array);
   
% get the original pressure source     
pressure_source_original = pressure_source;

if para.deconvolve_source
   pressure_source = 1;
end


if para.filter_signals
    
    % choose the cut-off frequencies for the filter, if the first component
    % is set zero, the filter will be low-pass. The upper bound of the
    % filter is chosen twice th maximum dfrequency supported by the grid
    % for the k-Wave simulation.
    cutoff_freq = [0, 2 * simulation_prop.f_max];
else
    cutoff_freq = nan;
end


switch para.shot_angles
    
    case 'raylinking'
        
        switch dim
            case 2
                
        % the method for ray linking
        raylinking_method = 'Regula-Falsi';
            case 3
           % the method for ray linking
        raylinking_method = 'Quasi-Newton';  
        end
        
        % the stopping threshold for ray linking
        raylinking_threshold = 1e-6;  
        
        % the initial angles for the rays
        initial_angles = [];
        
        % the query points for interpolation onto the receivers are not used,
        % so an empty variable
        interp_receiver = [];
        
    case 'even_spaced_angle'
        
        % This approach is only used as a benchmark for testing. Using this approach,
        % the rays are traced with equidistant initial angles, and ray linking is not
        % performed. First, by choosing an empty variable for initial angle of
        % the rays, the rays are initialised by initial angles aligning straight lines
        % from emitter to the receivers. By choosing the stopping criterion
        % infinity, the initial traced rays are accepted as the optimal
        % ray.
        
        % the method for ray linking
        raylinking_method = 'Secant';
        
        % the stopping threshold for ray linking
        raylinking_threshold = inf;
        
        % by setting the initial angles of the rays an empty variable (as an input for the
        % ray tracing function), the initial angles are chosen aligning straight lines
        % from each emitter to all the receivers.
        initial_angles = [];
        
        % the polar angle of the receivers (used as the the query points for interpolation
        % of the parameters from the end point of the rays onto the receivers)
        angle_receiver = (cart2pol(receiver.positions(1,:),...
            receiver.positions(2,:)))';
        
        % sort the query points in an ascending order, and their associated
        % indices
        [interp_receiver.angles, interp_receiver.indices] = sort(angle_receiver, 'ascend');
        
end

% Give an error message, if the shot angles are set 'even_spaced_angle', and the
% method for computing the geomterical portion of the attenuation is chosen 'raylinked'
if strcmp(para.shot_angles, 'even_spaced_angle') && strcmp(para.attenuation_geom_method, 'raylinked')
    error(['The method for computing the geometrical spreading portion of the attenuation'...
        ' must be set only auxiliary, if the rays are traced by evenly-spaced angles.']);
end


% get the Boolean controlling whether the auxiliary rays are traced, or not
switch para.attenuation_geom_method
    case 'auxiliary'
        trace_auxiliary_ray = true;
    case 'raylinked'
        trace_auxiliary_ray = false;
end


% display the used approach for computing the geometrical portion of the
% attenuation
disp(['The approach for computing the geometrical portion of the amplitude is: '...
     para.attenuation_geom_method])
 
 if strcmp(para.attenuation_geom_method, 'auxiliary')
disp(['The approach for computing the auxiliary rays is: '...
     para.auxiliary_method])
 end
 

% get the ray spacing [m]
ray_spacing = para.grid_spacing * para.raytogrid_spacing;


% a factor for the extension of the grid beyond the maxium position of the
% transducers along each Cartesian coordinate
switch para.gridtoray_interp
    case 'Bilinear'
        
        grid_expansion_coeff = 1 + 0.07 * 1000 * para.grid_spacing;
    case 'Bspline'
        
        % B-spline uses four adjacent grid points for interpolation, so the
        % grid is enlarged such that the adjacent grid points about the
        % target grid points do not exceed the edge of the computational
        % grid.
        grid_expansion_coeff = 1 + 0.07 * 1000 * para.grid_spacing;
        
end



%% ========================================================================
% FILTER THE MEASURED SIGNALS
%==========================================================================
% filter the data in frequency domain, if requested
if all(isfinite(cutoff_freq))
    
    cutoff_freq(~cutoff_freq) = [];
    
    % choose the method for filtering based on the chosen frequency intervals
    if length(cutoff_freq) < 2
        
        % if the first component is zero, apply the 'low-pass' filter
        filter_mode = 'low-pass';
        
    else
        
        % If the both components are finite, apply the 'band-pass' filter
        filter_mode = 'band-pass';
    end
    
    data_filter_args = {'Mode', filter_mode, 'nworker_pool', para.num_worker_pool};
    
    % filter the data using the given cut-off frequency
    data_object = filterPressureData(data_object, cutoff_freq, dt, data_filter_args{:});
    data_water = filterPressureData(data_water, cutoff_freq, dt, data_filter_args{:});
    
end


%% ========================================================================
% MAKE THE GRID FOR IMAGE RECONSTRUCTION
% =========================================================================
% the optional inputs for construction of the grid for image reconstruction
grid_args = {'grid_expansion', grid_expansion_coeff * ones(1, dim)};
if dim == 3
    grid_args = {grid_args{:}, 'reconstruction_geometry', 'real',...
        'z_pos_height', para.z_pos_height};
end

% make the grid for image reconstruction
recon_grid = makeReconstructionGrid(para.grid_spacing * ones(1, dim),...
    detec_radius, grid_args{:});

%% ========================================================================
% GET THE MASK FOR RAY TRACING
% =========================================================================
% update the mask for ray tracing
switch dim
    case 2
mask_raytracing = recon_grid.x.^2 + recon_grid.y.^2 < (para.mask_coeff * detec_radius)^2;
    case 3
mask_raytracing = recon_grid.x.^2 + recon_grid.y.^2 + recon_grid.z.^2 < (para.mask_coeff * detec_radius)^2;
end

%% ========================================================================
% INTERPOLATE THE MEDIUM'S PROPERTIES ONTO THE GRID FOR RAY TRACING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the sound speed, and absorption coefficient maps, and interpolate them
% from the grid for data simulation onto the grid for ray tracing
switch recon_grid.dim
    case 2
        
        % interpolate the sound speed from the grid for data simulation
        % onto the grid for ray tracing and computing the Green's function
        sound_speed_phantom = interpn(simulation_prop.x, simulation_prop.y,...
            simulation_prop.sound_speed, recon_grid.x, recon_grid.y);
        sound_speed_phantom(~isfinite(sound_speed_phantom)) = simulation_prop.sound_speed_ref;
        
        
        if isfield(simulation_prop, 'alpha_coeff')
            
            % interpolate the absorption coefficient from the grid for data simulation
            % onto the grid for ray tracing and computing the Green's function
            absorption_phantom = interpn(simulation_prop.x, simulation_prop.y,...
                simulation_prop.alpha_coeff, recon_grid.x, recon_grid.y);
            absorption_phantom(~isfinite(absorption_phantom)) = 0;
            
        else
            absorption_phantom = 0;
        end
        
        
    case 3
        
        % interpolate the sound speed from the grid for data simulation
        % onto the grid for ray tracing and computing the Green's function
        sound_speed_phantom = interpn(simulation_prop.x, simulation_prop.y,...
                simulation_prop.z - simulation_prop.z_offset, simulation_prop.sound_speed,...
            recon_grid.x, recon_grid.y, recon_grid.z);
        sound_speed_phantom(~isfinite(sound_speed_phantom)) = simulation_prop.sound_speed_ref;
        
        % interpolate the absorption coefficient
        if isfield(simulation_prop, 'alpha_coeff')
            
          
            % interpolate the absorption coefficient from the grid for data simulation
            % onto the grid for ray tracing and computing the Green's function
            absorption_phantom = interpn(simulation_prop.x, simulation_prop.y,...
                simulation_prop.z - simulation_prop.z_offset, simulation_prop.alpha_coeff,...
                recon_grid.x, recon_grid.y, recon_grid.z);
            absorption_phantom(~isfinite(absorption_phantom)) = 0;
            
            
        else
            absorption_phantom = 0;
        end
        
end

% compute the distance between the single emitter and receivers
distance_emitter_receivers = calculateDistanceEmitterReceiver(...
    emitter.positions, receiver.positions, []);

if para.remove_scattered_waves
   

% get the first arrival [s] of the excitation pulse
tof_excit_args = {'Method', 'Modified_AIC', 'Length_moving_windows',...
        [1.0, nan, 0.25, nan], 'Threshold', 0};
    
% get the first arrival time [s] of the pressure source
time_excit = timeOfFlightPickingEachSignal([emitter.pulse,...
    zeros(1, num_time-length(emitter.pulse))],...
        time_array, emitter.pulse_duration,...
        [dt, dt * find(abs(emitter.pulse)/max(abs(emitter.pulse))>0.2,...
        1, 'first')], tof_excit_args{:} );  
  

% display the excitation time
disp(['The excitation time is:' num2str(1e6 * time_excit) 'microseconds'])
    
% calculate the centre of the time window including the transmitted pulse
% using the calculated time-of-flights of the measured (simulated) signals
time_window_centre = time_excit + 1/simulation_prop.sound_speed_ref...
    * distance_emitter_receivers + tof_discrepancy(:);

% the edges for the time window [s] used for removing the scattered waves
time_window_interval = [-2, +6] * 1e-6;

end

% make the distances empty for 2D case, but keep it for computing Green's
% function for 3D case
if dim == 2
    distance_emitter_receivers = [];
end

%%=========================================================================
% CONVERT THE MEASURED (SIMULATED) SIGNALS TO CELL ARRAY
%==========================================================================
% convert the object-in-water and only-water data to cell arrays
data_object_cell = cell(num_emitter, 1);
data_water_cell = cell(num_emitter, 1);

for ind_emitter = 1 : num_emitter
    for ind_receiver = 1 : num_receiver
        
        % get the measured signal for the current emitter-receiver pair
        signal = data_object(ind_receiver, :, ind_emitter);
        
        if para.remove_scattered_waves
            
            % get a time window including the transmitted measured pulse
            time_window = ...
                time_window_centre((ind_emitter-1) * num_receiver + ind_receiver)...
                + time_window_interval;
            
            % make the signal ouside the time window zero
            signal(time_array < time_window(1) |...
                time_array > time_window(2)) = 0;
        else
            
            time_window = [];
        end
        
        % get the signal correpsonding to the object-in-water data after 
        % nulling the scattered waves
        data_object_cell{ind_emitter}(ind_receiver, :) = signal;
        
        % get the signal correpsonding to the object-in-water data after 
        % nulling the scattered waves
        data_water_cell{ind_emitter} = data_water(:, :, ind_emitter);
    end
end


%% ========================================================================
% TRACE THE RAYS AND COMPUTE THE PARAMETERS ON RAYS
%==========================================================================
% get the index of the chosen freequency
[~, frequency_index] = min(abs(1/(2*pi) * omega - para.frequency_single));

    % define the optional inputs for ray tracing
    % For the last optionl input, the acoustic properties should be smoothed
    % once, either smoothing is applied on the inputs of the function (as explained above),
    % or smoothing is applied inside the function.
    ray_args = {'nworker_pool',  para.num_worker_pool, 'interp_method', para.gridtoray_interp ,...
        'raylinking_method', raylinking_method, 'raytogrid_spacing', para.raytogrid_spacing,...
        'auxiliary_ray', trace_auxiliary_ray, 'auxiliary_method', para.auxiliary_method,...
        'reference_angle', para.angle_perturb_auxiliary, 'max_num_points_factor', 2,...
        'smoothing_window_size', para.smoothing_window_size,...
        'max_iter', para.max_raylinking_iter, 'varepsilon', raylinking_threshold,...
        'angular_frequency_centre', omega(frequency_index), 'absorption_power',...
        simulation_prop.alpha_power};
    
   
    
    % compute the trajectory of the rays and accumulated parameters along
    % the rays for the only-water
    [~, ~, ~, ray_position_water, ray_time_water, ~, ~,...
        ray_position_left_water, ray_position_right_water, ~, ~, ~] = ...
        computeRaysParameters(recon_grid, ones(recon_grid.size),...
        zeros(recon_grid.size), emitter.positions, receiver.positions,...
        initial_angles, mask_raytracing, [], 0, ray_args{:});
    
    % the acoustic absorption in water is set scalar zero.
    ray_absorption_water = 0;
    
    
    switch dim
        
        case 2
           
    % compute the trajectory of the rays and accumulated parameters along
    % the rays for the object-in-water
    [~, ~, ~, ray_position, ray_time, ray_absorption,...
        rayspacing_receiver, ray_position_left, ray_position_right,...
        adjoint_ray_position_left, adjoint_ray_position_right, ~] =...
        computeRaysParameters(recon_grid, simulation_prop.sound_speed_ref./sound_speed_phantom,...
        absorption_phantom, emitter.positions, receiver.positions, initial_angles,...
        mask_raytracing, [], 0, ray_args{:});
    
        case 3
            
            % get the number of levels for ray tracing in the 3D case
            num_level = 4;
            
            for ind_level = 1 : num_level + 1
                
                % get the sound speed for the current level
                sound_speed_level = simulation_prop.sound_speed_ref...
                    + (ind_level-1)/num_level * (sound_speed_phantom -...
                    simulation_prop.sound_speed_ref);
                
                
                % compute the trajectory of the rays and accumulated parameters along
                % the rays for the object-in-water
                [~, ~, initial_angles, ray_position, ray_time, ray_absorption,...
                    rayspacing_receiver, ray_position_left, ray_position_right,...
                    adjoint_ray_position_left, adjoint_ray_position_right, ~] =...
                    computeRaysParameters(recon_grid, simulation_prop.sound_speed_ref./sound_speed_level,...
                    absorption_phantom, emitter.positions, receiver.positions, initial_angles,...
                    mask_raytracing, [], 0, ray_args{:});
                
            end
           
    end
    
   
% get the number of grid points inside the mask
% num_gridpoints = nnz(mask_raytracing);

% get the Cartesian position of the grid points inside the binary mask

% x coordinate
grid_x = recon_grid.x(mask_raytracing);

% y coordinate
grid_y = recon_grid.y(mask_raytracing);

% x-y coordinates
grid_pos = [grid_x, grid_y];

% add a column for the z-coordinate 
if dim == 3
    grid_pos = [grid_pos, recon_grid.z(mask_raytracing)];
    distance_emitter_receivers = 1;
else
    distance_emitter_receivers = [];
end


%% ========================================================================
% DEFINE THE IMPLICIT FUNCTIONS FOR APPROXIMATING GREEN'S FUNCTION
% =========================================================================
% Define a handle function for approximating the pressure field using the
% Green's function
approximate_pressure = @(pressure_source, time_delays, geom_attenuation, absorption,...
    caustic_number, source_mode) approxPressureGreens(pressure_source, time_delays,...
    geom_attenuation, absorption, caustic_number, omega, simulation_prop.alpha_power,...
    distance_emitter_receivers, source_mode);

% Define a handle function for approximating the parameters of the Green's
% function on the grid points and receivers
calc_parameters = @(ray_position, ray_time, ray_absorption,...
    ray_position_left, ray_position_right) calcParametersGreens(....
    ray_position, ray_time, ray_absorption, ray_position_left,...
    ray_position_right, sound_speed_phantom/simulation_prop.sound_speed_ref,...
    grid_pos, ray_spacing, detec_radius, para.mask_coeff,...
    'forward', false, interp_receiver);

% Define a handle function for approximating the pressure field on the grid points and
% the grid points using the Green's function
calc_pressure = @(approx_pressure, parameters_grid, caustic_number, nan_grid_binary,...
    parameters_receiver, caustic_receiver, receiver_order, interp_receiver)...
    calcPressureGreensTest(approx_pressure, parameters_grid, caustic_number, nan_grid_binary,...
    pressure_source, 'normal', parameters_receiver, caustic_receiver, receiver_order,...
    length(omega), interp_receiver);



% convert scalar zero to cell containing zeros, if the acoustic absorption is
% not accounted for
if ~iscell(ray_absorption)
    ray_absorption = cell(1, num_emitter);
    ray_absorption(:) = {0};
end

% allocate zero matrices for the dicrepancy between the pressure fieled
% simulated by k-Wave and approximated using the Green's function for
% heterogeneous media.

% matrix including the absorption
relative_discrepancy_absorbing = cell(1, num_emitter);

% matrix excluding the absorption
relative_discrepancy_nonabsorbing = cell(1, num_emitter);

% allocate cell arrays for the measured signals (here simulated using k-Wave)
signal_measured = cell(1, num_emitter);
signal_measured_water = cell(1, num_emitter);

% allocate cell arrays for the signals approximated using the Green's solver
signal_approximated = cell(1, num_emitter);
signal_approximated_water = cell(1, num_emitter);
signal_approximated_nonabsorbing = cell(num_emitter, 1);

% allocate cell array for the parameters on the grid
parameters_grid = cell(1, num_emitter);

% allocate cell arrays for the binary masks
nan_grid_binary = cell(1, num_emitter);

% allocate cell arrays for the time delays on the receivers with repect to each emitter for
% the object inside water, and only water
ray_time_receiver = cell(num_emitter, 1);
ray_time_receiver_water = cell(num_emitter, 1);


parfor (ind_emitter = 1: num_emitter, para.num_worker_pool)
   %  for ind_emitter = 1: num_emitter
    
    % get the object-inside-water measured times series for the current emitter
    signal_measured_emitter = data_object_cell{ind_emitter};
    
    % get the only-water measured times series for the current emitter
    signal_measured_water_emitter = data_water_cell{ind_emitter};
 
    % get the parameters of the Green's function for the object inside
    % water on the grid points and receivers
    [parameters_grid{ind_emitter}, parameters_receiver, nan_grid_binary{ind_emitter}, caustic_number,...
        receiver_order, caustic_receiver, ~, interp_receiver_emitter] = calc_parameters(...
        ray_position{ind_emitter}, ray_time{ind_emitter}, ray_absorption{ind_emitter},...
        ray_position_left{ind_emitter}, ray_position_right{ind_emitter});
    
    % approximate the pressure field for the absorbing object inside water on the grid points and receivers
    [~, signal_approximated{ind_emitter}, ray_time_receiver{ind_emitter}] = calc_pressure(...
        approximate_pressure, parameters_grid{ind_emitter}, caustic_number, nan_grid_binary{ind_emitter},...
        parameters_receiver, caustic_receiver, receiver_order, interp_receiver_emitter);
    
    % get the parameters of the Green's function for only water on the grid points and
    % receivers
    [parameters_grid_water, parameters_receiver_water, nan_grid_binary_water, ~,...
        receiver_order_water, ~, ~, ~] = calc_parameters(...
        ray_position_water{ind_emitter}, ray_time_water{ind_emitter}, ray_absorption_water,...
        ray_position_left_water{ind_emitter}, ray_position_right_water{ind_emitter});
     
    % approximate the pressure field for only water on the grid points and receivers
    [~, signal_approximated_water{ind_emitter}, ray_time_receiver_water{ind_emitter}] = calc_pressure(...
        approximate_pressure, parameters_grid_water, 0, nan_grid_binary_water,...
        parameters_receiver_water, 0, receiver_order_water, interp_receiver_emitter);
    
    % discard the column corresponding to acoustic absorption from the matrix for the parameters
    % of the Green's function, and compute the pressure field by neglecting the acoustic
    % absorption as benchmark for evaluating the accuracy of including acoustic absorption
    % in computing the pressure field using the Green's function.
    % discard the column corresponding to the acoustic absorption from parameters
    if ~isscalar(ray_absorption{ind_emitter})
        parameters_receiver(:, end) = [];
        parameters_grid{ind_emitter} = parameters_grid{ind_emitter}(:, [1, 2, 4]);
    end
    
    % approximate the pressure field for the object inside water
    % on the grid points and receivers, when the absorption effects are
    % ignored.
    [~, signal_approximated_nonabsorbing{ind_emitter}, ~] = calc_pressure(...
        approximate_pressure, parameters_grid{ind_emitter}, caustic_number, nan_grid_binary{ind_emitter},...
        parameters_receiver, caustic_receiver, receiver_order, interp_receiver_emitter);
    

    if para.deconvolve_source
        
        % deconvolve the object-inside-water measured time traces from the pressure
        % source in the frequency domain
        signal_measured_emitter = deconvolve(emitter.pulse,...
            signal_measured_emitter, time_array, omega, deconvolution_parameter);
        
        % deconvolve the only-water measured time traces from the pressure
        % source in the frequency domain
        signal_measured_water_emitter = deconvolve(emitter.pulse,...
            signal_measured_water_emitter, time_array, omega, deconvolution_parameter);
        
    else
        
        % Compute the Discrete Fourier Transform (DFT) of the
        % the object-inside-water measured signals (actually simlulated by k-Wave)
        signal_measured_emitter = discreteFourierTransform(...
            signal_measured_emitter, omega, time_array);
        
        % Compute the Discrete Fourier Transform (DFT) of the
        % the only-water measured signal (actually simlulated by k-wave)
        signal_measured_water_emitter = discreteFourierTransform(...
            signal_measured_water_emitter, omega, time_array);
    end
    
    % get the approximated (Green's) and simulated (k-Wave) time traces
    % recorded on all the receivers for the pressure field propagated after
    % being produced by the current emitter
    signal_measured{ind_emitter} = signal_measured_emitter;
    signal_measured_water{ind_emitter}= signal_measured_water_emitter;
    
    
    % calculate the percentage relative discrepancy of the signals approximeted
    % using Green's function and signals simulated using k-Wave for the object inside water,
    % when the absorption effects are included.    
    relative_discrepancy_absorbing{ind_emitter} = ...
        vecnorm(signal_approximated{ind_emitter} - signal_measured{ind_emitter}, 2, 2) ./...
        vecnorm(signal_measured{ind_emitter}, 2, 2) * 100;
    
    % calculate the percentage relative discrepancy of the signals approximeted
    % using Green's function and signals simulated using k-Wave for the object inside water,
    % when the absorption effects are ignored. 
    relative_discrepancy_nonabsorbing{ind_emitter} =...
        vecnorm(signal_approximated_nonabsorbing{ind_emitter} - signal_measured{ind_emitter}, 2, 2) ./...
        vecnorm(signal_measured{ind_emitter}, 2, 2) * 100;
    
    % display the number of emitter under run
    disp (['The number of emitter is:' num2str(ind_emitter)]);
    
end


% convert the cells for the relative discrepancies to matrix
relative_discrepancy_absorbing = cell2mat(relative_discrepancy_absorbing);
relative_discrepancy_nonabsorbing = cell2mat(relative_discrepancy_nonabsorbing);

%% ========================================================================
% PLOT THE SIMULATED/APPROXIMATED SIGNALS FOR THE CHOSEN EMITTER-RECEIVER
% PAIRS
%==========================================================================
% get the maximum absolute of the pressure source
source_max = max(abs(pressure_source_original));

% plot the amplitude of the pressure field simulated by k-Wave and approximated
% by the Green's function only for the chosen emitter-receiver pairs
h5 = figure;
% Green's function for object inside water and including absorption
plot(1/(2*pi) * omega, 1/source_max * abs(signal_approximated{emitter_index}(receiver_index, :)), 'r');hold on;
% Green's for object inside water and ignoring absorption
plot(1/(2*pi) * omega, 1/source_max * abs(signal_approximated_nonabsorbing{emitter_index}(receiver_index, :)), 'c');hold on;
% Green's for only water
plot(1/(2*pi) * omega, 1/source_max * abs(signal_approximated_water{emitter_index}(receiver_index, :)), 'g-.');hold on;
% k-Wave for object inside water and including absorption
plot(1/(2*pi) * omega, 1/source_max * abs(signal_measured{emitter_index}(receiver_index, :)), 'b--');
xlabel('Frequency [Hz]'); ylabel('Amplitude [Pa s]');
legend('Greens, absorbing breast',...
    'Greens, nonabsorbing breast',...
    'Greens, water',...
    'k-Wave, absorbing breast');

% get the phase for the object-in-water approximated signal for the current emitter-receiver pair
phase_approximated_receiver = wrapToPi(angle(signal_approximated{emitter_index}(receiver_index, :)))/pi;

% get the phase for the the only-water approximated signal for the current emitter-receiver pair
phase_approximated_water_receiver = wrapToPi(angle(signal_approximated_water{emitter_index}(receiver_index, :)))/pi;

% get the phase for the object-in-water measured signal for the current emitter-receiver pair
phase_measured_receiver = wrapToPi(angle(signal_measured{emitter_index}(receiver_index, :)))/pi;


% get the receiver indices on which the discrepancy in phase is
% greater than 1 pi rad. Those discrepncies are because of the jumps at -pi and +pi.
% For these prticular receiver indices, we shift the phase of the measured signals 
% (simulated by k-Wave) by -2pi or +2pi, with signs opposite to the sign of the
% phase on th receivers.
discontinuity_indices = find(abs(phase_approximated_receiver - phase_measured_receiver)> 1);

% shift the phase of the measured signals (simulated by k-Wave) by -2pi or +2pi
phase_measured_receiver(discontinuity_indices) =...
   phase_measured_receiver(discontinuity_indices)...
    - 2 * sign(phase_measured_receiver(discontinuity_indices));

% plot the phase of the pressure field simulated by k-Wand and approximated
% by the Green's function only for the chosen emitter-receiver pair
h6 = figure;
% Green's for the object inside water and including absorption
plot(1/(2*pi) * omega, phase_approximated_receiver, 'r');hold on;
% Green's for only water
plot(1/(2*pi) * omega, phase_approximated_water_receiver, 'g-.');hold on;
% k-Wave for the object inside water and including absorption
plot(1/(2*pi) * omega, phase_measured_receiver, 'b--');
xlabel('Frequency [Hz]'); ylabel('Phase [\pi rad]');
legend('Greens, absorbing breast',...
    'Greens, water', 'k-Wave, absorbing breast');


if para.save_plots
    
    saveas(h5, [plot_directory, 'signal_breast_amplitude'...
        para.shot_angles    '.fig']);
    saveas(h5, [plot_directory, 'signal_breast_amplitude'...
        para.shot_angles    '.png']);
    saveas(h5, [plot_directory, 'signal_breast_amplitude'...
        para.shot_angles    '.tiff']);
    saveas(h5, [plot_directory, 'signal_breast_amplitude'...
        para.shot_angles    '.eps'], 'epsc'); hold off;
    
    saveas(h6, [plot_directory, 'signal_breast_phase'...
        para.shot_angles  '.fig']);
    saveas(h6, [plot_directory, 'signal_breast_phase'...
        para.shot_angles  '.png']);
    saveas(h6, [plot_directory, 'signal_breast_phase'...
        para.shot_angles  '.tiff']);
    saveas(h6, [plot_directory, 'signal_breast_phase'...
        para.shot_angles  '.eps'], 'epsc'); hold off;
    
end


%% ========================================================================
% PLOT THE AMPLITUDES FOR THE CHOSEN EMITTER
%==========================================================================
% Plot the amplitude of the pressure recorded on the receivers
% The amplitudes can be displayed as the maximum of the pressure signals ('maximum'),
% or a chosen single frequency ('single-frequency')
switch para.amplitude_plot_approach
    case 'maximum'
        
        % get the maximum absolute amplitudes
        maximum_approximated = max(abs(signal_approximated{emitter_index}), [], 2);
        maximum_approximated_water = max(abs(signal_approximated_water{emitter_index}), [], 2);
        maximum_approximated_nonabsorbing = max(abs(signal_approximated_nonabsorbing{emitter_index}), [], 2);
        maximum_measured = max(abs(signal_measured{emitter_index}), [], 2);
        
    case 'single-frequency'
        
        % get the amplitudes at the chosen single frequency.
        maximum_approximated = abs(signal_approximated{emitter_index}(:, frequency_index));
        maximum_approximated_water = abs(signal_approximated_water{emitter_index}(:, frequency_index));
        maximum_approximated_nonabsorbing = abs(signal_approximated_nonabsorbing{emitter_index}(:, frequency_index));
        maximum_measured = abs(signal_measured{emitter_index}(:, frequency_index));
        
end


h7 = figure;
% Green's for the object inside water and including absorption
plot(1:num_receiver, 1/source_max * maximum_approximated, 'r');hold on;
% Green's for the object inside water and ignoring absorption
plot(1:num_receiver, 1/source_max * maximum_approximated_nonabsorbing, 'c');hold on;
% Green's for only water
plot(1:num_receiver, 1/source_max * maximum_approximated_water, 'g-.');hold on;
% k-Wave for the object inside water and including absorption
plot(1:num_receiver, 1/source_max * maximum_measured, 'b--');
xlabel('receiver index'); ylabel('Amplitude [Pa s]');
legend('Greens, absorbing breast',...
    'Greens, nonabsorbing breast',...
    'Greens, water ',...
    'k-Wave, absorbing breast'); hold off;
switch para.shot_angles
    case 'raylinking' 
        if para.deconvolve_source
            axis([0 256 0 0.8]);
        else
            switch dim
                case 2
            axis([0 256 0 0.02]);
                case 3
            end
        end
    case 'even-spaced-angle'
        axis([0 256 0 0.05]);
end
set(gca, 'FontSize', 12);


if para.save_plots

    switch para.shot_angles

        case 'raylinking'

            saveas(h7, [plot_directory, 'Fig4b' '.fig' ]); % (or Fig4d.fig')
            saveas(h7, [plot_directory, 'Fig4b' '.png' ]);
            saveas(h7, [plot_directory, 'Fig4b' '.tiff']);
            saveas(h7, [plot_directory, 'Fig4b' '.eps' ], 'epsc');hold off;

        case 'even-spaced-angle'

            saveas(h7, [plot_directory, 'Fig4f' '.fig' ]);
            saveas(h7, [plot_directory, 'Fig4f' '.png' ]);
            saveas(h7, [plot_directory, 'Fig4f' '.tiff']);
            saveas(h7, [plot_directory, 'Fig4f' '.eps' ], 'epsc');hold off;

    end
end

%% ========================================================================
% PLOT THE PHASE FOR THE CHOSEN EMITTER
%==========================================================================
% Calculate the phases in terms of [pi rad]
% Green's for the object inside water and including absorption
phase_approximated_emitter = wrapToPi(angle(signal_approximated{ind_emitter}(:, frequency_index)))/pi;
% Green's for only water
phase_approximated_water_emitter = wrapToPi(angle(signal_approximated_water{ind_emitter}(:, frequency_index)))/pi;
% k-Wave for the object inside water and including absorption
phase_measured_emitter = wrapToPi(angle(signal_measured{ind_emitter}(:, frequency_index)))/pi;

% get the receiver indices on which the discrepancy in phase [pi rad] is
% greater than 1. Those discrepncies are because of the jumps at -pi and +pi.
% For these prticular receiver indices, we shift the phase of the measured signals 
% (simulated by k-Wave) by -2pi or +2pi, with sign opposite to the sign of the
% phase.
discontinuity_indices = find(abs(phase_approximated_emitter - phase_measured_emitter)> 1);

% shift the phase of the measured signals (simulated by k-Wave) by -2pi or +2pi
phase_measured_emitter(discontinuity_indices) =...
    phase_measured_emitter(discontinuity_indices)...
    - 2 * sign(phase_measured_emitter(discontinuity_indices));



h8 = figure;
% Green's for the object inside water and including absorption
plot(phase_approximated_emitter, 'r'); hold on;
% Green's for only water
plot(phase_approximated_water_emitter, 'g-.');hold on;
% Green's for the object inside water and including absorption
plot(phase_measured_emitter, 'b--');

xlabel('receiver index'); ylabel('Phase [\pi rad]');
legend({'Greens, absorbing breast',...
    'Greens, water', 'k-Wave, absorbing breast',...
    }, 'Location', 'northeast'); hold off;

 axis([0 256 -1.1  +1.1]); 
set(gca, 'FontSize', 12);


if para.save_plots
    switch para.shot_angles
        case 'raylinking'

            saveas(h8, [plot_directory, 'Fig4a'  '.fig' ]); % (or Fig4c.fig)
            saveas(h8, [plot_directory, 'Fig4a'  '.png' ]);
            saveas(h8, [plot_directory, 'Fig4a'  '.tiff']);
            saveas(h8, [plot_directory, 'Fig4a'  '.eps' ], 'epsc');hold off;

        case 'even-spaced-angle'

            saveas(h8, [plot_directory, 'Fig4e' '.fig' ]);
            saveas(h8, [plot_directory, 'Fig4e' '.png' ]);
            saveas(h8, [plot_directory, 'Fig4e' '.tiff']);
            saveas(h8, [plot_directory, 'Fig4e' '.eps' ], 'epsc');hold off;
    end
end


end
    