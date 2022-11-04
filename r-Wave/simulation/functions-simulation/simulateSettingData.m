function [data, data_ref, data_noisy, data_ref_noisy, emitter, receiver,...
    kgrid, medium, simulation_prop, data_simulation_time] =...
    simulateSettingData(do_data_sim, dim, scenario,...
    simulation_purpose, excit_pulse_name, half_grid_size, num_receiver,...
    ratio_receiver_to_emitter, do_absorption, oa_breast_path, machine_name,...
    result_path, local_result_path, varargin)
%SIMULATESETTINGDATA simulates an ultrasound setting and ultrasound data
%   
%
% DESCRIPTION:
%             simulateSettingData simultes an experimental setting for
%             ultrasound imaging and simulates the propagation of the
%             acoustic waves using k-Wave
%
% USAGE:
%     
%
% INPUTS:
%       do_data_sim       - Boolean controlling whether the ultrasound data
%                           are simulated (true), or already stored data
%                           corresponding to the inputs are loaded (false). 
%       dim               - the dimensions of the k-Wave simulation
%       scenario          - the scenario, in which the simulation is done
%                           for all emitter-receiver pairs ('standard'), or simulation
%                           is only done for a single emitter ('single_emitter')
%       simulation_purpose - the purpose of simulation of data. This can be 
%                            'raytracing_validation', or 'image_reconstruction'
%                            (Default: image reconstruction)
%       excit_pulse_name  - the name of the available excitation pulse. The
%                           pulse is 'pammoth'       
%       half_grid_size    - the half-size [m] of side of the computational
%                           grid
%       num_receiver      - the number of receiver
%       ratio_receiver_to_emitter - the ratio of the number of receivers to the
%                           number of emitters
%       do_absorption     - Boolean controlling whether absorption (and dispersion)
%                           are included in the k-Wave simulation or not.
%       oa_breast_path    - the path for the digital phantom
%       machine_name      - the name of the machine used for simulations
%       result_path       - the path for saving the results
%       local_result_path - the local path for saving the results 
%  
%      
 
% OPTIONAL INPUTS:
%      'num_worker_pool'  - the number of workers for parallel programming
%                           (Default: 8)
%      'num_worker_hdf5'  - the number of workers for creating HDF5 files
%                           (Default: 8)
%      'data_save_format' - the format for storing the ultrasound time
%                           series. This can be 'mat' or 'hdf5'.
%                           (default:'mat')
%      'do_calculate_tofs'- Boolean controlling whether the time-of-flights
%                           are computed from simulated (or loaded) data,
%                           or the already computed tofs are used.
%      'cfl'              - the cfl number, which adjusts the temporal
%                           sampling to the grid spacing
%      'grid_spacing_data'- the grid spacing [m] for data simulation 
%                           (Default: 1e-3)
%      'grid_spacing_reconstruction' - the grid spacing for ray tracing and
%                           image reconstruction
%      'z_pos_height'     - the maximum z position [m] of the grid points, if
%                           the reconstruction geometry is set 'real'
%                           (Default : 1e-2)
%      'transducer_geom'  - the geometry of the transducers. This can be
%                           'point' (2D, 3D), 'line' (2D), or 'disc' (3D)
%                           (Default: 'point')
%      'transducer_interp_approach' - the approach for interpolation from
%                            the grid points onto the transducers and vice versa
%                            This can be This can be 'nearest' or 'offgrid'. 
%                            (Default = 'offgrid')
%      'transducer_radius' - the radius [m] of the transducer
%      'do_low_filter'     - Boolean controlling whether a low-pas filter 
%                            is applied to the excitation pulse, or not
%      'low_filter_coeff'  - a coefficient multiplied by the minimum speed
%                            of sound for adjusting the low pass filter.
%                            The minimum sound speed of the medium may
%                            not sufficiently cut off the high
%                            frequencies, and these frequencies may lead to
%                            some oscilations in the simulated signals
%     'normalise_excitation_amplitude' - Boolean controlling whether the
%                            amplitude of the excitation pulse is normalised
%                            or not.
%     'time_downsampling_rate' - the downsampling rate applied to the second
%                            (time) dimension of the data simulated by k-Wave
%                           (Default: 2)
%     'emitter_downsampling_rate' - the downsampling rate applied to the third
%                            (emitter) dimension of the data simulated by k-Wave (Default: 1)
%     'receiver_downsampling_rate' - the downsampling rate applied to the first
%                            (receiver) dimension of the data simulated by k-Wave (Default: 1)
%     'noise_level'         - the SNR of each signal, (Default: 40)
%     'noise_mode'          - the mode for adding noise. This can be 'peak',
%                            or 'snr'. (default:'peak').
%     'code_type'          - the code type for the k-Wave simulation. This
%                            can be 'Matlab', or 'CUDA'
%     'save_to_disk'       - control what information is saved. This can
%                          - 0, 1, or 2
%                            0: For each excitation, an HDF5 file
%                            containing the inputs for each k-Wave
%                            simulations is saved, and the pressure
%                            data is also simulated and saved
%                            1: only the pressure data is simulated
%                            and saved, 2: only the HDF5 files containing
%                            the inputs for the k-Wave simulation is saved
%      'data_cast'          - the precision for data storage. This can be
%                             'single' or 'double' (Default: 'double' (2D),
%                             'single' (3D))
%      'gpu_index'          - the index of the used gpu. This can be 0,...,
%                             n-1, where n is the number of the available
%                             gpus
%      'single_emitter_receiver' - the choice of single-emitter pairs when
%                              the scenario for simulations is set
%                              'single_emitter'. If the first component is
%                              not infinite, the simulation is done for a
%                              single emitter and for all the receivers.
%      'plot_sim'           - Bolean controlling whether the plots, e.g.,
%                             exctatoon pulse in time or frequnecy
%                             domain,... are simulated, or not. 
%                             (Default : true)
%      
%      
%     
%
% OUTPUTS:
%       data                - the num_receiver x num_time x num_emitter
%                             ultrasound data for the object in water with
%                             num_receiver the number of receivers, num_time
%                             the length of the time, and num_emitter the
%                             number of emitters
%       data_ref            - the num_receiver x num_time x num_emitter
%                             ultrasound data for the only water with
%                             num_receiver the number of receivers, num_time
%                             the length of the time, and num_emitter the
%                             number of emitters
%       data_noisy          - the noisy version of object-in-water data
%       data_ref_noisy      - the noisy version of only-water data
%       emitter           - a struct which defines the properties of the
%                           excitation as follows: This includes the
%                           fields 'positions', 'pulse', 'pulse_duration',
%                           and 'shot_time'
%       'positions'       - 2/3 x num_emitter array of Cartesian position 
%                           of the center of the emitter objects
%       'pulse'           - a vector of size 1 x num_te with num_te <= num_time.
%       'shot_time'       - the first arrival of the excitation pulse, which
%                           can be negative [sec]
%       receiver          - a struct which includes the fields
%       'positions'       - 2/3 x N array of Cartesian position of the
%                           centre of the receiver objects. For rotating
%                           ultrasound systems, the position of the receivers
%                           is changed for different excitations. For the rotating case,
%                           this is a cell array with length num_angle (the number of rotation
%                           angles) with each cell a 2/3 x num_transducer array of
%                           Cartesian position of transducers for each
%                           angle. (num_receiver = num_angle x num_transducer)
%      kgrid               - a struct array including the data for the k-wave simulation 
%      medium              - a struct array for the acoustic maps used for
%                            the k-Wave simulation. This strut includes the
%                            fields:
%      'sound_speed'       - the sound spped distribution [m/s]
%      'alpha_coeff'       - the absorption coefficient [dBMHz^{-y}
%                            cm^{-1}]
%      'alpha_power'       - the exponent power of the attenuation used for
%                            the k-Wave simulation 
%      simulation_prop      - a struct array containing the properties of the
%                           simulation. This struct includs the fields: 
%      'PML'                - a  1 x dim vector of the number of PML layers
%                           along each dimension
%      'f_max'              - the maximum frequency [Hz] supported by the
%                           computational grid    
%      'detec_radius'       - the radius of the detection ring (surface)
%      data_simulation_time - the cpu time [s] for simulating ultrasound data
%                             for the object in water
%       
%    
%
% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 10.08.2022
%
% This script is part of the r-Wave Tool-box (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian and Ben Cox



%% ========================================================================
%  SET OPTIONAL PARAMETERS
%  ========================================================================
if any(strcmp(machine_name(1:end-1), {'kinsler','nyborg'}))
    
    % the number of workers for data simulation
    para.num_worker_pool = 16;
    
    % the number of workers for creating HDF5 files
    para.num_worker_hdf5 = 16;
    
else
    
    % give a message
    disp(['The number of workers are set' num2str(4) ', but the user could'...
        'adjust the number of workers for data simulation and image reconstruction.'])
    
    % the number of workers for data simulation
    para.num_worker_pool = 4;
    
    % the number of workers for creating HDF5 files
    para.num_worker_hdf5 = 4;
    
end


% save the pressure times series in the mat format
para.data_save_format = 'mat'; % or 'hdf5'

% check the number of dimenstions is set 2 or 3.
if ~ismember(dim, [2,3])
    error('dim must be either 2, or 3.')
end
    
switch dim
    case 2
        
        % the time-of-flights are computed from simulated or loaded data,
        % or already computed time-of-flights are used.
        para.do_calculate_tofs  = true;
        
        % the height [m] of the chest wall (volume of z > 0), not applied to the 2D case
        % the maximum
        para.z_pos_height = 0;
        
        % the grid spacing [m] for data simulation
        para.grid_spacing_data = 4000e-7;
        
        % the code type for data simulation
        para.code_type = 'Matlab';
        
        % the precision for data storage
        para.data_cast = 'double';
        
    case 3
        
        % the time-of-flights are computed from simulated or loaded data,
        % or already computed time-of-flights are used.
        para.do_calculate_tofs  = true;
        
        % the height [m] of the chest wall (volume of z > 0)
        para.z_pos_height = 5e-3;
        
        % the grid spacing [m] for data simulation
        para.grid_spacing_data = 5000e-7;
        
        % the code type for data simulation
        para.code_type = 'CUDA';
        
        % the precision for data storage
        para.data_cast = 'single';
end

% the cfl number, which adjusts the temporal sampling to the grid spacing
para.cfl = 0.1;

% the grid spacing [m] for image reconstruction
para.grid_spacing_reconstruction = 1e-3;

% the geometry of the transducers
para.transducer_geom = 'point';

% the appproch for interpolating the pressure from grid points to the
% transducers and vice versa
para.transducer_interp_approach = 'offgrid';

% Boolean controlling whether a low-pass filter is applied on the excitation
% pulse for removing the high frequencies which are not supported by the 
% computational grid or not.
para.do_low_filter = true;

% a coefficient multiplied by the minimum speed of sound for adjusting the
% low pass filter applied on the excitation pulse. This low pass filter is 
% applied on the excitation pulse to ensure that the maximum frequency of the
% excitation pulse be less than the maximum frequency supported by the
% computational grid
para.low_filter_coeff = 1;

% Boolean controlling whether the amplitude of the excitation pulse is
% normalised or not.
para.normalise_excitation_amplitude = true;

% the mode for adding noise
para.noise_mode = 'peak';

% the noise level
para.noise_level = 40;

% the downsampling rate of data along the time dimension. Different temporal sampling rate
% is used for the data simulated by  the k-Wave simulation and the data used for image 
% reconstruction in order to avoid an inverse crime in time sampling.
para.time_downsampling_rate = 2;  % (default : 2)

% the downsampling rate of emitters
para.emitter_downsampling_rate = 1; 

% the downsampling rate of receivers
para.receiver_downsampling_rate = 1; 


% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


if strcmp(para.transducer_geom, 'point')
    
    % the radius of the transducers
    para.transducer_radius = [];
else
    % the radius [m]  of the transducers, i.e., the half-length [m] of the line (for 2D) or
    % radius [m] of the disc (for 3D)
    para.transducer_radius = 1e-3;
end


% choose the smoothing window size based on the grid spacing [m]
if para.grid_spacing_reconstruction < 1e-3
    error(['For the sizes used for breast imaging, and using ray-based methods,'...
        'the grid spacing smaller than 1mm is not necessary.'])
elseif para.grid_spacing_reconstruction < 1.5e-3 
    para.smoothing_window_size = 7;
elseif para.grid_spacing_reconstruction < 2.01e-3  
    para.smoothing_window_size = 5;
else
    error('The grid spacing must not be larger than 2mm.')
end

    

% spacing of the k-Wave simulation
if do_data_sim
    
    % control what information is saved. This can
    % 0, 1, or 2
    % 0: For each excitation, an HDF5 file
    % containing the inputs for each k-Wave
    % simulations is saved, and the pressure
    % data is also simulated and saved
    % 1: only the pressure data is simulated
    % and saved
    % 2: only the HDF5 files containing the inputs for
    % the k-Wave simulation is saved
    para.save_to_disk = 1;
end

% the index of the used GPU
para.gpu_index = 4; 

% display the index of the used GPU
disp(['The index of the used gpu:' num2str(para.gpu_index)]);

    
switch scenario
    case 'standard'
        
        % allocate an empty variable
        para.single_emitter_receiver = [];
        
    case 'single_emitter'
        
        % get the indices of a single pair of receiver-emitter pairs used for
        % validating ray tracing as a 1 x 2 vector, where the first and
        % second components correpond to the index of receiver and emitter,
        % respectively. if the first component is set nan, all receivers
        % are chosen.
        
        switch dim
            case 2
        para.single_emitter_receiver = [nan, 1];
            case 3
        para.single_emitter_receiver = [nan, 1000]; 
        end
        
        if strcmp(simulation_purpose, 'image_reconstruction')
            error('A single emitter is not sufficient for image reconstruction.')
        end
        
end


% Boolean controlling whether the plots, e.g., exctitation pulse in time or frequency
% domain,... are simulated, or not.
para.plot_sim = true;


% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% check the sampling rates
if   para.time_downsampling_rate < 1  || rem(log2(para.time_downsampling_rate), 1) ~= 0
    error(['The time sampling rate must be a power of 2',...
        'with power an integer number greater than zero']);
end
if   para.emitter_downsampling_rate < 1  || rem(log2(para.emitter_downsampling_rate), 1) ~= 0
    error(['The emitter sampling rate must be a power of 2',...
        'with power an integer number greater than zero']);
end
if   para.receiver_downsampling_rate < 1  || rem(log2(para.receiver_downsampling_rate), 1) ~= 0
    error(['The receiver sampling rate must be a power of 2',...
        'with power an integer number greater than zero']);
end


if rem(para.smoothing_window_size, 2) == 0 && para.smoothing_window_size > 3
    error('The smoothing window size must not be set an even value, if it is larger than 3.')
end

if strcmp(scenario, 'single_emitter')
    
    % emiters are not downsmapled if the simulation is performed for a
    % single emitter.
    para.emitter_downsampling_rate = 1;
end

%% ========================================================================
% SET FIXED INPUT PARAMETERS
%==========================================================================
% get the grid spacings for the grid for data simulation
grid_spacing_data = para.grid_spacing_data;

% get the grid spacings for the grid for image reconstruction
grid_spacing_reconstruction = para.grid_spacing_reconstruction;

% display the grid spacing for k-Wave simulation, and image
% reconstruction(or validation of ray tracing)
disp(['The smoothing window size for image reconstruction (or validation of ray tracing) is:'...
    num2str(para.smoothing_window_size)])
disp(['The grid spacing for k-Wave simulation is:' num2str(grid_spacing_data, '%2.2e') 'm'])
disp(['The grid spacing for image reconstruction or validation of ray tracing is:'...
    num2str(grid_spacing_reconstruction, '%2.2e') 'm'])

% the mode for image reconstruction. This can be 'transmission' or 'reflection'
% At the moment, only the 'transmission' mode is included.
ultrasound_mode = 'transmission';

% choose the phantom. All experiments are tested with this phantom.
phantom = 'breast';

% Boolean controlling whether the phantom is created and saved, or loaded.
% If oa_breast_path is not given, it has been automatically set nan in startup.m 
if isnan(oa_breast_path)
    do_create_phantom = false;
else
    do_create_phantom = true;
end

    
% set the Boolean controlling whether the pressure data are simulated for
% the entire computational grid or not. The pressure time series are not required
% for the entire grid, because that will dramatically increase the
% computational memory,
% so this varibale is always set false.
simulate_on_grid = false;

% get the order of the gpu
% The gpu numbers are started from zero, but gpu orders are started from 1.
gpu_order = para.gpu_index + 1;

% Boolean controlling whether the initial pressure distribution is smoothed
% in space. This is always set false.
smooth_source = false;

% Boolean controlling whether the sound speed distribution is smoothed
% in space. This is always set false.
smooth_sound_speed = false;

% check the number of receiver is a power of 2, with power an integer
% number greater than zero.
if  num_receiver < 1  || rem(log2(num_receiver), 1) ~= 0
    error(['The number of receiver must be a power of 2',... 
        'with power an integer number greater than zero']);
end   
    

%% ========================================================================
%  DEFINE THE PHYSICAL AND COMPUTATIONAL SETTINGS FOR DATA SIMULATION
%  ========================================================================
settings_args = {'DataSim', do_data_sim, 'Scenario', scenario, 'Mode', ultrasound_mode, ...
    'Excit', excit_pulse_name, 'InterpType', para.transducer_interp_approach, 'TransGeom', para.transducer_geom, ...
    'singleER', para.single_emitter_receiver, 'Absorption', do_absorption, 'Low_Filter', para.do_low_filter,...
    'Low_Filter_Fac', para.low_filter_coeff, 'CodeVersion', para.code_type, 'CreatePhantom',...
    do_create_phantom, 'plot_sim', para.plot_sim};


[kgrid, medium , emitter, receiver, simulation_prop, data_paths] = ...
    makeDataSettings(dim, para.grid_spacing_data, half_grid_size, para.z_pos_height, num_receiver,...
    ratio_receiver_to_emitter, phantom, para.cfl, oa_breast_path, result_path,...
    local_result_path, settings_args{:});


if para.normalise_excitation_amplitude
    
    % normalise the amplitude of the excitation pulse, if required
    emitter.pulse = 1/max(abs(emitter.pulse)) * emitter.pulse;
end

%% ========================================================================
% SET THE SMOOTHING WINDOW SIZE APPLIED ON THE GRID FOR THE K-WAVE
% SIMULATION, IF THE PURPOSE OF SIMULATION IS VALIDATING RAY TRACING
%==========================================================================

if strcmp(simulation_purpose, 'raytracing_validation')
    
    % For validation of ray tracing, the acoustic properties of the medium
    % for the k-Wave simulation must be smoothed by a smoothing window
    % with size equivalent to the window size which will be used for ray tracing.
    
    % compute the smoothing window size for k-wave simulation
    if rem(grid_spacing_reconstruction, grid_spacing_data) == 0
        
        smoothing_window_size_data = grid_spacing_reconstruction/grid_spacing_data *...
            para.smoothing_window_size;
        
        % make the smoothing window size for simulation an odd number
        if rem(smoothing_window_size_data, 2) == 0
            smoothing_window_size_data = smoothing_window_size_data + 1;
        end
    else
        smoothing_window_size_data = [floor(grid_spacing_reconstruction/grid_spacing_data *...
             para.smoothing_window_size),...
            ceil(grid_spacing_reconstruction/grid_spacing_data * para.smoothing_window_size)];
        
        % remove the even component
        smoothing_window_size_data(rem(smoothing_window_size_data, 2)==0) = [];
    end
    
    % smooth the sound speed of the medium for data simulation
    medium.sound_speed = smoothField(medium.sound_speed,...
        smoothing_window_size_data, medium.sound_speed_ref);
    
    if isfield(medium, 'alpha_coeff')
        
    % smooth the absorption coefficient of the medium
    medium.alpha_coeff = smoothField(medium.alpha_coeff,...
        smoothing_window_size_data, 0);
    
    
    % make the nonnegative absorption coefficients zero
    medium.alpha_coeff(medium.alpha_coeff < 0) = 0;
    
    end
    
    
    % set the exponent power 2 
    % medium.alpha_power = 2;
    
end


% determine the phantom type (smmoth or nonsmooth) from the puropose of
% simulation
switch simulation_purpose
    case 'raytracing_validation'
        phantom_type = ['_smooth_', num2str(smoothing_window_size_data)];
    case 'image_reconstruction'
        phantom_type = '_nonsmooth';
end

if strcmp(scenario, 'single_emitter')
    phantom_type = [phantom_type, '_', num2str(para.single_emitter_receiver(2))];
end

% get the path for saving or loading the simulated pressure time series
data_path = [data_paths.main_directory, data_paths.directory, data_paths.name_data,...
        'data', num2str(1e4 * grid_spacing_data), phantom_type];


if do_data_sim
    
    %% ========================================================================
    %  DEFINE A MASK AND A SPARSE MATRIX FOR INTERPOLATION OF EMITTERS AND
    %  RECEIVERS
    %  ========================================================================
    interp_args = {'Trans_Geom', para.transducer_geom, 'Interp_Method', para.transducer_interp_approach,...
        'radius', para.transducer_radius, 'upsampling_rate', 16, 'simulate_on_grid', simulate_on_grid};
    
    % get the interpolation parameters for emitters
    [interp_params.emitter.mask, interp_params.emitter.mask_entire_volume,...
        interp_params.emitter.matrix , interp_params.emitter.elpased_time] =...
        calcInterpParameters(kgrid, emitter.positions, simulation_prop.z_offset, interp_args{:});
    
    % get the interpolation parameters for receivers
    [interp_params.receiver.mask, interp_params.receiver.mask_entire_volume,...
        interp_params.receiver.matrix, interp_params.receiver.elapsed_time] =...
        calcInterpParameters(kgrid, receiver.positions, simulation_prop.z_offset, interp_args{:});
    
    
    
    %% ====================================================================
    %  SIMULATE (OR LOAD) THE DATA
    %======================================================================
    % get the optional inputs for the k-Wave simulation
    sim_args =  {'savetodisk', para.save_to_disk, 'DataCast', para.data_cast,...
        'code_type', para.code_type, 'num_worker_data', para.num_worker_pool,...
        'num_worker_hdf5', para.num_worker_hdf5, 'smooth_source',...
        smooth_source, 'smooth_sos', smooth_sound_speed,...
        'PMLSize', simulation_prop.PML};
    if strcmp(para.code_type, 'CUDA')
        sim_args = {sim_args{:}, 'gpu_index', para.gpu_index, 'gpu_order', gpu_order};
    end
    
   
    % simulate the first data set for the object inside water
    disp('Running simulation for the object inside water:');
    [data, ~, data_paths, data_simulation_time] = kwaveSimulateDataUST(kgrid, medium, ...
        emitter, receiver, interp_params, data_paths, sim_args{:} );
     

        switch para.data_save_format
            case 'mat'
                
                % save the pressure time series and the hdf5 file in a mat
                % file
                save(data_path, 'data', 'sim_args', '-v7.3');
            case 'hdf5'
                if ~exist(data_path, '.h5')
                    
                    % create a hdf5 file and save the pressure time series
                    % inside that
                    h5create([data_path, '.h5'], '/data', size(data))
                    
                    % save the optinal inputs for the simulation and the
                    % hdf5 file in a mat file
                    save(data_path, '.mat', 'sim_args', '-v7.3');
                end 
             
        end
   
    
    if strcmp(ultrasound_mode, 'transmission')
        
        % define a homogeneous and nonabsorbing medium (water)
        medium_ref.sound_speed_ref = medium.sound_speed_ref;
        medium_ref.sound_speed = medium_ref.sound_speed_ref * ones(size(kgrid.x));
        medium_ref.density = 1000;
        
        disp('Running simulation for water:');
        % clear data medium
        [data_ref, ~, ~, data_simulation_time] = kwaveSimulateDataUST(kgrid, medium_ref, ...
            emitter, receiver, interp_params, data_paths, sim_args{:} );
        
            
            switch para.data_save_format
                case 'mat'
                    
                    % save the pressure time series and the hdf5 file in a mat
                    % file
                    save(data_path, 'data_ref', '-append');
                    
                case 'hdf5'
                    if ~exist(data_path, '.h5')
                        
                        % create a hdf5 file and save the pressure time series
                        % inside that
                        h5create([data_path, '.h5'], '/data_ref', size(data_ref))
                        
                    end
            end
            
     end
    
    
else
    
% load the data for the object inside water
         switch para.data_save_format
             case 'mat'
                 if dim == 2   || strcmp(simulation_purpose, 'raytracing_validation')
                     load([data_path, '.mat'], 'data', 'data_ref');
                 else
                     data = []; data_ref = [];
                 end
             case 'hdf5'
         end
   
    
    % set the time for data simulation an empty variable
     data_simulation_time = [];
    
end

% For dim = 3, transform the z position of the transducers from the centred-origin
% (used for k-Wave simulation) to real
if dim == 3
    
    % transform the z position of receivers
     receiver.positions(3,:) = receiver.positions(3,:) - simulation_prop.z_offset;
     
     % transform the z position of emitters
     emitter.positions(3,:) = emitter.positions(3,:) - simulation_prop.z_offset;  
     
     % remove the field 'z_offset' from the struct for simulation
     % properties
     % simulation_prop = rmfield(simulation_prop, 'z_offset');
     
end

%% ========================================================================
% CORRECT FOR THE MISSED TEMPORAL INTEGRAL 
%==========================================================================
% Boolean controlling whether the original exciation pulse, or its first
% derivative is given as the output. This is not optional parameter for this script,
% and must be always set true. 
derivative_emitter_pulse = true;

if derivative_emitter_pulse
    
    % the coefficient indicating to the time integration of the excitation pulse
    % included in the data simulation
    data_attenuation_factor =  1/2 * kgrid.dt * medium.sound_speed_ref/kgrid.dx;
    
    % calcuate numerically the derivative of the excitation pulse,
    excitation_pulse =  1/ data_attenuation_factor *  diff([emitter.pulse, 0]);
    
    if dim == 3
        excitation_pulse =  kgrid.dx * excitation_pulse;
    end
else
   
    % alternatively the user can integrate the exciation pulse in time before
    % giving it as an input to the k-Wave. For this case,
    % 'derivative_emitter_pulse' must be set false.
    % excitation_pulse = emitter.pulse;
    error('Not implemented yet!')  
    
end


if para.do_calculate_tofs
        
        %% ========================================================================
        % ADD NOISE TO THE TIME-SERIES
        %==========================================================================
        noise_args = {'NoiseType', para.noise_mode, 'SNR', para.noise_level};
        [data_noisy] = addUsctNoise(data, noise_args{:});
        [data_ref_noisy] = addUsctNoise(data_ref, noise_args{:});
        
        %%=================================================================
        % DOWNSMAPLE THE DATA
        %==================================================================
        % downsample in time
        % display the time downsampling rate
        disp(['The downsmapling rate is:' num2str(para.time_downsampling_rate)])
        
        % downsample the time array in time
        kgrid.t_array = kgrid.t_array(1:para.time_downsampling_rate:end);
        
        % downsample the exciation pulse in time, and update the emitter struct
        % based on the updated excitation pulse
        emitter.pulse = excitation_pulse(1:para.time_downsampling_rate:end);
        
        
        % downsample the emitters and the data, if required
        [data, data_ref, data_noisy, data_ref_noisy, emitter] = subsampleData(data,...
            data_ref, data_noisy, data_ref_noisy, emitter,...
            para.emitter_downsampling_rate, 'emitter');
        
        % downsample the receivers and the data, if required
        [data, data_ref, data_noisy, data_ref_noisy, receiver] = subsampleData(data,...
            data_ref, data_noisy, data_ref_noisy, receiver,...
            para.receiver_downsampling_rate, 'receiver');
        
        % downsample the data in time
        [data, data_ref, data_noisy, data_ref_noisy, ~] = subsampleData(data,...
            data_ref, data_noisy, data_ref_noisy, [],...
            para.time_downsampling_rate, 'time');
        
        else
        
        
        %%=================================================================
        % DOWNSMAPLE THE DATA
        %==================================================================
       
        % display the time downsampling rate
        disp(['The downsmapling rate is:' num2str(para.time_downsampling_rate)])
        
        % downsample the time array in time
        kgrid.t_array = kgrid.t_array(1:para.time_downsampling_rate:end);
        
        % downsample the exciation pulse in time, and update the emitter struct
        % based on the updated excitation pulse
        emitter.pulse = excitation_pulse(1:para.time_downsampling_rate:end);
        
        
        % downsample the emitters and the data, if required
        [data, data_ref, ~, ~, emitter] = subsampleData(data, data_ref,...
            [], [], emitter, para.emitter_downsampling_rate, 'emitter');
        
        % downsample the receivers and the data, if required
        [data, data_ref, ~, ~, receiver] = subsampleData(data, data_ref,...
            [], [], receiver, para.receiver_downsampling_rate, 'receiver');
        
         % downsample the data in time
        [data, data_ref, ~, ~, ~] = subsampleData(data, data_ref,...
            [], [], [], para.time_downsampling_rate, 'time');

        % The data simulations took one week on 8 GPUs, and may be
        % impractical for many users. 
        % For the 3D case, only time-of-flight imaging is done using the
        % already computed time-of-flights so that data simulation using
        % k-Wave is not used.
        data_noisy = [];
        data_ref_noisy = [];
end



end