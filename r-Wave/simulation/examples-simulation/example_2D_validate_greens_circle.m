% Example for evaluating accuracy of ray approximation to Green's function
% in approximating phases and amplitudes [1] (cf. section References)
%..........................................................................
% The scenario is based on section 6.2. ❝Numerical validation of the ray 
% approximation to the Green’s function.❞ in [1]
%
% The approach taken for computing the geomterical portion of the amplitudes
% is based on paraxial ray tracing. (cf. Algorithm 1 in [1])

%%=========================================================================
% SIMULATION (OR LOAD) OF UST DATA USING K-WAVE
%==========================================================================
% This example first simulates (or loads) two ultrasound data after excitation
% of a single emitter, one for only water and another for the breast in water.
% The datasets can be either simulated or downlaoded from the Zenodo link
% provided in the referene [3]. 
% The simulated data were recorded on all receivers. The simulations were
% performed on a grid with spacing 0.4 mm using the k-Wave toolbox [5]. For the
% breast-in-water data, the sound speed and absorption coefficient maps are 
% smoothed by an averaging window of size 17 grid points. The absorption and
% dispersion are accounted for based on a frequency-power law. The pressure
% time series simulated on the receivers are then decomposed into phase and
% amplitude using a Fourier transform operator. The details are provided in
% section 6.2. ❝Numerical validation of the ray approximation to the Green’s
% function.❞ in [1]


%%=========================================================================
% GREEN'S FORMULA IN HOMOGENEOUS MEDIA (ONLY WATER)
%==========================================================================
% The homogeneous Green's function will then be used for approximating pressure
% field on all receivers after excitation of a chosen single emitter and
% propagation of the pressure field in only water. The approximated pressure
% field on all receivers are compared to the pressure time series simulated by
% the k-Wave in terms of phase and amplitude at the single frequency 1MHz.

%%=========================================================================
% RAY APPROXIMATION TO GREEN'S FORMULA IN HETEROGENEOUS MEDIA (BREAST-IN-WATER)
%==========================================================================
% The ray approximation to heterogeneous Green's function using ray tracing
% will then be used for approximating pressure field on all receivers after
% excitation of a chosen single emitter and propagation of the pressure field
% across the breast in water. The pressure time series approximated on all
% receivers are compared with the pressure time series simulated by the k-Wave
% after decomposing into phases and amplitudes at the single frequency 1MHz.
% Using the Green's approach, the amplitudes on the receivers are plotted before
% and after incorporating the acoustic absorption relying on the frqeuency power law.
% Using the Green's approach, the dispersion is accounted for in computing the wrapped
% phases on receivers.
% It is highly recommended that papers [1] and [2], especially the details
% given in:
% section 6.2. ❝Numerical validation of the ray approximation to the Green’s
% function❞ are carefully read by the users.

%%=========================================================================
% DIGITAL BREAST PHANTOM AND K-WAVE
%==========================================================================
% The Breast is simulated using a digital breast phantom developed by
% Mark Anastasio group [4]. The phantom data must be downloaed via the link:...
% https://anastasio.bioengineering.illinois.edu/downloadable-content/oa-breast-database/ [3].
% The pressure fied  used as the benchmark is simulated using the k-Wave
% toolbox. www.k-Wave.org (v. 1.3. or 1.4.) [5].
% It was shown in [6] that the step for inclusion of the source in the k-Wave
% must be adapted to a source $s$ so that the simulated transmission synthetic
% data agree with the Green's formula. 
%%=========================================================================
% REFERENCES
%==========================================================================
% If you find the toolbox useful for your research, please consider citing these papers:...
% 1 - A. Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution
% quantitative ultrasound tomography❞, 2022, https://arxiv.org/abs/2211.00316.
% 2 - A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering
% for biomedical ultrasound tomography❞, Inverse Problems vol. 37, no.11, 115003, 2021.

% These data-base/toolbox were used in this project.
% 3 - A. Javaherian, 2023, ❝Transmission ultrasound data simulated using the k-Wave toolbox as a benchmark for biomedical quantitative ultrasound tomography using a ray approximation to Green's function❞ (1.1) [Data set]. Zenodo. 
% https://zenodo.org/records/8330926 
% 4- Y. Lou, W. Zhou, T. P. Matthews, C. M. Appleton and M. A. Anastasio, ❝Generation of anatomically realistic
% numerical phantoms for photoacoustic and ultrasonic breast imaging❞, J. Biomed. Opt., vol. 22, no. 4, pp. 041015,
% 2017.
% 5 - B. E. Treeby and B. T. Cox, ❝k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic
% wave fields❞, J. Biomed. Opt. vol. 15, no. 2, 021314, 2010.
%
% You may also think about citing this preprint:
% 6- A Javaherian, ❝Full-waveform Approximation of Finite-Sized Acoustic Apertures:
% Forward and Adjoint Wavefields❞, https://arxiv.org/abs/2212.04466.
%
%
% author: Ashkan Javaherian
% date:            - 11.04.2020
% last update:     - 05.10.2022
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian


% the three essential c's
clear all
close all
clc

% run the startup script for defining the paths
startup_simulation_ust;

%% ========================================================================
% OPTIONAL PARAMETERS
%==========================================================================
% get the Boolean controlling whether the k-Wave simulation is performed and
% the simulated synthetic UST data is stored in a directory (true), or alternatively,
% the already simulated pressure time series are loaded from the same directory (false).
% For setting 'do_data_sim = false', the data muat have been stored in the
% associated path.

% The directories for saving the UST simulated by k-Wave are:
% ['.../simulation/data/data_ust_kWave_transmission/2D/name_data/file_name',
% '.../simulation/data/data_ust_kWave_transmission/3D/name_data/file_name'

% The folder ''data_ust_kWave_transmission'', which includes the k-Wave simulated data
% used in the examples in this project, can be downloaded via:
% https://zenodo.org/records/8330926.
%
%
% name_data is set in the function 'makeDataSettings.m' using the line
% commands:
% data_paths.name_data = ['Pulse' para.Excit  '_dx' num2str(1e4 * dx) '_cfl' num2str(10 * cfl) '_Nr'...
% num2str(num_receiver) '_Ne' num2str(num_emitter)  '_Interp' para.InterpType...
% '_Transgeom' para.TransGeom '_Absorption' num2str(para.Absorption) '_Code'   para.CodeVersion  '/'];


% For example,
% 1- a file_name 'data_4nonsmooth.m' means data which is simulated using
% k-Wave on a grid of size 0.4 mm using all 64 emitters.
% For this simulation, the sound speed and absorption coefficient maps
% are not smoothed, because it will be used for image reconstruction.
% 2- a file_name 'data_4smooth_17_20.m' means data which is
% simulated using k-Wave on a grid of size 0.4 mm using emitter 20. For this
% k-Wave simulation, the fields are smoothed by a window of size 17 grid points,
% equivalent to 7 grid points on a grid of size 1 mm, and is used for validation of
% the Green's approach in water analytically and in the heterogeneous medium
% (breast in water).
% Please see Section 'Numerical results' in the preprint:
% 'A. Javaherian, Hessian-inversion-free ray-born inversion....'
do_data_sim = false;


% get the Boolean controlling whether the plots are saved or not.
% This Boolean determines whether the plots associated with validation
% of the Green's approach in homogeneous or heterogeneous medium is stored
% or not.
save_plots = false;

%%=========================================================================
% GET THE PURPOSE OF RUNNING THE SCRIPT
%==========================================================================
% get the number of dimensions (2 or 3)
dim = 2;

% get the numer of workers for parallel programming
% the user may want to change the number of workers.
num_worker_pool = 0;


% In general, we have two options
% 1) Set:
% scenario = 'standard'; purpose = 'image_reconstruction';
% The UST data is simulated for all emitter-receiver pairs, and will be used
% for image reconstruction. For simulation of data, the sound speed and
% absorption coefficient maps are not smoothed.
% 2) Set...
% scenario = 'single_emitter'; purpose = 'raytracing_validation';
% The UST data is simulated for a single emitter, and is used for testing
% and validation of the ray tracing algorithm for calculating amplitude and
% phase.
% For this example script, the purpose is validation of ray tracing using
% a single emitter

% get the scenario, which is 'single_emitter'
scenario = 'single_emitter';

% get the purpose, which is 'raytracing_validation'.
purpose = 'raytracing_validation';



%%=========================================================================
% DEFAULTS USED IN THE PAPER, BUT IT CAN BE CHANGED BY THE USER
%%=========================================================================
% get the emitter index for validation of the ray approximation to the
% Green's function for homogeneous and heterogeneus media (phase and
% amplitude plots for the pressure field on all the receivers.)
% For emitters 1 and 20, the stored data is available at zenodo.
emitter_index = 1;   

% get the receiver index for validation of the ray approximation to the
% Green's function for homogeneous and heterogeneus media (phase and amplitude plot
% for a single-emitter-receiver pair and on all the discretised frequencies.)
receiver_index = 100;
   
%%=========================================================================
% THE FIXED PARAMETERS, THE PARAMETERS WHICH MUST NOT BE CAHNGED BY THE
% USER
%==========================================================================
% get the grid spacing [m] for ray tracing
grid_spacing_raytracing = 1e-3;

% the name of the used excitation pulse in our study
excit_pulse_name = 'Pammoth_1';

% get the cfl number to be used for time spacing for the k-Wave simulations
cfl_number = 0.1;

% the temporal downsampling rate 
time_downsampling_rate = 2;

% get the Boolean controlling whether the pressure source is deconvolved
% from the pressure time series or not. The presure time series are not
% deconvolved from the source, when the purpose is validation of ray
% tracing, because the plan is compare the amplitudes simlated using the
% k-Wave and approximated using the Green's formula on the receivers,.
deconvolve_source = false;

% Boolean controlling whether the time-of-flights are computed from simulated
% (or loaded) data, or not. For this example, because the purpose is
% validation of ray tracing, not image reconstruction, this parameter must
% be set false.
do_calculate_tofs = false;

% get the number of receivers
num_receiver = 2^8;

% get the haLf size [m] of the computational grid
half_grid_size = 10e-2;

% get the grid spacing [m] for the k-Wave simulations
grid_spacing_data = 4e-4;   %   (Default : 4e-4)

% get the Boolean controlling whether the acoustic absorption and dispersion are
% included or not
do_absorption = true;

% the approach for interpolation of the pressure field from the
% emitters to grid points and from grid points to receivers.
% Using 'offgrid', the offgrid toolbox is used. 
transducer_interp_approach = 'offgrid';

% get the downsampling rate for the emitters
emitter_downsampling_rate = 1;

% get the downsampling rate for the receivers
receiver_downsampling_rate = 1;

% get the ratio of the number of receivers to emitters
% It is not used when scenario is set 'single_emitter'.
ratio_receiver_to_emitter = 4;

% the approach for computing the geometrical attenuation, auxiliary: using
% auxiliary rays, or 'raylinked' (linked rays without using auxiliary rays)
attenuation_geom_method = 'auxiliary';

% the approach for tracing auxiliary rays, which will be used for computing
% the geomterical attenuation, 'paraxial' or 'angle_perturbation'.
auxiliary_method = 'paraxial';


if save_plots
    
    % get the directory for saving the plots
    plot_directory = 'results/simulation/2D/greens_plots/';
    
    % make the directory for saving the plots, if it does not exist
    makeDirectory(plot_directory);
else
    
    % allocate an empty variable for the directory
    plot_directory = [];
    
end


%%=================================================================
% SIMULATE THE ONLY-WATER AND BREAST-IN-WATER PRESSURE DATA ON RECEIVERS
% USING K-WAVE
%==================================================================
% get the optional inputs for defining the simulation settings, or simulating
% (or loading) the pressure data sets. The simulations are performed using
% the k-Wave toolbox.
data_args = {'num_worker_pool', num_worker_pool, 'grid_spacing_data', grid_spacing_data,...
    'grid_spacing_reconstruction', grid_spacing_raytracing, 'cfl', cfl_number,...
    'transducer_interp_approach', transducer_interp_approach,...
    'single_emitter_receiver', [nan, emitter_index], 'emitter_downsampling_rate',...
    emitter_downsampling_rate, 'receiver_downsampling_rate',...
    receiver_downsampling_rate, 'time_downsampling_rate', time_downsampling_rate,...
    'do_calculate_tofs', do_calculate_tofs};


% if the puropose is validation of ray tracing, the noise-free data is used.
[data, data_ref, ~, ~, emitter, receiver, simulation_prop, data_simulation_time] =...
    simulateSettingData(do_data_sim, dim, scenario, purpose, excit_pulse_name,...
    half_grid_size, num_receiver, ratio_receiver_to_emitter, do_absorption,...
    oa_breast_path, machine_name, res_path, local_res_path, data_args{:});


%%==================================================================
% COMPARE THE GREEN'S FUNCTION AND K-WAVE FOR HOMOGENEOUS
% (ONLY-WATER) MEDIUM
%===================================================================
% get the optional inputs for homogeneous Green's function
greens_hom_args = {'emitter_index', 1, 'receiver_index', receiver_index,...
    'deconvolve_source', deconvolve_source, 'save_plots', save_plots};

% get the relative discrepancies for the pressure field approximated
% by the Green's function and k-Wave, on the receivers (or single receiver)
% after being produced by a single emitter and propgating through only water
[relative_discrepancy_signal_water_emitter, relative_discrepancy_signal_water_receiver] = ...
    validateGreensHomogeneous(data_ref, simulation_prop.t_array, emitter, receiver,...
    simulation_prop.sound_speed_ref, simulation_prop.f_max, plot_directory, greens_hom_args{:});

%%=================================================================
% COMPARE THE GREEN'S FUNCTION AND K-WAVE FOR HETEROGENEOUS
% (BREAST-IN-WATER) MEDIUM
%===================================================================
% get the optional inputs for ray approaximtion to heterogeneous Green's function
greens_het_args = {'num_worker_pool', num_worker_pool, 'emitter_index',...
    1, 'receiver_index', receiver_index,...
    'deconvolve_source', false, 'save_plots', save_plots, ...
    'auxiliary_method', auxiliary_method,...
    'attenuation_geom_method', attenuation_geom_method};


% compare the pressure field approximated by the Green's function and k-Wave
% on the receivers after being produced by the chosen single emitter and 
% propagation through the breast in water
% In addition, get the parameters on the sampled points along the forward and
% adjoint (backprojection) rays.
[ray_position, ray_time, ray_absorption, ray_position_left, ray_position_right,...
    adjoint_ray_position_left, adjoint_ray_position_right, ray_spacing,...
    rayspacing_receiver, relative_discrepancy_absorbing, relative_discrepancy_nonabsorbing] =...
    validateGreensHeterogeneous(data, data_ref, simulation_prop.t_array, emitter, receiver,...
    simulation_prop, [], plot_directory, greens_het_args{:});

%%=================================================================
% DISPLAY THE RAYS AFTER VALIDATION OF RAY TRACING
%==================================================================


% get the optional inputs
% emitter index must be set  zero, because emitter.positions now includes
% only the position of one emitter.
display_args = {'save_plots', save_plots, 'emitter_index', 1,...
    'receiver_index', receiver_index};

% display the forward ray initialised from a single emitter and
% adjoint rays initilialised from a single receiver
displayRays(emitter.positions, receiver.positions, ray_position, ray_time,...
    ray_absorption, ray_position_left, ray_position_right,...
    adjoint_ray_position_left, adjoint_ray_position_right, simulation_prop,...
    ray_spacing, rayspacing_receiver, plot_directory, display_args{:});



      