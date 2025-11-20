% Example for the image reconstruction of the sound speed from ultrasound
% data simulated on a number of rotating linear arrays (2D case)

%%=========================================================================
%  SIMULATION (OR LOAD) OF UST DATA USING K-WAVE AND INITIAL GUESS
%==========================================================================
% This example first simulates (or loads) two ultrasound data sets, one for
% breast inside water and another for only water.The datasets can be either
% simulated or downlaoded from the Zenodo link provided in the referene [2].
% 16 line arrays are placed on the periphery of a circle, with 16 transducers
% placed on each line array.
% The time-of-flights of the UST time series were then computed using the 
% Akaike-information-criterion approach. A sound speed image is then reconstructed
% from the time-of-flights using few iterations (linearizations) of a time-of-flight
% algorithm. Each ToF-based linearised subproblem is solved using 'SART' inner 
% iterations. 
% For construction of the system matrix associated with each linearisation in
% the time-of-flight-based image reconstruction approach, the ray tracing can
% be done using the four algorithms developed in this toolbox.

%%=========================================================================
% DIGITAL PHANTOM AND K-WAVE
%==========================================================================
% The Breast is simulated using a digital breast phantom developed by the group
% of Professor Mark Anastasio [3]. The phantom data must be downloaded via:'...
% https://anastasio.bioengineering.illinois.edu/downloadable-content/oa-breast-database/ [3].
% The synthetic UST pressure times series used in this study are simulated using the k-Wave
% toolbox. www.k-Wave.org (v. 1.3. or 1.4.) [4].

%
%%=========================================================================
% REFERENCES
%==========================================================================
% If you find the toolbox useful for your research, please consider citing these papers:...
% 1 - A. Javaherian et al., ❝Refraction-corrected ray-based inversion
% for three-dimensional ultrasound tomography of the breast❞, Inverse Problems, 36 125010. https://iopscience.iop.org/article/10.1088/1361-6420/abc0fc/  
% 
% These data-base/toolbox were used in this project.
% 2 - A. Javaherian, 2023, ❝Transmission ultrasound data simulated using the k-Wave toolbox as a benchmark for biomedical quantitative ultrasound tomography using a ray approximation to Green's function❞ (1.1) [Data set]. Zenodo. 
% https://zenodo.org/records/8330926 
% 3 - Y. Lou, W. Zhou, T. P. Matthews, C. M. Appleton and M. A. Anastasio, ❝Generation of anatomically realistic
% numerical phantoms for photoacoustic and ultrasonic breast imaging❞, J. Biomed. Opt., vol. 22, no. 4, pp. 041015,
% 2017.
% 4 - B. E. Treeby and B. T. Cox, ❝k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic
% wave fields❞, J. Biomed. Opt. vol. 15, no. 2, 021314, 2010.
%
% You may also think about citing this preprint:
% 4- A Javaherian, ❝Full-waveform Approximation of Finite-Sized Acoustic Apertures:
% Forward and Adjoint Wavefields❝, https://arxiv.org/abs/2212.04466.
%
%
% author: Ashkan Javaherian
% date:            - 17.02.2020
% last update:     - 05.07.2023
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian


% the three essential c's
clear all
close all
clc

% run the startup script for defining the paths
startup_simulation_ust;

%%=========================================================================
% THE PARAMETERS THAT CAN BE CHANGED BY THE USER BASED ON THEIR DESCRIPTIONS
%==========================================================================
% the signal-to-noise ratio of the simulated ultrasound data
% This can be set 40, 30, or 25
noise_level = 40;

% get the optimisation approach for image reconstruction
% using the ray approaximation to heterogeneous Green's function
% This can be set 'hessian', or 'backprojection'.
% 'backprojection' is based on the Green's approach proposed in [1].
% 'hessian' is based on the Green's approach proposed in [2].
greens_optimisation_approach = 'backprojection';

% the type of the absorption coefficient map used for image reconstruction
% using the ray approaximation to heterogeneous Green's function. This can
% be 'true', 'homogeneous', or 'none'
% Please read section Numerical results in [1]!
absorption_map = 'true';

% get the grid spacing [m] for image reconstruction, 1e-3 or 2e-3.
% (Default: 1e-3)
% The grid spacing 1e-3 m is used in the papers, and is recommended.
% If your access to memory is limited, you can also use 2e-3 m.
grid_spacing_reconstruction = 1e-3;

% get the numer of workers for parallel programming
% the user may want to change the number of workers.
num_worker_pool = 16;

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

% get the Boolean controlling whether the mat results are saved or not
% This Boolean determines whether the mat files associated with the images
% reconstructed using the time-of-flight-based or Green's approach are
% saved or not.
save_results = true;


%%=========================================================================
% DEFAULTS USED IN THE PAPER, BUT IT CAN BE CHANGED BY THE USER.
%%=========================================================================
% choose the approach used for ray tracing for image reconstruction using
% time-of-flight data. This can be set: 'Mixed-step', 'Dual-update',
% 'Characteristics', or 'Runge-kutta-2nd'.(Default: 'Runge-kutta-2nd')
raytracing_method_tof = 'Runge-kutta-2nd'; %'Runge-kutta-2nd';

%%=========================================================================
% THE FIXED PARAMETERS, THE PARAMETERS WHICH MUST NOT BE CAHNGED BY THE
% USER
%==========================================================================
% get the number of dimensions
dim = 2;

% In general, you have two options
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

% In this eaxample script, the purpose is image reconstruction (not validation of
% ray tracing), and the UST data simulated for all emitter-receiver pairs
% (not only a single emitter) must be used.

scenario = 'standard';

% the purpose is either 'image_reconstruction',
purpose = 'image_reconstruction';

% get the detection geomtery
detection_geom = 'plane';

% get the Boolean controlling whether the image reconstruction using the Green's
% approach is performed, or not.
do_greens_approach = false;

% the name of the used excitation pulse in our study
excit_pulse_name = 'Pammoth_1';

% get the cfl number to be used for time spacing for the k-Wave simulations
cfl_number = 0.1;

% the temporal downsampling rate used for (partially) avoiding an inverse crime
% in the time spacing
time_downsampling_rate = 2;

% get the Boolean controlling whether the pressure source is deconvolved
% from the pressure time series or not. If purpose = 'image_reconstruction',
% we set  deconvolve_source = true.
deconvolve_source = true;

% Boolean controlling whether the time-of-flights
% are computed from simulated (or loaded) data,
% or the already computed tofs are used.
do_calculate_tofs = true;

% get the ray-to-grid and grid-to-ray interpolation approach for
% time-of-flight-based image reconstruction approach
gridtoray_interp_tof = 'Bspline';

% get the number of receivers 
% the number of receivers must be the same as the number of the stored
% element positions. It is not an optional imput and must not be changed 
% by the user.
num_receiver = 2^8;

% get the half size [m] of the computational grid
half_grid_size = 11e-2;

% get the grid spacing [m] for the k-Wave simulations
grid_spacing_data = 4e-4;     %   (Default : 4e-4)

% get the Boolean controlling whether the acoustic absorption and dispersion are
% included or not
do_absorption = true;

% get the number of linearised problems (outer iterations) for
% reconstructing an image from the time of flights. (Note that the number
% of iterations for reconstructing the optimal image is 12-15, but for
% providing an initial guess for the Green's approach, only few iterations
% of the 'sart' algorithm (simulateneous algebraic reconstruction algorithm)
% is sufficient for providing an image accurate enough and with minimal artefact.
num_iterout_tof = 5;

% the approach for interpolation of the pressure field from the
% emitters to grid points and from grid points to receivers.
% Using 'offgrid', the offgrid toolbox is used. (c.f. reference [5]
% in the description of the script.)
transducer_interp_approach = 'offgrid';

% get the downsampling rate for the emitters
emitter_downsampling_rate = 1;

% get the downsampling rate for the receivers
receiver_downsampling_rate = 1;

% get the ratio of te number of receivers to the number of emitters
ratio_receiver_to_emitter = 4;

% get the approach for linearisation. This cab be 'absolute' or
% 'difference'. The latter is deprecated.
linearisation_approach = 'absolute';

% get the appproach for solving each linearised subproblem. Using a
% lineraisation using the 'absloute' approach, the arisng linearised
% subproblems can be solved using the 'sart' or 'conjuge_gradient'
% algorithm.
linear_subproblem_method = 'sart';

% the approach fpr computing the geomterical attenuation, auxiliary: using
% auxiliary rays, or 'raylinked' (linked rays without using auxiiary rays)
attenuation_geom_method = 'auxiliary';

% the approach for tracing auxiliary rays, which will be used for computing
% the geomterical attenuation, 'paraxial' or 'angle_perturbation'.
auxiliary_method = 'paraxial';


if save_results
    
    % get the directory for saving the mat results
    results_directory = 'results/simulation/2D/tof/';
    
    % make the directory, if it does not exist
    makeDirectory(results_directory);
    
else
    
    % allocate an empty variable for the directory
    results_directory = [];
    
end


%%=================================================================
% DEFINE THE SIMULATION SETTINGS AND SIMULATE (OR LOAD) THE PRESSURE UST DATA
% SETS USING K-WAVE
%==================================================================
% get the optional inputs for defining the settings for simulation using
% the k-Wave, or loading the pressure data sets already simulated using the k-Wave
% toolbox and with the same settings.
data_args = {'num_worker_pool', num_worker_pool, 'grid_spacing_data', grid_spacing_data,...
    'grid_spacing_reconstruction', grid_spacing_reconstruction, 'detection_geom', detection_geom,...
     'cfl', cfl_number, 'transducer_interp_approach', transducer_interp_approach,...
    'emitter_downsampling_rate', emitter_downsampling_rate,...
    'receiver_downsampling_rate', receiver_downsampling_rate,...
    'time_downsampling_rate', time_downsampling_rate,...
    'do_calculate_tofs', do_calculate_tofs, 'noise_level', noise_level,...
    'step_length', 1};

% if the purpose is image reconstruction, the noise-contaminated
% data is used.
[~, ~, data_noisy, data_ref_noisy, emitter, receiver,...
    simulation_prop, data_simulation_time] =...
    simulateSettingData(do_data_sim, dim, scenario,...
    purpose, excit_pulse_name, half_grid_size, num_receiver,...
    ratio_receiver_to_emitter, do_absorption, oa_breast_path, machine_name,...
    res_path, local_res_path, data_args{:});

%%=========================================================================
% RESHAPE THE data BASED ON ROTATION INDICES
%==========================================================================
% The full-wave simulation for each emitter per array has been done
% simultaneously for all receiver arrays. But the rays initialised at each emitter
% per array is traced for each receiver array separately. That means for
% each emitter per array, the ray tracing must be done num_array times.
% This is because given an emitter position, the ray linking problem can be
% solved using a single equation for the receivers, i.e., a single array.

% get the number of time instants in data
num_time = size(data_noisy, 2);

% get the number of receivers for each angular position
num_receiver = size(data_noisy, 1)/max(emitter.rotation_indices);

% get the number of emitters
num_emitter = length(emitter.rotation_indices);

% reshape the object-in-water data
data_noisy = permute(reshape(permute(data_noisy, [1,3,2]),...
    [num_receiver, num_emitter, num_time]), [1,3,2]);

% reshape the only-water data
data_ref_noisy = permute(reshape(permute(data_ref_noisy, [1,3,2]),...
    [num_receiver, num_emitter, num_time]), [1,3,2]);


%%=================================================================
% RECONSTRUCT THE SOUND SPEED IMAGE USING THE TIME-OF-FLIGHT-BASED
% APPROACH
%==================================================================
% get the optional inputs
reconst_args_tof = {'num_worker_pool', num_worker_pool,...
    'reconstruct_image', true, ...
    'matrix_construction_method', 'bent-ray',...
    'linear_subproblem_method', linear_subproblem_method,...
    'linearisation_approach', linearisation_approach,...
    'raytracing_method', raytracing_method_tof,...
    'gridtoray_interp', gridtoray_interp_tof,...
    'grid_spacing', grid_spacing_reconstruction,...
    'binaries_emitter_receiver', 'distances',...
    'num_iterout', num_iterout_tof,...
    'emitter_downsampling_rate', emitter_downsampling_rate,...
    'receiver_downsampling_rate', receiver_downsampling_rate};

% compute the time-of-flights from ultrasound data sets for the
% breast-in-water and only-water data, and reconstruct the sound
% speed image iteratively from the time-of-flight data
[img_tof, recon_grid, ~, ray_initial_angles, out_tof, para_tof] = ...
    reconstructTimeofFlightImage(data_noisy, data_ref_noisy,...
    simulation_prop.t_array, emitter, receiver, simulation_prop.sound_speed_ref,...
    simulation_prop, simulation_prop.data_path, reconst_args_tof{:});

% display the reconstructed image
figure; imagesc(img_tof); axis image; colorbar;

if save_results
    
    % save the results in the given path
    save([results_directory, 'results_', greens_optimisation_approach,  '_' absorption_map,...
        '_' num2str(noise_level) 'db'],...
        'img_tof', 'out_tof', 'para_tof', 'recon_grid', '-v7.3');
    
end


if do_greens_approach

% error (Not implemented yet!)  
error(['The Greens approach for the planar geometry has not been computed yet.'])
    
%%=================================================================
% RECONSTRUCT THE SOUND SPEED IMAGE USING THE RAY APPROAXIMATION
% TO HETEROGENEOUS GREENS FUNCTION ITERATIVELY
%==================================================================
% get the directories for storing the results associated with the image
% reconstruction using the Green's approach
directories_greens = [];
directories_greens.results_directory = results_directory;

% get the optional inputs
greens_args = {'num_worker_pool', num_worker_pool,...
    'optimisation_approach', greens_optimisation_approach,...
    'absorption_map', absorption_map,...
    'deconvolve_source', deconvolve_source,...
    'noise_level', noise_level};

% compute the sound speed image using the ray approximation to
% the heterogeneous Green's function. The image is
% reconstructed iteratively from low to high frequencies.
[img_greens, out_greens, para_greens] = reconstructGreensImage(data_noisy,...
    simulation_prop.t_array, recon_grid, emitter, receiver, simulation_prop.sound_speed_ref,...
    img_tof, ray_initial_angles, simulation_prop, directories_greens, greens_args{:});

if save_results
    
    % save the results in the given path
    save([results_directory 'results_' greens_optimisation_approach  '_' absorption_map...
        '_' num2str(noise_level) 'db'],...
        'img_greens', 'out_greens', 'para_greens', '-append');
    
end

end


