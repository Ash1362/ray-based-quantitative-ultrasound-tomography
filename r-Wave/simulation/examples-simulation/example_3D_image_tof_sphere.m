% Example for the image-reconstruction part of the scenarios in paper
% [1] (cf. section References)
%............................................................................................
% This work was sucessfully implemented on the Pammoth systemhttps (www.pammoth-2020.eu/),
% leading to significant contrast between using straight and bent rays.

%%=========================================================================
% SIMULATION (OR LOAD) OF UST DATA USING K-WAVE
%==========================================================================
% Two ultrasound datasets, one for breast inside water and another for only water,
% were simulated based on the scenarios in section 5.2-5.5 in [1]. The simulation of
% the two only-water and breast-in-water data sets using the k-Wave [4] on a grid with
% grid spacing 0.5mm took a week on 8 GPUs. Therefore, for the 3D case, the data were
% already simulated and time-of-flights were computed using the"Akaike information criterion"
% algorithm. The time-of-flight data are available at the link provided in the reference [3],
% i.e., https://zenodo.org/records/8330926.
%%=========================================================================
% RAY TRACING 
%==========================================================================
% This script will then compute a low-resolution image of the sound speed
% from time-of-flights of two simulated transmission ultrasoud datasets for
% only-water and breast-in-water cases.
% For construction of the system matrix associated with each linearization
% in the time-of-flight-based image reconstruction, ray tracing algorithms
% described in section 3 in [1] and ray linking algorithm described in section 4
% in [1] were used.
% For ray tracing, the user can choose between ray tracing algorithms:
% 1- 'Dual-update' (section 3.3.1 in [1])
% 2- 'Mixed-step' (section 3.3.2. i [1])
% 3- 'Charactersitics' (section 5.1, Eq. (54) in [2])
% 4 -'Rung-kutta-2nd' (section 5.1, Eq (54) using Algorithm 2 in [2])
% ([3] and [4] were added after publishing the paper [1].

%%=========================================================================
% RAY-TO-GRID/GRID-TO-RAY INTERPOLATION
%==========================================================================
% 1- 'Bilinear' (section 3.4. in [1])
% 2- 'Bspline' (section 5.3. in [2])

% 'Bspline' is more accurate than 'Bilinear', but it will make ray tracing slower.
% 'Bilinear' was used for implementation on the Pammoth system, but 'Bspline' was used
% as default in this simulation example.
%%=========================================================================
% RAY LINKING 
%==========================================================================
% For ray linking, the user is limited to the best algorithm in the paper.
% The chosen ray linking algorithm uses a Quasi-Newton algorithm based on a
% derivative-free ❝Broyden❞ approach. 

%%=========================================================================
% IMAGE RECONSTRUCTION
%==========================================================================
% For image reconstruction, the objective function, which is an L2 norm of
% the discrepancy of the time-of-flights associated with the breast-in-water
% and only-water data is iteratively linearised, and each linearised subproblem
% is solved using either 'sart', or 'conjugate_gradient' algorithms. These
% algorithms converge faster then the steepest descent, and require maximum 5
% inner iterations for soving each linearised subproblem.  (default: 'sart')
% The image reconstruction algorithm is terminated after 12 linearisations.


%%=========================================================================
% DIGITAL PHANTOM AND K-WAVE
%==========================================================================
% A digital breast phantom developed by team of Professor Mark Anastasio is used
% in this example script [4]. The phantom data must be downloaded via:'...
% https://anastasio.bioengineering.illinois.edu/downloadable-content/oa-breast-database/ [3].
% The folder 'Neg_47_Left' is for simulation of UST data for all studies in this toolbox.
% The folder 'Neg_47_Left' must be added to the directory:
% '.../simulation/data/phantom/OA-BREAST/Neg_47_Left...'
% The pressure field used as the benchmark is simulated using the k-Wave
% toolbox.  www.k-Wave.org (Version 1.3. or 1.4.) [5].
% 
%%=========================================================================
% REFERENCES
%==========================================================================
% If you find this example script useful for your research, please consider citing this paper:
% 1- A. Javaherian, F. Lucka and B. T. Cox, ❝Refraction-corrected ray-based inversion for three-dimensional
% ultrasound tomography of the breast❞, Inverse Problems, 36 125010.
% 2 - A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering
% for biomedical ultrasound tomography❞, Inverse Problems vol. 37, no.11, 115003, 2021.

% These data-base/toolbox were used in this project.
% 3 - A. Javaherian, 2023, ❝Transmission ultrasound data simulated using the k-Wave toolbox as a benchmark for biomedical quantitative ultrasound tomography using a ray approximation to Green's function❞ (1.1) [Data set]. Zenodo. 
% https://zenodo.org/records/8330926 
% 4- Y. Lou, W. Zhou, T. P. Matthews, C. M. Appleton and M. A. Anastasio, ❝Generation of anatomically realistic
% numerical phantoms for photoacoustic and ultrasonic breast imaging❞, J. Biomed. Opt., vol. 22, no. 4, pp. 041015,
% 2017.
% 5- B. E. Treeby and B. T. Cox, ❝k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic
% wave fields❞, J. Biomed. Opt. vol. 15, no. 2, 021314, 2010.
%
%
% author: Ashkan Javaherian
% date:            - 29.8.2019
% last update:     - 05.10.2022
%
% This script is part of the r-Wave Toolbox
% Copyright (c) 2022 Ashkan Javaherian


% the three essential c's
clear all
close all
clc

% run the startup script for defining the paths
startup_simulation_ust;

%%=========================================================================
% THE PARAMETERS THAT CAN BE CHANGED BY THE USER BASED on the description
%==========================================================================
% get the grid spacing [m] for image reconstruction
% The grid spacing can be either 1e-3 m, or 2e-3 m. (Default :2e-3 m)
% The grid spacing 2mm is good choice for 3D, and was used for implemenetation
% on the Pammoth system. Note that an extension of
% image reconstruction using the Green' approach to the 3D case is not
% completed yet. For 3D case, only a time-of-flight image is reconstructed
% using the paper mentioned above.
grid_spacing_reconstruction = 2e-3;

% the method for reconstruction of the system matrix after each
% linearisation of the objective function, which is an L2 norm of the
% misfit between the only-water and object-in-water time-of-flight data.
% This can be set 'straight-ray', or 'bent-ray'. Using 'straight-ray',
% the system matrix is computed once, but using 'bent-ray'
% approach, the system matrix is iteratively computed for each linearised
% subproblem, based on the last update of the refractive index.
matrix_construction_method = 'bent-ray';

% choose the approach used for ray tracing for image reconstruction using
% time-of-flight data. This can be set: 1)'Mixed-step', 2)'Dual-update',
% 3)'Characteristics', or 4)'Runge-kutta-2nd'. (default:'Mixed-step')
% The 'Mixed-step' algorithm is quite fast, and was used for implementation
% on the Pammoth system.
raytracing_method_tof = 'Runge-kutta-2nd';  

% get the numer of workers for parallel programming
% the user may want to change the number of workers.
num_worker_pool = 16;

% get the downsampling rate for the emitters. This can be set 1, 2, 4, 8.
% 8 gives the faster image reconstruction. (Default: 2)
emitter_downsampling_rate = 2;

% get the downsampling rate for the receivers. This can be set 1, 2, 4, 8.
% 8 gives the faster image reconstruction. (Default: 2)
receiver_downsampling_rate = 2;


%% ========================================================================
% OPTIONAL PARAMETERS
%==========================================================================
% get the Boolean controlling whether the mat results are saved or not
% This Boolean determines whether the mat files associated with the images
% reconstructed using the time-of-flight-based or Green's approach are
% saved or not.
save_results = true;

% get the ray-to-grid and grid-to-ray interpolation approach for
% time-of-flight-based image reconstruction approach. this can be set
% 'Bilinear' or 'Bspline'
gridtoray_interp_tof = 'Bspline';

% get the appproach for solving each linearised subproblem. Using a
% lineraisation using the 'absloute' approach, the arisng linearised
% subproblems can be solved using the 'sart' or 'conjuge_gradient'
% algorithm.
linear_subproblem_method = 'sart';

%%=========================================================================
% THE FIXED PARAMETERS, THE PARAMETERS WHICH MUST NOT BE CHANGED BY THE
% USER, BECAUSE THEY CORRESPOND TO (CONSISTENT WITH) THE STORED
% TIME-OF-FLIGHTS.
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
do_data_sim = false;

% Boolean controlling whether the time-of-flights are computed from simulated
% (or loaded) data, or the already computed tofs are used.
% By setting true, the time-of-flights are computed from the synthetic UST
% data and stored at the same directory for storing UST data. Please read the description
% for data_sim
% By setting false, the already computed time-of-flights, which are stored
% at the same directory for storing the UST data, are loaded and used for
% image reconstruction. For this example script, the already computed
% time-of-flights will be used.
do_calculate_tofs = false;

% get the number of dimensions (2 or 3)
dim = 3;

% get the Boolean controlling whether the image reconstruction (or validation
% of ray tracing) using the Green's approach is performed, or not.

% The Green's approach has not been extended to 3D case yet.
% For 3D case, only a time-of-flight image is reconstructed using
% the bent rays.
do_greens_approach = false;

% the signal-to-noise ratio of the simulated ultrasound data
% For this example script, the available time-of-flights correspond to a 40 dB SNR
% (fixed)
noise_level = 40;

% scenario 'standard' means the k-Wave simulation is performed or has been performed
% for all emitter-receiver pairs.'single_emitter' means simulation for a single emitter
% for testing purposes. For this example script, the scenario must be set 'standard'
scenario = 'standard';

% For this example script, the purpose is either 'image_reconstruction'.
purpose = 'image_reconstruction';

% the name of the used excitation pulse in our study
excit_pulse_name = 'Pammoth_1';

% get the cfl number to be used for time spacing for the k-Wave simulations
cfl_number = 0.1;

% the temporal downsampling rate used for (partially) avoiding an inverse crime
% in the time spacing. (Note that because this study uses two inherently different
% approaches for data simulation and image reconstruction, the inverse
% crime is not the case, but the author still considered these points (avoiding
% inverse crime for time sampling, grid psacing, realsitic pressure source,...) for
% making the study reliable as much as possible.)
time_downsampling_rate = 2;

% get the number of receivers
num_receiver = 2^12;

% get the half size [m] of the computational grid along the x-y
% plane
half_grid_size = 13e-2;

% get the grid spacng [m] for the k-Wave simulations
grid_spacing_data = 5e-4;

% get the Boolean controlling whether the acoustic absorption and dispersion are
% included or not
do_absorption = false;

% get the approach for interpolation of the pressure field from the
% emitters to grid points and from grid points to receivers.
% Using 'offgrid', the offgrid toolbox is used. (c.f. reference [5]
% in the description of the script.)
transducer_interp_approach = 'nearest';

% get the ratio of te number of receivers to the number of emitters
ratio_receiver_to_emitter = 4;

% get the approach for linearisation. This can only be set 'absolute'.
linearisation_approach = 'absolute';

if save_results
    
    % get the directory for saving the mat results
    results_directory = 'results/simulation/3D/tof/';
    
    % make the directory, if it does not exist
    makeDirectory(results_directory);
    
else
    
    % allocate an empty variable for the directory
    results_directory = [];
    
end

%%=================================================================
% SIMULATE PRESSURE UST DATA SETS USING K-WAVE OR GET THE SETTINGS
% OF THE SIMULATION
%==================================================================
% get the optional inputs for defining the settings which are already set
% for simulating the pressure data sets using the k-Wave. Because do_calculate_tofs
% is set false, this function is used for only reproducing the simulation settings
% corresponding to the stored Time-of-flight data.
data_args = {'num_worker_pool', num_worker_pool, 'grid_spacing_data', grid_spacing_data,...
    'grid_spacing_reconstruction', grid_spacing_reconstruction, 'cfl', cfl_number,...
    'transducer_interp_approach', transducer_interp_approach,...
    'emitter_downsampling_rate', emitter_downsampling_rate, 'receiver_downsampling_rate',...
    receiver_downsampling_rate, 'time_downsampling_rate', time_downsampling_rate,...
    'do_calculate_tofs', do_calculate_tofs, 'noise_level', noise_level};


% if the purpose is image reconstruction, the noise-contaminated
% data is used.
[~, ~, data_noisy, data_ref_noisy, emitter, receiver, simulation_prop,...
    data_simulation_time] =...
    simulateSettingData(do_data_sim, dim, scenario,...
    purpose, excit_pulse_name, half_grid_size, num_receiver,...
    ratio_receiver_to_emitter, do_absorption, oa_breast_path, machine_name,...
    res_path, local_res_path, data_args{:});


%%=================================================================
% RECONSTRUCT THE SOUND SPEED IMAGE USING THE TIME-OF-FLIGHT-BASED
% APPROACH
%==================================================================
% get the optional inputs
reconst_args_tof = {'num_worker_pool', num_worker_pool,...
    'reconstruct_image', true, ...
    'matrix_construction_method', matrix_construction_method,...
    'linear_subproblem_method', linear_subproblem_method,...
    'linearisation_approach', linearisation_approach,...
    'raytracing_method', raytracing_method_tof,...
    'gridtoray_interp', gridtoray_interp_tof,...
    'grid_spacing', grid_spacing_reconstruction,...
    'binaries_emitter_receiver', 'distances',...
    'emitter_downsampling_rate', emitter_downsampling_rate,...
    'receiver_downsampling_rate', receiver_downsampling_rate};


if ~do_calculate_tofs
    simulation_prop.t_array = [];
end


% load the time-of-flight data from the directory for storing data and
% reconstruct the sound speed image iteratively from the time-of-flight data
[img_tof, recon_grid, ~, ray_initial_angles, out_tof, para_tof] = ...
    reconstructTimeofFlightImage(data_noisy, data_ref_noisy,...
    simulation_prop.t_array, emitter, receiver, simulation_prop.sound_speed_ref,...
    simulation_prop, simulation_prop.data_path,...
    reconst_args_tof{:});

if save_results
    
    % save the results in the given path
    save([results_directory 'results_' greens_optimisation_approach  '_' absorption_map...
        '_' num2str(noise_level) 'db'],...
        'img_tof', 'out_tof', 'para_tof', 'recon_grid', '-v7.3');
    
end


        