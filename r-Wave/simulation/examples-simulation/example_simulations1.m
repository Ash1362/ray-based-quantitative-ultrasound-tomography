% Example for reproduction of the results for the image-reconstruction portion of the scenarios in paper:
% A. Javaherian and B. Cox, Ray-based inversion accounting for scattering for biomedical ultrasound tomography,
% Inverse Problems vol. 37, no.11, 115003, 2021.
%
% This example first simulates (or loads) two ultrasound data sets, one for
% breast inside water and another for only water based on the scenarios in
% paper [2]. The script will then computes an image from the time-of-flights of
% two simulated (or loaded) ultrasoud time series using the early iterations of a 
% time-of-flight-based algorithm.
% The reconstructed image using the time-of-flights of the ultrasound
% pressure time series will then be used as the initial guess for an
% inversion approach based on the ray approaximation to heterogeneous Green's
% function. Using the ray-based Green's inversion approach, the image
% reconstruction is performed from low to high frequencies.

% The Breast is simulated
% using a digital breast phantom developed by Mark Anastasio group available at:'...
% https://anastasio.bioengineering.illinois.edu/downloadable-content/oa-breast-database/ [3].
% The presssuer field used as the benchmark is simulated using the k-Wave
% toolbox. www.k-Wave.org [4].
% For our study in [2], the interpolation of the pressure filed from
% from emitters to grid points and from grid points to receivers are simulated using
% the off-grid toolbox [5].
% 'Therefore,i f you find the toolbox useful for your resarch: please cite these papers:'...
% 1 - A Javaherian, F Lucka and B T Cox, Refraction-corrected ray-based inversion for three-dimensional
% ultrasound tomography of the breast, Inverse Problems, 36 125010.
% 2 -  A. Javaherian and B. Cox, Ray-based inversion accounting for scattering for biomedical ultrasound tomography,
% Inverse Problems vol. 37, no.11, 115003, 2021.
% Please also add the other toolboxes used in our tudy.
% 3- Y. Lou, W. Zhou, T. P. Matthews, C. M. Appleton and M. A. Anastasio, Generation of anatomically realistic
% numerical phantoms for photoacoustic and ultrasonic breast imaging, J. Biomed. Opt., vol. 22, no. 4, pp. 041015,
% 2017.
% 4 - B. E. Treeby and B. T. Cox, k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic
% wave fields, J. Biomed. Opt. vol. 15, no. 2, 021314, 2010.
% 5 -E. S. Wise, B. T. Cox, J. Jaros, B. E. Treeby, Representing arbitrary acoustic source and sensor distributions in
% Fourier collocation methods, J. Acoust. Soc. of Am., vol. 146, no. 1, pp. 278-288, 2019.
% 
%
% author: Ashkan Javaherian
% date:            - 29.12.2019
% last update:     - 05.10.2022
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2020 Ashkan Javaherian


% the three essential c's
clear all
close all
clc

% run the startup script for defining the paths
startup_simulation_ust;

%%=========================================================================
% THE PARAMETERS THAT CAN BE CHANGED BY THE USER
%==========================================================================
% the type of the absorption coefficient map used for image reconstruction
% using the ray approaximation to heterogeneous Green's function. This can
% be 'true', 'homogeneous', or 'none'
absorption_map = 'true';

% the optimisation approach for image reconstruction
% using the ray approaximation to heterogeneous Green's function
% This can be set 'hessian', or 'backprojection'
greens_optimisation_approach = 'hessian';  

% the signal-to-noise ratio of the simulated ultrasound data
noise_level = 40;

%% ========================================================================
% THE FIXED INPUT PARAMETERS
%==========================================================================
% get the numer of workers for parallel programing
% the user may want to change the number of workers. 
num_worker_pool = 16;
% get the Boolean controlling whether the k-Wave simulation is performed and stored in
% a directory, or alternatively, the already simulated pressure time series are loaded
% from the same directory.
do_data_sim = true;

% get the Boolean controlling whether the mat results are saved or not
save_results = false;

% get the Boolean controlling whether the plots are saved or not
save_plots = false;

% get the nuber of dimensions (2 or 3)
dim = 2;

% get the scenario ('standard': the k-Wave simulation is performed for all
% emitter-receiver pairs, or 'single_emitter': a single emitter for testing
% purposes.)
scenario = 'single_emitter';            % 'standard';            %              

% the purpose is either 'image_reconstruction', or 'raytracing_validation'.
purpose = 'raytracing_validation';      % 'image_reconstruction';    %   

if strcmp(purpose, 'raytracing_validation')
        
        
    % get the emitter index for validation of the ray approximation to the 
    % Green's function for homogeneous and heterogeneus media (phase and
    % amplitude plots for the pressure field on all the receivers.)
    emitter_index = 1;
    
    % get the receiver index for validation of the ray approximation to the
    % Green's function for homogeneous and heterogeneus media (phase and amplitude plot
    % for a single-emitter-receiver pair and on all the discretised frequencies.)
    switch dim
        case 2
            receiver_index = 100;
        case 3
            receiver_index = 1000;
    end
    
else
    emitter_index = [];
            
end


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


% get the Boolean controlling whether the pressure source is deconvolved
% from the pressure time series or not.
switch purpose
    case 'raytracing_validation'
        
        deconvolve_source = false;
        
    case 'image_reconstruction'
        
        deconvolve_source = true;
        
end




switch dim
    case 2
        
        % Boolean controlling whether the time-of-flights
        % are computed from simulated (or loaded) data,
        % or the already computed tofs are used.
        switch purpose
            case 'raytracing_validation'
                
                do_calculate_tofs = false;
                
            case 'image_reconstruction'
                
                do_calculate_tofs = true;
                
        end
        

        % get the Boolean controlling whether the image reconstruction (or validation
        % of ray tracing) using the Green's approach is performed, or not.
        do_greens_approach = true;
        
        % get the ray-to-grid and grid-to-ray interpolation approach for
        % time-of-flight-based image reconstruction approach
        raytogrid_interp_tof = 'Bspline';
        
        % get the number of receivers
        num_receiver = 2^8;
        
        % get the haLf size [m] of the computational grid
        half_grid_size = 10e-2;
        
        % get the grid spacing [m] for the k-Wave simulations
        grid_spacing_data = 1e-3;   %   (Default : 4e-4) 
        
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
        
        % get the grid spacing [m] for image recosntruction
        grid_spacing_reconstruction = 1e-3;
        
    case 3
        
        
        % give an error if dim=3, and the purpose is raytracing_validation
        % if strcmp(purpose, 'raytracing_validation')
        %    error(['The extension of ray approximation to heterogeneous Greens'...
        %        'function to 3D is not completed yet.'])
        %end
           
        % Boolean controlling whether the time-of-flights
        % are computed from simulated (or loaded) data,
        % or the already computed tofs are used.
        do_calculate_tofs = false;
        
        % get the Boolean controlling whether the image reconstruction (or validation
        % of ray tracing) using the Green's approach is performed, or not.
        do_greens_approach = false;
        
        % get the ray-to-grid and grid-to-ray interpolation approach for
        % time-of-flight-based image reconstruction approach
        raytogrid_interp_tof = 'Bspline';
        
        % get the number of receivers
        num_receiver = 2^12;
        
        % get the half size [m] of the computational grid along the x-y
        % plane
        half_grid_size = 13e-2;
        
        % get the grid spacng [m] for the k-Wave simulations
        grid_spacing_data = 1e-3;
        
        % get the Boolean controlling whether the acoustic absorption and dispersion are
        % included or not
        do_absorption = false;
        
        % get the number of linearised problems (outer iterations) for
        % reconstructing an image from the time of flights. Note athat each 
        % linearsied subproblem is solved using a 'sart' algorithm.
        num_iterout_tof = 4;
        
        % get the approach for interpolation of the pressure field from the
        % emitters to grid points and from grid points to receivers.
        % Using 'offgrid', the offgrid toolbox is used. (c.f. reference [5]
        % in the description of the script.)
        transducer_interp_approach = 'offgrid';
        
        % get the downsampling rate for the emitters
        emitter_downsampling_rate = 2;
        
        % get the downsampling rate for the receivers
        receiver_downsampling_rate = 2;
        
         % get the grid spacing [m] for image recosntruction
        grid_spacing_reconstruction = 2e-3;
end

% get the ratio of te number of receivers to the number of emitters
ratio_receiver_to_emitter = 4;

% get the approach for linearisation. This cab be 'absolute' or
% 'difference'. The latter is deprecated.
linearisation_approach = 'absolute';  %     'difference';  % 

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


if save_plots
    
    % get the directory for saving the plots
    plot_directory = 'simulation/results/plot_greens/';
    
    % make the directory for saving the plots, if it does not exist
    makeDirectory(plot_directory);
else
    
    % allocate an empty variable for the directory
    plot_directory = [];
    
end


if save_results
    
    % get the directory for saving the mat results
    results_directory = 'simulation/results/mat_greens/';
    
    % make the directory, if it does not exist
    makeDirectory(results_directory);
    
else
    
    % allocate an empty variable for the directory
    results_directory = [];
    
end



% get the optional inputs for simulating the pressure data sets using the
% k-Wave toolbox
data_args = {'num_worker_pool', num_worker_pool, 'grid_spacing_data', grid_spacing_data,...
    'grid_spacing_reconstruction', grid_spacing_reconstruction, 'cfl', cfl_number,...
    'transducer_interp_approach', transducer_interp_approach,...
    'single_emitter_receiver', [nan, emitter_index], 'emitter_downsampling_rate',...
    emitter_downsampling_rate, 'receiver_downsampling_rate',...
    receiver_downsampling_rate, 'time_downsampling_rate', time_downsampling_rate,...
    'do_calculate_tofs', do_calculate_tofs, 'noise_level', noise_level};


switch purpose
    case 'raytracing_validation'

        %%=================================================================
        % SIMULATE PRESSURE UST DATA SETS USING k-WAVE
        %==================================================================
        % if the puropose is validation of ray tracing, the clean data (without noise) is
        % used.
        [data, data_ref, ~, ~, emitter, receiver, kgrid,...
            medium, simulation_prop, data_simulation_time] =...
            simulateSettingData(do_data_sim, dim, scenario,...
            purpose, excit_pulse_name, half_grid_size, num_receiver,...
            ratio_receiver_to_emitter, do_absorption, oa_breast_path, machine_name,...
            res_path, local_res_path, data_args{:});
        
       if do_greens_approach
       %%==================================================================
       % COMPARE THE GREEN'S FUNCTION AND K-WAVE FOR HOMOGENEOUS
       % (ONLY-WATER) MEDIUM
       %===================================================================
       % get the optional inputs for homogeneous Green's function
       greens_hom_args = {'emitter_index', emitter_index, 'receiver_index', receiver_index,...
           'deconvolve_source', deconvolve_source, 'save_plots', save_plots};
       
       % get the relative discrepancies for the pressure field approximated
       % by the Green's function and k-Wave, on the receivers (or single receiver)
       % after being produced by a single emitter and propgating through only water
       [relative_discrepancy_signal_water_emitter, relative_discrepancy_signal_water_receiver] = ...
           validateGreensHomogeneous(data_ref, kgrid.t_array, emitter, receiver,...
           medium.sound_speed_ref, simulation_prop.f_max, plot_directory, greens_hom_args{:});

        %%=================================================================
       % COMPARE THE GREEN'S FUNCTION AND K-WAVE FOR HETEROGENEOUS
       % (BREAST-IN-WATER) MEDIUM
       %===================================================================
        % get the optional inputs for ray approaximtion to heterogeneous Green's function
        greens_het_args = {'num_worker_pool', num_worker_pool, 'emitter_index',...
            emitter_index, 'receiver_index', receiver_index,...
            'deconvolve_source', false, 'save_plots', save_plots, ...
            'auxiliary_method', auxiliary_method,...
            'attenuation_geom_method', attenuation_geom_method};
        
        
        % get the relative discrepancies for the pressure field approximated
        % by the Green's function and k-Wave, on the receivers after being produced by
        % a single emitter and propagating through the breast in water,
        % when the acoustic absorption (and dispersion) of the breast is
        % included or ignored.
        % Also, get the parameters on the sampled points along the forward and
        % adjoint (backprojection) rays.
        [ray_position, ray_time, ray_absorption, ray_position_left, ray_position_right,...
            adjoint_ray_position_left, adjoint_ray_position_right, ray_direction, ray_spacing,...
            rayspacing_receiver, relative_discrepancy_absorbing, relative_discrepancy_nonabsorbing] =...
            validateGreensHeterogeneous(data, data_ref, kgrid.t_array, emitter, receiver,...
            kgrid, medium, simulation_prop, [], plot_directory, greens_het_args{:});
        
        %%=================================================================
        % DISPLAY THE RAYS AFTER RAY TRACING VALIDATION
        %==================================================================
        % get the optional inputs
        display_args = {'save_plots', false, 'emitter_index', emitter_index, 'receiver_index', receiver_index};
        
        % display the forward ray initialised from a single emitter and
        % adjoint rays initilialised from a single receiver
        displayRays(emitter.positions, receiver.positions, ray_position, ray_time,...
            ray_absorption, ray_position_left, ray_position_right,...
            adjoint_ray_position_left, adjoint_ray_position_right, ray_direction, ray_spacing,...
            rayspacing_receiver, plot_directory, display_args{:});
        
       end
       
       
    case 'image_reconstruction'
        
        %%=================================================================
        % SIMULATE PRESSURE UST DATA SETS USING K-WAVE
        %==================================================================
        % if the purpose is image reconstruction, the noise-contaminated
        % data is used.
        [~, ~, data_noisy, data_ref_noisy, emitter, receiver, kgrid,...
            medium, simulation_prop, data_simulation_time] =...
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
            'matrix_construction_method', 'bent-ray',...
            'linear_subproblem_method', linear_subproblem_method,...
            'linearisation_approach', linearisation_approach,...
            'raytogrid_interp', raytogrid_interp_tof,...
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
            kgrid.t_array, emitter, receiver, medium.sound_speed_ref,...
            simulation_prop.z_offset, kgrid, medium.sound_speed, reconst_args_tof{:});
        
        if save_results
            
            % save the results in the given path
            save([results_directory 'results_' greens_optimisation_approach  '_' absorption_map...
                '_' num2str(noise_level) 'db'],...
                'img_tof', 'out_tof', 'para_tof', 'recon_grid', '-v7.3');
            
        end
        
        %%=================================================================
        % RECONSTRUCT THE SOUND SPEED IMAGE USING THE RAY APPROAXIMATION
        % TO HETEROGENEOUS GREENS FUNCTION ITERATIVELY
        %==================================================================
        if do_greens_approach
            
            % get the directories for the image reconstruction using the Green's approach
            directories_greens = [];
            directories_greens.results_directory = results_directory;
            directories_greens.images_directory = plot_directory;
            
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
                recon_grid, emitter, receiver, medium, img_tof, grid_spacing_reconstruction,...
                ray_initial_angles, kgrid, directories_greens, greens_args{:});
            
            if save_results
                
                % save the results in the given path
                save([results_directory 'results_' greens_optimisation_approach  '_' absorption_map...
                    '_' num2str(noise_level) 'db'],...
                    'img_greens', 'out_greens', 'para_greens', '-append');
                
            end
            
            
        end
        
end