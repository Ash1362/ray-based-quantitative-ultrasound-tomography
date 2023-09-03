% Example for image reconstruction of the speed from experimental data 
% measured on a ring.
%..........................................................................
% 
%%=========================================================================
% EXPERIMENTAL SETTING
%==========================================================================
% This example reconstructs a low-resolution sound speed image from
% time-of-flight data obtained from ultrasound time series measured in 
% an experimental setting. In this setting, 16 line arrays are placed on the
% periphery of a crircle, and 16 transducers are placed on each line array.
% In contrast with 'example_2d_tof_real_data_circle.m', for each emitter,
% the ray linking for each receiver array is done seprately using its specific
% line equation.
% The experimental setting is described in:
% ''https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=10139828&tag=1''
%%=========================================================================
% LOAD OF UST DATA 
%==========================================================================
% The time-of-flight data used in this experiment can be loaded from 
% ''https://github.com/ucl-bug/ust-sart''
%
%
% author: Ashkan Javaherian
% date:            - 05.02.2023
% last update:     - 15.07.2023
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian

clear all
close all
clc

% run startup
startup_simulation_ust;

%%=========================================================================
% GET THE TASK
%==========================================================================

% get the Boolean controlling whether the results are saved or not
save_results = true;


% get the Boolean whether the plots are saved or not
save_plots = true;

%%=========================================================================
% THE PARAMETERS WHICH CAN BE CHANGED BY THE USER.
%%=========================================================================

% get the approach for construction of the system matrix ('straight-ray' or
% 'bent-ray')
matrix_construction_method = 'bent-ray';

% get the appproach for solving each linearised subproblem. Using a
% lineraisation using the 'absloute' approach, the arising linearised
% subproblems can be solved using the 'sart' or 'conjuge_gradient'
% algorithm.
linear_subproblem_method = 'sart';

% get the grid spacing [m] for image reconstruction, 1e-3 or 2e-3.
% (Default: 1e-3)
% The grid spacing 1e-3 m is used in the papers, and is recommended.
% If your access to memory is limited, you can also use 2e-3 m.
grid_spacing_reconstruction = 2e-3;

% get the numer of workers for parallel programming
% the user may want to change the number of workers.
num_worker_pool = 16;

% choose the approach used for ray tracing for image reconstruction using
% time-of-flight data. This can be set: 'Mixed-step', 'Dual-update',
% 'Characteristics', or 'Runge-kutta-2nd'.(Default: 'Runge-kutta-2nd')
raytracing_method_tof = 'Runge-kutta-2nd';

% get the ray-to-grid and grid-to-ray interpolation approach for
% time-of-flight-based image reconstruction approach
gridtoray_interp_tof = 'Bspline';

% get the number of linearised problems (outer iterations) for
% reconstructing an image from the time of flights. (Note that the number
% of iterations for reconstructing the optimal image is 12-15, but for
% providing an initial guess for the Green's approach, only few iterations
% of the 'sart' algorithm (simulateneous algebraic reconstruction algorithm)
% is sufficient for providing an image accurate enough and with minimal artefact.

switch linear_subproblem_method
    case 'conjugate_gradient'
        
        % get the number of outer iterations
        num_iterout_tof = 5;
        
        % get the step length
        step_length = 0.1;
        
        
    case 'sart'
        
        % get the number of outer iterations
        num_iterout_tof = 10;
        
        % get the step length
        step_length = 0.5;
end


% get the approach for linearisation. This can be 'absolute' or
% 'difference'. The latter is deprecated.
linearisation_approach = 'absolute';

%%=========================================================================
% GET THE DIRECTORY FOR DATA AND STORING RESULTS
%==========================================================================

% get the path for the data
data_directory = 'data/experimental_2D/experiment2/';

% load the position of the transducers 
load([data_directory, 'element_positions.mat']);

% get the path in which the tof sinogram has been stored
load([data_directory, 'delta_tof.mat']);

if save_results
    
% get the path for the result
results_directory = 'results/experimental_2D/experiment2/';

% make the directory for the results
makeDirectory(results_directory)

end

if save_plots
    
    % get the path for storing the plots results
    plot_directory = 'plots/experiment2/';
    
    % make the directory for storing the plots
    makeDirectory(plot_directory)
    
end

%%=========================================================================
% THE PARAMETERS WHICH ARE ALWAYS FIXED AND ARE SPECIFIC TO THIS EXPERIMENT
%==========================================================================
    
% get the number of dimensions 
dim = 2 ;

% get the number of linear arrays
num_array = 16;

% get the number of transducers at each linear array
num_transducer_per_array = 16;

% get the number of emitters 
num_excitation = num_array * num_transducer_per_array;

% get the sound speed [m/s] in water
sound_speed_water = 1480;


%%=========================================================================
% GET POSITION OF EMITTERS AND RECEIVERS
%==========================================================================
% get the matrix containg the position of emitters

% get the rotation indices as {1,...,num_array, 1,...., num_array,...
% 1,...,num_array}, which iterates num_array * num_transducer_per_array
% times
emitter.rotation_indices = vectorise(repmat((1: num_array).',...
   [1, num_array * num_transducer_per_array]));

% get the number of emitters 
num_emitter = length(emitter.rotation_indices);

% allocate a zero matrix for position of emitters
emitter.positions = zeros(dim, num_emitter);

% get the x Cartesian coordinate of the transducers via duplicating each emitter
% position num_array times. The reason is that the code reads each excition
% (emitter) as each pair (emitter element, receiver_array). Therefore, the number
% of emitters will be the number of transducer elements times number of arrays, in 
% which the number of transducers elements is number of arrays times the number
% of transducer per arrays.
emitter.positions(1, :) = vectorise(repmat(element_positions(:, 1).',...
    [num_array, 1])).';

% get the y Cartesian coordinate of the transducers in the same way
emitter.positions(2, :) = vectorise(repmat(element_positions(:, 2).',...
    [num_array, 1])).';

% get the cell array containing the position of receivers
receiver.positions = mat2cell(element_positions.', dim,...
    num_transducer_per_array * ones (1, num_array) );


%%=========================================================================
% GET THE DIFFERENCE TIME-OF-FLIGHT SINOGRAM
%==========================================================================
% make the Nans in the tof sinogram zero
delta_tof(~isfinite(delta_tof)) = 0;

% get the discrepancy of time-of-flight data
% make the lower triangle the same as upper triangle
tof_data = reshape((delta_tof + delta_tof').', [num_transducer_per_array, num_emitter]); 

figure; imagesc((delta_tof + delta_tof').'); colorbar;
figure; imagesc(tof_data); colorbar;


% get the path and file for the adapted time-of-flight sinogram
tof_path = [data_directory, 'data'];

% save the adapted tof sinogram.
% By convention, the name of mat file for the time-of-flight data must be
% ended by '_tof_sinogram'
save([tof_path, '_tof_sinogram.mat'], 'tof_data', '-v7.3')

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
    'step_length', step_length,...
    'grid_spacing', grid_spacing_reconstruction,...
    'binaries_emitter_receiver', 'distances',...
    'num_iterout', num_iterout_tof,...
    'minimum_distance_coeff', sqrt(2),...
    'emitter_downsampling_rate', 1,...
    'receiver_downsampling_rate', 1};

% compute the time-of-flights from ultrasound data sets for the
% breast-in-water and only-water data, and reconstruct the sound
% speed image iteratively from the time-of-flight data
[img_tof, recon_grid, ~, ray_initial_angles, out_tof, para_tof] = ...
    reconstructTimeofFlightImage([], [],...
    [], emitter, receiver, sound_speed_water,...
    [], tof_path, reconst_args_tof{:});


if save_results
    
    % save the results in the given path
    save([results_directory, 'results_', matrix_construction_method '_'...
        linear_subproblem_method '.mat'],...
        'img_tof', 'out_tof', 'para_tof', 'recon_grid', '-v7.3');
    
end


% display the reconstructed image
h1 = figure; imagesc(img_tof); axis image; axis off; colorbar;

if save_plots
    
    % store the figure of the reconstructed image in different formats
    switch matrix_construction_method
        
        % get the figure name
        case 'straight-ray'
            
            switch linear_subproblem_method
                case 'conjugate_gradient'
                    fig_name = 'fig_1b';
                case 'sart'
                    fig_name = 'fig_1c';
            end
            
        case 'bent-ray'
            
            switch linear_subproblem_method
                case 'conjugate_gradient'
                    fig_name = 'fig_1d';
                case 'sart'
                    fig_name = 'fig_1e';
            end
            
    end
    
    saveas(h1, [plot_directory, fig_name '.fig']);
    saveas(h1, [plot_directory, fig_name '.png']);
    saveas(h1, [plot_directory, fig_name '.tiff']);
    saveas(h1, [plot_directory, fig_name '.eps'], 'epsc');
end


