% Example for image reconstruction of the speed from experimental data 
% measured in an experimental setting on a ring.
%..........................................................................
% 
%%=========================================================================
% EXPERIMENTAL SETTING
%==========================================================================
% The transducers are placed on a ring with radius of 13cm. The time step is
% 40ns, each signal starts at t = 80us after acoustic pulse emission. A single
% emitter placed on the ring is excited by an extric pulse and the induced acoustic 
% waves are measured in time by receivers placed on an arc encompassing 120 degrees of the
% same ring in front of the emitter. The emitter is rotated and is excited and
% the receivers are correspondingly rotated on the ring until the
% excitaions are performed for 2pi Rad.

%%=========================================================================
% LOAD OF UST DATA 
%==========================================================================
% The UST data used in this experiment can be loaded from ''www.zenodo.org''.
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
save_results = false;

% get the Boolean whether the plots are saved or not
save_plots = false;

%%=========================================================================
% THE PARAMETERS WHICH CAN BE CHANGED BY THE USER.
%==========================================================================

% get the approach for construction of the system matrix ('straight-ray' or
% 'bent-ray')
matrix_construction_method = 'bent-ray';

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
num_iterout_tof = 5;

% get the approach for linearisation. This can be 'absolute' or
% 'difference'. The latter is deprecated.
linearisation_approach = 'absolute';

% get the appproach for solving each linearised subproblem. Using a
% lineraisation using the 'absloute' approach, the arisng linearised
% subproblems can be solved using the 'sart' or 'conjuge_gradient'
% algorithm.
linear_subproblem_method = 'sart';


% get the step length
switch linear_subproblem_method
    case 'sart'
        step_length = 0.3;
    case 'conjugate_gradient'
        step_length = 0.1;
end


% An assumption of the maximum sound speed difference [m/s] with water in the medium 
% For each emitter-receiver pair,the first arrivals are are picked within a time
% window $t_s + [d/c_max, d/c_min]$, where $t_s$ is the shot (excitation) time,
% $d$ is the distance between emitter and receiver, and 
% c_min = sound_speed_water - sound_speed_time_window;
% c_max = sound_speed_water + sound_speed_time_window;
sound_speed_time_window = 80;

% A scalar representing the fraction of the maximum absolute amplitude of
% the signal as the right edge of the time window.
% If chosen nonzero, the right edge of the time window is chosen as 
tof_peak_frac_threshold = 0.35;

%%=========================================================================
% GET THE DIRECTORY FOR DATA AND STORING RESULTS
%==========================================================================

% get the data path
data_path = 'data/experimental_2D/experiment1/';

% make the directory, if it does not exist
makeDirectory(data_path);

% get the path for storing the time-of-flight data
tof_path = [data_path, 'experiment1' ];

if save_results
    
% get the path for the result
results_directory = 'results/experimental_2D/experiment1/';

% make the directory for the results
makeDirectory(results_directory)

end

if save_plots
    
    % get the path for storing the plots results
    plot_directory = 'plots/experiment1/';
    
    % make the directory for storing the plots
    makeDirectory(plot_directory)
    
end

%%=========================================================================
% THE PARAMETERS WHICH ARE ALWAYS FIXED AND ARE SPECIFIC TO THIS EXPERIMENT
%==========================================================================
    
% get the number of dimensions 
dim = 2;

% get the radius of the detection ring
detec_radius = 13000e-5;

% get the sound speed [m/s] in water
sound_speed_water = 1480;

%%=========================================================================
% GET EMITTERS
%==========================================================================

% get the angular poition of the emitters
% the emitters are rotated for matching the final reconstrucrted image 
% with the photograph of the phantom. The pi/6 angular offset was found
% empirically.
emitters_angular_position = pi/6 + deg2rad(0:2:358);

% get the number of emitters
num_emitter = size(emitters_angular_position, 2);

% get an estimated the excitation pulse duration
emitter.pulse_duration = 1.5e-6;

% get the shot time
emitter.shot_time = -80e-6;

% allocate a zero matrix for position of emitters
emitter.positions = zeros(dim, num_emitter);

% get the Cartesian position of emitters
[emitter.positions(1,:), emitter.positions(2,:)] = pol2cart(emitters_angular_position,...
    detec_radius);

% get the rotation indices
emitter.rotation_indices = 1: num_emitter;

%%=========================================================================
% GET RECEIVERS
%==========================================================================

% get the angular position of all receivers
receivers_angle = deg2rad(-88:2:90);

% get the angular position of receivers for all emitters
receivers_angular_position = emitters_angular_position' + pi + receivers_angle;

% get the number of receivers per emitter (excitation)
num_receiver = size(receivers_angular_position, 2);

% get the struct for the rotating position of receivers
receiver.positions = cell(1, num_emitter);

% get the rotating position of receivers
for ind_emitter = 1: num_emitter
    receiver.positions{ind_emitter} = zeros(dim, num_receiver);
    [receiver.positions{ind_emitter}(1,:), receiver.positions{ind_emitter}(2,:)] =...
        pol2cart(receivers_angular_position(ind_emitter,:), detec_radius);
end


%%=========================================================================
% GET THE MEASURED DATA
%==========================================================================

% get the data measured for object in water
load ([data_path, 'raw_data.mat']);
data = permute(S, [2, 1, 3]);
data = data(:, :, 1:num_emitter);
clear S;


% get the data measured for only water
load([data_path, 'datam90p90step250avg.mat']);
S(:, end)= [];

% match the exitations for only water with object-in-water data
data_ref = repmat(S, [1, 1, num_emitter]);

% expand the object-in-water data to match it with only-water data (include
% the side tranducers for only-wter data into the object-in-water data)
% For expansion, the time series are shifted by 1 so that the time
% difference won't be absolute zero, and therefore, the corresponding rays
% are included in SART. Because using SART, a more coverge for tranducers 
% increase the accuracy. Remind that the emitter-receiver pair with
% absolute zero time difference will not be included in the system matrix
% and SART.
data = cat(1, data_ref(1:15,[1,1:end-1],:),data, data_ref(16+size(data, 1):end,[1,1:end-1],:));
clear S


%%=========================================================================
% GET TIME ARRAY
%==========================================================================

% get the time spacing
dt = 4000e-11;

% get the number of time instants in the time array used for measuring data
num_time = size(data, 2);

% get the time array used for measuring data
time_array = cumsum([0, dt * ones(1, num_time-1)]);

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
    'sound_speed_time_window', sound_speed_time_window,...
    'tof_frac_peak_threshold', tof_peak_frac_threshold,...
    'tof_outliers', true,...
    'step_length', step_length,...
    'num_iterout', num_iterout_tof,...
    'emitter_downsampling_rate', 1,...
    'receiver_downsampling_rate', 1};

% compute the time-of-flights from ultrasound data sets for the
% breast-in-water and only-water data, and reconstruct the sound
% speed image iteratively from the time-of-flight data
[img_tof, recon_grid, ~, ray_initial_angles, out_tof, para_tof] = ...
    reconstructTimeofFlightImage(data, data_ref,...
    time_array, emitter, receiver, sound_speed_water,...
    [], tof_path, reconst_args_tof{:});



if save_results
    
    % save the results in the given path
    save([results_directory, 'results_', matrix_construction_method '_'...
        linear_subproblem_method '.mat'],...
        'img_tof', 'out_tof', 'para_tof', 'recon_grid', '-v7.3');
    
end

% match the reconstructed image with the photograph of the phantom by
% flipping the columns of the reconstructed image
img_tof = flip(img_tof);

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

% display the sound speed values on the diagonal of the reconstructed image
h2 = figure; plot(diag(img_tof)); xlabel('Index of grid point');
ylabel('c [ms^{-1}]'); axis image;


% display the sound speed values on the diagonal of the transpose of the 
% reconstructed image
h3 = figure; plot(diag(img_tof')); xlabel('Index of grid point');
ylabel('c [ms^{-1}]'); axis image;




