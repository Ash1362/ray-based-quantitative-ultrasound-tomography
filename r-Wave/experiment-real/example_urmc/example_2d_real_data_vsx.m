% Example for image reconstruction of the speed from the in-vitro transmission
% ultrasound datasets measured on a ring. The results obtained
% from running this example script for the settings:
% "emitter_downsampling_rate = 1" and "receiver_downsampling_rate = 1"
% are reported in the reference [3]. 
% Two datasets measured from two slices of the assocaited phantom is freely
% available and should be downloaded from the link provided in the reference [5].
% 
%%=========================================================================
% EXPERIMENTAL SETTING
%==========================================================================
% The transducers are placed on a ring with a radius of approximately 11cm. 
% For further details, See: Section ''IV. RESULTS, C. In-Vitro Experiments''
% of reference [4].
%%=========================================================================
% IMAGE RECONSTRUCTION
%==========================================================================
% This script first loads the in-vitro dataset from the data_path, which is
% set 'data/experiment_rochester/' be default. The user must upload the
% dataset from the links provided in [5], and add these datasets:
% 'VSX_YezitronixPhantom1' and 'VSX_YezitronixPhantom2' to the data path
% before running this script. These datasets correspond to two slices of
% the associated in-vitro phantom. The script loads one of these data sets 
% depending on the setting experment_number 1 or 2 in the "GET THE TASK" Section
% of this manuscript.
% By running this script, a low-resolution and low-contrast image of the
% sound speed is first reconstructed from only three iterations (linearizations)
% of the time-of-flight-based algorithm.
% The time-of-flight-based reconstructed sound speed image is then used as 
% the initial guess for the main image reconstruction algorithm, which is 
% the Hessian-free ray-Born inversion. 
% The users are referred to the reference [1] for details of the
% Hessian-free ray-Born inversion, implemented by this script as main image
% reconstruction algorithm.
%%=========================================================================
% USER INPUT CHOICE
%==========================================================================
% By default, both emitter_downsampling_rate and receiver_downsampling_rate
% are chosen 2 in this script, i.e., the data for 128 emitters and
% 128 receivers are used for image reconstruction.
% However, the results in the associated paper are based on setting both
% sampling rates 1, i.e., the data for all 256 emitters and 256 receivers
% were used for image reconstruction.
% The user's choice is between 1 or 2. (A higher value leads to less
% computational cost, but lower accuracy.)
%%=========================================================================
% REFERENCES
%==========================================================================
% If you find the toolbox useful for your research, please consider citing these papers:...
% 1 - A. Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution
% quantitative ultrasound tomography❞, 2022, https://arxiv.org/abs/2211.00316.
% 2 - A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering
% for biomedical ultrasound tomography❞, Inverse Problems vol. 37, no.11, 115003, 2021.
% 3- A. Javaherian, ❝The first in vitro and in vivo validation of the hessian-free
% ray-Born inversion for quantitative ultarsound tomography❞, 2025.
% 4 - R. Ali et al., ❝2-D Slicewise Waveform Inversion of Sound Speed and
% Acoustic Attenuation for Ring Array Ultrasound Tomography Based on a 
% Block LU Solver,❞ in IEEE Transactions on Medical Imaging, 
% vol. 43, no. 8, pp. 2988-3000, Aug. 2024.

% The following data-base/toolbox were used in this project.
% 5 - R. Ali, GitHub link: https://github.com/rehmanali1994/WaveformInversionUST
% (DOI: https://zenodo.org/badge/latestdoi/684631232)

% You may also think about citing this preprint:
% 6- A. Javaherian, ❝Full-waveform Approximation of Finite-Sized Acoustic Apertures:
% Forward and Adjoint Wavefields❞, https://arxiv.org/abs/2212.04466.

% author: Ashkan Javaherian
% date:            - 07.02.2024
% last update:     - 15.07.2025
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

% get the experiment number (should be set 1 or 2 for downloading the
% 'VSX_YezitronixPhantom1' and 'VSX_YezitronixPhantom2' datasets, respectively.)
experiment_number = 2; % 1,2


%%=========================================================================
% THE USER CHOICE FOR EMITTER/RECEIVER DOWNSAMPLING RATE (1,2, OR 4)
%==========================================================================
% By default, both emitter_downsampling_rate and receiver_downsampling_rate
% are chosen 2 in this script, i.e., the data for 128 emitters and
% 128 receivers are used for image reconstruction.
% However, the results in the associated paper [3] are based on choosing both
% sampling rates 1, i.e., the data for all 256 emitters and 256 receivers
% were used for image reconstruction.
% The user's choice is between 1 or 2. (Higher value leads to less
% computational cost, but lower accuracy.)

% get the downsampling rate for the emitters (1 or 2)
emitter_downsampling_rate = 2;

% get the downsampling rate for the receivers (1 or 2)
receiver_downsampling_rate = 2;

%%=========================================================================
% GET THE DIRECTORY FOR DATA AND STORING RESULTS
%==========================================================================

% display the experiment number 
disp(['The experiment number is:' num2str(experiment_number)])

% get the data path
data_path = 'data/experiment_rochester/';

% make the directory, if it does not exist
makeDirectory(data_path);

if save_results

    % get the path for the result
    results_directory = ['results/experiment_rochester/experiments/'...
    'experiment' num2str(experiment_number) '/'];

    % make the directory for the results
    makeDirectory(results_directory)
else

    % allocate an empty variable
    results_directory = [];

end

%%=========================================================================
% THE PARAMETERS WHICH ARE SET FOR THIS SPECIFIC DATASET
%==========================================================================
% the type of the absorption coefficient map used for image reconstruction
% using the ray approaximation to heterogeneous Green's function. This can
% be 'true', 'homogeneous', or 'none'
% Please read section Numerical results in [1]!
absorption_map = 'none';

% An assumption of the maximum sound speed difference [m/s] with water in the medium
% For each emitter-receiver pair,the first arrivals are are picked within a time
% window $t_s + [d/c_max, d/c_min]$, where $t_s$ is the shot (excitation) time,
% $d$ is the distance between emitter and receiver, and
% c_min = sound_speed_water - sound_speed_time_window;
% c_max = sound_speed_water + sound_speed_time_window;
sound_speed_time_window = 40;

%%=========================================================================
% THE PARAMETERS WHICH ARE ALWAYS FIXED
%==========================================================================
% Boolean controlling whether the Green's approach is used, or only the
% image reconstructed using the time-of-flight-based approach is
% reconstructed.
do_greens_approach = true;

% get the approach for construction of the system matrix ('straight-ray' or
% 'bent-ray')
matrix_construction_method = 'bent-ray';

% get the grid spacing [m] for image reconstruction, 1e-3 or 2e-3.
% (Default: 1e-3)
% The grid spacing 1e-3 m is used in the papers, and is recommended.
grid_spacing_reconstruction = 1e-3; 

% get the smoothing window size for the ToF-based approach
smoothing_window_size_tof = 7;

% get the numer of workers for parallel programming
% the user may want to change the number of workers.
num_worker_pool = 12;

% set the approach used for ray tracing for image reconstruction using
% time-of-flight data. 
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
num_iterout_tof = 3; 

% get the approach for linearisation. This can be 'absolute' or
% 'difference'. The latter is deprecated.
linearisation_approach = 'absolute';

% get the appproach for solving each linearised subproblem. Using a
% linearisation using the 'absloute' approach, the arisng linearised
% subproblems can be solved using the 'sart' or 'conjugate_gradient'
% algorithm.
linear_subproblem_method = 'sart';

% get the step length for the time-of-flight-based image reconstruction.
% This is a scalar value which is multiplied by the update direction of the
% sound speed after solving each linearized subproblem. % Default: 1e-1 
step_length_tof = 1e-1; 

% get the step length for the Hessian-free ray-Born image reconstruction
% approach. This is a scalar value which is multiplied by the update direction
% of the sound speed after solving each linearized subproblem. % Default: 1.2e-1
step_length_greens = 1.2e-1;  

% get the relative radius of the binary mask for ray tracing
mask_coeff = 1.02; 

% A scalar representing the fraction of the maximum absolute amplitude of
% the signal as the right edge of the time window.
% If chosen nonzero, the right edge of the time window is chosen as
tof_peak_frac_threshold = 0.25;    

% get the downsampling rate for time
time_downsampling_rate = 1;

% get the optimisation approach for image reconstruction
% using the ray approaximation to heterogeneous Green's function
% This can be set 'hessian', or 'backprojection'.
% 'backprojection' is based on the Green's approach proposed in [1].
% 'hessian' is based on the Green's approach proposed in [2].
greens_optimisation_approach = 'backprojection';

% get the Boolean controlling whether the pressure source is deconvolved
% from the pressure time series or not. If purpose = 'image_reconstruction',
% we set deconvolve_source = true.
deconvolve_source = true;

% the minimum angular distance [Rad] between emitter-receiver pairs included in
% the image reconstruction. This is a coefficient in Rad and its
% multiplication by the radius [m] of the detection ring gives the permissible
% minimum distance for emitter-receiver pairs for including in the image
% reconstruction.
minimum_distance_angular = pi/4; 

% the fraction of the emitter-receiver pairs excluded from image reconstruction
% for a specific frequency, because their magnitudes are outliers.
% (default: 0.01)
outliers_fraction = 1e-2;

% Boolean controlling whether an exponential time window normalised (with maximum)
% at the arrival time of each signal is applied on each measured time trace or not.
multiply_data_window = true;


% The approch taken for projecting the real positions of transducers to the
% ring. (The ray-based image reconstruction approaches assume the transducers
% are positioned on the ring upto machine precision.)
% This can be set 'equidistant-angle' or 'real-angle'
project_transducer_approach = 'real-angle';


% the approach for computing the geometrical attenuation, auxiliary: using
% auxiliary rays, or 'raylinked' (linked rays without using auxiliary rays)
attenuation_geom_method = 'auxiliary';

% the approach for tracing auxiliary rays, which will be used for computing
% the geomterical attenuation, 'paraxial' or 'angle_perturbation'.
auxiliary_method = 'paraxial';

% get the number of dimensions
dim = 2;

% get the sound speed [m/s] in water
sound_speed_water = 1490;

%==========================================================================
% LOAD TIME ARRAY, TRANSDUCERS' POSITIONS AND DATA SETS
%==========================================================================
% get the file name       
filename = ['VSX_YezitronixPhantom' num2str(experiment_number)];

% load the data set itself and information about the experiment setting, i.e.,
% transducers' positions and time array
load([data_path, filename, '.mat'], ...
    'time', 'transducerPositionsXY', 'full_dataset');
full_dataset = single(full_dataset);


%%=========================================================================
% GET EMITTERS
%==========================================================================

% get an estimated the excitation pulse duration
emitter.pulse_duration = 4e-6;

% get the shot time
emitter.shot_time = 0;

% set the emitter pulse a Diract delta function at the origin with a
% normalised amplitude. The real amplitude of the source will then affect
% the scaling factor of the data set.
source_amplitude = 1;

% get the emitter indices included in the image reconstruction
emitter_indices = 1:emitter_downsampling_rate:size(transducerPositionsXY, 2);

% get the real position of emitters
emitter.positions_real = transducerPositionsXY(:, emitter_indices);

% get the number of emitters
num_emitter = length(emitter_indices);

% set the radius of the detection ring the mean distrance of the real
% position of emitters to the centre of the detection ring
detec_radius = mean(vecnorm(emitter.positions_real));

% get the real angular position of the transducers
[emitter_angles_real, ~] = cart2pol(emitter.positions_real(1,:),...
    emitter.positions_real(2,:));

switch project_transducer_approach

    case 'equidistant-angle'

        % get equidistant emitters' angles for ray tracing and image reconstruction
        emitter_angles = wrapToPi(emitter_angles_real(1): (-2*pi/num_emitter):...
            (emitter_angles_real(1)-2*pi/num_emitter*(num_emitter-1)));

    case 'real-angle'

        % keep the real angles of the transducers and project the positions
        % onto the ring aligned with radius coordinate
        emitter_angles = emitter_angles_real;


end

% project the emitters' real positions to the periphery of a circle with the mean radius
% of all real positions and equidistant angular positions on the periphery
% of that circle
[emitter.positions(1,:), emitter.positions(2,:)] = pol2cart(emitter_angles, detec_radius);


%[emitter.positions(1,:), emitter.positions(2,:)] = pol2cart(...
%   cart2pol(emitter.positions_real(1,:), emitter.positions_real(2,:)), detec_radius);




%%=========================================================================
% GET RECEIVERS
%==========================================================================

% get the receiver indices included in the image reconstruction
receiver_indices = (floor(receiver_downsampling_rate/2) + 1) :...
    receiver_downsampling_rate:size(transducerPositionsXY, 2);

% get the number of receivers
num_receiver = length(receiver_indices);

% get the real position of receivers
receiver.positions_real = transducerPositionsXY(:, receiver_indices);


% get the real angular position of the receivers
[receiver_angles_real, ~] = cart2pol(receiver.positions_real(1,:),...
    receiver.positions_real(2,:));

switch project_transducer_approach

    case 'equidistant-angle'

        % get equidistant emitters' angles for ray tracing and image reconstruction
        receiver_angles = wrapToPi(receiver_angles_real(1): (-2*pi/num_receiver):...
            (receiver_angles_real(1)-2*pi/num_receiver*(num_receiver-1)));

    case 'real-angle'

        % keep the real angles of the transducers and project the positions
        % onto the ring aligned with radius coordinate
        receiver_angles = receiver_angles_real;

end

% project the receivers' real positions to the periphery of a circle with
% the mean radius of all receivers and equidistant angular positions
% on the pheriphery of that circle
[receiver.positions(1,:), receiver.positions(2,:)] = pol2cart(receiver_angles, detec_radius);


%[receiver.positions(1,:), receiver.positions(2,:)] = pol2cart(...
%    cart2pol(receiver.positions_real(1,:), receiver.positions_real(2,:)), detec_radius);



%%=========================================================================
% PLOT THE DISPLACEMENTS OF TRANSDUCERS' POSITIONS IN THE CARTESIAN AND POLAR
% COORDINATES
%==========================================================================

% display the plot for the difference of polar position of real and
% projected emitters
figure; plot(1:num_emitter, emitter_angles-emitter_angles_real);
xlabel('Emitter index'); ylabel('Error in angle of emitters [m]');

% display the plot for the distance of real and projected emitters
figure; plot(1:num_emitter, vecnorm(emitter.positions-emitter.positions_real));
xlabel('Emitter index'); ylabel('Error in position of emitters [m]');

% display the plot for the difference of polar position of real and
% projected receivers
figure; plot(1:num_receiver, receiver_angles-receiver_angles_real);
xlabel('Receiver index'); ylabel('Error in angle of receivers [m]');

% display the plot for the distance of real and projected receiver
figure; plot(1:num_receiver, vecnorm(receiver.positions-receiver.positions_real));
xlabel('Receiver index'); ylabel('Error in position of receivers [m]');

%%=========================================================================
% GET THE MEASURED DATA
%==========================================================================

% get the data measured for object in water
full_dataset = permute(full_dataset, [2, 1, 3]);

% sample the data for receivers, time and emitters
data = full_dataset(receiver_indices, 1:time_downsampling_rate:end, emitter_indices);

% clear the original data set
clear full_dataset

%%=========================================================================
% GET TIME ARRAY
%==========================================================================

% get the time array used for measuring data
time_array = time(1:time_downsampling_rate:end);

% get the number of time instants in the time array used for measuring data
num_time = size(data, 2);

if ~(length(time_array)-num_time)
    disp('The time array in data and time array itself are consistent.')
end


%==========================================================================
% COMPUTE THE ONLY-WATER TOFS FROM REAL DISTANCES OF EMITTER-RECEIVER PAIRS
%==========================================================================
% compute the only-water TOFs [s] for the rael position of emitter-receiver
% pairs
tofs_water_real = calculateDistanceEmitterReceiver(emitter.positions_real,...
    receiver.positions_real, [])/sound_speed_water;



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
    'raytogrid_spacing', 1,...
    'gridtoray_interp', gridtoray_interp_tof,...
    'grid_spacing', grid_spacing_reconstruction,...
    'smoothing_window_size', smoothing_window_size_tof,...
    'binaries_emitter_receiver', 'distances',...
    'minimum_distance_coeff', 2 * minimum_distance_angular,...
    'mask_coeff', mask_coeff,...
    'sound_speed_time_window', sound_speed_time_window,...
    'tof_frac_peak_threshold', tof_peak_frac_threshold,...
    'tof_outliers', false,...
    'step_length', step_length_tof,...
    'num_iterout', num_iterout_tof,...
    'emitter_downsampling_rate', 1,...
    'receiver_downsampling_rate', 1};




% compute the time-of-flights from ultrasound data sets for the
% breast-in-water and only-water data, and reconstruct the sound
% speed image iteratively from the time-of-flight data
[img_tof, recon_grid, ~, ray_initial_angles, out_tof, para_tof] = ...
    reconstructTimeofFlightImage(data, tofs_water_real,...
    time_array, emitter, receiver, sound_speed_water,...
    [], [], reconst_args_tof{:});



if save_results

    % save the results in the given path
    save([results_directory, 'results_' greens_optimisation_approach, ...
        '_e' num2str(emitter_downsampling_rate)], 'data', 'time_array', 'recon_grid',...
        'emitter', 'receiver', 'sound_speed_water', 'img_tof', 'out_tof',...
        'para_tof', 'recon_grid', 'results_directory', 'greens_optimisation_approach',...
        'emitter_downsampling_rate', '-v7.3');

end


% display the reconstructed image
h1 = figure; imagesc(img_tof.'); axis image; axis off; colorbar;

% display the sound speed values on the diagonal of the reconstructed image
h2 = figure; plot(diag(img_tof)); xlabel('Index of grid point');
ylabel('c [ms^{-1}]'); axis image;


%%=========================================================================
% GET THE EMITTING PULSE
%==========================================================================

% allocate a zero vector for emitter pulse
% Using experimental data, the pulse should be given to the struct emitter
% after time-of-flight picking step is completed.
% get the emitter pulse as a Dirac delta function at time origin
emitter.pulse = [1, zeros(1, length(time_array)-1)];

%%=================================================================
% RECONSTRUCT THE SOUND SPEED IMAGE USING THE RAY APPROXIMATION
% TO HETEROGENEOUS GREENS FUNCTION ITERATIVELY
%=================================================================
if do_greens_approach

    % get the directories for storing the results associated with the image
    % reconstruction using the Green's approach
    directories_greens.results_directory = results_directory;

    % get the directories for storing the plots associated with the image
    % reconstruction using the Green's approach
    directories_greens.image_directory = [];


    % initialize the frqeuency levels by 1
    frequency_level = 1;

    % get the optional inputs
    greens_args = {'num_worker_pool', num_worker_pool,...
        'optimisation_approach', greens_optimisation_approach,...
        'max_raylinking_iter', 500,...
        'absorption_map', absorption_map,...
        'deconvolve_source', deconvolve_source,...
        'do_scaling', true,...
        'attenuation_geom_method', attenuation_geom_method,...
        'auxiliary_method', auxiliary_method,...
        'num_frequency', 140,... % Default:140
        'num_frequency_level', 70,... % Default:70
        'num_frequency_per_level', 2,...
        'step_length', step_length_greens,...
        'f_range', [0.35e6, 0.95e6],...% [0.35e6, 0.95e6]
        'mask_coeff', mask_coeff-0.01,...
        'minimum_distance_coeff', minimum_distance_angular,... % Default:pi/4
        'outliers_fraction', outliers_fraction,...
        'multiply_data_window', multiply_data_window,...
        'cut_off', 0.75,...
        'order', inf,...
        'frequency_level', frequency_level:frequency_level + 6};


    % initialize the output struct for all iterations
    out_greens_tot = [];
    out_greens_tot.gradient_obj = 0;
    out_greens_tot.ray_cpu_time = 0;
    out_greens_tot.update_norm = 0;
    out_greens_tot.gradient_cpu_time = 0;
    out_greens_tot.objective_function = 0;
    out_greens_tot.ray_initial_angles = ray_initial_angles;

    % get the initial guess for images reconstructed using the Green's
    % approach
    img_greens = img_tof;

    for levels = 1:10

        % compute the sound speed image using the ray approximation to
        % the heterogeneous Green's function. The image is
        % reconstructed iteratively from low to high frequencies.
        [img_greens, out_greens, para_greens] = reconstructGreensImage(data,...
            time_array, recon_grid, emitter, receiver, sound_speed_water,...
            img_greens, out_greens_tot.ray_initial_angles, [], greens_args{:});


        % update the output struct for all iterations
        out_greens_tot.gradient_obj = out_greens_tot.gradient_obj + out_greens.gradient_obj;
        out_greens_tot.ray_cpu_time = out_greens_tot.ray_cpu_time + out_greens.ray_cpu_time;
        out_greens_tot.update_norm = out_greens_tot.update_norm + out_greens.update_norm;
        out_greens_tot.gradient_cpu_time = out_greens_tot.gradient_cpu_time + out_greens.gradient_cpu_time;
        out_greens_tot.objective_function = out_greens_tot.objective_function + out_greens.objective_function;
        out_greens_tot.ray_initial_angles = out_greens.ray_initial_angles;


        % update the frequency levels
        greens_args{end} = greens_args{end} + 7;

        if save_results

            % save the results in the given path
            save([results_directory, 'results_' greens_optimisation_approach,...
                '_e' num2str(emitter_downsampling_rate) num2str(levels)], 'img_greens', 'out_greens_tot',...
                'para_greens', 'directories_greens', 'greens_args', 'emitter', '-v7.3');

        end


    end


end



