function [img, recon_grid, data_ids, out, para] = reconstructSoundspeedImagePammoth(...
     data_ids, results_directory, varargin)
%RECONSTRUCTSOUNDSPEEDIMAGEPAMMOTH reconstructs the sound speed image from
%the ultrasound portion of the Pammoth data.
%   
%
% DESCRIPTION:
%             reconstructSoundspeedImagePammoth loads the ultrasound
%             portion of the hybrid photoacoustic-ultrasound data collected
%             from the Pammoth system and computes the time-of-flight data,
%             or alternatively loads the already computed time-of-flights
%             from the ultrasound data. The sound speed image is then
%             reconstructed iteratively from the time-of-flight data by
%             iteratively linearising the objective function, which is the
%             misfit of the ray-based modelled and computed time-of-flights
%             in an L2 norm sense.
%            
%
% USAGE:
%     
%
% INPUTS:
%       data_ids          - a struct array containing the IDs (and information)
%                           about the data. This struct includes the fields:
%       'main_path'       - the main directory, which must be fixed for all 
%                           UST data
%       'object'            the type of object for which the ultrasound data is measured.
%                           This can be 'volunteers', 'patients', 'us_measurements',
%                           'breast_phantom', 'channel_phantom', or
%                           'water_phantom'.
%       'id'                the ID of the data
%       'date'              the date on which the data were measured
%       'measurement_number' the measurement number as a string
%       'cup_size'           the cup size, which can be integers 2,3,...,8.
%       'main_path_user_data' - the main path for storing the results, which can
%                            be set by the user. It is recommended that the
%                            default path given in the example is used.
%       results_directory -  a directory for saving the results on the local
%                            machine. If this input is given as an empty
%                            variable, it will be automatically set 
%                            ['SoSreconstructions/' data_ids.id]]
%  
%       
% OPTIONAL INPUTS:
%      'num_worker_pool'    - the number of workers for parallel programming
%                             on CPUs
%      'do_calculate_tofs'  - Boolean controlling whether the time-of-flights
%                             (TOFs) are calculated, or the already calculated
%                             TOFs are loaded. This Boolean will be automatically
%                             set true, if the TOFs do not exist in the
%                             defined path.
%     'matrix_construction_method' - the approach for reconstruction of the
%                             sound speed. This can be set 'straight-ray',
%                             or 'bent-ray'.
%      'linear_subproblem_method' - the method for solving the linearised
%                           subproblems, which can be either 'sart',
%                           'conjugate_gradient', or 'steepest_descent'.
%                           (Default : 'sart')
%     'linearisation_approach' - the approach for linearisation of the
%                             misfit function. This can be either 'difference',
%                             or 'absolute'. (default: 'absolute')
%     'choose_last_us_position' - A parameter determining how to choose the last
%                                 US position, which can be 'early' (chosen by the user),
%                                 or automatically set all the angular US positions - 1 ('all'). 
%                                 The choice for the latter is based on the fact that 
%                                 the last position is the same as the first position,
%                                 so the last position is removed. (Default = 'all')
%      'mask_coeff'        - a coefficient which defines the normalised radius of
%                            the circular binary mask with radius = mask_coef * detec_radius
%                            for image reconstruction. Here, detec_radius
%                            is the radius of the detection surface, which
%                            must be chosen greater 1. (Default : 1.03)
%      'do_remove_bottom_emitters' - Boolean controlling whether the difference
%                               time-of-flights associated with the emitters 
%                               at the bottom of the cup are nulled. 
%                               All the rays which are intialised from those
%                               emitters and intercepted by all the receivers
%                               travel only across the water.
%     'transducer_height_position' - the nimimum permissible z position [m] for
%                                the emitters
%     'z_pos_height'        - the maximum z position [m] of the grid points, if
%                             the reconstruction geometry is set 'real'
%                             (Default: 1e-2)
%     'gridtoray_interp'    - the method for interpolation of the refractive
%                             index from the grid points to arbitrary off-grid
%                             points and vice versa. This can be set
%                             'Bilinear', or 'Bspline'.
%     'reconstruct_image'   - Boolean controlling whether the image is
%                             reconstructed or not. If not, only the
%                             time-of-flights are computed and saved.
%     'save_tofs'           - Boolean controlling whether the computed 
%                             time-of-flights are saved, or not.
%     'save_results'        - Boolean controlling whether the results are
%                             saved, or not.
%     'display_transducers' - Boolean controlling whether the transducers during
%                             the rotations are visualised or not (Default = false)
%                          
%      
%       
%
% OUTPUTS:
%      img                    - the reconstructed image of the sound speed
%      recon_grid             - the grid for image reconstruction
%      data_ids               - a struct array containing the IDs (and information)
%                               about the data. Compared to the same
%                               data_ids given as an input to the function,
%                               this struct array includes additional
%                               field:
%      'only_water_measurement_date' - the date on which the only-water data
%                               is measured.
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
%      'imgs'                 - all the reconstructed images of the sound speed
%                               after solving the linearised subproblems
%                               If the sound_speed_ground_truth is given, 
%                               this struct can also includes the fields:
%      'relative_error'       - the L2 norm of the discrepancy of the
%                               updates of the sound speed and the ground
%                               truth over the L2 norm of the disrepancy of
%                               water and the ground truth time 100.
%      'rmse'                 - the rmse of the reconstructed sound speed images
%                               with respect to the ground truth.
%       para                  - all the information about the approaches
%                               used for time-of-flight picking and image
%                               reconstruction.
%
%
%
% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 10.08.2022
%
% This script is part of the r-Wave Toolbox (http://www.r-Wave.org).
% Copyright (c) 2020 Ashkan Javaherian and Ben Cox

para.num_worker_pool = 16;
para.do_calculate_tofs = true; 
para.matrix_construction_method = 'bent-ray';
para.linearisation_approach = 'absolute';
para.choose_last_us_position = 'all';
para.do_remove_bottom_emitters = true;
para.transducer_height_position = -0.10;
para.z_pos_height = 1e-2; 
para.gridtoray_interp = 'Bilinear';
para.reconstruct_image = true;
para.save_results = false;
para.display_transducers = false;
para.save_tofs = true;




% get the taken approach for solving the arising linearised subproblems
switch para.linearisation_approach
    case 'difference'

error('The difference approach for linearisation has been deprectaed.')

% get the approach for solving the linearised subproblems        
linear_subproblem_method = 'steepest_descent'; 
% get the smoothing window size
para.smoothing_window_size = 13;
% get the normalised radius of the mask for image reconstruction
para.mask_coeff = 0.87;

    case 'absolute'
        
% get the approach for solving the linearised subproblems            
para.linear_subproblem_method = 'sart';

% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% get the smoothing window size
para.smoothing_window_size = 5;
% get the normalised radius of the mask for image reconstruction
para.mask_coeff = 1.06;

end


% get the number of iterations
if strcmp(para.linear_subproblem_method, 'sart')
    para.num_iterout = 8;
else
    para.num_iterout = 6;
end


% display the approach for solving the linearised subproblems
disp(['The chosen approach for solving the linear subproblems is:'...
    para.linear_subproblem_method])

switch para.linear_subproblem_method
    case 'steepest_descent'
        
        % get the number of iterations for solving the linearised
        % subproblems
        para.num_iterin = 10;
        
        % get the step length
        para.step_length = 1;
        
    case 'conjugate_gradient'
        
        para.num_iterin = 5;
        
        % choose a step length, a factor that is multiplied
        % by the update direction solved from each linearised
        % subproblem in order to ensure the convergence
        para.step_length = 0.2;
        
    case 'sart'
        
        para.num_iterin = 10;
        
        % choose a step length, a factor that is multplied
        % by the update direction solved from each linerised
        % subproblem in order to ensure the convergence
        para.step_length = 0.2; 
end


%% ========================================================================
% CHOOSE THE LAST US POSITION
%==========================================================================
switch para.choose_last_us_position
    case 'early'
        
        % If set 'early', the index of the last angular position
        % included in the image reconstruction is set manually few number
        % of angular positions. If a value larger than all the US position
        % is chosen, the last us position is chosen 3 by default. Three
        % angular positions is equivalent to an approximate
        % full rotation over all the emitters for a data including 180
        % excitations per position.
        para.last_angular_position = 3;
 
        
    case 'all'
        if strcmp(data_ids.object, 'us_measurements')
            
            % 'us_measurements' is so large, so half of the positions
            % can be excluded by default.
            para.last_angular_position = 31;
        else
            
            % If set 'all', only the last position is removed.
            para.last_angular_position = 100;
        end
end

% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

if strcmp(para.choose_last_us_position, 'early')
    
    % The results are not saved, if the last us position is chosen
    % 'early'
    para.save_results = false;
end
    
%if ~any(strcmp(machine_name(1:end-1), {'kinsler','nyborg'}))
    
    % if the codes are not run on the UCL's servers, 1) the TOFs must be
    % available and loaded. 2) the last US position must be st all.
    % para.do_calculate_tofs = false;
   %  para.choose_last_us_position = 'all';

% end



if para.do_calculate_tofs
    
    % if do_calculate_tofs is set true, choose the data set for only-water data
    only_water_measurement_date = '05-2021';
else
    
    % if do_calculate_tofs is set false, set a nan value
    only_water_measurement_date = nan;
end

%% ========================================================================
% DEFINE THE BOOLEANs FOR REQUESTING A TASK FROM THE SCRIPT
%==========================================================================
% Boolean controlling whether the images are saved on a path for the images or not
% save_images = false;


% Boolean controlling whether the computed TOFs are postprocessed or not.
% If true, the bad tofs are determined, and replaced by the median values
% from the neighboring transducers. Because the distance between the
% transducers are large for the Pammoth system, this is set false.
% (Default = true)
para.postprocess_tofs = false;

if para.postprocess_tofs
    
% Choose a factor between zero and one for the min/max tolerance for the outliers
% in the tof picking. For bad data, for example because of movements, smaller 
% tolerance (smaller value) could be used. (Default = 3.00)
para.tofs_fac = 5;

if para.tofs_fac < 0.5
    
    % Give a warning, if the tofs_fac is chosen smaller than 0.5.
    warning(['A small tolerance for accepting the TOFs may affect the true TOFs.']);
end

else
    
    % this factor is not used when postprocess_tofs is set false.
    para.tofs_fac = nan;

end

% Boolean controlling whether the TOFs from bad transducers are removed or
% not. If this is set true, the bad transducers based are removed from 
% the difference TOF map, and are set zero. The bad transducers found by
% eye inspection from the difference TOF map are consistent with those in 
% the pdf file provided by the PA imaging company. (default = true)
para.manual_remove_bad_transducers = true;

% Boolean controlling whether the boundary of the cone for accepting the
% emitter-receiver pairs are smoothed or not. Setting this parameter true
% will considerably reduce the artefact in water, and reduces erroneous
% refractions, but it may affect the image inside the cup. 
para.smooth_cone_boundary = false;

% A parameter determining how the included  emitter-receiver pairs are chosen,
% based on 'distances'or 'open_angle'. (default: 'open_angle')
para.choose_emitter_receiver = 'open_angle';

% Boolean controlling whether the only-water data is used as a reference,
% or the image reconstruction is done only using the object-in-water data.
% (Default : true)
para.include_only_water = true;

% Boolean controlling whether the emitters at the bottom of the cup are excluded from
% image reconstruction or not. All the rays which are intialised from those
% emitters and intercepted by all the receivers travel only across the
% water.
% do_remove_bottom_emitters = true;

% Boolean controlling whether a Binary mask is used for image reconstruction or not
% apply_mask_cup = false;


% The transducers below the chosen height value [m] are ignored in time-of-flight
% picking, and image reconstruction. It is recommended that the user determines 
% the best height value for each cup size, and creates a look-up table.
%if do_remove_bottom_emitters
%    transducer_height_position = - 0.10;
%end

% a factor for the extension of the grid beyond the maxium position of the
% transducers along each Cartesian coordinate
%switch para.gridtoray_interp
%    case 'Bilinear'
        
    %    grid_expansion_coeff = 1.05;
        
%    case 'Bspline-Denis'
        
        % B-spline uses four adjacent grid points for interpolation, so the
        % grid is enlarged such that the adjacent grid points about the
        % grid points which are close to the edge of grid do not exceed
        % the edge of the computational grid
        grid_expansion_coeff = 1.10;
        
%end

% A minimum distance threshold for choosing emitter-receiver pairs included
% in the image reconstruction. Only the rays linking emitter-receiver
% pairs that hit the breast cup contain information for Radon-like transmission
% image reconstruction.
minimum_distance = 0.18;

if minimum_distance > 0.20
    error(['The minimum distance must not be chosen larger than' num2str(20) 'cm']);
end

% A range of sound speed for choosing a time window for maximum possible
% varitation of first-arrival times for the object-in-water data and only-water data.
% This time window is chosen based on assuming a minimum-maximum
% homogeneous sound speed inside the bowl. (Default = [-50, 50] for the
% object-in-water data, and [-20, 20] for the only water data.) 
% In general, the better the information about the shot time of the
% excitation pulse, namely US shot time, be available, the narrower the ranges 
% for choosing a time window for the first-arrival of the signals could be. 
para.sound_speed_ranges_diff = 50 * [-1, +1];    %  [m/s]

% The angle of view for displaying the transducers during the rotation
% This can be set 'horizontal' or 'oblique', here is set always 'oblique'
% for better visualisation.
% view_angle = 'oblique';


%% ========================================================================
% GET THE ALWAYS FIXED PARAMETERS
%==========================================================================
% get the sampling rate [Hz] for the excitation pulse
time_info.emit_frequency = 5.0e7;

% get the sampling rate [Hz] for the measured ultrasound signals 
time_info.us_frequency = 2.5e7;

% get the time index [a.u.] corresponding to the shot time of the signal,
% according to the pdf file provided by the PA imaging company
time_info.shot_time_index = 637;

% get the radius [m] of the bowl, the hemi-spherical surface on which the 
% transducers are positioned.
radius_bowl = 13e-2;

% get the grid spacing [m]
grid_spacing = 2e-3;

% The number of data segments for loading data and time-of-flight picking.
% This will be set the number of the US angular positions, if it is chosen larger than
% the number of US positions (default : 100). The data are loaded
% part-by-part in order to reduce the memory required for the computations.
num_us_segment = 100;

%==========================================================================
% CHOOSE THE EMITTER-RECEIVER PAIRS
%%=========================================================================
% The 0/1 open angle of a binary cone for each emitter. The axis of the cone
% connects the emitter (on the vetex of the cone) to the center of the bowl.
% For each emitter, only the receivers inside the cone are included. This cone
% is chosen based on the fact that the normal vector to all the transducers(emitters)
% face to the centre of the bowl. Alternatively, the distances can be used.
% Note that both approches, 'distances' or 'open_angle' are the same by
% choosing equivalent prameters. The open angle is chosen based on the cup size.
para.open_angle = max(pi/6, pi/4 - 3 * (8 - data_ids.cup_size) * pi/180);  % pi/4;  %

% Display the cup size
disp(['The cup size is:' num2str(data_ids.cup_size)])

% Display the chosen open angle
disp(['The chosen open angle is:' num2str(para.open_angle * 180/pi) ' degrees'])


%% ========================================================================
% GET THE DATA PATH AND THE TABLE
%==========================================================================
% get the path associated with the object-in-water data, and the table for
% the corresponding angular positions, the temperatures for the water 
% encompassing the object during the rotations, and the cup size.
[data_paths_object, table_array_position, water_object_temperatures,...
    ~] = setPathGetInfos(data_ids, para.last_angular_position);

%%=========================================================================
%  DEFINE THE PATH FOR STORING AND LOADING THE DATA FOR THE OBJECT ON THE LOCAL
%  MACHINE
%==========================================================================
% define a path for the user on the local machine
user_data_paths = [];
user_data_paths.main_directory = data_ids.main_path_user_data;
user_data_paths.directory = data_paths_object.directory;
user_data_paths.data_name = data_paths_object.data_name;

% convert the user data path from struct to a string
user_data_path = [user_data_paths.main_directory, user_data_paths.directory,...
    user_data_paths.data_name '/'];


% get the default directory for saving the results
if isempty(results_directory)
    
    % get the directory for the results
    results_directory = [data_ids.main_path_user_data, 'SoSreconstructions/', data_ids.id];
    
    % make the directory, if it does not exist.
    makeDirectory(results_directory)
    
end

% get the initial characters of the file name for saving the results
switch para.matrix_construction_method
    case 'straight-ray'
        initial_file_name = 'st_';
    case  'bent-ray'
        initial_file_name = 'ss_';
end

%% ========================================================================
% GET POSITION OF EMITTERS/RECEIVERS
%==========================================================================
% get the position of emitters and receivers
position_args = {'do_calculate_tofs', para.do_calculate_tofs};
[emitter, receiver, emit_shot_index, us_shot_index, emit_firing_transducer_index,...
    rotation_indices] =...
     getPositions(data_paths_object, table_array_position, radius_bowl,...
     user_data_paths, position_args{:});

% get the number of dimensions
dim = size(emitter.positions, 1);
 
% get the number of emitters
num_emitter = size(emitter.positions, 2);

% get the number of receivers
num_receiver = size(receiver.positions{1}, 2);

% get the number of excitations
% num_excitation = length(emit_shot_index);

%% ========================================================================
% GET THE ONLY-WATER DATA
%==========================================================================
% laod the only-water data, and the sound speed for only-water data.
% This data is availabe for only one full rotation over all emitters, because
% the tofs in only water are independent of the angular positions, and are only
% affeced by the distance of emitter-receiver pairs, so the data for one full
% rotation over all emitters can be extended to all excitations.
[data_water, ~, ~, sound_speed_only_water] = loadDataOnlyWater(data_paths_object,...
    only_water_measurement_date, num_receiver);


%% ========================================================================
% GET THE INFORMATION ABOUT THE WATER FOR OBJECT-IN-WATER DATA
%==========================================================================
% get the information about the water encompassing the object 
[sound_speed_object_correction_coeff, refractive_background_object_alldata,...
    sound_speed_object_water] = getWaterObject(water_object_temperatures, rotation_indices,...
    sound_speed_only_water, true, num_emitter, num_receiver);

%% ========================================================================
% COMPUTE (OR LOAD) THE TIME-OF-FLIGHTS
%==========================================================================

% the optional inputs for computing the Tofs
tof_args = {'num_worker_pool', para.num_worker_pool,...
'do_calculate_tofs', para.do_calculate_tofs,...
'include_only_water', para.include_only_water,...
'choose_last_us_position', para.choose_last_us_position,...
'choose_emitter_receiver', para.choose_emitter_receiver,...
'open_angle', para.open_angle,...
'num_us_segment', num_us_segment,...
'postprocess_tofs', para.postprocess_tofs,...
'tofs_fac', para.tofs_fac,...
'smooth_cone_boundary', para.smooth_cone_boundary,...
'manual_remove_bad_transducers', para.manual_remove_bad_transducers,...
'sound_speed_ranges_diff', para.sound_speed_ranges_diff,...
'minimum_distance', minimum_distance,...
'save_results', para.save_tofs,...
'absorption_frequencies', nan};

% compute the time-of-flights
[tof_data, weight_ray] = calcLoadTofs(data_paths_object, data_water, emitter,...
    receiver, time_info, emit_shot_index, us_shot_index, emit_firing_transducer_index,...
    rotation_indices, sound_speed_object_correction_coeff, refractive_background_object_alldata,...
    sound_speed_object_water, sound_speed_only_water, user_data_path, tof_args{:});

 if para.reconstruct_image
%% ========================================================================
% REMOVE THE BOTTOM EMITTERS
% =========================================================================

% remove the bottom emitters
if para.do_remove_bottom_emitters
    
    % get the 1 x num_emitter binaries for exciations with emitters with z 
    % Cartesian coordinate above the transducer_height_position
    top_emitter_binaries  = emitter.positions(3, :) > para.transducer_height_position;
    
    % update the emitter positions 
    emitter.positions = emitter.positions(:, top_emitter_binaries);
    
    % update the rotataion indices
    rotation_indices = rotation_indices(top_emitter_binaries);
    
    % get the num_receiver x num_emitter linear indices for the binaries
    % for exciations with emitters with z 
    % Cartesian coordinate above the transducer_height_position
    top_emitter_binaries_linear = repmat(top_emitter_binaries, [num_receiver, 1]);
    
    % update the vector of time-of-flight differences
    tof_data = tof_data(top_emitter_binaries_linear);
    
    % update the refractive index of the water encompassing the object
    % during the roatation of the bowl for all emitter-receiver pairs
    refractive_background_object_alldata = ...
        refractive_background_object_alldata(top_emitter_binaries_linear);
   
    % update the computed weights for the rays. The weights are set 1 for
    % all the rays.
    weight_ray = weight_ray(top_emitter_binaries_linear);
    
end

% Display the position of the top emitter
disp(['The position of the top transducer is:'...
    num2str(1000 * max(emitter.positions(3,:))) 'mm'])

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
    recon_grid = makeReconstructionGrid(grid_spacing * ones(1, dim),...
        radius_bowl, grid_args{:});
    
    % define an ROI (binary mask) for the first iteration
    mask_raytracing  = recon_grid.x.^2 + recon_grid.y.^2 + recon_grid.z.^2 <= ...
        (para.mask_coeff * radius_bowl)^2  & recon_grid.z < para.z_pos_height- 2 * recon_grid.dx;  
          
   % get the binary mask for the image reconstruction. 
   % switch para.linearisation_approach
   %     case 'absolute'
   %         mask_reconst  = recon_grid.x.^2 + recon_grid.y.^2 + recon_grid.z.^2 <= ...
   %             ((para.mask_coeff-0.15) * radius_bowl)^2  &   recon_grid.z <  1/4 * para.z_pos_height;  
   %     case 'difference'
   %         mask_reconst  = recon_grid.x.^2 + recon_grid.y.^2 + recon_grid.z.^2 <= ...
   %             ((para.mask_coeff-0.05) * radius_bowl)^2  &   recon_grid.z <   1/4 * para.z_pos_height;  
   % end
       
    % The binary mask is chosen a digital cup larger than the physical cup
    % such that all the possible movements of the physical cup is included
    % in the binary mask.
    % get a binary mask for gradient of the refractive index 
      mask_reconst = getBreastCupMaskModified(recon_grid, data_ids.cup_size, -4e-3);
      mask_reconst(recon_grid.z > 1/4 * para.z_pos_height) = false;
    % mask_reconst = true(recon_grid.size);
   
    %% ====================================================================
    % RECONSTRUCT THE IMAGE USING THE ITERATIVE ALGORITHM
    %======================================================================
    % get the sound speed in only water
    water = [];
    water.refractive_object_water = refractive_background_object_alldata;
    water.sound_speed_only_water = sound_speed_only_water;
    water.sound_speed_object_water = sound_speed_object_water;
    
    % get the optional parameters
    reconst_args = {'matrix_construction_method', para.matrix_construction_method,...
        'linearisation_approach', para.linearisation_approach,...
        'linear_subproblem_method', para.linear_subproblem_method,...
        'gridtoray_interp', para.gridtoray_interp,...
        'num_iterout', para.num_iterout, 'num_iterin', para.num_iterin,...
        'step_length', para.step_length, 'binaries_emitter_receiver',...
        para.choose_emitter_receiver, 'open_angle', para.open_angle,...
        'smoothing_window_size', para.smoothing_window_size};
   
    % reconstruct the image iteratively
    [img, out, ~] = reconstructIterativeSoundspeedPammoth(weight_ray .* tof_data, recon_grid,...
       emitter.positions, receiver.positions, water, mask_raytracing, mask_reconst,...
       rotation_indices, [], reconst_args{:});

    % include the time and date of the image reconstruction in the output
    % struct
    out.time_stamp = datetime;
    
    % include the only-water sound speed in the output struct
    out.sound_speed_only_water = sound_speed_only_water;
    
    % include the mean sound speed of the water encompassing the object in the output struct
    out.sound_speed_object_water = sound_speed_object_water;
    
    % include the date on which the only water data are measured in the data_ids
    data_ids.only_water_measurement_date = only_water_measurement_date;
    
    %% ====================================================================
    % SAVE THE RESULTS
    %======================================================================
    if para.save_results
        
        % save the results in the specified directory
        save([results_directory, initial_file_name, num2str(data_ids.measurement_number)...
            '_' '1_' para.linear_subproblem_method  '.mat'], 'img', 'recon_grid', 'data_ids',...
            'out', 'para', '-v7.3');
        
    end
    
 else
     
     img = [];
     recon_grid = [];
     out = [];
     
 end

    
end