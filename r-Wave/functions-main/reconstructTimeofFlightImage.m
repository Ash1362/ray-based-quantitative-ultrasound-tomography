function [img, recon_grid, tof_discrepancy, ray_initial_angles, out, para] =...
    reconstructTimeofFlightImage(data_object, data_water, time_array,...
    emitter, receiver, sound_speed_water, simulation_prop, tof_path, varargin)
%RECONSTRUCTTIMEOFFLIGHTIMAGE computes a sound speed image using time-of-flights of
%the ultarsound time series.
%
%
% DESCRIPTION:
%             reconstructTimeofFlightImage computes the image of speed of
%             sound from ultrasound time series
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
%                           and num_receiver the number of receivers
%       data_water        - the ultrasound time series measured from water.
%                           This is a matrix of size num_receiver x num_time
%                           x num_emitter with num_emitter the number of emitters,
%                           num_time the number of time instants for measurement
%                           and num_receiver the number of receivers
%       time_array        - a time array [s] of size 1 x num_time on which the 
%                           the UST time series are measured. 
%       emitter           - a struct which defines the properties of the
%                           excitation as follows: This includes the
%                           fields 'positions', 'pulse', 'pulse_duration',
%                           and 'shot_time'
%       emitter.positions - 2/3 x num_emitter array of Cartesian position
%                           of the center of the emitter objects
%       emitter.rotation_indices  - a vector of size num_emitter x 1 and
%                           includes the rotation indices of the excitations.
%                           If this field is not given, the transducers are
%                           assumed fixed (not rotated)
%       pulse             - a vector of size 1 x num_te with num_te <= num_time.
%       shot_time         - the first arrival of the excitation pulse, which
%                           can be negative [sec]
%       receiver          - a struct which includes the fields
%       positions         - 2/3 x N array of Cartesian position of the
%                           centre of the receiver objects. For rotating
%                           ultrasound systems, the position of the receivers
%                           is changed for different excitations. For the rotating case,
%                           this is a cell array with length num_angle (the number of rotation
%                           angles) with each cell a 2/3 x num_transducer array of
%                           Cartesian position of transducers for each
%                           angle. (num_receiver = num_angle x num_transducer)
%      sound_speed_water   - the sound speed [m/s] in only water
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
%      'z_offset'           - the z offset [m] between the grid for k-Wave simulation
%                             and the grid for image reconstruction. The image 
%                             reconstruction is done on a detection surface
%                             (ring) with centre at the origin.
%      tof_path             - the data path for storing the TOFs
%
%
% OPTIONAL INPUTS:
%      'num_worker_pool'   - the number of workers for parallel programming
%      'reconstruct_image' - Boolean controlling whether the image is
%                            reconstructed, or only time-of-flights are
%                            computed and are given as the input of the function.
%      'emittersegments_to_nworkers' - The number of segments for dividing emitters to the number of
%                            cpu workers.  For avoiding memory leakage during the
%                            parallel programming, the emitters are first
%                            divided to a number of segments, and then
%                            the segment are consecutively run over the threads.
%      'matrix_construction_method' - the method for construction of the
%                           system matrix and reconstructing the
%                           time-of-flight image. (default: 'bent-ray')
%      'linearisation_approach' - the approach for linearisation of the
%                           misfit function. This can be either 'difference',
%                           or 'absolute'. (default: 'difference')
%      'linear_subproblem_method' - the method for solving the linearised
%                           subproblems, which can be either 'conjugate_gradient',
%                           'sart', or 'steepest_descent'.
%                           (Default : 'sart')
%      'grid_spacing'      - the grid spacing for image reconstruction
%                            (default: 1e-3 [m])
%      'z_pos_height'      - the maximum z position [m] of the grid points, if
%                            the reconstruction geometry is set 'real'
%                            (Default: 1e-2)
%      'mask_coeff'        - a coefficient which defines the normalised radius of
%                            the circular binary mask with radius = mask_coef * detec_radius
%                            for image reconstruction. Here, detec_radius
%                            is the radius of the detection surface, which
%                            must be chosen greater 1. (Default : 1.03)
%      'tof_frac_peak_threshold' - the fraction of the peak amplitude for
%                            chosing the right edge of the time window
%                            for time-of-flight picking. The left edge of the 
%                            the time window is chosen based on an assumption
%                            of a maximum sound speed for the entire medium.
%                            If it is set zero, the right edge of the time 
%                            window is chosen based on an assumption of a
%                            homogeneous minimum sound speed. (Default: 0.5)
%      'tof_outliers'       - the Boolean controlling whether the outliers in the
%                            difference time-of-flight matrix are replaced
%                            using an interpolation approach or not.                         
%      'binaries_emitter_receiver' - the method for choosing the
%                               emitter-receiver pairs. This can be
%                               'open_angle', i.e., the angle between
%                               a geometrical vactor connecting the
%                               emitter to the centre of the detection suface (ring),
%                               and another vector connecting the emitter to
%                               receivers. The emitter-receiver pairs with open angles
%                               larger than a threshold are excluded. This can also be set
%                               'distances', i.e., emitters and receivers with
%                               distances smaller than a specific threshold
%                               are excluded. Bothe approaches are
%                               equivalent.
%      'minimum_distance_coeff' - a scalar coefficient indicating the minimum distance
%                            of the transducers for including in the image
%                            reconstruction. The emitter-receiver pairs with distances
%                            less than minimum_distance = minimum_distance_coeff * detec_radius
%                            The measured signals associated with close emitter-receiver
%                            pairs may be deteriorated by the effects of directivity
%                            of the transducers. In addition, rays linking the
%                            close emitter-receiver pairs only travel inside water.
%        'open_angle'          - the maximum angle between two vecors, one
%                                connecting the emitter to the centre of the
%                                detection surface (ring), and another connecting
%                                the emitter to the receivers. For each emitter,
%                                the included receivers will be inside a cone
%                                with axis the vector connecting the emitter to
%                                the centre of the detection surface.
%      'smoothing_window_size' - the size of the smoothing window applied
%                              on the updates of the refractive index.
%      'raylinking_method'     - the method for ray linking. This can be either
%                              'Regula-Falsi' or 'Secant' for 2D case, or
%                              'Quasi-Newton' for 3D case.
%      'raytracing_method'  -  the method for ray tracing, which can be
%                              'Mixed-step', 'Dual-update', 'Characteristics',
%                              or 'Runge-kutta-2nd'. (Default : 'Mixed-step')
%      'gridtoray_interp'      - the method for interplating from grid points to
%                                rays and vice versa. This can be 'Bilinear', or
%                               'Bspline'.
%      'max_raylinking_iter'   - the maximum number of iterations for ray
%                                linking
%      'raylinking_threshold'  - the threshold distance [m] for stopping the
%                               ray linking iterations. If the distance
%                               between the interception point of the ray
%                               by the detection surface (ring) and the centre
%                               of the receiver becomes less than
%                               raylinking_threshod, the ray linking
%                               algorithm will be stopped.
%      'raytogrid_spacing'    - the spacing of the sampled points along the ray
%                               to the grid spacing (Default: 1/2 (2D case),
%                               1 (3D case))
%      'sound_speed_time_window' - an assumption of homogeneous min/max for
%                               the sound speed. The distance between each
%                               emitter-receiver pair over these min/max
%                               values is used for choosing a time window
%                               for AIC first-arrival picking for the
%                               corresponding ultrasound signal. (Defult:
%                               50)
%     'emitter_downsampling_rate' - the downsampling rate applied to the third
%                            (emitter) dimension of the data simulated by k-Wave (Default: 1)
%                             (This is only used if data is empty, and therefore, the time-of-flights are loaded.
%                              For simulation studies, it must be set the
%                              same as in 'simulateSettingData.m')
%     'receiver_downsampling_rate' - the downsampling rate applied to the first
%                            (receiver) dimension of the data simulated by k-Wave (Default: 1)
%                             (This is only used if data is empty, and therefore, the time-of-flights are loaded.
%                              For simulation studies, it must be set the
%                              same as in 'simulateSettingData.m')
%
% OUTPUTS:
%      img                    - the reconstructed image of the sound speed
%      recon_grid             - the computational grid for image reconstruction
%      tof_discrepancy        - the num_reciver x num_emitter matrix of the
%                               the discrepancy of the time-of-flights
%                               for the object-in-water data and only-water
%                               data
%      ray_initial_angles     - the initial angle of the linked (optimal) rays for all
%                               emitter-receiver pairs after ray linking applied to the
%                               last update of the sound speed
%      out                    - a struct containing the stored data of
%                               the iterative reconstruction. This
%                               struct includes the fields:
%      'norm_res'             - the L2 norm of res updates for the cg algorithm
%                               and the L2 norm of update directions of the sound speed
%                               for the steepest descent algorithm at each inner iteration
%      'residual_norm'        - the first component of norm_res
%      'num_rays'             - the number of rays for ray linking between
%                               the emitter-receiver pairs
% cartesian_position_end_point' - the cartesian position of the end point of
%                               the linked rays for all iterations
%      'system_matrix_time'   - the cpu time for construction of the system
%                               matrix
%      'update_time'          - the cpu time for solving each linear subproblem
%       'imgs'                - all the reconstructed images of the sound speed
%                               after solving the linearised subproblems
%                               If the sound_speed_ground_truth is given,
%                               this struct can also includes the fields:
%       'relative_error'      - the L2 norm of the discrepancy of the
%                               updates of the sound speed and the ground
%                               truth over the L2 norm of the disrepancy of
%                               water and the ground truth time 100.
%        'rmse'               - the rmse of the reconstructed images with
%                               respect to the ground truth.
%        para                 - the ptional inputs used for the image
%                               reconstruction

%
%
% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 05.10.2022
%
% This script is part of the r-Wave Toolbox
% Copyright (c) 2022 Ashkan Javaherian


para.num_worker_pool = 16;
para.reconstruct_image = true;
para.matrix_construction_method  = 'bent-ray';
para.raytracing_method = 'Mixed-step';
para.linearisation_approach = 'absolute';
para.linear_subproblem_method = 'sart';
para.grid_spacing = 2e-3;
para.binaries_emitter_receiver = 'distances';
para.minimum_distance_coeff = 0.7368;
para.raylinking_threshold = 1e-6;
para.sound_speed_time_window = 50;
para.emitter_downsampling_rate = 1;
para.emitter_downsampling_rate = 1;
para.tof_frac_peak_threshold = 0.5;
para.tof_outliers = false;

% get the number of the dimensions
dim = size(emitter.positions, 1);

switch dim
    case 2
        
        para.raylinking_method ='Secant';
        para.gridtoray_interp = 'Bspline';
        para.max_raylinking_iter = 1000;
        para.raytogrid_spacing = 1/2;
        para.emittersegments_to_nworkers = 4;
        para.open_angle = pi/3;
        para.z_pos_height = 0;
        
    case 3
        
        para.raylinking_method ='Quasi-Newton';
        para.gridtoray_interp = 'Bilinear';
        para.max_raylinking_iter = 500;
        para.raytogrid_spacing = 1;
        para.emittersegments_to_nworkers = 32;
        para.open_angle = pi/4;
        para.z_pos_height = 1.5e-2;
        
end

% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


% give error for inconsistent optional inputs
if strcmp(para.linearisation_approach, 'absolute') && ...
        strcmp(para.linear_subproblem_method, 'steepest_descent')
    error(['Using the absolute approach for linearisation, the steepest descent method'...
        'converges very slowly. The user must use ''conjuage_gradient'' or ''sart''.'])
end

%if strcmp(para.linearisation_approach, 'difference')
%    error(['This linearisation approach has been deprecated, and will be omitted.'...
%       'Please use only the absolute approach for linearisation!'])
%end
%    ~strcmp(para.linear_subproblem_method, 'steepest_descent')
%    error(['If the approach taken for linearisation is ''difference'', the approach for'...
%    'solving the linearised subproblem must be set ''steepest_descent''.'])
%end


% choose the smoothing window size for the supported  grid spacing [m]
% Give an error if the grid spacing is not set 2mm, 1mm or 0.5 mm
if  para.grid_spacing < (5e-4)-(1e-10)
     error('The grid spacing smaller than 0.5mm is not necessary for a ray-based method.')
elseif  abs(para.grid_spacing - (5e-4))< 1e-10
    para.smoothing_window_size = 13;
elseif abs(para.grid_spacing - (7.5e-4))< 1e-10
    para.smoothing_window_size = 9;  
elseif abs(para.grid_spacing - (1e-3))< 1e-10
    para.smoothing_window_size = 7;
elseif abs(para.grid_spacing - (2e-3))< 1e-10
    para.smoothing_window_size = 5;
elseif para.grid_spacing > (2e-3)+(1e-10)
    error('The grid spacing must not be larger than 2mm.')
else
    error('The supported grid spacings are 2mm, 1mm, 0.75mm or 0.5mm.')
end


if rem(para.smoothing_window_size, 2) == 0 && para.smoothing_window_size > 3
    error('The smoothing window size must not be an even value, if it is larger than 3.')
end


% get the parameters for image reconstruction
switch dim
    case 2
        switch para.raylinking_method
            case 'Regula-Falsi'
                para.num_iterout = 5;
            case 'Secant'
                
                switch para.matrix_construction_method
                    case 'straight-ray'
                        switch para.linear_subproblem_method
                            case 'conjugate_gradient'
                                para.num_iterout = 6;
                            case 'steepest_descent'
                                para.num_iterout = 6;
                            case 'sart'
                                para.num_iterout = 6;
                        end
                    case 'bent-ray'
                        switch para.linear_subproblem_method
                            case 'conjugate_gradient'
                                para.num_iterout = 6;  %(Default : 6)
                            case 'steepest_descent'
                                para.num_iterout = 3;
                            case 'sart'
                                para.num_iterout = 5;  %(default : 5)
                        end
                end
                
        end
    case 3
        switch para.linear_subproblem_method
            case 'conjugate_gradient'
                para.num_iterout = 12;
            case 'steepest_descent'
                para.num_iterout = 4;
            case 'sart'
                para.num_iterout = 12;
        end
end

% the number of inner iterations
switch dim
    case 2
        switch   para.linear_subproblem_method
            case 'steepest_descent'
                
                % get the number of iterations for solving each
                % linearised subproblem
                para.num_iterin = 300;
                
                % get the step length
                para.step_length = 1;
                
            case 'conjugate_gradient'
                
                % get the number of iterations for solving each
                % linearised subproblem
                para.num_iterin = 5;
                
                % choose a step length, a factor that is multiplied
                % by the update direction solved from each linerised
                % subproblem such that the convergence is ensured
                para.step_length = 0.2;
                
            case 'sart'
                
                % get the number of iterations for solving each
                % linearised subproblem
                para.num_iterin = 5;
                
                % choose a step length, a factor that is multiplied
                % by the update direction solved from each linerised
                % subproblem such that the convergence is ensured
                para.step_length = 1;
                
        end
    case 3
        switch   para.linear_subproblem_method
            
            case 'steepest_descent'
                
                % get the number of iterations for solving each
                % linearised subproblem
                para.num_iterin = 300;   % (default:400)
                
                % get the step length
                para.step_length = 1;
                
            case 'conjugate_gradient'
                
                % get the number of iterations for solving each
                % linearised subproblem
                para.num_iterin = 5;
                
                % choose a step length, a factor that is multiplied
                % by the update direction solved from each linerised
                % subproblem such that the convergence is ensured
                para.step_length = 0.2;
                
            case 'sart'
                
                % get the number of iterations for solving each
                % linearised subproblem
                para.num_iterin = 5;
                
                % choose a step length, a factor that is multiplied
                % by the update direction solved from each linerised
                % subproblem such that the convergence is ensured
                para.step_length = 1;
        end
end

% get the normalised radius of the mask for image reconstruction
switch para.linearisation_approach
    case 'absolute'
        para.mask_coeff = 1.03;
    case 'difference'
        error('The difference approach for image reconstruction is deprecated.')
        para.mask_coeff = 0.87;
end

% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% error messages
if para.mask_coeff < 1  &&  strcmp(para.linearisation_approach, 'absolute')
    error(['Using the absolute approach, the mask_coeff, which is the normalised radius'...
        'of the binary mask for image reconstruction, must be larger than the radius of'...
        'the detection surface (ring).'])
end


if nnz([isempty(data_object), isempty(data_water), isempty(time_array)]) > 0
    
    if nnz([isempty(data_object), isempty(data_water), isempty(time_array)]) < 3
        error(['The object-in-water and only-water data and the time array must all be empty, or'...
            ' all nonempty.'])
    end
    
else
    
    if any([size(data_object,1)-size(data_water,1);...
            size(data_object, ndims(data_object))-...
            size(data_water, ndims(data_water))])
        error('The sizes of object-in-water and only-water data sets must be consistent.')
    end
    
end



% give an empty variable to the field rotation_indices, if it does not
% exist, i.e., the position of receivers are fixed with changes in
% excitations.
if ~isfield(emitter, 'rotation_indices')
    
    % add an empty field 'rotation_indices' to struct emitter.
    emitter.rotation_indices = [];
    
    % change the number of receivers from matrix to cell array
    receiver.positions = mat2cell(receiver.positions,...
        size(receiver.positions, 1),...
        size(receiver.positions, 2));
    
end

if size(receiver.positions{1}, 2) < 3
    error('The number of receivers for each excitation must be at least 3.')
end

% get the radius of the detection surface (if circle or hemi-sphere)
radius1 = norm(receiver.positions{1}(:, 1));
radius2 = norm(receiver.positions{1}(:, 2));
radius3 = norm(receiver.positions{1}(:, 3));


% check if the detection surface is circle or hemi-sphere
if abs(radius1 - radius2)< 1e-10  && abs(radius1 - radius3)< 1e-10
    
    % the detection geometry is circle (2D) or sphere (3D)
    detec_geom = 'sphere';
    
    detec_radius = radius1;
    
    % get the minimum distance [m] for accepting the emitter-receiver pairs
    % included in the image reconstruction
    minimum_distance = para.minimum_distance_coeff * detec_radius;
    
    if dim == 3
        minimum_distance = 1.5 * minimum_distance;
    end
    
    
    
else
    
    if dim == 3
        
        % get the radius of the detection surface in x-y plane
        radius1 = norm(receiver.positions{1}(1:2, 1));
        radius2 = norm(receiver.positions{1}(1:2, 2));
        radius3 = norm(receiver.positions{1}(1:2, 3));
        
    end
    
    
    if dim == 3  && abs(radius1 - radius2)< 1e-10   &&   abs(radius1 - radius3)< 1e-10
        
        
        % the detection geometry is cylinder (only 3D)
        detec_geom = 'cylinder';
        
        error('Not implemented yet!')
        
    else
        
        % get the slope of the detection surface (if linear or planar)
        slope1 = (receiver.positions{1}(2, 2)- receiver.positions{1}(2, 1))./...
            (receiver.positions{1}(1, 2)- receiver.positions{1}(1, 1));
        
        slope2 = (receiver.positions{1}(2, end)-receiver.positions{1}(2, end-1))./...
            (receiver.positions{1}(1, end)- receiver.positions{1}(1, end-1));
        
        if abs(slope1-slope2)< 1e-10   || all(isinf([slope1, slope2]))
            
            % get the number of receivers per excitation
            num_receiver_per_emitter = size(receiver.positions{1}(:, 2));
            
            % the detection geometry is planar (2D or 3D)
            detec_geom = 'plane';
            
            % the assumption is that two planes, one for emitters and
            % another for receivers, face each other with sufficient
            % distance and rotate together around a point(2D) or line (3D),
            % so all emitter-receiver pairs are accepted.
            % However, in practical settings, the user may like to add a
            % maximum permissible distance between the emitter-receiver pairs
            % to ommit the emitter-receiver pairs which are far and do not
            % face each other.
            
            % get the index or indices for computing the radius [m] of rotation
            detec_radius = mean(vecnorm(receiver.positions{1}(:, unique([floor((num_receiver_per_emitter + 1)/2),...
                floor((num_receiver_per_emitter + 1)/2) + rem(num_receiver_per_emitter + 1, 2)])) ) );
            
            % get the minimum distance [m] for accepting the emitter-receiver pairs
            % included in the image reconstruction
            minimum_distance = para.minimum_distance_coeff * detec_radius;
    
        end
        
    end
    
end

if isempty(emitter.rotation_indices)  
    
    % change the number of receivers from cell array to matrix
    receiver.positions = cell2mat(receiver.positions);
end

% get the Boolean controlling whether the time-of-flights are computed or
% loaded.
do_calculate_tofs = ~isempty(time_array);


disp(['The step length is:' num2str(para.step_length)])
disp(['The method for solving the linearised subproblem is:'...
    para.linear_subproblem_method])

% a factor for the extension of the grid beyond the maxium position of the
% transducers along each Cartesian coordinate so that there are enough
% neigboring points for interpolation. the larger interpolation (B-spline)
% was considered.
grid_expansion_coeff = 1 + 0.07 * 1000 * para.grid_spacing;  % 1.07

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
% INTERPOLATE THE MEDIUM'S PROPERTIES TO GRID FOR IMAGE RECONSTRUCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is only used for measuring the error of the iteratively
% reconstructed sound speed images.

% interpolate the sound speed phantom to the grid for image
% reconstruction

if ~isempty(simulation_prop)
switch recon_grid.dim
    case 2
        
        sound_speed_ground_truth = interpn(simulation_prop.x, simulation_prop.y,...
            simulation_prop.sound_speed, recon_grid.x , recon_grid.y);
        sound_speed_ground_truth(~isfinite(sound_speed_ground_truth)) = sound_speed_water;
        
        % display the sound speed phantom (ground truth)
        figure;imagesc(sound_speed_ground_truth); axis image; colorbar;
        
    case 3
        
        sound_speed_ground_truth = interpn(simulation_prop.x, simulation_prop.y,...
            simulation_prop.z - simulation_prop.z_offset, simulation_prop.sound_speed, recon_grid.x, recon_grid.y,...
            recon_grid.z, 'spline');
        
        sound_speed_ground_truth(~isfinite(sound_speed_ground_truth)) = sound_speed_water;
        
        % display the sound speed phantom (ground truth)
        scrollView(sound_speed_ground_truth, 1);
end

else
    
    % allocate an empty varibale, if the information about simulation is
    % not given, for example when experimental data is used.
    sound_speed_ground_truth = [];
end

%% ========================================================================
% DEFINE A REGION-OF-INTEREST FOR IMAGE RECONSTRUCTION
% =========================================================================
switch dim
    
    case 2
        
        switch detec_geom
            
            case 'sphere'
                
                
                % get a binary mask for ray tracing
                mask_raytracing = recon_grid.x.^2 + recon_grid.y.^2 < (para.mask_coeff * detec_radius)^2;
               
                % get a binary mask for image reconstruction
                mask_reconst = recon_grid.x.^2 + recon_grid.y.^2 < ((para.mask_coeff-0.05) * detec_radius)^2;
                
                
            case 'plane'
                
                % get a binary mask for ray tracing
                mask_raytracing = true(recon_grid.size);
                
                % get a binary mask for image reconstruction
                mask_reconst = recon_grid.x.^2 + recon_grid.y.^2 < ((para.mask_coeff-0.05) * detec_radius)^2;
                
            case 'cylinder'
                
                % give an error
                error('The 2D case cannot have cylindrical geomtery for position of receivers.')
                
        end
        
    case 3
        
        switch detec_geom
            
            case 'sphere'
                
                % get the maximum z coordinate for the mask for ray tracing and
                % image reconstruction
                z_mask_max = 8e-3;
                
                % get a binary mask for ray tracing
                mask_raytracing = recon_grid.x.^2 + recon_grid.y.^2 + recon_grid.z.^2<= ...
                    (para.mask_coeff * detec_radius)^2   &   recon_grid.z < z_mask_max;
                
                switch para.linearisation_approach
                    case 'absolute'
                        
                        % get a binary mask for image reconstruction
                        mask_reconst = recon_grid.x.^2 + recon_grid.y.^2 + recon_grid.z.^2<= ...
                            ((para.mask_coeff-0.15) * detec_radius)^2  &   recon_grid.z < z_mask_max - 6e-3;
                        
                    case 'difference'
                        
                        % get a binary mask for image reconstruction
                        mask_reconst = recon_grid.x.^2 + recon_grid.y.^2 + recon_grid.z.^2<= ...
                            ((para.mask_coeff-0.05) * detec_radius)^2  &   recon_grid.z < z_mask_max - 6e-3;
                end
                
            case 'plane'
                
                error('Not implemented yet!')
                
            case 'cylinder'
                
                error('Not implemented yet!')
                
        end
        
end

if do_calculate_tofs
    
    % the optional inputs for calculating the tim-of-flights from the measured data
    tof_args = {'nWorkerPool', para.num_worker_pool, 'Method', 'Modified_AIC',...
         'SoundSpeedRanges', sound_speed_water + para.sound_speed_time_window * [-1, 1],...
         'binaries_emitter_receiver', para.binaries_emitter_receiver, 'minimum_distance', minimum_distance,...
         'Threshold', para.tof_frac_peak_threshold};
    
    % compute the time-of-flights [s] for the object-in-water data
    [~, tof_het, ~] = timeOfFlightPicking(data_object,...
        [], emitter, receiver, time_array, emitter.rotation_indices, tof_args{:});
    
    % the optional inputs for calculating the time-of-flights from the measured data
    % tof_args = {'nWorkerPool', para.num_worker_pool, 'Method', 'Modified_AIC',...
    %   'SoundSpeedRanges', sound_speed_water + 2/5 * para.sound_speed_time_window * [-1, 1],...
    %    'binaries_emitter_receiver', 'distances',...
    %    'minimum_distance', minimum_distance, 'Threshold', para.tof_frac_peak_threshold};
    
    % clear the object-in-water data for saving memory
    clear data_object
    

    if ~(length(time_array)-size(data_water, 2))
        
    % compute the time-of-flights [sec] for the only-water data
    [~, ~, tof_hom] = timeOfFlightPicking([],...
        data_water, emitter, receiver, time_array, emitter.rotation_indices, tof_args{:});

    else

    % set tof_hom the TOF for only-water data 
    tof_hom = data_water;

    % set the only-water TOFs for close emitter-receiver pairs zero
    tof_hom(~tof_het)=0;
    
    end

    % clear the only-water data for saving memory
    clear data_water
    
    % compute the discrepancy of time-of-flights [sec]
    tof_discrepancy = tof_het - tof_hom;
    
    % display the difference time-of-flight data
    figure; imagesc(tof_discrepancy); colorbar;
    
    if para.tof_outliers
        
        switch dim 
            case 2
                
        % if requested, remove outliers and replace by a Modified Akima cubic
        % Hermite interpolation
        % see function filloutliers in Matlab
        tof_discrepancy = filloutliers(tof_discrepancy, 'makima', 'mean');
        
        % enforce a minimum bound on the difference TOFs
        tof_discrepancy(tof_discrepancy<-1e-5) = -1e-5;
        
        % enforce a maximum bound on the difference TOFs
        tof_discrepancy(tof_discrepancy>+1e-5) = +1e-5;

        % display the difference time-of-flight data
        figure; imagesc(tof_discrepancy); colorbar;
    
            case 3
                
                error('Not implemented yet.')
        end
        
    end
        

    % get the tof data
    tof_data = vectorise(tof_discrepancy);
    
    if para.emitter_downsampling_rate == 1   &&   para.receiver_downsampling_rate == 1  && ...
           ~isempty(tof_path)
        
        % save the stack vector of time-of-flight discrepancies, if the emitters and receivers
        % are not downsampled.
        save([tof_path, '_tof_sinogram.mat'], 'tof_data');
    end
    
else
    
    % get the number of receivers
    if iscell(receiver.positions)
        num_receiver = size(receiver.positions{1}, 2);
    else
        num_receiver = size(receiver.positions, 2);
    end
    
    % get the number of emitters
    num_emitter = size(emitter.positions, 2);
   
    % load the map of time-of-flight discrepancies
    load([tof_path, '_tof_sinogram.mat'], 'tof_data');
  
    tof_discrepancy = reshape(tof_data,...
        [para.receiver_downsampling_rate * num_receiver,...
        para.emitter_downsampling_rate * num_emitter]);
    
    tof_discrepancy = tof_discrepancy(para.receiver_downsampling_rate:para.receiver_downsampling_rate:end,...
        para.emitter_downsampling_rate:para.emitter_downsampling_rate:end);
end

% display the size of the TOF discrepancy map
disp(['The size of the TOF discrepnacy map is:' num2str(size(tof_discrepancy))])
% display the norm of the TOF difference data
disp(['The norm of the TOF difference data is:' num2str(norm(tof_discrepancy), '%2.5e') 'sec'])

if para.reconstruct_image
    %% ========================================================================
    % RECONSTRUCT THE IMAGE USING THE ITERATIVE ALGORITHM
    %==========================================================================
    % get the sound speed in only water
    water = [];
    water.sound_speed_only_water = sound_speed_water;
    
    % get the optional parameters
    reconst_args = {'matrix_construction_method', para.matrix_construction_method,...
        'linearisation_approach', para.linearisation_approach,...
        'linear_subproblem_method', para.linear_subproblem_method,...
        'raytracing_method', para.raytracing_method,...
        'gridtoray_interp', para.gridtoray_interp,...
        'raytogrid_spacing', para.raytogrid_spacing,...
        'smoothing_window_size', para.smoothing_window_size,...
        'num_iterout', para.num_iterout, 'num_iterin', para.num_iterin,...
        'step_length', para.step_length, 'binaries_emitter_receiver',...
        para.binaries_emitter_receiver,...
        'minimum_distance', minimum_distance};

    % reconstruct the image iteratively
    [img, out, ray_initial_angles] = reconstructIterativeSoundspeed(...
        tof_discrepancy(:), recon_grid, emitter.positions,...
        receiver.positions, water, mask_raytracing, mask_reconst,...
        emitter.rotation_indices, sound_speed_ground_truth,...
        reconst_args{:});
    
else
    
    img = [];
    out = [];
    ray_initial_angles = [];
end

end