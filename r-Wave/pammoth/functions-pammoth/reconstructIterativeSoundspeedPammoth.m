function [img, out, initial_angles] = reconstructIterativeSoundspeedPammoth(tof_discrepancy,...
    recon_grid, emitter_positions, receiver_positions, water, mask_raytracing,...
    mask_reconst, rotation_indices, sound_speed_true, varargin)
%RECONSTRUCTITERATIVESOUNDSPEED reconstructs the image of sound speed iteratively
%
%
% DESCRIPTION:
%             reconstructIterativeSoundSpeed includes the iterative portion of
%             the rconstruction of the sound speed from time-of-flight
%             data. The objective function is iteratively linearised, and 
%             and the arising linear subproblems are solved  iteratively using
%             'sart', 'conjugate_gradient', or 'steepest_descent'.
%             
%      
%
% USAGE:
%     
%
% INPUTS:
%       tof_discrepancy    - the vector of the discrepancy of time-of-flights
%                            between the object-in-water and only-water data
%       recon_grid         - the grid for image reconstruction
%       emitter_positions  - the dim x num_emitter Cartsian position of
%                            emitters
%       receiver_positions - the dim x num_receiver Cartesian position of
%                            receivers
%       water              - a struct containing the sound speed (or refractive index)
%                            in water. This includes the fields: 
%       'sound_speed_only_water' - the sound speed in only water
%       'sound_speed_object_water' - the mean sound speed of water
%                            encompassing the object during the rotation of
%                            the detection surface
%       'refractive_object_water' - the vector of the refractive indices in 
%                            water encompassing the object during the rotation of
%                            the detection surface.
%       mask_raytracing    - a binary mask for ray tracing. Outside this binary
%                            mask, the gradient of the refractive index is
%                            set zero. Therefore, outside this binary mask,
%                            the rays are staright.
%       mask_reconst       - a binary mask for image reconstruction. Outside this 
%                            binary mask, the sound speed is set the sound speed
%                            in water. For the 'static' setting, i.e., when the 
%                            position of the transducers are fixed, it will be
%                            set the sound speed in only water. For the
%                            'rotating' setting, it will be set the mean of
%                            the sound speed ecncompassing the object
%                            during rotation of the detection surface.
%      rotation_indices    - 1 x num_emitter vector of indices of angluar
%                            positions for excitions
%      sound_speed_true    - the true sound speed map used as the ground
%                            truth for simulation studies
%  
%       
% OPTIONAL INPUTS:
%      'num_worker_pool'   - the number of workers for parallel programming
%      'emittersegments_to_nworkers' - The number of segments for dividing emitters
%                            to the number of cpu workers. For avoiding memory 
%                            leakage during the parallel programming, the 
%                            emitters are first divided to a number of segments,
%                            and then the segment are consecutively run over the threads.                        
%      'matrix_construction_method' - the method for construction of the
%                           system matrix and reconstructing the
%                           time-of-flight image. (default: 'bent-ray')
%      'linear_subproblem_method' - the method for solving the linearised
%                           subproblems, which can be either 'sart',
%                           'conjugate_gradient', or 'steepest_descent'.
%                           (Default : 'sart')
%      'linearisation_approach' - the approach for linearisation of the
%                           misfit function. This can be either 'difference',
%                           or 'absolute'. (default: 'difference')
%      'smoothing_window_size' - the size of the smoothing window applied
%                            on the updates of the refractive index.
%      'raylinking_method' - the method for ray linking. This can be either
%                            'Regula-Falsi' or 'Secant' for 2D case, or 
%                            'Quasi-Newton' for 3D case. 
%      'gridtoray_interp' - the method for interplating from grid points to
%                           rays and vice versa. This can be 'Bilinear', or
%                           'Bspline'.
%      'max_raylinking_iter' - the maximum number of iterations for ray
%                             linking
%      'raylinking_threshold' - the threshold distance [m] for stopping the
%                               ray linking iterations. If the distance
%                               between the interception point of the ray
%                               by the detection surface (ring) and the centre 
%                               of the receiver becomes less than
%                               raylinking_threshod, the ray linking
%                               algorithm will be stopped.
%      'raytogrid_spacing'    - the spacing of the sampled points along the ray
%                               to the grid spacing (Default: 1/2 (2D case),
%                               1 (3D case))
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
%       'minimum_distance'     - the mimimum distance for emitter-receiver
%                                pairs
%        'open_angle'          - the maximum angle between two vecors, one
%                                connecting the emitter to the centre of the
%                                detection surface (ring), and another connecting
%                                the emitter to the receivers. For each emitter, 
%                                the included receivers will be inside a cone
%                                with axis the vector connecting the emitter to
%                                the centre of the detection surface.
%        'num_iterout'         - the outer number of iterations
%        'num_iterin'          - the inner number of iterations
%        'step_length'         - the step length enforced on the updates of
%                                the sound speed after each linerisation
%
% OUTPUTS:
%      img                    - the reconstructed image of the sound speed
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
%       'relative_error'      - the L2 norm of the discrepancy of the
%                               updates of the sound speed and the ground
%                               truth over the L2 norm of the disrepancy of
%                               water and the ground truth, times 100.
%        'rmse'               - the rmse of the reconstructed images with
%                               respect to the ground truth.
%       initial_angles        - the initial angle of the linked (optimal) rays for all
%                               emitter-receiver pairs after ray linking
%                               applied to the last update of the sound
%                               speed
%                              
%    
%
% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 10.08.2022
%
% This script is part of the r-Wave Toolbox (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian and Ben Cox


para.num_worker_pool = 16;
para.matrix_construction_method  = 'bent-ray';
para.linearisation_approach = 'absolute';
para.linear_subproblem_method = 'conjugate_gradient';  
para.raylinking_threshold = 1e-6;
para.binaries_emitter_receiver = 'distances';


% get the dimensions
dim = size(emitter_positions, 1);

switch dim
    case 2
        
        para.raylinking_method ='Secant';
        para.gridtoray_interp = 'Bspline';
        para.max_raylinking_iter = 1000;
        para.raytogrid_spacing = 1/2;
        para.emittersegments_to_nworkers = 4;
        para.minimum_distance = 0.07;
        para.open_angle = pi/3;
        
    case 3
        
        para.raylinking_method ='Quasi-Newton';
        para.gridtoray_interp = 'Bilinear';
        para.max_raylinking_iter = 500;
        para.raytogrid_spacing = 1;
        para.emittersegments_to_nworkers = 32;
        para.minimum_distance = 0.16;
        para.open_angle = pi/4;
        
end


% choose the smoothing window size based on the grid spacing [m]
if recon_grid.dx < 1e-3
    error('The grid spacing smaller than 1mm is not necessary for a ray-based method.')
elseif recon_grid.dx < 1.5e-3
    para.smoothing_window_size = 7;
elseif recon_grid.dx < 2.01e-3
    para.smoothing_window_size = 5;
else
    error('The grid spacing must not be larger than 2mm.')
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
    
%if strcmp(para.linearisation_approach, 'difference')  && ... 
%    ~strcmp(para.linear_subproblem_method, 'steepest_descent')
%    error(['If the approach taken for linearisation is ''difference'', the approach for'...
%    'solving the linearised subproblem must be set ''steepest_descent''.'])
%end



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
                    para.num_iterin = 10;

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
                    
               case 'preconditioned_conjugate_gradient'    
                
                    % get the number of iterations for solving each
                    % linearised subproblem
                    para.num_iterin = 5;
                    
                    % choose a step length, a factor that is multiplied
                    % by the update direction solved from each linerised
                    % subproblem such that the convergence is ensured
                    para.step_length = 0.2;
                    
                    % get the number of iterations for solving the linearised problem for
                    % inverting the preconditioning matrix
                    para.preconditioning_iter = 100;
                    
                    % get the stooping tolerance for solving the linearised problem for
                    % inverting the preconditioning matrix
                    para.preconditioning_tol = 1e-3;
                    
                    % get the smoothing parameter for total-variation preconditioning
                    para.smoothing_parameter = 1e-3;
                    
                    % the scalar parameter added to the preconditioning matrix for ensuring
                    % the positive definiteness
                    para.gamma = 1e-2;   
                    
            end
        case 3
            switch   para.linear_subproblem_method
                
                case 'steepest_descent'
                    
                    % get the number of iterations for solving each
                    % linearised subproblem
                    para.num_iterin = 10;
                    
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
                    para.step_length = 0.2;
                    
                    
                case 'preconditioned_conjugate_gradient'
                    
                    % get the number of iterations for solving each
                    % linearised subproblem
                    para.num_iterin = 5;
                    
                    % choose a step length, a factor that is multiplied
                    % by the update direction solved from each linerised
                    % subproblem such that the convergence is ensured
                    para.step_length = 0.2;
                    
                    % get the number of iterations for solving the linearised problem for
                    % inverting the preconditioning matrix
                    para.preconditioning_iter = 100;
                    
                    % get the stooping tolerance for solving the linearised problem for
                    % inverting the preconditioning matrix
                    para.preconditioning_tol = 1e-3;
                    
                    % get the smoothing parameter for total-variation preconditioning
                    para.smoothing_parameter = 1e-3;
                    
                    % the scalar parameter added to the preconditioning matrix for ensuring
                    % the positive definiteness
                    para.gamma = 1e-2;
                       
            end
    end


if rem(para.smoothing_window_size, 2) == 0 && para.smoothing_window_size > 3
    error('The smoothing window size must not be an even value, if it is larger than 3.')
end

% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end
    


if iscell(receiver_positions)
    
    % if receiver_positions is given as a struct array, the transducers
    % rotate during the ust measurement.
    ust_setting = 'rotational';
      
else
    
    % if receiver_positions is given as a matrix array, the transducers
    % are fixed during the ust measurement.
    ust_setting = 'static';
    
end

% get the parameters for image reconstruction as scalars 
% get the outer number of iterations
num_iterout = para.num_iterout;

% get the inner number of iterations
num_iterin = para.num_iterin;

% get the step length
step_length = para.step_length;


% compute the vector of acoustic lengths in water between emitter-receiver pairs
% as the distance between emitter-receiver pairs
if strcmp(para.linearisation_approach, 'absolute')
acoustic_length_water = vectorise(calculateDistanceEmitterReceiver(...
    emitter_positions, receiver_positions, rotation_indices));
end

% display the approach for linearisaton of the objective function, and
% solving the linearised subproblems
disp(['The approach for linearisaton of the objective function is:' para.linearisation_approach])
disp(['The approach for solving the linearised paroblems is:' para.linear_subproblem_method])

if strcmp(ust_setting, 'static')
    
    % get the number of emitters
    num_emitter = size(emitter_positions, 2);
    
    % get the number of receivers
    num_receiver = size(receiver_positions, 2);

    % if the detection surface does not rotate
    water.sound_speed_object_water = water.sound_speed_only_water;
    water.refractive_object_water = ones(num_emitter * num_receiver, 1);
end

    if strcmp(para.linearisation_approach, 'absolute')
        
    % correct the in-water acoustic lengths for the changes in sound speed
    % of water during the rotation of the detection surface
    acoustic_length_water = water.refractive_object_water .*...
        acoustic_length_water;
    end

     % get the mean refractive index of water
    refractive_water_mean = water.sound_speed_only_water/water.sound_speed_object_water;
 
    % get the difference accoustic lengths
    % (here, the reference sound speed used for image reconstruction is used the
    % the same as the reference sound speed used for k-Wave simulation, but
    % they can be different.)
    acoustic_length = water.sound_speed_only_water * tof_discrepancy(:);
    
    % remove nans and zeros from difference TOFs
    acoustic_length(~isfinite(acoustic_length))= 0;
    zero_indices = ~acoustic_length;
    acoustic_length(zero_indices) = []; 
    
    switch para.linearisation_approach
        case 'absolute'
        acoustic_length_water(zero_indices) = [];
        case 'difference'
        water.refractive_object_water(zero_indices) = []; 
    end
    
    % get the ray spacing
    ray_spacing = para.raytogrid_spacing * recon_grid.dx;
    
    % allocate an struct for storing the outputs
    out = [];
      
    % a cell array containing the number of rays for each emitter-receiver
    out.num_rays = cell(1, num_iterout);
    
    % compute the system matrix using ray tracing for the homogeneous water,
    % which is the initial guess for TOF-based image reconstruction
    ray_args = {'nworker_pool', para.num_worker_pool,...
        'emittersegments_to_nworkers', para.emittersegments_to_nworkers,...
        'interp_method', para.gridtoray_interp,...
        'max_iter', para.max_raylinking_iter,...
        'raylinking_method', para.raylinking_method,...
        'raytogrid_spacing', para.raytogrid_spacing,...
        'binaries_emitter_receiver', para.binaries_emitter_receiver ,...
        'minimum_distance', para.minimum_distance,...
        'open_angle', para.open_angle,...
        'smoothing_window_size', para.smoothing_window_size,...
        'emittersegments_to_nworkers', para.emittersegments_to_nworkers,...
        'varepsilon', para.raylinking_threshold,...
        'refractive_background', refractive_water_mean};
    
    [system_matrix, ~, out.num_rays{1}, initial_angles, system_matrix_time(1)] = ...
        calcTransmissionMatrix(recon_grid, refractive_water_mean * ones(recon_grid.size), emitter_positions,...
        receiver_positions, [], mask_raytracing, rotation_indices, 0, ray_args{:});
    
    
    % remove the rows corresponding to the zeros in the difference TOFs
    system_matrix(zero_indices ,:) = [];
    
    % remove the zero indices from the vector for number of traced rays for ray
    % linking
    out.num_rays{1}(zero_indices) = [];
    
    % calculate the inverse coefficients for the rows (emitter-receiver pairs)
    coeff_res_inv = sum(system_matrix, 2);
    coeff_res_inv(~isfinite(coeff_res_inv)) = false;
    
    % remove the grid points which travel outside the binary mask for
    % ray tracing
    system_matrix(:, ~mask_raytracing(:)') = [];
    
    % get the binary vector for the rows which must be nulled
    res_null = coeff_res_inv < ray_spacing; % ~coeff_res_inv; %  
    
    % confine the rows of the system matrix to the emitter-receiver pairs with
    % rays having coefficients (for emitter-receiver pairs closer than a
    % minimum distance, the row of the system matrix is filled with full zeros.)
    system_matrix(res_null, :) = [];
    
    % get the absolute acoustic length
    switch para.linearisation_approach
        case 'difference'
            acoustic_length_absolute = acoustic_length(~res_null)...
                + water.refractive_object_water(~res_null) .*...
                (system_matrix * ones(size(system_matrix, 2), 1));
        case 'absolute'
            acoustic_length_absolute = acoustic_length + acoustic_length_water;
            acoustic_length_absolute(res_null) = [];
    end
    
    % parameters for image reconstruction using SART
    if strcmp(para.linear_subproblem_method, 'sart')
        
        % calculate the inverse coefficients for the columns (grid points)
        coeff_gridpoints_inv = sum(system_matrix, 1);
        
        % remove the grid points with zero inverse coefficient from the system matrix 
        system_matrix = system_matrix(:, coeff_gridpoints_inv > 0);
        
        % calculate the inverse coefficients for the grid points
        coeff_gridpoints = 1./coeff_gridpoints_inv(coeff_gridpoints_inv > 0); 
        
        % calculate the inverse coefficients for the residual
        coeff_res_inv(res_null) = [];
        coeff_res = 1./coeff_res_inv;
    
        step_length = struct('step_length', para.step_length, 'coeff_gridpoints',...
            coeff_gridpoints, 'coeff_res', coeff_res);
    end
    
    % initialise the refractive index distribution
    refractive = ones(recon_grid.size);
    
    % allocate variables for the results
    % residual norm
    out.residual_norm = zeros(1, num_iterout);
    
    % residual norm during solving the linear subproblem
    out.norm_res = zeros(num_iterout, num_iterin);
    
    if ~isempty(sound_speed_true)
        
    % relative error of the reconstructed images with respect to the ground
    % truth
    out.relative_error = zeros(1, num_iterout);
    
    % RMSE of the reconstructed images with respect to the ground truth
    out.rmse = zeros(1, num_iterout);
    
    end
    
    % a cell array containing the reconstructed images after each UST iteration
    out.imgs = cell(1, num_iterout);
    
    % a vector containing the cpu time for computation of the
    % system matrix
    out.system_matrix_time = [system_matrix_time, zeros(1, num_iterout-1)];
    
    % a vector containing the the cpu time for solving the linearised
    % subproblems
    out.update_time = zeros(1, num_iterout);
    
    % get the initial guess for the solution for the entire grid
    x_all = refractive_water_mean * ones(prod(recon_grid.size), 1);
    
    
    % get the previous norm of the residual
    norm_residual_previous = inf;
    
    
    
    for iterout = 1 : num_iterout
        
        if strcmp(para.raylinking_method, 'Regula-Falsi')
            
            % Using the Regula-Falsi approach, the initial guess for
            % initial unit direction of the ray is determined inside the 
            % function for ray tracing.
            initial_angles = [];
        end
        
        
        
        if strcmp(para.linear_subproblem_method, 'priorconditioned_conjugate_gradient')
            
            % define the handle function for computing the action of the
            % priorconditioning matrix on the residual vector
            priorconditioner = @(field) calcSpatialVariationMatrix(field,...
                refractive, mask_raytracing, recon_grid.size, para.smoothing_parameter,...
                para.gamma, 1);
            
        else
            
           priorconditioner = [];

        end
        
               
        if (strcmp(para.matrix_construction_method, 'bent-ray') && iterout > 1)...
                
        % update the system Matrix using the update of refractive index
        [system_matrix, ~, out.num_rays{iterout}, initial_angles,...
            out.system_matrix_time(iterout)] = calcTransmissionMatrix(recon_grid,...
            refractive, emitter_positions, receiver_positions, initial_angles,...
            mask_raytracing, rotation_indices, 0, ray_args{:});
        
        
        % remove the rows corresponding to the zeros in the difference TOFs
        system_matrix(zero_indices,:) = [];
        
        % remove the zero indices from the vector for number of traced rays for ray
        % linking
        out.num_rays{iterout}(zero_indices) = [];
    
        % calculate the inverse coefficients for the rows (emitter-receiver
        % pairs)
        coeff_res_inv = sum(system_matrix, 2);
        coeff_res_inv(~isfinite(coeff_res_inv)) = false;
        
        % confine the columns of the system matrix to grid points inside the
        % binary mask for ray tracing
        system_matrix(:, ~mask_raytracing(:)') = [];
        
        
        % get the binary vector for the rows which must be nulled
        res_null = coeff_res_inv < ray_spacing; % ~coeff_res_inv;  
        
        % confine the rows of the system matrix to the emitter-receiver pairs with
        % rays having nonzero coefficients, because for emitter-receiver pairs closer than a
        % minimum distance, the row of the system matrix is filled entirely with zeros.
        system_matrix(res_null, :)= [];
        
        % update the absolute acoustic length
        switch para.linearisation_approach
            case 'difference'
                acoustic_length_absolute = acoustic_length(~res_null)...
                    + water.refractive_object_water(~res_null) .*...
                    (system_matrix * ones(size(system_matrix, 2), 1));
            case 'absolute'
                acoustic_length_absolute = acoustic_length + acoustic_length_water;
                acoustic_length_absolute(res_null) = [];
        end

        % parameters for image reconstruction using SART
        if strcmp(para.linear_subproblem_method, 'sart')
            
            % calculate the inverse coefficients for the columns (grid points)
            coeff_gridpoints_inv = sum(system_matrix, 1);
            
            % remove the grid points with zero inverse coefficient from the system matrix
            system_matrix = system_matrix(:, coeff_gridpoints_inv > 0);
            
            % calculate the inverse coefficients for the grid points
            coeff_gridpoints = 1./coeff_gridpoints_inv(coeff_gridpoints_inv > 0); % (maskal)';
            %%
 
            % claculate the inverse coefficients for the residual
            coeff_res_inv(res_null) = [];
            coeff_res = 1./coeff_res_inv;
            
            % make an struct for the sart parameters
            step_length.coeff_gridpoints = coeff_gridpoints;
            step_length.coeff_res = coeff_res;
            
        end
        
        end
        
        
        if strcmp(para.linear_subproblem_method, 'sart')
            
            % get the update inside the mask for ray tracing
            x1 = refractive_water_mean * x_all(mask_raytracing(:));
            
            % get the update inside the mask and for grid points and nonzero
            % coefficients
            x = x1(coeff_gridpoints_inv > 0);
   
        else
            
            % get the update inside the mask for ray tracing
            x = x_all(mask_raytracing(:));
        end
        
  
        % get the norm of residual
        norm_residual = norm(acoustic_length_absolute - system_matrix * x);
          
        
        % display the system matrix
        disp(['The size of system matrix is:' num2str(size(system_matrix))])
        
        % display the norm of residual
        disp(['The norm of residual is:' num2str(norm_residual)])
        
       % check if the norm of residual decreases. Otherwise stop the algorithm. 
       if norm_residual < norm_residual_previous
           
        % solve the linearised problem
        [x, out.norm_res(iterout, :), out.residual_norm(iterout), out.update_time(iterout)] ...
                = solveTransmissionLinearSystem(system_matrix,...
                acoustic_length_absolute, x , num_iterin, step_length, para.linear_subproblem_method);
       else
           
           break;
           
       end
             
        
        % store the norm of residual for the next iteration
        norm_residual_previous = norm_residual;      
             
  
        % update the refractive index values for the grid points inside the
        % binary mask for ray tracing
        if strcmp(para.linear_subproblem_method, 'sart')
            
        x1(coeff_gridpoints_inv > 0) = x;
        
        % update the refractive index for the entire computational grid
        x_all(mask_raytracing(:)) = x1;
        
        else
            
         % update the refractive index for the entire computational grid
        x_all(mask_raytracing(:)) = x;
        
        end
        
        % make the refractive index ouside the smaller binary mask for
        % image recontruction 1
        x_all(~mask_reconst(:)) = refractive_water_mean;
        
        % update the refractive index
        refractive = reshape(x_all, recon_grid.size);
        
        % update the sound speed
        img = water.sound_speed_only_water./refractive;
        
        % save the sound speed update for the current linear subproblem
        out.imgs{iterout}.img = img;
        
        if ~isempty(sound_speed_true)
            
            % calculate the relative error
            out.relative_error(iterout)= 100*norm((img(mask_reconst)) - sound_speed_true(mask_reconst))/...
                norm(sound_speed_true(mask_reconst) - water.sound_speed_only_water);
            
            % calculate the RMSE
            out.rmse(iterout) = sqrt(sum((img(mask_reconst)...
                - sound_speed_true(mask_reconst)).^2) /nnz(mask_reconst) );

            % display the relative error for the current linear subproblem
            disp(['Relative Error:' num2str(out.relative_error(iterout), '%1.3f')]);
               
        end
        
        % display the number of linear subproblem
        disp(['The number of iteration is:'  num2str(iterout)]);
        
    end



end