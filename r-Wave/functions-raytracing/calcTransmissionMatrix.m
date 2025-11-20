function [system_matrix, cartesian_position_endpoint, num_rays, optimal_polar_initial_direction,...
    matrix_construction_time] = calcTransmissionMatrix(recon_grid, refractive, cartesian_position_allemitters,...
    cartesian_position_allreceivers, polar_initial_direction, mask, rotation_indices, z_offset, varargin)
% CALCTRANSMISSIONMATRIX performs ray tracing on a given refractive index map, and
% constructs a sparse system matrix containing the contributions of the grid points
% from the traced rays.
% DESCRIPTION:
%       calcTransmissionMatrix performs ray tracing on a given refractive index map, and
%       constructs a sparse system matrix which contains num_emitter x
%       num_receiver rows and num_gridpoints columns. The matrix contains contains
%       the contributions of the grid points from the traced rays.
%
% USAGE:
%
%
% INPUTS:
%       recon_grid                     - the grid used for image reconstruction
%       refractive                     - the refractive index distribution
%       cartesian_position_allemitters - a dim x num_emitter matrix for
%                                        the Cartesian poition of all emitters
%       cartesian_position_allreceivers - a dim x num_receiver matrix for
%                                        the Cartesian position of all
%                                        receivers, if the data is
%                                        measured, if the geomtery for
%                                        poistion of the transducers is
%                                        'fixed'. Alternatively, a 1 x num_rotation
%                                        cell array, each containing the dim x
%                                        num_receiver cartesian position of
%                                        the receivers for each individual
%                                        rotation, if the geomtery for
%                                        position of the transducers is
%                                        set 'rotational'.
%       polar_initial_direction        - a 1 x num_emitter cell ray
%                                        containing the dim-1 x num_receiver
%                                        polar initial direction of the first
%                                        initial guess for the ray linking
%       mask                           - a binary mask for ray tracing. outside
%                                        this mask is assumed homogeneous
%                                        and water. Therefore, the rays travelling
%                                        outside the binary mask is traced as straight lines.
%      rotation_indices                - a 1 x num_emitter vector containing
%                                        the index of rotation for each emitter
%      z_offset                        - a scalar representing the discrepancy of
%                                        z-axis between the grid for image reconstruction
%                                        and an 'orgin-cenred' grid used
%                                        for k-Wave simulations, if
%                                        applied.
%
%
% OPTIONAL INPUTS:
%        'refrcative_background'- the reference refractive index
%        'binaries_emitter_receiver' - the method for choosing the emitter-receiver
%                                 pairs incorporated into the image reconstruction.
%                                 This can be set 'distances' or
%                                 'open_angle'. Using 'distances', the
%                                 emitter-receiver pairs with distance smaller than a threshold
%                                 are excluded. Using 'open_angle', the emitter-receiver pairs with
%                                 the angle between geometrical vectors emitter-centre
%                                 and emitter-receiver smaller than a threshold are excluded.
%                                 (useful for non-omnidirectional transducers)
%                                 (default = 'distances')
%        'minimum_distance'     - the minimum distance between emitter-receiver
%                                 pairs for which the rays are traced and are included
%                                 in the image reconstruction
%                                 (used if binaries_emitter_receiver =
%                                 'distances')
%        'open_angle'           - the minimum open angle of the emitter-receiver
%                                 pairs for being included in the image
%                                 reconstruction (used if binaries_emitter_receiver =
%                                 'open_angle')
%        'raytogrid_spacing'    - the ray-to-grid spacing, (default = 1
%                                 (2D), 1/2 (3D))
%       'raylink_method'     -  method for ray linking. For 2D case, this can be
%                              'Regula-Falsi' or 'Secant', and for 3D case, this
%                              can be 'Quasi-Newton'. 'Regula-Falsi' converges
%                              well with initial guess far from true, but it
%                              converges solwly. 'Secant' and 'Quasi-Newton'
%                              are fast, but converge badly for initial guesses
%                              far from true, and are therefore used
%                              through iteratively reconstruction of the
%                              sound speed, where the linked ray for each
%                              iteration is used as initial guess for ray
%                              linking for the next iteration. (Default:
%                              'Regula-Falsi (2D), and 'Quasi-Newton' (3D))
%       'raytracing_method'  - the method for ray tracing, which can be
%                              'Mixed-step', 'Dual-update',
%                              'Characteristics', or 'Runge-kutta-2nd'.
%                              (Default : 'Mixed-step')
%       'interp_method'      - method for interpolation, which can be
%                              'Bilinear' or 'Bspline'. (Default: 'Bilinear')
%        'max_iter'             - the maxium permissble number of iterations for ray
%                                 linking (default = 5000 (2D), 500 (3D))
%        'stopping_tolerance'   - the stopping tolerance for ray linking
%                                 (default = 1e-6)
%        'nWorkerPool'          - the number of used workers for paralel
%                                 programming (default = 8)
% OUTPUTS
%
%      system_matrix            - a sparse matrix of size with num_receiver x
%                                 num_emitter rows correponding the rays
%                                 linking all emitter-receiver pairs and
%                                 num_grid_points columns, where
%                                 num_grid_points is the number of the grid
%                                 points. Each row correponds to the ray
%                                 corresponding to a pair of
%                                 emitter-receivers and represents the
%                                 contribution of each grid point to the
%                                 ray.
%      cartesian_position_endpoint - the cartesian position of the end point of
%                                 the rays after ray linking. It should match the
%                                 position of the corresponding receiver.
%      num_rays                  - the number of traced rays for solving
%                                 the ray linking problem corresponding to
%                                 each emitter-receiver pair
%      optimal_polar_initial_direction - the optimal direction of the ray
%                                 in the polar coordinates after ray
%                                 linking
%      matrix_construction_time  - the cpu time for constructio of the
%                                 system matrix


% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 10.10.2022
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian

para = [];
para.nworker_pool = 16;
para.interp_method = 'Bilinear';
para.emittersegments_to_nworkers = 4;
para.smoothing_window_size = 7;
para.binaries_emitter_receiver = 'open_angle';
para.refractive_background = 1;
para.max_iter = 500;
para.raytracing_method = 'Mixed-step';

switch para.binaries_emitter_receiver
    case 'distances'
        para.minimum_distance = 0.08;
    case 'open_angle'
        para.open_angle = pi/4;
end

% get the dimension of the medium
dim = recon_grid.dim;

switch dim
    case 2
        para.raytogrid_spacing = 1/2;
        para.raylinking_method = 'Regula-Falsi';
    case 3
        para.raytogrid_spacing = 1;
        para.raylinking_method = 'Quasi-Newton';
end

para.varepsilon = 1e-6;


if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% the setting for position of the transducers: fixed or rotational
if iscell(cartesian_position_allreceivers)
    if isempty(rotation_indices)
        error('For a rotation setting, the rotation indices of the emitters must be specified.');
    end
    transducer_positions_setting = 'rotational';
else
    transducer_positions_setting = 'fixed';
end


% start the run time
start_time = tic;

%%=========================================================================
% GET THE DETECTION GEOMTERY
%==========================================================================
% allocate an empty variable for detection geometry
detec_geom = [];


if strcmp(transducer_positions_setting, 'fixed')
    
    % change the number of receivers from matrix to cell array
    cartesian_position_allreceivers = mat2cell(cartesian_position_allreceivers,...
        size(cartesian_position_allreceivers, 1),...
        size(cartesian_position_allreceivers, 2));
    
end

% get the number of the receivers
num_receiver = size(cartesian_position_allreceivers{1}, 2);

if num_receiver < 3
    error('The number of receivers for each position must be at least 3.')
end

% get the radius of the detection surface (if circle or hemi-sphere)
radius1 = norm(cartesian_position_allreceivers{1}(:, 1));
radius2 = norm(cartesian_position_allreceivers{1}(:, 2));
radius3 = norm(cartesian_position_allreceivers{1}(:, 3));

% check if the detection surface is circle or hemi-sphere
if abs(radius1 - radius2)< 1e-10  && abs(radius1 - radius3)< 1e-10
    
    % get the radius [m] of the detection surface, if it is circular or
    % hemi-spaherical
    detec_geom.radius_circle = radius1;
    
else
    
    if dim == 3
        
        % get the radius of the detection surface (if circle or hemi-sphere)
        radius1 = norm(cartesian_position_allreceivers{1}(1:2, 1));
        radius2 = norm(cartesian_position_allreceivers{1}(1:2, 2));
        radius3 = norm(cartesian_position_allreceivers{1}(1:2, 3));
        
    end
    
    if dim == 3  && abs(radius1 - radius2)< 1e-10   &&   abs(radius1 - radius3)< 1e-10
        
        % get the radius [m] of the detection surface, if it is cylinder
        % (2.5D)
        detec_geom.radius_cylinder = radius1;
        
    else
        
        % get the slope of the detection surface (if linear or planar)
        slope1 = (cartesian_position_allreceivers{1}(2, 2)- cartesian_position_allreceivers{1}(2, 1))./...
            (cartesian_position_allreceivers{1}(1, 2)- cartesian_position_allreceivers{1}(1, 1));
        
        slope2 = (cartesian_position_allreceivers{1}(2, end)- cartesian_position_allreceivers{1}(2, end-1))./...
            (cartesian_position_allreceivers{1}(1, end)- cartesian_position_allreceivers{1}(1, end-1));
        
        
        if abs(slope1-slope2) < 1e-10 || all(isinf([slope1, slope2]))
            
            
            % get the number of angular positions
            if isempty(rotation_indices)
                num_pos = 1;
            else
                num_pos = length(unique(rotation_indices));
            end
            
            % get the coefficients [a,b,c] for the line ax+by=c, where
            % c= -slope*x_0+y_0 with (x_0, y_0) a point on the line
            detec_geom = cell(num_pos, 1);
            

            for ind_pos = 1: num_pos
                
                % get the slope of the detection surface (if linear or planar)
                slope1 = (cartesian_position_allreceivers{ind_pos}(2, 2)- cartesian_position_allreceivers{ind_pos}(2, 1))./...
                    (cartesian_position_allreceivers{ind_pos}(1, 2)- cartesian_position_allreceivers{ind_pos}(1, 1));
                
                if isinf(slope1)
                    
                    % get the coefficient of the linear or planar
                    % detection surface along x-y plane in the form of
                    % ax+by=c with a=1, y=0, c=x_0
                    detec_geom{ind_pos}.line_coeff = [1, 0,...
                        cartesian_position_allreceivers{ind_pos}(1, 1)];
                    
                else
                    
                    % get the coefficient of the linear or planar
                    % detection surface along x-y plane in the form of
                    % ax + by = c with a = -slope, b = 1 and c=
                    % -slope * x_0 + y_0 with (x_0, y_0) a point on the line
                    detec_geom{ind_pos}.line_coeff = [-slope1, 1,...
                        -slope1 * cartesian_position_allreceivers{ind_pos}(1, 1) + ...
                        cartesian_position_allreceivers{ind_pos}(2, 1)];
                end
                
            end
            
        end
        
    end
    
end

if strcmp(transducer_positions_setting, 'fixed')
    
    % change the number of receivers from cell array to matrix
    cartesian_position_allreceivers = cell2mat(cartesian_position_allreceivers);
    
    % change detec_geom from cell array to matrix
    if iscell(detec_geom)
        detec_geom = cell2mat(detec_geom);
    end
    
end

% get the grid spacing
grid_spacing = recon_grid.dx;

% the number of pixels along the Carstesian coordinates
grid_size = [recon_grid.Nx , recon_grid.Ny];

% the vector for the Cartesian coordinates
xvec = recon_grid.x_vec;
yvec = recon_grid.y_vec;

% the cartesian position of starting and the ending grid points
pos_grid_first  = [xvec(1); yvec(1)];
pos_grid_end = [xvec(end); yvec(end)];


% extend the grid size and Cartesian coordinates to 3D
if dim == 3
    grid_size = [grid_size, recon_grid.Nz];
    zvec = recon_grid.z_vec - z_offset;
    pos_grid_first = [pos_grid_first; zvec(1)];
    pos_grid_end  = [pos_grid_end ; zvec(end)];
else
    zvec = [];
end


% the number of grid points
num_gridpoints = prod(grid_size);

if length(unique(refractive(:)))>1
    
    %  smooth the refractive index distribution, if the medium is heterogeneous
    if para.smoothing_window_size > 3.99
        refractive = smoothField(refractive,...
            para.smoothing_window_size, para.refractive_background);
    else
        
        switch dim
            case 2
                refractive = imgaussfilt(refractive,...
                    para.smoothing_window_size);
            case 3
                refractive = imgaussfilt3(refractive,...
                    para.smoothing_window_size);
        end
    end
end

% calculate the directional gradient fields along each Cartesian coordinate
switch dim
    case 2
        [refractive_gradient_y, refractive_gradient_x] = gradient(refractive, grid_spacing, grid_spacing);
    case 3
        [refractive_gradient_y, refractive_gradient_x, refractive_gradient_z] = ...
            gradient(refractive, grid_spacing, grid_spacing, grid_spacing);
end


% set zero the gradient of refraction outside the support of binary mask
% so the ray will be straight line ouside the binary mask
refractive_gradient_x = mask .* refractive_gradient_x;
refractive_gradient_y = mask .* refractive_gradient_y;


if dim == 3
    refractive_gradient_z =  mask.* refractive_gradient_z;
else
    refractive_gradient_z = [];
end


% display the ray-to-grid spacing
disp(['The ray-to-grid spacing is:' num2str(para.raytogrid_spacing)]);

% choose a step size for solving the rays
% the grid spacing is set the same for all Cartesian coordinates.
ray_spacing = para.raytogrid_spacing * grid_spacing;

% allocate a struct for interpolation parameters
ray_interp_coeffs = [];

switch para.interp_method
    
    
    case 'Bilinear'
        
        % get the directional gradients of the refractive index
        ray_interp_coeffs.refractive_gradient_x = refractive_gradient_x;
        ray_interp_coeffs.refractive_gradient_y = refractive_gradient_y;
        ray_interp_coeffs.refractive_gradient_z = refractive_gradient_z;
        
    case 'Bspline'
        
        % get the indices of the grid points included in the
        % Bspline interpolation with respect to the ray's point
        % (off-grid point)
        indices_vec = [-1; 0; 1; 2];
        
        switch dim
            case 2
                ray_interp_coeffs.raytogrid_indices_x = vectorise(repmat(indices_vec, [1, 4]));
                ray_interp_coeffs.raytogrid_indices_y = vectorise(repmat(indices_vec.', [4, 1]));
                ray_interp_coeffs.raytogrid_indices_z = [];
            case 3
                ray_interp_coeffs.raytogrid_indices_x = vectorise(repmat(indices_vec, [1, 4, 4]));
                ray_interp_coeffs.raytogrid_indices_y = vectorise(repmat(indices_vec.', [4, 1, 4]));
                ray_interp_coeffs.raytogrid_indices_z = vectorise(repmat(permute(indices_vec, [2 3 1]),...
                    [4, 4, 1]));
        end
        
        
        % get the coefficient of the polynomials for interpolation of the refractive
        % index
        ray_interp_coeffs.raytogrid_coeff_matrix = 1/6 * [-1, 3,-3, 1;...
            3,-6, 0, 4;...
            -3, 3, 3, 1;...
            1, 0, 0, 0];
        
        % get the coefficient of the polynomials for interpolation of the first-oder
        % gradient of the refractive index
        ray_interp_coeffs.raytogrid_coeff_derivative_matrix = 1/(6*grid_spacing) * [-3, 6, -3;...
            9,-12, 0;...
            -9,  6, 3;...
            3,  0, 0];
        
end


% define a handle function for ray linking and construction of the system
% matrix. The system matrix is a sparse matrix of the ray-to-grid
% interpolation coefficients. The function is run for each emitter separately.
calc_interp_coefficients = @(cartesian_position_emitter, polar_direction_allreceivers,...
    polar_initial_direction_allreceivers, position_receivers, detec_geom) rayLink(ray_interp_coeffs,...
    refractive, cartesian_position_emitter, position_receivers, polar_direction_allreceivers,...
    polar_initial_direction_allreceivers, xvec, yvec, zvec, pos_grid_first, pos_grid_end,...
    grid_spacing, ray_spacing, grid_size, dim, detec_geom, mask, para);


% the number of emitters
num_emitter = size(cartesian_position_allemitters, 2);

% specify a number of segment as multiply of number of workers for parallel
% programming
num_emitter_segments = para.emittersegments_to_nworkers * para.nworker_pool;

% the number of emitters for each segment, except for the last segment
emitter_per_segment = floor(num_emitter/num_emitter_segments);

% the number of remaining emitters
num_emitter_rem = rem(num_emitter, num_emitter_segments);

% allocate a cell array for the system matrix
system_matrix = cell(num_emitter_segments, 1);

% allocate a cell array for the cartesian position of the end point of the
% linked (optimal) rays
cartesian_position_endpoint = cell(num_emitter_segments, 1);

if ~strcmp(para.raylinking_method, 'Regula-Falsi')
    % allocate a sub-cell array for the polar direction of geometerical vector
    % from the emitters to the end points of the optimal rays
    optimal_polar_initial_direction = cell(num_emitter_segments, 1);
end

% allocate a cell array for the maximum number of rays for ray linking
num_rays = cell(num_emitter_segments, 1);


parfor (ind_emitter_segment = 1 : num_emitter_segments, para.nworker_pool)
   %  for ind_emitter_segment = 1 : num_emitter_segments
    
    % extract the indices of emitters for the segment
    if ind_emitter_segment < num_emitter_rem + 1
        indices_emitters_segment = (ind_emitter_segment - 1) * (emitter_per_segment + 1) + 1 :...
            ind_emitter_segment * (emitter_per_segment + 1);
    else
        indices_emitters_segment = ( num_emitter_rem * (emitter_per_segment + 1))  + ...
            ((ind_emitter_segment - 1 - num_emitter_rem) * emitter_per_segment + 1 :...
            (ind_emitter_segment - num_emitter_rem) * emitter_per_segment);
    end
    
    % choose the cartesian position of emitters for the segment
    cartesian_position_emitter_segment = cartesian_position_allemitters(:,...
        indices_emitters_segment);
    
    % the number of emitters for the segment
    num_emitter_per_segment = size(cartesian_position_emitter_segment, 2);
    
    % allocate a sparse matrix for the current segment
    system_matrix_segment = sparse(num_emitter_per_segment * num_receiver, num_gridpoints);
    
    % allocate a sub-cell array for the carstesian position of the end point of
    % the optimal rays
    cartesian_position_endpoint_segment = cell(num_emitter_per_segment, 1);
    
    if ~strcmp(para.raylinking_method, 'Regula-Falsi')
        % allocate a sub-cell array for the polar direction of geometrical vector
        % from the emitters to the end points of the optimal rays
        optimal_polar_initial_direction_segment = cell(num_emitter_per_segment, 1);
    end
    % allocate a sub-cell array for the number of rays for reaching
    % sufficiently close to the receivers after ray linking
    num_rays_segment = cell(num_emitter_per_segment, 1);
    
    if strcmp(transducer_positions_setting, 'rotational')
        rotation_indices_segment = rotation_indices(indices_emitters_segment);
    end
    
    % get the polar initial directions for the current segment, if given
    if ~isempty(polar_initial_direction)
        polar_initial_direction_segment = polar_initial_direction(indices_emitters_segment);
    end
    
    
    for ind_emitter = 1:num_emitter_per_segment
        
        disp( ['Number of emitter:'  num2str(indices_emitters_segment(ind_emitter))] )
        
        % the cartesian position of the emitter
        cartesian_position_emitter = cartesian_position_emitter_segment(:, ind_emitter);
        
        
        if strcmp(transducer_positions_setting, 'rotational')
            
            % if the detection geometry is rotational, the matrix 'cartesian_position_allreceivers'
            % is a cell array containing the position of receivers for each rotatation angle
            
            % get the index of rotation for getting the position of receivers
            ind_rot = rotation_indices_segment(ind_emitter);
            
            % if the setting for position of the transducers is 'rotational', the matrix of position
            % of receivers is separate for each 'ind_rot' (index of rotation)
            cartesian_position_allreceivers_singleemitter = cartesian_position_allreceivers{ind_rot};
            
        else
            
            % if the setting for position of the transducers is 'fixed', the matrix of position of
            % receivers is a dim x num_receiver matrix, and is fixed for all emitters
            cartesian_position_allreceivers_singleemitter = cartesian_position_allreceivers;
            
        end
        
  
        % calculate the cartesian direction of geometrical vectors from the
        % current emitter to the receivers
        cartesian_direction_allreceivers = cartesian_position_allreceivers_singleemitter...
            - cartesian_position_emitter;
        

        switch dim
            case 2
                
                % calculate the polar direction from emitter to the receivers
                switch para.raylinking_method
                    case {'Secant'}
                        [polar_direction_allreceivers, ~] = cart2pol(cartesian_direction_allreceivers(1,:),...
                            cartesian_direction_allreceivers(2,:));
                    case 'Regula-Falsi'
                        
                        % using 'Regula Falsi' method, the angles are
                        % calculated with respect to a reference geomterical vector from
                        % emitter to the centre of the detection circle (the origin of
                        % the Cartesian coordinates)
                        polar_direction_allreceivers = zeros(dim-1, num_receiver);
                        for ind_receiver = 1:num_receiver
                            polar_direction_allreceivers(:, ind_receiver) = ...
                                calcDirectionalAngle([-cartesian_position_emitter;  0],...
                                [cartesian_direction_allreceivers(:, ind_receiver); 0]);
                        end
                end
                
                
            case 3
                switch para.raylinking_method
                    case {'Newton','Quasi-Newton'}
                        [azimuthal_angle, elevation_angle, ~ ] = cart2sph(cartesian_direction_allreceivers(1,:),...
                            cartesian_direction_allreceivers(2,:), cartesian_direction_allreceivers(3,:));
                        polar_direction_allreceivers = [azimuthal_angle; elevation_angle];
                    otherwise
                end
                
            otherwise
                
                error (' The dimension is not correct.');
                
        end
        
        % get the row indices of the sparse matrix for this segment
        row_indices = (ind_emitter-1)* num_receiver + 1: ind_emitter * num_receiver;
        
        % get the detec_geom, which defines the stopping criteria for the rays
        % based on the position of receivers, for the current emitter
        if iscell(detec_geom)
            
            % get the receiver geometry for the current emitter
            detec_geom_pos = detec_geom{ind_rot};
            
        else
            
            % get the fixed receiver geometry
            % note that the setting can be rotational, but detec_geom
            % be fixed.
            detec_geom_pos = detec_geom;
        end
        
        switch para.raylinking_method
            case 'Regula-Falsi'
                
                % calculate the system matrix
                [system_matrix_segment(row_indices, :), ~, cartesian_position_endpoint_segment{ind_emitter},...
                    num_rays_segment{ind_emitter}] = calc_interp_coefficients(...
                    cartesian_position_emitter, polar_direction_allreceivers, [],...
                    cartesian_position_allreceivers_singleemitter, detec_geom_pos);
            case {'Secant','Newton','Quasi-Newton'}
                
                % calculate the polar initial direction of the initial guess
                % for the rays
                if isempty(polar_initial_direction)
                    polar_initial_direction_allreceivers = polar_direction_allreceivers;
                else
                    polar_initial_direction_allreceivers = polar_initial_direction_segment{ind_emitter};
                end
                % calculate the system matrix
                [system_matrix_segment(row_indices, :), optimal_polar_initial_direction_segment{ind_emitter},...
                    cartesian_position_endpoint_segment{ind_emitter}, num_rays_segment{ind_emitter}] = calc_interp_coefficients(...
                    cartesian_position_emitter, polar_direction_allreceivers,....
                    polar_initial_direction_allreceivers, cartesian_position_allreceivers_singleemitter, detec_geom_pos);
                
        end
        
        % display the number of receivers for which the end point of the ray
        % does not match the position of the receivers using the maximum
        % permissible number of iteration.
        % note that the end point of a bad linked ray may be very close to the
        % receiver, also note that emitter and receiver are assumed emission and reception points,
        % respectively.
        disp(['The number of bad linkings:'  num2str(nnz(num_rays_segment{ind_emitter} > para.max_iter-1))])
    end
    
    %% ========================================================================
    % GET THE OUTPUTS
    % =========================================================================
    % get the system matrix for the current segment
    system_matrix{ind_emitter_segment} = system_matrix_segment;
    
    % get the cartesian position of the end points for the current segment
    cartesian_position_endpoint{ind_emitter_segment} = cartesian_position_endpoint_segment;
    
    if ~strcmp(para.raylinking_method, 'Regula-Falsi')
        % get the polar direction of geometerical vector from the emitters to the end points
        % of the optimal rays for the current segment
        optimal_polar_initial_direction{ind_emitter_segment} = optimal_polar_initial_direction_segment;
    end
    
    % get the number of the linked rays for the current segment
    num_rays{ind_emitter_segment} = num_rays_segment;
    
    
end

% convert the system matrix as a cell with each element for each segment
% to a single sparse matrix
system_matrix = cat(1, system_matrix{:});

if strcmp(para.raylinking_method, 'Regula-Falsi')
    optimal_polar_initial_direction = [];
else
    
    % convert the polar direction of geometrical vectors from emitters to the end points of the linked rays
    % as a cell with each element for each segment to a cell with each element for each emitter
    optimal_polar_initial_direction = cat(1, optimal_polar_initial_direction{:});
end

% convert the Cartesian position of the end points of the linked rays as a cell with each element for each segment
% to a cell with each element for each emitter
cartesian_position_endpoint = cat(1, cartesian_position_endpoint{:});

% convert the number of rays as a cell with each element for each segment
% to a cell with each element for each emitter
num_rays = cat(1, num_rays{:});
% convert the cell for each emitter to a stack vector
num_rays = cat(1, num_rays{:});


if length(unique(refractive(:))) > 1
    disp('The system matrix was reconstructed using two-point ray tracing (ray linking)')
    disp(['The percentage fraction of bad raylinkings is'...
        num2str(nnz(num_rays>para.max_iter-1)/nnz(num_rays>1)*100, '%1.5e') '%'])
else
    disp('The system matrix was reconstructed using straight rays')
end

% the whole run time for construction of the system matrix
matrix_construction_time = toc(start_time);

end