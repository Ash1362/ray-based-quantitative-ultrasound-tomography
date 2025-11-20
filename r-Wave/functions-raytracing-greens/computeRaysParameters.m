function [cartesian_position_endpoint, cell_num_rays, optimal_polar_initial_direction,...
    cartesian_position_raypoints, ray_times, ray_absorption, rayspacing_receivers,...
    cartesian_position_auxiliary_left, cartesian_position_auxiliary_right,...
    adjoint_cartesian_position_auxiliary_left, adjoint_cartesian_position_auxiliary_right,...
    matrix_construction_time] =...
    computeRaysParameters(recon_grid, refractive, absorption_coeff, cartesian_position_allemitters,...
    cartesian_position_allreceivers, polar_initial_direction, mask, rotation_indices, z_offset, varargin)
%COMPUTERAYSPARAMETERS traces rays using ray linking and
%calculates the parameters required for approximating the pressure field
%using heterogeneous Green's funtion along the rays
%
% DESCRIPTION:
%       computeRaysParameters links the rays between emitter-receivers, and
%       calculates the parameters required for approximating the pressure field
%       using heterogeneous Green's function along the rays and on the receivers
%
% USAGE:
%
%
% INPUTS:
%       recon_grid                     - the grid used for image reconstruction
%       refractive                     - the refractive index [a.u.] distribution
%       absorption_coeff               - the absorption coefficient
%                                        [nepers (rad/s)^(-y) m^(-1)]
%       cartesian_position_allemitters - a dim x num_emitter matrix for
%                                        the Cartesian poition of all emitters
%       cartesian_position_allreceivers - a dim x num_receiver matrix for
%                                        the Cartesian position of all
%                                        receivers, if the data is
%                                        simulated in a 'fixed' geomotery
%                                        for position of the transducers.
%                                        a 1 x num_rotation cell array each
%                                        containing the dim x num_receiver
%                                        Cartesian position of the receivers
%                                        for each individual rotation, if
%                                        the goemetry is 'rotational'
%       polar_initial_direction        - a 1 x num_emitter cell ray each
%                                        containing the dim-1 x num_receiver
%                                        polar initial direction of the first
%                                        initial guess for the initial
%                                        angle of the rays for ray linking
%       mask                           - a binary mask for ray tracing.
%                                        Outside the mask is assumed homgenoeus
%                                        water, and therefore, the rays are
%                                        assumed as straight lines. This may
%                                        be used for improving the efficiency
%                                        of ray tracing
%      rotation_indices                - a 1 x num_emitter vector containing
%                                         the index of rotation for each
%                                         emitter. An empty variable for the
%                                         fixed setting
%      z_offset                        - a scalar representing the discrepancy of
%                                        z-axis between the grid for image
%                                        reconstruction and an 'orgin-centred' grid, eg.
%                                        kgrid in k-Wave
%
%
% OPTIONAL INPUTS:
%        'raytogrid_spacing'    - the ray-to-grid spacing, (Default = 1
%                                 (2D), 1/2 (3D))
%        'interp_method         - the method for ray to grid interpolation,
%                                 (Default = 'Bilinear')
%        'raylinking_method'    - the method for ray linking,
%                                 (Default = 'Regula-Falsi' (2D),
%                                 'Quasi-Newton' (3D))
%        'raylinking_initialisation' - the method for initialisation of the
%                                      ray linking (Default = 'Local')
%        'auxiliary_ray'        - Boolean controlling whether two auxiliary ray
%                                 is traced (Default = true)
%        'auxiliary_method'     - the approach for tracing auxiliary rays, which will
%                                 be used for computing the geomterical attenuation.
%                                 This can be set 'paraxial' or 'angle_perturbation'
%                                 (default = 'paraxial').
%        'reference_angle'      - the angle between main and auxiliary
%                                 rays
%        'max_iter'             - the maxium permissble number of iterations for ray
%                                 linking (Default = 500)
%        'varepsilon'            - the stopping tolerance for ray linking
%                                 (Default = 1e-6)
%        'smoothing_window_size' - the smoothing window size applied on the
%                                  refractive index for ray tracing
%                                  For grid with 1mm grid spacing,
%                                  (Default = 7 (2D), Default = 3 (3D))
%                                  Our ray linking solver developed for 3D
%                                  case can deal with sharp gradients even
%                                  better than protype 2D ray linking
%                                  approaches.
%        'reference_sound_speed' - the sound speed in water
%        'dim'                   - the dimension of eth grid
%        'angular_frequency_centre' - the frequency for including the
%                                   dispersion effects in ray tracing
%        'absorption_power'       - a saclar representing the absorption
%                                   exponent factor used for including the
%                                   dispersion effects in ray tracing
%         'minimum_distance'      - the minimum distance [m] between emitter-receiver
%                                   pairs for which the rays are traced and are included
%                                   in the image reconstruction
%        'max_num_points_factor'  - max_num_points_factor is factor gretaer than 1
%                                   to accunt for refraction effects on an
%                                   increase in the number of rays' points.
%                                   the size of the vector (or matrix) used
%                                   for storing the accumulated information
%                                   along the rays will then be
%                                   guesses by"
%                                   max_num_points_factor x the expected
%                                   number of points along an equivalent
%                                   straight rays. This approach is more efficient 
%                                   than increasing the size of the
%                                   vector/matrix array for every points
%                                   along the ray.
%                                   
%        'nworkerPool'          - the number of used workers for paralel
%                                 programming (Default = 16)
% OUTPUTS
%
%      system_matrix            - a sparse matrix linking the refrcative
%                                 index on the grid points to the accumulated
%                                 acoustic length along the rays linking
%                                 the emitter-receiver pairs
%      cartesian_position_endpoint - the Cartesian position of the end
%                                    point of rays
%      cell_num_rays               - a num_emitter x 1 cell array each containing
%                                    1 x num_receiver vector containing the number of rays
%                                    for ray linking between each emitter-receiver pair
%      optimal_polar_initial_direction - a num_emitter x 1 cell array each containing
%                                    dim-1 x num_receiver matrix containing
%                                    the polar initial direction of the
%                                    linked (optimal) rays
%      cartesian_position_raypoints - the dim x num_emitter Cartesian position
%                                     of rays' points
%      ray_times                    - the accumulated time delays along the
%                                     linked rays
%      ray_absorption               - the accumulated acoustic absorption
%                                     along the linked rays
%      rayspacing_receivers         - the ray spacing for the end point of
%                                     the linked ray on the receiver
%     cartesian_position_auxiliary_left - the Cartesian position of the
%                                      left auxiliary ray for the forward field
%     cartesian_position_auxiliary_right - the Cartesian position of the
%                                      right auxiliary ray for the forward field
%     adjoint_cartesian_position_auxiliary_left  - the Cartesian position of the
%                                      left auxiliary ray for the adjoint field
%     adjoint_cartesian_position_auxiliary_right - the Cartesian position of the
%                                      right auxiliary ray for the adjoint field

%     matrix_construction_time    - the time of running the m-file function
%
%
%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian

para = [];
para.nworker_pool = 16;
para.interp_method = 'Bspline';
para.smoothing_window_size = 7;
para.reference_sound_speed = 1500;
para.do_perturb_initial_position = false;
para.auxiliary_ray = true;
para.max_num_points_factor = 1.3;
para.auxiliary_method = 'paraxial';
para.max_iter = 500;
para.absorption_power = 1.4;
para.sound_speed_reference = 1500;
para.angular_frequency_centre = 2 * pi * 1e6;
para.varepsilon = 1e-6;
para.smooth_medium = true;
para.refractive_background = 1;
para.do_direction_angle = true;
para.minimum_distance = 0.08;
% get the number of dimensions
dim = recon_grid.dim;

switch dim
    case 2
        
        para.raytogrid_spacing = 1/2;
        para.raylinking_method = 'Regula-Falsi';
        para.reference_angle = pi/(2*180);
        
    case 3
        
        para.raytogrid_spacing = 1;
        para.raylinking_method = 'Quasi-Newton';
        
end

if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% start the run time
start_time = tic;



if ~(para.auxiliary_ray  && strcmp(para.auxiliary_method, 'paraxial'))

    % the perturbed rays' directions in the polar coordinates because of
    % perturbation to the initial positions are not computed, if the
    % auxiliary rays are not computed, or the approach for tracing auxiliary
    % rays is not set 'paraxial'.
    para.do_direction_angle = false;
    
end

% the number of emitters
num_emitter = size(cartesian_position_allemitters, 2);

% the setting for position of the transducers: fixed or rotational
if iscell(cartesian_position_allreceivers)
    if isempty(rotation_indices)
        error('For a rotation setting, the rotation indices of the emitters must be specified.');
    end
    transducer_positions_setting = 'rotational';
else
    transducer_positions_setting = 'fixed';
end





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
pos_grid_first = [xvec(1); yvec(1)];
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

% calculate the accumulated acosutic absorption, if the absorption
% coefficient map is nonzero
do_absorption = nnz(absorption_coeff) > 0;

if do_absorption
    
    % convert dB MHz^{-y} cm^{-1} to nepers (rad/s)^{-y} m^{-1}
    absorption_coeff = db2neper(absorption_coeff, para.absorption_power);
    
    % get the nonsmoothed refractive index
    refractive_nonsmoothed = refractive .* (1 + absorption_coeff .* para.sound_speed_reference./refractive...
        .* tan(pi * para.absorption_power/2) * para.angular_frequency_centre^(para.absorption_power-1) );
    
else
    
    refractive_nonsmoothed = refractive;
    
end



if length(unique(refractive_nonsmoothed(:)))>1
    
    %  smooth the refractive index distribution, if the medium is heterogeneous
    if para.smoothing_window_size > 3.99
        refractive = smoothField(refractive_nonsmoothed,...
            para.smoothing_window_size, para.refractive_background);
    else
        
        switch dim
            case 2
                refractive = imgaussfilt(refractive_nonsmoothed,...
                    para.smoothing_window_size);
            case 3
                refractive = imgaussfilt3(refractive_nonsmoothed,...
                    para.smoothing_window_size);
        end
    end
    
  %  refractive = ones(grid_size);
end



% calculate the directional gradient fields along each Cartesian
% coordinate
switch dim
    case 2
        [refractive_gradient_y, refractive_gradient_x] = gradient(refractive,...
            grid_spacing, grid_spacing);
        
    case 3
        [refractive_gradient_y, refractive_gradient_x, refractive_gradient_z] = ...
            gradient(refractive, grid_spacing, grid_spacing, grid_spacing);
end


% set zero the gradient of refraction outside the support of binary mask
% so the ray will be straight line ouside the binary mask
refractive_gradient_x = mask.* refractive_gradient_x;
refractive_gradient_y = mask.* refractive_gradient_y;


if dim == 3
    refractive_gradient_z = mask .* refractive_gradient_z;
else
    refractive_gradient_z = [];
end


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
        
        
        % get the coefficient of the polynomials for interpolation of the first-order
        % gradient of the refractive index
        ray_interp_coeffs .raytogrid_coeff_derivative_matrix = 1/(6*grid_spacing) * [-3, 6, -3;...
            9,-12, 0;...
            -9,  6, 3;...
            3,  0, 0];
        
        
        % get the coefficient of the polynomials for interpolation of the second-order
        % gradient of the refractive index
        ray_interp_coeffs.raytogrid_coeff_second_derivative_matrix = 1/(6*grid_spacing^2) *[-6, 6;...
            18,-12;...
            -18, 6;...
            6, 0];
        
end

% estimate the maximum number of rays' points for forming the matrices for
% storing the accumulated parameters along the rays
if isfield(detec_geom, 'radius_circle')
    num_points = ceil(2 * para.max_num_points_factor * detec_geom.radius_circle/ray_spacing) ;
elseif isfield(detec_geom, 'radius_cylinder')
    num_points = ceil(2 * para.max_num_points_factor * detec_geom.radius_cylinder/ray_spacing);
else
    error('Not implemented yet!')
end
    

% allocate a cell array for the optimal polar initial directions, which are
% the polar initial directions of the linked rays
optimal_polar_initial_direction = cell(1, num_emitter);

% allocate a cell array for the cartesian position of the end point of the
% linked (optimal) rays
cell_cartesian_position_endpoint = cell(1, num_emitter);

% allocate a cell array for the maximum number of rays for ray linking
cell_num_rays = cell(1, num_emitter);

% allocate a cell array for the Cartesian position of the points along the
% rays
cartesian_position_raypoints = cell(1, num_emitter);

% allocate a cell array for the time delays along the rays
ray_times = cell(1, num_emitter);

% allocate a cell array for the accumulated absorption along the rays
ray_absorption = cell(1, num_emitter);

% allocate a cell array for the Cartesian position of the forward auxiliary
% rays
% left
cartesian_position_auxiliary_left = cell(1, num_emitter);
%right
cartesian_position_auxiliary_right = cell(1, num_emitter);

% allocate a cell array for the Cartesian position of the adjoint auxiliary
% rays
% left
adjoint_cartesian_position_auxiliary_left = cell(1, num_emitter);
% right
adjoint_cartesian_position_auxiliary_right = cell(1, num_emitter);

% allocate a cell array for the last spacing along the ray
rayspacing_receivers = cell(1, num_emitter);

% define ahandle function for ray linking and construction of the system
% matrix for each single emitter
calc_parameters = @(cartesian_position_emitter, polar_direction_allreceivers,...
    polar_initial_direction_allreceivers, position_receivers, detec_geom) rayLinkFullWave(ray_interp_coeffs,...
    refractive, refractive_nonsmoothed, absorption_coeff, cartesian_position_emitter,...
    position_receivers, polar_direction_allreceivers,...
    polar_initial_direction_allreceivers, xvec, yvec, zvec, pos_grid_first, pos_grid_end,...
    grid_spacing, ray_spacing, grid_size, dim, detec_geom, mask, num_points, num_emitter, para);


 parfor (ind_emitter = 1 : num_emitter, para.nworker_pool)
       % for ind_emitter = 1:num_emitter % (for test)
       
    disp( ['Number of emitter:'  num2str(ind_emitter)] )
    
    % the cartesian position of the emitter
    cartesian_position_emitter = cartesian_position_allemitters(:, ind_emitter);
    
    if iscell(cartesian_position_allreceivers)
        
        % if the detection geometry is rotational, the matrix 'cartesian_position_allreceivers'
        % is a cell array containing the position of receivers for each individual emitter,
        % a dim x num_receiver matrix for each individual emitter
        
        if ~isempty(rotation_indices)
            
            % get the index of rotation associated with the emitter
            ind_rot = rotation_indices(ind_emitter);
        else
            error('For a rotation setting, the rotation indices of the emitters must be specified.');
        end
        
        % if the detection geometry is 'rotational', the matrix of position
        % of receivers is separate for each 'ind_rot' (index of rotation)
        cartesian_position_allreceivers_singleemitter = cartesian_position_allreceivers{ind_rot};
        
    else
        
        % if the detection geometry is 'fixed', the matrix of position of
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
                    
                    % using 'Rregula Falsi' method, the angles are
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
    
    
    switch para.raylinking_method
        case 'Regula-Falsi'
            
            polar_initial_direction_allreceivers = [];
            
        case {'Secant','Newton','Quasi-Newton'}
            % calculate the polar initial direction of the initial guess
            % for the rays
            if isempty(polar_initial_direction)
                polar_initial_direction_allreceivers = polar_direction_allreceivers;
            else
                polar_initial_direction_allreceivers = polar_initial_direction{ind_emitter};
            end
            
    end
    
    
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
    
    
    % get the rays' parameters
    [optimal_polar_initial_direction{ind_emitter},...
        cell_cartesian_position_endpoint{ind_emitter}, cell_num_rays{ind_emitter},...
        cartesian_position_raypoints{ind_emitter}, acoustic_length_emitter, absorption_emitter, rayspacing_receivers{ind_emitter},...
        cartesian_position_auxiliary_left{ind_emitter}, cartesian_position_auxiliary_right{ind_emitter},...
        adjoint_cartesian_position_auxiliary_left{ind_emitter}, adjoint_cartesian_position_auxiliary_right{ind_emitter}]...
        = feval(calc_parameters, cartesian_position_emitter, polar_direction_allreceivers,....
        polar_initial_direction_allreceivers, cartesian_position_allreceivers_singleemitter,...
        detec_geom_pos);
    
    % calculate the time delays along the ray from the acoustic
    % lengths divided by the reference (water) sound speed
    ray_times{ind_emitter} = acoustic_length_emitter/para.reference_sound_speed;
    
    if do_absorption
        
        % get the accumulated acoustic absorption along the ray
        ray_absorption{ind_emitter} = absorption_emitter;
        
    else
        
        % set the accumulated absorption on the rays zero
        ray_absorption{ind_emitter} = 0;
        
    end
    
    % display the number of receivers for which the end point of the ray
    % does not match the position of the receivers uisng the maximum permissible numer of iteration
    % note that the end point of a bad linked ray may be very close to the
    % receiver. also note that emitter and receiver are assumed emission and reception points,
    % respectively.
    disp(['The number of bad linkings:'  num2str(nnz(cell_num_rays{ind_emitter} > para.max_iter-1))])
    
end


if strcmp(para.raylinking_method, 'Regula-Falsi')
    optimal_polar_initial_direction = [];
end


cartesian_position_endpoint = cell2mat(cell_cartesian_position_endpoint);
clear cell_cartesian_position_endpoint


% the whole run time for construction of the system matrix
matrix_construction_time = toc(start_time);

end