function [cpl_mask, C] = constructSensorArrayUSCT(kgrid, pos, sen_type, interp_type, varargin)
%CONSTRUCTSENSORARRAYUSCT creates a description of the acoustic sensor array
%
% DESCRIPTION:
%       constructSensorArrayUSCT is a modified version of
%       constructSensorArray in Pammoth toolbox, and constructs the objects needed to model 
%       the signals recorded at a spatially extened sensor array. It 
%       constructs a binary mask that describes at which simulation voxels 
%       the pressure field needs to be recorded and a matrix C that maps 
%       the pressure time courses at these voxels to the pressure time courses
%       recorded at the physical sensor, e.g. by spatial averaging. 
%
% USAGE:
%       [cpl_mask, C] = constructSensorArray(pos, setting)
%       [cpl_mask, C] = constructSensorArray(pos, setting, sen_type)
%       [cpl_mask, C] = constructSensorArray(pos, setting, sen_type, ...
%                                                             interp_type)
%       [cpl_mask, C] = constructSensorArray(pos, setting, sen_type, ...
%                                                       interp_type, para)
%
% INPUTS:
%       kgrid       - a struct in k-Wave format for properties of the
%                     computational grid
%       pos         - 2/3 x N array of cartesian points containing the centers of
%                     the transducer objects
%       sen_type    - type of sensor used: 'point' for idealized point
%                     sensors, 'line' for line-like sensors (only 2D) whose normals
%                     point (0,0), 'disk' for disk-like sensors (in 3D!) that model spatially
%                     extended sensors (the disks' normals point to (0,0,0))
%
%       interp_type - type of spatial interpolation used: 'nearest' for
%                     simple nearst-neighbour interpolation, 'linear' for
%                     linear interpolation based on a triangulation of the
%                     domain and 'fourier' to use the off-grid source
%                     technique
%
% OPTIONAL INPUTS:    
%      'radius'          - radius of the disk for sen_type 'disk'
%       'tol'            - threshold tolerance for weights for interp_type
%                         'fourier' (default: 10^-10)
%      'N_mask'          - size of to spatial window to consider for 
%                         interp_type 'fourier' (default: 7)
%      'upsampling_rate' - upsampling factor when constructing
%                         'disk' type sensors
%
% OUTPUTS:
%      cpl_mask      - binary mask to be used in the k-Wave simulation. This is nonzro
%                      at a union of all grid points associated with all the transducers 
%         C          - a sparse matrix that maps the time pressure series
%                      f recorded at bmask to the physical transducers
%                      by g = C * f
%
% ABOUT:
%       author          - Felix Lucka, the 'fourier' interpolation is based 
%                         on code by Elliott Wise
%       date            - 15.12.2019
%       last update     - 15.12.2019
%       
%
% 
% This function is part of the r-Wave Toolbox (http://www.r-wave.org)
% Copyright (C) 2020 Ashkan Javaherian and Ben Cox
%%

para.radius          = 1e-3;
para.upsampling_rate = 4;
para.N_mask          = 7;
para.tol             = 10^-10;

% read additional arguments or overwrite default ones (no sanity check!)
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% check user defined value for sen_type, otherwise assign default value
if nargin < 3 || isempty(sen_type)
    sen_type = 'point';
end

% check user defined value for interp_type, otherwise assign default value
if nargin < 4 || isempty(interp_type)
    interp_type = 'nearest';
end

% check user defined value for para, otherwise assign default value
if nargin < 5 || isempty(para)
    para.tol = 10^-10;
    para.N_mask = 7;
    para.radius = 1e-3;
    para.upsampling_rate = 4;
end





% the number of grid points and grid spacing 
% along cartesian coordinates
Nx  = kgrid.Nx;
Ny  = kgrid.Ny;
dx  = kgrid.dx;
dim = kgrid.dim;


switch dim
    case 2
N    =  Nx*Ny;
Nxyz = [Nx,Ny];          
    case 3
Nz = kgrid.Nz;
N    = Nx*Ny*Nz;
Nxyz = [Nx,Ny,Nz];
end



% construct spatial grids and a fast function that returns the nearest grid
% point to a given point

node_numbers = reshape(1:N, Nxyz);
switch dim
    case 2
        cpl_mask    = false(Nx, Ny);
        [x, y] = ndgrid(kgrid.x_vec, kgrid.y_vec);
        grid_interp = griddedInterpolant(x , y , node_numbers, 'nearest');
    case 3
        cpl_mask    = false( Nx, Ny, Nz);
        [x, y, z] = ndgrid(kgrid.x_vec, kgrid.y_vec, kgrid.z_vec);
        grid_interp = griddedInterpolant(x, y, z, node_numbers, 'nearest');
end


% check for invalid inputs
switch sen_type
    case 'line'
        if(dim == 3)
            error('sen_type ''line''  only valid for dim = 2')  
        end 
    case 'disk'
        if(dim == 2)
            error('sen_type ''disk'' only valid for dim = 3')  
        end
end


n_sensors = size(pos, 2);
% construct differnt sensor types
switch sen_type
    case 'point'        
        switch interp_type
            case 'nearest'

                % construct cpl_mask
                indices           = grid_interp(pos');
                cpl_mask(indices) = true;
                cpl_mask_indices  = find(cpl_mask(:));
                
                % construct C
                [~, loc] = ismember(indices, cpl_mask_indices);
                C = sparse(1:n_sensors, loc, ones(1, n_sensors), n_sensors, ...
                            nnz(cpl_mask));
                
            case 'linear'
                                
                % extract the nearest nodes to each point
                indices       = grid_interp(pos');
                mask          = false( Nxyz);
                mask(indices) = true;
                % expand to 9/27 neighbourhood
                mask          = convn(mask, true(3,3,3), 'same') > 0;
                mask_indices  = find(mask(:));
                
                % compute triangulation and interpolation weights
                switch dim
                    case 2
                        [tri, bc] = gridDataFast2D(x(mask), y(mask), ...
                            pos(1,:), pos(2,:));
                    case 3
                        [tri, bc] = gridDataFast3D(x(mask), y(mask), ...
                            z(mask), pos(1,:), pos(2,:), pos(3,:));
                end
                tri                   = mask_indices(tri);
                cpl_mask(unique(tri)) = true;
                cpl_mask_indices      = find(cpl_mask(:));
                [~, loc]              = ismember(tri, cpl_mask_indices);
                C_i                   = repmat(1:n_sensors, dim + 1, 1);
                C_j                   = loc';
                C_v                   = bc';
                % construct C
                C = sparse(C_i(:), C_j(:), C_v(:), n_sensors, nnz(cpl_mask));
                
%             case 'fourier_cmp'
%                 
%                 tol = para.tol;
%                 
%                 % construct complete mask
%                 for i_sen = 1:n_sensors
%                     mask_i = offGridPoints(setting.kgrid, ...
%                         pos(:,i_sen) - setting.kgrid0Ind, 1);
%                     mask_i_support = abs(mask_i) > tol;
%                     cpl_mask = cpl_mask | mask_i_support;
%                 end
%                 
%                 cpl_mask_indices = find(cpl_mask(:));
%                 C_i = cell(n_sensors, 1);
%                 C_j = cell(n_sensors, 1);
%                 C_v = cell(n_sensors, 1);
%                 
%                 for i_sen = 1:n_sensors
%                     mask_i = offGridPoints(setting.kgrid, ...
%                         pos(:,i_sen) - setting.kgrid0Ind, 1);
%                     mask_i_support = abs(mask_i) > tol;
%                     
%                     C_i{i_sen} = i_sen * ones(nnz(mask_i_support),1);
%                     
%                     [tf, loc] = ismember(find(mask_i_support), cpl_mask_indices);
%                     if(any(not(tf)))
%                         error('single sensor mask support not part of the total sensor mask support')
%                     end
%                     C_j{i_sen} = loc;
%                     C_v{i_sen} = mask_i(mask_i_support);
%                 end
%                 
%                 C = sparse(cell2mat(C_i), cell2mat(C_j), cell2mat(C_v), n_sensors, nnz(cpl_mask));
                
            case 'fourier'
                
                tol     = para.tol;
                N_mask  = para.N_mask;
                win_1D  = -N_mask:N_mask;
                
                % extract the nearest nodes to each point and construct
                % spatial window functions
                switch dim
                    case 2
                        x_vec           = kgrid.x_vec;
                        y_vec           = kgrid.y_vec;
                        indices         = grid_interp(pos');
                        [ind_x, ind_y]  = ind2sub(Nxyz, indices);
                        [win_x, win_y]  = ndgrid(-N_mask:N_mask, -N_mask:N_mask);
                        win_ind         = [win_x(:), win_y(:)];
                        window_bool     = logical(abs(win_x.*win_y) <= N_mask);
                    case 3
                        x_vec                 = kgrid.x_vec;
                        y_vec                 = kgrid.y_vec;
                        z_vec                 = kgrid.z_vec;
                        indices               = grid_interp(pos');
                        [ind_x, ind_y, ind_z] = ind2sub(Nxyz, indices);
                        dum_z                 = zeros(1, 1, 2*N_mask+1);
                        [win_x, win_y, win_z] = ndgrid(-N_mask:N_mask, ...
                                            -N_mask:N_mask, -N_mask:N_mask);
                        win_ind     = [win_x(:), win_y(:), win_z(:)];
                        window_bool = logical(abs(win_x.*win_y.*win_z) <= N_mask);
                end
                
                cpl_mask(indices) = true;
                % expand mask by convolution with spatial window mask
                cpl_mask = convn(cpl_mask, window_bool, 'same') > 0;
                
                cpl_mask_indices = find(cpl_mask(:));
                
                % P implements a fast way to compute the ismember function
                % for searching for the indicies of a sub-set in a larger 
                % but fixed set 
                P                   = zeros(1,numel(x));
                P(cpl_mask_indices) = 1:length(cpl_mask_indices);
                
                % construct C
                C_i = cell(n_sensors, 1);
                C_j = cell(n_sensors, 1);
                C_v = cell(n_sensors, 1);
                
                for i_sen = 1:n_sensors
              
                    ind_x_i_sen = ind_x(i_sen) + win_1D;
                    ind_y_i_sen = ind_y(i_sen) + win_1D;
                    
                    % compute x part of the BLI
                    mask_i = sinc(pi*(x_vec(ind_x_i_sen) - pos(1,i_sen)) ./ dx);
                    
                    switch dim
                        case 2
                            
                            % compute y part of the BLI
                            mask_i         = mask_i * sinc(pi*(y_vec(ind_y_i_sen)' - pos(2,i_sen)) ./ dx);
                            
                            % extract the voxel indices where BLI is largest 
                            mask_i_support = window_bool & (abs(mask_i) > tol);
                            mask_indicies  = win_ind + [ind_x(i_sen), ind_y(i_sen)];
                            mask_indicies  = mask_indicies(mask_i_support,:);
                            mask_indicies  = sub2ind(Nxyz, mask_indicies(:,1), mask_indicies(:,2));
                            
                        case 3
                            
                            % compute y and z part of the BLI
                            ind_z_i_sen   = ind_z(i_sen) + win_1D;
                            mask_i        = mask_i * sinc(pi*(y_vec(ind_y_i_sen)' - pos(2,i_sen)) ./ dx);
                            dum_z(1:end)  = sinc(pi*(z_vec(ind_z_i_sen) - pos(3,i_sen)) ./ dx);
                            mask_i        = bsxfun(@times, mask_i, dum_z);
                            
                            % extract the voxel indices where BLI is largest 
                            mask_i_support = window_bool & (abs(mask_i) > tol);
                            mask_indicies  = win_ind + [ind_x(i_sen), ind_y(i_sen), ind_z(i_sen)];
                            mask_indicies  = mask_indicies(mask_i_support(:),:);
                            mask_indicies  = sub2ind(Nxyz, mask_indicies(:,1), mask_indicies(:,2), mask_indicies(:,3));
                    end
                    
                    
                    C_i{i_sen} = i_sen * ones(length(mask_indicies),1);
                    
                    %                     [tf, loc] = ismember(mask_indicies, cpl_mask_indices);
                    %                     if(any(not(tf)))
                    %                         error('single sensor mask support not part of the total sensor mask support')
                    %                     end
                    C_j{i_sen} =  P(mask_indicies)';
                    C_v{i_sen} = mask_i(mask_i_support);
                    
                end
                
                C_i = cell2mat(C_i);
                C_j = cell2mat(C_j);
                C_v = cell2mat(C_v);
                C   = sparse(C_i, C_j, C_v, n_sensors, nnz(cpl_mask));
                clear C_i C_j C_v
                
                % remove all-zero columns
                non_zeros = abs(full(sum(abs(C),1))) > 0;
                cpl_mask_indices = cpl_mask_indices(non_zeros);
                C = C(:, non_zeros);
                cpl_mask = false(Nxyz);
                cpl_mask(cpl_mask_indices) = true;
                
        end
    case {'disk', 'line'}
        
        switch sen_type
            case 'disk'
                % the disks are assumed to point to (0,0,0)!
                [ele_pts_tmpl, ~] = makeSpiral(kgrid, para.radius, para.upsampling_rate);
            case 'line'
                line_dx = dx / para.upsampling_rate;
                ele_pts_tmpl = -para.radius:line_dx:para.radius;
                ele_pts_tmpl = [ele_pts_tmpl; zeros(1, length(ele_pts_tmpl))];
        end
        n_ele_pts         = size(ele_pts_tmpl, 2);

        
        switch interp_type
            case 'nearest'

                % construct cpl_mask
                point_indices = zeros(n_ele_pts, n_sensors);
                for i_sen = 1:n_sensors
                    R = computeRotation(pos(:, i_sen));
                    ele_points = R * ele_pts_tmpl;
                    ele_points = bsxfun(@plus, ele_points, pos(:, i_sen));
                    point_indices(:, i_sen) = grid_interp(ele_points');
                end
                
                [cpl_mask_indices, ~, C_j] = unique(point_indices);
                cpl_mask(cpl_mask_indices) = true;
                
                C_i = repmat(1:n_sensors, n_ele_pts, 1);
                C_v = 1/n_ele_pts * ones(n_ele_pts * n_sensors,1);
                
                % construct C
                C = sparse(C_i, C_j, C_v, n_sensors, nnz(cpl_mask));
                
            case 'linear'
                
                % extract the nearest nodes to each point
                point_indices = zeros(n_ele_pts, n_sensors);
                for i_sen = 1:n_sensors
                    R = computeRotation(pos(:, i_sen));
                    ele_points = R * ele_pts_tmpl;
                    ele_points = bsxfun(@plus, ele_points, pos(:, i_sen));
                    point_indices(:, i_sen) = grid_interp(ele_points');
                end
                
                cpl_mask = false(Nxyz);
                cpl_mask(unique(point_indices)) = true;
                % expand to 9/27 neighbourhood
                cpl_mask = convn(cpl_mask, true(3,3,3), 'same') > 0;
                cpl_mask_indices = find(cpl_mask(:));
                % P implements a fast way to compute the ismember function
                % for searching for the indicies of a sub-set in a larger 
                % but fixed set 
                P = zeros(1,numel(x));
                P(cpl_mask_indices) = 1:length(cpl_mask_indices);
                
                % compute triangulation and interpolation weights
                switch dim
                    case 2
                        tri = delaunayTriangulation(x(cpl_mask), y(cpl_mask)); 
                    case 3
                        tri = delaunayTriangulation(x(cpl_mask), y(cpl_mask), z(cpl_mask)); 
                end
                
                % construct C
                C_j = zeros(n_ele_pts * (dim + 1), n_sensors);
                C_v = zeros(n_ele_pts * (dim + 1), n_sensors);
                
                %find the nearest triangle and the corresponding Barycentric coordinates
                for i_sen = 1:n_sensors
                    R = computeRotation(pos(:, i_sen));
                    ele_points = R * ele_pts_tmpl;
                    ele_points = bsxfun(@plus, ele_points, pos(:, i_sen));
                    [t, bc] = pointLocation(tri,ele_points');
                    loc = P(cpl_mask_indices(tri(t,:)));
                    C_j(:, i_sen) = loc(:);
                    C_v(:, i_sen) = bc(:) / n_ele_pts; 
                end
                C_i = repmat(1:n_sensors, n_ele_pts * (dim + 1), 1);
                
                
                C = sparse(C_i(:), C_j(:), C_v(:), n_sensors, nnz(cpl_mask));
                
                % remove all-zero columns
                non_zeros = abs(full(sum(abs(C),1))) > 0;
                cpl_mask_indices = cpl_mask_indices(non_zeros);
                C = C(:, non_zeros);
                cpl_mask = false(Nxyz);
                cpl_mask(cpl_mask_indices) = true;
            
            case 'fourier'
                notImpErr
        end 
    otherwise
        notImpErr
end