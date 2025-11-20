function [binary_mask, mask_entire_volume, interp_matrix, elapsed_time] = calcInterpParameters(kgrid,...
    transducer_position, z_offset, varargin)
%CALCINTERPPARAMETERS calculates the the interpolation parameters
%
% DESCRIPTION:
%       calcInterpParameters constructs the objects needed for mapping 
%       the pressure field from emitters to the grid points, and conversely from 
%       grid points to the receivers. It constructs a binary mask that 
%       describes at which grid points the pressure field needs to be recorded
%       and a matrix interpolation_matrix which contains the weights associated with the binary
%       mask
%
%
% USAGE:
%       
%
% INPUTS:
%       kgrid    - a k-Wave struct 
%       tr_pos   - 2/3 x N array of cartesian points containing the centers of
%                  the transducers
%       z_offset - the off-set of the grid points in kgrid along
%                  z-directions with respect to the struct setting.
%                  For example, the z-axis (-13,0) cm in setting is equlivalent
%                  to the axis(-6.5,+6.5) cm in the kgrid, because the
%                  coordinates in kgrid must be symmetric with respect to
%                  the origin by convention
%
% OPTIONAL INPUTS:
%       'Trans_Geom'      - the geomtery for the transducers. This can
%                           'point' (2D or 3D), 'line'(2D) or 'disk' (3D) (Default:'point')
%       'Interp_Method'   - the method used for interpolation. This can 
%                           be 'nearest', 'linear', 'fourier' or 'offgrid' (Default:'offgrid')
%       'simulate_on_grid' - a boolean controlling whether the binary mask is
%                              chosen the entire volume (area) inside the detection
%                              surface, or chosen only a volume
%                              encompassing the transducers
%       'radius_mask_extension' - the radius of the mask for storing the
%                                 pressure time series on the simulation grid,
%                                 if requested, (if simulate_on_grid =
%                                 true)
%       'radius'          - radius of the disk for Trans_Geom = 'disk'
%       'upsampling_rate' - upsampling factor when constructing transducers
%                           in the form 'disk' 
%       'N_mask'          - size of the spatial window to consider for 
%                           Interp_Method 'fourier' (Default = 7)
%       'tol'             - threshold tolerance for weights for interp_type
%                          'fourier' (Default = 10^-10)
%
% OUTPUTS:
%        binary_mask       - binary mask to be used in the k-Wave simulation. Using it,
%                           k-Wave will record the pressure time series as a union of
%                           all grid points associated with all the transducers
%        mask_entire_volume - binary mask to be used in the k-Wave
%                             simulation, which is the entire volume (area) 
%                             inside the detection surface (line)
%                             
%        interp_matrix    - a sparse matrix that maps the time pressure series
%                           recorded at binary_mask to the physical transducers
%                           by g = interp_matrix * f
%        elapsed_time     - elapsed time for calculation of the interpolation 
%                           parameters for all transducers
%
% ABOUT:
%       author           - Ashkan Javaherian
%       date             - 15.12.2019
%       last update      - 15.12.2019
%       
%
% This function is part of the r-Wave Toolbox
% Copyright (C) 2022 Ashkan Javaherian

para.Trans_Geom      = 'point';
para.Interp_Method   = 'offgrid';
para.simulate_on_grid = false;
para.radius_mask_extension = 1.05;
para.radius          = 1e-3;
para.upsampling_rate = 4;
para.N_mask          = 7;
para.tol             = 10^-10;

% read additional arguments or use default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% set the starting point for measuring the running time
ts = tic; 
switch para.Interp_Method
    case {'nearest'}
      
        %interp_args = {};
        % define the radius and upsampling rate for trasducers with finite area (volume)
        %if ~strcmp(para.Trans_Geom,'point')
        %    interp_args = {interp_args{:}, 'radius', para.radius, 'upsampling_rate', para.upsampling_rate};
        %end
        %if strcmp(para.Interp_Method,'fourier')
        %    interp_args = {interp_args{:}, 'N_mask', para.N_mask, 'tol', para.tol};
        %end
        [binary_mask, interp_matrix] = interpNeighborUST(kgrid, transducer_position);
    case{'offgrid'}
        interp_args = {'Trans_Geom',para.Trans_Geom, 'Interp_Method', para.Interp_Method,...
            'radius', para.radius, 'upsampling_rate', para.upsampling_rate};
        [binary_mask, interp_matrix] = calcOffGridParams(kgrid, transducer_position, z_offset, interp_args{:});       
end
        
if para.simulate_on_grid

switch kgrid.dim
    case 2
        mask_entire_volume = kgrid.x.^2 + kgrid.y.^2 <=...
            (para.radius_mask_extension * norm(transducer_position(:,1)))^2;
        % find the grid indices of the binary_mask
        indices = find(binary_mask);
        % find the grid indices of the mask inside the detection surface 
        indices_entire_volume = find(mask_entire_volume);
        loc = ismember(indices_entire_volume, indices);
        if nnz(loc) < length(indices)
            error(['the mask used for interpolation does not reside completely'...
                'on the mask used for storing the simulation on the grid.'...
                'The user must increase the radius of the mask used for siumlation on the grid'...
                'using the optional input for radius_mask_extension']);
        end
    case 3
        mask_entire_volume = [];
end

else
        mask_entire_volume = [];
end


% set the terminating point for measuring the running time        
elapsed_time = toc(ts);        
        
end