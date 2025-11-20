function [mask, C] = calcOffGridParams(kgrid, tr_pos, z_offset, varargin)
%CALCOFFGRIDPARAMS calculates the interpolation parameters using
%off-grid sources toolbox
%
% DESCRIPTION:
%       calcOffGridParameters uses the off-Grid toolboc for constructing
%       the objects needed for mapping the pressure field from emitters to the grid points,
%       and conversely from grid points to the receiver. It constructs a binary mask that 
%       describes at which grid points the pressure field needs to be recorded
%       and a matrix C which contains the weights associated with the binary
%       mask 
%
% USAGE:
%       
%
% INPUTS:
%       kgrid    - k-Wave struct 
%       tr_pos   - 2/3 x N array of cartesian points containing the centers 
%                  of the transducers
%       z_offset - the off-set of the grid points in kgrid along
%                  z-directions with respect to the struct setting. For example,
%                  (-13,0)cm in setting is converted to (-6.5cm,+6.5cm) in the kgrid, because the
%                  coordinates in kgrid must be symmetric with respect to
%                  the origin by convention
%
% OPTIONAL INPUTS:
%       'Trans_Geom'      - the geomtery for the transducers. This can
%                           'point' (2D or 3D), 'line'(2D) or 'disk' (3D) 
%                           (default: 'point')
%       'Interp_Method'   - the method used for interpolation. This can 
%                           be 'nearest', 'linear', 'fourier' or 'offgrid'
%                           (default: 'offgrid')
%       'radius'          - radius of the disk for sen_type 'disc'
%       'bli_tolerance'   - the band-limitted interpolant threshold
%       'upsampling_rate' - upsampling factor when constructing
%                           'disc' type sensors

%
% OUTPUTS:
%        mask            - binary mask, whic is used in the k-Wave simulation 
%                           Using this, k-Wave will record the pressure time 
%                           series at a union of all grid points associated
%                           with all the transducers
%         C               - a sparse matrix that maps the time pressure series
%                           recorded at mask to the physical transducers
%                          
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
para.radius          = 1e-3;
para.upsampling_rate = 16;
para.disc_type = 'typical';


% get the number of dimensions
dim = size(tr_pos, 1);

switch dim
    case 2
        para.bli_tolerance = 0.001;
    case 3
        para.bli_tolerance = 0.1;
end




% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


% the number of transducers
num_trans = size(tr_pos, 2);

switch para.Trans_Geom
    case 'point'
        input_args_offGrid = {};
        tr.pos       =  tr_pos;
    case 'line'
        input_args_offGrid  = {'upsampling_rate', para.upsampling_rate};
        [tr.start_pos, tr.end_pos, ~] = lineEdgesCircularDetection(tr_pos, para.radius);
    case 'disk'
        input_args_offGrid  = {'upsampling_rate', para.upsampling_rate, 'disc_type', para.disc_type,...
            'UseSpiral', false};

        tr.pos         = tr_pos;
        tr.radius      = para.radius;
        tr.focus       = [0, 0, z_offset].';
    otherwise
        
end

% add the chosen BLI threshold to the optional inputs
input_args_offGrid  = {input_args_offGrid{:}, 'bli_tolerance', para.bli_tolerance};

[weight, ~, indices_in_union, union_indices] = calculateMask(kgrid, tr,...
    para.Trans_Geom, input_args_offGrid{:});

mask = zeros(size(kgrid.x));
mask(union_indices) = 1;

C_i = [];
for i = 1:num_trans
    C_i = [C_i; i*ones(length(weight{i}),1)];
end
C_j = cell2mat(indices_in_union);
C_v = cell2mat(weight);

% Construct a sparse matrix that includes the weights
C = sparse (C_i, C_j, C_v, num_trans, nnz(union_indices));

end