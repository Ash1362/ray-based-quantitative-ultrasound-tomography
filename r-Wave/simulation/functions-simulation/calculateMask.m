function [weight,point_indices,indices_in_union,union_indices] = calculateMask(kgrid, tr_pos, trans_geom, varargin)
% CALCULATEMASK calculates a non-binary mask and the corresponding
% indices for each of the transducers

% DESCRIPTION:
%       calculateMask calculates a non-binary mask and its associated
%       weights for the tranducers. These are used for mapping the pressure
%       field from the emitters to the grid points, and conversely from the
%       the grid points to the receivers

%
% USAGE:
%
%
%
% INPUTS:
%       kgrid             - a struct representing the computational grid in
%                           the k_Wave format
%       tr_pos            - For 'line' geometry, this is a struct arrray
%                           containing two fields 'start_pos' and 'end_pos',
%                           which are 2/3 x N arrays of Cartesian points
%                           containing the edges of the transducers.
%                           For disc geometry, this is a struct array
%                           containing two fields 'pos', 'radius', and 'focus',
%                           'pos' 2/3 x N arrays of Cartesian points containing
%                           the centre of the transducers, 'radius' a scalar value
%                           containing the radius of the discs, and 'focus'
%                           is the position of a point at which the discs are focused.
%                           For 'point' geometry', this is a struct array
%                           including the field 'pos', 2/3 x N arrays of
%                           Cartesian points containing the centre of
%                           the transducers,
%      trans_geom         - the geometry of transducer

%
% OPTIONAL INPUTS:
%       'upsampling_rate' - the rate of unsampling of the grid points for
%                           generation of off-Grid points
%
% OUTPUTS:
%       weight              - A cell array containing the weights allocated to
%                             the grid points for each transducer
%       point_indices       - A cell array containing the indices of the grid
%                             points associated with each transducer
%       indices_in_union    - A cell array containing the indices of the grid
%                             points associated with each transducer within a
%                             union of the grid points associated with all transducers
%       union_indices       - a union of indices of grid points associated with
%                             all transducers within the entire computational grid
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 15.12.2019
%       last update     - 15.12.2019
%
% This function is part of the r-Wave Toolbox
% Copyright (C) 2022 Ashkan Javaherian



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



% calculate the non-binary mask using off-grid toolbox
switch trans_geom
    case 'point'

        % get the optional input for computing BLI coefficients
        offgridpoint_args = {'BLITolerance', para.bli_tolerance};

        num_trans = size(tr_pos.pos, 2);
        point_indices = cell(num_trans, 1);
        weight        = cell(num_trans, 1);
        for kt = 1:num_trans
            nb_mask = offGridPoints(kgrid, tr_pos.pos(:,kt), 1, offgridpoint_args{:});
            point_indices{kt} = find(nb_mask);
            weight{kt} = nb_mask(point_indices{kt});
        end


    case 'line'


        num_trans = size(tr_pos.start_pos, 2);
        point_indices = cell(num_trans, 1);
        weight        = cell(num_trans, 1);
        % karray = kWaveArray;
        for kt = 1:num_trans
            %karray.addLineElement(tr_pos.start_pos(:, kt), tr_pos.end_pos(:, kt));
            %nb_mask = getArrayGridWeights(karray, kgrid);
            nb_mask = offGridLine(kgrid, tr_pos.start_pos(:, kt), tr_pos.end_pos(:, kt));
            point_indices{kt} = find(nb_mask);
            weight{kt} = nb_mask(point_indices{kt});
        end
    case 'disk'

        switch para.disc_type
            case 'typical'

                num_trans = size(tr_pos.pos, 2);
                point_indices = cell(num_trans, 1);
                weight        = cell(num_trans, 1);
                karray = kWaveArray('BLITolerance', para.bli_tolerance,...
                    'UpsamplingRate', para.upsampling_rate);
                for kt = 1:num_trans
                    karray.addDiscElement(tr_pos.pos(:, kt), 2*tr_pos.radius , tr.focus);
                    nb_mask = getArrayGridWeights(karray, kgrid);
                    point_indices{kt} = find(nb_mask);
                    weight{kt} = nb_mask(point_indices{kt});
                end

            case 'triangulation'
               
                num_trans = size(tr_pos.pos, 2);
                point_indices = cell(num_trans, 1);
                weight        = cell(num_trans, 1);
                input_disc = {'UpsamplingRate', para.upsampling_rate};
                karray = kWaveArray('BLITolerance', para.bli_tolerance);

                for kt = 1:num_trans

                    % get the integartion points over the disc surface and its corresponding
                    % scaling factors
                    [disc_points, scale_coeff] = makeDiscPoints(tr_pos.pos(:, kt), tr_pos.radius, ...
                        kgrid.dx, input_disc{:});

                    karray.addCustomElement(disc_points, scale_coeff, 2, 'Disc');
                    nb_mask = getArrayGridWeights(karray, kgrid);
                    point_indices{kt} = find(nb_mask);
                    weight{kt} = nb_mask(point_indices{kt});
                end

        end



end

% find indices of the union of grid points with nonzero weights using all the
% transducers
cpl_mask = false(size(kgrid.x));
cpl_mask(unique(cell2mat(point_indices))) = true;
union_indices = find(cpl_mask(:));

% allocatre a new index to the grid points within the chosen grid points
P = zeros(1, numel(kgrid.x));
P(union_indices) = 1:length(union_indices);


indices_in_union = cell(num_trans, 1);
for kt = 1:num_trans
    % allocate each point for each line an associated index stored in P
    indices_in_union{kt} = P(point_indices{kt})';
end


end