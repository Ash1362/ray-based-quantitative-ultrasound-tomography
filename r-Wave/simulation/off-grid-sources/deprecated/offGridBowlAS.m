function mask = offGridBowlAS(kgrid, bowl_axial_pos, radius, diameter, varargin)
%OFFGRIDBOWLAS Generate a non-binary mask for a bowl source on an axisymmetric (2D) grid.
%
% DESCRIPTION:
%     offGridBowlAS computes a non-binary mask for implementing a bowl
%     source in axisymmetric simulations. It evenly samples the bowl by
%     generating uniform points along an arc (including the negative radial
%     region), which then cover the bowl due to the axisymmetric
%     assumption. For each sample point, a band-limited interpolant is
%     computed corresponding to a point source at that location. These
%     point source responses are summed and scaled to give the source mask.
%     The beam axis is parallel to the x-coordinate axis (axial dimension),
%     and points in the positive direction.
%
% USAGE:
%     mask = offGridBowlAS(kgrid, bowl_axial_pos, radius, diameter)
%     mask = offGridBowlAS(kgrid, bowl_axial_pos, radius, diameter, ...)
%
% INPUTS:
%     kgrid             - Object of the kWaveGrid class defining the
%                         Cartesian and k-space grid fields.
%     bowl_axial_pos    - Axial position of rear surface of the bowl [m].
%     radius            - Radius of curvature of the bowl [m].
%     diameter          - Aperture diameter of the bowl [m].
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%  
%     'BLITolerance'    - Scalar value controlling where the spatial extent
%                         of the BLI at each point is trunctated as a
%                         portion of the maximum value (default = 0.1).
%     'BLIType'         - String controlling the BLI expression that is
%                         used for each point source, either 'sinc' or
%                         'exact' (default = 'sinc'). BLITolerance is
%                         ignored if 'exact' is specified.
%     'Plot'            - Boolean controlling whether the bowl is plotted
%                         (default = false).
%     'UpsamplingRate'  - Oversampling used to distribute the off-grid
%                         points compared to the equivalent number of
%                         on-grid points (default = 4).
%
% OUTPUTS:
%     mask              - Axisymmetric non-binary source mask for a bowl.
%
% ABOUT:
%     author            - Elliott Wise and Bradley Treeby
%     date              - 31st October 2017
%     last update       - 22nd October 2018
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2018 Elliott Wise and Bradley Treeby
%
% See also makeArc, makeBowl, makeCartArc, offGridArc, offGridBowl,
% offGridPoints

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
%
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
% more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% define defaults
plot_bowl       = false;
bli_tolerance   = 0.1;
bli_type        = 'sinc';
upsampling_rate = 4;

% assign optional input parameters if provided
if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'BLITolerance'
                
                % assign input
                bli_tolerance = varargin{input_index + 1};
                
                % check parameter range
                if ~isnumeric(bli_tolerance) || (bli_tolerance < 0) || (bli_tolerance > 1)
                    error('Optional input ''BLITolerance'' must be between 0 and 1.');
                end
                
            case 'BLIType'
                
                % assign input
                bli_type = varargin{input_index + 1};
                
                % check parameter
                if ~any(strcmp(bli_type, {'sinc', 'exact'}))
                    error('Optional input ''BLIType'' must be either ''sinc'' or ''exact''.')
                end
                            
            case 'Plot'
                
                % assign input
                plot_bowl = varargin{input_index + 1};
                
                % check value
                if ~islogical(plot_bowl)
                    error('Optional input ''Plot'' must be Boolean.');
                end
                
            case 'UpsamplingRate'
                
                % assign input
                upsampling_rate = varargin{input_index + 1};
                
            otherwise
                error('Unknown optional input.');
        end
    end
end

% check the bowl position is scalar
if numel(bowl_axial_pos) > 1
    error('Input bowl_axial_pos should be scalar');
end

% mirror the kWaveGrid about the symmetry axis
kgrid2 = kWaveGrid(kgrid.Nx, kgrid.dx, 2 * kgrid.Ny, kgrid.dy);

% assign non-uniform transformations if needed
if kgrid.nonuniform
    kgrid2.setNUGrid(1, kgrid.xn_vec, kgrid.dxudxn, kgrid.xn_vec_sgx, kgrid.dxudxn_sgx);
    yn_vec = kgrid.yn_vec;
    yn_vec_sgy = kgrid.yn_vec_sgy;
    yn_vec_m = mirror_y_vec(yn_vec);
    yn_vec_sgy_m = mirror_y_vec(yn_vec_sgy);
    yn_vec2 = [0; yn_vec_m(1:end-1); 1+yn_vec]/2;
    yn_vec_sgy2 = [yn_vec_sgy_m; 1+yn_vec_sgy]/2;
    dyudyn2 = [1, fliplr(kgrid.dyudyn(2:end)), kgrid.dyudyn];
    dyudyn_sgy2 = [fliplr(kgrid.dyudyn_sgy), kgrid.dyudyn_sgy];
    kgrid2.setNUGrid(2, yn_vec2, dyudyn2, yn_vec_sgy2, dyudyn_sgy2);
end

% set the beam vector to point in the positive x direction
focus_pos = bowl_axial_pos + 1;

% generate an off-grid line or arc source in the mirrored domain
if isinf(radius)
    start_point = [bowl_axial_pos, -diameter/2];
    end_point = [bowl_axial_pos, diameter/2];
    arc2 = offGridLine(kgrid2, start_point, end_point, ...
        'UpsamplingRate', upsampling_rate, 'BLITolerance', bli_tolerance, 'BLIType', bli_type, 'Plot', plot_bowl);
else
    arc2 = offGridArc(kgrid2, [bowl_axial_pos, 0], radius, diameter, [focus_pos, 0], ...
        'UpsamplingRate', upsampling_rate, 'BLITolerance', bli_tolerance, 'BLIType', bli_type, 'Plot', plot_bowl);
end

% use the plotting functionality of offGridArc, then truncate the radial
% axis
if plot_bowl
    xlim([0, kgrid.y_size]);
end

% keep points in the positive y domain
mask = arc2(:, kgrid.Ny + 1:end);
    
end

function ym = mirror_y_vec(y)
% Subfunction to mirror the input vector about the left side, i.e.,
% assuming boundaries are WSWA where left-endpoint of domain is included,
% right is not.

Ny = length(y);
s = (0:1:Ny - 1).' / Ny;
ss = (1:1:Ny).' / Ny;
ym = -flipud(y - s) + ss;

end