function mask = offGridArc(kgrid, arc_pos, radius, diameter, focus_pos, varargin)
%OFFGRIDARC Generate a non-binary mask for a arc source on a 2D grid.
%
% DESCRIPTION:
%     offGridArc computes a non-binary mask for implementing an arc source
%     in 2D simulations. It evenly samples the arc, and for each sample
%     point computes a band-limited interpolant corresponding to a point
%     source at that location. These point source responses are summed and
%     scaled to give the source mask.
%
% USAGE:
%     mask = offGridArc(kgrid, arc_pos, radius, diameter, focus_pos)
%     mask = offGridArc(kgrid, arc_pos, radius, diameter, focus_pos, plot_arc)
%
% INPUTS:
%     kgrid             - Object of the kWaveGrid class defining the
%                         Cartesian and k-space grid fields.
%     arc_pos           - Centre of the rear surface of the arc (midpoint)
%                         given as a two element vector [bx, by] [m].
%     radius            - Radius of curvature of the bowl [m].
%     diameter          - Aperture diameter (length of line connecting arc
%                         endpoints) [m].
%     focus_pos         - Any point on the beam axis of the arc given as a
%                         two element vector [fx, fy] [m].
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
%     'Plot'            - Boolean controlling whether the arc sampling
%                         points are plotted (default = false).
%     'UpsamplingRate'  - Oversampling used to distribute the off-grid
%                         points compared to the equivalent number of
%                         on-grid points (default = 2).
%
% OUTPUTS:
%     mask              - 2D non-binary source mask for an arc.
%
% ABOUT:
%     author            - Elliott Wise and Bradley Treeby
%     date              - 15th May 2017
%     last update       - 2nd May 2018
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2018 Elliott Wise and Bradley Treeby
%
% See also makeArc, makeCartArc

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
plot_arc        = false;
bli_tolerance   = 0.1;
bli_type        = 'sinc';
upsampling_rate = 2;

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
                plot_arc = varargin{input_index + 1};
                
                % check value
                if ~islogical(plot_arc)
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
    
% compute arc angle from chord: https://en.wikipedia.org/wiki/Chord_(geometry)
varphi_max = asin(diameter ./ (2 * radius));
    
% compute the number of source points and subsequent scaling factor (note,
% this assumes dx = dy; this should be updated to project the area
% onto the kgrid properly)
arc_length = 2 * radius * varphi_max;
N_ongrid  = arc_length ./ kgrid.dx;

% compute the number of off-grid sources using the upsampling rate
N_offgrid = ceil(N_ongrid .* upsampling_rate);

% compute scaling factor
scale = N_ongrid ./ N_offgrid;    
    
% compute a set of uniformly spaced Cartesian points covering the arc
arc_points = makeCartArc(arc_pos, radius, diameter, focus_pos, N_offgrid);
    
% remove points which are outside the simulation
arc_points = trimCartPoints(kgrid, arc_points);

% generate a source mask using BLIs centered on the Cartesian bowl points
mask = offGridPoints(kgrid, arc_points, scale, ...
    'BLITolerance', bli_tolerance, 'BLIType', bli_type);
    
% plot source points if required
if plot_arc
    
    % plot source points
    figure;
    plot(arc_points(2, :), arc_points(1, :), 'k.')

    % plot beam axis
    dp = focus_pos - arc_pos;
    hold on;
    quiver(arc_pos(2), arc_pos(1), dp(2), dp(1), 0, 'k')

    % adjust axes
    axis equal;
    box on;
    grid on;
    axis(0.5 * [kgrid.y_size .* [-1, 1], kgrid.x_size .* [-1, 1]]);
    xlabel('y-position [m]');
    ylabel('x-position [m]');
    set(gca, 'YDir', 'Reverse');
    title(sprintf('%d source points', size(arc_points, 2)));
    
end