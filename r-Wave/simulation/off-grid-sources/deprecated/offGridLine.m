function mask = offGridLine(kgrid, start_point, end_point, varargin)
%OFFGRIDLINE Create a non-binary mask for a line source.
%
% DESCRIPTION:
%     offGridLine computes a non-binary mask for implementing a line source
%     in simulations with any dimensionality. It evenly samples the line,
%     and for each sample point computes a band-limited interpolant
%     corresponding to a point source at that location. These point source
%     responses are summed and scaled to give the source mask.
%
% USAGE:
%     mask = offGridLine(kgrid, start_point, end_point)
%     mask = offGridLine(kgrid, start_point, end_point, plot_line)
%
% INPUTS:
%     kgrid             - Object of the kWaveGrid class defining the
%                         Cartesian and k-space grid fields.
%     start_point       - Start coordinate for the line given as a one
%                         (1D), two (2D), or three (3D) element vector [m].
%     end_point         - End coordinate for the line given as a one (1D),
%                         two (2D), or three (3D) element vector [m]. 
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
%     'Plot'            - Boolean controlling whether the line sampling
%                         points are plotted (default = false).
%     'UpsamplingRate'  - Oversampling used to distribute the off-grid
%                         points compared to the equivalent number of
%                         on-grid points (default = 2).
%
% OUTPUTS:
%     mask              - Non-binary source mask for a line.
%
% ABOUT:
%     author            - Elliott Wise and Bradley Treeby
%     date              - 29th March 2017
%     last update       - 2nd May 2018
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2018 Elliott Wise and Bradley Treeby
%
% See also makeLine, offGridArc

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
plot_line       = false;
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
                plot_line = varargin{input_index + 1};
                
                % check value
                if ~islogical(plot_line)
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

% compute the number of source points and subsequent scaling factor (note,
% this assumes dx = dy = dz)
line_length = norm(end_point - start_point);
N_ongrid  = line_length ./ kgrid.dx;

% compute the number of off-grid sources using the upsampling rate
N_offgrid = ceil(N_ongrid .* upsampling_rate);

% compute scaling factor
scale = N_ongrid ./ N_offgrid;

% distance between points in each dimension
d = (end_point - start_point) ./ N_offgrid;

% compute a set of uniformly spaced Cartesian points covering the line
switch kgrid.dim
    case 1
        line_points = linspace(start_point + d(1)/2, end_point - d(1)/2, N_offgrid);
    case 2
        px = linspace(start_point(1) + d(1)/2, end_point(1) - d(1)/2, N_offgrid);
        py = linspace(start_point(2) + d(2)/2, end_point(2) - d(2)/2, N_offgrid);
        line_points = [px; py];
    case 3
        px = linspace(start_point(1) + d(1)/2, end_point(1) - d(1)/2, N_offgrid);
        py = linspace(start_point(2) + d(2)/2, end_point(2) - d(2)/2, N_offgrid);
        pz = linspace(start_point(3) + d(3)/2, end_point(3) - d(3)/2, N_offgrid);
        line_points = [px; py; pz];
end

% generate a source mask using BLIs centered on the Cartesian line points
mask = offGridPoints(kgrid, line_points, scale, ...
    'BLITolerance', bli_tolerance, 'BLIType', bli_type);

% plot source points if required
if plot_line
    figure;
    switch kgrid.dim
        case 1
            plot(line_points, zeros(size(line_points)), 'k.');
            xlim(0.5 * kgrid.x_size * [-1, 1]);
            ylim([-1, 1]);
            set(gca, 'YTick', [], 'YTickLabel', {});
            xlabel('x-position [m]');
        case 2
            plot(line_points(2, :), line_points(1, :), 'k.');
            axis(0.5 * [kgrid.y_size .* [-1, 1], kgrid.x_size .* [-1, 1]]);
            set(gca, 'YDir', 'Reverse');
            xlabel('y-position [m]');
            ylabel('x-position [m]');
        case 3
            plot3(line_points(2, :), line_points(1, :), line_points(3, :), 'k.');
            axis(0.5 * [kgrid.y_size .* [-1, 1], kgrid.x_size .* [-1, 1], kgrid.z_size .* [-1, 1]]);
            xlabel('y-position [m]');
            ylabel('x-position [m]');
            zlabel('z-position [m]');
    end
    title(sprintf('%d source points', N_offgrid));
end