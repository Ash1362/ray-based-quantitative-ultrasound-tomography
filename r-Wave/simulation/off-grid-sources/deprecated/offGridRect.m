function mask = offGridRect(kgrid, rect_pos, Lx, Ly, theta, varargin)
%OFFGRIDRECT Generate a non-binary mask for a rectangular source on a 2D/3D grid.
%
% DESCRIPTION:
%     offGridRect computes a non-binary mask for implementing a rectangular
%     source in 2D or 3D simulations. It evenly samples the rectangle, and
%     for each sample point computes a band-limited interpolant
%     corresponding to a point source at that location. These point source 
%     responses are summed and scaled to give the source mask.
%
% USAGE:
%     mask = offGridRect(kgrid, rect_pos, Lx, Ly)
%     mask = offGridRect(kgrid, rect_pos, Lx, Ly, theta)
%     mask = offGridRect(kgrid, rect_pos, Lx, Ly, [], ...)
%     mask = offGridRect(kgrid, rect_pos, Lx, Ly, theta, ...)
%
% INPUTS:
%     kgrid             - Object of the kWaveGrid class defining the
%                         Cartesian and k-space grid fields.
%     rect_pos          - Cartesian position of the centre of the rectangle
%                         given as a two (2D) or three (3D) element vector
%                         [m].
%     Lx                - Height of the rectangle (along the x-coordinate
%                         before rotation) [m].
%     Ly                - Width of the rectangle (along the y-coordinate
%                         before rotation) [m].
%     theta             - Either a scalar (2D) or three element vector 
%                         [tx, ty, tz] specifying the orientation of the
%                         rectangle [deg]. In 3D, the orientation is
%                         specified by yaw, pitch, and roll (intrinsic
%                         rotation axes).
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
%     'Plot'            - Boolean controlling whether the disc sampling
%                         points are plotted (default = false).
%     'UpsamplingRate'  - Oversampling used to distribute the off-grid
%                         points compared to the equivalent number of
%                         on-grid points (default = 4).
%
% OUTPUTS:
%     mask              - 2D/3D non-binary source mask for a rectangle
%
% ABOUT:
%     author            - Elliott Wise
%     date              - 30th April 2018
%     last update       - 2nd July 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018 Elliott Wise
%
% See also makeCartRect

% define defaults
plot_rect       = false;
bli_tolerance   = 0.1;
bli_type = 'sinc';
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
                plot_rect = varargin{input_index + 1};
                
                % check value
                if ~islogical(plot_rect)
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

% check for theta input
if nargin < 5
    theta = [];
end

% compute the number of source points and subsequent scaling factor (note,
% this assumes dx = dy = dz)
area = Lx .* Ly;
N_ongrid  = area ./ (kgrid.dx.^2);
    
% compute the number of off-grid sources using the upsampling rate
N_offgrid = ceil(N_ongrid .* upsampling_rate);

% compute a set of uniformly spaced Cartesian points covering the disc
if isempty(theta)
    rect_points = makeCartRect(rect_pos, Lx, Ly, [], N_offgrid);
else
    rect_points = makeCartRect(rect_pos, Lx, Ly, theta, N_offgrid);
end

% note adjusted number of Cartesian points, and compute scaling factor
N_offgrid = size(rect_points, 2);
scale = N_ongrid ./ N_offgrid;

% remove rectangle points which are outside the simulation
rect_points = trimCartPoints(kgrid, rect_points);

% generate a source mask using BLIs centered on the Cartesian disc points
mask = offGridPoints(kgrid, rect_points, scale, ...
    'BLITolerance', bli_tolerance, 'BLIType', bli_type);  

% plot source points if required
if plot_rect
    figure;
    switch kgrid.dim
        case 2
            
            % plot source points
            plot(rect_points(2, :), rect_points(1, :), 'k.');
            axis(0.5 .* [kgrid.y_size .* [-1, 1], kgrid.x_size .* [-1, 1]]);
            axis equal;
            xlabel('y-position [m]');
            ylabel('x-position [m]');
            
        case 3
            
            % plot source points
            plot3(rect_points(2, :), rect_points(1, :), rect_points(3, :), 'k.')

            % adjust axes
            view(45, 30);
            axis equal;
            box on;
            grid on;
            axis(0.5 .* [kgrid.y_size .* [-1, 1], kgrid.x_size .* [-1, 1], kgrid.z_size .* [-1, 1]]);
            xlabel('y-position [m]');
            ylabel('x-position [m]');
            zlabel('z-position [m]');
            
    end
    title(sprintf('%d source points', size(rect_points, 2)))
end