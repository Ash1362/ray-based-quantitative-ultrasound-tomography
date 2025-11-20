function mask = offGridDisc(kgrid, disc_pos, diameter, focus_pos, varargin)
%OFFGRIDDISC Generate a non-binary mask for a disc source on a 2D/3D grid.
%
% DESCRIPTION:
%     offGridDisc computes a non-binary mask for implementing a disc source
%     in 2D or 3D simulations. It evenly samples the disc using Fermat's
%     spiral, and for each sample point computes a band-limited interpolant
%     corresponding to a point source at that location. These point source
%     responses are summed and scaled to give the source mask.
%
% USAGE:
%     mask = offGridDisc(kgrid, disc_pos, diameter)
%     mask = offGridDisc(kgrid, disc_pos, diameter, focus_pos)
%     mask = offGridDisc(kgrid, disc_pos, diameter, [], ...)
%     mask = offGridDisc(kgrid, disc_pos, diameter, focus_pos, ...)
%
% INPUTS:
%     kgrid             - Object of the kWaveGrid class defining the
%                         Cartesian and k-space grid fields.
%     disc_pos          - Cartesian position of the centre of the disc
%                         given as a two (2D) or three (3D) element vector
%                         [m]. 
%     diameter          - Diameter of the disc [m].
%     focus_pos         - Any point on the beam axis of the disc given as a
%                         three element vector [m] (used for discs defined
%                         in 3D only).
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
%     'UseSpiral'       - Boolean controlling whether the disc sampling
%                         points are chosen using a spiral sampling pattern
%                         (default = false). Concentric sampling is used by
%                         default.
%
% OUTPUTS:
%     mask              - 2D/3D non-binary source mask for a disc
%
% ABOUT:
%     author            - Elliott Wise and Bradley Treeby
%     date              - 30th March 2017
%     last update       - 2nd July 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2019 Elliott Wise and Bradley Trebey
%
% See also makeDisc, makeCartDisc, offGridBowl

% define defaults
plot_disc       = false;
bli_tolerance   = 0.1;
bli_type        = 'sinc';
upsampling_rate = 4;
use_spiral      = false;

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
                plot_disc = varargin{input_index + 1};
                
                % check value
                if ~islogical(plot_disc)
                    error('Optional input ''Plot'' must be Boolean.');
                end
                
            case 'UpsamplingRate'
                
                % assign input
                upsampling_rate = varargin{input_index + 1};
                
            case 'UseSpiral'
                
                % assign input
                use_spiral = varargin{input_index + 1};
                
                % check value
                if ~islogical(use_spiral)
                    error('Optional input ''UseSpiral'' must be Boolean.');
                end
                
            otherwise
                error('Unknown optional input.');
        end
    end
end

% check for focus_pos input
if (kgrid.dim == 2) || (nargin < 4)
    focus_pos = [];
end

% compute the number of source points and subsequent scaling factor (note,
% this assumes dx = dy = dz)
area = pi .* (diameter / 2).^2;
N_ongrid  = area ./ (kgrid.dx.^2);
    
% compute the number of off-grid sources using the upsampling rate
N_offgrid = ceil(N_ongrid .* upsampling_rate);

% compute scaling factor
scale = N_ongrid ./ N_offgrid;

% compute a set of uniformly spaced Cartesian points covering the disc
if isempty(focus_pos)
    disc_points = makeCartDisc(disc_pos, diameter/2, [], N_offgrid, false, use_spiral);
else
    disc_points = makeCartDisc(disc_pos, diameter/2, focus_pos, N_offgrid, false, use_spiral);
end

% recompute number of points if concentric sampling was used
if ~use_spiral
    N_offgrid = size(disc_points, 2);
    scale = N_ongrid ./ N_offgrid;
end

% remove disc points which are outside the simulation
disc_points = trimCartPoints(kgrid, disc_points);

% generate a source mask using BLIs centered on the Cartesian disc points
mask = offGridPoints(kgrid, disc_points, scale, ...
    'BLITolerance', bli_tolerance, 'BLIType', bli_type);   

% plot source points if required
if plot_disc
    figure;
    switch kgrid.dim
        case 2
            
            % plot source points
            plot(disc_points(2, :), disc_points(1, :), 'k.');
            axis(0.5 .* [kgrid.y_size .* [-1, 1], kgrid.x_size .* [-1, 1]]);
            axis equal;
            xlabel('y-position [m]');
            ylabel('x-position [m]');
            
        case 3
            
            % plot source points
            plot3(disc_points(2, :), disc_points(1, :), disc_points(3, :), 'k.')
            
            % plot beam axis
            dp = focus_pos - disc_pos;
            hold on;
            quiver3(disc_pos(2), disc_pos(1), disc_pos(3), ...
                dp(2), dp(1), dp(3), 0, 'k');

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
    title(sprintf('%d source points', size(disc_points, 2)))
end