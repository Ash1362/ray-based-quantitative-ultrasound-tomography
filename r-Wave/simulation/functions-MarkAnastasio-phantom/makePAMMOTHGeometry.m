function setting = makePAMMOTHGeometry(radius, dx, radius_to_height, z_pos_height, sensor_dist, plot_setting)
%MAKEPAMMOTHGEOMETRY Create geometrical properties of the bowl-shaped 
% PAMMOTH scanner 
%
% DESCRIPTION:
%       makePAMMOTHGeometry returns a struct whose fields provide the
%       geometrical description of the bowl-shaped PAMMOTH scanner. The
%       bowl is centered at (0,0,0) and described by z <= 0. The input 
%       z_pos_height can be used to extend volume above z > 0. The spatial
%       grids are chosen such that they contain (0,0,0) but have an even
%       number of grid points.
%       Note that the corrdinates of the kgrid structure returned are 
%       shifted in z-direction such as the original coordinates used are 
%       not centered around z
%
% USAGE:
%       setting = makePAMMOTHGeometry(radius, dx)
%       setting = makePAMMOTHGeometry(radius, dx, z_pos_height)
%       setting = makePAMMOTHGeometry(radius, dx, z_pos_height, sensor_dist)
%       setting = makePAMMOTHGeometry(radius, dx, z_pos_height, sensor_dist, plot_setting)
%
% INPUTS:
%       radius - radius of the bowl [m]
%       dx     - spatial resolution of the model [m]
%       radius_to_height - the radius of the detection surface (or ring)
%                         to height (z size) of the grid
%       z_pos_height - heigth of volume included above z = 0 (default = 0)
%       sensor_dist  - additional volume included to guarantee a certain 
%                      distance of the sensors from the volume boundary
%                      (default = 0)
%       plot_setting - Boolean controlling whether the setting is
%                      plotted (default = false)
%
% OUTPUTS:
%       setting - struct describing the geometry with the fields
%
%       This has the following properties:
%
%       setting.Nx        - voxel in x/x/z direction
%                          (Nx, Nz accordinly)
%       setting.Nxyz      - [Nx, Ny, Nz] as an array;
%       setting.N         - total number of voxels Nx*Ny*Nz
%       setting.x_vec     - x grid locations (y_vec, z_vec accordingly)
%       setting.kgrid     - kgrid structure for acoustic simulation with 
%                           kWave
%       setting.kgrid0Ind - sub indices of center (0,0,0) in kgrid in the
%                           grid described by x_vec, y_vec, z_vec
%       setting.mask      - uint8 volume with lables: 
%                               0: below bowl surface z direction
%                               1: bowl interior
%                               2: bowl surface
%                               3: z > 0 volume
%       setting.radius    - radius or bowl
%
% ABOUT:
%       author          - Felix Lucka
%       date            - 18.02.2017
%       last update     - 10.07.2017

% check user defined value for z_pos_height, otherwise assign default value
if nargin < 3 || isempty(z_pos_height)
    z_pos_height = 0;
end

% check user defined value for sensor_dist, otherwise assign default value
if nargin < 4 || isempty(sensor_dist)
    sensor_dist = 0;
end

% check user defined value for plot_setting, otherwise assign default value
if nargin < 5 || isempty(plot_setting)
    plot_setting = false;
end

% number of voxels in +/- x/y direction
NxyOneDir = ceil((radius + sensor_dist) / dx);

% construct x and y grid locations such that 0 is in them and that they
% have an even length. For this, make them too large in negative direction
x_vec = dx * (-(NxyOneDir + 1):NxyOneDir);
Nx    = length(x_vec);
y_vec = x_vec;
Ny    = Nx;

% construct z grid, again such that 0 is in it and that it has even length
NzNeg   = ceil((radius/radius_to_height + sensor_dist) / dx);
NzPos   = ceil(max(z_pos_height, sensor_dist) / dx);
NzPos   = NzPos + mod(NzPos+NzNeg,2);
z_vec   = dx * [(-(NzNeg+1):1:0),(1:NzPos)];
Nz      = length(z_vec);

% construct kgrid
kgrid = kWaveGrid(Nx,dx,Ny,dx,Nz,dx);

% determine the index of (0,0,0) in the kgrid in the normal grid
kgrid0Ind    = zeros(3,1);
kgrid0Ind(1) = x_vec(kgrid.x_vec == 0);
kgrid0Ind(2) = y_vec(kgrid.y_vec == 0);
kgrid0Ind(3) = z_vec(kgrid.z_vec == 0);

% make extended spatials grids to construct segementation
x_vec_ext = [x_vec(1) - dx, x_vec, x_vec(end) + dx];
z_vec_ext = [z_vec(1) - dx, z_vec];
[X_ext, Y_ext, Z_ext] = ndgrid(x_vec_ext, x_vec_ext, z_vec_ext);

% define bow interior and boundary
bowl_mask_ext = (sqrt(X_ext.^2 + Y_ext.^2 + Z_ext.^2) <= radius & Z_ext <= 0);
bowl_surf_mask_ext   = extractBoundary(not(bowl_mask_ext));

% build setting mask
mask = zeros([Nx,Ny,Nz],'uint8');
mask(:, :, z_vec == 0) = 2;
mask(bowl_mask_ext(2:end-1, 2:end-1, 2:end))      = 1;
mask(bowl_surf_mask_ext(2:end-1, 2:end-1, 2:end)) = 2;
mask(:, :, z_vec > 0) = 3;



% gather results in struct
setting = [];
setting.dx         = dx;
setting.Nx         = Nx;
setting.Ny         = Ny;
setting.Nz         = Nz;
setting.Nxyz       = [Nx, Ny, Nz];
setting.N          = prod([Nx, Ny, Nz]);
setting.x_vec      = x_vec(:); 
setting.y_vec      = y_vec(:);
setting.z_vec      = z_vec(:);
setting.kgrid      = kgrid;
setting.kgrid0Ind  = kgrid0Ind;
setting.mask       = mask;
setting.radius     = radius;

setting = orderfields(setting);

% plot results
if plot_setting
    figure();
    imagesc(squeeze(double(setting.mask(:, Ny/2, :))));
    box(gca, 'on');
    set(gca, 'DataAspectRatio', [1 1 1], 'Layer', 'top');
end
