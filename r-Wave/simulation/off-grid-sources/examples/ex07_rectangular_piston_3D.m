% Script to plot the acoustic response from a rectangle using kWaveArray
% and the acoustic field propagator.
%
% author: Bradley Treeby
% date: 14th November 2018
% last update: 5th December 2019

clearvars;

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% select which k-Wave code to run
%   1: MATLAB CPU code
%   2: C++ code
model       = 1;

% define element properties
f           = 1e6;          % [Hz]
Lx          = 1e-3;         % [m]
Ly          = 8e-3;         % [m]
theta       = [0, 0, 0];    % [deg]

% ppw for calculation
ppw         = 3;

% domain size for calculation
x_size      = 5e-3;         % [m]
y_size      = 20e-3;        % [m]
z_size      = 80e-3;       % [m]

% medium properties
c0          = 1500;         % [m/s]

% offset for off-grid source
z_offset    = 50;

% =========================================================================
% RUN CALCULATION
% =========================================================================

% calculate the grid spacing based on the PPW and F0
dx = c0 / (ppw * f);

% compute the size of the grid
Nx = roundEven(x_size / dx);
Ny = roundEven(y_size / dx);
Nz = roundEven(z_size / dx);

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% create empty kWaveArray
karray = kWaveArray;

% add rectangular element
position = [kgrid.x_vec(Nx/2), kgrid.y_vec(Ny/2), kgrid.z_vec(z_offset)];
karray.addRectElement(position, Lx, Ly, theta)

% define source
amp_in = karray.getArrayGridWeights(kgrid);
phase_in = 0;

% calculate the field using the acoustic field propagator
switch model
    case 1
        [amp_out, phase_out] = acousticFieldPropagator(amp_in, phase_in, dx, f, c0);
    case 2
        [amp_out, phase_out] = acousticFieldPropagatorC(amp_in, phase_in, dx, f, c0);
end

% trim off the source off-set
amp_out = amp_out(:, :, z_offset:end);

% =========================================================================
% PLOT
% =========================================================================

% get plot axis
x_vec = 1e3 * kgrid.x_vec;
y_vec = 1e3 * kgrid.y_vec;
z_vec = 1e3 * kgrid.z_vec(z_offset:end);
z_vec = z_vec - z_vec(1);

% extract planes
xz = squeeze(amp_out(:, round(end/2), :));
yz = squeeze(amp_out(round(end/2), :, :));

% normalise
xz_norm = xz ./ max(xz(:));
yz_norm = yz ./ max(yz(:));

% display central planes
figure;

subplot(2, 1, 1);
imagesc(z_vec, x_vec, xz);
title('x-z plane');
axis image;
colorbar

subplot(2, 1, 2);
imagesc(z_vec, y_vec, yz);
title('y-z plane');
axis image;
colormap(parula(512)); 
colorbar

% display central planes
figure;
plot_scale = [-30, 0];

subplot(2, 1, 1);
imagesc(z_vec, x_vec, 20*log10(xz_norm), plot_scale);
title('x-z plane');
axis image;
colorbar

subplot(2, 1, 2);
imagesc(z_vec, y_vec, 20*log10(yz_norm), plot_scale);
title('y-z plane');
axis image;
colormap(parula(512)); 
colorbar

% transmit receive sensitivity
sens = yz .* flip(yz, 2);
sens = sens ./ max(sens(:));

% plot the transmit receive sensitivity
figure;
imagesc(z_vec, y_vec, sens);
axis image;
colorbar;

% extract the profile in the middle
fwhm(sens(:, round(end/2)), y_vec, true);

% take log
sens_log = 20*log10(sens);

% get peak value from the middle
mid_val = max(sens_log(:, round(end/2)));

% plot the transmit receive sensitivity
figure;
imagesc(z_vec, y_vec, sens_log, [mid_val - 6, mid_val]);
axis image;
colorbar;