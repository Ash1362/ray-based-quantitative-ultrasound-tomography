% Script to demonstrate using kWaveArray to create an array with ten
% arc-shaped transducer elements. These elements are then used to detect a
% photoacoustic signal (initial pressure distribution).
%
% author: Bradley Treeby
% date: 4th September 2018
% last update: 6th July 2019

clearvars;

% =========================================================================
% DEFINE KWAVEARRAY
% =========================================================================

% create empty array
karray = kWaveArray;

% set the properties for the arc shaped elements and the ring geometry in
% which they're placed
radius_of_curv = 100e-3;
diameter       = 8e-3;
ring_radius    = 50e-3;
num_elements   = 20;

% orient all elements towards the centre of the grid
focus_pos = [0, 0];

% generate the centre position for each element in Cartesian space using
% makeCartCircle (these positions could also be defined manually, etc)
elem_pos = makeCartCircle(ring_radius, num_elements, [0, 0]);

% add elements to the array
for ind = 1:num_elements
    karray.addArcElement(elem_pos(:, ind), radius_of_curv, diameter, focus_pos);
end

% =========================================================================
% DEFINE GRID PROPERTIES
% =========================================================================

% grid properties
Nx = 256;
dx = 0.5e-3;
Ny = 256;
dy = 0.5e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% medium properties
medium.sound_speed = 1500;

% time array
kgrid.makeTime(medium.sound_speed);

% =========================================================================
% CONVENTIONAL SIMULATION
% =========================================================================

% set source as a square (directional) and a disc
source.p0 = makeDisc(Nx, Ny, Nx/4 + 20, Ny/4, 3);
source.p0(100:120, 50:200) = 1;

% assign Cartesian points
sensor.mask = elem_pos;

% run the k-Wave simulation using point detectors in the normal way
sensor_data_point = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% =========================================================================
% KWAVEARRAY SIMULATION
% =========================================================================

% assign binary mask from karray to the sensor mask
sensor.mask = karray.getArrayBinaryMask(kgrid);

% run k-Wave simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% combine data to give one trace per physical array element
combined_sensor_data = karray.combineSensorData(kgrid, sensor_data);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot binary source mask (note it is non-local)
figure;
imagesc(kgrid.y_vec, kgrid.x_vec, double(sensor.mask | source.p0));
axis image;
colormap(flipud(gray));

% overlay the physical source positions
hold on;
karray.plotArray(false);

% plot recorded sensor data
figure;
subplot(2, 1, 1)
imagesc(kgrid.t_array * 1e6, 1:num_elements, sensor_data_point);
xlabel('Time [\mus]');
ylabel('Detector Number');
title('Cartesian point detectors');
colorbar;

subplot(2, 1, 2);
imagesc(kgrid.t_array * 1e6, 1:num_elements, combined_sensor_data);
xlabel('Time [\mus]');
ylabel('Detector Number');
title('Arc detectors');
colorbar;