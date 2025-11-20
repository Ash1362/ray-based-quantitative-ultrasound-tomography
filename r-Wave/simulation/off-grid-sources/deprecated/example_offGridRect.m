% http://www.medata.cz/_docs/en_ultrasonix_transducers.pdf
% C5-2/60 Convex Transducer

% grid properties
PPW = 4;
C0 = 1500; % [m/s]
BYPASS_POINTS = false;

% element properties
ELEM_PITCH = 0.47e-3; % [m]
ELEM_LENGTH = 14e-3; % [m]
NUM_ELEM = 128;

% transducer properties
FOV   = deg2rad(56); % [rad]
ROC   = 61e-3; %[m]
SOURCE_FREQ = 2e6;%5e6; %[Hz]

% derived trandsucer properties
% assume elem_kerf = 0.5 * elem_width;
elem_width = ELEM_PITCH ./ 1.5;
elem_kerf = 0.5*elem_width;
diameter = 2*ROC*sin(FOV/2);

% construct a computational grid
lambda = C0 ./ SOURCE_FREQ;
z_size = 1.5 * diameter;
Nz = ceil(PPW .* z_size ./ lambda);
dz = z_size ./ Nz;
dx = dz;
dy = dz;
Nx = ceil(1.5 * ELEM_LENGTH ./ dx);
Ny = ceil(1.5 * diameter ./ dy);
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% generate an arc of points for the element centres
back_pos = [0, 0, -diameter/2];
rect_centres = zeros(3, NUM_ELEM);
rect_centres(2:3, :) = makeCartArc(back_pos(2:3), ROC, diameter, back_pos(2:3) + [0, -1], NUM_ELEM);

% compute source mask or source points
if BYPASS_POINTS
    amp_in = zeros(Nx, Ny, Nz);
else
    points = [];
    upsampling_rate = 4;
    area = ELEM_LENGTH .* elem_width;
    N_ongrid  = area ./ (dx.^2);
    N_offgrid = ceil(N_ongrid .* upsampling_rate);
end
for i = 1:NUM_ELEM
    i
    rect_pos = rect_centres(:, i).';
    theta = [atan2(rect_pos(2), rect_pos(3)), 0, 0];
    if BYPASS_POINTS
        rect = offGridRect(kgrid, rect_pos, ELEM_LENGTH, elem_width, theta);
        amp_in = amp_in + rect;
    else
        points = [points, makeCartRect(rect_pos, ELEM_LENGTH, elem_width, N_offgrid, theta)];
    end
end

% plot a projection of the source mask or the source points
figure
subplot(2, 2, 1)
if BYPASS_POINTS
    
    imagesc(squeeze(max(abs(amp_in), [], 1)))
    
else
    
    plot3(points(2, :), points(1, :), points(3, :), 'k.')
    
    % adjust axes
    view(45, 30);
    axis equal;
    box on;
    grid on;
    axis(0.5 .* [kgrid.y_size .* [-1, 1], kgrid.x_size .* [-1, 1], kgrid.z_size .* [-1, 1]]);
    xlabel('y-position [m]');
    ylabel('x-position [m]');
    zlabel('z-position [m]');
    
    % generate a source mask
    N_offgrid = size(points, 2);
    scale = N_ongrid ./ N_offgrid;
    amp_in = offGridPoints(kgrid, points, scale, 'WaitBar', true);
    
end

% compute the acoustic field
[pressure, amp_out, ~] = acousticFieldPropagatorC(amp_in, 0, dx, SOURCE_FREQ, C0);

% plot the results
[~, scale, prefix] = scaleSI(max([kgrid.x_vec; kgrid.y_vec; kgrid.z_vec]));

% plot maximum projection of field amplitude
subplot(2,2,2)
p_plot = squeeze(max(abs(amp_out), [], 1));
imagesc(kgrid.z_vec * scale, kgrid.y_vec * scale, p_plot);
axis image

% plot central plane of field amplitude
subplot(2,2,3)
p_plot = squeeze(amp_out(round(Nx/2), :, :));
imagesc(kgrid.z_vec * scale, kgrid.y_vec * scale, p_plot);
axis image

% plot central plane of pressure wave
subplot(2,2,4)
p_plot = squeeze(real(pressure(round(Nx/2), :, :)));
imagesc(kgrid.z_vec * scale, kgrid.y_vec * scale, p_plot);
axis image
