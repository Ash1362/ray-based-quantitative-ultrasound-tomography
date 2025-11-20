% Example for quantifying accuracy of ray tracing algorithms using a
% "Maxell-fish-eye lens" phantom. This example script reproduces 
% the numerical results in the reference [3].
%...........................................................................
% By running this example script, the accuracy of ray tracing algorithms are
% quantified by implementing them on a refractive index field for which
% rays' trajectoris are known, and are arcs of analytic circular or
% spherical paths for 2D and 3D cases, respectively.  
% The associated phantom is called "Maxwell’s Fish-eye lens", which is a
% refractive index field defined by $n(x) = n_0/(1+(x/a)^2)$, where $a$ is a
% scaling scalar, and was chosen 1 here.

%%=========================================================================
% CRITERIA (MEASURES)
%==========================================================================
% The accuracy of four ray tracing algorithms developed in this toolbox are
% quantified in terms of two criteria:
% "Criterion 1, Mean Acoustic Length Deviation:" Mean deviation of traced
% rays from the expected analytic acoustic lengths.

% "Criterion 2, Mean Radius Deviation:" Mean deviation of trajectories of the
% traced rays from the expected analytic circular or spherical paths for
% the 2D and 3D cases, respectively.

% For details, the users should read the "Numerical Results" section in the
% Reference [3].

%%=========================================================================
% RAY TRACING
%==========================================================================
% This script uses the ray tracing algorithms:
% 1- 'Dual-update' (section 3.3.1 in [1])
% 2- 'Mixed-step' (section 3.3.2. i [1])
% 3- 'Charactersitics' (section 5.1, Eq. (54) in [2])
% 4 -'Rung-kutta-2nd' (section 5.1, Eq (54) using Algorithm 2 in [2])
% ([3] and [4] were added after publishing the paper [1].

%%=========================================================================
% RAY-TO-GRID/GRID-TO-RAY INTERPOLATION
%==========================================================================
% The user can choose between the two interpolation approaches.
% 1- 'Bilinear' (section 3.4. in [1])
% 2- 'Bspline' (section 5.3. in [2])

%%%=========================================================================
% REFERENCES
%==========================================================================
% If you find the toolbox useful for your research, please consider citing these papers:...
% 1- A Javaherian, F. Lucka and B. Cox, ❝Refraction-corrected ray-based inversion
% for three-dimensional ultrasound tomography of the breast❞, Inverse Problems 36 125010.
% 2 - A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering
% for biomedical ultrasound tomography❞, Inverse Problems vol. 37, no.11, 115003, 2021.
% 3- A. Javaherian, ❝An open-source MATLAB package for quantitative ultrasound tomography via
% ray-Born inversion with in vitro and in vivo validaton❞, 2025.

% author: Ashkan Javaherian
% date:            - 21.08.2019
% last update:     - 07.09.2023
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian

% run the startup script for defining the paths
clear all
close all
clc

% run the startup
startup_simulation_ust;

%%=========================================================================
% THE INPUTS WHICH CAN BE CHOSEN BY THE USER
%==========================================================================
% get the number of dimensions of the grid, which can be set 2 or 3
dim = 2;

% choose the criterion for testing ray tracing algorithms. this criterion
% can be set '1' or '2'
% criterion '1' for Mean Acoustic Length Deviation
% criterion '2' for Mean Radius Deviation
criterion = '2';

% Choose the approach for interpolation from rays' smapled points onto the grid
% and vice versa
% This can be set 'Bilinear', or 'Bspline'
interp_method = 'Bspline';

% Get the ray tracing algorithm for displaying
% This can be set 'Mixed-step', 'Dual-update', 'Characteristics', 'Runge-Kutta-2nd'
raytracing_approach_display = 'Runge-Kutta-2nd';

% get a sampling rate for visualization of rays' paths (for visulazation purpose)
ray_sampling_rate = 4;

% get the Boolean controlling whether the figures are stored or not
save_image = true;

%==========================================================================

% display the number of dimensions
disp(['The number of dimensions is:' num2str(dim)])

% display the number of criterion for measuring the accuracy of ray tracing
disp(['The criterion being tested is:' criterion])

% display the approach for interpolation
disp(['The approach for ray-to-grid (or grid-to-ray) interpolation is:'...
    interp_method])




% get the grid-to-ray spacing power, for which the ray trajectories
% and error plots will be displayed. this must be chosen a scalar value
% from the vector gridtoray_spacing_power = 4.5:-0.5:-3;
% The greeter value gives smaller spacing for the sampled points along the
% ray.
gridtoray_spacing_power_display = 0;


if save_image

    % get the directory for storing the figures
    images_directory = 'results/simulation/fish_eye/';

    % make the directory for storing the images, if it does not exist
    makeDirectory(images_directory);
end




%%=========================================================================
% THE INPUTS WHICH WHICH MUST BE KEPT FIXED
%==========================================================================

% get the cell array of ray trcaing approaches
raytracing_approaches = {'Mixed-step', 'Dual-update', 'Characteristics',...
    'Runge-Kutta-2nd'};

% get the vector of the grid-to-ray spacing powers
gridtoray_spacing_power = 4.5:-0.5:-3;

% get the struct for the computational grid
test_grid.dim = dim;

% get the number of rays used for measuring the accuracy of ray tracing
switch dim
    case 2
        num_rays = 101;
    case 3
        num_rays = 21;
end


% get the maximum of the Cartesian coordinates
switch criterion
    case '1'
        switch test_grid.dim
            case 2
                max_coord = 2.0;
            case 3
                max_coord = 4.0;
        end
    case '2'
        max_coord = 4.0;
end

%  get the background refractive index
n_0 = 1;

% get the scalar $a$
a = 1;

% get the grid spacing [m] equal to 1 degrees on a cricle with unit radius
test_grid.dx = 2*pi*a/360 ;

% get the vector of the x Cartesian position of the grid points
test_grid.x_vec = (-max_coord: test_grid.dx : +max_coord)';

% get the vector of the y Cartesian position of the grid points
test_grid.y_vec = (-max_coord: test_grid.dx : +max_coord)';

% get the x-y Cartesian position of the first grid point
pos_grid_first = [test_grid.x_vec(1);  test_grid.y_vec(1)];

% get the x-y Cartesian position of the end grid point
pos_grid_end = [test_grid.x_vec(end); test_grid.y_vec(end)];

% get the number of the grid points along the x coordinate
test_grid.Nx = length(test_grid.x_vec);

% get the number of the grid points along the y coordinate
test_grid.Ny = length(test_grid.y_vec);

% get the number of grid points along the x-y coordinates
test_grid.size = [test_grid.Nx, test_grid.Ny];


switch test_grid.dim

    case 2


        % set the vector of z cartesian position of the grid points
        % empty
        test_grid.z_vec = [];

        % get the matrix of the x Cartesian position of the grid points
        test_grid.x  = repmat(test_grid.x_vec, [1 test_grid.Ny]);

        % get the matrix of the y Cartesian position of the grid points
        test_grid.y = repmat(test_grid.y_vec', [test_grid.Nx  1]);

        % get the map of distance of the grid points to the origin
        test_grid.r = sqrt ((test_grid.x).^2 + (test_grid.y).^2);

    case 3


        % get the vector of the Cartesian position of the grid points
        test_grid.z_vec = (-max_coord: test_grid.dx: +max_coord)';

        % get the z Cartesian position of the first grid point
        pos_grid_first = [pos_grid_first; test_grid.z_vec(1)];

        % get the z Cartesian position of the first grid point
        pos_grid_end = [pos_grid_end; test_grid.z_vec(end)];

        % get the length of the vector of the Cartesian position of the grid points
        test_grid.Nz = length(test_grid.z_vec);

        % get the number of the grid points along the z coordinates
        test_grid.size = [test_grid.size, test_grid.Nz];

        % get the matrix of the x Cartesian position of the grid points
        test_grid.x = repmat(test_grid.x_vec, [1, test_grid.Ny, test_grid.Nz]);

        % get the matrix of the y Cartesian position of the grid points
        test_grid.y = repmat(test_grid.y_vec.', [test_grid.Nx, 1, test_grid.Nz]);

        % get the matrix of the z Cartesian position of the grid points
        test_grid.z = repmat(permute(test_grid.z_vec, [2 3 1]), [test_grid.Nx, test_grid.Ny, 1]);

        % get the map of distance of the grid points to the origin
        test_grid.r = sqrt ((test_grid.x).^2 + (test_grid.y).^2 +...
            (test_grid.z).^2 );
end

% get the refractive index map as a function of the scalar $a$
n = n_0./(1 + (test_grid.r/a).^2) ;

% compute the directional gradients of the refractive index along the Cartesian
% coordinates
switch dim
    case 2
        [n_y, n_x] = gradient(n, test_grid.dx, test_grid.dx);
        n_z = [];
    case 3
        [n_y, n_x, n_z] = ...
            gradient(n, test_grid.dx, test_grid.dx, test_grid.dx);
end


switch test_grid.dim


    case 2

        % get the initial point of the ray on the circle with radius $a$
        initial_position = [0,a]';

        switch criterion
            case '1'

                % For the 2D case, given the initial position $p$, only one ray
                % can be traced along a circle with a given centre.
                % Instead, we trace a ray bundle whose centre is origin (c=o).
                % So for an initial position $p$ and origin $o$, the initial
                % direction of the central ray is $p-o$, and the initial
                % direction of other rays will be computed by a rotation of
                % $p-o$ with rotation angles in the interval [-pi/3,pi/3].
                % The user is referred to the obtained figure for the trajectory
                % (ray-to-grid coefficients) of rays after running this example

                % set the centre of trajectory of the ray bundle
                % origin (c = o)
                centre_position = [0,0]';

                % get the rotation angles for defining the initial direction
                % of the rays around the geometrical vector p-o
                directions_rotation_angles = linspace(-pi/3, pi/3, num_rays);

                % get the geometrical vector of the initial
                % direction at angle 0. The initial direction of the
                % rays are obtained by rotation of this geometrical
                % vector.
                directions_unrotated = centre_position - initial_position;

                % get the analytic (exact) acoustic length along a ray initialised
                % from a point and reaching its mirror with res
                acoustic_length_true = 2 * a * atan(1);


            case '2'

                % For the 2D case, given the initial position $p$, only one ray
                % can be traced along a circle with a given centre.
                centre_position = [a,0]';

                % only one ray along normal to the geomterical vector
                % connecting the initial position of the rays to the
                % centre of trajectory of rays (the vector p-c) is traced.
                directions_rotation_angles = 0;

                % the initial direction tangent to a circle is obtained by
                % the pi/2 rotation of the radius
                directions_unrotated = cat(1, [cos(pi/2), -sin(pi/2)],...
                    [sin(pi/2), cos(pi/2)]) * (centre_position  - initial_position);

                % get the true radius for rotation of the ray tangent to an
                % associated circle
                radius_true = norm(centre_position - initial_position);

                % set the number of rays 1
                num_rays = 1;

                disp('For 2D case and criterion 2, only one ray is traced.');

        end



    case 3


        % get the initial point of the ray on the circle with radius $a$
        initial_position = [0,0,a]';

        % let the point $p$ be the initial position of the rays and
        % $c$ the centre of the trajectory of the rays initialised
        % at $p$. the rays are traced with initial directions forming a
        % Rodriguez rotation along a geometrical vector orthognal
        % to the geometrical vector p-c.
        centre_position = [a,a,0]';

        % get the rotation angles for defining the initial direction
        % of the rays around a vector orthogonal to p-o
        directions_rotation_angles = linspace(0, (num_rays-1)/num_rays * (2*pi),...
            num_rays);

        % get the geometrical vector of the initial
        % direction at angle 0. The initial direction of the
        % rays are obtained by rotation of this goemetrical
        % vector.
        % (Note that this is the derivation by assuming a=1.)
        directions_unrotated = -[1/sqrt(2),1/sqrt(2),sqrt(2)]';

        % get the axis of the Rodriguez rotation
        directions_axis = centre_position - initial_position;

        switch criterion
            case '1'

                % get the analytic (exact) acoustic length along a ray initialised
                % from a point, traveling along a sphere and reaching the
                % initial point
                acoustic_length_true = 4 * a * atan(1);

            case '2'

                % get the true radius for rotation of the ray tangent to an
                % associated circle
                radius_true = norm(centre_position - initial_position);
        end

end

% get the number of the chosen ray spacings
num_ray_spacing = length(gridtoray_spacing_power);

switch criterion
    case '1'

        % get the cell array for the error in acoustic length using
        % different ray tracing (or interpolation) approaches with
        % respect to an analytically computed acoustic length
        acoustic_length_error = cell(length(raytracing_approaches), 1);
    case '2'

        % get the cell array for the error in rays' trajectory
        % using different ray tracing (or interpolation) approaches with
        % respect to an analytically computed acoustic length
        radius_error = cell(length(raytracing_approaches), 1);

end



if dim == 2 && strcmp(criterion,'2')

    hm = figure();
    % reduce the maximum coordinates for improving visualization
    max_coord_visual = 1;
    % set the binaries for the reduced coordinates
    x_binary = test_grid.x_vec>=- max_coord_visual  &  test_grid.x_vec<= max_coord_visual;
    y_binary = test_grid.y_vec>=- max_coord_visual  &  test_grid.y_vec<= max_coord_visual;
    % show the refractive index image
    imagesc(test_grid.x_vec(x_binary), test_grid.y_vec(y_binary),...
        n(x_binary, y_binary)); axis image;
    % add axis lables
    xticks(max_coord_visual * (-1:0.5:+1));
    yticks(max_coord_visual * (-1:0.5:+1));
    % add colorbars
    cl = colorbar;
    cl.Ticks =  0.1:0.1:1;  % set tick positions
    clim([0.01 1]);
    ylabel(cl, 'n [a.u.]', 'FontSize', 14);
    set(gca, 'FontSize', 14);

    % save the image of the Maxwell Fish-eye phantom
    saveas(hm, [images_directory, 'Dim' num2str(test_grid.dim)  '_Maxwell' '.fig']);
    saveas(hm, [images_directory, 'Dim' num2str(test_grid.dim)  '_Maxwell' '.tiff']);
    saveas(hm, [images_directory, 'Dim' num2str(test_grid.dim)  '_Maxwell' '.png']);
    saveas(hm, [images_directory, 'Dim' num2str(test_grid.dim)  '_Maxwell' '.eps'], 'epsc');

end





% get an empty figure for plotting criteria
hg = figure;

for approach_index = 1: length(raytracing_approaches)


    % allocate a variable for parameters reuired for ray tracing
    ray_gradients_interpcoeffs = [];

    switch interp_method


        case 'Bilinear'

            % get the directional gradients of the refrcative index
            ray_gradients_interpcoeffs.refractive_gradient_x = n_x;
            ray_gradients_interpcoeffs.refractive_gradient_y = n_y;
            ray_gradients_interpcoeffs.refractive_gradient_z = n_z;

        case 'Bspline'

            indices_vec = [-1; 0; 1; 2];
            switch dim
                case 2

                    % get the indices of the grid points included in the
                    % Bspline interpolation
                    ray_gradients_interpcoeffs.raytogrid_indices_x = vectorise(repmat(indices_vec, [1, 4]));
                    ray_gradients_interpcoeffs.raytogrid_indices_y = vectorise(repmat(indices_vec.', [4, 1]));
                    ray_gradients_interpcoeffs.raytogrid_indices_z = [];
                case 3

                    % get the indices of the grid points included in the
                    % Bspline interpolation
                    ray_gradients_interpcoeffs.raytogrid_indices_x = vectorise(repmat(indices_vec, [1, 4, 4]));
                    ray_gradients_interpcoeffs.raytogrid_indices_y = vectorise(repmat(indices_vec.', [4, 1, 4]));
                    ray_gradients_interpcoeffs.raytogrid_indices_z = vectorise(repmat(permute(indices_vec, [2 3 1]),...
                        [4, 4, 1]));
            end

            % get the coefficient matrix for interpolation of the refractive
            % index
            ray_gradients_interpcoeffs.raytogrid_coeff_matrix = 1/6 * [-1, 3,-3, 1;...
                3,-6, 0, 4;...
                -3, 3, 3, 1;...
                1, 0, 0, 0];

            % get the coefficient matrix for interpolation of the fisrt gradient
            % of the refractive index
            ray_gradients_interpcoeffs.raytogrid_coeff_derivative_matrix = 1/(6 * test_grid.dx) * [-3, 6, -3;...
                9,-12, 0;...
                -9, 6, 3;...
                3,  0, 0];


    end


    switch criterion
        case '1'

            % allocate zero matrices for the acoustic length and its
            % error for the current approach
            acoustic_length_error{approach_index} = zeros(num_rays, num_ray_spacing);
            acoustic_length_rays = zeros(num_rays, num_ray_spacing);

        case '2'

            % allocate a zero matrix for the deviation of the radius of
            % the rays' trajectories for the current apprpach
            radius_error{approach_index} = zeros(num_rays, num_ray_spacing);
    end


    for ray_spacing_index = 1:num_ray_spacing

        % get the grid to ray spacing
        gridtoray_spacing = 2.^gridtoray_spacing_power(ray_spacing_index);

        % get the ray spacing
        ray_spacing = test_grid.dx./gridtoray_spacing;


        % get the the handle function for ray tracing
        switch approach_index

            case 1

                % choose the Mixed step algorithm as the approach for ray tracing.
                calc_acoustic_length = @(cartesian_initial_direction) calcRayMixedStep(...
                    n, ray_gradients_interpcoeffs, initial_position, cartesian_initial_direction,...
                    test_grid.x_vec, test_grid.y_vec, test_grid.z_vec, pos_grid_first, pos_grid_end,...
                    test_grid.dx, ray_spacing, test_grid.size, test_grid.dim, criterion,...
                    true(test_grid.size), true);

                % get the plot properties
                plot_prop = 'b-o';

            case 2

                % choose the Dual update algorithm as the approach for ray tracing.
                calc_acoustic_length = @(cartesian_initial_direction) calcRayDualUpdate(...
                    n, ray_gradients_interpcoeffs, initial_position, cartesian_initial_direction,...
                    test_grid.x_vec, test_grid.y_vec, test_grid.z_vec, pos_grid_first, pos_grid_end,...
                    test_grid.dx, ray_spacing, test_grid.size, test_grid.dim, criterion,...
                    true(test_grid.size), true);

                % get the plot properties
                plot_prop = 'r--s';

            case 3

                % choose the method of Characteristics as the approach for ray tracing.
                calc_acoustic_length = @(cartesian_initial_direction) calcRayCharacteristics(...
                    n, ray_gradients_interpcoeffs, initial_position, cartesian_initial_direction,...
                    test_grid.x_vec, test_grid.y_vec, test_grid.z_vec, pos_grid_first, pos_grid_end,...
                    test_grid.dx, ray_spacing, test_grid.size, test_grid.dim, criterion,...
                    true(test_grid.size), true);

                % get the plot properties
                plot_prop = 'g:^';
            case 4

                % choose the second order Runge Kutta secheme as the approach for ray tracing.
                calc_acoustic_length = @(cartesian_initial_direction) calcRayRungeKutta2nd(...
                    n, ray_gradients_interpcoeffs, initial_position, cartesian_initial_direction,...
                    test_grid.x_vec, test_grid.y_vec, test_grid.z_vec, pos_grid_first, pos_grid_end,...
                    test_grid.dx, ray_spacing, test_grid.size, test_grid.dim, criterion ,...
                    true(test_grid.size), true);

                % get the plot properties
                plot_prop = 'k-.d';
        end

        % allocate a cell array of length num_rays for the position of
        % points along the rays
        ray_points_position = cell(num_rays, 1);


        % if sctrcom(criterion, '1')

        % allocate a zero vector for the grid-to-ray coefficients
        % or points on the ray
        %rays_coeff = zeros(prod(test_grid.size), 1);

        % end

        % start the time
        ts = tic;
        for ray_index = 1:num_rays

            % get the rays' initial directions via rotation along the
            % chosen axis
            switch test_grid.dim
                case 2
                    [ray_initial_direction] = rotateRodriguez(directions_unrotated,...
                        directions_rotation_angles(ray_index));
                case 3
                    [ray_initial_direction] = rotateRodriguez(directions_unrotated,...
                        directions_rotation_angles(ray_index), directions_axis);
            end



            switch criterion
                case '1'

                    % get the ray-to-grid interpolation coefficients and the Cartesian position of the
                    % sampled points along the ray, and the cartesian end position
                    [end_position, coeff_vec, single_ray_position] = calc_acoustic_length(ray_initial_direction);

                    % get the number of points along the ray excluding the initial
                    % point, which is obviously on the circular trajectory
                    num_ray_points = size(single_ray_position, 2);

                    % get the acoustic length of the ray
                    acoustic_length_rays(ray_index, ray_spacing_index) = coeff_vec' * n(:);

                    % get the relative error in the acoustic length
                    acoustic_length_error{approach_index}(ray_index, ray_spacing_index) =...
                        abs(acoustic_length_rays(ray_index, ray_spacing_index)-acoustic_length_true)/...
                        acoustic_length_true * 100;

                    % add the ray-to-grid coefficients
                    % rays_coeff = rays_coeff + coeff_vec(:);


                case '2'


                    % get the Cartesian position of the sampled points along the
                    % ray and the Cartesian end position
                    [end_position, ~, single_ray_position] = calc_acoustic_length(ray_initial_direction);

                    % get the number of points along the ray excluding the initial
                    % point, which is obviously on the circular trajectory
                    num_ray_points = size(single_ray_position, 2);


                    % get the distance of the ray points to the centre of the ray trajectory [m]
                    radius_ray_points = vecnorm(single_ray_position - centre_position).';

                    % get the relative error in the position of the ray in
                    % percentage
                    radius_error{approach_index}(ray_index, ray_spacing_index) = mean(1/...
                        radius_true * abs(radius_ray_points(2:end) - radius_true)) * 100;


            end



            % the number of sampled points on the ray
            disp(['The number of the sampled points on the ray is:'...
                num2str(num_ray_points)]);

            % add the positions as rows to the matrix
            ray_points_position{ray_index} = single_ray_position;

        end


        if abs(gridtoray_spacing_power(ray_spacing_index)- gridtoray_spacing_power_display)< eps


            switch criterion

                case '1'

                    % plot the analytical acoustic lengths
                    if approach_index == 1
                        h_all(1) = plot(directions_rotation_angles, acoustic_length_true * ones(1, num_rays),'m',...
                            'MarkerSize', 3); hold on;
                    end

                    % plot the acoustic lenths of rays' points associated
                    % with the current ray tracing approach
                    h = plot(directions_rotation_angles, acoustic_length_rays(:, ray_spacing_index), plot_prop,'MarkerSize', 3); hold on;

                    % add information for the current plot to the all plots
                    h_all = [h_all, h];

                    if approach_index == 4

                        switch dim
                            case 2
                                xlabel('$\theta\ \ \mathrm{[rad]}$', 'Interpreter', 'latex');
                                ylabel('$L_{(\theta;p,\tilde{p})} \ \mathrm{[m]}$', 'Interpreter', 'latex');
                            case 3
                                xlabel('$\varphi\ \ \mathrm{[rad]}$', 'Interpreter', 'latex');
                                ylabel('$L_{(\varphi;p,p)} \ \mathrm{[m]}$', 'Interpreter', 'latex');
                        end

                        % Add the legend and hold off the figure
                        legend(['Analytical', raytracing_approaches], 'Location', 'northwest'); hold off;

                    end



                    %if  strcmp(raytracing_approaches{approach_index}, raytracing_approach_display)

                    % display the ray trajectories
                    %rays_coeff = reshape(rays_coeff, test_grid.size);

                    %switch test_grid.dim
                    %   case 2
                    %       rays_coeff((test_grid.x - initial_position(1)).^2+...
                    %           (test_grid.y - initial_position(2)).^2 < 4 * test_grid.dx^2 |...
                    %           (test_grid.x - end_position(1)).^2+...
                    %           (test_grid.y - end_position(2)).^2 < 4 * test_grid.dx^2) = 0;
                    %       figure; imagesc(rays_coeff); axis image;

                    %   case 3
                    %       rays_coeff((test_grid.x - initial_position(1)).^2+...
                    %           (test_grid.y - initial_position(2)).^2 +....
                    %           (test_grid.z - initial_position(3)).^2 < 4 * test_grid.dx^2) = 0;
                    %       figure; scrollView(rays_coeff, 1);

                    %end

                    %end



                case '2'

                    % plot the acoustic lengths
                    if approach_index == 1
                        h_all = plot(1:ray_sampling_rate:num_ray_points,...
                            radius_true * ones(size(1:ray_sampling_rate:num_ray_points)),...
                            'm', 'MarkerSize', 1); hold on;
                    end

                    % plot distances of the rays's points associated with the current
                    % ray tracing approach to the center
                    h = plot(1:ray_sampling_rate:num_ray_points,...
                        radius_ray_points(1:ray_sampling_rate:num_ray_points),...
                        plot_prop, 'MarkerSize', 1); hold on;

                    % add information for the current plot to the all plots
                    h_all = [h_all, h];

                    if approach_index == 4

                        xlabel('$m$', 'Interpreter', 'latex'); 

                        switch dim
                            case 2
                                
                                % add ylabel
                                ylabel('$\left\| \mathbf{x}_{(s_m;\theta=\pi/2;p,p^*)}-\mathbf{x}_c\right\| \ \mathrm{[m]}$', 'Interpreter', 'latex');
                                % Add the legend and hold off the figure
                                legend(['Analytical', raytracing_approaches], 'Location', 'northeast');
                            case 3
                                % add ylabel
                                ylabel('$\left\| \mathbf{x}_{(s_m;\varphi_{\mathrm{max}};p,p^*)}-\mathbf{x}_c\right\| \ \mathrm{[m]}$', 'Interpreter', 'latex');
                                % Add the legend and hold off the figure
                                legend(['Analytical', raytracing_approaches], 'Location', 'northwest');
                        end
                        hold off;

                    end


            end

            % Set font size for current axes
            set(gca, 'FontSize', 14);

        end



        % get the CPU time
        cpu_time = toc(ts);

        switch criterion
            case '1'


                % disp the mean deviation of the rays' accumulated acoustic length
                % from the analytical acoustic length
                disp(['The mean deviation of the  is:' num2str(mean(mean(acoustic_length_error{approach_index})), '%1.5f')...
                    '[m]']);

            case '2'
                % disp the mean of the deviation of the rays' points from
                % the centre of the rays' trajectories
                disp(['The mean deviation of the ray points from the centre of '...
                    'the ray trajectory is:' num2str(mean(mean(radius_error{approach_index})), '%1.5f')...
                    '[m]']);
        end


        if strcmp(raytracing_approaches{approach_index}, raytracing_approach_display) &&...
                abs(gridtoray_spacing_power(ray_spacing_index)- gridtoray_spacing_power_display)< eps

            % make an empty figure
            hs = figure;
            switch test_grid.dim

                case 2


                    % plot the percentage relative error in the acoustic length in the base-10
                    % logarithmic scale for different ray tracing algorithms
                    for ray_index = 1:num_rays

                        % get the number of sampled points for
                        % the current ray
                        num_ray_points = size(ray_points_position{ray_index}, 2);

                        % plot the trajectory of rays
                        scatter(ray_points_position{ray_index}(1,1:ray_sampling_rate:num_ray_points),...
                            ray_points_position{ray_index}(2,1:ray_sampling_rate:num_ray_points), 5, 'k');
                        hold on;
                    end

                    % add a pointer to the initial position
                    b = initial_position + [0.5, 0].';
                    dx_q = initial_position - b;
                    quiver(b(1), b(2), dx_q(1), dx_q(2),...
                        'Color', 'b', 'LineWidth', 1, 'MaxHeadSize', 2); hold on;
                    scatter(initial_position(1), initial_position(2), 30, 'b', 'filled');
                    % Red text label
                    text(b(1)+ 0.4, b(2), '$\mathbf{x}_p = [0,1]$', ...
                        'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
                        'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold'); hold on;

                    % add a pointer to the mirror position
                    if strcmp(criterion, '1')
                        b = -initial_position + [0.5, 0].';
                        dx_q = -initial_position - b;
                        quiver(b(1), b(2), dx_q(1), dx_q(2),...
                            'Color', 'b', 'LineWidth', 1, 'MaxHeadSize', 2); hold on;
                        scatter(-initial_position(1), -initial_position(2), 30, 'b', 'filled');
                        % Red text label
                        text(b(1) + 0.4, b(2), '$\mathbf{x}_{\tilde{p}} = [0,-1]$', ...
                            'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
                            'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold'); hold on;
                    end

                    % add a pointer to the center position
                    b = centre_position + [0.5, 0].';
                    dx_q = centre_position - b;
                    quiver(b(1), b(2), dx_q(1), dx_q(2), 'Color', 'r',...
                        'LineWidth', 1, 'MaxHeadSize',2); hold on;
                    scatter(centre_position(1), centre_position(2), 30, 'r', 'filled');
                    % Red text label
                    switch criterion
                        case '1'
                            text(b(1) + 0.4, b(2), '$\mathbf{x}_c = [0,0]$', ...
                                'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
                                'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold'); hold on;
                        case '2'
                            text(b(1) + 0.4, b(2), '$\mathbf{x}_c = [1,0]$', ...
                                'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
                                'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold');
                    end

                    axis image; xlabel('x [m]'); ylabel('y [m]');
                    hold off;

                case 3



                    for ray_index = 1:num_rays

                        % get the number of sampled points for
                        % the current ray
                        num_ray_points = size(ray_points_position{ray_index}, 2);

                        % plot the trajectory of rays
                        scatter3(ray_points_position{ray_index}(1,1:ray_sampling_rate:num_ray_points),...
                            ray_points_position{ray_index}(2,1:ray_sampling_rate:num_ray_points),...
                            ray_points_position{ray_index}(3,1:ray_sampling_rate:num_ray_points), 5, 'k');
                        hold on;
                    end

                    % add a pointer to the initial position
                    b = initial_position + [-1, -1, 1].';
                    dx_q = initial_position - b;
                    quiver3(b(1), b(2), b(3), dx_q(1), dx_q(2), dx_q(3),...
                        'Color', 'b', 'LineWidth', 1, 'MaxHeadSize',2); hold on;
                    scatter3(initial_position(1), initial_position(2), initial_position(3), 30, 'b', 'filled');
                    % Red text label
                    text(b(1), b(2), b(3) + 0.1, '$\mathbf{x}_p = [0,0,1]$', ...
                        'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
                        'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold'); hold on;

                    % add a pointer to the center position
                    b = centre_position + [0.5, 0.5, 0.5].';
                    dx_q = centre_position - b;
                    quiver3(b(1), b(2), b(3), dx_q(1), dx_q(2), dx_q(3),...
                        'Color', 'r','LineWidth', 1, 'MaxHeadSize',2); hold on;
                    scatter3(centre_position(1), centre_position(2), centre_position(3), 30, 'r', 'filled');
                    % Red text label
                    text(b(1), b(2), b(3) + 0.1, '$\mathbf{x}_c = [1,1,0]$', ...
                        'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
                        'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold'); hold on;
                    axis image; xlabel('x [m]'); ylabel('y [m]');zlabel('z [m]'); view(-15, 15);
                    hold off;

            end


            % Set font size for current axes
            set(gca, 'FontSize', 14);

            if save_image

                % save the plot in the resuts directory
                % Criterion 1: Acoustic length error (%)
                % Criterion 2: Distances error (%)
                saveas(hs, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Points' '.fig']);
                saveas(hs, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Points' '.tiff']);
                saveas(hs, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Points' '.png']);
                saveas(hs, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Points' '.eps'], 'epsc');

            end

        end

        % display the ray-to-grid spcing
        disp(['The ray-to-grid spacing is:' num2str(1/gridtoray_spacing, '%2.5e')])


    end

    disp(['The ray tracing using the approach '...
        raytracing_approaches{approach_index} ' was completed.'])

end


if save_image

    % save the plot in the resuts directory
    % Criterion 1: Acoustic length error (%) against ray-to-grid spacing for
    % all ray tracing approaches
    % Criterion 2: Distances error (%) against ray-to-grid spacing for
    % all ray tracing approaches
    saveas(hg, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Plot1' '.fig']);
    saveas(hg, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Plot1' '.tiff']);
    saveas(hg, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Plot1' '.png']);
    saveas(hg, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Plot1' '.eps'], 'epsc');

end

% make an empty figure
hp = figure;

switch criterion
    case '1'

        % plot the percentage relative error in the acoustic length in the base-10
        % logarithmic scale for different ray tracing algorithms
        for approach_index = 1:length(raytracing_approaches)

            % get the plot properties
            switch approach_index
                case 1
                    plot_prop = 'b-o';
                case 2
                    plot_prop = 'r--s';
                case 3
                    plot_prop = 'g:^';
                case 4
                    plot_prop = 'k-.d';
            end


            loglog(2.^-gridtoray_spacing_power,...
                mean(acoustic_length_error{approach_index}), plot_prop);
            hold on;
        end
        hold off;
        xlabel('\Delta s / \Delta x');
        ylabel('RE_{al} (%)');
        legend(raytracing_approaches, 'Location', 'northwest');

    case '2'

        % plot the percentage relative error in the deviation of the radius
        % of the trajectory in the base-10 logarithmic scale for
        % different ray tracing algorithms
        for approach_index = 1:length(raytracing_approaches)

            % get the plot properties
            switch approach_index
                case 1
                    plot_prop = 'b-o';
                case 2
                    plot_prop = 'r--s';
                case 3
                    plot_prop = 'g:^';
                case 4
                    plot_prop = 'k-.d';
            end

            loglog(2.^-gridtoray_spacing_power, mean(radius_error{approach_index}, 1),...
                plot_prop);
            hold on;
        end

        hold off;

        xlabel('\Delta s / \Delta x');
        ylabel('RE_{rd} (%)');
        legend(raytracing_approaches, 'Location', 'northwest');

end


% Set font size for current axes
set(gca, 'FontSize', 14);


if save_image

    % save the plot in the resuts directory
    % Criterion 1: Mean Acoustic length error (%) against ray-to-grid spacing for
    % all ray tracing approaches
    % Criterion 2: Mean Distances error (%) against ray-to-grid spacing for
    % all ray tracing approaches
    saveas(hp, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Plotrg' '.fig']);
    saveas(hp, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Plotrg' '.tiff']);
    saveas(hp, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Plotrg' '.png']);
    saveas(hp, [images_directory, 'Dim' num2str(test_grid.dim) '_Criterion' criterion '_Plotrg' '.eps'], 'epsc');

end


