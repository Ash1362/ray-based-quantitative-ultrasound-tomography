% Example for visualizing the results in the image-reconstruction part of the
% scenarios in the reference [1] (cf. section References), given the
% simulated and saved results after running the script:
% 'example_2D_image_tof_greens_circle.m'
% By setting script_task: 'reconstructed_image' or 'line_plot', the
% reconstructed images or line plots in the reference [1] are visulaized,
% respectively.
% Before running the currect m-file script, the example script:
% 'example_2D_image_tof_greens_circle.m'
% must be run for all the 27 possible cases:
%..........................................................................
% 1) The signal-to-noise ratio of the simulated ultrasound data
% This can be set 40, 30, or 25.
% 2) The optimisation approach for image reconstruction
% using the ray approaximation to heterogeneous Green's function.
% This can be set 'hessian' or 'backprojection'.
% 'hessian' itself can be run by setting the variable: do_hom_greens, which
% can be set true or false. By setting do_hom_greens = true, the Green's functions
% are approximated by the assumption of homogeneous water (only), and the
% updated wavenumber fields are incorporated only into the scattering
% potential maps. Therefore, setting do_hom_greens = true is equivalent to
% implementing a prototype image recosntruction using the Born
% approximation, which is used as a benchmark for evaluating the accuracy
% of the image reconstruction approaches using ray approximaton to the Green's
% function. By setting do_hom_greens = false, the updated wavenumber
% fields are incorporated into both the Green's functions and the scattering
% potential maps. Therefore, setting do_hom_greens = false is equivalent to
% an image reconstruction approach using a ray approximation to
% the Green's function.
% For the optimization approach 'backprojection', do_hom_greens is always
% set false.
%  For do_hom_greens = false,
% 'backprojection' is based on the Green's approach proposed in [1].
% 'hessian' is based on the Green's approach proposed in [2].
% 3) The type of the absorption coefficient map used for image reconstruction
% using the ray approaximation to heterogeneous Green's function. This can
% be set: 'true', 'homogeneous', or 'none'.
% By setting do_hom_greens = false, the chosen spatially-varying absorption
% coefficient map \(alpha_0\) and its associated exponent power \(y\),
% which is assumed spatially constant, are incorporated into both the Green's
% functions and scattering potentials.
% By setting do_hom_greens = false, the chosen spatially-varying absorption
% coefficient map \(alpha_0\) and its associated exponent power \(y\),
% which is assumed spatially constant, are incorporated only in the
% scattering potential. This is equivalent to a Born approximation using
% a scattering potential map in terms of a complex wavenumber.
% Please read section Numerical results in [1].
%%=========================================================================
% DIGITAL PHANTOM AND K-WAVE
%==========================================================================
% The Breast is simulated using a digital breast phantom developed by the group
% of Professor Mark Anastasio [3]. The phantom data must be downloaded via:'...
% https://anastasio.bioengineering.illinois.edu/downloadable-content/oa-breast-database/ [3].
% The synthetic UST pressure times series used in this study are simulated using the k-Wave
% toolbox. www.k-Wave.org (v. 1.3. or 1.4.) [4].

%
%%=========================================================================
% REFERENCES
%==========================================================================
% If you find the toolbox useful for your research, please consider citing these papers:...
% 1 - A. Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution
% quantitative ultrasound tomography❞, 2022, https://arxiv.org/abs/2211.00316.
% 2 - A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering
% for biomedical ultrasound tomography❞, Inverse Problems vol. 37, no.11, 115003, 2021.

% These data-base/toolbox were used in this project.
% 3- Y. Lou, W. Zhou, T. P. Matthews, C. M. Appleton and M. A. Anastasio, ❝Generation of anatomically realistic
% numerical phantoms for photoacoustic and ultrasonic breast imaging❞, J. Biomed. Opt., vol. 22, no. 4, pp. 041015,
% 2017.
% 4 - B. E. Treeby and B. T. Cox, ❝k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic
% wave fields❞, J. Biomed. Opt. vol. 15, no. 2, 021314, 2010.
%
% You may also think about citing this preprint:
% 5- A Javaherian, ❝Full-waveform Approximation of Finite-Sized Acoustic Apertures:
% Forward and Adjoint Wavefields❞, https://arxiv.org/abs/2212.04466.
%
%
%
% author: Ashkan Javaherian
% date:            - 17.02.2020
% last update:     - 07.11.2024
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian


% the three essential c's
clear all; close all; clc;

% run the startup script for defining the paths
startup_simulation_ust;

% get the requested task, which is done by running the Matlab script
% This can be set 'reconstructed_image' or 'line_plot'.
script_task = 'line_plot'; %'reconstructed_image'

% get the path for loading the results stored in the workspace
results_directory = 'results/simulation/2D/greens5/';

% make the directory for the results in the workspace, if it does not exist
makeDirectory(results_directory);

% get the path for saving figures
images_directory = 'results/simulation/2D/images5/';

% make the directory for saving figures, if it does not exist
makeDirectory(images_directory);

% get all noise levels
noise_levels = {'40','30','25'};

% get all asumptions for the absorption map
absorption_maps = {'true','homogeneous','none'};

% get all Greens-based optimisation approaches
greens_optimisation_approaches = {'hessian','hessian','backprojection'};

% get the color range for displaying
display_color_range = [1440, 1600];


switch script_task

    case 'reconstructed_image'



        % get the subfigures indices
        subfigures = 'abcdefghijk';



        for i = 1:length(noise_levels)
            for j = 1:length(absorption_maps)
                for k = 1:length(greens_optimisation_approaches)


                    if j==1 && k==1

                        % load the image reconstructed using the Green's approach for the current cases
                        load([results_directory, 'results_' greens_optimisation_approaches{k}, '_' absorption_maps{j},...
                            '_' num2str(k==1) '_' noise_levels{i} 'db'], 'img_tof');


                        % display the image reconstructed using the time-of-flight-based approach
                        ht = figure;imagesc(img_tof, display_color_range);
                        colormap(gray);
                        c = colorbar;
                        c.Label.String = 'c [m s^{-1}]'; c.YTick = display_color_range(1):20:display_color_range(end);
                        c.FontSize = 12; axis image; axis off;


                        % get the file name for the image reconstructed using the time-of-flight-based approach
                        file_name_tof = [num2str(i+4), subfigures(2)];


                        % save the image reconstructed using the time-of-flight-based approach
                        saveas(ht, [images_directory, 'Fig' file_name_tof '.fig']);
                        saveas(ht, [images_directory, 'Fig' file_name_tof '.png']);
                        saveas(ht, [images_directory, 'Fig' file_name_tof '.tiff']);
                        saveas(ht, [images_directory, 'Fig' file_name_tof '.eps'], 'epsc');



                    end

                    % load the image reconstructed using the Green's approach for the current cases
                    load([results_directory, 'results_' greens_optimisation_approaches{k},  '_' absorption_maps{j},...
                        '_' num2str(k==1) '_' noise_levels{i} 'db'], 'img_greens');

                    % display the image reconstructed using the Green's approach
                    hg = figure;imagesc(img_greens, display_color_range);
                    colormap(gray);
                    c = colorbar;
                    c.Label.String = 'c [m s^{-1}]'; c.YTick = display_color_range(1):20:display_color_range(end);
                    c.FontSize = 12; axis image; axis off;

                    % get the file name for the image reconstructed using the Green's approach
                    file_name_greens = [num2str(i+4), subfigures(3*j+k-1)];

                    % save the image reconstructed using the Green's approach
                    saveas(hg, [images_directory, 'Fig' file_name_greens '.fig']);
                    saveas(hg, [images_directory, 'Fig' file_name_greens '.png']);
                    saveas(hg, [images_directory, 'Fig' file_name_greens '.tiff']);
                    saveas(hg, [images_directory, 'Fig' file_name_greens '.eps'], 'epsc');

                end
            end
        end

    case 'line_plot'

        % get the subfigures indices
        subfigures = 'abcdefghijkl';

        % get the radius of the detection ring
        detec_radius = 0.0950;


        % load the sound speed phantom
        load('results/simulation/2D', 'simulation_prop');

        % load the image reconstructed using the Green's approach for the current cases
        load([results_directory, 'results_' greens_optimisation_approaches{1},  '_' absorption_maps{2},...
            '_1'  '_' noise_levels{1} 'db'], 'recon_grid');


        % get the x vector of the phantom grid
        x_vec_phantom = simulation_prop.x(:,1);

        % get the y vector of the phantom grid
        y_vec_phantom = simulation_prop.y(1,:);

        % get the x vector of the image reconstruction grid
        x_vec_image = recon_grid.x(:,1);

        % get the y vector of the image reconstruction grid
        y_vec_image = recon_grid.y(1,:);




        for plot_number = 1:4

            switch plot_number

                case 1

                    % get the index in the phantom grid for the x=0 line
                    [~, ix_phantom] = min(abs(x_vec_phantom));

                    % get the index in the image reconstruction grid for the x=0 line
                    [~, ix_image] = min(abs(x_vec_image));

                    % get the y coordinate of the sampled points on the data-simulation grid
                    coord_phantom = y_vec_phantom;

                    % get the y coordinate of the sampled points on the image-reconstruction grid
                    coord_image = y_vec_image;

                    % get a binary vector for the y coordinate of the data-simulation grid
                    binary_phantom = coord_phantom >= -0.1 & coord_phantom <= +0.1;

                    % get a binary vector for the y coordiante of the image-reconstruction grid
                    binary_image = coord_image >= -0.1 & coord_image <= +0.1;

                    % get the values of the sound speed phantom along the x=0 line
                    sound_speed_phantom = simulation_prop.sound_speed(ix_phantom, :);

                case 2

                    % get the index in the phantom grid for the y=0 line
                    [~, iy_phantom] = min(abs(y_vec_phantom));

                    % get the index in the image reconstruction grid for the y=0 line
                    [~, iy_image] = min(abs(y_vec_image));

                    % get the y coordinate of the sampled points on the data-simulation grid
                    coord_phantom = y_vec_phantom;

                    % get the y coordinate of the sampled points on the image-reconstruction grid
                    coord_image = y_vec_image;

                    % get a binary vector for the y coordinate of the data-simulation grid
                    binary_phantom = coord_phantom >= -0.1 & coord_phantom <= +0.1;

                    % get a binary vector for the y coordiante of the image-reconstruction grid
                    binary_image = coord_image >= -0.1 & coord_image <= +0.1;

                    % get the values of the sound speed phantom along the y=0 line
                    sound_speed_phantom = simulation_prop.sound_speed(:, iy_phantom);


                case {3,4}

                    %  get the radial coordinate of the sampled points on the first diagonal of the
                    %  grid for data simulation
                    coord_phantom = sign(x_vec_phantom) .* sqrt(x_vec_phantom.^2 + y_vec_phantom.'.^2);

                    %  get the radial coordinate of the sampled points on the first diagonal of the
                    %  grid for image reconstruction
                    coord_image = sign(x_vec_image) .* sqrt(x_vec_image.^2 + y_vec_image.'.^2);

                    % get a binary vector for the radial coordinate of the data-simulation grid
                    binary_phantom = coord_phantom >= -0.1 & coord_phantom <= +0.1;

                    % get a binary vector for the radial coordiante of the image-reconstruction grid
                    binary_image = coord_image >= -0.1 & coord_image <= +0.1;

                    if plot_number == 3

                        % get the values of the sound speed phantom along the first diagonal
                        sound_speed_phantom = diag(simulation_prop.sound_speed);

                    else

                        % get the values of the sound speed phantom along the second diagonal
                        sound_speed_phantom = diag(fliplr(simulation_prop.sound_speed));

                    end

            end


            for i = 1:length(noise_levels)


                hl = figure(); hl0 = plot(100*coord_phantom(binary_phantom),...
                    sound_speed_phantom(binary_phantom), 'k-');
                hl0.LineWidth = 1;
                hold on;


                for k = 1:3


                    if k==1

                        % For the current noise level, load the image reconstructed using
                        % the Hessian-based Green's approach
                        load([results_directory, 'results_' greens_optimisation_approaches{k},  '_' absorption_maps{2},...
                            '_' num2str(1) '_' noise_levels{i} 'db'], 'img_tof');

                        % get the line sound speed values of the image reconstructed
                        % using the time-of-flight-based approach
                        switch plot_number
                            case 1
                                % along the x=0 axis
                                img = img_tof(ix_image,:);
                            case 2
                                % along the y=0 axis
                                img = img_tof(:,iy_image);
                            case 3
                                % along the fisrt diagonal
                                img = diag(img_tof);
                            case 4
                                % along the second digonal
                                img = diag(fliplr(img_tof));
                        end

                        % plot the sound speed values
                        hl1 = plot(100*coord_image(binary_image), img(binary_image), 'g:');
                        hl1.LineWidth = 1;
                        hold on;

                    end


                    % For the current noise level, load the image reconstructed using
                    % the Hessian-based Green's approach
                    load([results_directory, 'results_' greens_optimisation_approaches{k},  '_' absorption_maps{2},...
                        '_' num2str(k==1) '_' noise_levels{i} 'db'], 'img_greens');


                    % get the line sound speed values of the image reconstructed using
                    % the Green's approach
                    switch plot_number
                        case 1
                            % along the x=0 axis
                            img = img_greens(ix_image,:);
                        case 2
                            % along the y=0 axis
                            img = img_greens(:,iy_image);
                        case 3
                            % along the fisrt diagonal
                            img = diag(img_greens);
                        case 4
                            % along the second digonal
                            img = diag(fliplr(img_greens));
                    end

                    switch k
                        case 1
                            plot_prop = 'm--';
                        case 2
                            plot_prop = 'b-.';
                        case 3
                            plot_prop = 'r-';
                    end


                    % plot the sound speed values
                    hl2 = plot(100*coord_image(binary_image), img(binary_image), plot_prop);
                    hl2.LineWidth = 1;
                    hold on;

                end


                lg = legend('Ground truth', 'TOF-based', 'Born', 'Hessian-based', 'Hessian-free');
                lg.FontSize = 14;
                lg.Location = 'north';


                % get the file name for the image reconstructed using the Green's approach
                file_name_plot = ['8', subfigures(3*(plot_number-1)+i)];
                
                % get the x axis label
                switch plot_number
                    case 1
                        xlabel('y [cm]');
                    case 2
                        xlabel('x [cm]');
                    case 3
                        xlabel('r [cm]');
                    case 4
                        xlabel('r [cm]');
                end

                % get the y axis label
                ylabel(' c [ms^{-1}]');

                set(gca,'FontSize', 14)


                % save the image reconstructed using the Green's approach
                saveas(hl, [images_directory, 'Fig' file_name_plot '.fig']);
                saveas(hl, [images_directory, 'Fig' file_name_plot '.png']);
                saveas(hl, [images_directory, 'Fig' file_name_plot '.tiff']);
                saveas(hl, [images_directory, 'Fig' file_name_plot '.eps'], 'epsc');
                hold off;

            end

        end

end






