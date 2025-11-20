% Example for visualizing the expermental results from the Rochester
% datasets
%
%%=========================================================================
% This script loads the results associated with Rochester datasets from the
% results_directory, and visualize them in the same format as in the
% paper.
% Depending on the chosen experiment_number in the ''GET THE TASK'' Section
% in this script, i.e., 1,2,3 or 4, the following example scripts  must be
% run be run before running this script, with setting save_results = true
% in ''GET THE TASK'' Section in the associated main example scripts so
% that the results are stored in the defined results_directory.
% -example_2d_real_data_vsx.m (with experiment_number 1 or 2),
% -example_2d_real_data_benigncyst.m (with experiment_number 3),
% -example_2d_real_data_malignancy.m (with experiment_number 4),
%
% The variable emitter_downsampling_rate should be set the same as chosen
% in the associated main example scripts.

% author: Ashkan Javaherian
% date:            - 07.02.2024
% last update:     - 15.07.2025
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian


% run 3 c's
clear all
close all
clc

% run the startup
startup_simulation_ust;


%%=========================================================================
% GET THE TASK
%==========================================================================

% set the Boolean whether the created be saved or not
save_figures = true;

% get the experiment number
experiment_number = 4;

%%=========================================================================
% The parameter which should be set the same as chosen in the associated main
% example script
%==========================================================================
% get the downsampling rate for the emitters
emitter_downsampling_rate = 4; %1,2,4

%%=========================================================================
% THE PARAMETERS WHICH ARE ALWAYS FIXED
%==========================================================================

% get the frquency level for plotting
if experiment_number == 4
    level_plot = 7;
else
    level_plot = 10;
end

% get the final frequency level
level_final = 10;




switch experiment_number
    case {1, 2}

        % get the frequency range
        frequency_range = [0.35, 0.95] * 1e6;
        % get the number of freqnecy
        num_frequency = 70;
    case 3

        % get the frequency range
        frequency_range = [0.30, 1.25] * 1e6;
        % get the number of freqnecy
        num_frequency = 100;
    case 4

        % get the frequency range
        frequency_range = [0.20, 1.25] * 1e6;
        % get the number of freqnecy
        num_frequency = 100;
end


% get the color range for displaying the sound speed map
switch experiment_number

    case {1,2}

        % get range of colorbar for displaying the reconstructed images
        color_range = [1450, 1550];
        % get spacing of colorbar for displaying the reconstructed images
        color_spacing = 20;
    case {3,4}

        % get range of colorbar for displaying the reconstructed images
        color_range = [1350, 1600];
        % get spacing of colorbar for displaying the reconstructed images
        color_spacing = 50;
end


% get the lateral and axial axes for displaying sound speed profiles
switch experiment_number
    case 1
        x_slice = -0.015;
        y_slice= 0.01;
    case 2
        x_slice = -0.025;
        y_slice= 0.015;
    case 3
        x_slice = 0.02;
        y_slice = 0.02;
    case 4
        x_slice = 0.01;
        y_slice= 0;
end

% get the optimisation approach used for reconstructing images using the ray
% approaximation to the heterogeneous Green's function
greens_optimisation_approach = 'backprojection';

%%=========================================================================
% GET THE DIRECTORIES FOR LOADING THE RESULTS AND STORING THE CREATED
% FIGURES
%==========================================================================
% get the path for the results
results_directory = ['results/experiment_rochester/experiments/'...
    'experiment' num2str(experiment_number) '/'];

% get the path for storing the figures
figures_directory = 'results/experiment_rochester/figures/';

% make the directory for storing the figure
makeDirectory(figures_directory);

% get the path for the file name of the ultrasound data
switch experiment_number
    case 1
        filename = 'VSX_YezitronixPhantom1';
    case 2
        filename = 'VSX_YezitronixPhantom2';
    case 3
        filename = 'BenignCyst';
    case 4
        filename = 'Malignancy';
end

% get the path for the image reconstruced using a  frequncy-domain
% full-wave inversion approach based on solving the Helmholtz wave
% equation
% This image is not an actual ground truth, but is used as a benchmark for
% comparison with image reconstructed using our ray-Born inversion approach.
switch experiment_number
    case 1
        fullwave_directory = 'results/experiment_rochester/full-wave/YezitronixPhantom1FinalSoS';
    case 2
        fullwave_directory = 'results/experiment_rochester/full-wave/YezitronixPhantom2FinalSoS';
    case 3
        fullwave_directory = 'results/experiment_rochester/full-wave/BenignCystFinalSoS';
    case 4
        fullwave_directory = 'results/experiment_rochester/full-wave/MalignancyFinalSoS';
end

% get the data path
data_path = 'data/experiment_rochester/';

%%=========================================================================
% LOAD THE TRANSDUCERS' POSITIONS
%==========================================================================
% load the transducers' positions
load([data_path, filename, '.mat'], ...
    'transducerPositionsXY');

% get the emitters' positions
emitter.positions_real = transducerPositionsXY;

% get the receivers' positions
receiver.positions_real = transducerPositionsXY;

%%=========================================================================
% LOAD AND DISPLAY THE IMAGE RECONSTRUCTED USING THE FULL-WAVE INVERSION
%==========================================================================
% load the image reconstrycted using the full-wave approach from its
% directory
if  ~strcmp(fullwave_directory(end-8), '1')
    load(fullwave_directory, 'VEL_ESTIM', 'xi','yi');
    VEL_ESTIM_FINAL = VEL_ESTIM;
else
    load(fullwave_directory, 'VEL_ESTIM_FINAL', 'xi','yi');
end
img_fullwave = VEL_ESTIM_FINAL;

% display the image reconstructed using the ToF-based approach
ha = figure; imagesc(xi, yi, img_fullwave, color_range);
colormap(gray); c = colorbar;
c.Label.String = 'c [m s^{-1}]';
c.YTick = color_range(1):color_spacing:color_range(end);
c.FontSize = 12; axis image; % axis off;
hold on;
switch experiment_number
    case 1
        text(-0.025, -0.04, '1', 'color', 'y', 'FontSize', 12);
        text(0.005,  -0.03, '2', 'color', 'y', 'FontSize', 12);
        text(0,  0.02, '3', 'color', 'y', 'FontSize', 12);
    case 2
        text(-0.0475, -0.04, '4', 'color', 'y', 'FontSize', 12);
        text(-0.03,  -0.0375, '5', 'color', 'y', 'FontSize', 12);
        text(-0.005,  -0.045, '6', 'color', 'y', 'FontSize', 12);
        text(0.0175,  -0.04, '7', 'color', 'y', 'FontSize', 12);
        text(-0.005,  0.015, '8', 'color', 'y', 'FontSize', 12);
    otherwise
end
hold on;
scatter(emitter.positions_real(1, :), emitter.positions_real(2, :), 1, 'white');
scatter(receiver.positions_real(1, :), receiver.positions_real(2, :), 1, 'white');
xticks(1/100* [-10,-5,0,5,10]); yticks(1/100 *[-10,-5,0,5,10]); axis on;
xlabel('Lateral [m]'); ylabel('Axial [m]'); set(gca, 'Fontsize', 12);
axis image;

if save_figures

    % store the image reconstructed using the full-wave inversion approach
    saveas(ha, [figures_directory, 'Fig' num2str(experiment_number) ...
        'a' '.fig']);
    saveas(ha, [figures_directory, 'Fig' num2str(experiment_number) ...
        'a' '.png']);
    saveas(ha, [figures_directory, 'Fig' num2str(experiment_number) ...
        'a' '.tiff']);
    saveas(ha, [figures_directory, 'Fig' num2str(experiment_number) ...
        'a' '.eps'], 'epsc');
end
%%=========================================================================
% LOAD AND DISPLAY THE IMAGE RECONSTRUCTED USING THE TOF-BASED APPROACH
%==========================================================================
% load the time-of-flight-based image, i.e., initial guess for image
% reconstruction approach using ray approximation to the Green's function
load([results_directory, 'results_backprojection_e'...
    num2str(emitter_downsampling_rate) '.mat'] , 'img_tof', 'recon_grid');

% get the grid points inside the detection ring
msk = sqrt(recon_grid.x.^2 + recon_grid.y.^2) < mean(vecnorm(emitter.positions_real));

if ~all(experiment_number -[1,2])

    % set the sound speed outside the detection ring the sound speed of water
    img_tof(img_tof==1490) = 1480;
end

% align the image reconstructed using the ToF-based and ray-Born approach
% with the image reconstructed using the full-wave approach
img_tof_t = img_tof.';


% display the image reconstructed using the ToF-based approach
hb = figure; imagesc(recon_grid.x_vec, recon_grid.y_vec, img_tof_t, color_range);
colormap(gray); c = colorbar;
c.Label.String = 'c [m s^{-1}]';
c.YTick = color_range(1):color_spacing:color_range(end);
c.FontSize = 12; axis image; % axis off;
hold on;
scatter(emitter.positions_real(1, :), emitter.positions_real(2, :), 1, 'white');
scatter(receiver.positions_real(1, :), receiver.positions_real(2, :), 1, 'white');
xticks(1/100* [-10,-5,0,5,10]); yticks(1/100 *[-10,-5,0,5,10]); axis on;
xlabel('Lateral [m]'); ylabel('Axial [m]'); set(gca, 'Fontsize', 12);
axis image;

if save_figures

    % store the image reconstructed using the ToF-based approach
    saveas(hb, [figures_directory, 'Fig' num2str(experiment_number) ...
        'b' '.fig']);
    saveas(hb, [figures_directory, 'Fig' num2str(experiment_number) ...
        'b' '.png']);
    saveas(hb, [figures_directory, 'Fig' num2str(experiment_number) ...
        'b' '.tiff']);
    saveas(hb, [figures_directory, 'Fig' num2str(experiment_number) ...
        'b' '.eps'], 'epsc');
end
%%=========================================================================
% LOAD AND DISPLAY THE IMAGE RECONSTRUCTED USING THE RAY-BORN APPROACH
%==========================================================================

% get an empty cell for the images reconstructed after every 10 frequency iterations
img_greens_fv = cell(level_final, 1);

for  levels = 1:level_final

    % load the results after the current 10 frequency iterations from the given
    % results directory
    load([results_directory, 'results_' greens_optimisation_approach,...
        '_e' num2str(emitter_downsampling_rate) num2str(levels)], 'img_greens');

    if ~all(experiment_number -[1,2])

        % set the sound speed outside the detection ring the sound speed of water
        img_greens(img_greens==1490) = 1480;

    end

    % align the image reconstructed using the ToF-based and ray-Born approach
    % with the image reconstructed using the full-wave approach
    img_greens_t = img_greens.';

    % get the reconstructed image after every 10 frequency iterations
    hc = figure; imagesc(recon_grid.x_vec, recon_grid.y_vec, img_greens_t, color_range);
    img_greens_fv{levels} = img_greens_t;
    colormap(gray); c = colorbar;
    c.Label.String = 'c [m s^{-1}]';
    c.YTick = color_range(1):color_spacing:color_range(end);
    c.FontSize = 12; axis image; % axis off;
    hold on;
    switch experiment_number
        case 1
            text(-0.025, -0.04, '1', 'color', 'y', 'FontSize', 12);
            text(0.005,  -0.03, '2', 'color', 'y', 'FontSize', 12);
            text(0,  0.02, '3', 'color', 'y', 'FontSize', 12);
        case 2
            text(-0.0475, -0.04, '4', 'color', 'y', 'FontSize', 12);
            text(-0.03,  -0.0375, '5', 'color', 'y', 'FontSize', 12);
            text(-0.005,  -0.045, '6', 'color', 'y', 'FontSize', 12);
            text(0.0175,  -0.04, '7', 'color', 'y', 'FontSize', 12);
            text(-0.005,  0.015, '8', 'color', 'y', 'FontSize', 12);
        otherwise
    end
    hold on;
    scatter(emitter.positions_real(1, :), emitter.positions_real(2, :), 1, 'white');
    scatter(receiver.positions_real(1, :), receiver.positions_real(2, :), 1, 'white');
    xticks(1/100* [-10,-5,0,5,10]); yticks(1/100 *[-10,-5,0,5,10]); axis on;
    xlabel('Lateral [m]'); ylabel('Axial [m]'); set(gca, 'Fontsize', 12);
    axis image;
    if save_figures && levels == level_plot

        % save the images reconstructed by the Ray-Born image
        % reconstruction approach
        saveas(hc, [figures_directory, 'Fig' num2str(experiment_number) 'c' '.fig'] );
        saveas(hc, [figures_directory, 'Fig' num2str(experiment_number) 'c' '.png'] );
        saveas(hc, [figures_directory, 'Fig' num2str(experiment_number) 'c' '.tiff'] );
        saveas(hc, [figures_directory, 'Fig' num2str(experiment_number) 'c' '.eps'] , 'epsc');
    end

end


% the vector of frequencies [MHz]
freq = 1e-6 * linspace(frequency_range(1),...
    frequency_range(2), num_frequency);

% get the index of the last frequency for plotting the L2 norm of
% the update direction
freq_plot_index = level_plot/level_final * num_frequency;

% truncate the frequencies used for plotting the L2 norm of the update
% directions
freq = freq(1:freq_plot_index);

% For updating the L2 norm of the update direction
% load the outputs for the last frequency level
load([results_directory, 'results_' greens_optimisation_approach,...
    '_e' num2str(emitter_downsampling_rate) num2str(levels)], 'out_greens_tot');

% get the norm of update direction for all frequency levels
gn = out_greens_tot.update_norm(1:freq_plot_index);
%gn = out_greens_tot.gradient_norm;


hd = figure;plot(freq, gn, 'k');
xlim([freq(1), 1.1*freq(end)]);
xticks(freq(1):0.1:1.1*freq(end));
xlabel('Frequency [MHz]'); ylabel('L2 norm of update direction');
set(gca, 'Fontsize', 12);

if save_figures

    % store the L2 norm of the update directions for each frequnecy level
    saveas(hd, [figures_directory, 'Fig' num2str(experiment_number) 'd' '.fig'] );
    saveas(hd, [figures_directory, 'Fig' num2str(experiment_number) 'd' '.png'] );
    saveas(hd, [figures_directory, 'Fig' num2str(experiment_number) 'd' '.tiff'] );
    saveas(hd, [figures_directory, 'Fig' num2str(experiment_number) 'd' '.eps'] , 'epsc');
end


%%=========================================================================
% CREATE THE LINE PLOTS OF THE RECONSTRUCTE SOUND SPEED IMAGES
%==========================================================================
% get the reconstruction grid used for the full-wave image reconstruction approach
recon_grid_full = [];
[recon_grid_full.x, recon_grid_full.y] = ndgrid(xi, yi);
% get the x vector and y vector for the reconstruction grid used for
% full-wave image reconstruction
recon_grid_full.x_vec = xi;
recon_grid_full.y_vec = yi;

% interpolate the image reconstructed using the full-wave approach onto the
% grid used for the ray-Born image reconstruction approach
img_fullwave_coarse = interpn(recon_grid_full.x, recon_grid_full.y, img_fullwave, ...
    recon_grid.x, recon_grid.y);

% get the imag reconstructed at the frequency level level_plot, i.e., the 
% freqency level chosen for visualizing the reconstructed images and plotting 
img_greens_t = img_greens_fv{level_plot};

%  get the signed radial axis for the grid used for the full-Wave image
%  reconstruction
r_axis_full = sign(recon_grid_full.x_vec) .* sqrt(recon_grid_full.x_vec.^2 + recon_grid_full.y_vec.^2);

%  get the signed radial axis for the grid used for the Ray-Born image
%  reconstruction
r_axis = sign(recon_grid.x_vec) .* sqrt(recon_grid.x_vec.^2 + recon_grid.y_vec.^2);

% plot the image reconstructed using the full-wave inversion approach along the first diagonal
he = figure; plot(r_axis_full, diag(img_fullwave), 'b', 'LineWidth', 2); hold on;
% plot the image reconstructed using the ToF-based inversion approach along the first diagonal
plot(r_axis,diag(img_tof_t), 'g', 'LineWidth', 2);
% plot the image reconstructed using the Ray-Born inversion approach along the first diagonal
plot(r_axis,diag(img_greens_t), 'r', 'LineWidth', 2);hold on;
% Show the transducer positions in the plot
xl1 = xline(-0.11,'-.k','Transducer');
xl1.LabelVerticalAlignment = 'bottom';
xl1.LabelHorizontalAlignment = 'center';
xl2 = xline(0.11,'-.k','Transducer');
xl2.LabelVerticalAlignment = 'bottom';
xl2.LabelHorizontalAlignment = 'center';
yl1 = yline(1480,'--','            water');
yl1.LabelVerticalAlignment = 'middle';
yl1.LabelHorizontalAlignment = 'left';
xlim([-0.13, 0.13]);
legend('Full-wave','ToF-based','Ray-Born');
xlabel('Radial [m]'); ylabel('Sound speed [m s^{-1}]');
set(gca, 'Fontsize', 12);

if save_figures

    % store the line plot (along the first diagonal) of the recontructed sound speed image
    saveas(he, [figures_directory, 'Fig'  num2str(experiment_number) 'e' '.fig'] );
    saveas(he, [figures_directory, 'Fig'  num2str(experiment_number) 'e' '.png'] );
    saveas(he, [figures_directory, 'Fig'  num2str(experiment_number) 'e' '.tiff'] );
    saveas(he, [figures_directory, 'Fig'  num2str(experiment_number) 'e''.eps'] , 'epsc');
end


% plot the image reconstructed using the full-wave inversion approach along the second diagonal
hf = figure; plot(r_axis_full, fliplr(diag(fliplr(img_fullwave))), 'b', 'LineWidth', 2); hold on;
% plot the image reconstructed using the ToF-based inversion approach along the second diagonal
plot(r_axis,fliplr(diag(fliplr(img_tof_t))), 'g', 'LineWidth', 2);
% plot the image reconstructed using the Ray-Born inversion approach along the second diagonal
plot(r_axis,fliplr(diag(fliplr(img_greens_t))), 'r', 'LineWidth', 2); hold on;
% Show the transducer positions within the plot
xl1 = xline(-0.11,'-.k','Transducer');
xl1.LabelVerticalAlignment = 'bottom';
xl1.LabelHorizontalAlignment = 'center';
xl2 = xline(0.11,'-.k','Transducer');
xl2.LabelVerticalAlignment = 'bottom';
xl2.LabelHorizontalAlignment = 'center';
yl1 = yline(1480,'--','            water');
yl1.LabelVerticalAlignment = 'middle';
yl1.LabelHorizontalAlignment = 'left';
xlim([-0.13, 0.13]);
legend('Full-wave','ToF-based','Ray-Born');
xlabel('Radial [m]'); ylabel('Sound speed [m s^{-1}]');
set(gca, 'Fontsize', 12);

if save_figures

    % store the line plot (along the second diagonal) of the recontructed sound speed image
    saveas(hf, [figures_directory, 'Fig'  num2str(experiment_number) 'f' '.fig'] );
    saveas(hf, [figures_directory, 'Fig'  num2str(experiment_number) 'f' '.png'] );
    saveas(hf, [figures_directory, 'Fig'  num2str(experiment_number) 'f' '.tiff'] );
    saveas(hf, [figures_directory, 'Fig'  num2str(experiment_number) 'f''.eps'] , 'epsc');
end

% get the index of x=x_slice axis for the grid used for the full_Wave image
% reconstruction
[~, ixf] = min(abs(recon_grid_full.x_vec - x_slice));

% ix = 118;
% get the index of x=x_slice axis for the grid used for the Ray-Born image
% reconstruction
[~, ix] = min(abs(recon_grid.x_vec - x_slice));

% get the y axis
y_axis = recon_grid.y_vec;
yf_axis = recon_grid_full.y_vec;

% plot the image reconstructed using the full-wave inversion approach along
% x=0 axis
hg = figure; plot(yf_axis, img_fullwave(ixf,:), 'b', 'LineWidth', 2); hold on;
% plot the image reconstructed using the ToF-based inversion approach along
% x=0 axis
plot(y_axis,img_tof_t(ix,:),'g', 'LineWidth', 2);
% plot the image reconstructed using the ray-Born inversion approach along
% x=0 axis
plot(y_axis,img_greens_t(ix,:), 'r', 'LineWidth', 2);hold on;
% Show the transducer positions within the plot
xl1 = xline(-0.11,'-.k','Transducer');
xl1.LabelVerticalAlignment = 'bottom';
xl1.LabelHorizontalAlignment = 'center';
xl2 = xline(0.11,'-.k','Transducer');
xl2.LabelVerticalAlignment = 'bottom';
xl2.LabelHorizontalAlignment = 'center';
yl1 = yline(1480,'--','            water');
yl1.LabelVerticalAlignment = 'middle';
yl1.LabelHorizontalAlignment = 'left';
xlim([-0.13, 0.13]);

legend('Full-wave','ToF-based','Ray-Born');
xlabel('Lateral [m]'); ylabel('Sound speed [m s^{-1}]');
set(gca, 'Fontsize', 12);

if save_figures

    % store the line plot (along x=x_slice axis) of the recontructed sound speed image
    saveas(hg, [figures_directory, 'Fig'  num2str(experiment_number) 'g' '.fig'] );
    saveas(hg, [figures_directory, 'Fig'  num2str(experiment_number) 'g' '.png'] );
    saveas(hg, [figures_directory, 'Fig'  num2str(experiment_number) 'g' '.tiff'] );
    saveas(hg, [figures_directory, 'Fig'  num2str(experiment_number) 'g''.eps'] , 'epsc');
end

% iy = 118;

% get the index of y=y_slice axis for the grid used for the full-Wave image
% reconstruction
[~, iyf] = min(abs(recon_grid_full.y_vec - y_slice));

% iy = 118;

% get the index of y=y_slice axis for the grid used for the ray-Born image
% reconstruction
[~, iy] = min(abs(recon_grid.y_vec - y_slice));

% get the x axis
x_axis = recon_grid.x_vec;
xf_axis = recon_grid_full.x_vec;

% plot the image reconstructed using the full-wave inversion approach along
% y=0 axis
hg = figure; plot(xf_axis, img_fullwave(:, iyf), 'b', 'LineWidth', 2); hold on;
% plot the image reconstructed using the ToF-based inversion approach along
% y=0 axis
plot(x_axis,img_tof_t(:,iy),'g', 'LineWidth', 2);
% plot the image reconstructed using the ray-Born inversion approach along
% y=0 axis
plot(x_axis,img_greens_t(:,iy), 'r', 'LineWidth', 2);hold on;
% Show the transducer positions within the plot
xl1 = xline(-0.11,'-.k','Transducer');
xl1.LabelVerticalAlignment = 'bottom';
xl1.LabelHorizontalAlignment = 'center';
xl2 = xline(0.11,'-.k','Transducer');
xl2.LabelVerticalAlignment = 'bottom';
xl2.LabelHorizontalAlignment = 'center';
yl1 = yline(1480,'--','            water');
yl1.LabelVerticalAlignment = 'middle';
yl1.LabelHorizontalAlignment = 'left';
xlim([-0.13, 0.13]);
legend('Full-wave','ToF-based','Ray-Born');
xlabel('Axial [m]'); ylabel('Sound speed [m s^{-1}]');
set(gca, 'Fontsize', 12);

if save_figures

    % store the line plot (along y=y_slice axis) of the recontructed sound speed image
    saveas(hg, [figures_directory, 'Fig'  num2str(experiment_number) 'h' '.fig'] );
    saveas(hg, [figures_directory, 'Fig'  num2str(experiment_number) 'h' '.png'] );
    saveas(hg, [figures_directory, 'Fig'  num2str(experiment_number) 'h' '.tiff'] );
    saveas(hg, [figures_directory, 'Fig'  num2str(experiment_number) 'h''.eps'] , 'epsc');
end

% get the relative discrepancty between the ToF-based and full-wave reconstructed images
relative_discrepancy_tof = norm(img_tof_t(msk)-img_fullwave_coarse(msk))/norm(1480-img_fullwave_coarse(msk))*100;

% display the relative error
disp(['The percantage relative discrepancty between the tof-based and full-wave reconstructed image is:' num2str(relative_discrepancy_tof, '%2.5f')])


for levels = 1:level_final

    % get the reconstructed image after the current 10-frequency levels
    img_greens_t = img_greens_fv{levels};

    % get the relative discrepancty between the scattering-corrected ray-based and full-wave reconstructed images
    relative_discrepancy_greens(levels) = norm(img_greens_t(msk)-img_fullwave_coarse(msk))/norm(1480-img_fullwave_coarse(msk))*100;

    % display the relative error
    disp(['The percentage relative discrepancty between the scattering-corrected ray-based and full-wave reconstructed image' 'at frequency level' ...
        'is:' num2str(relative_discrepancy_greens(levels), '%2.5f')])

end






