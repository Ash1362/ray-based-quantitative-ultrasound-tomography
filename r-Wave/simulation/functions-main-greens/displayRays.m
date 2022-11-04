function displayRays(emitter_positions, receiver_positions, ray_position, ray_time,...
    ray_absorption, ray_position_left, ray_position_right,...   
adjoint_ray_position_left, adjoint_ray_position_right, ray_direction, ray_spacing,...
rayspacing_receiver, plot_directory, varargin)
%DISPLAYRAYS displays the sampled points along the forward and adjoint rays
%   
%
% DESCRIPTION:
%             displayRays displays the sampled points along the forward and adjoint rays
%      
%
% USAGE:
% 
% INPUTS:
%      emitter_positions                       - the dim x num_emitter
%                                                position of the emitter
%      receiver_positions                      - the dim x num_receiver 
%                                                position of the receivers
%      receiver_position
%      ray_position                            - the num_emitter x 1 cell array
%                                                 for the position of the sampled
%                                                 points along the rays
%      ray_time                                 - the num_emitter x 1 cell array
%                                                 for the accumulated time delays
%                                                 on the sampled points along
%                                                 the ray
%      ray_absorption                           - the num_emitter x 1 cell array
%                                                 for the accumulated acoustic 
%                                                 absorption on the sampled
%                                                 points along the rays
%      ray_position_left                        - the position of the
%                                                 sampled points along the left
%                                                 auxiliary ray
%      ray_position_right                       - the position of the
%                                                 sampled points along the right auxiliary
%                                                 ray
%      adjoint_ray_position_left                - the position of the
%                                                 sampled points along the left auxiliary
%                                                 adjoint ray
%      adjoint_ray_position_right               - the position of the
%                                                 sampled points along the right auxiliary
%                                                 adjoint ray
%      ray_direction                            - the perturbed direction
%                                                 of the forward rays
%      ray_spacing                              - a scalar value representing the 
%                                                 spacing [m] of the
%                                                 sampled points along the
%                                                 rays
%      ray_spacing_receiver                     - the spacing [m] of the
%                                                 last two points along the linked rays
%      plot_directory                           - the directory for saving the plots     
%
% OPTIONAL INPUTS:
%      'save_plots'        - Boolean controlling whether the plots are saved
%                            or not. (default:true)
%      'emitter_index'     - the index of emitter on which the diplayed forward rays
%                            are initialised. (default: 1)
%      'receiver_index'    - the index of receiver on which the displayed
%                            adjoint rays are initialised. (default:100)
%      'sample_coeff'      - the coefficint for sampling the points along
%                            the rays
para.save_plots = false;
para.emitter_index = 1;
para.receiver_index = 100;
para.sample_coeff = 8;

  
% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% get the chosen emitter and receiver indices
emitter_index = para.emitter_index;
receiver_index = para.receiver_index;

if rem(para.sample_coeff, 2) > 0  ||  para.sample_coeff<=0 
     error('The coefficient for sampling the rays must be even.')
end
  
% check if any of the auxiliary varibales are empty
if isempty(ray_position_left)
  
    % compute the ray's Jacobian with respect to the initial angle 
    % only using the linked (main) rays, not auxiliary rays
    attenuation_geom_method = 'raylinked';
    
    % set the auxiliary method nan
    auxiliary_method = nan;

else
    
    % compute the ray's Jacobian with repect to the initial angle 
    % using the auxiliary rays 
    attenuation_geom_method = 'auxiliary';
    
    if isempty(ray_position_right{1})
        
        % the method used for tracing auxiliary rays is 'paraxial', the right
        % auxiliary rays are not given 
        auxiliary_method = 'paraxial';
    else
        auxiliary_method = 'angle_perturbation';
    end
    
end

% get the x Cartesian position of the sampled points on the rays for the chosen emitter
ray_position_x_emitter = ray_position{emitter_index}(1:2*para.sample_coeff:end, :);

% get the y Cartesian position of the sampled points on the rays for the chosen emitter
ray_position_y_emitter = ray_position{emitter_index}(2:2*para.sample_coeff:end, :);


if strcmp(attenuation_geom_method, 'auxiliary')
    

% remove the nans from the x position of the forward rays initilaised from
% the chosen emitter
nan_binary = isnan(ray_position_x_emitter);
ray_position_x_emitter = ray_position_x_emitter(~nan_binary);

% remove the nans from the y position of the forward rays initilaised from
% the chosen emitter
ray_position_y_emitter = ray_position_y_emitter(~nan_binary);


if strcmp(attenuation_geom_method, 'auxiliary')
ray_position_x_left_emitter = ray_position_left{emitter_index}(1:2*para.sample_coeff:end, :);
ray_position_y_left_emitter = ray_position_left{emitter_index}(2:2*para.sample_coeff:end, :);
ray_position_x_left_emitter = ray_position_x_left_emitter(~nan_binary);
ray_position_y_left_emitter = ray_position_y_left_emitter(~nan_binary);
if strcmp(auxiliary_method, 'angle_perturbation')
ray_position_x_right_emitter = ray_position_right{emitter_index}(1:2*para.sample_coeff:end, :);
ray_position_y_right_emitter = ray_position_right{emitter_index}(2:2*para.sample_coeff:end, :);
ray_position_x_right_emitter = ray_position_x_right_emitter(~nan_binary);
ray_position_y_right_emitter = ray_position_y_right_emitter(~nan_binary);
end
end


% get the time delays with respect to the position of the chosen emitter
ray_time_emitter = ray_time{emitter_index};
ray_time_emitter = ray_time_emitter(1:para.sample_coeff:end, :);

% remove the nans from the accumulated time delays along the forward rays
% initilaised from the chosen emitter
ray_time_emitter = ray_time_emitter(~nan_binary);


% compute the adjoint rays by reversing the forward rays
[ray_position_adjoint, ray_time_adjoint, ~,...
    ray_position_left_adjoint, ray_position_right_adjoint] = calcRaysAdjoint(...
    ray_position, ray_time, ray_absorption, adjoint_ray_position_left,...
    adjoint_ray_position_right, ray_direction, ray_spacing, rayspacing_receiver);


% get the parameters on the adjoint rays for the chosen receiver

% get the x-coordinate poistion for the chosen receiver
ray_position_x_receiver = ray_position_adjoint{receiver_index}(1:2 * para.sample_coeff:end, :);

% get the y-coordinate position for the chosen receiver
ray_position_y_receiver = ray_position_adjoint{receiver_index}(2:2*para.sample_coeff:end, :);

% remove the nans the nans from the x position of the adjoint rays 
nan_binary_adjoint = isnan(ray_position_x_receiver);
ray_position_x_receiver = ray_position_x_receiver(~nan_binary_adjoint);

% remove the nans the nans from the y position of the adjoint rays 
ray_position_y_receiver = ray_position_y_receiver(~nan_binary_adjoint);

if strcmp(attenuation_geom_method, 'auxiliary')

ray_position_x_left_receiver = ray_position_left_adjoint{receiver_index}(1:2*para.sample_coeff:end, :);
ray_position_y_left_receiver = ray_position_left_adjoint{receiver_index}(2:2*para.sample_coeff:end, :);

% remove the nans
ray_position_x_left_receiver = ray_position_x_left_receiver(~nan_binary_adjoint);
ray_position_y_left_receiver = ray_position_y_left_receiver(~nan_binary_adjoint);

switch auxiliary_method
    case 'angle_perturbation'
ray_position_x_right_receiver = ray_position_right_adjoint{receiver_index}(1:2*para.sample_coeff:end, :);
ray_position_y_right_receiver = ray_position_right_adjoint{receiver_index}(2:2*para.sample_coeff:end, :);

% remove the nans
ray_position_x_right_receiver = ray_position_x_right_receiver(~nan_binary_adjoint);
ray_position_y_right_receiver = ray_position_y_right_receiver(~nan_binary_adjoint);

    case 'paraxial'
        
        % convert the perturbed positions to an implicit position such that the
        % paraxial rays can be displayed, if the method for tracing
        % auxiliary rays is set 'paraxial'.
        
        % emitter
        ray_position_x_left_emitter = ray_position_x_emitter + 0.01 *...
            ray_position_x_left_emitter;
        ray_position_y_left_emitter = ray_position_y_emitter + 0.01 *...
            ray_position_y_left_emitter;
        
        % receiver
        ray_position_x_left_receiver = ray_position_x_receiver + 0.01 *...
            ray_position_x_left_receiver;
        ray_position_y_left_receiver = ray_position_y_receiver + 0.01 *...
            ray_position_y_left_receiver;

end

end

% modify the perturbation to the position of the rays such that the
% auxiliary rays can be displayed, if the method for tracing auxiliary ray
% is set 'paraxial'.


   

% get the accumulated time delays on the sampled points the adjoint rays 
% initialised from the chosen receiver.
ray_time_receiver = ray_time_adjoint{receiver_index}(1:para.sample_coeff:end, :);

% remove the nans from the time delays along the adjoint rays
ray_time_receiver = ray_time_receiver(~nan_binary_adjoint);


%%=========================================================================
% DISPLAY THE TIME DELAYS ALONG THE FORWARD RAYS
%==========================================================================
h1= figure();

% display the sampled points on the main rays for the chosen emitter
scatter(vectorise(ray_position_x_emitter),...
    vectorise(ray_position_y_emitter), 10, vectorise(ray_time_emitter));

hold on;
% display the emitter
scatter(emitter_positions(1, emitter_index),...
    emitter_positions(2, emitter_index), 50, 'g', 'filled'); hold on;
if strcmp(attenuation_geom_method, 'auxiliary')
hold on;
scatter(vectorise(ray_position_x_left_emitter),...
    vectorise(ray_position_y_left_emitter), 5, 'g');
hold on;
if strcmp(auxiliary_method, 'angle_perturbation')
scatter(vectorise(ray_position_x_right_emitter),...
    vectorise(ray_position_y_right_emitter), 5, 'b');
end
end
axis image;

title('The sampled main forward rays. The auxiliary ray is shown by green.')
set(gca, 'Fontsize', 12); colorbar;

% add colorbar 
a = colorbar;
ylabel(a, 'Time [s]');

if para.save_plots
saveas(h1, [plot_directory, 'emitter_time_delays' '.fig']);
saveas(h1, [plot_directory, 'emitter_time_delays' '.png']);
saveas(h1, [plot_directory, 'emitter_time_delays' '.tiff']);
saveas(h1, [plot_directory, 'emitter_time_delays' '.eps'], 'epsc');
end

%%=========================================================================
% DISPLAY THE TIME DELAYS ALONG THE ADJOINT RAYS
%==========================================================================

h2 = figure();

% display the sampled points on the main adjoint rays for the chosen receiver
scatter(vectorise(ray_position_x_receiver),...
    vectorise(ray_position_y_receiver), 10, vectorise(ray_time_receiver));

% display the sampled points on the auxiliary adjoint rays for the chosen receiver
if strcmp(attenuation_geom_method, 'auxiliary')
hold on;
scatter(vectorise(ray_position_x_left_receiver),...
    vectorise(ray_position_y_left_receiver), 5, 'g');
hold on;
if strcmp(auxiliary_method, 'angle_perturbation')
scatter(vectorise(ray_position_x_right_receiver),...
    vectorise(ray_position_y_right_receiver), 5, 'b');
end
end
hold on;


% display the receiver
%scatter(receiver.positions(1,:), receiver.positions(2, :), 10, 'k', 'filled');

% display the receiver
scatter(receiver_positions(1, receiver_index),...
    receiver_positions(2, receiver_index), 50, 'r', 'filled'); hold on;
axis image;

title('The sampled main and auxiliary adjoint rays. The auxiliary ray is shown by green.')

% add colorbar 
a = colorbar;
ylabel(a, 'Time [s]');

set(gca, 'Fontsize', 12);


if para.save_plots
saveas(h2, [plot_directory, 'receiver_time_delays' '.fig']);
saveas(h2, [plot_directory, 'receiver_time_delays' '.png']);
saveas(h2, [plot_directory, 'receiver_time_delays' '.tiff']);
saveas(h2, [plot_directory, 'receiver_time_delays' '.eps'], 'epsc');
end  
        

end