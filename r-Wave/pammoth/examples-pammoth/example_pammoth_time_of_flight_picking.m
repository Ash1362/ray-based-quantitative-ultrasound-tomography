% Example:  Time-of-flight picking from pressure time series measured 
% in a 3D transmission ultrasound setting
%
% This example computes the first-arrival time of the pressure time series
% measured in a 3D transmission ultrasound setting
%
%
% author: Ashkan Javaherian
% date:            - 04.08.2020
% last update:     - 21.10.2021
%
%
% This script is part of the r-Wave Toolbox.
% Copyright (C) 2021 Ashkan Javaherian



clear all
close all
clc

% define the paths
startup_pammoth_ust;

%% ========================================================================
% GET THE DATA IDS
%==========================================================================
data_ids = [];

% the main directory, which must be fixed for all UST data
data_ids.main_path = main_path;

% get the object, the type of object for which the UST data have been
% measured. (check the m-file function 'setPathGetInfos')
data_ids.object = 'us_measurements';

% get the ID of the data
data_ids.id = 'us_measurements/';

% get the date on which the data are measured
data_ids.date = '20201123_';
      
% get the measurement number as a string
data_ids.measurement_number = '00'; 

% get the cup size
data_ids.cup_size = 8;

% get the main path for saving the user's data. This should be first defined
% by the user, and will be fixed for all the reconstructions.
data_ids.main_path_user_data = 'pammoth/data/';

% the number of workers
num_worker_pool = 16; 

% get the z ccordinate of emitter. The emitter closest to the chosen
% z coordinate will be chosen
emitter_coord_z = -0.01;

% get the distance of emitter-receiver pair for analysing the
% time-of-flight picking. The receiver closest to the chosen will be
% selected.
receiver_distance = 0.22;

% get the directory for storing the plots
directory_figures = 'plots/tof_picking/';

% make the directory for the plots
makeDirectory(directory_figures);

% Boolean cpontrolling whether the figures are saved or not
save_figures = true;

%% ========================================================================
% GET THE ALWAYS FIXED PARAMETERS
%==========================================================================
% get the sampling rate [Hz] for the excitation pulse
time_info.emit_frequency = 5.0e7;

% get the sampling rate [Hz] for the measured ultrasound signals 
time_info.us_frequency = 2.5e7;

% get the time index [a.u.] corresponding to the shot time of the signal,
% according to the pdf file provided by the PA imaging company
time_info.shot_time_index = 637;

% choose the data set for only-water data
only_water_measurement_date = '05-2021';

% the amplitude threshold for chooing the right edge of the time window
% For some US data sets including this, the threshold is set zero. If set
% zero, the both edges  will be chosen based
tof_peak_frac_threshold = 0;

% get the total angular positions for the stored data
angular_position_total = 3;

% get the last angular position for time-of-flight picking
angular_position_last = 1;

% get the radius [m] of the bowl, the hemi-spherical surface on which the 
% transducers are positioned.
radius_bowl = 13e-2;

% a minimum distance for the emitter-receiver pairs for which time-of-flights
% are computed.
minimum_distance = 18e-2;

% get the index of the time series used for time-of-flight picking 
time_index_end = 6000;

% A range of sound speed for choosing a time window for maximum possible
% varitation of first-arrival times for the object-in-water data and only-water data.
% This time window is chosen based on assuming a minimum-maximum
% homogeneous sound speed inside the bowl. (Default = [-50, 50] for the
% object-in-water data, and [-20, 20] for the only water data.) 
% In general, the better the information about the shot time of the
% excitation pulse, namely US shot time, be available, the narrower the ranges 
% for choosing a time window for the first-arrival of the signals could be. 
sound_speed_ranges_diff = 50 * [-1, +1];    %  [m/s]

%% ========================================================================
% GET THE DATA PATH AND THE TABLE
%==========================================================================
% get the path associated with the object-in-water data, and the table for
% the corresponding angular positions, the temperatures for the water 
% encompassing the object during the rotations, and the cup size.
[data_paths_object, table_array_position, water_object_temperatures,...
    ~] = setPathGetInfos(data_ids, angular_position_total);

%%=========================================================================
%  DEFINE THE PATH FOR STORING AND LOADING THE DATA FOR THE OBJECT ON THE LOCAL
%  MACHINE
%==========================================================================
% define a path for the user on the local machine
user_data_paths = [];
user_data_paths.main_directory = data_ids.main_path_user_data;
user_data_paths.directory = data_paths_object.directory;
user_data_paths.data_name = data_paths_object.data_name;

% convert the user data path from struct to a string
user_data_path = [user_data_paths.main_directory, user_data_paths.directory,...
    user_data_paths.data_name '/'];

%% ========================================================================
% GET POSITION OF EMITTERS/RECEIVERS
%==========================================================================

% get the position of emitter-receiver pairs
position_args = {'do_calculate_tofs', false};
[emitter, receiver, emit_shot_index, us_shot_index, emit_firing_transducer_index,...
    rotation_indices] =...
     getPositions(data_paths_object, table_array_position, radius_bowl,...
     user_data_paths, position_args{:});

% get the emitter positions for the first angular position
emitter.positions = emitter.positions(:, rotation_indices <= angular_position_last);

% get the receiver positions for the first angular position and convert to
% a matrix
receiver.positions = cell2mat(receiver.positions(1:angular_position_last));

% get the the indices of excitations for the first angular position
rotation_indices = rotation_indices(rotation_indices <= angular_position_last);

% get the number of dimensions
dim = size(emitter.positions, 1);
 
% get the number of emitters
num_emitter = size(emitter.positions, 2);

% get the number of receivers
num_receiver = size(receiver.positions, 2);


%% ========================================================================
% GET THE ONLY-WATER DATA
%==========================================================================
% load the only-water data, and the sound speed for only-water data.
% This data is availabe for only one full rotation over all emitters, because
% the tofs in only water are independent of the angular positions, and are only
% affeced by the distance of emitter-receiver pairs, so the data for one full
% rotation over all emitters can be extended to all excitations.
[data_water, ~, ~, sound_speed_only_water] = loadDataOnlyWater(data_paths_object,...
    only_water_measurement_date, num_receiver);

% truncate the end part of the time series using last_time_index
data_water = data_water(:, 1: time_index_end, :);

%% ========================================================================
% GET THE INFORMATION ABOUT THE WATER FOR OBJECT-IN-WATER DATA
%==========================================================================
% get the information about the water encompassing the object 
[sound_speed_object_correction_coeff, refractive_background_object_alldata,...
    sound_speed_object_water] = getWaterObject(water_object_temperatures, rotation_indices,...
    sound_speed_only_water, true, num_emitter, num_receiver);

%%=========================================================================
% GET THE INFORMATION OF TIME ARRAY
%==========================================================================
% get the time array
time_array = cumsum([0, (1/time_info.us_frequency) * ones(1, size(data_water, 2) - 1)]);

% get the us shot time based on the pdf provided by the PA imaging company.
% The shot time is fixed for all excitatons, and agrees well with the
% calculated TOF differences.
emitter.shot_time = (time_info.shot_time_index - 1) * (1/time_info.emit_frequency);

% get an approximate pulse duration.
% This is not used when the AIC approach is used for time-of-flight
% picking.
emitter.pulse_duration = 3e-6;


%%=========================================================================
% GET THE OBJECT-IN-WATER DATA
%==========================================================================


% Get the optional inputs for time-of-flight picking function, which are not
% actually optional and are set in the main script
tof_args = {'nWorkerPool', num_worker_pool,...
    'SoundSpeedRanges', 2/5 * sound_speed_ranges_diff + sound_speed_object_water,...
    'SoundSpeedRef', sound_speed_only_water, 'Threshold', tof_peak_frac_threshold,...
    'minimum_distance', minimum_distance};

% compute the time-of-flight sinogram for only-water (reference) data for
% the first angular position
disp('computing time-of-flight sinogram for the only-water data:');
[~, ~, tof_hom_first] = timeOfFlightPicking([], data_water,...
    emitter, receiver, time_array, [], tof_args{:});

% enlrage the the sound speed ranges for time-of-flight picking for the 
% object-in-water data
tof_args{4} = sound_speed_ranges_diff + sound_speed_object_water;

% load the object-in-water data by reading the associated HDF5 file, which has 
% been stored on the defined user's data path
data_object = h5read([user_data_path, 'time_series_early.h5'], '/data_breast',...
    [1, 1, 1], [inf, time_index_end, num_emitter]);

disp('computing time-of-flight sinogram for the object-in-water data:');
% compute the time-of-flight sinogram for object-in-water data for
% the first angular position
[~, tof_het, ~] = timeOfFlightPicking(data_object, [], emitter,...
    receiver, time_array, [], tof_args{:});

% compute the time-of-flights for the object-in-water data by subtracting
% the shot time from first-arrival times
tof_het = tof_het - emitter.shot_time;

% compute the time-of-flights for the only-water data by subtracting
% the shot time from first-arrival times and then correct for the
% relative difference with the water encompassing the object during the
% rotation of the bowl
tof_hom_first = sound_speed_object_correction_coeff .* ...
    (tof_hom_first - emitter.shot_time);

% compute the difference TOFs, the discrepancy of the picked TOFs between
% the object-in-water and only-water data
tof_data = tof_het - tof_hom_first;

% compute the distance between emitter-receiver pairs
distance_emitters_receivers = calculateDistanceEmitterReceiver(...
    emitter.positions, receiver.positions, []);

% calculate the true TOFs between the emitters and receivers for the water
% encompassing the object
tof_true = 1/sound_speed_only_water * ...
    distance_emitters_receivers .*...
    reshape(refractive_background_object_alldata, num_receiver, num_emitter);

%%=========================================================================
% DISPLAY THE PLOTS
%==========================================================================

% get the index of the chosen emitter for plotting the time-of-flights
[~, emitter_index] = min(abs(emitter.positions(3,:) - emitter_coord_z));

% get the distance of the chosen emitter and all receivers and sort based 
% on distances
[distance_sorted, indices_distance_sorted] =...
    sort(distance_emitters_receivers(:, emitter_index));


% display the plot
h1 = figure; plot(1e2 * distance_sorted, 1e6 * tof_true(indices_distance_sorted, emitter_index), 'k'); hold on; 
plot(1e2 * distance_sorted, 1e6 * tof_hom_first(indices_distance_sorted, emitter_index), 'r--'); hold on;
plot(1e2 * distance_sorted, 1e6 * tof_het(indices_distance_sorted, emitter_index), 'b'); 
xlabel('Distances [cm]'); ylabel('TOFs [\mu s]');
axis([1e2 * (minimum_distance - 0.01),...
    1e2 * (2 * radius_bowl + 0.01),...
    1e6/sound_speed_only_water * (minimum_distance - 0.01),...
    1e6/sound_speed_only_water * (2 * radius_bowl + 0.01)]);
legend('T_w (true)', 'T_w',...
    'T_{obj}', 'Location', 'Northwest');

if save_figures
    saveas(h1, [directory_figures, 'fig_1b.fig']);
    saveas(h1, [directory_figures, 'fig_1b.tiff']);
    saveas(h1, [directory_figures, 'fig_1b.png']);
    saveas(h1, [directory_figures, 'fig_1b.eps'], 'epsc');
end

% get the index of receiver, which is closest to the distance to
% the chosen emitter
[~, receiver_index] = min(abs(distance_emitters_receivers(:, emitter_index)-...
receiver_distance ));

% get the included receivers
receiver_indices = distance_emitters_receivers(:, emitter_index)> minimum_distance;

% display the emitters and receivers
h2 = figure; scatter3(100 * emitter.positions(1,:), 100 * emitter.positions(2,:),...
100 * emitter.positions(3,:), 10, 'k'); hold on;
scatter3(100 * emitter.positions(1, receiver_indices), 100 * emitter.positions(2, receiver_indices),...
100 * emitter.positions(3, receiver_indices), 10, 'b'); hold on;
scatter3(100 * emitter.positions(1, emitter_index), 100 * emitter.positions(2, emitter_index),...
100 * emitter.positions(3, emitter_index), 50, 'g','filled'); hold on;
scatter3(100 * receiver.positions(1, receiver_index), 100 * receiver.positions(2, receiver_index),...
100 * receiver.positions(3, receiver_index), 50, 'r','filled'); axis image; view(60, 30); axis image;
xlim([-14, 14]); ylim([-14, 14]); zlim([-14, 0]);
xticks([-10,0,10]); yticks([-10,0,10]); zticks([-10,0,10]);
xlabel('[cm]'); ylabel('[cm]'); zlabel('[cm]');
legend('All transducers', 'Included receivers', 'Emitter', 'Receiver');

if save_figures
    saveas(h2, [directory_figures, 'fig_1a.fig']);
    saveas(h2, [directory_figures, 'fig_1a.tiff']);
    saveas(h2, [directory_figures, 'fig_1a.png']);
    saveas(h2, [directory_figures, 'fig_1a.eps'], 'epsc');
end


%%=========================================================================
% ANALYSE THE STEPS FOR TIME-OF-FLIGHT PICKING FOR THE TIME TRACE
% ASSOCIATED WITH THE CHOSEN EMITTER-RECEIVER PAIR
%==========================================================================

% get the pressure signal for the chosen emitter-receiver pair
data_water_receiver = data_water(receiver_index, :, emitter_index);

% get the maximum absolute amplitude
amplitude_max = max(abs(data_water_receiver));

% plot the original signal for the chosen emitter-receiver pair
h3 = figure; plot(1e6 * time_array, data_water_receiver);
xlabel('Time [\mu s]'); ylabel('Amplitude [a.u.]');
axis([0, 1.1 * 1e6 * time_array(end), amplitude_max * [-1, 1]]);

if save_figures
    saveas(h3, [directory_figures, 'fig_1c.fig']);
    saveas(h3, [directory_figures, 'fig_1c.tiff']);
    saveas(h3, [directory_figures, 'fig_1c.png']);
    saveas(h3, [directory_figures, 'fig_1c.eps'], 'epsc');
end

% get the time spacing [s]
dt = time_array(2)-time_array(1);

% calculate the minimum and maximum time (in the time_array) for the time
% window. The time window may be later made tighter.
t_min = (distance_emitters_receivers(receiver_index, emitter_index) ./...
    (sound_speed_ranges_diff(2) + sound_speed_object_water)) + emitter.shot_time;
t_max = (distance_emitters_receivers(receiver_index, emitter_index) ./...
    (sound_speed_ranges_diff(1) + sound_speed_object_water)) + emitter.shot_time;

% get a Tukey window
tukey = getWin(length(time_array), 'Tukey', 'Param', 0.25).';

% remove the DC offset from the reference (in only water) data for the current
% emitter, and then apply a Tukey window on the time-domain signal.
signal = bsxfun(@times, tukey, removeDataDc(data_water_receiver) );

% get the maximum absolute amplitude
amplitude_max = max(abs(signal));

% plot the original signal for the chosen emitter-receiver pair
%h4 = figure; plot(1e6 * time_array, signal);
%xlabel('time [\mu s]'); ylabel('amplitude [a.u.]');
%axis([0, 1.1 * 1e6 * time_array(end), amplitude_max * [-1, 1]]);

%if save_figures
%    save(h4, [directory_figures, 'fig2.mat']);
%    save(h4, [directory_figures, 'fig2.tiff']);
%    save(h4, [directory_figures, 'fig2.png']);
%    save(h4, [directory_figures, 'fig2.eps'], 'epsc');
%end


% get the window length factors
window_length_factor =  [1.0, nan, 0.25, nan];

% calculate the integer length of the windows for computing the
% first-arrivals.
window_length = round(window_length_factor * emitter.pulse_duration/dt);

%  length of the signal
length_signal = length(signal);

% length of the backward window
length_back = window_length(1);

% length of the long window (if applied)
length_long = window_length(2);

% length of the short window (if applied)
if isfinite(window_length(3))
    length_short = window_length(3);
end

% get the minimum and maximum integer indices for the time window
time_index_start = round(t_min/dt) + 1;
time_index_end = round(t_max/dt) + 1;

if  tof_peak_frac_threshold
    
    % calculate the envelope of the signal
    signal_envelope = abs(signal);
    
    % truncate the envelope of the signal using the starting and ending indices
    % which are calculated above
    signal_envelope = signal_envelope(time_index_start : time_index_end);
    
    % find the first index of the time windowed envelope after which the
    % signal is larger the chosen fraction of the maximum absolute amplitude
    ixm = find(signal_envelope/max(signal_envelope) > tof_peak_frac_threshold, 1, 'first');
    
    % shift the index for getting the time index with respect to time origin
    ixm = ixm + time_index_start - 1;
    
end

% normalise the amplitude of the signal
signal_normalised = signal/max(abs(signal));

if tof_peak_frac_threshold
    
    % find the minimum and maximum indices for the time window for finding
    % the minimum AIC. The end point is ixm, and the starting point is
    % ixm - length_back.
    time_index_start = max (1, ixm - length_back);
    time_index_end = ixm;
    
end

% plot the zero-DC signal for the chosen emitter-receiver pair
%h4 = figure; plot(1e6 * time_array, signal); hold on;
%xline(1e6 * dt * (time_index_start - 1), '--r',{'t_l'}); hold on;
%xline(1e6 * dt * (time_index_end - 1), '--r',{'t_r'});
%xlabel('time [\mu s]'); ylabel('amplitude [a.u.]');
%axis([0, 1.1 * 1e6 * time_array(end), amplitude_max * [-1, 1]]);

%if save_figures
%    save(h4, [directory_figures, 'fig3.mat']);
%    save(h4, [directory_figures, 'fig3.tiff']);
%    save(h4, [directory_figures, 'fig3.png']);
%    save(h4, [directory_figures, 'fig3.eps'], 'epsc');
%end


% truncate the time array for the chosen window
t_array_window = time_array(time_index_start:time_index_end);

% allocate a vector for Akaike-Information-Criterion (AIC)
AIC = zeros(size(t_array_window));

% get the length of the window for computing the AIC valuse
length_window = length(t_array_window);

for i = 1:length_window
    k = i + time_index_start - 1;
    variance_noise = var(signal_normalised(1:k));
    variance_signal = var(signal_normalised(k+1:length_signal));
    AIC(i) = k * log(variance_noise) + (length_signal-1-k) * log(variance_signal);
end

% calculate the minimal AIC value within the time window
[AICm, ixm] = min(AIC);

h4 = figure; plot(1e6 * t_array_window, AIC); hold on;
xline(1e6 * (dt * (ixm - 1) + t_array_window(1)), '--k', 'AIC_{min}',...
    'LabelHorizontalAlignment', 'center'); hold on;
xlabel('Time [\mu s]'); ylabel('AIC [a.u.]');

if save_figures
    saveas(h4, [directory_figures, 'fig_1d.fig']);
    saveas(h4, [directory_figures, 'fig_1d.tiff']);
    saveas(h4, [directory_figures, 'fig_1d.png']);
    saveas(h4, [directory_figures, 'fig_1d.eps'], 'epsc');
end



% select the short window for picking the first arrival
indices = max(1, ixm - length_short): min(length_window, ixm + length_short);

% compute the exponential of discrepancy of the AIC and the calculated minimum AIC
% in the chosen short time window about the minimum AIC
exponal_discrepancy_akaike = exp(-(AIC(indices)-AICm)/2);

% calculate the Akaike weights for each data sample within the time window
weights = exponal_discrepancy_akaike/sum(exponal_discrepancy_akaike);

% calculate a weighted average TOF using the Akaike weights
tof = weights * t_array_window(indices)';

% plot the zero-DC signal for the chosen emitter-receiver pair
h5 = figure; plot(1e6 * time_array, signal_normalised); hold on;
xline(1e6 * tof, '--k',{'t_{fa}'},...
     'LabelHorizontalAlignment', 'center'); 
xline(1e6 * dt * (time_index_start - 1), '--r',{'t_l'},...
    'LabelHorizontalAlignment', 'center',...
    'LabelVerticalAlignment', 'bottom'); hold on;
xline(1e6 * dt * (time_index_end - 1), '--r',{'t_r'},...
    'LabelHorizontalAlignment', 'center',...
    'LabelVerticalAlignment', 'bottom'); hold on;
xlabel('Time [\mu s]'); ylabel('Amplitude [a.u.]');
axis([0, 1.1 * 1e6 * time_array(end), [-1, 1]]);

if save_figures
    saveas(h5, [directory_figures, 'fig_1e.fig']);
    saveas(h5, [directory_figures, 'fig_1e.tiff']);
    saveas(h5, [directory_figures, 'fig_1e.png']);
    saveas(h5, [directory_figures, 'fig_1e.eps'], 'epsc');
end

% plot the zero-DC signal for the inside the time window
h6 = figure; plot(1e6 * time_array, signal_normalised); hold on;
xline(1e6 * tof, '--k',{'t_{fa}'}, 'LabelHorizontalAlignment', 'center'); 
xlabel('Time [\mu s]'); ylabel('Amplitude [a.u.]');
axis([1e6 * t_array_window([1, end]), 0.25 * [-1,+1]]);

if save_figures
    saveas(h6, [directory_figures, 'fig_1f.fig']);
    saveas(h6, [directory_figures, 'fig_1f.tiff']);
    saveas(h6, [directory_figures, 'fig_1f.png']);
    saveas(h6, [directory_figures, 'fig_1f.eps'], 'epsc');
end

