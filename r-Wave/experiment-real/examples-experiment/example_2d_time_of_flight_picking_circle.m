% Example for validating time-of-flight picking algorithm using
% ultrasound data measured in an experimental setting.
%..........................................................................
% 
%%=========================================================================
% EXPERIMENTAL SETTING
%==========================================================================
% The transducers are placed on a ring with radius of 13cm. The time step is
% 40ns, each signal starts at t = 80us after acoustic pulse emission. A single
% emitter placed on the ring is excited by an extric pulse and the induced acoustic 
% waves are measured in time by receivers placed on an arc encompassing 120 degrees of the
% same ring in front of the emitter. The emitter is rotated and is excited and
% the receivers are correspondingly rotated on the ring until the
% excitaions are performed for 2pi Rad.

%%=========================================================================
% LOAD OF UST DATA 
%==========================================================================
% The UST data used in this experiment can be loaded from ''www.zenodo.org''.
%
%
% author: Ashkan Javaherian
% date:            - 05.02.2023
% last update:     - 15.07.2023
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian

clear all
close all
clc

% run startup
startup_simulation_ust;

% get the data path
data_path = 'data/experimental_2D/experiment1/';

% make the directory, if it does not exist
makeDirectory(data_path);

% get the number of dimensions
dim = 2;

% get the radius of the detection ring
detec_radius = 13000e-5;

% get the sound speed [m/s] in water
sound_speed_water = 1480;

%%=========================================================================
% THE PARAMETERS WHICH CAN BE CHANGED BY THE USER.
%%=========================================================================

% get the index of emitter
emitter_index = 10;

% get the index of receiver for analysing the steps of time-of-flight
% picking algorithms
receiver_index = 5;

% An assumption of the maximum sound speed difference [m/s] with water in the medium
% For each emitter-receiver pair,the first arrivals are are picked within a time
% window $t_s + [d/c_max, d/c_min]$, where $t_s$ is the shot (excitation) time,
% $d$ is the distance between emitter and receiver, and
% c_min = sound_speed_water - sound_speed_time_window;
% c_max = sound_speed_water + sound_speed_time_window;
sound_speed_time_window = 80;

% A scalar representing the fraction of the maximum absolute amplitude of
% the signal as the right edge of the time window.
% If chosen nonzero, the right edge of the time window is chosen as
tof_peak_frac_threshold = 0.35;

%%=========================================================================
% GET EMITTERS
%==========================================================================

% get the angular poition of the emitters
emitters_angular_position = deg2rad(0:2:358);

% get the number of meitters
num_emitter = size(emitters_angular_position, 2);

% get an estimated the excitation pulse duration
emitter.pulse_duration =  5.0e-6;

% get the shot time
emitter.shot_time = -80e-6;

% allocate a zero matrix for position of emitters
emitter.positions = zeros(dim, num_emitter);

% get the Cartesian position of emitters
[emitter.positions(1,:), emitter.positions(2,:)] = pol2cart(emitters_angular_position,...
    detec_radius);

% get the rotation indices
emitter.rotation_indices = 1:num_emitter;

%%=========================================================================
% GET RECEIVERS
%==========================================================================

% choose the receivers included in the time-of-flight-based image
% reconstruction. The receivers which face the emitter with bad directivity
% are excluded.
indices = 6:55;

% get the angular position of all receivers
receivers_angle = deg2rad(-58:2:60);

% get the angular position of the used receivers
receivers_angle = receivers_angle(indices);

% get the angular position of receivers for all emitters
receivers_angular_position = emitters_angular_position' + pi + receivers_angle;

% get the number of receivers per emitter (excitation)
num_receiver = size(receivers_angular_position, 2);

% get the struct for the rotating position of receivers
receiver.positions = cell(1, num_emitter);

% get the rotating position of receivers
for ind_emitter = 1: num_emitter
    receiver.positions{ind_emitter} = zeros(dim, num_receiver);
    [receiver.positions{ind_emitter}(1,:), receiver.positions{ind_emitter}(2,:)] =...
        pol2cart(receivers_angular_position(ind_emitter,:), detec_radius);
end

%%=========================================================================
% GET THE MEASURED DATA
%==========================================================================

% get the data measured for object in water
load ([data_path, 'raw_data.mat']);
data = permute(S(:,indices,:), [2, 1, 3]);
data = data(:, :, 1:num_emitter);
clear S;


% get the data measured for only water
load([data_path, 'datam90p90step250avg.mat']);
S(:, end)= [];

% match the exitations for only water with
data_water = repmat(S(15 + indices, :), [1, 1, num_emitter]);
clear S


%%=========================================================================
% GET THE TIME ARRAY
%==========================================================================

% get the time spacing
dt = 4000e-11;

% get the number of time instants in the time array used for measuring data
num_time = size(data, 2);

% get the time array used for measuring data
time_array = cumsum([0, dt * ones(1, num_time-1)]);


%%=========================================================================
% COMPUTE THE TIME-OF-FLIGHTS FOR THE OBJECT-IN-WATER DATA FOR THE CHOSEN EMITTER
% AND FOR ALL RECEIVERS USING THE MAIN FUNCTION
%==========================================================================

% get the pressure data for the chosen emitter
data_water_emitter = data_water(:,:, emitter_index);

% the optional inputs for calculating the tim-of-flights from the measured data
tof_args = {'nWorkerPool', 0, 'SoundSpeedRanges',...
    sound_speed_water + sound_speed_time_window * [-1, 1],...
    'binaries_emitter_receiver', 'distances',...
    'Threshold', tof_peak_frac_threshold};

emitter_single = emitter;
emitter_single.positions = emitter.positions(:, emitter_index);
emitter_single.rotation_indices = emitter.rotation_indices(emitter_index);

% compute the time-of-flights [s] for the only-water data
[~, ~, tof_water] = timeOfFlightPicking([],...
    data_water_emitter, emitter_single, receiver, time_array,...
    emitter_single.rotation_indices, tof_args{:});

% correct the only-water time-of-flights for the excitation time
tof_water = tof_water - emitter.shot_time;

% compute the true only-water TOFs
tof_water_true = 1/sound_speed_water * calculateDistanceEmitterReceiver(...
    emitter_single.positions,...
    receiver.positions, emitter_single.rotation_indices);

% get the index of receivers based on distances to emitter in an ascending
% order
[~, indices_ascend] = sort(tof_water_true);

% plot the only-water TOFs with respect to dis
figure;  hold on;
plot(1:num_receiver, tof_water_true(indices_ascend), 'k'); hold on;
plot(1:num_receiver, tof_water(indices_ascend), 'r');
legend('True TOFs', 'Computed TOFs');


%%=========================================================================
% ANALYSE THE STEPS FOR TIME-OF-FLIGHT PICKING FOR THE TIME TRACE
% ASSOCIATED WITH THE CHOSEN EMITTER-RECEIVER PAIR
%==========================================================================


% get the pressure signal for the chosen emitter-receiver pair
data_water_receiver = data_water_emitter(receiver_index, :);

% plot the original signal for the chosen emitter-receiver pair
figure; plot(1e6 * time_array, data_water_receiver);
xlabel('time [\mu s]'); ylabel('amplitude [a.u.]');

% calculate the distances between the chosen emitter and all receivers
distance_emitter_receivers = calculateDistanceEmitterReceiver(emitter_single.positions,...
    receiver.positions, emitter_single.rotation_indices);

% calculate the minimum and maximum time (in the time_array) for the time
% window. The time window may be later made tighter.
t_min = (distance_emitter_receivers(receiver_index) ./...
    (sound_speed_water + sound_speed_time_window)) + emitter.shot_time;
t_max = (distance_emitter_receivers(receiver_index) ./...
    (sound_speed_water - sound_speed_time_window)) + emitter.shot_time;

% get a Tukey window
tukey = getWin(length(time_array), 'Tukey', 'Param', 0.25).';

% remove the DC offset from the reference (in only water) data for the current
% emitter, and then apply a Tukey window on the time-domain signal.
signal = bsxfun(@times, tukey, removeDataDc(data_water_receiver) );

% get the maximum absolute amplitude
amplitude_max = max(abs(signal));

% plot the zero-DC signal for the chosen emitter-receiver pair
figure; plot(1e6 * time_array, signal);
xlabel('time [\mu s]'); ylabel('amplitude [a.u.]');
axis([0, 1.1 * 1e6 * time_array(end), amplitude_max * [-1, 1]]);

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
figure; plot(1e6 * time_array, signal); hold on;
xline(1e6 * dt * (time_index_start - 1), '--r',{'t_l'}); hold on;
xline(1e6 * dt * (time_index_end - 1), '--r',{'t_r'});
xlabel('time [\mu s]'); ylabel('amplitude [a.u.]');
axis([0, 1.1 * 1e6 * time_array(end), amplitude_max * [-1, 1]]);

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

figure; plot(1e6 * t_array_window, AIC); hold on;
xline(1e6 * (dt * (ixm - 1) + t_array_window(1)), '--r', 'AIC_{min}'); hold on;
xlabel('time [\mu s]'); ylabel('AIC [a.u.]');

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
figure; plot(1e6 * time_array, signal); hold on;
xline(1e6 * tof, '--k',{'first arrival'}); 
xlabel('time [\mu s]'); ylabel('amplitude [a.u.]');
axis([0, 1.1 * 1e6 * time_array(end), amplitude_max * [-1, 1]]);
