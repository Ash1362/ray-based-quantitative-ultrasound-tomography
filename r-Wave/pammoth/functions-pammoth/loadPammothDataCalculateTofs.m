function [tof_het, tof_hom, amplitude_het, amplitude_hom, elapsed_time] =...
    loadPammothDataCalculateTofs(emitter, receiver, rotation_indices,...
    us_shot_index, emit_shot_index, time_array, data_water, num_us_segment,...
    data_paths, user_data_path, varargin)
%LOADPAMMOTHDATACALCULATETOFS calculates the discrepancy of the first-arrival times
%between two sets of data collected by the Pammoth system
%
%
% DESCRIPTION:
%       loadPammothDataCalculateTofs calculates the discrepancy of the time-of-flights
%       between two data sets, one is measured from only water, and another is measured
%       from the object inside water. The data is loaded, and processed
%       segment-by-segment in order to save memory.
%
% USAGE:
%
%
% INPUTS:
%      emitter                     - a struct that defines the properties of the
%                                    excitation. This includes the
%                                    fields 'positions', 'pulse' and 'pulse_duration'
%      emitter.positions           - 2/3 x N array of Cartesian position of
%                                    the centre of the emitter objects
%      emitter.pulse               - the excitation pulse
%                                    pulses for different emitters are different
%      emitter.pulse_duration      - a scalar indicating the the approximate
%                                    time duration of the excitation pulses,
%                                    which will be used for calculation of the
%                                    first-arrival of the signals
%      emitter.shot_time           - the shot time of the excitation
%                                    pulse (or the ultrasound shot time),
%                                    which can be negative [sec].
%                                    If negative, the time series have been
%                                    started after pulse excitation.
%      receiver.positions          - 2/3 x N array of cartesian position of
%                                    the centre of the receiver objects (a
%                                    fixed measuremnt setting)
%                                    a cell array of length num_rot with
%                                    num_rot the number of angular
%                                    positions containg /3 x N array of cartesian position of
%                                    the centre of the receiver objects for
%                                    each angular poistion (a rotational
%                                    measuremnt setting)
%     rotation_indices             - a vector containing the index of
%                                    angular position for each excitation
%     us_shot_index                - the index of US shot indices (See
%                                    Pammoth measurement protocol)
%     time_array                   - the time array used for data measurement
%                                    [s]
%     data_water                   - the data from water
%     num_us_segment               - the number of data segments for loading
%                                  - the number of data segments for time-of-flight picking
%                                    (for saving memory)
%     data_paths                   - path for loading data
%     user_data_path               - local (user) path for storing data
%
% OPTIONAL INPUTS:
%        'nWorkerPool'          - the number of used workers
%        'Method'               - the method for calculation of the first-arrival
%                                 of the signals. The best method for
%                                 medical ultrasiund is 'Modified_AIC'
%        'SoundSpeedRangesDiff' - a prior assumption for changes in the sound
%                                 speed of water, compared to a reference
%                                 sound speed [m/s]
%        'SoundSpeedRef_water'  - the sound speed for only water [m/s]
%        'SoundSpeedRef_object' - the mean sound speed of water encompassing
%                                 the object during the rotations [m/s]
%        'minimum_distance'     - the minimum distance between
%                                  emitter-receiver for TOF picking [m]
%        'time_index_end'       - the last time index of the time array for including
%                                 in the TOF picking
%      'choose_last_us_position' - A parameter determining how to choose the last
%                           US position, which can 'early' (chosen by the user),
%                           or automatically set all the angular US positions - 1 ('all'). 
%                           The choice for the latter is based on the fact that 
%                           the last position is the same as the first position,
%                           so the last position is removed. (Default = 'all')
%        'binaries_emitter_receiver' - the method for choosing the
%                                 emitter-receiver pair for TOF picking
%        'open_angle'           - an open angle for a cone with central
%                                 axis a line segment from emitter to the
%                                 centre of the bowl. The TOFs are calculated
%                                 for receivers inside the cone
%        'Cutoff_freq'          - the cut-off frequncy of the filter [Hz]
%                                 This can be set nan (no filtering), has one
%                                 component (low-pass filtering), or  two
%                                 components (band-pass filtering)
%        'absorption_frequencies' - the frequencies at which the absorption
%                                   is reconstructed (Hz)
%        'save_data'              - the Boolean controlling whether the data is
%                                 saved or not. 
%                                 
% OUTPUTS:
%        tof_dif               - the time discrepancy of the time-of-flights
%                                for the object and only water[s]
%        tof_het               - a matrix of size num_receiver x num_emitter containing
%                                the time-of-flights [s] for the object data
%        tof_hom               - a matrix of size num_receiver x num_emitter containing
%                                the time-of-flights [s] for the water data
%        amplitud_het          - a struct array each containing a matrix of size
%                                num_receiver x num_emitter containing the
%                                amplitudes [a.u.] for the object data at
%                                the chosen freqiencies.
%       amplitude_hom          - a struct array each containing a matrix of size
%                                num_receiver x num_emitter containing the
%                                amplitudes [a.u.] for the water data at
%                                the chosen freqiencies.
%        elapsed_time          - the computational time [s] for running this m-file
%                                function
%
%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 05.08.2020
%       last update     - 12.12.2020
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2020 Ashkan Javaherian

% optinal inputs
para = [];
para.nWorkerPool = 16;
para.Method = 'Modified_AIC';
para.SoundSpeedRangesDiff = [-50, 50];
para.SoundSpeedRef_water = 1500;
para.SoundSpeedRef_object = 1500;
para.minimum_distance = 0.16;
para.time_index_end = 6000;
para.choose_last_us_position = 'all'; 
para.binaries_emitter_receiver = 'open_angle';
para.Cutoff_freq = nan;
para.open_angle = pi/4;
para.absorption_frequencies = nan;
para.save_data = true;



% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


% start the time
time_start = tic;

% Get the optional inputs for picking the first-arrival time [s] of the measured
% time traces
tof_args = {'nWorkerPool', para.nWorkerPool, 'Method', para.Method,...
    'SoundSpeedRanges', 2/5 * para.SoundSpeedRangesDiff + para.SoundSpeedRef_water,...
    'SoundSpeedRef', para.SoundSpeedRef_water, 'Threshold', 0,...
    'minimum_distance', para.minimum_distance,...
    'Plot', false', 'binaries_emitter_receiver', para.binaries_emitter_receiver,...
    'open_angle', para.open_angle, 'Cutoff_freq', para.Cutoff_freq};


% Get the optional input for computing amplitudes [a.u.] of the measured
% time traces
amp_args = {'nWorkerPool', para.nWorkerPool, ...
    'SoundSpeedRanges', [4 * para.SoundSpeedRangesDiff(1), para.SoundSpeedRangesDiff(2)] + para.SoundSpeedRef_water,...
    'SoundSpeedRef', para.SoundSpeedRef_water, 'minimum_distance', para.minimum_distance,...
    'Plot', false', 'binaries_emitter_receiver', para.binaries_emitter_receiver,...
    'open_angle', para.open_angle};


% Get the number of receivers
num_receiver = size(receiver.positions{1}, 2);

% Get the number of emitters
num_emitter = size(emitter.positions, 2);

% Get the number of angular US positions
num_us_position = max(rotation_indices);

% Get the number of angular positions per data segment
num_us_position_per_segment = floor(num_us_position/num_us_segment);

% the number of remaining emitters
num_emitter_rem = rem(num_us_position, num_us_segment);

if num_us_position_per_segment == 0
    
    % if the number of US positions are less than the chosen number of
    % data segments
    num_us_segment = num_us_position;
end

% Calculate the TOFs for only-water data for the first num_receiver
% excitations
% Note that num_receiver, the number of receivers, also equals the number of
% transducers, ie. the number of emitters for one full rotation over all the
% transducers.
emitter_ind = emitter;


% extract the position of emitters, receivers and the rotations for the
% first num_receiver exciations
emitter_ind.positions = emitter.positions(:, 1:num_receiver);
receiver_ind.positions = receiver.positions(1:rotation_indices(num_receiver));
rotation_indices_ind = rotation_indices(1: num_receiver);

% Calculate the time-of-flight sinogram for only-water (reference) data for
% the first angular position
disp('computing time-of-flight sinogram for the only-water data:');
[~, ~, tof_hom_first] = timeOfFlightPicking([], data_water(:, 1:para.time_index_end,:),...
    emitter_ind, receiver_ind, time_array(1:para.time_index_end), rotation_indices_ind, tof_args{:});

if all(isfinite(para.absorption_frequencies))
    disp('computing amplitudes for the only-water data:');
    
    % compute the amplitudes for reference (only water) data for the first angular position
    [~, amplitude_hom_first] = calcAmplitude([], data_water(:, 1:para.time_index_end,:),...
        emitter_ind, receiver_ind, time_array(1:para.time_index_end), rotation_indices_ind,...
        para.absorption_frequencies, amp_args{:});
    
    % convert the amplitudes for the only-water data from matrix to cell array
    % each containing the amplitudes for all emitter-receiver pairs at the chosen
    % frequencies. Note that for the only-water data, the number of emitters and receivers are
    % both equal the number of transducers.
    amplitude_hom_first = mat2cell(amplitude_hom_first,...
        num_receiver * num_receiver, ones(1, length(para.absorption_frequencies)));
    
    % Expand the only-water amplitudes for making them consistent with the object-in-water
    % amplitudes for all the rotations
    amplitude_hom = cell(size(amplitude_hom_first));
    
    for ind_freq = 1:length(para.absorption_frequencies)
        
        % reshape the vector for amplitudes for all emitter-receiver pairs at
        % the specific frequency to a matrix of size num_receiver x
        % num_receiver
        amplitude_hom_first{ind_freq} = ...
            reshape(amplitude_hom_first{ind_freq}, [num_receiver, num_receiver]);
        
        % Expand the columns from num_receiver to num_emitter by iterating
        % over every num_receiver excitations.
        
        % The resulting matrix will be of size num_receiver x num_emitter,
        % and contains the only-water amplitudes for all emitter-receiver pairs
        % at the specific frequency
        amplitude_hom{ind_freq} =...
            [repmat(amplitude_hom_first{ind_freq}, [1, floor(num_emitter/num_receiver)]),...
            amplitude_hom_first{ind_freq}(:, 1: rem(num_emitter, num_receiver))];
    end
else
    
    % allocate an empty varibale for the only-water amplitudes, if
    % para.absorption_frequencies is set nan.
    amplitude_hom = [];
    
end

% Expand the only-water TOFs for making them consistent with the TOFs
% for the object-in-water data for all rotations
tof_hom = [repmat(tof_hom_first, [1, floor(num_emitter/num_receiver)]),...
    tof_hom_first(:, 1: rem(num_emitter, num_receiver))];


% get the nan binary in the US shot index
[notnans_binary, ~] = ismember(us_shot_index, emit_shot_index);


% get the rotation indices for the US shots from the rotation indices for the emit shots
rotation_indices_us_shots = zeros(size(us_shot_index));

% rotation_indices_us_shots(num_receiver + 1 :num_receiver + 1: end) = nan;
rotation_indices_us_shots(~notnans_binary) = nan; 
rotation_indices_us_shots(~isnan(rotation_indices_us_shots)) = rotation_indices;


% Get the optional inputs for time-of-flight picking function, which are not
% actually optional and are set in the main script
tof_args = {'nWorkerPool', para.nWorkerPool, 'Method', para.Method,...
    'SoundSpeedRanges', para.SoundSpeedRangesDiff + para.SoundSpeedRef_object,...
    'SoundSpeedRef', para.SoundSpeedRef_object, 'Threshold', 0,...
    'minimum_distance', para.minimum_distance,...
    'Plot', false, 'binaries_emitter_receiver', para.binaries_emitter_receiver,...
    'open_angle', para.open_angle, 'Cutoff_freq', para.Cutoff_freq};

% Get the optional input for computing amplitudes [a.u.] of the measured
% time traces
amp_args = {'nWorkerPool', para.nWorkerPool, ...
    'SoundSpeedRanges',[4 * para.SoundSpeedRangesDiff(1), para.SoundSpeedRangesDiff(2)] + para.SoundSpeedRef_object,...
    'SoundSpeedRef', para.SoundSpeedRef_object, 'minimum_distance', para.minimum_distance,...
    'Plot', false', 'binaries_emitter_receiver', para.binaries_emitter_receiver,...
    'open_angle', para.open_angle};


% Allocate a zero matrix for the TOFs for the object-in-water data
tof_het = zeros(num_receiver, num_emitter);

if all(isfinite(para.absorption_frequencies))
    
    % Allocate a struct array containing zero matrices for the amplitudes
    % for the object-in-water data at the chosen frequencies
    amplitude_het = cell(length(para.absorption_frequencies), 1);
    
    for ind_freq = 1:length(para.absorption_frequencies)
        amplitude_het{ind_freq} = zeros(num_receiver, num_emitter);
    end
else
    
    % allocate an empty variable for the object-in-water amplitudes,
    % para.absorption_frequencies is set nan.
    amplitude_het = [];
end

% get an empty variable
if para.save_data && strcmp(para.choose_last_us_position, 'early') 
data_breast = [];
end

for ind_segment = 1 : num_us_segment
    
    if ind_segment < num_emitter_rem + 1
        
        % Get the binary indices in the 'emit_shot_index' for the current data
        % segment
        [emitter_binaries_segment, ~] = ismember(rotation_indices, (ind_segment - 1) * ...
            (num_us_position_per_segment + 1) + 1: ind_segment * (num_us_position_per_segment + 1));
        
        % get the binary indices in the 'us_shot_index' for the current data
        % segment
        [us_binaries_segment, ~] = ismember(rotation_indices_us_shots, (ind_segment - 1) * ...
            (num_us_position_per_segment + 1) + 1: ind_segment * (num_us_position_per_segment + 1));
        
    else
        
        % Get the binary indices in the 'emit_shot_index' for the current data
        % segment
        [emitter_binaries_segment, ~] = ismember(rotation_indices, num_emitter_rem * (num_us_position_per_segment + 1) + ...
            (ind_segment - num_emitter_rem - 1) * num_us_position_per_segment + 1: (ind_segment - num_emitter_rem) * num_us_position_per_segment);
        
        % Get the binary indices in the 'us_shot_index' for the current data
        % segment
        [us_binaries_segment, ~] = ismember(rotation_indices_us_shots, num_emitter_rem * (num_us_position_per_segment + 1) + ...
            (ind_segment - num_emitter_rem - 1) * num_us_position_per_segment + 1: (ind_segment - num_emitter_rem) * num_us_position_per_segment);
        
    end
    
    % get the number of emitters (excitations) for the current data segment
    num_emitter_ind = nnz(emitter_binaries_segment);
    
    % Get the position of emitters associated with the current data segment
    emitter_ind.positions = emitter.positions(:, emitter_binaries_segment);
    
    % Get the position of the receivers for the current data segment
    if ind_segment < num_emitter_rem + 1
        receiver_ind.positions = receiver.positions( (ind_segment - 1) * ...
            (num_us_position_per_segment + 1) + 1: ind_segment * (num_us_position_per_segment + 1) );
    else
        receiver_ind.positions = receiver.positions( num_emitter_rem * (num_us_position_per_segment + 1) + ...
            (ind_segment - num_emitter_rem - 1) * num_us_position_per_segment + 1: (ind_segment - num_emitter_rem) * num_us_position_per_segment  );
    end
    
    % Get the rotation indices for the current data
    us_indices_segment_first = find(us_binaries_segment, 1, 'first');
    us_indices_segment_last = find(us_binaries_segment, 1, 'last');
    
    % Load object data for the object data indexed by ind_segment
    disp(['loading the data from object:' num2str(ind_segment)]);
    data_object = h5read([data_paths.main_directory,...
        data_paths.directory, data_paths.data_name '.hdf5'], '/us_data',...
        [1, 1, us_indices_segment_first],...
        [para.time_index_end, inf, us_indices_segment_last-(us_indices_segment_first-1)]);
    
    % Remove the nan index from the data, convert it to a single array,
    % and exchange the first and second dimensions of data
    % for making it consistent with the ray tracing codes
    data_object = permute(single(data_object(:,:,us_binaries_segment...
        (us_indices_segment_first:us_indices_segment_last) )), [2, 1, 3]);
    
    % Get the 1-based rotation indices for the current data segment
    rotation_indices_ind = rotation_indices(emitter_binaries_segment)...
        - min(rotation_indices(emitter_binaries_segment & rotation_indices > 0)) + 1;
    
    if para.save_data && strcmp(para.choose_last_us_position, 'early') 
    data_breast = cat(3, data_breast, data_object);
    end
      
    % Calculate the time-of-flights sinogram for the object data for
    % the current data segment
    disp(['computing time-of-flight sinogram for the object-in-water data' num2str(ind_segment)]);
    [~, tof_het(:, emitter_binaries_segment) , ~] = timeOfFlightPicking(data_object, [], emitter_ind,...
        receiver_ind, time_array(1:para.time_index_end), rotation_indices_ind, tof_args{:});
    
    if all(isfinite(para.absorption_frequencies))
        
        disp(['computing amplitudes for the object-in-water data:' num2str(ind_segment)]);
        % Compute the amplitudes for the object-in-water data for the
        % current data segment
        [amplitude_het_segment, ~] = calcAmplitude(data_object, [], ...
            emitter_ind, receiver_ind, time_array(1:para.time_index_end),...
            rotation_indices_ind, para.absorption_frequencies, amp_args{:});
        
        % convert the amplitudes from matrix to a cell array each containing the
        % amplitudes for all emitter-receiver pairs at the chosen frequencies
        amplitude_het_segment = mat2cell(amplitude_het_segment,...
            num_receiver * num_emitter_ind, ones(1, length(para.absorption_frequencies)));
        
        
        for ind_freq = 1:length(para.absorption_frequencies)
            
            % Put the object-in-water amplitudes for the current data
            % segment in the matrix for the object-in-water amplitudes for
            % all the emitters
            amplitude_het{ind_freq}(:, emitter_binaries_segment) = ...
                reshape(amplitude_het_segment{ind_freq}, [num_receiver, num_emitter_ind]);
        end
        
    end
    
end

if para.save_data && strcmp(para.choose_last_us_position, 'early') 

if ~exist([user_data_path 'time_series_early.h5' ])
    
    % create the HDF5 paths
    h5create([user_data_path, 'time_series_early.h5'], '/data_breast', size(data_breast));
end

% write the HDF5 file on the defined user data path
h5write([user_data_path, 'time_series_early.h5'], '/data_breast', data_breast);

end

% calculate the elapsed time
elapsed_time = toc(time_start);

end