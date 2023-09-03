function [data_first_position, elapsed_time] = loadPammothFirstPositionData(...
    rotation_indices, us_shot_index, num_us_segment, num_receiver, data_paths, varargin)
%LOADPAMMOTHFIRSTPOSITIONDATA loads the data for the first position 
%
%
% DESCRIPTION:
%       loadPammothFirstPositionData loads the data for the first position
%       and selects the position of the corresponding emitter-receiver pairs 
%
% USAGE:
%
%
% INPUTS:

%      emitter.pulse               - the excitation pulse
%                                    pulses for different emitters are different
%      emitter.pulse_duration      - a scalar indicating the the approximate
%                                    time duration of the excitation pulses,
%                                    which will be used for calculation of the
%                                    first-arrival of the signals
%      emitter.pulse_first_arrival - the first arrival of the excitation
%                                    pulse (or the ultrasound shot time),
%                                     which can be negative [sec]
%                                    if negative, the time series has been started
%                                    after pulse excitation
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
%                                    [sec]
%     data_water                   - the data from water
%     num_us_segment               - the number of data segments for loading
%                                  - the number of data segments for time-of-flight picking
%                                    (for saving memeory)
%     data_paths                   - path for loading data
%
% OPTIONAL INPUTS:
%     'nWorkerPool'                - the number of used workers


% OUTPUTS:
%     data_first_position    - the object-in-water data for the first position
%
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
para.time_index_end = 6000;




% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


% start the time
time_start = tic;

% Get the number of angular US positions
num_us_position = max(rotation_indices);

% Get the number of angular positions per data segment
num_us_position_per_segment = floor(num_us_position/num_us_segment);

% the number of remaining emitters
num_emitter_rem = rem(num_us_position, num_us_segment);


if num_us_position_per_segment == 0
    
    % if the number of US positions are less than the chosen number of
    % data segments
    num_us_segment = num_emitter_rem;
end



% Convert the rotation indices for the emit shots to the rotation
% indices for the emit shots
rotation_indices_us_shots = zeros(size(us_shot_index));
rotation_indices_us_shots(num_receiver + 1 :num_receiver + 1: end) = nan;
rotation_indices_us_shots(~isnan(rotation_indices_us_shots)) = rotation_indices;


% allocate an empty variable to the data for the first position
data_first_position = [];


for ind_segment = 1 : num_us_segment
    
    if ind_segment < num_emitter_rem + 1
        
        % the binary indices in the 'us_shot_index' for the current data
        % segment
        [us_binaries_segment, ~] = ismember(rotation_indices_us_shots, (ind_segment - 1) * ...
            (num_us_position_per_segment + 1) + 1: ind_segment * (num_us_position_per_segment + 1));
        
    else
        
       
        % the binary indices in the 'us_shot_index' for the current data
        % segment
        [us_binaries_segment, ~] = ismember(rotation_indices_us_shots, num_emitter_rem * (num_us_position_per_segment + 1) + ...
            (ind_segment - num_emitter_rem - 1) * num_us_position_per_segment + 1: (ind_segment - num_emitter_rem) * num_us_position_per_segment);
        
    end
    
    
 
    
    % Get the rotation indices for the current data
    us_indices_segment_first = find(us_binaries_segment , 1, 'first');
    us_indices_segment_last = find(us_binaries_segment , 1, 'last');
    % Load object data for the object data indexed by ind_segment
    disp(['loading the data from object:' num2str(ind_segment)]);
    data_object = h5read([data_paths.main_directory,...
        data_paths.directory, data_paths.data_name '.hdf5'], '/us_data',...
        [1, 1, us_indices_segment_first],...
        [para.time_index_end, inf, us_indices_segment_last-(us_indices_segment_first-1)]);
    
    % Remove nan index from the data, convert it to a single array,
    % and exchange the first and second dimensions of data
    % in order to make it consistent with the ray tracing codes
    data_object = permute(single(data_object(:,:,us_binaries_segment...
        (us_indices_segment_first:us_indices_segment_last) )), [2, 1, 3]);
    
    % add the data for the current data segment to the data for the first
    % position
    data_first_position = cat(3, data_first_position, data_object);
    
    
end

% calculate the elapsed time
elapsed_time = toc(time_start);

end