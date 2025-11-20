function [amplitude_het, amplitude_hom] = calcAmplitude(data_het,...
    data_hom, emitter, receiver, time_array, rotation_indices, frequencies, varargin)
%calcAmplitude calculates the amplitudes, or the logarithmic relative change in the amplitudes between the
%data measured from the object inside water and only water
%
%
% DESCRIPTION:
%       calcAmplitude calculates the amplitudes or logarithmic relative change [dB] in the amplitudes between the
%       the data measured from the object inside water and only water at specific frequencies.
%
% USAGE:
%
%
% INPUTS:
%      data_het                   -  a matrix of size num_receiver x num_time x num_emitter,
%                                    and represents the data set that is measured
%                                    from the object inside the water (main data)
%      data_hom                    - a matrix of size num_receiver x num_time x num_emitter,
%                                    and represents the data set that is measured
%                                    from a homogeneous liquid like water (reference data)
%      emitter                     - a struct that defines the properties of the
%                                    excitation as follows: This includes the
%                                    fields 'positions', 'pulse' and 'pulse_duration'
%      emitter.positions           - 2/3 x N array of cartesian points containing
%                                    the centers of the emitter objects
%      emitter.pulse               - a vector of size 1 x num_te
%                                    with num_te <= num_time.
%                                    For real data, this may be a matrix of size
%                                    num_emitter x num_te, because the excitation
%                                    pulses for different emitters are different
%      emitter.pulse_duration      - a scalar indicating the approximate
%                                    time duration of the main lobe of the
%                                    excitation pulses, which will be used
%                                    for calculation of the first-arrival of
%                                    the signals
%      emitter.shot_time           - the first arrival of the excitation
%                                    pulse, which can be negative [sec]
%                                    if negative, the time after pulse excitation
%                                    the data measurement has been started [sec]
%      receiver.positions          - 2/3 x N array of cartesian points containing
%                                    the centers of the receiver objects
%                                    For the cases at which the position of
%                                    the receivers is changed with position of
%                                    the emitter, the array will be a cell
%                                    array with size the number of rotation
%                                    angles with each cell a 2/3 x N array of
%                                    cartesian points for the position o
%                                    transducers for each angle
%       time_array                 - 1 x N_t vector of measurement times
%                                    [sec]
%       'frequencies'              - the specific frequencies at which the
%                                   attenuation is computed
% OPTIONAL INPUTS:
%        'nWorkerPool'          - the number of used workers
%        'SoundSpeedRanges'     - a vector of size 1 x 2 containing a
%                                 minimal and maximal for the sound speed
%                                 [m/sec]
%        'binaries_emitter_receiver' - the method for choosing the
%                                 emitter-receiver pairs. This can be 'open_angle'
%                                 (only the directivity of the emitters are included)
%                                 or 'distances' (emitters and receivers with
%                                 distances larger than a specific threshold
%                                 are included, so the directivity of both emitters
%                                 and receivers are included.)
%        'minimum_distance'      - the mimimum distance for emitter-receiver
%                                  pair being included
%        'open_angle'            - the maximum angle between two vecors, the
%                                 first of which connecting the emitter to
%                                 the centre, and the second connecting the
%                                 the emitter to the receiver. For each
%                                 emitter, the included receivers will be
%                                 inside a cone with axis the vector
%                                 connecting the emitter to all receivers.
%        'Plot'                 - boolean for plotting the calculated time
%                                 differences as a sinogram
% OUTPUTS:
%        amplitude              - the logarithmic relative change in the amplitudes
%                                 at specific frequencies
%        amplitude_het          - the amplitudes at the chosen frequencies for the data
%                                 measured from the object inside water
%        amplitude_hom          - the amplitudes at the chosen frequencies for the data
%                                 measured from the water
%        elapsed_time           - the computational time of running this m-file
%                                 function
%
%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian

para.nWorkerPool = 16;
para.SoundSpeedRanges = [1400, 1600];
para.SoundSpeedRef = 1500;
para.minimum_distance = 0.08;
para.binaries_emitter_receiver = 'open_angle';
para.Order = 8;
if strcmp(para.binaries_emitter_receiver, 'open_angle')
    para.open_angle = pi/4;
end
para.minimum_distance = 0.08;

% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end



% the number of emitters
num_emitter = size(emitter.positions, 2);

% the number of receivers
if  iscell(receiver.positions)
    num_receiver = size(receiver.positions{1}, 2);
else
    num_receiver = size(receiver.positions, 2);
end


% calculate the time spacing of the measurement
dt = time_array(2) - time_array(1);

% get a Tukey window
tukey = getWin(length(time_array), 'Tukey', 'Param', 0.25).';

% set the optional inpus for FFTs
f_args = {'FFTLength', 2 * length(time_array)};

% calculate the distances between each pair of emitters and receivers
distances = calculateDistanceEmitterReceiver(emitter.positions,...
    receiver.positions, rotation_indices);

switch  para.binaries_emitter_receiver
    case 'distances'
        binaries_emitter_receiver = mat2cell(distances > para.minimum_distance,...
            num_receiver, ones(1, num_emitter) );
    case 'open_angle'
        binaries_emitter_receiver = mat2cell(calcAngleEmitterReceiver(emitter.positions,...
            receiver.positions, rotation_indices) < para.open_angle, ...
            num_receiver, ones(1, num_emitter) );
end


if isfield(emitter, 'pulse')
    
    % if the excitation pulse is given, calcuate the first_arrival of the
    % excitation pulse (This must be modified if we have separate
    % exciataion pulse for each emitter)
    excitation_pulse = zeros(1, length(time_array));
    excitation_pulse(1:length(emitter.pulse)) = emitter.pulse;
    M = para.Order;
    time_excit = araic(excitation_pulse, dt, time_array, emitter.pulse_duration,...
        10, 1, 0.75, M, time_array(M + 1), time_array(end - M));
    
    % display the excitation time
    disp(['The excitation time is:' num2str(time_excit) '[s]'])
    
elseif isfield(emitter, 'shot_time')
    
    % if the excitation pulse is not given, but the first arrival of
    % the excitation pulse is given
    time_excit = emitter.shot_time;
    
end

% display the excitation time
disp(['The time of excitation in the time array is:' num2str(1e6*time_excit,'%3.5f') ':'  'microseconds']);

% calculate the minimum and maximum time (in time_array) for the occurrence of the first
% arrival in the main data
tof_min = mat2cell(distances./para.SoundSpeedRanges(2) + time_excit,...
    num_receiver, ones(1, num_emitter) );
tof_max = mat2cell(distances./para.SoundSpeedRanges(1) + time_excit,...
    num_receiver, ones(1, num_emitter) );

if ~isempty(data_hom)
    
    % convert the reference data matrix to a cell array
    data_hom =  mat2cell(data_hom, num_receiver, size(data_hom, 2),...
        ones(num_emitter, 1) );
    
    % allocate a cell array for tofs for the reference data
    amplitude_hom_cell = cell(num_emitter, 1);
    
    % calcuation of first-arrivals for the reference data
   parfor (ind_emitter = 1:num_emitter, para.nWorkerPool)
      %   for ind_emitter = 1:num_emitter % (for test)
        
        % get the reference data for the current emitter
        data_hom_emitter = bsxfun(@times, tukey, removeDataDc(data_hom{ind_emitter}) );
        
        % get the binaries for the emitter
        binaries_emitter = binaries_emitter_receiver{ind_emitter};
        
        % get the minimum expected TOFs for the emitter
        tof_min_emitter = tof_min{ind_emitter};
        
        % get the maximum expected TOFs for the emitter
        tof_max_emitter = tof_max{ind_emitter};
        
        % allocate zero matrix
        amplitude_hom = zeros(num_receiver, length(frequencies));
          
        if any(frequencies) > 0
            
            
            for ind_receiver = 1:num_receiver
                %  disp(ind_receiver);
                
                if binaries_emitter(ind_receiver)
                    
                    
                    time_binary = time_array >= tof_min_emitter(ind_receiver) & ...
                        time_array <= tof_max_emitter(ind_receiver);
                    
                    if any(time_binary)
                        
                    % get the FFT of the signals
                    [f, amplitude_hom_receiver, ~] = spect(data_hom_emitter(ind_receiver, time_binary).', 1/dt, f_args{:});
                    
                    % get the indices for the interested frequencies
                    [~, f_indices] = min(abs(frequencies - f));
                    
                    % get the amplitudes at the chosen frequency indices
                    amplitude_hom(ind_receiver, :) = amplitude_hom_receiver(f_indices).';
                    
                    else
                        
                    amplitude_hom(ind_receiver, :) = 1 ;
                    
                    end
                    
                else
                    
                    % get the maximum absolute amplitudes
                    amplitude_hom(ind_receiver, :) = 1;
                    
                end
                
            end
            
        else
            
            for ind_receiver = 1:num_receiver
                %  disp(ind_receiver);
                
                if binaries_emitter(ind_receiver)
                    
                   time_binary = time_array >= tof_min_emitter(ind_receiver) & ...
                        time_array <= tof_max_emitter(ind_receiver);
                    
                    if any(time_binary)
                        
                    % get the maximum absolute amplitudes
                    amplitude_hom(ind_receiver) = max(abs(data_hom_emitter(ind_receiver,...
                        time_binary) ) );
                    else
                        
                    amplitude_hom(ind_receiver) = 1;
                    end
                    
                else
                    
                    amplitude_hom(ind_receiver) = 1;
                    
                end
                
            end
        end
        

        amplitude_hom_cell{ind_emitter} = amplitude_hom;
        
        disp(['Computing the amplitudes for the only-water data, Elapsed Time Percentage:'...
            num2str(ind_emitter/num_emitter * 100)]);
        
        
        
    end
    
    % clear the reference data for saving memory
    clear data_hom
    
else
    
    amplitude_hom_cell = [];
    
end



if ~isempty(data_het)
    
    % convert the reference data matrix to a cell array
    data_het = mat2cell(data_het, num_receiver, size(data_het, 2),...
        ones(num_emitter, 1) );
    
    % allocate a cell array for the tofs for the reference data
    amplitude_het_cell = cell(num_emitter, 1);
    
    % calcuation of first-arrivals for the reference data
    parfor (ind_emitter = 1:num_emitter, para.nWorkerPool)
        % for ind_emitter = 1:num_emitter % (for test)
        
        % get the reference data for the current emitter
        data_het_emitter = bsxfun(@times, tukey, removeDataDc(data_het{ind_emitter}) );
        
        % get the binaries for the emitter
        binaries_emitter = binaries_emitter_receiver{ind_emitter};
        
        % get the minimum expected TOFs for the emitter
        tof_min_emitter = tof_min{ind_emitter};
        
        % get the maximum expected TOFs for the emitter
        tof_max_emitter = tof_max{ind_emitter};
        
        % allocate zero matrix
        amplitude_het = zeros(num_receiver, length(frequencies));
          
        if any(frequencies) > 0
            
            
            for ind_receiver = 1:num_receiver
                %  disp(ind_receiver);
                
                if binaries_emitter(ind_receiver)
                    
                    
                    time_binary = time_array >= tof_min_emitter(ind_receiver) & ...
                        time_array <= tof_max_emitter(ind_receiver);
                    
                    if any(time_binary)
                    % get the FFT of the signals
                    [f, amplitude_het_receiver, ~] = spect(data_het_emitter(ind_receiver, time_array >= tof_min_emitter(ind_receiver) & ...
                        time_array <= tof_max_emitter(ind_receiver)).', 1/dt, f_args{:});
                    
                    % get the indices for the interested frequencies
                    [~, f_indices] = min(abs(frequencies - f));
                    
                    % get the amplitudes at the chosen frequency indices
                    amplitude_het(ind_receiver, :) = amplitude_het_receiver(f_indices).';
                    
                    
                    else
                        
                      % get the maximum absolute amplitudes
                    amplitude_het(ind_receiver, :) = 1;  
                        
                    end
                    
                else
                    
                    % get the maximum absolute amplitudes
                    amplitude_het(ind_receiver, :) = 1;
                    
                end
                
            end
            
        else
            
            for ind_receiver = 1:num_receiver
                %  disp(ind_receiver);
                
                if binaries_emitter(ind_receiver)
                    
                     time_binary = time_array >= tof_min_emitter(ind_receiver) & ...
                        time_array <= tof_max_emitter(ind_receiver);
                    
                    if any(time_binary)
                    
                    % get the maximum absolute amplitudes
                    amplitude_het(ind_receiver) = max(abs(data_het_emitter(ind_receiver,...
                        time_array >= tof_min_emitter(ind_receiver) & ...
                        time_array <= tof_max_emitter(ind_receiver)) ) );
                    
                    else
                        
                        amplitude_het(ind_receiver) = 1;
                        
                    end
                    
                else
                    
                    amplitude_het(ind_receiver) = 1;
                    
                end
                
            end
        end
        

        amplitude_het_cell{ind_emitter} = amplitude_het;
        
        disp(['Computing the amplitudes for the only-water data, Elapsed Time Percentage:'...
            num2str(ind_emitter/num_emitter * 100)]);
        
    end
    
    % clear the reference data for saving memory
    clear data_het
    
else
    
    amplitude_het_cell = [];
    
end







% convert the cells to a stacked vector for all emitters-receivers
if ~isempty(amplitude_hom_cell)
    amplitude_hom = cell2mat(amplitude_hom_cell);
else
    amplitude_hom = [];
end

if ~isempty(amplitude_het_cell)
    amplitude_het = cell2mat(amplitude_het_cell);
else
    amplitude_het = [];
end



end