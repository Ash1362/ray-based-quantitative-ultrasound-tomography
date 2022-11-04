function [tof, tof_het, tof_hom] = timeOfFlightPicking(data_het, data_hom,...
    emitter, receiver, time_array, rotation_indices, varargin)
%TIMEOFFLIGHTPICKING calculates the discrepancy of the first-arrival times
%between two sets of data
%
%
% DESCRIPTION:
%       timeOfFlightPicking calculates the discrepancy of the first-arrival times between
%       two data sets, the first of which is measured from the phantom
%       inside water, and the second is measured from only water.
%       This will be used for calculation of discrepancy of time-of-flight
%       of the (first-arrival) acoustic waves from an emitter to a receiver
%       between the phantom inside water and only water.
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
%      emitter.shot_time           - the shot time [s] of the excitation
%                                    pulse, which can be negative.
%                                    Negative value indicates that the data 
%                                    measurement has been started after
%                                    excitation time.
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
%
% OPTIONAL INPUTS:
%        'nWorkerPool'          - the number of used workers
%        'Method'               - the method for calculation of the first-arrival
%                                 of the signals. This must be set
%                                 'Modified_AIC'. Other methods
%                                 'Short_Time_Average/Long_Time_Average',
%                                 'Modified_Energy_Ratio',
%                                 'Modified_Coppens'
%                                 'Modified_Short_Term_Average/Long_Term_Average'
%                                  are for the evaluation purposes, and are
%                                  used only by the code developer.
%        'SoundSpeedRanges'     - a vector of size 1 x 2 containing a
%                                 minimal and maximal for the sound speed
%                                 [m/sec]
%        'DistanceThreshold'    - For distances between emitter and receiver
%                                 smaller than this threshold, the
%                                 time discrepancy of the first-arrivals is
%                                 not calculated, and set zero [m]
%      'binaries_emitter_receiver' - the method for choosing the
%                               emitter-receiver pairs. This can be
%                               'open_angle', i.e., the angle between
%                               a geometrical vactor connecting the
%                               emitter to the centre of the detection suface (ring),
%                               and another vector connecting the emitter to
%                               receivers. The emitter-receiver pairs with open angles
%                               larger than a threshold are excluded. This can also be set
%                               'distances', i.e., emitters and receivers with
%                               distances smaller than a specific threshold
%                               are excluded. Bothe approaches are
%                               equivalent.
%       'minimum_distance'     - the mimimum distance for emitter-receiver
%                                pairs
%        'open_angle'          - the maximum angle between two vecors, one
%                                connecting the emitter to the centre of the
%                                detection surface (ring), and another connecting
%                                the emitter to the receivers. For each emitter, 
%                                the included receivers will be inside a cone
%                                with axis the vector connecting the emitter to
%                                the centre of the detection surface.
%        'Cutoff_freq'           - the cut-off frequency of the filter, if
%                                 one componet (low-pass filter), if two component
%                                 (band-pass filter)
%         'Order'                - the order of the regression equation, if 'AR_AIC'
%                                  approach is used.
%         'Beta'                 - a scalar used for smoothing the
%                                  short-time-average to long-time-average
%                                 (this is only used for 'Modified
%                                 Coppens' method)
%        'Plot'                 - boolean for plotting the calculated time
%                                 differences as a sinogram
% OUTPUTS:
%        tof                    - the time discrepancy of the first-arrivals
%                                 (time-of-flights) [sec]
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
% Copyright (c) 2020 Ashkan Javaherian

% optinal inputs
para = [];
para.nWorkerPool = 8;
para.Method = 'Modified_AIC';
para.SoundSpeedRanges = [1400, 1600];
para.SoundSpeedRef = 1500;
para.Beta = 1e-1;
para.Order = 8;
para.Plot = true;
para.binaries_emitter_receiver = 'open_angle';

if strcmp(para.binaries_emitter_receiver, 'open_angle')
    para.open_angle = pi/4;
end
para.minimum_distance = 0.08;
para.Cutoff_freq = nan; % [0.5, 2]*1e6;


% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


% the method used for calculation of the first arrival
picking_method = para.Method;

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

    
% get the transition width
trans_width_MHz = 10e6;

% get the frequency
fs = 1/dt;


trans_width_prop = trans_width_MHz / fs;


if length(para.Cutoff_freq) > 1
    error('cut-off frequency must be a scalar value or nan.')
    
end


% calculate the distances between each pair of emitters and receivers
distances = calculateDistanceEmitterReceiver(emitter.positions,...
    receiver.positions, rotation_indices);

switch  para.binaries_emitter_receiver
    case 'distances'
        binaries_emitter_receiver  = distances > para.minimum_distance;
    case 'open_angle'
        binaries_emitter_receiver = calcAngleEmitterReceiver(emitter.positions,...
            receiver.positions, rotation_indices) < para.open_angle;
end


% an estimate for the pulse duration
T = emitter.pulse_duration;

if isfield(emitter, 'pulse')
    
    % if the excitation pulse is given, calcuate the first_arrival of the
    % excitation pulse (This must be modified if we have separate
    % exciataion pulse for each emitter)
    excitation_pulse = zeros(1, length(time_array));
    excitation_pulse(1:length(emitter.pulse)) = emitter.pulse;
    
    % M = para.Order;
    excitation_time_max = dt * find(abs(excitation_pulse)/max(abs(excitation_pulse))>0.2, 1, 'first');
    tof_args ={'Method', picking_method, 'Length_moving_windows',...
        [1.0, nan, 0.25, nan], 'Threshold', 0};
    time_excit = timeOfFlightPickingEachSignal(excitation_pulse,...
        time_array, T, [dt, excitation_time_max], tof_args{:} );
    
    % the line written below are deprecated.
    %time_excit = araic(excitation_pulse, dt, time_array, T, 10, 1, 0.75, M,...
    %  time_array(M + 1), time_array(end - M));
    
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
TOF_min = (distances./para.SoundSpeedRanges(2)) + time_excit;
TOF_max = (distances./para.SoundSpeedRanges(1)) + time_excit;

% calculate the minimum and maximum time (in time_array) for the occurrence of the first
% arrival in the reference data
TOF_min_ref  = (distances./para.SoundSpeedRanges(2)) + time_excit;
TOF_max_ref  = (distances./para.SoundSpeedRanges(1)) + time_excit;


% Apply a band-pass filter on the data sets in the frequency domain
% get a Tukey window
tukey = getWin(length(time_array), 'Tukey', 'Param', 0.25).';

% define the optional inputs
% define a vector indicating the length of the
% moving windows:
% 1 - backward window, 2 - long window,
% 3 - short window, 4 - EPS window


% Table of prefererred parameters for different methods using
% the excitation pulse 'Pammoth_1'
% ---------------------------------------------------
%  Method:     Back    Long  Short   EPS    Threshold
%----------------------------------------------------
%  STA/LTA      1.0     2.5   0.5    nan      0.5
%  mSTA/LTA     1.0     2.0   0.5    0.25     0.5
%  MER          1.0     1.0   nan    nan      0.5
%  mCoppens     1.0     nan   0.5    0.25     0.5
%  mAIC         1.0     nan   0.25   nan      0.0

switch picking_method
    case 'Short_Time_Average/Long_Time_Average'
        tof_args ={'Method', picking_method, 'Length_moving_windows',...
            [1.0, 2.5, 0.5, nan], 'Threshold', 0.5};
    case 'Modified_Short_Time_Average/Long_Time_Average'
        tof_args ={'Method', picking_method, 'Length_moving_windows',...
            [1.0, 2.0, 0.5, 0.25], 'Threshold', 0.5};
    case 'Modified_Energy_Ratio'
        tof_args ={'Method', picking_method, 'Length_moving_windows',...
            [1.0, 1.0, nan, nan], 'Threshold', 0.5};
    case 'Modified_Coppens'
        tof_args ={'Method', picking_method, 'Length_moving_windows',...
            [1.0, nan, 0.5, 0.25], 'Threshold', 0.5, 'Beta', 1e-1};
    case 'Modified_AIC'
        tof_args ={'Method', picking_method, 'Length_moving_windows',...
            [1.0, nan, 0.25, nan], 'Threshold', 0.5};
    case 'AR_AIC'
        tof_args ={'Method', picking_method, 'Length_moving_windows',...
            [1.0, nan, 0.25, nan], 'Threshold', 0.5, 'Order', 9};
end


if ~isempty(data_hom)
    
    
    % convert the reference data matrix to a cell array
    data_hom =  mat2cell(data_hom, num_receiver, size(data_hom, 2), ones(num_emitter, 1) );
    
    % allocate a cell array for tofs for the reference data
    tof_hom_cell = cell(num_emitter, 1);
    
    % calcuation of first-arrivals for the reference data
    parfor (ind_emitter = 1:num_emitter, para.nWorkerPool)
        % for ind_emitter = 1:num_emitter % (for test)
        
        
        % remove the DC offset from the reference (in only water) data for the current
        % emitter, and then apply a Tukey window on the time-domain signal.
        data_hom_emitter = bsxfun(@times, tukey, removeDataDc(data_hom{ind_emitter}) );
        
        % allocate zero values
        tof_hom = zeros(num_receiver, 1);
        
        for ind_receiver = 1:num_receiver
            %  disp(ind_receiver);
            if binaries_emitter_receiver(ind_receiver, ind_emitter)
                
                if isfinite(para.Cutoff_freq)
                    
                    %apply a filter in the frequency domain, not
                    %recommended
                    data_hom_emitter(ind_receiver, :) = applyFilter(data_hom_emitter(ind_receiver, :),...
                        fs, para.Cutoff_freq, 'LowPass', 'ZeroPhase', true, 'TransitionWidth', trans_width_prop);
                end
                
                [tof_hom(ind_receiver)] = timeOfFlightPickingEachSignal(data_hom_emitter(ind_receiver,:),...
                    time_array, T, [TOF_min_ref(ind_receiver, ind_emitter),...
                    TOF_max_ref(ind_receiver, ind_emitter)], tof_args{:});
            else
                tof_hom(ind_receiver) = 0;
            end
            
        end
        
        disp(['Computing the first-arrivals for the reference data, Elapsed Time Percentage:' num2str(ind_emitter/num_emitter*100)]);
        
        tof_hom_cell{ind_emitter} = tof_hom;
        
        
    end
    
    % clear the reference data for saving memory
    clear data_hom
    
else
    
    tof_hom_cell = [];
    
end




if ~isempty(data_het)
    
    
    % convert the object data matrix to a cell array
    data_het = mat2cell(data_het, num_receiver, size(data_het, 2), ones(num_emitter, 1) );
    
    % allocate a cell array for tofs for the object data
    tof_het_cell = cell(num_emitter, 1);
    
    % calcuation of first-arrivals for the reference data
    parfor (ind_emitter = 1:num_emitter, para.nWorkerPool)
        % for ind_emitter = 1:num_emitter % (for test)
        
        % remove the DC offset from the object data for the current
        % emitter, and then apply a Tukey window on the time-domain signal.
        data_het_emitter = bsxfun(@times, tukey, removeDataDc(data_het{ind_emitter}) );
        
        tof_het = zeros(num_receiver, 1);
        
        for ind_receiver = 1:num_receiver
            %  disp(ind_receiver);
            if binaries_emitter_receiver(ind_receiver, ind_emitter)
                
                if isfinite(para.Cutoff_freq)
                    
                    %apply a filter in the frequency domain, not
                    %recommended
                    data_het_emitter(ind_receiver, :) = applyFilter(data_het_emitter(ind_receiver, :),...
                        fs, para.Cutoff_freq, 'LowPass', 'ZeroPhase', true, 'TransitionWidth', trans_width_prop);
                end
                
                [tof_het(ind_receiver)] = timeOfFlightPickingEachSignal(data_het_emitter(ind_receiver,:),...
                    time_array, T, [TOF_min(ind_receiver, ind_emitter),...
                    TOF_max(ind_receiver, ind_emitter)], tof_args{:} );
                
            else
                tof_het(ind_receiver) = 0;
            end
            
        end
        
        disp(['Computing the first-arrivals for the phantom data, Elapsed Time Percentage:' num2str(ind_emitter/num_emitter*100)]);
        
        tof_het_cell{ind_emitter} = tof_het;
        
        
        if ~isempty(tof_hom_cell)
            tof_cell{ind_emitter} = tof_het_cell{ind_emitter} - tof_hom_cell{ind_emitter};
        end
        
        
    end
    % clear the main data for saving memory
    clear data_het
    
else
    
    tof_het_cell = [];
    
end


% convert the cells to a stacked vector for all emitters-receivers
if ~isempty(tof_hom_cell)
    tof_hom = cell2mat(tof_hom_cell);
    tof_hom = reshape(tof_hom, [num_receiver, num_emitter]);
else
    tof_hom = [];
end

if ~isempty(tof_het_cell)
    tof_het = cell2mat(tof_het_cell);
    tof_het = reshape(tof_het, [num_receiver, num_emitter]);
else
    tof_het = [];
end

if ~(isempty(tof_hom_cell) || isempty(tof_het_cell))
    tof = cell2mat(tof_cell);
    tof = reshape(tof, [num_receiver, num_emitter]);
else
    tof = [];
end



end