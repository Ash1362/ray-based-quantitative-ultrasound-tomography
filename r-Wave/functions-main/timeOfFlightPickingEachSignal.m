function [tof] = timeOfFlightPickingEachSignal(signal,t_array, T, time_window_in, varargin)
%TIMEOFFLIGHTPICKINGEACHSIGNAL calculates the first-arrival of a signal
%
%
% DESCRIPTION:
%       selectTimeWindow selects a minimum and maximum time for the first-arrival
%       of the signal
%
%
% USAGE:
%
%
% INPUTS:
%       signal            - a 1D signal (the pressure data that is excited
%                           through an emitter, and is  measured by a
%                           receiver)
%       t_array                - measurement time_array [sec]
%       T                 - the approximate pulse duration of the excitation
%                           pulse [sec]
%      time_window_in     - the minimum and maximum times for the input
%                           window
%
% OPTIONAL INPUTS:
%        'Method'         - the method for calculation of the first-arrival
%                           of the signals. This must be set:
%                          'Modified_AIC'. Other choices are
%                           included, but those are for testing purposes
%                           for developing codes, and the user will get error
%                           for those cases.
%                           'Short_Time_Average/Long_Time_Average',
%                           'Modified_Energy_Ratio',
%                           'Modified_Coppens'
%                           'Modified_Short_Term_Average/Long_Term_Average'
%       'Length_moving_window'  - a vector indicating the length of the
%                                 moving windows:
%                                 1 - backward window, 2 - long window,
%                                 3 - short window, 4 - EPS window
%       'Threshold'             - a threshold used for selection of the
%                                 end point of a fixed window for selection
%                                 of the first-arrival
%       'Beta'                  - a parameter used for smoothing the relative
%                                 energy, using 'Modified Coppens' method
%

% OUTPUTS:
%       tof               - the time of the first-arrival of the signal
%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 


para.Method = 'Modified_AIC';
para.Length_moving_windows = [1, nan, 0.25, nan];
para.Threshold = 0.5;
para.Beta = 1e-1;
para.Order = 8;

% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

if strcmp(para.Method, 'AR_AIC')
    regression_order = para.Order;
end



window_threshold     = para.Threshold;
window_length_factor = para.Length_moving_windows;

% calculate the time spacing of the measurement
dt = t_array(2) - t_array(1);
%% ========================================================================
% CALCULATE THE INTEGER LENGTH OF MOVINNG WINDOWS
% =========================================================================
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

% calculate the length of the window(s) for calculation of the
% first-arrivals. To do that, for each method, a corresponding
% window_length_factor must be used.
window_length = round(window_length_factor*T/dt);



%  length of the signal
length_signal = length(signal);
% length of the backward window
length_back = window_length(1);
% length of the long window
length_long = window_length(2);
% length of the short window if applied)
if isfinite(window_length(3))
    length_short = window_length(3);
end

if isfinite(window_length(4))
    length_eps = window_length(4);
    if round(mod(length_eps, 2)) == 0
        length_eps = length_eps + 1;
    end
end





%% ========================================================================
% SELECT THE TIME WINDOW
% =========================================================================
% get the minimum and maximum index for the time window
time_index_start = round(time_window_in(1)/dt)+1;
time_index_end   = round(time_window_in(2)/dt)+1;

if window_threshold
    
    % calculate the envelope of the signal
    signal_envelope = abs(signal); 
    
    % truncate the envelope of the signal using the starting and ending indices
    % which are calculated above
    signal_envelope = signal_envelope(time_index_start : time_index_end);
    
    % find the first index of the time windowed envelope after which the
    % signal is larger the chosen fraction of the maximum absolute amplitude 
    ixm = find(signal_envelope/max(signal_envelope)>window_threshold, 1, 'first');
    
     % shift the index for getting the time index with respect to time origin
    ixm = ixm + time_index_start - 1;
    
end

% normalise the amplitude of the signal
signal = signal/max(abs(signal));

switch para.Method
    
    case 'Short_Time_Average/Long_Time_Average'
        
        % This method should not be used by the users.
        Error(['The user must use the modified AIC approach'...
            'for time-of-flight picking.'])
        
        if window_threshold
            
            % find a minimum and maximum index for a window for the first-arrival
            % the end point is ixm, and for the starting point, move backward
            % for length_back
            time_index_start  = max(length_long + 1, ixm - length_back);
            time_index_end = ixm;
        end
        
        if time_index_start > time_index_end
            error('The window sizes are too large');
        end
        % truncate the signal
        signal = signal(time_index_start - length_long : time_index_end);
        
        
        % calculate the energy of the signal within the short and
        % long windows
        energy_short_window = zeros(size(signal));
        energy_long_window  = zeros(size(signal));
        
        for i = length_long + 1 : length(signal)
            indices_short = i - length_short : i;
            indices_long = i - length_long : i;
            
            energy_short_window(i) = signal(indices_short)*conj(signal(indices_short))';
            energy_long_window(i) = signal(indices_long)*conj(signal(indices_long))';
        end
        
        energy_short_window = 1/length_short*energy_short_window;
        energy_long_window = 1/length_long*energy_long_window;
        
        attribute = energy_short_window./energy_long_window;
        attribute = attribute(length_long + 1 : end);
        
        
        % find the time index with maximum derivative within this window
        [~,ix] = max(diff(attribute));
        ix = time_index_start - 1 +  ix ;
        % use ix-1, because the first element in the time array is zero
        % (ix = 2 is equivalent to t = dt)
        tof = dt*(ix-1);
        
    case 'Modified_Short_Time_Average/Long_Time_Average'
        
        
        % This method should not be used by the users.
        Error(['The user must use the modified AIC approach'...
            'for time-of-flight picking.'])
        
        if window_threshold
            
            % find a minimum and maximum index for a window for the first-arrival
            % the end point is ixm, and for the starting point, move backward
            % for length_back
            time_index_start  = max (length_long + length_eps, ixm - length_back);
            time_index_end = min (ixm , length_signal - length_eps + 1);
        end
        % truncate the signal
        signal = signal(time_index_start - (length_long + length_eps) + 1 : time_index_end + length_eps - 1);
        
        % calculate the energy of the signal within the short and
        % long windows
        energy_short_window = zeros(size(signal));
        energy_long_window  = zeros(size(signal));
        
        for i = length_long + 1 : length(signal)
            indices_short = i - length_short : i;
            indices_long = i - length_long : i;
            
            energy_short_window(i) = signal(indices_short)*conj(signal(indices_short))';
            energy_long_window(i) = signal(indices_long)*conj(signal(indices_long))';
        end
        
        energy_short_window = 1/length_short*energy_short_window;
        energy_long_window = 1/length_long*energy_long_window;
        
        attribute = energy_short_window./energy_long_window;
        attribute = attribute(length_long + 1 : end);
        
        
        % apply EPS filter in order to enhance sharp changes
        [ attribute_eps ] = epsFilter(attribute, length_eps);
        
        % find the time index with maximum derivative within this window
        [~,ix] = max(diff(attribute_eps));
        ix = time_index_start - 1 + ix;
        
        % use ix-1, because the first element in the time array is zero
        % (ix = 2 is equivalent to t = dt )
        tof = dt*(ix-1);
        
    case 'Modified_Energy_Ratio'
        
        % This method should not be used by the users.
        Error(['The user must use the modified AIC approach'...
            'for time-of-flight picking.'])
 
        if window_threshold
            
            % find a minimum and maximum index for a window for the first-arrival
            % the end point is ixm, and for the starting point, move backward
            % for length_back
            time_index_start  = max (length_long + 1 , ixm - length_back);
            time_index_end = min (ixm , length_signal - length_long);
        end
        
        signal =   signal(time_index_start - length_long : time_index_end + length_long);
        
        % calculate the energy of the signal within the short window
        energy_forward_window = zeros(size(signal));
        energy_backward_window = zeros(size(signal));
        for i = length_long + 1 : length(signal) - length_long
            
            indices_forward = i : i + length_long;
            indices_backward = i - length_long : i;
            
            energy_forward_window(i) = signal(indices_forward)*conj(signal(indices_forward))';
            energy_backward_window(i) = signal(indices_backward)*conj(signal(indices_backward))';
            
        end
        
        
        energy_relative = energy_forward_window./energy_backward_window;
        
        attribute = (abs(signal(length_long + 1 : length(signal) - length_long)).*...
            energy_relative(length_long + 1 : length(signal) - length_long)).^3;
        
        % find the time index of the maximum of the signal within this window
        [ ~ , ix ] = max( attribute );
        ix = ix + time_index_start - 1;
        
        % use ix-1, because the first element in the time array is zero
        % (ix = 2 is equivalent to t = dt )
        tof = dt*(ix-1);
        
        
    case 'Modified_Coppens'
        
        % This method should not be used by the users.
        Error(['The user must use the modified AIC approach'...
            'for time-of-flight picking.'])
        
        if window_threshold
            
            % find a minimum and maximum index for a window for the first-arrival
            % the end point is ixm, and for the starting point, move backward
            % for length_back
            time_index_start  = max (length_short + length_eps - 1, ixm - length_back );
            time_index_end = min (ixm, length_signal - length_eps + 1);
        end
        
        % calculate the energy of the signal within the long window
        energy_long_window = cumsum(signal.*conj(signal));
        
        % truncate the signal and the long window
        signal = signal(time_index_start - (length_short+length_eps) + 2 : time_index_end + length_eps - 1);
        energy_long_window = energy_long_window(time_index_start - (length_short+length_eps) + 2 : time_index_end + length_eps - 1);
        
        
        energy_short_window = zeros(size(signal));
        for i = length_short:length(signal)
            indices  =  i - length_short + 1 : i;
            energy_short_window(i) = signal(indices)*conj(signal(indices))';
        end
        
        
        attribute = energy_short_window./(energy_long_window + para.Beta);
        
        % apply EPS filter in order to enhance sharp changes
        attribute_eps = epsFilter(attribute(length_short:length(signal)), length_eps);
        
        
        % find the time index with maximum derivative within this window
        [~,ix] = max( diff( attribute_eps ));
        ix = ix + time_index_start - 1;
        
        tof = dt*( ix - 1 );
        
        
    case 'Modified_AIC'
        
        
        if window_threshold
            
            % find the minimum and maximum indices for the time window for finding
            % the minimum AIC. The end point is ixm, and the starting point is
            % ixm - length_back.
            time_index_start = max (1, ixm - length_back);
            time_index_end = ixm;
            
        end
        
        % truncate the time array for the chosen window
        t_array_window = t_array(time_index_start:time_index_end);
        
        % allocate a vector for Akaike-Information-Criterion (AIC)
        AIC = zeros(size(t_array_window));
        
        % get the length of the window for computing the AIC valuse
        length_window = length(t_array_window);
        
        for i = 1:length_window
            k = i + time_index_start - 1;
            variance_noise = var(signal(1:k));
            variance_signal = var(signal(k+1:length_signal));
            AIC(i) = k*log(variance_noise) + (length_signal-1-k)*log(variance_signal);
        end
        
        % calculate the minimum AIC within the time window
        [AICm, ixm] = min(AIC);
        
        % select the short window for picking the first arrival
        indices =  max(1, ixm - length_short): min(length_window, ixm + length_short);
        
        % compute the exponential of discrepancy of the AIC and the calculated minimum AIC
        % in the chosen short time window about the minimum AIC
        exponal_discrepancy_akaike = exp(-(AIC(indices)-AICm)/2);
        
        % calculate the Akaike weights for each data sample within the time window
        weights = exponal_discrepancy_akaike/sum(exponal_discrepancy_akaike);
        
        % calculate a weighted average TOF using the Akaike weights
        tof = weights * t_array_window(indices)';
        
    case 'AR_AIC'
        
        % This method should not be used by the users.
        Error(['The user must use the modified AIC approach'...
            'for time-of-flight picking.'])
        
        if window_threshold
            % find a minimum and maximum index for a window for the first-arrival
            % the end point is ixm, and for the starting point, move backward
            % for length_back
            time_index_start  = max(length_long + 1, ixm - length_back);
            time_index_end = ixm;
        end
        
        if time_index_start > time_index_end
            error('The window sizes are too large');
        end
        
        
        
        % truncate the time array using the selected window
        t_array_window = t_array(time_index_start:time_index_end);
        
        % allocate a vector for AIC values
        AIC = zeros(size(t_array_window));
        length_window = length(t_array_window);
        
        
        for i = 1:length_window
            k = i + time_index_start - 1;
            [ ~ , e_w1 ] = aryule( signal(1:k), regression_order);
            [ ~ , e_w2 ] = aryule( flip (signal(k+1:length_signal)), regression_order);
            AIC(i) = (k - regression_order)*log(e_w1) + (length_signal - regression_order - k)*log(e_w2);
        end
        
        % get the length of the window for computing the AIC valuse
        length_window_short = length(t_array_window);
        
        % calculate the minimal AIC value within the time window
        [AICm, ixm] = min(AIC);
        
        % select the short window for section of first arribval
        indices =  max(1, ixm - length_short): min(length_window_short, ixm + length_short);
        
        % calculate the exponal discrepancy of the AIC for the data
        % within the time window indicated by indices and the calculated minimum AIC
        exponal_discrepancy_akaike = exp(-(AIC(indices)-AICm)/2);
        
        % calculate the Akaike weights for each data sample within the time window
        weights = exponal_discrepancy_akaike/sum(exponal_discrepancy_akaike);
        
        % calculate a weighted average TOF using the Akaike weights
        tof = weights*t_array_window(indices)';
        
end



end