function signal_frequency = discreteFourierTransform1(signal, omega, time_array)
%DISCRETEFOURIERTRANSFORM calculates the Discrete Fourier transform (dft) of
% the recorded time-series data

% DESCRIPTION:
% discreteFourierTransform calculates the Discrete Fourier transform (dft) of
% the pressure time series on the given discretised angular frequencies
%
% USAGE:
%
% INPUTS:
%       signal             - the input signals in the time domain
%       omega              - a vector of angular frequencies [rad/s]
%       time_array         - time array [s]
% 
%
% OPTIONAL INPUTS:
%
%
%
% OUTPUTS:
%       signal_frequency    - the signals in the frequency domain
%
% ABOUT:
%       author            - Ashkan Javaherian
%       date              - 25.07.2020
%       last update       - 25.07.2020
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian


% get the length of time array
num_time = length(time_array);

if size(signal, 2) < num_time
    
    % make the length of the signal and time array consistent
    signal = [signal, zeros(1, num_time-length(signal))];
end

% calculate the DFT of the signals
signal_frequency =  signal * exp(-1i * (omega' * time_array)' );