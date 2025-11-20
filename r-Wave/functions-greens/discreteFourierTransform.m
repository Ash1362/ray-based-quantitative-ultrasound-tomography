function signal_frequency = discreteFourierTransform(signal, omega, time_array, do_window)
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
%       do_window          - Boolean controlling whether the windowing is
%                            applied for mitigating the frequency leakage
%                            or not. 
%
% OPTIONAL INPUTS:
%
%
%
% OUTPUTS:
%       signal_frequency    - the signals in the frequency domain

% ABOUT:
%       author            - Ashkan Javaherian
%       date              - 25.07.2020
%       last update       - 25.07.2020
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian


if nargin < 4
    
    % Do apply window, if the last input is not given
    do_window = false;
end

% get the length of time array
num_time = length(time_array);

% get the time spacing [s]
time_spacing = time_array(2) - time_array(1);

if size(signal, 2) < num_time
    
    % make the length of the signal and time array consistent
    signal = [signal, zeros(1, num_time-length(signal))];
    
end


if do_window
    
% get the time window
% time_window = getWin(num_time, 'Tukey', 'Param', 0.25).';    
time_window = chamming(num_time).';
 
% compute $e^(i \omega t)$
omegat = exp(1i * (omega' * time_array)');  

% apply the window on the amplitude, and compute the signal in the
% frequency domain  
signal_frequency = 2 * abs(2/num_time * (time_window .* signal) * omegat).*...
exp(1i*angle(2/num_time * signal * omegat));

%signal_frequency = 2 * abs(time_spacing * (time_window .* signal) * omegat).*...
%exp(1i*angle(time_spacing  * signal * omegat));



else

% compute the signal in the frequency domain
signal_frequency = 2/num_time * signal * exp(1i * (omega' * time_array)' );
%signal_frequency = time_spacing * signal * exp(1i * (omega' * time_array)' );

end

function w = chamming(n)
%return the n-point Hamming window.
%
if n > 1
    w = .54 - .46*cos(2*pi*(0:n-1)'/(n-1));
else
    error('n must be greater than 1.')
end
return

function w = channing(n)
%return the n-point Hanning window in a column vector.

w = .5*(1 - cos(2*pi*(1:n)'/(n+1)));
return