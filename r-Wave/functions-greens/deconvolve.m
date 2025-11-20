function signal_frequency = deconvolve(signal_in, signal_out, time_array,...
    omega, deconvolution_parameter)
%DECONVOLVE computes the deconvolution of two time-domain signals in the
%frequency domain
% DESCRIPTION:
% deconvolve deconvolves the signal_in from signal_out in the frequency domain.
% USAGE:
%
% INPUTS:
%       signal_in            - the 1 x num_time input signal in the time domain
%       signal_out           - the N x num_time output signal in the time domain
%       time_array           - time array [s]
%       omega                - the angular frequency array on which the
%                              deconvolution is done
%       deconvolution_parameter - A small scalar representing the parameter
%                                  for regularisation. The zeros in the frequency
%                                  domain version of signal_in will make the
%                                  deconvolution operator ill-conditioned,
%                                  so regularisation is required for avoiding
%                                  ill-conditioning.
%
%
%
% OUTPUTS:
%       signal_frequency    - the deconvolved signal in the frequency domain

% ABOUT:
%       author            - Ashkan Javaherian
%       date              - 25.07.2020
%       last update       - 25.07.2020
%
% This script is part of the r-Wave Tool-box.
% Copyright (c) 2022 Ashkan Javaherian

if size(signal_in,1) > 1 
    error(['The code currently accepts only a single pressure source for all'...
        'the emitters.']);
end
    

% get the length of time array
num_time = length(time_array);

if size(signal_in, 2) < num_time
    
    % make the length of the signal and time array consistent
    signal_in = [signal_in, zeros(1, num_time - length(signal_in))];
    
end


% compute the signal_in in the frequency domain
signal_in_frequency =  2/num_time * signal_in * exp(1i * (omega' * time_array)');

% compute the frequency domain version of signal_out and deconvolve it from
% signal_in (all the given traces simultaneously)
signal_frequency = 2/num_time *(signal_out * exp(1i * (omega' * time_array)')) .*...
    (conj(signal_in_frequency)./(abs(signal_in_frequency).^2 + ...
    (deconvolution_parameter * max(abs(signal_in_frequency))).^2));


end