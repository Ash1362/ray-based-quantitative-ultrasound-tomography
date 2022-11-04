function signal_frequency = fourierTransformSampledWindowed(signal, f, time_array, windowing_method)
% Discrete Fourier transform (dft) on the specified frequency, not fft
%
%Input:     data:   recorded time-series data
%           f:      frequency (Hz) at which one wants to do the spectral
%                   analysis
%           delt_t: time increment (s, second) between two neighbouring data
%                   point. delt_t is related to the sampling rate, e.g. if
%                   the total data collection time is T, and the length of
%                   data is n, then delt_t = T/(n - 1).
%           wyn:    'y' -- windowing on data
%                   'n' -- no windowing applied
%Output:    amp:    amplitude of the signal
%           pha:    phase (radian) of the signal
%
% by Hongxue Cai (h-cai@northwestern.edu)

n = length(time_array);
% Hemming window is preferred.
switch windowing_method
    case 'Hamming'
    window = chamming(n);   %Hamming windowing
    case 'Hanning'
    window = channing(n);   %Hanning windowing
    otherwise
    window = 1;
end


signal_frequency = window .* signal * exp(-1i * (f' * time_array)' );
function w = chamming(n)
%Return the n-point Hamming window.
%
if n > 1
    w = .54 - .46*cos(2*pi*(0:n-1)/(n-1));
else
    error('n must be greater than 1.')
end
return
%
%
function w = channing(n)
%Return the n-point Hanning window in a column vector.

w = .5*(1 - cos(2*pi*(1:n)/(n+1)));
return