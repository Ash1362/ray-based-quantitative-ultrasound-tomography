function [data] = removeDataDc(data)
% REMOVEDATADC removes the DC part of pressure time series

% DESCRIPTION:
% removeDataDc removes the DC part of pressure time series, an is extracted
% from the function 'removeDCoffset.m' written by Brad Treeby at UCL)
%
% USAGE:
%
% INPUTS:
%       data             - the num_receiver x num_time matrix of
%                          ultrasound pressure time series
% OUTPUTS:
%       data_processed    - the processed pressure data

%
% ABOUT:
%       author          - Brad Treeby, Ashkan Javaherian
%       date            - 30.06.2019
%       last update     - 30.06.2019
%
%
% See also 



% get the length of the time array
num_time = size(data, 2);

% get the filter for removing the DC offset
filter_vector = ones(1, num_time);
if rem(num_time, 2)
    filter_vector((end-1)/2+1) = 0;
else
    filter_vector(end/2+1) = 0;
end


% apply the filter on the data, and remove the DC offset
data = real(ifft(bsxfun(@times, ifftshift(filter_vector), fft(data, [], 2)), [], 2));

end