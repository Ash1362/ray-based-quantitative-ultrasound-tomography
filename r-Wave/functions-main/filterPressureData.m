function [data_filtered] = filterPressureData(data, cutoff_frequencies, dt, varargin)
%FILTERPRESSUREDATA applies a frequency filter on the data set
%
%
% DESCRIPTION:
%  filterPressureData applies a frequency filter on the data set
%
% USAGE:
%
% INPUTS:
%       signal             - the input signals in the time domain
%       cutoff_frequencies - cut-off frequencies
%       dt                 - time spacing [s]
%
%
% OPTIONAL INPUTS:
%      Mode                - the mode for frequency filtering, which can be
%                           'low-pass' or 'band-pass'
%      'nworker_pool'      - the number of workers
%
%
% OUTPUTS:
%       data_filtered    - the filtered signals

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 18.03.2020
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2020 Ashkan Javaherian

para = [];
para.Mode = 'low-pass';
para.nworker_pool = 16;


if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% get the transition width
trans_width_MHz = 10e6;

% get the frequency
fs = 1/dt;


trans_width_prop = trans_width_MHz / fs;


% number of receivers
num_receiver = size(data, 1);

% number of emitters
num_emitter = size(data, 3);

% number of time indices
num_time = size(data, 2);


if num_emitter < para.nworker_pool
    para.nworker_pool = 0;
end

% convert the data matrix to a cell array
data = mat2cell(data, num_receiver, num_time, ones(num_emitter, 1) );


% allocate a cell array for the filtered data
data_filtered = cell(1, num_emitter);


parfor (ind_emitter = 1 : num_emitter, para.nworker_pool)
    
    switch para.Mode
        case 'low-pass'
            
            % apply the low-pass filter on the time traces
             data_filtered{ind_emitter} = lowpass(data{ind_emitter}', cutoff_frequencies, fs);
            %   for ind_receiver = 1:num_receiver
            %data_filtered{ind_emitter}(ind_receiver, :) = applyFilter(data{ind_emitter}(ind_receiver, :),...
            %   fs, freq_pass, 'LowPass', 'ZeroPhase', true, 'TransitionWidth', trans_width_prop);
            
            
            
            %   end
        case 'band-pass'
            
            % apply the band-pass filter on the time traces
             data_filtered{ind_emitter} = bandpass(data{ind_emitter}', cutoff_frequencies, fs);
            % for ind_receiver = 1:num_receiver
            %data_filtered{ind_emitter}(ind_receiver, :) = applyFilter(data{ind_emitter}(ind_receiver, :),...
             %   fs, freq_pass, 'BandPass', 'ZeroPhase', true, 'TransitionWidth', trans_width_prop);
            % end
    end
    
    disp(['The emitter for which data is being filtered is :' num2str(ind_emitter)])
    
end


% convert the cells to a matrix
data_filtered = permute(reshape(cell2mat(data_filtered),...
    [num_time, num_receiver, num_emitter]), [2, 1, 3]);


end