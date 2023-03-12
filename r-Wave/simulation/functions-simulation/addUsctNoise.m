function [data_noisy] = addUsctNoise(data, varargin)
%ADDUSCTNOISE add noise to each trace of the USCT simulated data

% DESCRIPTION:
%       AddUsctNoise add noise to each trace of the USCT simulated data
%      
% USAGE:
%       [data_noisy] = AddUsctNoise(data,SNR,varargin)
%
% INPUTS:
%       data - the simulated data

%       
% OPTIONAL INPUTS:
%       Type - the kWave mode for adding noise, 'peak' or 'rms'
%       SNR  - a desired signal to noise ratio for the simulated data    
% OUTPUTS:
%       data_noisy - the noisy data with the sime size as the input data
% ABOUT:
%       author            - Ashkan Javaherian
%       date              - 30.12.2019
%       last update       - 30.12.2019
%
% See also addnoise.m
% ABOUT:
%       author           - Ashkan Javaherian
%       date             - 15.12.2019
%       last update      - 15.12.2019
%       
%
% This function is part of the r-Wave Toolbox
% Copyright (C) 2022 Ashkan Javaherian




% set kWave default parameter
para      = [];
para.NoiseType = 'peak';
para.SNR = 40;
para.nworker_pool = 16;
% read additional arguments or overwrite default ones (no sanity check!)
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% get the number of emitters   
num_emitter = size(data,3);

% get the number of receivers
num_receiver = size(data,1);

% get the number of time smaples
num_time = size(data, 2);

% convert the original noise-free data to a cell array
data = mat2cell(data, num_receiver, num_time, ones(num_emitter, 1));

% get the cell array for the noise-contaminated data
data_noisy = cell(size(data));

parfor (ind_emitter = 1:num_emitter, para.nworker_pool)
    
disp(['adding noise:' num2str(ind_emitter)])
    for ind_receiver = 1:num_receiver
    
    % add noise to the time trace using k-Wave    
    data_noisy{ind_emitter}(ind_receiver, :) = addNoise(data{ind_emitter}(ind_receiver, :),...
        para.SNR, para.NoiseType);

    end
end

% convert the cell array for the noise-contaminated data to a mat array
data_noisy = cell2mat(data_noisy);    


end