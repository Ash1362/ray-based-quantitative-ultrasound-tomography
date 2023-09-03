function [data, data_ref, data_noisy, data_ref_noisy,...
    transducer] = subsampleData(data, data_ref,...
    data_noisy, data_ref_noisy, transducer, sampling_factor, mode)
%SUBSMAPLEDATA subsamples the transducers, and the corresponding data
%sets
%
%
% DESCRIPTION:
% subsampleData is used for image reconstruction with a subsample of
% excitations (emitters). The emitters are subsmapled homogeneously by the
% sampling_factor
%
% USAGE:
%
%
% INPUTS:
%       data                   - the clean object-in-water data of size
%                                num_receiver x num_time x num_emitter
%       data_ref               - the clean only-water data of size
%                                num_receiver x num_time x num_emitter
%       data_noisy             - the noise-contaminated object-in-water data of size
%                                num_receiver x num_time x num_emitter
%       data_ref_noisy         - the noise-contaminated only-water data of
%                                size num_receiver x num_time x num_emitter
%       transducer             - an object defining the transducers to be
%                                subsampled
%       sampling_factor        - a natural number for subsampling the
%                                excitations
%       mode                   - the mode of subsampling, which can be
%                                'receiver','time', or 'emitter'

% OPTIONAL INPUTS:
%

% OUTPUTS:
%       data                   - the downsmapled clean object-in-water data 
%       data_ref               - the downsmapled clean only-water data 
%       data_noisy             - the downsmapled noise-contaminated object-in-water data
%       data_ref_noisy         - the downsmapled noise-contaminated only-water data
%       transducer             - the sampled transducers
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 29.12.2019
%       last update     - 29.12.2019
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2020 Ashkan Javaherian



if any(strcmp(mode,{'emitter', 'receiver'})) && isempty(transducer)
    error(['For subsmapling emitters or receivers, the object for transducers'...
        'must be given as an input.'])
end


if iscell(transducer)
    error(['The downsampling of the receivers in a rotational setting is not'...
        'supported yet.'])
end


switch mode
    case 'receiver'
        
        % subsample the noise contaminated data (used for image reconstruction)
        % along the receivers
        if ~isempty(data_noisy)
        data_noisy = data_noisy(sampling_factor:sampling_factor:end, :, :);
        end
        if ~isempty(data_ref_noisy)
        data_ref_noisy = data_ref_noisy(sampling_factor:sampling_factor:end, :, :);
        end
        
        if ~isempty(transducer)
            
        % subsample the receivers
        transducer.positions = transducer.positions(:, sampling_factor:sampling_factor:end);
        end
        
        % subsample the original clean data along the receivers
        if ~isempty(data)
            data = data(sampling_factor:sampling_factor:end, :, :);
        end
        if ~isempty(data_ref)
            data_ref = data_ref(sampling_factor:sampling_factor:end, :, :);
        end
        
    case 'emitter'
        
        % subsample the noise contaminated data (used for image reconstruction)
        % along the emitters
        if ~isempty(data_noisy)
            data_noisy = data_noisy(:, :, sampling_factor:sampling_factor:end);
        end
        if ~isempty(data_ref_noisy)
            data_ref_noisy = data_ref_noisy(:, :, sampling_factor:sampling_factor:end);
        end
        
        if ~isempty(transducer)
            
        % subsample the emitters
        transducer.positions = transducer.positions(:, sampling_factor:sampling_factor:end);
        end
        
        % subsample the original clean data along the emitters
        if ~isempty(data)
            data = data(:, :, sampling_factor:sampling_factor:end);
        end
        if ~isempty(data_ref)
            data_ref = data_ref(:, :, sampling_factor:sampling_factor:end);
        end
    case 'time'
        
        % subsample the noise contaminated data (used for image reconstruction)
        % in time
        if ~isempty(data_noisy)
        data_noisy = data_noisy(:, 1:sampling_factor:end, :);
        end
        if ~isempty(data_ref_noisy)
        data_ref_noisy = data_ref_noisy(:, 1:sampling_factor:end, :);
        end
        
        % subsample the original clean data in time
        if ~isempty(data)
            data = data(:, 1:sampling_factor:end, :);
        end
        if ~isempty(data_ref)
            data_ref = data_ref(:, 1:sampling_factor:end, :);
        end
        
end


end