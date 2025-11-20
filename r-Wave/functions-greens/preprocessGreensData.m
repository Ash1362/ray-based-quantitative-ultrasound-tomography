function [data_freq, source_freq, binary_data] = preprocessGreensData(...
    data_pressure, time_array, omega, emitter, receiver, sound_speed_water,...
    minimum_distance, varargin)
%PREPROCESSGREENSDATA PREPROCESS THE ULTRASOUND DATA AND TRANSFORM IT TO
%A GIVEN DISCRETISED FREQUENCY DOMIAN
%
% DESCRIPTION:
% PREPROCESSGreensData preprocess the ultrasoud data and excition pulse and
% transform it from a given discretised time domain  to a given discretised
% frequency domain using a discrete Fourier transform operator

% USAGE:
%        [img, out, para] = reconstructGreensImage(data, time_array, recon_grid,...
%        emitter, receiver, sound_speed_water, img, ray_initial_angles, simulation_prop,...
%        directories, varargin)
% INPUTS:
%       data_pressure      - a num_receiver x num_time x num_emitter
%                            matrix representing the ultrasound data
%       time_array         - a 1 x num_time vector representing the time
%                            array
%       omega              - a 1 x num_freq matrix representing the
%                            frequency array, i.e., frequencies on which
%                            the image is reconstructed.
%       emitter            - a struct containing the fields:
%       'pulse'            - a time-varying signal representing the source
%       'positions_real'     the real position of emitters
%       'positions'          the projected position of emitters to a circle
%                            with fixed radius
%       receiver           - a struct containing the fields
%       'positions_real'     the real position of receivers
%       'positions'          the projected position of receivers to a circle
%                            with fixed radius
%       sound_speed_water  - the sound speed [m/s] in water
%       minimum_distance   - the minimum distance of the transducers for including
%                            in the image reconstruction. The emitter-receiver pairs with
%                            distances less than minimum_distance are excluded from
%                            image reconstruction.
%                            The measured signals associated with close emitter-receiver
%                            pairs may be deteriorated by the effects of directivity
%                            of the transducers. In addition, rays linking the
%                            close emitter-receiver pairs only travel inside water.
%
%
% OPTIONAL INPUTS:
%      'num_workers'       - The number of workers for parallel programming
%      'greens_domain'     - the domain for approximating the pressure field
%                            using Green's function, this can be set 'time',
%                            or 'frequency'
%      'greens_method'     - the method for calculating the pressure
%                            field, This can be 'analytic',...
%      'deconvolve_source' - Boolean controlling whether the measured signals
%                            are deconvolved from the pressure source or
%                            not. (default:true)
%      'deconvolution_parameter' - get the regularisation parameter for
%                            deconvolution (default : 0)                      
%      'filter_data'        - Boolean controlling whether the data is
%                             filtered in the frequency domain or not.
%                             (default : false)
%      'outliers_fraction'  - the fraction of the emitter-receiver pairs
%                             excluded from image reconstruction for a specific
%                             frequency, because their magnitudes are outliers.
%                             (default: 0.01) 
%      'multiply_data_window'- Boolean controlling whether the time series 
%                              are multiplied by an expoentially-varying
%                              time window or not (default : false)
%
% OUTPUTS:
%       data_freq          - a num_emitter x num_freq x num_receiver matrix
%                            representing the frequency-domain ultrsasound
%                            data                          
%       source_freq        - a 1 x num_freq signal reepresenting the source
%                            in the frequency domain
%       binary_data        - a cell array of length num_emitter, each
%                           containing a num_receiver x num_freq matrix 
%                           of accepted data for including in the image
%                           reconstruction.

% ABOUT:
%       author            - Ashkan Javaherian
%       date              - 25.07.2020
%       last update       - 25.07.2020
%
% This script is part of the r-Wave Tool-box
% Copyright (c) 2022 Ashkan Javaherian

para = [];
para.num_workers = 16;
para.greens_domain = 'frequency';
para.greens_method = 'analytic';
para.deconvolve_source = false;
para.deconvolution_parameter = 0;
para.filter_data = false;
para.outliers_fraction = 0.01;
para.multiply_data_window = false;


if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% initialise the run time
ts = tic;

% get the length of the time arrays
num_time = length(time_array);

% get the length of the frequency array
num_freq = length(omega);

% get the number of emitters 
num_emitter = size(emitter.positions, 2);

% get the number of receivers
num_receiver = size(receiver.positions, 2);


% If the real position of transducers are given, the time delays due to the
% real position of the transducers are corrected.
do_correct_time_delays = isfield(emitter, 'positions_real');


%%=========================================================================
% COMPUTE THE TIME CORRECTIONS ASSOCIATED WITH PROJECTING TRANSDUCERS'
% POSITIONS
%==========================================================================
if do_correct_time_delays

% compute the only-water distances [m] for the real position of emitter-receiver
% pairs
distances_real = calculateDistanceEmitterReceiver(emitter.positions_real,...
    receiver.positions_real, []);
end

% compute the distances between the projected position of emitter-receiver
% pairs onto the detection ring (surface)
distances = calculateDistanceEmitterReceiver(emitter.positions,...
    receiver.positions, []);


if length(size(data_pressure))< 3

    % if the data is given for one emitter, expand the dimensions of data
    data_pressure = reshape(data_pressure, [size(data_pressure),1]);

end

% check the size of the given forward pressure data
if any(size(data_pressure) ~= [num_receiver, num_time, num_emitter])
    error('The size of the given forward pressure data is not consistent.')
end



%==========================================================================
% APPLY FILTER ON THE ULTARSOUND SIGNALS
%==========================================================================

if para.filter_data

    % choose the cut-off frequencies for the filter
    cutoff_freq = [0, 4] * 1e6;
else

    cutoff_freq = nan;

end



% filter the data in the frequency domain, if requested
if all(isfinite(cutoff_freq))

    cutoff_freq(~cutoff_freq) = [];

    % choose the method for filtering based on the chosen frequency intervals
    if length(cutoff_freq) < 2

        % if the first component is zero, apply the 'low-pass' filter
        filter_mode = 'low-pass';

    else
        % if the both components are finite, apply the 'band-pass' filter
        filter_mode = 'band-pass';

    end

    % set the optional inputs
    data_filter_args = {'Mode', filter_mode, 'nworker_pool', para.num_worker_pool};

    % apply the filter on the data
    data_pressure = filterPressureData(data_pressure, cutoff_freq, time_array, data_filter_args{:});


end



%%=========================================================================
% MULTIPLY A TIME WINDOW BY THE ULTRASOUND TIME SERIES 
%==========================================================================
if para.multiply_data_window


    % get the time-window edge for the times before the arrival times
    if do_correct_time_delays
        time_window_pre = 0.05/(sound_speed_water + 60) * max(distances_real(:));
    else
        time_window_pre = 0.05/(sound_speed_water + 60) * max(distances(:));
    end

% get the coefficient for a time window for the times after the arrival
% times
time_window_post = inf;


for ind_receiver = 1:num_receiver

    % get the tofs for the current receiver
    if do_correct_time_delays
        tofs_receiver = distances_real(ind_receiver,:)/(sound_speed_water + 60);
    else
        tofs_receiver = distances(ind_receiver,:)/(sound_speed_water + 60);
    end

    % make a mesh grid from TOFs associated with the current receiver
    [tofs_receiver, time_matrix] = meshgrid(tofs_receiver, time_array);

    % get an exponential time window 
    time_window = exp(-(1/2)*(subplus(time_matrix-tofs_receiver)/time_window_post + ...
       subplus(tofs_receiver- time_matrix)/time_window_pre).^2);
   

    % multiply the time window by the signals for the current receiver
    data_pressure(ind_receiver,:,:) = time_window .* squeeze(data_pressure(ind_receiver,:,:));
end

end

%%=========================================================================
% TRANSFORM THE TIME SERIES TO THE FREQUENCY DOMAIN
%==========================================================================

% get the excitation pulse
switch para.greens_domain

    case 'time'
        error('Not implemented yet.')
    case 'frequency'

        switch para.greens_method


            case  'analytic'

                if para.deconvolve_source

                    % deconvolve the pressure time series from the
                    % excitation pulse
                    transform_frequency = @(data) deconvolve(emitter.pulse,...
                        data, time_array, omega, para.deconvolution_parameter);

                    % make the source 1, if the data is deconvolved
                    % from the source
                    source_freq = 1;


                else

                    % Define a handle function for computing the Discrete
                    % Fourier Transform (DFT) of a set of time series data
                    transform_frequency = @(data) discreteFourierTransform(...
                        data, omega, time_array);

                    % Get the DFT of the pressure source
                    source_freq = transform_frequency([emiter.pulse,...
                        zeros(1, num_time-length(emitter.pulse))]);


                end

        end




        % allocate a 3D matrix for the frequency-domain ultrasound data
        data_freq = zeros(num_receiver, num_freq, num_emitter);

        for ind_emitter = 1: num_emitter

            % transform the data from the given time domain to the
            % frequency domain
            data_freq(:,:, ind_emitter) = transform_frequency(data_pressure(:,:, ind_emitter));

        end

end

%%=========================================================================
% APPLY TIME CORRECTIONS
%==========================================================================
if do_correct_time_delays

% get the discrepancy of the computed TOfs [s] between the projected and real 
% position of emitter-receiver pairs
tofs_water_correction = 1/(sound_speed_water + 60)  * (distances - distances_real);

% display the corrected TOFs
figure;imagesc(tofs_water_correction); colorbar;

% apply the time delay corrections to the time series in the frequency
% domain
for ind_freq = 1: num_freq
data_freq(:, ind_freq, :) = reshape(exp(1i * omega(ind_freq) * tofs_water_correction) .*...
    squeeze(data_freq(:, ind_freq, :)), [num_receiver, 1, num_emitter]);
end

end

%%=========================================================================
% GET A MAP FOR THE ACCEPTED TIME SERIES AND REMOVING OUTLIERS 
%==========================================================================

% get a cell array of length num_emitter, each a num_receiver x 
% num_freq matrix representing the accepted time traces for each frequency

% get a binary map based on a minimum distance for the accepted emitter-receiver pairs
% The close emitter-receiver pairs usually do not provide good signals
% because of directivity effects. This binary map is the same for all
% frequencies.
binary_data_allfrequencies = distances > minimum_distance;

% get the number of outliers
num_outliers = round(para.outliers_fraction * num_receiver * num_emitter);

% get a binary map for accepting the time series 
% This binary map is determined at each freqency separately.
binary_data = true(size(data_freq));

for ind_freq = 1: num_freq
   
  % get the magnitudes for the current emitter
   magnitudes_current = binary_data_allfrequencies .* abs(squeeze(data_freq(:, ind_freq, :)));

   % detemine the ouliers as a portion (here 1%) of greater magnitudes
   [~, outliers_indices] = maxk(magnitudes_current(:), num_outliers);
   
   % make the outliers zero
   magnitudes_current(outliers_indices) = 0;

   % get the data binary map for the current frequency
   binary_data(:,ind_freq,:) = reshape(magnitudes_current > 0,...
       [num_receiver, 1, num_emitter]);

end

% For testing the 'equidistant angle' approach for projecting transducers
% to the fixed-radius ring, the last emitter and receiver must be removed
% because of large error caused in projection.
% binary_data(end,:,:) = false;
% binary_data(:,:,end) = false;

% make the binary matrix for accepting time traces a cell array of length
% num-emitter, each containg a matrix of size num_receiver x num_freq
binary_data = squeeze(mat2cell(binary_data, num_receiver, num_freq, ones(num_emitter, 1)));

% compute the whole run time for preprocessing the ultrasound data
time_run = toc(ts);

disp(['The time for preprocessing the data was' num2str(time_run, '%.2f') 'sec.'])

end