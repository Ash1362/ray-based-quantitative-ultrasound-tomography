function [tof_data, weight_ray] = calcLoadTofs(data_paths_object, data_water, emitter, receiver,...
    time_info, emit_shot_index, us_shot_index, emit_firing_transducer_index,...
    rotation_indices, sound_speed_object_correction_coeff, refractive_background_object_alldata,...
    sound_speed_object_water, sound_speed_only_water, user_data_path, varargin)
%CALCLOADTOFS computes the time-of-flights (TOFs)
%   
%
% DESCRIPTION:
%             calcLoadTofs computes the time-of-flights, or alternatively 
%             loads the time-of-flights from the user data path.
%      
%
% USAGE:
%     
%
% INPUTS:
%       data_paths_object - the path in which the object-in-water data are
%                           stored.
%       data_water      - the num_receiver x num_time x num_receiver 
%                         only-water data. This data contain a single full rotation
%                         of excitations over all the transducers.
%       emitter         - a struct with field
%       'positions'     - the position of emitters
%       'shot_time'     - the US shot time, which is used as the time
%                         origing for computing the time-of-flights,
%                         i.e., tof = first arrival_time - us shot time,
%                         where us shot time is fixed for all emitters.
%       receiver        - a struct with field
%       'positions'     - the position of receivers
%       time_info       - the information about time
%    emit_shot_index    - the vector including the indices of excitions
%    us_shot_index         - the vector including the indices of US shots 
%    rotation_indices      - the vector including the indices of angular
%                            positions for exitations
%    sound_speed_object_correction_coeff - the 1 x num_emitter vector of
%                             correction coeffficients applied to the
%                             computed time-of-flights for the only-water
%                             data for correcting the temperature changes
%    refractive_background_object_alldata - the (num_emitter *
%                            num_receiver) x 1 vector of the refractive index of water
%                            encomapassing the object for all exciations/receivers
%                            (all time series)
%    sound_speed_object_water - A scalar value representing the mean sound
%                               speed of water encompassing the object
%    sound_speed_only_water   - A scalar value representing the sound speed [m/s]
%                               of water for the only-water data
%    user_data_path      - the information of the UST data, together with the 
%                           computed TOFs, are saved in the user_data_path. 
%                           So the user can reconstruct the image without
%                           the need for the main data by setting
%                           do_calculate_tofs false.
%
%       
% OPTIONAL INPUTS:
%      num_worker_pool    - the number of workers for parallel programming
%                           on CPUs 
%      do_calculate_tofs  - Boolean controlling whether the time-of-flights
%                          (TOFs) are calculated, or the already calculated
%                           TOFs are loaded. This Boolean will be automatically
%                           set true, if the TOFs do not exist in the
%                           defined path. (Default:true)
%      include_only_water  - Boolean controlling whether the only-water
%                            data is included in the image reconstruction or not. 
%      choose_last_us_position - A parameter determining how to choose the last
%                           US position, which can 'early' (chosen by the user),
%                           or automatically set all the angular US positions - 1 ('all'). 
%                           The choice for the latter is based on the fact that 
%                           the last position is the same as the first position,
%                           so the last position is removed. (Default = 'all')
%      choose_emitter_receiver - the approach for choosing the
%                           emitter-receiver pairs. This can be set
%                           'open_angle', or 'distances'. Both approaches
%                           are equivalent. (Default: 'open_angle')
%      open_angle           - the maximum permissible open angle for
%                            including the emitter-receiver pairs in the
%                            image reconstruction
%      num_us_segment       - the object-in-water data is divided to
%                            to num_us_segment segments, and  the time-of-flights
%                            are computed segment by segment
%      postprocess_tofs     - Boolean controlling whether the computed TOfs
%                            are postprocessed or not. (default : true)
%      tofs_fac             - a factor between zero and one for the min/max 
%                            tolerance for the outliers in the tof picking.
%                            For data with lower quality, for example because of
%                            movements, smaller tolerance (smaller value) must
%                            be used. ( Default : 3.0)
%      smooth_cone_boundary - Boolean controlling whether the boundary of the cone
%                             for accepting the emitter-receiver pairs are smoothed
%                             or not. Setting this parameter true will considerably
%                             reduce the artefact in water, and reduces erroneous
%                             refractions. (default : true)
%      manual_remove_bad_transducers - Boolean controlling whether the TOFs 
%                            from bad transducers are removed or not. If this 
%                            is set true, the bad transducers are removed from 
%                            the difference TOF map, and are set zero.
%                            The bad transducers found by eye inspection from
%                            the difference TOF map are consistent with those in 
%                            the pdf file provided by the PA imaging company. (default = true)
%      sound_speed_ranges_diff - An assumption of a min/max homogeneous sound speed
%                            inside the bowl minus the only water soud speed.
%                            (Default = [-50, 50] for the object-in-water data) 
%      absorption_frequencies - the vector of frequencies at which the
%                             relative amplitudes are computed. If it is set zero,
%                             the maximum relative amplitude is computed. If it set
%                             nan, the relative amplitudes are not computed. 
%         
%
% OUTPUTS:
%      tof_data              - the num_receiver x num_emitter matrix of the
%                              difference time-of-flight data [s]
%      weight_ray            - the coefficient enforced on the
%                              computed time-of-flight for each
%                              emitter-receiver pair
%      tof_emitters_receivers - the vector of time-of-flights [s] for the water 
%                               encompassing the object and computed using the
%                               given temperatue of water during the rotation 
%                               of the bowl.
%    
%    
%
% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 10.08.2022
%
% This script is part of the r-Wave Tool-box (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian and Ben Cox


para.num_worker_pool = 16;
para.include_only_water = true;
para.do_calculate_tofs = true;
para.choose_last_us_position = 'all';
para.choose_emitter_receiver = 'open_angle';
para.open_angle = pi/4;
para.num_us_segment = 100;
para.postprocess_tofs = true;
para.tofs_fac = 3;
para.smooth_cone_boundary = false;
para.manual_remove_bad_transducers = true;
para.sound_speed_ranges_diff = 50 *[-1, +1];
para.save_results = true;
para.absorption_frequencies = nan;


% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end


if ~para.do_calculate_tofs
    para.save_results = false;
end


if para.postprocess_tofs
    
    % Choose a factor between zero and one for the min/max tolerance for the outliers
    % in the tof picking. For bad data, for example because of movements, smaller
    % tolerance (smaller value) could be used. (Default = 0.8)
    para.tofs_fac = 3;
    
    if para.tofs_fac < 0.5
        
        % Give a warning, if the tofs_fac is chosen smaller than 0.5.
        warning(['A small tolerance for accepting the TOFs may affect the true TOFs.']);
    end
    
else
    
    % this factor is not used when postprocess_tofs is set false.
    para.tofs_fac = nan;
    
end
 
% get the number of receivers (transducers)
num_receiver = size(receiver.positions{1}, 2);

% get the number of emitters
num_emitter = size(emitter.positions, 2);


switch para.choose_emitter_receiver
    
    case 'distances'
        
        % calculate the distance between emitter-receiver pairs
        distance_emitters_receivers = vectorise(calculateDistanceEmitterReceiver(emitter.positions,...
            receiver.positions, rotation_indices));
        
        % reshape the resulting vector to a matrix of size num_receiver x
        % num_emitter
        distance_emitters_receivers_matrix = reshape(distance_emitters_receivers,...
            [num_receiver, num_emitter]);
        
        % get a vector of binaries for the emitters and receivers included
        % in the image reconstruction
        binaries_emitter_receiver = distance_emitters_receivers_matrix > para.minimum_distance;
        
        % getthe weights for all the rays 1
        weight_ray = one(num_receiver, num_emitter);
        
    case 'open_angle'
        
        % calculate the angles between a geomterical vector normal to surface of emitters
        % (from emitters to the centre of the bowl) and a geomterical vector from
        % emitters to the receivers
        angles_emitters_receivers = calcAngleEmitterReceiver(emitter.positions,...
            receiver.positions, rotation_indices);
        
        % choose the emitters and receivers for which the system matrix is
        % computed
        binaries_emitter_receiver = angles_emitters_receivers < para.open_angle;
        
        if para.smooth_cone_boundary
            
        % smooth_cone_boundary reduces the artefact in water by enforcing weights
        % smaller than 1 to the rays travelling only in water, and also may slightly
        % affect the image for large cups, so it is not recomended for
        % large cups, say (6, 7 , 8).
            
        % get the weights as a function of angles_emitter_receivers.  The weights
        % will be used for smoothing the interface between the volumes for
        % which the rays are included in the image reconstruction, and volumes for
        % which the rays are neglected. The rays form a cone with
        % vetex the emitter and the axis a vector connecting the emitter to the centre
        % of the bowl. All the rays intercepted by receivers inside the cone are
        % included in the image reconstruction. However, the discontinuity in the
        % boundary of this cone will appear as a large artefact. To remove this
        % continuity, the boundary of the cone is smoothed using a tukey window
        % as function of angles_emitters_receivers
        
        % make a tukey window for equidistant angles. (One index is added for zero angle,
        % which is aligned by the axis of the cone.) The coefficient are specified to the rays
        % such that the contribution of the rays close to the boundary of the cone
        % will gradually decay, and thus the discontinuity on the cone's boundary
        % is removed.
        
        % get the angles 
        angles_sampled = (0:round(para.open_angle*180/pi)) * pi/180;
        
        % apply the tukey window on angles
        coeff_angles = getWin(length(angles_sampled), 'Tukey', 'Param', 0.25); 
        
        % remove the filter from the small angles (the rays close to the axis of
        % the cone)
        ix = find(coeff_angles == 1, 1, 'first');
        coeff_angles(1:ix) = 1;
        
        % get the interpolation operator for the asending angles and their
        % specified coefficients.
        interpolant_operator = griddedInterpolant(angles_sampled, coeff_angles);

        % allocate a zero matrix for ray coefficients
        weight_ray = zeros(num_receiver, num_emitter);
        
        for ind_emitter = 1:num_emitter
            
            % get the ascending angles between the geometerical vectors connecting
            % the emitter and all receivers (transducers) and the axis of the cone
            % associated with the emitter
            [angles_emitter_ascend, ~] = ...
                sort(angles_emitters_receivers(:, ind_emitter));
            
            % get the coefficients for the ascending angles
            coeffs_interpolated_emitter = interpolant_operator(angles_emitter_ascend);
            
            % get the location of the original angles in the sorted angles
            [~, loc] = ismember(angles_emitters_receivers(:, ind_emitter), angles_emitter_ascend);
            
            % get the interpolated coefficients for the current emitter
            weight_ray(:, ind_emitter)= coeffs_interpolated_emitter(loc);
            
        end
        
        else
        
        % make all the weighting coefficients 1
         weight_ray = ones(num_receiver, num_emitter);
         
        end
end


% display the percentage of the rays included in the image reconstruction
disp(['The percentage of the used emitter-receiver pairs for image reconstruction is:'...
    num2str(nnz(binaries_emitter_receiver)/(num_emitter * num_receiver)*100) '%'])

%% ========================================================================
% CALCULATE OR LOAD THE TOFS
% =========================================================================

    % compute the distances between emitter-receiver pairs
    distance_emitters_receivers = vectorise(calculateDistanceEmitterReceiver(emitter.positions,...
        receiver.positions, rotation_indices));
    
    % compute the true TOFs between emitters and receivers via including the
    % sound speed of water encompassing the object
    tof_emitters_receivers = 1/sound_speed_only_water * ...
        distance_emitters_receivers .* refractive_background_object_alldata;

    
if  para.do_calculate_tofs  
    
    % get the time array
    time_array = cumsum([0, (1/time_info.us_frequency) * ones(1, size(data_water, 2) - 1)]);
    
    % get the us shot time based on the pdf provided by the PA imaging company. 
    % The shot time is fixed for all excitatons, and agrees well with the 
    % calculated TOF differences.
    emitter.shot_time = (time_info.shot_time_index - 1) * (1/time_info.emit_frequency);
    
    % get an approximate pulse duration.
    % This is not used when the AIC approach is used for time-of-flight
    % picking.
    emitter.pulse_duration = 3e-6;
    
   
    % get the optional parameters for loading the data and time-of-flight
    % picking
    data_args = {'nWorkerPool', para.num_worker_pool, 'Method', 'Modified_AIC',...
        'SoundSpeedRangesDiff', para.sound_speed_ranges_diff, 'SoundSpeedRef_water',...
        sound_speed_only_water, 'SoundSpeedRef_object', sound_speed_object_water,...
        'minimum_distance', para.minimum_distance, 'binaries_emitter_receiver',...
        para.choose_emitter_receiver, 'open_angle', para.open_angle,...
        'Cutoff_freq', nan, 'absorption_frequencies', para.absorption_frequencies,...
        'choose_last_us_position', para.choose_last_us_position, 'save_data', para.save_results};
    
    % load the ultrasound data part by part, and compute the
    % time-of-flights
    [tof_het, tof_hom, amplitude_het, amplitude_hom, ~] = ...
        loadPammothDataCalculateTofs(...
        emitter, receiver, rotation_indices, us_shot_index, emit_shot_index, time_array,...
        data_water, para.num_us_segment, data_paths_object, user_data_path, data_args{:});
    
    if strcmp(para.choose_last_us_position, 'all')  &&    para.save_results
            
            % get the emitter shot time [s]
            emitter_shot_time = emitter.shot_time;

            % save the time-of-flights and amplitudes for the object-in-water and only-water data
            save([user_data_path 'tofs.mat'], 'tof_hom', 'tof_het', 'emitter_shot_time',...
                'sound_speed_object_correction_coeff', 'tof_emitters_receivers', '-v7.3');
            
    end
    
else
    
    if strcmp(para.choose_last_us_position, 'all')
        
        % load the time-of-flights and amplitudes for the object-in-water and only-water data
        load([user_data_path 'tofs.mat'], 'tof_hom', 'tof_het', 'emitter_shot_time',...
                'sound_speed_object_correction_coeff', 'tof_emitters_receivers');
        
        % add the shot time to the emitter struct
        emitter.shot_time = emitter_shot_time;
        
    else
        
        % give an error if only the early US angular positions are included.  
        error(['If only the early angular positions are included, the saved time-of-flights'...
        'are not available.'])
    
    end
    
end

% compute the time-of-flights for the object-in-water data by subtracting
% the shot time from first-arrival times
tof_het = tof_het - emitter.shot_time;

% calculate the discrepancy of time-of-flights between the object-in-water data
% and the only-water data
if para.include_only_water
    
    % compute the time-of-flights for the only-water data by subtracting
    % the shot time from first-arrival times and then correct for the
    % relative difference with the water encompassing the object during the
    % rotation of the bowl
    tof_hom = sound_speed_object_correction_coeff .* ...
        (tof_hom - emitter.shot_time);
    
    % compute the difference TOFs, the discrepancy of the picked TOFs between
    % the object-in-water and only-water data
    tof_data = tof_het - tof_hom;
    
else
    
    % calculate the discrepancy of the picked TOFs between object in water
    % and the true TOfs
    tof_data = tof_het - tof_emitters_receivers;
end


if para.manual_remove_bad_transducers
    
    % get the bad transducers based on Laurens' PDF file
    bad_transducers_indices_sub = [1,  1,  3, 3, 4, 6, 12, 9;...
        22, 23, 30,31,31,1, 23, 6]' + 1;
    bad_transducers_indices = sub2ind([32 16], bad_transducers_indices_sub(:,2),...
        bad_transducers_indices_sub(:,1) );
    
    % get the ugly transducers based on Lauren's PDF file
    ugly_transducers_indices_sub = [13, 12,  8, 6, 10;...
        28, 31, 31, 3, 1]' + 1;
    
    ugly_transducers_indices = sub2ind([32 16], ugly_transducers_indices_sub(:,2),...
        ugly_transducers_indices_sub(:,1) );
    
    % get the remaining bad transducers by eye inspection from the difference TOF data
    badtof_tranducer_indices =  [18; 123; 186; 231; 351; 407; 414];   %   [18, 123, 186, 351, 392, 414]';  % [123; 351];   % 
    
    % get all the bad indices (bad and ugly)
    bad_transducers_indices = [bad_transducers_indices; ugly_transducers_indices;...
         badtof_tranducer_indices];
    
    % remove the bad and ugly transducers from the TOFs
    [bad_transducer_binaries, ~] = ismember(emit_firing_transducer_index + 1,...
       bad_transducers_indices);
    
    % set the TOF for bad receivers zero
    tof_data(bad_transducers_indices, :) = 0;
    
    % set the TOF for bad emitters zero
    tof_data(:, bad_transducer_binaries) = 0;
    
end

% remove the unwanted TOF differences using the user defined binaries
% the noninfinite TOFs correspond to emitter-receiver pairs for which the
% same tranducer is used for both emission and reception
tof_data(~binaries_emitter_receiver | ~isfinite(tof_data)) = 0;


% get the logarithmic relative amplitudes for the chosen frequencies
if all(isfinite(para.absorption_frequencies))
    
    amplitude_data = cell(length(para.absorption_frequencies), 1);
    for ind_freq = 1: length(para.absorption_frequencies)
        
        % Calculate the logarithm of the only-water amplitudes to the object-in-water
        % amplitudes
        amplitude_data{ind_freq} = 0.2 * log10(amplitude_hom{ind_freq} ./...
            amplitude_het{ind_freq});

    end
    
end

disp(['The percentage of bad pickings for nonzero tofs are:'...
    num2str(nnz(abs(tof_data(:)) > 4e-6)/nnz(abs(tof_data(:))>0)*100)]);

if para.postprocess_tofs
    % choose the neighboring transducers for each transducer using a
    % triangulation
    [connected_receiver_indices, ~] = findConnectedPoints(receiver.positions{1});
    
    % the mean of the nonzero difference tofs
    tofs_mean = mean(tof_data(abs(tof_data)>0));
    
    % the standard deviation of the nonzero difference tofs
    tofs_std = std(tof_data(abs(tof_data)>0));
    
    % minimum/maximum tolerance
    tols = [tofs_mean - para.tofs_fac * tofs_std, tofs_mean + para.tofs_fac * tofs_std];
    
    % postprocess the calculated difference tofs for each emitter aligning the receivers
    % with nonzeros tof difference
    tof_data = postProcessTOFsAlignedReceivers(tof_data, ...
        connected_receiver_indices, tols);
    
    % postprocess the difference tofs again for each receiver aligning the
    % emitters with nonzeros tof difference
    tof_data = postProcessTOFsAlignedEmitters(tof_data,...
        connected_receiver_indices, emit_firing_transducer_index, rotation_indices,...
        tols);

end

% display the TOF matrix for the first full rotation of the excitation over
% all the meitters
figure; imagesc(tof_data);colorbar;

% get the vector of TOF differences as the input for image reconstruction
tof_data = tof_data(:);

end