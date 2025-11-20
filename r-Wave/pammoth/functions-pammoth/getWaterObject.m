function [sound_speed_object_correction_coeff,...
    refractive_background_object_alldata,...
    sound_speed_object_water] = getWaterObject(water_object_temperatures, rotation_indices,...
    sound_speed_only_water, correct_temperatures, num_emitter, num_receiver)
%GETWATEROBJECT gets the information about the water encompassing the object
%   
%
% DESCRIPTION:
%             getWaterObject gets the information about water encompassing
%             the object during the rotation of the bowl
%      
%
% USAGE:
%     
%
% INPUTS:
%        water_object_temperatures - the last_us_position x 2 marix array of 
%                           temperatures of the water encompassing the object
%                           for the object-in-water data during the rotation
%                           of the bowl.
%        rotation_indices - the vector including the indices of angular
%                            positions for exitations
%        sound_speed_only_water - the scalar value for the sound speed
%                           [m/s] for only-water data
%        correct_temperatures - Boolean controlling whether the changes of
%                            the sound speed in water encompassing the object
%                            during the rotation of the bowl are accounted
%                            for or not.
%        num_emitter       - the number of emitters
%        num_receiver      - the number of receivers (transducers)
%  
%       
% OPTIONAL INPUTS:
%        sound_speed_object_correction_coeff - the num_emitter x 1 vector of 
%        sound speed correction coefficient for all excitations during the rotation
%        of the bowl.
%        refractive_background_object_alldata - the (num_emitter *
%        num_receiver) * 1 vector of the refractive index of water
%        encomapassing the object for all exciations/receivers (all time series)
%        sound_speed_object_water - A scalar value representing the mean sound
%                                  speed of water encompassing the object
%                          
%      
%
% OUTPUTS:
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
%% ========================================================================
% GET THE SOUND SPEED IN WATER
%==========================================================================
% Get the sound speed in water encompassing the object in rotations
sound_speed_background_object = waterSoundSpeed(1/100 * mean(water_object_temperatures, 2));

% get the mean sound speed of the water encomapssing the object.
% This is used for the sound speed outside the chosen Region-Of-Interest (ROI)
% for image reconstruction
sound_speed_object_water = mean(sound_speed_background_object);

if correct_temperatures
    
    % get the sound speed of water encompassing the object for all excitations during the rotation
    sound_speed_background_object_allexcitations = sound_speed_background_object(rotation_indices)';
    
    % get the sound speed correction coefficient for all excitations during the rotation
    % of the bowl. This is used for making the only-water data consistent with the
    % object-in-water data for the purpose of calculating the matrix of difference
    % time-of-flights (TOFs)
    sound_speed_object_correction_coeff = sound_speed_only_water./...
        sound_speed_background_object_allexcitations;
    
    % get the refractive index of water encomapassing the object for all
    % exciations/receivers (all time series) during the rotation of the bowl
    refractive_background_object_alldata = vectorise(...
        repmat(sound_speed_object_correction_coeff, [num_receiver, 1]));
    
else
    
    % the background sound speed (the sound speed of water encompassing the object)
    % for all excitations during the rotation
    % sound_speed_background_object_allexcitations = sound_speed_only_water* ones(1, num_emitter);
    
    % get the sound speed correction coefficient for all excitations during the rotation
    % of the bowl.
    sound_speed_object_correction_coeff = ones(1, num_emitter);
    
    % get the refractive index of water encompassing the object for all
    % exciations/receivers (all time series) during the rotation of the bowl
    refractive_background_object_alldata = ones(num_emitter * num_receiver, 1);
    
end


        
end