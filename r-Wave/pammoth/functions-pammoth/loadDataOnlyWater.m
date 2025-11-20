function [data_water, emit_shot_index_water, us_shot_index_water,...
    sound_speed_only_water] = loadDataOnlyWater(data_paths_object,...
    only_water_measurement_date, num_receiver)
%LOADDATAONLYWATER LOADS THE ONLY-WATER DATA
%
%
% DESCRIPTION:
%             loadDataOnlyWater loads the only-water data
%
%
% USAGE:
%
%
% INPUTS:
%       data_paths_object - a struct array containing the paths associated
%                           with the object-in-water data
%       only_water_measurement_date - the date on which the only-water data
%                           are measured.
%       num_receiver      - the number of receivers or transducers
%
%
% OPTIONAL INPUTS:
%
%
%
%
% OUTPUTS:
%      data_water         - the num_transduccer * num_time * num_transducer
%                           data associated with one rotation
%      emit_shot_index_water - the vector of emit_shots for the only-water
%                           data. (See the pdf file provided by the PA
%                           imaging company.)
%      us_shot_index_water - the vector of us_shots for the only-water
%                           data. (See the pdf file provided by the PA
%                           imaging company.)
%      sound_speed_only_water - a scalar value for the sound speed of the
%                           only water
%
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

switch only_water_measurement_date
    case '05-2021'
        
        % Get the path for the only-water data 
        data_paths_water.directory = '20210528_primer_elevated_water_level/';
        data_paths_water.data_name = '20210528_00.hdf5';
        
        % get the data in water from the corresponding HDF5 file
        data_water = h5read([data_paths_object.main_directory,...
            data_paths_water.directory, data_paths_water.data_name], '/us_data');
        
        % convert the data to the standard format with dimensions
        % num_receiver x num_time x num_emitter
        data_water = permute(double(data_water), [2, 1 ,3]);
        
        % get the only-water data for single excitation from all the
        % transducers.
        % This will then be extended to all excitations.
        data_water = data_water(:, :, 1: num_receiver);

    otherwise
        
        % set an empty variable
        data_water = [];
end


% get the emit shot index for only-water data from the corresponding HDF5 file
% emit_shot_index_water = h5read([data_paths_object.main_directory,...
%    data_paths_water.directory, data_paths_water.data_name], '/emit_shot_index');

% get the us shot index for only-water data from the corresponding HDF5 file
% us_shot_index_water = h5read([data_paths_object.main_directory,...
%    data_paths_water.directory, data_paths_water.data_name], '/us_shot_index');
% These are available in the object-in-water data, and are redundant.
emit_shot_index_water = [];
us_shot_index_water = [];

% get the table for temperature for the only-water data (reference)
info_table_water = importdata([data_paths_object.main_directory,...
    'tables/water_aux_tbl.txt']);


% get the sound speed [m/s] for the only-water data (reference sound speed)
% from the the mean of the temperature of the stationary data for
% only-water data
sound_speed_only_water = waterSoundSpeed(1/100* mean(info_table_water.data(1:2)));


end