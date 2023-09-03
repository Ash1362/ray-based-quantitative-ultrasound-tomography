function [data_paths_object, table_array_position, water_object_temperatures,...
    cup_size] = setPathGetInfos(data_ids, last_angular_position)
%SETPATHGETINFOS gets the path for data, and upload the table
%
%
% DESCRIPTION:
%             setPathGetInfos gives the paths associated with the data, and
%             upload the information included in the table in the associated path
%
%
% USAGE:
%
%
% INPUTS:
%       data_ids          - a struct array containing the IDs (and information)
%                           about the data. This struct includes the fields:
%       'object'            - the type of object for which the ultrasound data is measured.
%                           This can be 'volunteers', 'patients', 'us_measurements',
%                           'breast_phantom', 'channel_phantom', or 'water_phantom'
%       'main_path'         - the main path for all the stored ultrasound
%                           data and information. It can be adjusted in the script
%                           'startup_pammoth.m'
%       'id'                the ID of the data
%       'date'              the date on which the data were measured
%       'measurement_number'- the measurement number as a string
%       'cup_size'          the cup size, which can be integers 2,3,...,8.
%       'main_path_user_data' - the main path for storing the results, which can
%                            be set by the user. It is recommended that the
%                            default path given in the example is used.
%       last_angular_position - the last included angular position

%
%
% OPTIONAL INPUTS:
%
%
%

%
%
% OUTPUTS:
%      data_paths_object  - the paths associated with the data
%      table_array_position - the table array associated with the angular
%                             position of the transducers
%      water_object_temperatures - the last_us_position x 2 marix array of
%                           temperatures of the water encompassing the object
%                           for the object-in-water data during the rotation
%                           of the bowl.
%      cup_size           - the cup size, which can be 2,3,...,8.
%
%
% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 10.08.2022
%
% This script is part of the r-Wave Tool-box (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian and Ben Cox


% get the object for which the ultarsound data were measured
object = data_ids.object;

% get the ID of the data
id = data_ids.id;

% get the date on which the date were measured
date = data_ids.date;

% get the cup size
cup_size = data_ids.cup_size;

% get the measurement number
measurement_number = data_ids.measurement_number;

% get the main path for the data (should be fixed for all data)
data_paths_object.main_directory = data_ids.main_path;

% get the main path for the data (should be fixed for all data)
data_paths_object.table_directory = [data_ids.main_path_user_data, 'tables/'];

%% ========================================================================
% GET THE PATH FOR THE DATA FROM OBJECT AND THEIR CORRESPONDING INFORMATION
%==========================================================================
% Get the path for the data from object
% overwrite the already known information
switch object
    case 'US_phantom_PAA_filled'
        
        % These data sets are deprecated.
        Error('These data sets are deprecated.')
        
    case 'volunteers'
        
        
        % get the data path for the object
        data_paths_object.directory = [date, id];
        
        switch id(end-3:end-1)
            
            case {'027'}
                
                % get the name of the data for the object
                data_paths_object.data_name = [date, measurement_number,...
                    '_averaged'];
                
                % get the table for the data for object
                info_table = importdata([data_paths_object.table_directory,...
                    data_paths_object.directory, date, measurement_number...
                    '_aux_tbl' '.txt']);
                
            case '086'
                
                % get the name of the data for the object
                data_paths_object.data_name = [date, measurement_number];
                
                % get the table for the data for the object
                info_table = importdata([data_paths_object.table_directory,...
                    data_paths_object.directory, 'Aux_', data_paths_object.data_name...
                    '.txt']);
                
            otherwise
                
                %      case {'016', '049', '058'}
                 
                % get the name of the data for the object
                data_paths_object.data_name = [date, measurement_number];
                
                % display the path for table
                disp([data_paths_object.table_directory,...
                    data_paths_object.directory, data_paths_object.data_name,...
                    '_aux_tbl' '.txt'])
                
                % get the name of the data for the object
                data_paths_object.data_name = [date, measurement_number];
                
                % get the table for the data for object
                info_table = importdata([data_paths_object.table_directory,...
                    data_paths_object.directory, data_paths_object.data_name,...
                    '_aux_tbl' '.txt']);
        end
        
        
    case 'us_measurements'
        
        % Get the data path for the object
        data_paths_object.directory = [date, id];
        
        % Get the data name for the object
        data_paths_object.data_name = [date, measurement_number...
            '_US_only'];
        
        % display the path for table
        disp([data_paths_object.table_directory,...
            data_paths_object.directory, 'Aux_', date, measurement_number...
            '.txt'])
        
        % get the table for the data for the object
        info_table = importdata([data_paths_object.table_directory,...
            data_paths_object.directory, 'Aux_', date, measurement_number...
            '.txt']);
        
    case {'patients', 'breast_phantom', 'water_phantom', 'phantoms'}
        
        % get the data path for the object
        data_paths_object.directory = [date, id];
        
        % get the name of the data for the object
        data_paths_object.data_name = [date, measurement_number...
            ];
        
        % display the path for table
        disp([data_paths_object.table_directory,...
            data_paths_object.directory, data_paths_object.data_name,...
            '_aux_tbl' '.txt'])
        
        % get the table for the data for object
        info_table = importdata([data_paths_object.table_directory,...
            data_paths_object.directory, data_paths_object.data_name,...
            '_aux_tbl' '.txt']);
        
end

% get the table for the angular position of us shots, and convert it to an
% array
info_table_textdata = info_table.textdata(2:end-2, :);
table_array_position = zeros(size(info_table_textdata));

for i = 1: size(table_array_position, 1)
    for j = 1: size(table_array_position, 2)
        table_array_position(i, j) = str2double(cell2mat(info_table_textdata(i, j)));
    end
end

% the last angular position have been chosen based on the assumption of
% 101 total angular positions, if total angular positions are more than
% 101, and the approach for choosing the last angular position is not
% set 'early', the last angular position must be updated.
if size(table_array_position, 1) > 101 && last_angular_position > 10
    last_angular_position = size(table_array_position, 1) - 1;
end

% reduce the size of the angular positions based on the last position
table_array_position = table_array_position(1: last_angular_position, :);

% get the temperature of water for us shots during the measurements
water_object_temperatures = info_table.data(1: last_angular_position, 1:2);


end