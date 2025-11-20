function [emitter, receiver, emit_shot_index, us_shot_index,...
    emit_firing_transducer_index, rotation_indices] = getPositions(data_paths_object,...
    table_array_position, radius_bowl, user_data_paths, varargin)
%LOADPAMMOTHUSTDATA loads the UST portion of the data from the pammoth system,
% and computes the time-of-flights (TOFs)
%   
%
% DESCRIPTION:
%             loadProcessPammothUstTofData loads the ultrasound portion of
%             the Pammoth data part by part, and computes the TOFs of the
%             ultrasound data using the AIC first-arrival picking algorithm.
%      
%
% USAGE:
%     
%
% INPUTS:

%       data_paths_object - a struct array containing the paths associated
%                           with the object-in-water data
%       table_array_position - the table array associated with the angular
%                           position of the transducers
%       radius_bowl       - the radius [m] of the bowl
%       user_data_paths    - a struct containing the path for storing the information
%                           associated with the UST data, together with the 
%                           computed TOFs. These are saved in the user_data_paths. 
%                           Therefore, the user can reconstruct the image without
%                           the need for the main data by setting do_calculate_tofs
%                           false, if the associated information and computed
%                           time-of-flights are stored on the path.
%  
%       
% OPTIONAL INPUTS:
%      do_calculate_tofs  - Boolean controlling whether the time-of-flights
%                          (TOFs) are calculated, or the already calculated
%                           TOFs are loaded. This Boolean will be automatically
%                           set true, if the TOFs do not exist in the
%                           defined path.
%     'choose_last_us_position' - A parameter determining how to choose the last
%                                 US position, which can be 'early' (chosen by the user),
%                                 or automatically set all the angular US positions - 1 ('all'). 
%                                 The choice for the latter is based on the fact that 
%                                 the last position is the same as the first position,
%                                 so the last position is removed. (Default = 'all')
%
% OUTPUTS:
%    emitter               - a struct with field
%    'positions'           - the position of emitters
%    receiver              - a struct with field
%    'positions'           - the position of receivers
%    emitter_shot_index    - the vector including the indices of excitions
%    us_shot_index         - the vector including the indices of US shots 
%    emit_firing_transducer_index - the zero-based (0,1,..., 511) transducer
%                            indices corresponding to the emitters, i.e.,
%                            the index of transducers in the excitation mode
%    rotation_indices      - the vector including the indices of angular
%                            positions for exitations
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

para.do_calculate_tofs = true;
para.choose_last_us_position = 'all';

% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

%%=========================================================================
% SET DIRECTORIES
%==========================================================================
% convert the user data path from struct to a string
user_data_path = [user_data_paths.main_directory, user_data_paths.directory,...
    user_data_paths.data_name '/'];

% make the user (local) directory, if it does not exist
makeDirectory(user_data_path);

if strcmp(para.choose_last_us_position, 'early')
   
    % make the user data path for data for early iterations
    % for saving in the directory, if requested.
    user_data_path = [user_data_path, 'early_'];
end

%% ========================================================================
% DEFINE THE SETTINGS FOR THE PAMMOTH DATA
%==========================================================================
% Get the number of the US shot indices for each measurement set. The indices are
% zero-based, i.e., 0,...,569 (fixed for all PAMMOTH measurements, and must
% not be changed.)
% num_us_shot_per_measurement = 570;

% get the number of the angular positions
num_us_position = size(table_array_position, 1);

% get the number of US shot indices per each position
% (equal to the first shot index for the second position (row))
% each row of array 'table_array_position' corresponds to an angular position
num_us_shot_per_position = table_array_position(2, 1);

% display the path for data
disp([data_paths_object.main_directory,...
        data_paths_object.directory, data_paths_object.data_name, '.hdf5'])
    
% get the shot_counter vector matching the 'us_data', which contains 1, 2, ..., 9, 11,...
% a vector with size num_us_shots x 1, i. e. 18090 x 1
if  para.do_calculate_tofs
    
    % get the us shot indics for the object-in-water data
    us_shot_index = double(h5read([data_paths_object.main_directory,...
        data_paths_object.directory, data_paths_object.data_name, '.hdf5'], '/us_shot_index'));
    
    % get the emit shot indices (the us shot indices without nans)
    emit_shot_index = double(h5read([data_paths_object.main_directory,...
        data_paths_object.directory, data_paths_object.data_name, '.hdf5'], '/emit_shot_index'));
    
    % get the zero-based (0,1,..., 511) transducer indices corresponding to the
    % emitters, ie. the index of transducers in the excitation mode
    emit_firing_transducer_index = h5read([data_paths_object.main_directory,...
        data_paths_object.directory, data_paths_object.data_name, '.hdf5'],...
        '/emit_firing_transducer_index');
    
    % save the US indices on the user's private path
    % save([user_data_path 'us_indices.mat'], 'us_shot_index', 'emit_shot_index',...
    %  'emit_firing_transducer_index', '-v7.3');
    if ~exist([user_data_path 'us_data.h5' ])
        % create the HDF5 paths
        h5create([user_data_path, 'us_data.h5'], '/us_shot_index', size(us_shot_index))
        h5create([user_data_path, 'us_data.h5'], '/emit_shot_index', size(emit_shot_index))
        h5create([user_data_path, 'us_data.h5'], '/emit_firing_transducer_index', size(emit_firing_transducer_index))
    end
    
    % write the HDF5 file on the defined user data path
    h5write([user_data_path, 'us_data.h5'], '/us_shot_index', us_shot_index )
    h5write([user_data_path, 'us_data.h5'], '/emit_shot_index', emit_shot_index)
    h5write([user_data_path, 'us_data.h5'], '/emit_firing_transducer_index', emit_firing_transducer_index)
    
else
    
    
    % load the US indices from the user's private path
    % load([user_data_path 'us_indices.mat'], 'us_shot_index', 'emit_shot_index',...
    %   'emit_firing_transducer_index');
    
    % read the HDF5 file for getting the information in data
    us_shot_index = h5read([user_data_path, 'us_data.h5'], '/us_shot_index');
    emit_shot_index = h5read([user_data_path, 'us_data.h5'], '/emit_shot_index');
    emit_firing_transducer_index = h5read([user_data_path, 'us_data.h5'], '/emit_firing_transducer_index');
    
    
end

    
% get the number of US shots
% num_us_shot = us_shot_index(end) + 1;

% get the nan binary in the US shot index
[notnans_binary, ~] = ismember(us_shot_index, emit_shot_index);

% find the indices with nans
nans_index = find(~notnans_binary);

% get the index of each angular position for each excitation in the emit_shot_index
% vector. The rotation indices are 1-based ({1,...,num_us_position})
rotation_indices = floor(emit_shot_index / num_us_shot_per_position) + 1;

% Get the rotation indices for each excitation in the us_shot_index
rotation_indices_us = ones(size(us_shot_index));
rotation_indices_us(nans_index) = 0;
rotation_indices_us(rotation_indices_us > 0) = rotation_indices;

% determine the last index in the emit_shot_index vector for the last angular position index.
num_excitation = find(rotation_indices == num_us_position, 1, 'last');

% replaced by last_angular_position

% determine the last index in the us_shot_index vector for the last angular position index.
last_rotation_index_us = find(rotation_indices_us == num_us_position, 1, 'last');
% replaced by last_angular_position

% reduce the vector of rotation_indices based on the chosen last ultrasound position
rotation_indices = rotation_indices(1 : num_excitation);

% reduce the vector of emit shot indices based on the chosen last ultrasound position
emit_shot_index = emit_shot_index(1 : num_excitation);

% reduce the vector of zero-based transducer indices corresponding to the
% emitters, ie. the index of transducers in the excitation mode
emit_firing_transducer_index = emit_firing_transducer_index(1 : num_excitation);

% reduce the vector of us shot indices based on the chosen last US
% position.
us_shot_index = us_shot_index(1 : last_rotation_index_us + 1);

if rotation_indices_us(last_rotation_index_us + 1)
    
    % remove the last element from us_shot_index, if it is nonzero, ie. if it
    % does not corresponds to nans
    us_shot_index(end) = [];
end


% get the nan binary in the US shot index 
[notnans_binary_reduced, ~] = ismember(us_shot_index, emit_shot_index);

% find the indices with nans
nans_index_reduced = find(~notnans_binary_reduced);

% get the US shot indices without nans (for checking)
us_shot_index_notnan_reduced = us_shot_index;
us_shot_index_notnan_reduced(nans_index_reduced) = [];

% check if us_shot_index_notnan is the same as the given
% emit_shot_index, and display a confirmation message
if norm(emit_shot_index - us_shot_index_notnan_reduced) == 0
    disp('The data indexing has been done correctly so far.');
else
    error('The data reading process has some issues.');
end

% Get the number of emitters
%num_emitter = size(emit_shot_index, 1);

%%=========================================================================
% LOAD THE TABLE FOR POSITIONS
%==========================================================================
% get the table for positions
transducer_positions_txt = table2array(readtable([user_data_paths.main_directory, 'tables/'...
    'transducer_channel_and_position_data.txt'], 'ReadVariableNames', false));

% convert the table to array
table_transducer_positions_fixed = convertTable2Array(transducer_positions_txt(2:end, :) );

% get the zero-based orders
element_number = 32 * table_transducer_positions_fixed(:, 10) ...
    + table_transducer_positions_fixed(:, 11);

% order the rows in the table using the 1-based element numbers
table_transducer_positions_fixed(element_number + 1, :) = table_transducer_positions_fixed;

%% ========================================================================
% GET THE POSITION OF RECEIVERS
%==========================================================================
% modify the phi positions according to the pdf provided by PA imaging
% company

% get and modify the azimuthal angle
phi_rad_fixed = - table_transducer_positions_fixed(:, 4) + 2.259;

% get the polar angle
theta_rad_fixed = table_transducer_positions_fixed(:, 3);

% get the increment in the azimuthal angle (phi)
phi_inc = table_array_position(3, 2) - table_array_position(2, 2);

% get the Cartesian position of the transducers in a fixed setting (first
% rotational angle)
position_fixed = radius_bowl * [cos(phi_rad_fixed) .* sin(theta_rad_fixed), ...
    sin(phi_rad_fixed) .* sin(theta_rad_fixed), cos(theta_rad_fixed)]';

% Get the number of receivers
num_receiver = size(position_fixed, 2);

% Allocate vectors for the azimuthal and polar angles of the
% receivers during the rotation
theta_rad = vectorise(repmat(theta_rad_fixed, 1, num_us_position));
phi_rad = vectorise(repmat(phi_rad_fixed, 1, num_us_position));

% desired rotation angle increment in [degrees]
phi_inc_rad = (phi_inc / 180) * pi;

% get the rotation angles
rotation_angles  = (0:num_us_position-1) * phi_inc_rad;

% get the azimuthal component of the angular positions
phi_rad = phi_rad + vectorise(repmat(rotation_angles, num_receiver, 1));

% get the position of receivers
receiver_positions = radius_bowl * [cos(phi_rad) .* sin(theta_rad),...
    sin(phi_rad) .* sin(theta_rad), cos(theta_rad)]';

% get the Cartesian position of the transducers in a rotational setting (first
% rotational angle)
receiver.positions = mat2cell(receiver_positions, size(receiver_positions, 1), ...
    num_receiver * ones(1, num_us_position));

%%=========================================================================
% GET THE CARTESIAN POSITION OF EMITTERS
% =========================================================================
% get the fixed position of emitters
phi_fixed_emitters = phi_rad_fixed(emit_firing_transducer_index + 1);

% get the azimuthal angle for all US shots in the rotation
% This is a matrix of size num_us_shot_per_position x num_us_position
phi_rotate_emitters = repmat(rotation_angles, [num_us_shot_per_position, 1]);

% sort the azimuthal angles for all the excitations (emitters)
phi_rotate_emitters = phi_rotate_emitters(emit_shot_index + 1);

% get the azimuthal angle of the emitters in the rotation
phi_emitters = phi_fixed_emitters + phi_rotate_emitters;

% get the polar angle of emitters during the rotations
theta_emitters = theta_rad_fixed(emit_firing_transducer_index + 1);

% get the position of emitters
emitter.positions = radius_bowl * [cos(phi_emitters) .* sin(theta_emitters),...
    sin(phi_emitters) .* sin(theta_emitters), cos(theta_emitters)]';


% ensure the number of emitter positions is the same as the number of
% excitations
if size(emitter.positions, 2) ~= num_excitation
error(['The number of emitter positions is not'...
    'consistent with the number of excitations.'])
end

end