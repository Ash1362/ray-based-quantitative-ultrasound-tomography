function [tof_matrix_processed] = postProcessTOFsAlignedEmitters(tof_matrix,...
   connected_emitter_indices, emit_firing_transducer_index, rotation_indices,...
    tol)
%POSTPROCESSTOFSALIGNEDEMITTERS repalces the ouliers in the calculated
%difference TOF matrix with the median values for the adjacent emitters.

% DESCRIPTION:
%     postProcessTOFsAlignedEmitters replaces the outliers in the
%     calculated difference time-of-flight (TOF) matrix with the median values for the
%     adjacent emitters.
%
% USAGE:
%     
%
% INPUTS:
%     tof_matrix               - the num_emitter x num_receiver matrix of the
%                               difference TOFs [sec]
%     connected_emitter_indices - a cell array containg the index of adjacent 
%                                transducers to each emitter
%     emit_firing_transducer_index - a vector containging the index of transducers
%                                    for each excitation
%     rotation_indices         - a vector containging the index of transducers
%                                  for each excitation
%     tol                      - the lower/upper tolerance for detecting the
%                                outliers in the TOF map  
%
% OUTPUTS:
%     tof_matrix_processed  - the processed matrix of the difference TOFs
%                            [sec]
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 05.08.2020
%       last update     - 12.12.2020
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2020 Ashkan Javaherian
%


% the total number of emitters
num_emitter = size(tof_matrix, 2);

if num_emitter ~= length(rotation_indices)
    error(['The length of the mat vector for the rotation indices must be'...
        'the same as the number of columns in the TOF matrix']);
end

% the number of receivers
num_receiver = size(tof_matrix, 1);
if num_receiver ~= length(connected_emitter_indices)
    error(['The length of the cell array for the connected receivers must be'...
        'the same as the number of rows in the TOF matrix']);
end

% add one to the transducer indices for the emitters so that they are started
% from 1, instead of zero
emit_firing_transducer_index = emit_firing_transducer_index + 1;

% the number of rotational positions for the transducers
num_rot_positions = rotation_indices(end);

% allocate a full-zero matrix for the processed difference tofs
tof_matrix_processed = zeros(size(tof_matrix));


for ind_position = 1 : num_rot_positions
    
    
    disp(['Replacing outliers for each receiver aligning the emitters,'...
        'The number of current position is:' num2str(ind_position)]);
    
    % get the binaries for the current position
    binary_rot_position = ~(rotation_indices - ind_position);
    
    % get the tof matrix for the current rotational position of the
    % transducers
    tof_matrix_position = tof_matrix(:, binary_rot_position);
    
    % the number of emitters for the current position
    num_emitter_position = size(tof_matrix_position, 2);
    
    % get the rotational indices for the current position
    emit_firing_transducer_index_position = emit_firing_transducer_index(binary_rot_position);
    
    % allocate a full-zero matrix for the processed difference tofs for the
    % current rotational position
    tof_matrix_position_processed = zeros(num_receiver, num_emitter_position);
    
    
    
    for ind_receiver = 1:num_receiver
        
        
        % get the difference tofs for the current receiver
        tof_receiver = tof_matrix_position(ind_receiver, :);
        
        % allocate a vector for the processed difference tofs for the current
        % receiver
        tof_receiver_processed = zeros(1, num_emitter_position);
        
        for ind_emitter = 1:num_emitter_position
            
            % process the tofs for the current position and emitter using the nozeros difference tofs
            % for all the emitters for that receiver
            if tof_receiver(ind_emitter)
                
                % get the transducer index
                transducer_index = emit_firing_transducer_index_position(ind_emitter);
                
                % get the indices of close transducers to the current emitter (for this
                % rotational position)
                neighbor_transducers = connected_emitter_indices{transducer_index};
                
                % find the the neigboring transducers which are used for
                % thr current poisition
                [neighbor_binaries,~] = ismember(emit_firing_transducer_index_position,...
                    neighbor_transducers);
                
                % get the difference tofs of the close emitters to the current receiver
                tof_current = tof_receiver(neighbor_binaries);
                % remove the emitters with zero values
                tof_current = tof_current(abs(tof_current) > 0);
                
                % if at least the tof for one of the neighboring transducers are nonzero,
                % calculate the median value
                if isempty(tof_current)
                    
                    tof_receiver_processed(ind_emitter) = tof_receiver(ind_emitter);
                    
                else
                    % calculate the median value among the nonzeros tofs for the close
                    % emitters
                    median_current = median(tof_current);
                    
                    % replace the outliers with the median value
                    if tof_receiver(ind_emitter)-median_current < tol(1) || tof_receiver(ind_emitter)-median_current > tol(2)
                        tof_receiver_processed(ind_emitter) = median_current;
                    else
                        tof_receiver_processed(ind_emitter) = tof_receiver(ind_emitter);
                    end
                    
                end
                
            end
            
            
        end
        
        % fill the allocated matrix for the current position with the processed difference tofs for the current receiver
        tof_matrix_position_processed(ind_receiver, :) = tof_receiver_processed;
        
    end
    
    % fill the allocated matrix with the processed difference tofs for the current position
    tof_matrix_processed(:, binary_rot_position) = tof_matrix_position_processed;
    
end

end

