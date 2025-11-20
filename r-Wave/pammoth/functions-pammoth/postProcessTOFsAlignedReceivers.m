function [tof_matrix_processed] = postProcessTOFsAlignedReceivers(tof_matrix, connected_receiver_indices,...
       tol)
%POSTPROCESSTOFSALIGNEDRECEIVERS repalces the ouliers in the calculated
%difference TOF matrix with the median values for the adjacent receivers.

% DESCRIPTION:
%     postProcessTOFsAlignedReceivers replaces the outliers in the
%     calculated difference time-of-flight (TOF) matrix with the median values for the
%     adjacent receivers.
%
% USAGE:
%     
%
% INPUTS:
%     tof_matrix               - the num_emitter x num_receiver matrix of the
%                               difference TOFs [sec]
%     connected_receiver_indices - a cell array containg the index of adjacent 
%                                transducers to each receiver
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

% the number of emitters
num_emitter = size(tof_matrix, 2);
% the number of receivers
num_receiver = size(tof_matrix, 1);
if length(connected_receiver_indices) ~= num_receiver
    error('The connected receivers are not consistent with the size of the TOF vector');
end

% allocate a full-zero matrix for the processed difference tofs
tof_matrix_processed = zeros(size(tof_matrix));


for ind_emitter = 1:num_emitter
 

disp(['Replacing outliers for each emitter aligning the receivers,'...
    'The number of current emitter is:' num2str(ind_emitter)]);    
    
    
% get the difference tofs for the current emitter    
tof_emitter = tof_matrix(:, ind_emitter);
    

% allocate a vector for the processed difference tofs for the current emitter
tof_emitter_processed = zeros(num_receiver, 1);

for ind_receiver = 1:num_receiver
    
    % preprocess the tofs for the current emitter using the nozero difference tofs
    % for all the receivers for that emitter
    if tof_emitter(ind_receiver)
    
    % get the indices of close receivers to the current receiver
    connected_receivers_current = connected_receiver_indices{ind_receiver};
    
    % get the difference tofs of the close receivers to the current receiver
    tof_current = tof_emitter(connected_receivers_current);
    
    % remove the receivers with zero values
    tof_current = tof_current(abs(tof_current)>0);
    
    % if at least the tof for one of the neighboring transducers are nonzero,
    % calculate the median value
    if isempty(tof_current)
        
    tof_emitter_processed(ind_receiver) = tof_emitter(ind_receiver);     
        
    else
        
    % calculate the median value among the nonzeros tofs for the close
    % receivers
    median_current = median(tof_current);
    
    % replace the outliers with the median value
    if tof_emitter(ind_receiver)-median_current< tol(1) || tof_emitter(ind_receiver)-median_current> tol(2)
    tof_emitter_processed(ind_receiver) = median_current;    
    else
    tof_emitter_processed(ind_receiver) = tof_emitter(ind_receiver); 
    end
    
    end
     
    end
   

end


tof_matrix_processed(:, ind_emitter) = tof_emitter_processed;

end





end

