function [distances] = calculateDistanceEmitterReceiver(emitter_positions,...
    receiver_positions, rotation_indices)
%CALCULATEDISTANCEEMITTERRECEIVER calculates the distance between each
%pair of emitters and receivers
%
%
% DESCRIPTION:
%       calculateDistanceEmitterReceiver calculates the distance between
%       each pair of emitters and receivers using the position of emitters
%       and receivers
%
% USAGE:
%     
%
% INPUTS:
%       emitter_positions       - 2/3 x N array of cartesian points containing 
%                                 the centers of the emitter objects                             
%       receiver_positions      - 2/3 x N array of cartesian points containing 
%                                 the centers of the receiver objects
%                                 For the cases at which the position of
%                                 the receivers is changed with position of
%                                 the emitter, this is a cell array
%                                 consisting of  2/3 x N array of cartesian
%                                 points for each position of emitters
%       rotation_indices        - the indices of angular positions for each
%                                 emitter
%                                
%
% OPTIONAL INPUTS:
%       
%                                
%                                 
% OUTPUTS:
%        distances             - the distances between all pairs of
%                                emitters and receivers
%                                                                              
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2020 Ashkan Javaherian 


num_emitter = size(emitter_positions, 2);
% the number of receivers
if  iscell(receiver_positions)
    
    num_receiver = size(receiver_positions{1}, 2);
    distances = zeros(num_receiver, num_emitter);
    
    % if the position of the receivers is changed with the position of
    % emitter
    for ind_emitter = 1 : num_emitter

            distances(:, ind_emitter) ...
                = vecnorm(receiver_positions{rotation_indices(ind_emitter)}...
                - emitter_positions(:, ind_emitter))';

    end
    
else
    num_receiver = size(receiver_positions, 2);
    distances = zeros(num_receiver, num_emitter);
    
    % if the position of the receivers is fixed with respect to the
    % position of emitter
    for ind_emitter = 1 : num_emitter
        for ind_receiver = 1 : num_receiver
            distances(ind_receiver, ind_emitter) ...
                = norm(receiver_positions(:, ind_receiver)...
                - emitter_positions(:, ind_emitter));
        end
    end
    
end
     
end