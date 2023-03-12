function [angles] = calcAngleEmitterReceiver(emitter_positions,...
    receiver_positions, rotation_indices)
%CALCANGLESEEMITTERRECEIVER calculates the angle between vectors normal to
% emitters and the vectors from the emitters to receivers
%
%
% DESCRIPTION:
%       calcAngleEmittersReceivers calculates the angle between the geometerical vectors
%       normal to the emitters and the angles from the geometrical vectors from emitters
%       to the receivers
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
%
%
% OPTIONAL INPUTS:
%
%
%
% OUTPUTS:
%        angles                 - the angle between vectors normal to
%                                 emitters and the vectors from the emitters
%                                 to the receivers
%
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 


num_emitter = size(emitter_positions, 2);


% the number of receivers
if  iscell(receiver_positions)   
    num_receiver = size(receiver_positions{1}, 2);
else
    num_receiver = size(receiver_positions, 2); 
end

    angles = zeros(num_receiver, num_emitter);
    
    % if the position of the receivers is changed with the position of
    % emitter
    for ind_emitter = 1 : num_emitter
        
        % calculate the Cartesian vectors from the emitter to the receivers
        if iscell(receiver_positions)
            cartesian_direction_emitter_receiver = ...
                receiver_positions{rotation_indices(ind_emitter)} - emitter_positions(:, ind_emitter);
        else
            cartesian_direction_emitter_receiver = ...
                receiver_positions - emitter_positions(:, ind_emitter);
        end
        
        for ind_receiver = 1: num_receiver
            
            angles(ind_receiver, ind_emitter) = atan2(norm(cross(-emitter_positions(:, ind_emitter),...
                cartesian_direction_emitter_receiver(:, ind_receiver))),...
                dot(-emitter_positions(:, ind_emitter),...
                cartesian_direction_emitter_receiver(:, ind_receiver)));

        end
    end
    
end
