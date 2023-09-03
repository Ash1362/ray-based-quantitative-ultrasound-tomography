function displayTransducerRotations(emitter, receiver, rotation_indices,...
    num_us_position, view_angle)
%DISPLAYTRANSDUCERSROTATIONS displays the transducers during the
%rotation of the bowl
%   
%
% DESCRIPTION:
%             displayTransducersRotations displays the transducers
%             during the rotation of the bowl
%      
%
% USAGE:
%     
%
% INPUTS:
%    emitter               - a struct with field
%    'positions'           - the position of emitters
%    receiver              - a struct with field
%    'positions'           - the position of receivers
%    rotation_indices      - the vector including the indices of angular
%                            positions for exitations
%    num_us_position       - the number of angular positions
%    view_angle           - the angle of view
%  
%       
% OPTIONAL INPUTS:
%    
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



    % visualise the receivers continously position by position
    figure;
    for ind_position = 1: num_us_position
        
        emitter_indices = find(~(rotation_indices-ind_position));
        
        % display emitters and receivers
        for i = ( (ind_position-1) * length(emitter_indices) + 1 ) : ind_position * length(emitter_indices) 
        scatter3(emitter.positions(1, i), ...
            emitter.positions(2, i),...
            emitter.positions(3, i),...
            50, 'g');
        hold on;
        scatter3(receiver.positions{ind_position}(1,:), ...
            receiver.positions{ind_position}(2,:),...
            receiver.positions{ind_position}(3,:),...
            10, 'k');
        pause(0.025)
        hold off;
        end
        if strcmp(view_angle, 'horizontal')
            view(0,90)
        end
        legend('emitters', 'receivers');
        pause(0.25)
        
    end
    
        
end