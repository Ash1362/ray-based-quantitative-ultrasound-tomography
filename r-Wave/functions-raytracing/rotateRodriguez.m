function [u_rotated] = rotateRodriguez(u, theta, u_ref)
%ROTATERESPECTEMITTER rotates a geomterical vector 
%
% DESCRIPTION:
% rotaterespectemitter applies a Rodriguez rotation to a geometrical vector
% connecting an emitter to a receiver such that the position of emitter is 
% the centre of rotation.
% 
% INPUTS:
%     u                     - the geomterical vector 
%     theta                 - the angle of rotation
%     u_ref                 - the axis of rotation

% OUTPUTS:
%      u_rotated            - the rotated geometrical vector

% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 30.12.2019
%       last update     - 30.12.2019
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian 


if any(u ~= 0)
    u = u./ norm(u);
end

if nargin < 3
        
% rotate counterclock-wise
rot = [cos(theta) -sin(theta);  sin(theta) cos(theta)];

u_rotated  = rot * u;
else
    
u_ref = u_ref./ norm(u_ref);  

A = [0, -u_ref(3), u_ref(2); u_ref(3), 0, -u_ref(1); -u_ref(2), u_ref(1), 0];
rot = cos(theta) * eye(3) + sin(theta) * A + (1 - cos(theta)) * (u_ref * u_ref.');
u_rotated  = rot * u;
end



end

