function makeDirectory(info_path)
%MAKEDIRECTORY makes a directory for storing the information
%
% DESCRIPTION: 
%   makeDirectory check if a directory does not exist, and if so creates a direcctory
%   for storing the information 
%
%  USAGE:
%       
%
%  INPUTS:
%       info_path          - the path for storing the information

%
%  OUTPUTS:

%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 05.12.2019
%       last update     - 05.12.2019
%
% See also
%
% This function is part of the r-Wave Toolbox.
% Copyright (c) 2022 Ashkan Javaherian
% 

if ~exist(info_path, 'dir')

    disp(['Creating new directory: ' info_path])
    [status, msg, msg_id] = mkdir(info_path);
    
    if ~status
        warning(msg, msg_id)
    end
    
else
    
    disp('The requested directory for making already exists.')
end

end