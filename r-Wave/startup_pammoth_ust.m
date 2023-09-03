% startup.m file for 'Ray-based toolbox'
%
% author: Ashkan Javaherian
% date: 30th June 2020
% last update: 30th June 2020

disp(['***Ray-based Toolbox:'...
    'Ray-based ultrasound imaging for biomedical tissues***']);

% set global variables
global root_path

% get the root path
root_path = pwd;

addpath(genpath(root_path))


% Get the main directory for the Pammoth data
% Get the machine name
% [~, machine_name] = system('hostname');

% set the machine name a local machine (You can use it unchanged, if you
% are using your local pc.
machine_name = 'ashkan';


if any(strcmp(machine_name(1:end-1), {'kinsler','nyborg'}))
    
    % For BUG servers, choose the following path
    % specify the main path for the data 
    main_path = '/mnt/Pammoth/Data/PAM3_measurements/';
    
    % display the machine name
    disp(['The machine name is:' machine_name(1:end-1)]);
else
    
    % the main path, where the information of the US data and time-of-flights are stored.
    % (default:data/pammoth/)
    % specify the main path for the data 
    main_path = 'data/pammoth/';
    
end