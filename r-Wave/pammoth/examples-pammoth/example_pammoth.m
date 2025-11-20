% Example:  Transmission UST image reconstruction from Pammoth data
%
% This example reconstructs the sound speed distribution from the UST portion of
% the Pammoth data. 
%
%
% author: Ashkan Javaherian
% date:            - 04.08.2020
% last update:     - 04.08.2020
%
%
% This script is part of the r-Wave Toolbox.
% Copyright (C) 2021 Ashkan Javaherian

clear all
close all
clc

% define the paths
startup_pammoth_ust;
%% ========================================================================
% OPTIONAL INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boolean controlling whether the time-of-flights are computed or not.
do_calculate_tofs = false; % (Default:true)

% Boolean controlling whether the computed time-of-flights are saved or not.
save_tofs = false; % (Default:true)

% Boolean controlling whether the image is reconstructed or not.
reconstruct_image = true;

% Boolean controlling whether the results (reconstructed images) are saved or not.
save_results = false; % (Default:true)

% get the approach for solving the linearised subproblem
linear_subproblem_method = 'sart';

%% ========================================================================
% GET THE DATA IDS
%==========================================================================
data_ids = [];

% the main directory, which must be fixed for all UST data
% (see startup_pammoth_ust.)
data_ids.main_path = main_path;

% get the object, the type of object for which the UST data have been
% measured. (check the m-file function 'setPathGetInfos')
data_ids.object = 'volunteers';

% get the ID of the data
data_ids.id = 'PAMMOTH-1-20-027/';

% get the date on which the data are measured
data_ids.date = '20201202_';
      
% get the measurement number as a string
data_ids.measurement_number = '07'; 

% get the cup size
data_ids.cup_size = 6;

% get the main path for saving the user's data. This should be first defined
% by the user, and will be fixed for all the reconstructions.
data_ids.main_path_user_data = 'data/pammoth/';


if save_results
    
    % get the directory for storing the results
    results_directory = [data_ids.main_path_user_data  'SoSreconstructions/', data_ids.id];
    
    % make the directory, if it does not exist.
    makeDirectory(results_directory)
    
else
    
    results_directory = [];
end

% get the optional inputs for loading ultrasound data and time-of-flight
% picking
reconst_args = {'do_calculate_tofs', do_calculate_tofs, 'reconstruct_image',...
    reconstruct_image, 'save_tofs', save_tofs, 'save_results', save_results,...
    'linear_subproblem_method', linear_subproblem_method};

% load the ultrasound data and compute (or load) the time-of-flights, and
% reconstruct the sound speed image
[img, recon_grid, data_ids, out, para] = reconstructSoundspeedImagePammoth(...
     data_ids, results_directory, reconst_args{:});

