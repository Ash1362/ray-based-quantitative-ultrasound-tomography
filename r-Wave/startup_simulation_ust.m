% startup.m file for ray-based toolbox
%
% author: Ashkan Javaherian
% date: 30th June 2020
% last update: 20 November 2025

% display the Tool-box title and the last version
disp(['***r-Wave: An open-source MATLAB package for quantitative ultrasound tomography' ...
    ' via ray-Born inversion with in vitro and in vivo validation***']);
disp('Version 1.1, Released 3rd September 2025');
disp('Version 1.2, Released 20 November 2025');

% set global variables
global root_path

% get the root path
root_path = pwd;

% add the subdirectories
addpath(genpath(root_path))



% get the main directory for saving and later loading the data
% get the machine name
[~, machine_name] = system('hostname');

% for testing if the codes work with any other machines.
% the user doesn't need to change it.
machine_name = 'ashkan';


% 'kinsler' and 'nyborg' are the machine names used by the code's author.
if any(strcmp(machine_name(1:end-1), {'kinsler','nyborg'}))
    
    % display the machine name
    disp(['The machine name is:' machine_name(1:end-1)]);
    
    % add the directory for the results
    res_path = '/mnt/Pammoth/user/Ashkan/data/';
    
    % get the k-Wave directory
    % kwave_path  = '/mnt/Pammoth/Software/k-wave-toolbox-version-1.3/k-Wave';
    % get the k-wave path
    kwave_path = 'k-wave-toolbox-version-1.3/k-Wave';
    % get the path for Mark Anastasio's breast phantom 
    oa_breast_path = '/mnt/Pammoth/Software/OA-BREAST/';
    

    addpath(oa_breast_path)
    
else
    
    % A directory for saving and later loading the data must be specified by the user.  
    % For example, the user could create a folder data in the current folder
    res_path = 'data/simulation/';
    
    % make the directory for saving and later loading the data 
    makeDirectory(res_path);
    
    
    if isfolder('k-wave-toolbox-version-1.3')
        % get the k-wave path
        kwave_path = 'k-wave-toolbox-version-1.3/k-Wave';
        
    else
        error(['The user must specify the path for the k-wave toolbox (kwave_path) in the next line.'...
            'The user may see how the author of the codes specifies kwave_path on his own machine '...
            'by checking the equivalent line command written above.']);
    end
    
    % by setting oa_breast_path nan, the breast phanstom is loaded from the
    % local directory
    % oa_breast_path = nan;
    oa_breast_path = 'data/simulation/phantom/OA-BREAST/';
    
end

% add the paths
addpath(kwave_path)
% get the path for saving the data on the local machine. res_path will be used
% for saving and later loading the portion of data, which require larger memory,
% but local_res_path is used for saving and loading the portion of data, which require
% smaller memory, e.g., phantom,... . By default, the same directory name is used
% for res_path and local_res_path, if the user is not using UCL's servers,
% 'nyborg' or 'kinsler', but the user may like to adjust the paths.
local_res_path = 'data/simulation/';

% make the directory for saving the data on the local machine, or
% verify if the local directory already exists
makeDirectory(local_res_path);

% messages informing the users about the external software resources used for developing the ray-based toolbox

% the digital breast phantom developed by Mark Anastasio's group
disp(['This toolbox uses the digital oa-breast-database phantom:'...
    'Please download it via the link:'...
    'https://anastasio.bioengineering.illinois.edu/downloadable-content/oa-breast-database/' ...
    'and add it to the path:...simulation/data/phantom/OA-BREAST/Neg_47_Left...'...
    'Lou Y, Zhou W, Matthews T P, Appleton C M and Anastasio M A 2017,'...
    'Generation of anatomically realistic numerical phantoms for photoacoustic and ultrasonic breast'...
    'imaging, J. Biomed. Opt. 22041015.']);

% the k-Wave toolbox
disp(['This toolbox uses k-wave toolbox for wave simulations. Please read: http://www.k-wave.org/license.php'...
    'B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic'...
    'wave-fields," J. Biomed. Opt., vol. 15, no. 2, p. 021314, 2010.'])

% the off-grid toolbox
% disp(['This toolbox uses off-grid toolbox. Please cite:'...
% 'E. S. Wise, B. T. Cox, J. Jaros, B. E. Treeby, Representing arbitrary acoustic source and sensor distributions in'...
% 'Fourier collocation methods, J. Acoust. Soc. of Am., vol. 146, no. 1, pp. 278-288, 2019.'])



