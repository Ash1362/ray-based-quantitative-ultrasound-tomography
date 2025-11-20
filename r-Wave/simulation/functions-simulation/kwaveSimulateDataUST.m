function [pressure_out, pressure_out_binary, data_paths, elapsed_time] = ...
    kwaveSimulateDataUST(kgrid, medium, emitter, receiver_positions, interp,...
    data_paths, varargin)
%KWAVESIMULATEDATA is a wrapper for calling kWave for simulation of data for
%UST imaging
%
% DESCRIPTION:
%        kwaveSimulateDataUSCT implements the k-Wave pseudo-spectral time domain
%        method using the k-Wave toolbox for simulating the pressure time series
%        used for UST imaging. k-Wave simulations are sequentially run by excitation
%        of emitters. If the code type is set 'Matlab', a number of simulations are
%        parallelized on different threads by a single call of function. The function
%        automatically shares the excitations (emitters) over the threads using
%        a parloop. However, if the simulations are run using high
%        performance codes, i.e. 'C++' or 'CUDA', a simultaneous simulation 
%        for different emitters using a single call of this function will not 
%        be manageable. 
%        Using 'CUDA', it is possible that the user share the emitters to 
%        the number of available threads, and run a script for each 
%        thread. The difference between the run scripts will be only
%        1- the chosen index for the GPU (This will be an optional input
%           for k-Wave mfile function kspaceFirstOrder3DG)
%        2- the chosen columns of emitter.positions
%        3- the chosen rows of emitter_interp_C (the chosen rows of 
%           emitter_interp_C must be the same as the chosen columns of
%           emitter.positions 

%        
%
% USAGE:
%       
%
% INPUTS:
%       kgrid              - a struct array defining the computational 
%                            grid in the k-Wave format
%       medium             - a struct array defining the medium's properties
%                            in the k-Wave format
%       emitter.positions  - 2/3 x N array of Cartesian points containing the centers of
%                            the emitter objects
%       emitter.pulse      - the excitational pulse
%       receiver_positions - 2/3 x N array of Cartesian points containing the centers of
%                            the receiver objects
%       interp.emitters    - a struct containing the interpolation optional parameters
%                            corresponding to the emitters. This contains
%                            the fields 'mask' (binary mask), 'matrix'
%                            (a sparse interpolation matrix)
%       interp.receivers    - a struct containing the interpolation optional parameters
%                            corresponding to the receivers. This contains
%                            the fields 'mask', 'matrix'
%                            (a sparse interpolation matrix), and maybe
%                            'mask_entire_volume' for recording the presure 
%                             time series on the grid points
%       data_paths         - a struct containing fields, 'main_directory', 
%                            'directory' and 'name_data'. 
%                            The 'HDF5' files will be saved on (or
%                            loaded from) the same directory.      
%
% OPTIONAL INPUTS:
%       'savetodisk'    - % control what information is saved. This can 
                          % 0, 1, or 2
                          % 0: the HDF5 file containing the simulation inputs is saved,
                          % and the pressure data is simulated and saved
                          % 1: only the pressure data is simulated and saved
                          % 2: only the HDF5 file containing the simulation inputs
                          % is saved without running the simulations
%       'data_cast'      - a matlab class for the output (Default ='single')
%       'code_type'      - the type of simulations supported by k-Wave 
%                         (Default= 'Matlab')
%       'num_worker_hdf5'   - the number of workers for creation and saving HDF5 files
%       'num_worker_data' - the numer of worker on CPUs or the total number of
%                         GPUs for executing k-Wave on different example
%                         scripts 
%       'gpu_index'     - the index of a single gpu for saving the data
%       'gpu_order'     - the order of the single GPU in use
%       'smooth_source'  - Boolean controlling whether the pressure sources
%                         should be smoothed in space 
%                         (Default = false)
%       'smooth_sos'     - boolean controlling whether the speed of sound
%                         map should be smoothed in the simulation
%                         (Default = true)
%       'PML_size'       - the number of PML layers along the different
%                         Cartesian coordinates
%
% OUTPUTS:
%       pressure_out    - A matrix representing the pressure data measured by
%                         the receivers. This matrix is of size
%                         num_receiver*kgrid.Nt*num_emitter with 
%                         num_receiver the number of receivers,
%                         Nt the length of k-Wave time array, and 
%                         num_emitter the number of the emitters
%       pressure_out_binary - A matrix representing the pressure data measured
%                         on a mask for the grid points.
%       data_paths      - a struct containing fields, 'main_directory', 
%                         'directory' and 'name_data', and also 'HDF5_path'
%                         if 'savetodisk' set 'true'.
%                         
%                         
%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 5.12.2019
%       last update     - 5.15.2019
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2020 Ashkan Javaherian


% =========================================================================
% CHECK INPUT STRUCTURES AND OPTIONAL INPUTS
% =========================================================================
% get the number of dimensions
dim = kgrid.dim;

% set k-Wave default optional_paramsmeter
optional_params = [];
optional_params.savetodisk = 1;
optional_params.data_cast  = 'single';
optional_params.code_type  = 'Matlab';
optional_params.num_worker_hdf5  = 4;
optional_params.num_worker_data  = 4;
optional_params.gpu_index = 0;
optional_params.gpu_order = 1;
optional_params.smooth_source = false;
optional_params.smooth_sos    = false;
optional_params.PML_size      = 10 + (dim == 2) * 10;

  
% read additional arguments or overwrite default ones (no sanity check!)
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to optional_paramsmeter struct (or overwrite fields)
        optional_params.(varargin{input_index}) = varargin{input_index + 1};
    end
end


% remove the fields that are used in parloop 
% from optional_params
code_type = optional_params.code_type;
gpu_index = optional_params.gpu_index;
data_cast = optional_params.data_cast;
rm_fields = {'code_type', 'gpu_index', 'data_cast'};
optional_params= rmfield(optional_params, rm_fields); 

% starting point for measuring simulation times
ts = tic;

% get the number of emitters 
num_emitter = size(emitter.positions, 2);

% get the number of receivers
num_receiver = size(receiver_positions, 2);

% get the size of the grid
grid_size = size(kgrid.x);

% get indices associated with the binary mask for emitters
emitter_indices = find(interp.emitter.mask);

% get the interpolation matrix for emitters
if strcmp(code_type, 'CUDA')
    
    % the chosen gpu order (1-based number of gpu) cannot be larger than the total number of
    % workers (gpus) 
    optional_params.num_worker_data = max(optional_params.num_worker_data,...
        optional_params.gpu_order);
         
    % if CUDA is used, the emitters are shared between gpus using multiple
    % main scripts. All the main scripts are the same, except for the 
    % input of gpu_index, which indicates the number of a single gpu
    % for running the associated k-Wave simulations. The input variable gpu_index 
    % allocates a gpu to a portion of emitters.
    
    
    % get the number of emitters per worker
    num_emitter_per_worker = num_emitter/optional_params.num_worker_data;
    
    if num_emitter_per_worker >= 1   &&   rem(num_emitter, optional_params.num_worker_data) ~= 0   
        error(['The number of emitters must be a natural factor of the number of the computational threads,'...
            'if the number of emitters are larger than the number of threads.'])
    else
    
    if num_emitter_per_worker < 1
        split_emitters_indices = (1:num_emitter)';
    else
        % choose the emitters using the index of the used gpu
        total_emitters_matrix = reshape(1:num_emitter,...
            [num_emitter_per_worker, optional_params.num_worker_data]);
        
        % the n-th column corresponds to the n-th gpu
        split_emitters_indices = total_emitters_matrix(:, optional_params.gpu_order);
    end
   
    % choose the position of emitters corresponding the gpu order
    emitter.positions = emitter.positions(:, split_emitters_indices);
    
    % the interpolation matrix for the chosen emitters
    emitter_matrix = interp.emitter.matrix(split_emitters_indices, :);
    
    % get the number of chosen emitters
    num_emitter = size(emitter.positions, 2);
    
    end
else
    % the interpolation matrix for the emitters
    emitter_matrix = interp.emitter.matrix;
end



if ~isempty(interp.receiver.mask_entire_volume)
mask_entire_volume = interp.receiver.mask_entire_volume;
indices = find(interp.receiver.mask);
indices_entire_volume = find(mask_entire_volume);
common_indices = ismember(indices_entire_volume, indices);
% set the sensor mask the entire inside the detection surface
sensor.mask = mask_entire_volume;
else
mask_entire_volume = [];    
% get the sensor mask
sensor.mask = interp.receiver.mask;
common_indices = [];
end
% get interpolation optional_parameters for receivers (fixed for all excitations)
receiver_matrix  = interp.receiver.matrix; 
% if the pressure data must be simulated and saved
if ismember(optional_params.savetodisk , [0, 1])
    pressure_out = zeros ([num_receiver, kgrid.Nt, num_emitter], data_cast );
    if ~isempty(interp.receiver.mask_entire_volume)
        pressure_out_binary = zeros ([nnz(mask_entire_volume), kgrid.Nt, num_emitter], ...
            'single');
    else
        pressure_out_binary = [];
    end
end

clear interp

% get the excitation pulse
pressure_in = emitter.pulse;

% if the HDF5 file must be saved
if ismember(optional_params.savetodisk , [0, 2])
    data_paths.HDF5_path = 'kwave_input_data/';
    HDF5_path = [data_paths.main_directory data_paths.directory data_paths.name_data data_paths.HDF5_path];
    makeDir(HDF5_path);
end

if optional_params.savetodisk == 2
    pressure_out = [];
    pressure_out_binary = [];
end 

% get the optional inputs for k-Wave simulation(s)
input_args = {'DataCast', data_cast , 'PMLSize', optional_params.PML_size,...
    'Smooth',[optional_params.smooth_source, optional_params.smooth_sos, false] ,...
    'PMLInside', false, 'PlotPML', false, 'PlotSim', false};

% =========================================================================
% SAVE INPUT DATA USING HDF5 FILES
% =========================================================================
if ismember(optional_params.savetodisk, [0, 2])
    
    parfor (emitter_index = 1:num_emitter, optional_params.num_worker_hdf5)
        % display the progess for saving the HDF5 files
        disp(emitter_index/num_emitter*100)
        % define a name containing the number of emitter for the HDF5 file
        input_args_HDF5 = {input_args{:}, 'savetodisk',...
            [HDF5_path 'kwave_input_data_'  int2str(emitter_index) '.h5']};
        % define a time-varying source (an input for k-Wave)
        source = [];
        source.p_mode = 'additive';
        source.p_mask = false(grid_size);
        source.p_mask(emitter_indices(find(emitter_matrix(emitter_index,:)))) = true;
        source.p = nonzeros(emitter_matrix(emitter_index,:)) * pressure_in;
        switch dim
            case 2
                kspaceFirstOrder2D(kgrid, medium, ...
                    source, sensor, input_args_HDF5{:});
            case 3
                kspaceFirstOrder3D(kgrid, medium, ...
                    source, sensor, input_args_HDF5{:});
        end
    end
end
clear input_args_HDF5
if strcmp(code_type, 'CUDA')
    input_args = {input_args{:}, 'DeviceNum', gpu_index};
end

% =========================================================================
% SIMULATE PRESSURE DATA
% =========================================================================
if ismember(optional_params.savetodisk, [0, 1])
    
    % if high performing codes are used, the user must share the excitations using multiple
    % main scripts, not this function
    if strcmp(code_type, 'CUDA') || strcmp(code_type, 'C++')
        optional_params.num_worker_data = 0;
    end
    
    parfor (emitter_index = 1:num_emitter, optional_params.num_worker_data)
        %   for emitter_index = 1: num_emitter
        
        % display the progess of k-Wave simulations
        disp(emitter_index/num_emitter*100)
        % define a time-varying source (an input for k-Wave)
        source = [];
        source.p_mode = 'additive';
        source.p_mask = false(grid_size);
        source.p_mask(emitter_indices(find(emitter_matrix(emitter_index,:)))) = true;
        source.p = nonzeros(emitter_matrix(emitter_index,:)) * pressure_in;
        
        switch dim
            case 2
                switch code_type
                    case 'Matlab'
                        pressure_on_mask = kspaceFirstOrder2D(kgrid, medium, ...
                            source, sensor, input_args{:});
                    case 'C++'
                        pressure_on_mask = kspaceFirstOrder2DC(kgrid, medium,...
                            source, sensor, input_args{:});
                    case 'CUDA'
                        Error('A 2D simulation using k-Wave has not been supported by CUDA yet.')
                    otherwise
                end
            case 3
                switch code_type
                    case 'Matlab'
                        pressure_on_mask = kspaceFirstOrder3D(kgrid, medium, ...
                            source, sensor, input_args{:} );
                    case 'C++'
                        pressure_on_mask = kspaceFirstOrder3DC(kgrid, medium, ...
                            source, sensor, input_args{:} );
                    case 'CUDA'
                        
                        pressure_on_mask = kspaceFirstOrder3DG(kgrid, medium, ...
                            source, sensor, input_args{:}  );
                    otherwise
                end
        end
        
        if ~isempty(mask_entire_volume)
            pressure_out_binary(:, : , emitter_index) = single(pressure_on_mask);
            % confine the pressure field on the grid points to those
            % corresponding to the receivers
            pressure_on_mask = pressure_on_mask(common_indices, :);
        end
        
        switch data_cast
            case 'single'
                pressure_out(:, :, emitter_index) = single(receiver_matrix * double(pressure_on_mask));
            case'double'
                pressure_out(:, :, emitter_index) = receiver_matrix * double(pressure_on_mask);
        end
       
        
        
    end
end

elapsed_time = toc(ts);

end