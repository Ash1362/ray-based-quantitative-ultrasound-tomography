function [kgrid, medium, emitter, receiver, simulation_prop, data_paths] = ...
    makeDataSettings(dim, dx, radius, z_pos_height, num_receiver, ratio_rtoe,...
    phantom, cfl, oa_breast_path, res_path, local_res_path, varargin)
%MAKEDATASETTINGS provides physical and computational settings for USCT
%data simulation

% DESCRIPTION:
% makeDataSettings is used for definition of the computational grid, emitters, receivers and intervening medium's
% acoustic properties for USCT data simulation
%
% USAGE:
% 
%
% INPUTS:
%       dim                    - the dimension of the grid
%       dx                     - the grid spacing  [m]
%       radius                 - the radius [m] of the bowl (hemi-spherical detection surface) 
%       z_pos_height           - the thickness [m] of the chest wall (This
%                                is zero for 2D case.)
%       num_receiver           - the number of receivers
%       ratio_rtoe             - the ratio of the number of receivers to emitters
%                               (this must be set 1,2, or 4)  
%       phantom                - the sound speed phantom (Default =
%                                'breast')
%       cfl                    - the CFL number
%       oa_breast_path         - the path for the MRI Breast phantom 
%                                (Mark Anastasio's phantom)
%       res_path               - the path for the results
%       local_res_path         - the local path for saving the required
%                                information

% OPTIONAL INPUTS:
%      DataSim                 - boolean controlling whether the data is 
%                                simulated and saved, or an already simulated and
%                                saved data set is loaded  
%      Scenario                - the scenario for imaging, 'standard' (for imaging),
%                                'single_emitter' (can be used for validating ray tracing)
%                                (default: 'standrad') 
%      Mode                    - the mode for imaging ('transmission'
%                                or 'reflection') (Default = 'transmission')
%      Excit                   - excitation pulse, which can be either 
%                               'Impulse', 'Pammoth_1' (Defualt =
%                               'Pammoth_1')
%      InterpType              - the method for mapping the pressure field from
%                                the grid points to the transducers and vice
%                                versa ('nearest', 'linear' or 'offGrid')
%                                (default = 'offGrid')
%      TransGeom               - the geomtery of the transducers ('point' for 2D or 3D,
%                               'line' for 2D, or 'disc' for 3D) (default = 'point')
%      singleER                - a vector of size 1 x 2 for choosing a pair of emitters
%                               and receivers, when Scenario = 'single_emitter'
%                               (default = [n_receiver/2, 1])                   
%      Absorption              - boolean controlling whether the absorption
%                                is included in simulations (default:'true')
%      Low_Filter              - boolean controlling whether a low-pass filter 
%                                is used for cutting-off the frequencies 
%                                of the excitation pulse that are not supported
%                                by the computational grid (default = 'true')
%      Low_Filter_Fac          - a factor multiplied by the minimum sound
%                                speed for adjusting the sound speed 
%                                for the low pass filter
%                                (the minimum sound speed of the medium may
%                                not sufficiently cut off the high
%                                frequencies, and these frequencies may lead to
%                                some oscilations in the measured signals
%      CodeVersion             - the k-Wave code version used for
%                                simulation (Default = 'Matlab' (2D case),
%                                'CUDA' (3D case))
%      CreatePhantom           - Boolean controlling whether the phantom is
%                                created or loaded. (default = true)
%      Plot                    - Boolean controlling whether the simulated settings are
%                                plotted (default = true)
% OUTPUTS:
%       kgrid                   - the k-Wave computational grid
%       medium                  - a k-wave struct containg the medium's
%                                 proprties:
%       medium.sound_speed      -  a matrix with the same size as the computational grid 
%                                  representing the sound speed [m/s] 
%       medium.sound_speed_ref  - a scalar representing the reference sound speed [m/s] 
%       medium.density          - a matrix with the same size as the computational grid 
%                                 representing the density [kg/m^3]
%       medium.alpha_power      - a scalar representing the power law absorption exponent
%       medium.alpha_coeff      - a matrix with the same size as the computational grid 
%                                 representing the absorption coefficient
%                                 [dB/(MHz^y cm)]
%       emitter                 - a struct defining the excitations                 
%       emitter.positions       - a matrix representing the position [m] of the 
%                                 centre of emitters in Cartesian coordinates
%                                 This has a size dim x num_emitter with num_emitter
%                                 the number of emitters                               
%       emitter.pulse           - a vector of size 1 x num_time_exc
%                                 with num_time_exc <= num_time_array. 
%                                 Here, num_time_exc and num_time_array are the length of
%                                 the excitation pulse and time array,
%                                 respectively. 
%       emitter.pulse_duration  - a scalar indicating the the approximate 
%                                 time duration of the excitation pulses, 
%                                 whichmay be used for processing of 
%                                 the simulated data
%       receiver.positions      - a matrix representing the position [m] of the
%                                 centre of receivers in Cartesian coordinates
%                                 This has a size dim x num_receiver with num_receiver
%                                 the number of receivers 
%       simulation_prop         - a struct containing some information
%                                 about the simulation:
%       simulation_prop.PML     - a vector of size 1 x dim representing
%                                 the number of PML layers added to each side of
%                                 the grid along the Cartesian coordinates
%       simulation_prop.f_max   - a scalar representing the maximum frequency 
%                                 of the excitations pulse supported by the computional grid
%       simulation_prop.z_offset - the off-set position of the transducers along 
%                                 z axis in 'kgrid' compared to the real case, because the Cartesian 
%                                 coordinates in 'kgrid' are symmetric with respect to the origin 
%      data_paths               - a struct containing the path for saving (or loading) the data
%      data_paths.main_directory- the main path for data storage
%      data_paths.directory     - the path for data storage 
%      data_paths.name_data     - a name for the folder containing the saved data 

%
% ABOUT:
%       author          - Ashkan Javaherian
%       date            - 29.12.2019
%       last update     - 29.12.2019
%
% This script is part of the r-Wave Tool-box (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian and Ben Cox

para.DataSim        = false;
para.Scenario       = 'standard';
para.Mode           = 'transmission';
para.Excit          = 'Pammoth_1'; 
para.InterpType     = 'offGrid';
para.TransGeom      = 'point';
para.Absorption     = false;
para.singleER       = [num_receiver/2, 1];
para.Low_Filter     = true;
para.Low_Filter_Fac = 1;
para.CodeVersion    = 'Matlab';
para.CreatePhantom = true;
para.Plot           = true;

% add to parameter struct 
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% get the number of emitters
num_emitter = num_receiver/ratio_rtoe;

% get the size [m] of the domain
diameter = 2 * radius; % [m]

% get the radius of the detection ring (surface)
detec_radius = 0.95 * radius;

%% ========================================================================
%  DEFINE A PATH FOR STORING THE DATA
%  ========================================================================
data_paths = [];
data_paths.main_directory = res_path;

data_paths.directory = ['study_' para.Mode '/' num2str(dim) 'D'  '/'];
    

%if strcmp(data_paths.main_directory(1:10), 'simulation' )
%data_paths.name_data = 'data1/';
%else
data_paths.name_data = ['Pulse' para.Excit  '_dx' num2str(dx,'%1.1e') '_Nr'...
num2str(num_receiver) '_Ne' num2str(num_emitter)  '_Interp' para.InterpType...
'_Transgeom' para.TransGeom '_Absorption' num2str(para.Absorption) '_Code'   para.CodeVersion  '/'];
%end

% Make the directory for the data
makeDirectory([data_paths.main_directory, data_paths.directory, data_paths.name_data]);

% display the directory for the data
if para.DataSim
    disp(['The data is simulated and will be stored at the directory:'...
        data_paths.main_directory, data_paths.directory, data_paths.name_data])
else
    disp(['The already simulated data is loaded from the directory:'...
        data_paths.main_directory, data_paths.directory, data_paths.name_data])
end

%% ========================================================================
%  CONSTRUCT THE COMPUTATIONAL GRID
%  ========================================================================
switch phantom
    case {'circles', 'circles_test'}
        Nx = intManipulation(ceil(diameter/dx),'evenUp');
        
        switch dim
            case 2
                kgrid  = kWaveGrid(Nx, dx, Nx, dx);
            case 3
                kgrid  = kWaveGrid(Nx, dx, Nx, dx, Nx, dx);
        end
        
    case 'breast'
        phantom_grid = makePAMMOTHGeometry(radius, dx, z_pos_height, 0, 0);
        
        switch dim
            case 2
                Nx    = phantom_grid.Nx;
                kgrid = kWaveGrid(Nx, dx, Nx, dx);
            case 3
                kgrid = phantom_grid.kgrid;
        end
        
end

%% ========================================================================
%  CHOOSE THE NUMBER OF THE PML LAYERS
%  ========================================================================
switch dim
    case 2
         simulation_prop.PML = kWave_getOptimalPMLSize([kgrid.Nx, kgrid.Ny], [20,30]);
    case 3
         simulation_prop.PML = kWave_getOptimalPMLSize([kgrid.Nx, kgrid.Ny, kgrid.Nz], [20,30]);
end              
        
% display the number of Grid points along the Carstesian coordinates with and without PMLs      
disp(['Nx = ' int2str(kgrid.Nx) ', PML = ' int2str(simulation_prop.PML(1))])
disp(['factors(Nx + 2 * PML): ' int2str( factor( kgrid.Nx + 2 * simulation_prop.PML(1)) ) ])
if dim == 3
disp(  ['factors(kgrid.Nz + 2 * PML: ' int2str( factor( kgrid.Nz + 2 * simulation_prop.PML(3) )  )]   )
end


%% ========================================================================
%  DEFINE THE POSITION OF TRANSDUCERS
%  ========================================================================

switch dim
    
    case 2
        
        % The emitters and receivers are chosen so that they are not
        % overlapped 
        % get the Cartesian position of the receivers
        receiver.positions = makeCartCircle(detec_radius, 2*num_receiver);
        receiver.positions = receiver.positions(:, 2:2:end);
        % get the Cartesian position of the emitters  
        emitter.positions = makeCartCircle(detec_radius, round(2*num_receiver/ratio_rtoe));
        emitter.positions = emitter.positions(:, 1:2:end);
        
        
    case 3
        
        % The emitters and receivers are overlapped for an equal number of
        % emitter and receivers
        % get the Cartesian position of the receivers 
        receiver.positions = makeCartHalfSphere(detec_radius, num_receiver, [0,0,0], false);
        % get the Cartesian position of the emitters
        emitter.positions = makeCartHalfSphere(detec_radius, round( num_receiver/ratio_rtoe ), [0,0,0], false);
        
end



switch dim
    case 2
        z_offset = [];
    case 3
        
        % correct the off-set position of the transducers along z axis 
        % For example, z=(-13,0)cm in 'setting' will give
        % z=(-6.5,6.5)cm in 'kgrid', because the Cartesian 
        % coordinates in 'kgrid' are symmetric with respect to the origin 
        if strcmp (phantom, 'breast')
            z_offset = kgrid.z_vec(1)-phantom_grid.z_vec(1);
            % transform the z positions from those in setting to those in kgrid 
            receiver.positions(3,:) = receiver.positions(3,:) + z_offset;
            emitter.positions(3,:) = emitter.positions(3,:) + z_offset;
        else
            z_offset = 0;
        end
end



distance_min = inf;
count_small_distances = 0;
for emitter_index = 1:num_emitter
    for receiver_index = 1:num_receiver
        cur_distance = norm(receiver.positions(:, receiver_index)-emitter.positions(:, emitter_index));
        % update the minimum ditance between emitters and receivers
        distance_min = min(distance_min, cur_distance);
        % count the number of emitter-receiver pairs that are 
        % mapped to the same grid point, if 'nearest' interpolation is
        % used, or they are closer than the half grid spacing.
        % mapping the transducers to the same grid points does not cause any
        % problem for data generation, and is counted only for information
        if cur_distance < 1/2*sqrt(2)*kgrid.dx
            count_small_distances = count_small_distances+1;
        end
    end
end

disp(['The minimum distance between emitters and receivers are:'....
    num2str(distance_min,'% 5.2f' )] );
disp(['The number of distances between emitters and receivers of smaller than grid spacing are:'....
    num2str(count_small_distances)] );

% Choose a single emitter and all receivers, or an emitter-receiver pair,
% if the scenario is set 'single_emitter'. 
if strcmp(para.Scenario, 'single_emitter')
    if isfinite(para.singleER(1))
        if para.singleER(1) <= num_receiver
            receiver.positions = receiver.positions(:, para.singleER(1));
        else
            error(['If the index of receiver is finite, it must not be larger than the number'...
                'of receivers.'])
        end
    end
    if isfinite(para.singleER(2))
        if para.singleER(2) <= num_emitter
            emitter.positions = emitter.positions(:, para.singleER(2));
        else
            error(['The index of emitter must be finite and not larger than the number'...
                'of emitters.'])
        end
    else
        error(['The index of emitter must be finite and not larger than the number'...
            'of emitters.'])
    end
end
     
    
%% ========================================================================
%  DEFINE THE MEDIUM'S SOUND SPEED
%  ========================================================================
% only the breast phantom is tested.
switch phantom
    
    case 'breast'
         % get the local directory for saving the phantom
            local_phantom_path = [local_res_path, 'phantom/' phantom '/'];
            
            % make the directory
            makeDirectory(local_phantom_path)
            
            % get the reference sound speed
            medium.sound_speed_ref = 1500; % [m/s]
            
        if para.CreatePhantom
            
            % get the phantom
            phan = getBreastPhantom('Neg_47_Left', oa_breast_path, phantom_grid, false);
            
            % get the modified acoustic properties
            medium_acoustic = acousticPropertiesModified(phantom_grid, phan, para.Absorption, false);
            
            % save the struct for the created phantom in the path
            save([local_phantom_path,  num2str(1e5 * dx), '.mat'], 'medium_acoustic', 'phantom_grid', '-v7.3')
            
        else
            
            % load the struct for the already created phantom from the path
            load([local_phantom_path,  num2str(1e5 * dx), '.mat'], 'medium_acoustic')
             
        end
            switch dim
                case 2
                    
                    % choose a horizontal slice of the 3D phantom for 2D imaging
                    medium.sound_speed = medium_acoustic.sound_speed(: ,: ,round(3/4 * size(medium_acoustic.sound_speed, 3)));
                    
                    if para.Absorption
                    medium.alpha_coeff = medium_acoustic.alpha_coeff(: ,: ,round(3/4 * size(medium_acoustic.sound_speed, 3)));
                    medium.alpha_power = 1.4;
                    
                 
                    end
                    
                case 3
                    
                    medium.sound_speed = medium_acoustic.sound_speed;
                    
                    if para.Absorption
                    medium.alpha_coeff = medium_acoustic.alpha_coeff;
                    medium.alpha_power = 1.4;
                  
                    
                    end
                    
                otherwise
                    error('The dimension must be set 2 or 3.');
            end
            if para.Absorption
            disp(['The exponent power is:' num2str(medium.alpha_power)]); 
            end
            clear phan medium_acoustic
            

        
    case 'simple'  
        if para.DataSim 
            medium.sound_speed_ref = 1500; % [m/s]
            switch dim
                case 2
                    notImpErr
                    medium.sound_speed = medium.sound_speed_ref * ones(kgrid.Nx,kgrid.Ny);
                    medium.sound_speed(5*(kgrid.x - 0.025).^2 + 2*(kgrid.y - 0.025).^2 < 0.20*detec_radius^2) = medium.sound_speed_ref + 48;   % [m/s]
                    medium.sound_speed(3*(kgrid.x + 0.025).^2 + 3*(kgrid.y + 0.025).^2 < 0.10*detec_radius^2) = medium.sound_speed_ref - 40;   % [m/s]
                    medium.sound_speed(3*(kgrid.x + 0.030).^2 + 3*(kgrid.y - 0.030).^2 < 0.15*detec_radius^2) = medium.sound_speed_ref + 80;   % [m/s]
                    
                case 3
                    notImpErr
            end
        else
            notImpErr
        end
        
    case 'circles'  
        notImpErr
        medium.sound_speed_ref = 1500; % [m/s]
        switch dim
            case 2
            medium.sound_speed_ref  = 1500; % [m/s]
            medium.sound_speed = medium.sound_speed_ref * ones(kgrid.Nx, kgrid.Ny);   
            xvec = linspace(-0.6 * radius, +0.6 * radius, 5);
            yvec = linspace(-0.6 * radius, +0.6 * radius, 5);
            for ii = 1:5
                for jj = 1:5
                    if rem(ii+jj, 2) == 0
                        medium.sound_speed((kgrid.x-xvec(ii)).^2 + (kgrid.y-yvec(jj)).^2 < (0.05*detec_radius)^2) = medium.sound_speed_ref + 50;   % [m/s]
                    else
                        medium.sound_speed((kgrid.x-xvec(ii)).^2 + (kgrid.y-yvec(jj)).^2 < (0.05*detec_radius)^2) = medium.sound_speed_ref - 50;   % [m/s]
                    end
                end
            end 
            case 3
            notImpErr
            otherwise
             notImpErr
        end    
    otherwise
        notImpErr     
end
% add a homogeneous density distribution
medium.density = 1000; % *ones(size(kgrid.x));
%% ========================================================================
%  DEFINE THE MEDIUM'S ABSORPTION
%  ========================================================================
if para.Absorption
    switch phantom
        case 'simple'
    switch dim
        case 2
    notImpErr
    alpha_coeff_ref = 0.002;
    medium.alpha_coeff = alpha_coeff_ref * ones(kgrid.Nx, kgrid.Ny);
    medium.alpha_coeff(5*(kgrid.x-0.025).^2 + 2*(kgrid.y+0.025).^2 < 0.20*detec_radius^2) = 0.5;     % dB MHz^{-1} {cm}^{-1}
    medium.alpha_coeff(3*(kgrid.x+0.040).^2 + 3*(kgrid.y+0.040).^2 < 0.10*detec_radius^2) = 0.75;    % dB MHz^{-1} {cm}^{-1}
    medium.alpha_coeff(3*(kgrid.x+0.010).^2 + 3*(kgrid.y-0.010).^2 < 0.15*detec_radius^2) = 1;       % dB MHz^{-1} {cm}^{-1}
    medium.alpha_power = 1.5;
        case 3
            notImpErr
    end
        case{'breast'}
         % included in acousticPropertiesModified   
        otherwise
            notImpErr
    end
end
            
           
%% ========================================================================
%  DEFINE THE TIME SERIES FOR K_WAVE SIMULATION
%  ========================================================================
switch para.Mode
    case 'transmission'
         t_end =  1.2 * diameter / min(medium.sound_speed(:));   % [s]
    case 'reflection'
         t_end =  2.3 * diameter / min(medium.sound_speed(:));   % [s]
end

% get the time array associated with the chosen CFL
kgrid.makeTime(medium.sound_speed, cfl, t_end);

%if strcmp(para.Excit, 'Pammoth_1')
%    kgrid.dt = 5000e-11;
%    kgrid.t_array = 0:kgrid.dt:t_end;
%end
    
%% ========================================================================
%  DEFINE (OR LOAD) THE EXCITATION PULSE
%  ========================================================================
% pulse_prop is a struct for properties of the excitation pulse
% only 'Pammpth_1' is tested.
switch para.Excit
    
    case 'Impulse'
        notImpErr
        % For an 'Impulse' excitation pulse, the defined parameters will be not
        % used for generation of the pulse, but will be later used as properties
        % of the signal for data pre-processing.
        pulse_prop.T  = 1.6e-6;    % the time duration of the excitation pulse [s]
    case {'Impulse-Dirichlet','Impulse-Gaussian'}
        notImpErr
        % define a coefficient for adjusting the sharpness of the Dirac delta
        % function, reciprocal to the temporal width of the excitation pulse
    
        pulse_max_freq = 1.6e6;           
        pulse_prop.T  = 1/pulse_max_freq;
        % the time peak of the approximated Dirac delta function [s]
        pulse_prop.t0 = 5e-6;      
    case 'Three-cycle-sinusoid'
        notImpErr
        pulse_prop.fc = 0.8e6;     % carrier frequency [Hz]
        pulse_prop.b  = 0.75e-6;   % temporal standard deviation of a Gaussian window [s]
        pulse_prop.T  = 3.2e-6;    % temporal mean of a Gaussian window [s]
        
    case 'Chirp'
        notImpErr
        %  'These parameters are used in: N.V. Ruiter, M. Zapf, R. Dapp, T. Hopp, W.A. Kaiser and H. Gemmeke, First Results of a Clinical Study with
        %  3D Ultrasound Computer Tomography, Joint UFFC, EFTF and PFM Symposium, 2013.)
        pulse_prop.fc = 2.5e6;     % central frequency [Hz]
        pulse_prop.b  = 1.0e6;     % frequecy Band-width [Hz]
        pulse_prop.T  = 12.8e-6;   % pulse duration [s]
        
    case 'Pammoth_1'
        pulse_prop.T  = 3.0e-6;   % the time duration of the excitation pulse [s]
        % For a realistic excitation pulse, e.g. 'Pammoth_1',
        % the corresponding sampling rate must be given
        pulse_prop.dt = 5000e-11;
        
end

% store the maximum supported frequency
simulation_prop.f_max = kgrid.k_max*min(medium.sound_speed(:))/(2*pi);  %   [Hz]
disp(['Maximal frequency supported by the computational grid:' num2str(simulation_prop.f_max, '%1.1e') '[Hz]'] );
% give as an output the off-set z (the z off-set between the z coordinates
% between the k-Wave and the rela z coordinate
simulation_prop.z_offset = z_offset;

% give as an output the radius of the dection surface (ring)
simulation_prop.detec_radius = detec_radius;

% simulate the excitation pulse
exc_args = {'Excit', para.Excit, 'Low_Filter', para.Low_Filter, 'Low_Filter_Fac',...
    para.Low_Filter_Fac, 'Plot', para.Plot, 'Plot_frequ_maxaxis', simulation_prop.f_max };
emitter.pulse = excitUSCTPulse(kgrid, medium.sound_speed, pulse_prop, exc_args{:});
emitter.pulse_duration = pulse_prop.T;

% the peak of the exciation pulse in time, if 'Impulse-Dirichlet' or 'Impulse-Gaussian'
% is used (approximated delta functions)
if isfield (pulse_prop, 't0')
    emitter.pulse_t0 = pulse_prop.t0;
end


%% ========================================================================
%  PLOT THE RESULTS
%  ========================================================================
if (para.Plot)
    switch dim
        case 2
            figure();
            imagesc(medium.sound_speed);axis image;colorbar;
            
            if para.Absorption
            figure();
            imagesc(medium.alpha_coeff);axis image;colorbar;
            end
            
            figure();
            scatter(emitter.positions(1,:), emitter.positions(2,:),'b');
            hold on;
            scatter(receiver.positions(1,:), receiver.positions(2,:),'r');
            axis image; legend('emitters', 'receivers'); hold off;
        case 3
            scrollView(medium.sound_speed, 3, []);colorbar;
            figure();
            scatter3(emitter.positions(1,:), emitter.positions(2,:),...
                emitter.positions(3,:), 'b');hold on;
            scatter3(receiver.positions(1,:), receiver.positions(2,:),...
                receiver.positions(3,:), 'r');
            axis image; legend('emitters', 'receivers'); hold off;
    end
end




end