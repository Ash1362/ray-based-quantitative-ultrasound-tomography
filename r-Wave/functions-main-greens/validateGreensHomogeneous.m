function [relative_discrepancy_signal_water_emitter, relative_discrepancy_signal_water_receiver] =...
    validateGreensHomogeneous(data_water, time_array, emitter, receiver,...
    sound_speed_water, f_max, plot_directory, varargin)
%VALIDATEGREENSHOMOGENEOUS compares the acoustic waves approximated by 
% the Green's function and those simulated by k-Wave in homogeneous medium
%   
%
% DESCRIPTION:
%             validateGreensHomogeneous compares the acoustic waves
%             approximated by the Greens' function and those simulated by
%             k-Wave in homogeneous medium.
%      
%
% USAGE:
%     
%
% INPUTS:
%       data_water        - the ultrasound time series measured from water.
%                           This is a matrix of size num_receiver x num_time
%                           x num_emitter with num_emitter the number of emitters,
%                           num_time the number of time instants for measurement
%                           and num_receiver the number of receivers
%       time_array        - a time array of size 1 x num_time 
%       emitter           - a struct which defines the properties of the
%                           excitation as follows: This includes the
%                           fields 'positions', 'pulse', 'pulse_duration',
%                           and 'shot_time'
%       emitter.positions - 2/3 x num_emitter array of Cartesian position 
%                           of the center of the emitter objects
%       emitter.pulse     - a vector of size 1 x num_te with num_te <= num_time.
%       emitter.shot_time - the first arrival of the excitation pulse, which
%                           can be negative [sec]
%       receiver          - a struct which includes the fields
%       receiver.positions - 2/3 x N array of Cartesian position of the
%                           centre of the receiver objects. For rotating
%                           ultrasound systems, the position of the receivers
%                           is changed for different excitations. For the rotating case,
%                           this is a cell array with length num_angle (the number of rotation
%                           angles) with each cell a 2/3 x num_transducer array of
%                           Cartesian position of transducers for each
%                           angle. (num_receiver = num_angle x num_transducer)
%      sound_speed_water   - the sound speed [m/s] in only water
%      f_max               - the maximum frequency supported by the grid
%                           for the k-Wave simulation
%      plot_directory      - the drectory for saving the plots 
%  
%      
 
% OPTIONAL INPUTS:
%      'save_plots'        - Boolean controlling whether the plots are saved
%                            or not.
%      'deconvolve_source' - Boolean controlling whether the measured signals
%                            are deconvolved from the pressure source or
%                            not.
%      'deconvolution_parameter' - get the regularisation parameter for
%                            deconvolution
%      'choose_receiver'   - the approach for choosing the receiver for
%                           plotting. The receiver can be chosen by index 
%                           ('index'), or alternatively by distancves to
%                           the chosen emitter ('distance')
%      'num_discretised_frequencies' - the number of discretised frequencies
%                          within the interval
%      'emitter_index'    - the index of a chosen emitter for analysis
%      'receiver_index'   - the index of a chosen receiver for analysis
%       
%
% OUTPUTS:
%      relative_discrepancy_signal_water_emitter -the percentage relative
%                                                 discrepancy between the
%                                                 signals simulated by
%                                                 k-Wave and the signals
%                                                 approximated by the
%                                                 Green's function in
%                                                 homogeneous water on all
%                                                 the receivers after being
%                                                 produced by a single
%                                                 chosen emitter.
%      relative_discrepancy_signal_water_receiver - the percentage relative
%                                                 discrepancy between the
%                                                 signals simulated by
%                                                 k-Waveand the signals
%                                                 approximated by the
%                                                 Green's function in
%                                                 homogeneous water on a
%                                                 single chosen receiver
%                                                 after being produced by a
%                                                 single chosen emitter.       
%    
%
% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 10.08.2022
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian 


para.save_plots = false;
para.deconvolve_source = false;
para.deconvolution_parameter = 1e-6;
para.choose_receiver = 'index';
para.num_discretised_frequencies = 50;
para.emitter_index = 1;

% get the number of the receivers
num_receiver = size(receiver.positions, 2);

% get the index of the receiver
para.receiver_index = round(num_receiver/2);

% get the number of dimensions
dim = size(emitter.positions, 1);

% get the preferred distance of the chosen emitter-receiver pairs in terms
% of the radius of the detection ring
switch dim
    case 2
        para.distance_range = 2 * [0.93, 0.95];
    case 3
        para.distance_range = 2 * [0.85, 0.87];
end
   
% read additional arguments or overwrite default ones
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        % add to parameter struct (or overwrite fields)
        para.(varargin{input_index}) = varargin{input_index + 1};
    end
end

% display the number of emitters and receivers
disp(['The number of the chosen emitter is:' num2str(para.emitter_index)]);
disp(['The number of the chosen receiver is:' num2str(para.receiver_index)]);

% get the number of emitters
% num_emitter = size(emitter.positions, 2);

% the radius [m] of the detection surface (ring for 2D)
% detec_radius = norm(emitter.positions(:, 1));

% choose a range for the angular frequency [rad/s]
frequency_range = 2 * pi * [0, f_max];

% get the length of the time array
% num_time = length(time_array);
        
% get the time spacing
dt = time_array(2) - time_array(1);
        
% get the sampling frequency rate
% fs = 1/dt;
     
% get the angular frequency array supported by the grid for data
% simulation
omega = linspace(0, 2*pi * f_max, para.num_discretised_frequencies);
        
% get the angular frequencies in the chosen frequency range
omega = omega(omega > frequency_range(1) &...
            omega <= frequency_range(2));
        
% get the pressure source in the frequency domain
pressure_source = discreteFourierTransform(emitter.pulse,...
            omega, time_array);
   
% get the original pressure source     
pressure_source_original = pressure_source;

if para.deconvolve_source
   pressure_source = 1;
end

% compute the distance between the single emitter and receivers
distance_emitter_receivers = calculateDistanceEmitterReceiver(...
    emitter.positions(:, para.emitter_index), receiver.positions, []);

% get the time delays
time_delays = 1/sound_speed_water * distance_emitter_receivers;

% make the distances empty for 2D case, but keep it for computing Green's
% function for 3D case
if dim == 2
    distance_emitter_receivers = [];
end

%%=========================================================================
% GET THE SIMULATED AND APPROAXIMTED SIGNALS FOR THE CHOSEN EMITTER
%==========================================================================

% get the simulated pressure time series for the chosen emitter and on all
% the receivers
signal_measured_water_emitter = data_water(:, :, para.emitter_index);

% get the approximated pressure time series for the chosen emitter and on
% all the receivers. The pressure times series are approximated using the
% analytical Green's function.
signal_approximated_water_emitter = approxPressureGreens(pressure_source,...
    time_delays, 1, 0, 0,  omega, 1.4, distance_emitter_receivers, 'normal');

if para.deconvolve_source
    
    % deconvolve the measured time traces from the pressure
    % source in the frequency domain
    signal_measured_water_emitter = deconvolve(emitter.pulse, signal_measured_water_emitter,...
        time_array, omega, para.deconvolution_parameter);
    
else
    
    % calcualate the measured signals in the frequency domain (used as benchmark for testing the
    % anlytic Green's function)
    signal_measured_water_emitter = discreteFourierTransform(signal_measured_water_emitter,...
        omega, time_array);
end

%%=========================================================================
% GET THE COMPARISON RESULTS FOR A SINGLE EMITTER-RECEIVER PAIR
%==========================================================================
% get the signal simulated using the k-wave for the chosen receiver
signal_measured_water_receiver = signal_measured_water_emitter(para.receiver_index, :);

% get the signal approximated using the Green's function for the chosen receiver
signal_approximated_water_receiver = signal_approximated_water_emitter(para.receiver_index, :);

% calculate the relative discrepancy of the signals approximated using Green's
% function and the signals simulated using the k-Wave for the chosen
% emitter and all the receivers
relative_discrepancy_signal_water_emitter = vecnorm(signal_approximated_water_emitter...
    - signal_measured_water_emitter, 2, 2) ./...
    vecnorm(signal_measured_water_emitter, 2, 2) * 100;

% compute the relative error of the signals approximated using Green's
% function and the signals simulated using the k-Wave for a single
% emitter-receiver pair
relative_discrepancy_signal_water_receiver = norm(signal_approximated_water_receiver...
    - signal_measured_water_receiver)/norm(signal_measured_water_receiver) * 100;

% display the mean relative error for the signals for the chosen emitter
disp(['The mean relative error for the chosen emitter is:'...
    num2str(mean(relative_discrepancy_signal_water_emitter)) '%' ]);
% display the mean relative error for the signals for the chosen
% emitter-receiver pair
disp(['The relative error of the signal for the chosen emitter-receiver pair is:'...
    num2str(relative_discrepancy_signal_water_receiver) '%']);


%% ========================================================================
% PLOT THE SIGNALS
%==========================================================================
% normalise the excitation pulse
excitation_pulse_normalised = 1/max(abs(emitter.pulse)) * emitter.pulse;

% get the maximum of the pressure source
pressure_source_max = max(abs(pressure_source_original));

% display the normalised pressure source in the time domain
h1 = figure; plot(time_array(1:length(emitter.pulse)),...
    excitation_pulse_normalised);
xlabel('Time [s]');  ylabel('Amplitude [Pa m^{-2}]'); axis([0 1.2e-5 -1.1 1.1]);
set(gca, 'FontSize', 12);

% display the normalised amplitude of the pressure source in the frequency domain
h2 = figure;
subplot(2, 1, 1); plot(1/(2*pi) * omega, abs(1/ pressure_source_max * pressure_source_original));
xline(f_max,'--r', 'f_{max}','LabelVerticalAlignment', 'middle',...
    'LabelHorizontalAlignment', 'center');
xlabel('Frequency [Hz]'); ylabel('Amplitude [Pa m^{-2} s]');
axis([0, 1.1 * f_max, 0, +1.1]); set(gca, 'FontSize', 12);

% display the phase of the pressure source in the frequency domain
subplot(2, 1, 2); plot(1/(2*pi) * omega, 1/pi * angle(pressure_source_original) );
xline(f_max,'--r', 'f_{max}','LabelVerticalAlignment', 'middle',...
    'LabelHorizontalAlignment', 'center');
xlabel('Frequency [Hz]'); ylabel('Phase [\pi rad]');
axis([0, 1.1 * f_max, -1.1, +1.1]); set(gca, 'FontSize', 12);

% save the figure, if requested by the user
if para.save_plots 
    
    saveas(h1, [plot_directory, 'Fig_3a' '.fig' ]);
    saveas(h1, [plot_directory, 'Fig_3a' '.png' ]);
    saveas(h1, [plot_directory, 'Fig_3a' '.tiff' ]);
    saveas(h1, [plot_directory, 'Fig_3a' '.eps' ], 'epsc');
    
    saveas(h2, [plot_directory, 'Fig_3b' '.fig' ]);
    saveas(h2, [plot_directory, 'Fig_3b' '.png' ]);
    saveas(h2, [plot_directory, 'Fig_3b' '.tiff' ]);
    saveas(h2, [plot_directory, 'Fig_3b' '.eps' ], 'epsc');
    
end



%% ========================================================================
% PLOT THE SIMULATED/APPROXIMATED SIGNALS FOR THE CHOSEN EMITTER-RECEIVER
% PAIRS
%==========================================================================
if length(omega) > 1
    
    % display the amplitude of the chosen signals in water in the frequency
    % domain, simulated using k-Wave (green) and approximated using the analytic Green's
    % function (blue)
    h3 = figure;plot(1/(2*pi) * omega, 1/pressure_source_max * abs(signal_approximated_water_receiver), 'g-.',...
        1/(2*pi) * omega, 1/pressure_source_max * abs(signal_measured_water_receiver), 'b--');
    
    axis([0, 1.1 * f_max, 1.1/pressure_source_max * min(abs(signal_approximated_water_receiver)),...
        1.1/pressure_source_max * max(abs(signal_approximated_water_receiver))]);
    xlabel('Frequency [Hz]'); ylabel('Amplitude [Pa s]');
    legend('Greens', 'k-Wave');set(gca, 'FontSize', 12);
    
    
    ph1 = 1/pi * angle(signal_approximated_water_receiver);
    ph2 = 1/pi * angle(signal_measured_water_receiver);
    
    % find the frequencies at which the discrepancies include a jump at -pi
    % or +pi
    jump_indices =  find(abs(ph1-ph2)> 1);
    
    % correct the jumps at -pi or +pi
    ph2(jump_indices) = ph2(jump_indices) - 2 * sign(ph2(jump_indices));
    
    % display the phase of the chosen signals in water in the frequency
    % domain, simulated using k-wave (green) and approximated using the analytic
    % Green's function (blue)
    h4 = figure;plot(1/(2*pi) * omega, ph1, 'g-.',...
        1/(2*pi)* omega, ph2, 'b--');
    axis([0, 1.1 * f_max, -1.1, 1.1]);
    xlabel('Frequency [Hz]'); ylabel('Phase [\pi rad]');
    legend('Greens', 'k-Wave');set(gca, 'FontSize', 12);
   
    if 0 % para.save_plots
        
        saveas(h3, [plot_directory, 'Fig_3c' '.fig' ]);
        saveas(h3, [plot_directory, 'Fig_3c' '.png' ]);
        saveas(h3, [plot_directory, 'Fig_3c' '.tiff' ]);
        saveas(h3, [plot_directory, 'Fig_3c' '.eps' ], 'epsc');
        
        saveas(h4, [plot_directory, 'Fig_3d' '.fig' ]);
        saveas(h4, [plot_directory, 'Fig_3d' '.png' ]);
        saveas(h4, [plot_directory, 'Fig_3d' '.tiff' ]);
        saveas(h4, [plot_directory, 'Fig_3d' '.eps' ], 'epsc');
        
    end
    
    
end


end