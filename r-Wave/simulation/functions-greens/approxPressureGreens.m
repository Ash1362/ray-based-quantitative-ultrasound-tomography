function [signal_delayed] = approxPressureGreens(signal, time_delays, attenuation,...
    absorption, caustic_number, omega, y, greens_method, source_mode)
%APPROXPRESSUREGREENS approximates the pressure field using an analytic Green's
%function
%
% DESCRIPTION:
% approxPressureGreens approximates the pressure field using Green's function
% given a time-varying pressure source

% USAGE:
%
% INPUTS:
%       signal             - a time varying pressure source(for the 'forward' field),
%                            or a time series data in the residual (for the 'adjoint'
%                            field)
%       time_delays        - a vector of the time delays on the grid points
%                            with respect to the emission point 
%       attenuation        - a vector of the refraction-induced portion of 
%                            attenuation due to geometrical spreading on the
%                            grid points with repect to the emission point 
%       absorption         - a vector of the accumulated acoustic
%                            absorption
%       caustic_number     - a vector of the integer number of caustic on the
%                            grid points along the rays 
%       omega              - angular frequency array [rad/s]
%       y                  - the power law exponent factor
%       greens_method      - the method for apprpoximating the pressure
%                            field using Green's function
%       source_mode        - the source mode for computing the pressure
%                            field. This can be 'normal' or 'time-reversal'
% 
%
% OPTIONAL INPUTS:
%
%
%
% OUTPUTS:
%       signal_delayed     - the delayed signals

% ABOUT:
%       author            - Ashkan Javaherian
%       date              - 25.07.2020
%       last update       - 25.07.2020
%


% calculate the product of time delays by the angular frequencies
time_delays = time_delays * omega;
switch greens_method
    
    case 'analytic'
    
    switch source_mode
        
        case 'normal'
            
            % apply the geometrical spreading, exponential decay due to the physical
            % absorption, time delay and pi/2 phase shifts due to caustics
            signal_delayed =  attenuation .* exp ( - absorption * abs(omega).^y ) .* ...
                (1/4) .* sqrt( 2./ (pi * time_delays ) ) .*...
                exp(1j * (time_delays + pi/4 - pi/2 * caustic_number) ) .* signal;
            
        case 'time-reversal'
            
            % apply the geometrical spreading, exponential decay due to the physical
            % absorption, time delay and pi/2 phase shifts due to caustics
            signal_delayed =  1./ ( attenuation .* exp ( - absorption * abs(omega).^y ) .* ...
                (1/4) .* sqrt( 2./ (pi * time_delays ))  ) .*...
                exp(-1j * (time_delays + pi/4 - pi/2 * caustic_number) ) .* signal;
            
    end
        
        
    otherwise
        
        error(['The chosen method for approximating the pressure using Greens'...
            'function is not supported yet, or deprecated.']);
        
end


end
