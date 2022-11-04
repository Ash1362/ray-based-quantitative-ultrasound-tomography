function [tof , AIC ] = araic( sig , dt ,  t_array , T , nw_fac ,   nw_fac_back , ther , M, T_min, T_max )
%   ARAIC calculates the first arrival of the signal using the ARAIC method
%   
%
% DESCRIPTION:
%      araic Calculates the first-arrival of the signal using the ARAIC method
%
% USAGE:
%     
%
% INPUTS:
%       sig               - a 1D signal associated with a pair of emitter
%                          and receivers  
%       dt                - measurement separation in time
%       T                 - the temporal width of the excitation pulse  
%       beta              - a smoothing parameter, which is used for
%                           mitigating the rapid fluctuations of the
%                           calculated attribute                      
%      T_min              - a minumum sound speed for the medium
%      T_max              - a maximum sound speed for the medium
 
% OPTIONAL INPUTS:
%
 
% OUTPUTS:
%       tof               - the calculated first-break arrival in time [s]
%       ER                - the attribute function
%       ER_EPS            - the atribute function after EPS filtering
%       above
% % ABOUT:
%       author          - Ashkan Javaherian
%       date            - 18.03.2020
%       last update     - 18.03.2020
%
% This script is part of the r-Wave Tool-box (http://www.r-wave.org).
% Copyright (c) 2020 Ashkan Javaherian and Ben Cox
 
nw =  round ( nw_fac*T/dt); 

  
%% normalise the signal with respect to the maximum absolute amplitude
sig  =  sig / max( abs(sig) );
N = length( t_array );
%% truncate the normalised signal with a min and max speed of sound
int_st   = round( T_min/dt );
int_end  = round( T_max/dt );


%% calculate the envelope of the normalised truncated signal
env_sig = abs(hilbert( sig )) ;
env_sig = env_sig(int_st : int_end);

%% find the first index exceeding the threshold
[~ , ixm ] = find( env_sig/max( env_sig ) > ther , 1 , 'first' );
%% shift the index for accounting for the truncated part
ixm = ixm + int_st - 1;


%% updtae the window using the envekope function
int_st  = max ( M + 1 , ixm - round( nw_fac_back*T/dt )); 
int_end = min ( ixm , N - M ); 
  
  
  
  
% select the part of the signal within the window 
  %sig      =   sig  (int_st:int_end);
  t_arrayc = t_array( int_st : int_end ); 
  
  
  
 
 
AIC = zeros( size( t_arrayc ) );
Nt  = length ( t_arrayc );
 
for i = 1 : Nt
    k = i + int_st - 1;
 [ ~ , e_w1 ] = aryule(        sig( 1 : k      )    , M  );
 [ ~ , e_w2 ] = aryule( flip ( sig( k + 1 : N  ) )  , M  );
AIC(i) = ( k - M )*log( e_w1 ) + ( N - M - k )*log( e_w2 );
end
 
%%
[ ~, ix_last]  = find ( diff(AIC) > 0 , 1 , 'last' );  
t_arrayc = t_arrayc ( 1 : ix_last );
AIC = AIC ( 1 : ix_last );
Nt = length ( t_arrayc );



%% calculate the minimal AIC value within the time window 
[ AICm , ixm ] = min ( AIC ); 


 IXl =  max( 1 , ixm - nw ): min( Nt , ixm + nw );

%% calculate the discrepancy of AIC for the data samples within the
%% selected window and the minimal value
Delta = AIC( IXl ) - AICm;

%% calculate the Akaike weight for each data sample within the time window
 Wv = exp ( - Delta/2 );
 W = Wv/ sum( Wv );
 
%% calculate a weighted average TOF using the Akaike weights
 tof = W*t_arrayc( IXl )';
 
 
 
end
