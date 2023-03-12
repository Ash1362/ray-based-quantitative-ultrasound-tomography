function [s] = excitUSCTPulse(kgrid, sound_speed, pulse_prop, varargin)
%EXITUSCTPULSE simulates an excitation pulse for USCT imaging
%
% DESCRIPTION:
%    exitUSCTPulse simulates an excitation pulse for USCT imaging.
%    The supported pulses are 'Impulse', 'three-cycle_sinusoid' or 'Chirp',
%    'Pammoth_1'
%    A 'Three-cycle-sinusoid' pulse is frequently used for full-wave invesrion for USCT.
%    (See T.P. Matthews and M.A. Anastasio 2017 Joint reconstruction of the initial pressure and speed of sound distributions from combined photoacoustic and ultrasound tomography
%    measurements, Inverse Problems 33 124002, section 4.1., Eq. 17)
%    A 'Chirp' excitation pulse is used by Karlsruhe Institute (KIT)
%    group because of a good compressibily after cross-correlation with a delayed-scaled pulse of the same shape.
%    Cross-correlation is used in the so-called adapted-matched filtering method for reflectivity imaging.
%    See N.V. Ruiter, M. Zapf, R. Dapp, T. Hopp, W.A. Kaiser and H. Gemmeke, First Results of a Clinical Study with
%    3D Ultrasound Computer Tomography, Joint UFFC, EFTF and PFM Symposium,
%    2013, DOI: 10.1109/ULTSYM.2013.0168.
%    See also Improvement of 3D Ultrasound Computer Tomography Images by Signal
%    Pre-Processing, N.V. Ruiter, G.F. Schwarzenberg, M. Zapf, H. Gemmeke, IEEE International Ultrasonics Symposium Proceedings,
%    2008, DOI: 10.1109/ULTSYM.2008.0205.

% USAGE:
%     
       
% INPUTS:
%       kgrid             - the computational grid for k-Wave simulation
%       sound_speed       - sound speed [m/s] distribution
%       pulse_prop.fc     - central frquency of the excitation pulse
%       pulse_prop.b      - the frequency bandwidth for the 'Chirp' pulse,
%                           or the temporal standard deviation of a Gaussian window
%                           for the 'Three-cycle-sinusoid' pulse
%       pulse_prop.T      - pulse duration for the 'Chirp' pulse and the temporal mean of 
%                           a Gaussian window for the 'Three-cycle-sinusoid' pulse
%       pulse_prop.t0     - the time for the peak of the approximated delta
%                           function (used if the excitation pulse is
%                           either 'Impulse-Dirichlet' or
%                           'Impulse-Gaussian')
% 
% OPTIONAL INPUTS:
%       Excit             - the excitation pulse that can be 'Impulse', 'Three-cycle_sinusoid' or 'Chirp',
%                          'Pammoth_1' (Default = 'Pammoth_1')
%       Low_Filter        - boolean controlling whether a low-pass filter 
%                           is used for cutting-off the higher frequencies,
%                           which are not supported by the computational
%                           grid (Default = 'true')
%       Plot              - Boolean controlling whether the simulated settings are
%                           plotted (default = true)

% OUTPUTS:
%        s                - excitation pulse in time 


% ABOUT:
%       author            - Ashkan Javaherian
%       date              - 29.12.2019
%       last update       - 29.12.2019
%
% This script is part of the r-Wave Tool-box 
% Copyright (c) 2022 Ashkan Javaherian 



para.Excit          = 'Pammoth_1'; 
para.Low_Filter     = true;
para.Low_Filter_Fac = 1;
para.Plot           = true;
para.Plot_frequ_maxaxis = 2e6;

% add to parameter struct (or overwrite fields)
if(~isempty(varargin))
    for input_index = 1:2:length(varargin)
        para.(varargin{input_index}) = varargin{input_index+1};
    end
end


% get the time spacing
dt = kgrid.dt;

switch para.Excit
case{'Impulse-Dirichlet','Impulse-Gaussian'}
   t = [-flip(kgrid.t_array(2:end)), kgrid.t_array];
    otherwise
   t = kgrid.t_array;
end



if isfield(pulse_prop,'T')
    T = pulse_prop.T;
end
if isfield(pulse_prop,'b')
    b = pulse_prop.b;
end
if isfield(pulse_prop,'fc')
    fc = pulse_prop.fc;
end
if isfield(pulse_prop,'t0')
    t0 = pulse_prop.t0;
end
    

switch para.Excit
    case 'Impulse'
        s = zeros(1, kgrid.Nt);
        s(1) = 1;
        
    case 'Impulse-Dirichlet'
        s = 1./(pi.*(t-t0)).*sin((t-t0)/T);
        
    case 'Impulse-Gaussian'
         s = 1/(sqrt(pi)*T) * exp(-((t-t0).^2)/(T^2));
        
    case 'Three-cycle-sinusoid'
        s = exp(-((t-T).^2)./(2 * b^2)) .* sin(2*pi*fc*t);
        
    case 'Chirp'
     f = fc - (b/2) + (b/(2*T)) .* t;
     s = sin(pi .* t/T) .* sin(2*pi.*f .* t);
   % s( t > T ) = 0; % Should zero padded to the signal
 
    case 'Pammoth_1'
s = [ -0.001818578
-0.001403717
-0.00091776
-0.000617935
-0.000653785
-0.001029553
-0.001587733
-0.001966936
-0.001894365
-0.001461619
-0.000901847
-0.000496246
-0.000483088
-0.000959491
-0.001631733
-0.002051989
-0.002018017
-0.001539803
-0.00082472
-0.000236848
-0.000285433
-0.000871016
-0.001640567
-0.002199404
-0.002239828
-0.001653313
-0.000681353
-2.01679E-05
-2.49324E-05
-0.000673031
-0.001660513
-0.002497338
-0.002616752
-0.001747195
-0.000551586
0.000341658
0.000470012
-0.000319664
-0.00176166
-0.00292875
-0.002999513
-0.001987487
-0.000392542
0.0010129
0.001415494
0.000209518
-0.001907909
-0.003581567
-0.003930429
-0.002658335
-0.000166528
0.002364409
0.002826418
0.000994152
-0.002201128
-0.005204129
-0.006374284
-0.004536855
0.000193841
0.004545696
0.00589809
0.00316279
-0.003010431
-0.01019169
-0.013774041
-0.009382517
0.000814728
0.011974517
0.017919905
0.01300357
-0.005789137
-0.032832913
-0.048228462
-0.034250435
0.022777286
0.128840644
0.279139256
0.450973318
0.594032432
0.670002309
0.652385585
0.528866149
0.303760885
-0.001839288
-0.332128247
-0.628158524
-0.852642432
-0.980447485
-1
-0.913944443
-0.728249662
-0.495177813
-0.244807084
-0.003435365
0.206395927
0.366952596
0.457732251
0.478788438
0.445049864
0.370558877
0.27272264
0.170450382
0.083455212
0.031188804
0.01096647
0.016704541
0.038527807
0.064764921
0.083997217
0.084061752
0.066236199
0.035601195
-0.001660854
-0.038921532
-0.070137706
-0.089422167
-0.094779575
-0.088669707
-0.073193584
-0.051179012
-0.026023009
-0.001747409
0.015311012
0.022336489
0.018001181
0.002881106
-0.020557545
-0.048083679
-0.072433371
-0.087106914
-0.089720492
-0.079998193
-0.059644521
-0.031994366
-0.001736077
0.025275459
0.046377975
0.06049814
0.067541265
0.06816233
0.063472208
0.054526338
0.043892633
0.032975856
0.0230974
0.015373757
0.010588246
0.009310668
0.010136152
0.011317247
0.011374643
0.009155384
0.004048388
-0.003896186
-0.013156
-0.021667088
-0.028091788
-0.031555105
-0.031715211
-0.02879299
-0.023310488
-0.01692382
-0.010441937
-0.004477575
0.000515696
0.004225487
0.006234369
0.006568864
0.00558325
0.003678212
0.001396909
-0.000632351
-0.001697418
-0.001121965
0.000805049
0.003683724
0.00695496
0.009999817
0.012242226
0.0130274
0.012567989
0.011171994
0.009179313
0.006937482
0.004751685
0.002885298
0.001445657
0.000385132
-0.00032951
-0.000746059
-0.000925935
-0.000946473
-0.000951653
-0.001095047
-0.001464848
-0.002091372
-0.002935824
-0.003883521
-0.004671224
-0.005040107
-0.004879205
-0.004174229
-0.003009098
-0.001553517
-6.2472E-05
0.0011707
0.002065522
0.002637461
0.002962651
0.003149945
0.003308575
0.003499056
0.003669873
0.003748121
0.00367532
0.003423466
0.003006492
0.002487166
0.001979741
0.001533458
0.001157123
0.000836307
0.000544207
0.000258417
-7.30871E-06
-0.000218213
-0.000347711
-0.000387048
-0.00034767
-0.000260885
-0.000173468
-0.000105252
-4.40859E-05
2.70415E-05
0.00011801
0.000221584
0.000288296
0.000244243
6.58388E-05
-0.000230164
-0.000588576
-0.000924934
-0.00112833
-0.001094837
-0.000870795
-0.000521359
-0.000128132
0.000228669
0.000487458
0.000602286
0.000617608
0.00056051
0.000446423
0.000284636
8.21291E-05
-0.000143238
-0.000330979
-0.000422107
-0.000377074
-0.000185634
0.000127216
0.000494725
0.000782101
0.000906422
0.000837157
0.000580556
0.000171978
-0.000333877
-0.000874464
-0.001395425
-0.001874599
-0.002292679
-0.002630277
-0.002867998
-0.002978186
-0.002970005
-0.002899556
-0.002817528
-0.002770897
-0.002793189
-0.002895589
-0.003050508
-0.003191459
-0.003281763
-0.003302064
-0.003249457
-0.003135072
-0.00297144
-0.002772769
-0.002539461
-0.002272757
-0.001982024
-0.001686732
-0.001421769
-0.001225241
-0.001088246
-0.000984611
-0.000877049
-0.000726902
-0.000505063
-0.000199963
0.000130854
0.000447363
0.000721209
0.000939086
0.001104033
0.001228949
0.001327532
0.001390983
0.001401401
0.001344388
0.001216206
0.001031451
0.000842744
0.000696949
0.000609477
0.000573673
0.000565092
0.000547951
0.000480109
0.000372275
0.000252652
0.00015013
8.59664E-05
6.69408E-05
7.91539E-05
8.61389E-05
7.42307E-05
5.09369E-05
3.55812E-05
5.23578E-05
0.000122791
0.00024343
0.000369224
0.000468038
0.000521054
0.000526807
0.000503185
0.000491984
0.000521793
0.000587773
0.000670868
0.000746483
0.00079265
0.000798709
0.000787881
0.00078894
0.000813085
0.000856488
0.000901723
0.00092036
0.00087895
0.000794842
0.000693924
0.000603727
0.000545216
0.000525676
0.00053011
0.000520766
0.000482145
0.00041682
0.000340136
0.000275152
0.00024737
0.000259132
0.000282483
0.00029482
0.000282324
0.000243972
0.000194331
0.00017025
0.000186077
0.000234203
0.000295038
0.000344635
0.000362256
0.000338512
0.00030115
0.000277748
0.000281604
0.000312213
0.000355576
0.000384694
0.00037285
0.000330527
0.000277909
0.000236906
0.000223963
0.000243128
0.000275276
0.000286028
0.000262946
0.000209475
0.000141401
8.22916E-05
5.92366E-05
7.05518E-05
9.1092E-05
9.87179E-05
7.89588E-05
2.95614E-05]';
   
% get the upsampling rate for the chosen CFL
 upsampling_rate = pulse_prop.dt/dt;

% get the time array associated with the given excitation pulse
 t_array_pulse = cumsum([0, pulse_prop.dt * ones(1, length(s)-1)]);

% get the time array for the upsampled excitation pulse using the determined
% upsampling rate
 t_array_pulse_upsampled = cumsum([0, dt * ones(1, ceil(upsampling_rate * length(s))-1)]);

% interpolate the excitation pulse to the upsampled time array 
 s = interp1(t_array_pulse, s, t_array_pulse_upsampled, 'spline'); 


    otherwise       
end

 

% Apply a low-pass filter on the excitation pulse such the maximum frequency 
% of the excitation pulse is supported by the computational grid.
if para.Low_Filter
     s = filterTimeSeries(kgrid, struct('sound_speed', para.Low_Filter_Fac * sound_speed), s, 'ZeroPhase', true);
end

% get the pulse in frequency domain
[freq, freq_spectrum] = spect(s, 1/dt);

switch para.Excit
case{'Impulse-Dirichlet','Impulse-Gaussian'}
   s(t<0) = [];
    otherwise
end


if para.Plot
    figure;subplot(2,1,1); plot(kgrid.t_array(1:size(s,2)), s); xlabel('Time [s]');
    ylabel('Pressure amplitude');title('Excitation pulse');
    subplot(2,1,2); plot(freq, abs(freq_spectrum)); xlabel('Frequency [Hz]');
    ylabel('Amplitude spectrum');title('Frequency spectrum');
    axis([0, para.Plot_frequ_maxaxis, 0, max(abs(freq_spectrum))]);
end





end

