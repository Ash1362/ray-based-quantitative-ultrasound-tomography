____________________________________________________________________________

                                  r-Wave

                    A MATLAB toolbox for high-resolution quantitative ultrasound
                     using two-point ray tracing 
____________________________________________________________________________

VERSION INFORMATION
____________________________________________________________________________

Version 1.0, Released 10th March 2023
Written by Ashkan Javaherian

Please report bugs and suggestions to email: a.javaherian@ucl.ac.uk, ashkan.javaherian@yahoo.com
% You will receive a response in 7 days.
The toolbox may be downloaded from https://github.com/Ash1362/ray-based-quantitative-ultrasound-tomography/

I would like to thank Dr. Yixuan Wu at university of Johns Hopkins for running the codes, and his useful feedbacks

____________________________________________________________________________

PRODUCT OVERVIEW
____________________________________________________________________________


____________________________________________________________________________

INSTALLATION INSTRUCTIONS
____________________________________________________________________________
run startup_simulation_ust.m
 
____________________________________________________________________________

RUNNING INSTRUCTIONS
____________________________________________________________________________
RUN ../simulation/examples-simulation.m

......................
NUMBER OF WORKERS
Set the number of workers for parallel programming by setting the varibale:
num_worker_pool = 16; (Default)

....................
DATA SIMULATION USING K-WAVE
For the first run, the user must simulate data using k-Wave by setting the variable:
do_data_sim = true;

For the next run, the user can load the data by setting:
do_data_sim = false;

Alternatively, the user can download the simulated UST data from ''www.zenodo.org''

.....................
PURPOSE OF RUNNING THE SCRIPT
For validation of ray tracing, the user must set the variales:
scenario = 'single_emitter';                            
purpose =  'raytracing_validation';        

For image reconstruction, the user must set the variales:
scenario = 'standard';                  
purpose = 'image_reconstruction'; 
......................
The scenarios for image reconstruction can be set by the variables:

The optimisation approach for image reconstruction
using the ray approaximation to heterogeneous Green's function
This can be set:
1-'backprojection', which was explained in:
Ashkan Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution quantitative ultrasound tomography❞, 
https://arxiv.org/abs/2211.00316, 2 Nov 2022.

2-'hessian', which was explained in:
A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering for biomedical ultrasound tomography❞,
Inverse Problems vol. 37, no.11, 115003, 2021.
.........................
The type of the absorption coefficient map used for image reconstruction
using the ray approaximation to heterogeneous Green's function. This can
be set 'true', 'homogeneous', or 'none'

....................
  
The signal-to-noise ratio of the simulated ultrasound data
This can be set 40 dB, 30 dB, or 25 dB.
noise_level = 40; (default)


Running the script gives the results in:
Ashkan Javaherian, Hessian-inversion-free ray-born inversion for high-resolution quantitative ultrasound tomography, 
https://arxiv.org/abs/2211.00316, 2 Nov 2022.

Please read the section ❝Numerical results❞ in that preprint for better understanding of the scenarios in the paper.
    
____________________________________________________________________________

RELEASE NOTES
____________________________________________________________________________


_______________________________________________________________

LICENSE
____________________________________________________________________________

r-Wave (c) 2022 Ashkan Javaherian
____________________________________________________________________________
