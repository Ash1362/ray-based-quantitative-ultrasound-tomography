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

For this version, the codes include the scenarios in the papers:
[1] Ashkan Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution quantitative ultrasound tomography❞, 
https://arxiv.org/abs/2211.00316, 2 Nov 2022.
[2] A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering for biomedical ultrasound tomography❞,
Inverse Problems vol. 37, no.11, 115003, 2021.
[3] A. Javaherian, F. Lucka and B. Cox, ❝Refraction-corrected ray-based inversion for three-dimensional ultrasound
tomography of the breast❞, Inverse Problems, vol. 36, 125010, 2020.

The studies [1] and [2] correspond to high-resolution quantitative ultrasound tomography of the breast from ultrasound data measured
using the transducers placed on a ring, and the image reconstruction is done in 2D. The studies [1] and [2] have not been tested on experimental UST data yet. The plan is to collaborate with University of Johns Hopkins for testing and implentation on experimental UST data. 

The study [3] corresponds to a full-3D refraction-corrected ray-based algorithm for reconstructing a low-resolution image of the sound speed from time-of-flight data. This study has been implemented on ultrasound data measured using the full-3D pammoth system. The application of [3] 
on experimental data shows significant contrast between the bent-ray approach proposed in [3] and straight-ray approach for quantitativce reconstruction of the sound speed from time-of-flight-data.

....................................................................
Acknowledgement:
I would like to thank Dr. Yixuan Wu at university of Johns Hopkins for useful comments for improving the usability of the codes, especially for implementation on experimental data. (https://yixuanwu.page/) 
Some comments will be applied on the next version (1.1), which will be commited in few weeks.
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
run ../simulation/examples-simulation.m

......................
NUMBER OF WORKERS
Set the number of workers for parallel programming by setting the varibale:
num_worker_pool = 16; (Default)

....................
DATA SIMULATION USING K-WAVE
For the first run, the user can simulate data using k-Wave by setting the variable:
do_data_sim = true;

For the next run, the user can load the data by setting:
do_data_sim = false;

Alternatively, the user can download the simulated UST data from ''www.zenodo.org'',
and always set do_data_sim = false;

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
1) greens_optimisation-approach, (Default: 'backprojection')
2) absorption_map, (Default: 'true')
3) noise _level, (Default: 40)


The optimisation approach for image reconstruction
using the ray approaximation to heterogeneous Green's function
This can be set:
1) greens_optimisation-approach = 'backprojection';

This image reconstruction approach was explained in:
Ashkan Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution quantitative ultrasound tomography❞, 
https://arxiv.org/abs/2211.00316, 2 Nov 2022.

2) greens_optimisation-approach = 'hessian';

This image reconstruction aproach was explained in:
A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering for biomedical ultrasound tomography❞,
Inverse Problems vol. 37, no.11, 115003, 2021.
.........................
The type of the absorption coefficient map used for image reconstruction
using the ray approaximation to heterogeneous Green's function. This can
be set:
1) aborption_map = 'true'; 
2) absorption_map ='homogeneous';
3) absorption_map = 'none';

....................
  
The signal-to-noise ratio of the simulated ultrasound data
This can be set
1) noise_level = 40 dB;
2) noise_level = 30 dB;
3) noise_level = 25 dB;


Running the script gives the results in the last paper::

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
