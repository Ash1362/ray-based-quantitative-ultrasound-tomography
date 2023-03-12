____________________________________________________________________________

                                  r-Wave

                    A MATLAB toolbox for high-resolution and quantitative ulrasound tomography
                     using two-point ray tracing 
____________________________________________________________________________

VERSION INFORMATION
____________________________________________________________________________

Version 1.0, Released 10th March 2023
Written by Ashkan Javaherian

Please report bugs and suggestions to email: a.javaherian@ucl.ac.uk, ashkan.javaherian@yahoo.com
You will receive a response in 7 days.
The toolbox may be downloaded from https://github.com/Ash1362/ray-based-quantitative-ultrasound-tomography/ 

____________________________________________________________________________

Acknowledgment
____________________________________________________________________________
The project was (only financially) supported by :
1) European Commission: PAMMOTH - Photoacoustic/Ultrasound Mammoscopy for evaluating screening-detected abnormalities in the breast (732411)
2) Research Councils UK: WHOLE-BODY, HIGH RESOLUTION, 3D, SMALL ANIMAL PHOTOACOUSTIC AND ULTRASOUND COMPUTED TOMOGRAPHY SYSTEM (EP/T014369/1)

Regarding technical points of view, since 2021, when the code developer raised some issues noticed, the contribution of the code developer to these projects was completely stopped by the associated department at University College London, and the code developer has been working in complete isolation (only few contacts via email) until the end of his appointment with UCL.  
____________________________________________________________________________
The examples include the scenarios in the papers:
1 - A. Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution quantitative ultrasound tomography❞, 2022, https://arxiv.org/abs/2211.00316/ .
2 - A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering for biomedical ultrasound tomography❞, Inverse Problems vol. 37, no.11, 115003, 2021.  https://iopscience.iop.org/article/10.1088/1361-6420/ac28ed/
3- A. Javaherian, F. Lucka and B. T. Cox, ❝Refraction-corrected ray-based inversion for three-dimensional ultrasound tomography of the breast❞, Inverse Problems, 36 125010.  https://iopscience.iop.org/article/10.1088/1361-6420/abc0fc/  

The studies [1] and [2] correspond to high-resolution and quantitative reconstruction of the sound speed of the breast from synthetic
ultrasound data simulated using the k-Wave toolbox (www.k-wave.org/). The transducers are assumed points which are placed on a 2D ring,
and the image reconstruction is done in 2D, but the approaches can be extended to 3D.

The study [3] corresponds to a full-3D and refraction-corrected ray-based algorithm for reconstructing a qunatitatively accurate and
low-resolution image of the sound speed from time-of-flight data. This study was sucessfully implemented on in-vivo ultrasound data
measured from the full-3D pammoth system. The application of [3] on Pammoth system shows significant contrast between the bent-ray approach
proposed in [3] and straight-ray approach for a fast and quantitative reconstruction of the sound speed from time-of-flight-data.

_________________________________________________________________

PRODUCT OVERVIEW
____________________________________________________________________________


____________________________________________________________________________

Getting started:
____________________________________________________________________________

1) Get access to the UST data
Download the folder ''data_ust_kWave_transmission'' via the link:
https://doi.org/10.5281/zenodo.7717290/ [7] 

Add the folder ''data_ust_kWave_transmission'' to the path:
.../r-Wave/simulation/data/...

This folder includes the ultrasound data simulated by the k-Wave toolbox. Alternatively, the user can set ''data_sim =true;'' in the
example scripts. By setting this parameter true, the ultrasound data are simulated and are stored in the associated path.

By setting ''data_sim=true;'', the following point must be considered. 
For the results in the papers [1] and [2], the transducers are assumed off-grid points which are placed on a 2D ring.
In the k-Wave, for interpolation of the simulated pressure field between the grid points and an off-grid point representing a tranducer,
an optional prameter 'BLITolerance', which is in the k-Wave by default 0.1 or 0.05, must be set. This parameter adjusts the portion of the
grid points included in the interpolation to an off-grid point, and is used for reducing the computational cost of the interpolation.
However, because using a k-space pseudospectal approach, the gradient at each point is computed using the information on all grid points,
all the grid points must be contributed to an interpolation to an off-grid point, when the same approach [6] is taken for interpolation.
But because including all the grid points for an interpolation to an off-grid point is very costly, 'BLITTolerance' is used for confining 
the interpolation to the grid points close to the off-grid point by setting the deafult 'BLITolerance' 0.1 or 0.05. Based on my experince,
these values for 'BLITolerance' are large and will deteriorate the first arrival of the signals.
For preserving the information in the simulated  pressure field including the first arrival of the signals, 'BLITolerance' was here set 0.001.
Note that smaller 'BLITolerance' means inclusion of larger number of grid points, and therefore increases the required memory for the k-Wave simulation, but it will make the k-Wave simulations accurate.

It is reminded that that in [1] and [2], early iterations of a time-of-flight-based image reconstruction algorithm using the first arrival of the signals is used for providing initial guess for the Green's inversion approaches. In study [3], which fully corresponds to an image reconstruction using time-of-flight data in full-3D geomtery, the interpolation is done using a neighboring approach.


2) A digital breast phantom developed by Mark Anastasio's group is used in this project. This phantom must be downloaded via the link:'...
https://anastasio.bioengineering.illinois.edu/downloadable-content/oa-breast-database/ [4].
The folder 'Neg_47_Left' is used for simulation of UST data for all studies in this toolbox.
The folder 'Neg_47_Left' must be added to the directory:
'.../simulation/data/phantom/OA-BREAST/Neg_47_Left/...'

3) Run ''startup_simulation_ust.m''

4) Run examples in the path ''...r-Wave/simulation/examples-simulation/....''
____________________________________________________________________________

RUNNING EXAMPLES
____________________________________________________________________________
1) .../simulation/example_2D_validate_greens.m  (paper [1])
This example validates ray approximation to heterogeneous Green's function in approximating phase and amplitude via
a comparison with the phase and amplitude simulated by the k-Wave.

2) .../simulation/example_2D_image_tof_greens.m (paper [1])
This example reconstructs image of the sound speed using inversion approaches based on ray approaximation to heterogeneous Green's function [1,2]. Early iterations of an image reconstruction algorithm using time-of-flight data will be used for providing initial guess for the Green's inversion approaches. Note that the Green's inversion approach used in [2] was used in [1] as the benchmark.

3) .../simulation/example_3D_image_tof.m (paper [3])
This example reconstructs a quantitatively accurate, low resolution and volumetric (full-3D) image of the sound speed
from time-of-flight of ultrasound data simulated on a hemispherical detection surface. The user can compare the images
reconstructed using the bent-ray inversion approach [3] with the same image reconstructed using straight rays.  
('matrix_construction_method= 'bent-ray'), or ('matrix_construction_method= 'straight-ray')

4) .../simulation/example_fish_eye.m (paper [3])
This example tests ray tracing algorithms on a ❝Maxwell's Fish-eye lens❞ phantom.
All the used ray tracing algorithms in this examplem can be used for time-of-flight-based reconstruction of the sound speed,
but only 'Rung-kutta-2nd' can be used for ray tracing for the Green's inversion approach. 

Please read the description of examples for better understanding the examples!
____________________________________________________________________________

Simulation of ultrasound data using K-WAVE
____________________________________________________________________________
The pressure field used as the benchmark is simulated using the k-Wave toolbox.  www.k-Wave.org (v. 1.3.) [5].
The k-Wave toolbox and the functions for an off-grid interpolation is available in this project. 
..............................................................................
Notation: It was noticed that correction steps described in [8] must be applied on the k-Wave, either v. 1.3 or v. 1.4,
for getting accurate signals matching the analytic Greens' formula in homogenous media (only water).
Here, instead of applying the required corrections to the k-Wave, the pressure time series simulated by the k-Wave are kept
uncahnged, and an inverse of the required correction steps was enforced on the emission pulse which is used as the input to
the Green's formula in homogeneous and heterogeneous media. Taking both approaches are the same, and it will affect only our 
assumption about the emission pulse. (please contact me if you need more information!)

____________________________________________________________________________

REFERENCES
__________________________________________________________________________

1 - A. Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution quantitative ultrasound tomography❞, 2022. https://arxiv.org/abs/2211.00316/ .
2 - A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering for biomedical ultrasound tomography❞, Inverse Problems vol. 37, no.11, 115003, 2021. https://iopscience.iop.org/article/10.1088/1361-6420/ac28ed/ 
3- A. Javaherian, F. Lucka and B. T. Cox, ❝Refraction-corrected ray-based inversion for three-dimensional ultrasound tomography of the breast❞, Inverse Problems, 36 125010. https://iopscience.iop.org/article/10.1088/1361-6420/abc0fc/  
4- Y. Lou, W. Zhou, T. P. Matthews, C. M. Appleton and M. A. Anastasio, ❝Generation of anatomically realistic numerical phantoms for photoacoustic and ultrasonic breast imaging❞, J. Biomed. Opt., vol. 22, no. 4, pp. 041015, 2017. 
https://anastasio.bioengineering.illinois.edu/downloadable-content/oa-breast-database/
5- B. E. Treeby and B. T. Cox, ❝k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave fields❞, J. Biomed. Opt. vol. 15, no. 2, 021314, 2010. http://www.k-wave.org/ 
6- E. S. Wise, B. T. Cox, J. Jaros, B. E. Treeby, ❝Representing arbitrary acoustic source and sensor distributions in Fourier collocation methods❞, J. Acoust. Soc. of Am., vol. 146, no. 1, pp. 278-288, 2019.
7- A. Javaherian, 2023, ❝Transmission ultrasound data simulated using the k-Wave toolbox as a benchmark for biomedical quantitative ultrasound tomography using a ray approximation to Green's function❞ (1.0) [Data set]. Zenodo. 
https://doi.org/10.5281/zenodo.7717290.
8- A Javaherian, ❝A note on an open-source toolbox for simulation of acoustic waves: inclusion of time-varying source❝, https://arxiv.org/abs/2212.04466.



RELEASE NOTES
____________________________________________________________________________


_______________________________________________________________

LICENSE
____________________________________________________________________________

r-Wave (c) 2022 Ashkan Javaherian
____________________________________________________________________________
