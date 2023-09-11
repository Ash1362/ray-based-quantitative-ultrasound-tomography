____________________________________________________________________________

                                  r-Wave

                    A MATLAB toolbox for quantitative ulrasound tomography
                    using two-point ray tracing 
____________________________________________________________________________

VERSION INFORMATION
____________________________________________________________________________

Version 1.1, Released 3rd September 2023
Written by Ashkan Javaherian

Please report bugs and suggestions to email: ashkan.javaherian@yahoo.com
The toolbox may be downloaded from https://github.com/Ash1362/ray-based-quantitative-ultrasound-tomography/ 

____________________________________________________________________________

ACKNOWLEDGEMENT
____________________________________________________________________________
_______________________________________________________________________
The examples include the scenarios in the papers:
1 - A. Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution quantitative ultrasound tomography❞, 2022,
https://arxiv.org/abs/2211.00316/.
2 - A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering for biomedical ultrasound tomography❞, Inverse Problems vol. 37, no.11, 115003, 2021.
https://iopscience.iop.org/article/10.1088/1361-6420/ac28ed/.
3- A. Javaherian, F. Lucka and B. T. Cox, ❝Refraction-corrected ray-based inversion for three-dimensional ultrasound tomography of the breast❞, Inverse Problems, 36 125010.
https://iopscience.iop.org/article/10.1088/1361-6420/abc0fc/.  

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

GETTING STARTED 
____________________________________________________________________________
____________________________________________________________________________

A) SIMULATION STUDIES
____________________________________________________________________________

1) Get access to the UST data
Download the folder ''simulation.zip'' via the link:
https://zenodo.org/record/8330926 [7] 

Add the folder ''data_ust_kWave_transmission'' to the path:
.../r-Wave/data/simulation/...

This folder includes the ultrasound data simulated by the k-Wave toolbox. Alternatively, the user can set ''data_sim = true;'' in the
example scripts. By setting this parameter true, the ultrasound data are simulated and are stored in the associated path.

By setting ''data_sim=true;'', the following point must be considered. 
For the results in the papers [1] and [2], the transducers are assumed off-grid points which are placed on a 2D ring.
In the k-Wave, for interpolation of the simulated pressure field between the grid points and an off-grid point representing a tranducer,
an optional prameter 'BLITolerance', which is in the k-Wave by default 0.1 or 0.05, must be set. Based on my experince,
these values for 'BLITolerance' are large and will deteriorate the first arrival of the signals. For preserving the information in the simulated
pressure field including the first arrival of the signals, 'BLITolerance' was here set 0.001 for the 2D case.
Note that smaller 'BLITolerance' means inclusion of larger number of grid points for interpolation onto the trasducers, and therefore an increase
in the required memory for the k-Wave simulation, but it will make the k-Wave simulations accurate.

It must be reminded that that in [1] and [2], early iterations of a time-of-flight-based image reconstruction algorithm using the first arrival of the
signals is used for providing initial guess for the Green's inversion approaches. In study [3], which fully corresponds to an image reconstruction using
time-of-flight data in full-3D geomtery, the interpolation is done using a neighboring approach.


2) A digital breast phantom developed by Mark Anastasio's group is used in this project. This phantom must be downloaded via the link:'...
https://anastasio.bioengineering.illinois.edu/downloadable-content/oa-breast-database/ [4].
The folder 'Neg_47_Left' is used for simulation of UST data for all studies in this toolbox.
The folder 'Neg_47_Left' must be added to the directory:
'.../data/simulation/phantom/OA-BREAST/Neg_47_Left/...'

3) Run ''startup_simulation_ust.m''

4) Run examples in the path ''...r-Wave/simulation/examples-simulation/....''

____________________________________________________________________________

RUNNING EXAMPLES-SIMULATION
____________________________________________________________________________
1) .../simulation/examples-simulation/example_2D_validate_greens_circle.m  (paper [1])
This example validates ray approximation to heterogeneous Green's function in approximating phase and amplitude via
a comparison with the phase and amplitude simulated by the k-Wave. The transducers are placed on a circular ring.

2) .../simulation/examples-simulation/example_2D_image_tof_greens_circle.m (paper [1])
This example reconstructs image of the sound speed using inversion approaches based on ray approaximation to heterogeneous Green's function [1,2].
Early iterations of an image reconstruction algorithm using time-of-flight data will be used for providing initial guess for the Green's inversion approaches.
Note that the Green's inversion approach used in [2] was used in [1] as the benchmark. The transducers are placed on a circular ring.

3) .../simulation/examples-simulation/example_3D_image_tof_sphere.m (paper [3])
This example reconstructs a quantitatively accurate, low resolution and volumetric (full-3D) image of the sound speed
from time-of-flight of ultrasound data simulated on a hemispherical detection surface. The user can compare the images
reconstructed using the bent-ray inversion approach [3] with the same image reconstructed using straight rays.  
('matrix_construction_method= 'bent-ray'), or ('matrix_construction_method= 'straight-ray')

4) .../simulation/examples-simulation/example_fish_eye.m (paper [3])
This example tests ray tracing algorithms on a ❝Maxwell's Fish-eye lens❞ phantom.
All the used ray tracing algorithms in this example can be used for time-of-flight-based reconstruction of the sound speed,
but only 'Runge-kutta-2nd' can be used for ray tracing for the Green's inversion approach. 

5) .../simulation/examples-simulation/example_2D_image_tof_plane.m (paper [1])
This example reconstructs image of the sound speed using inversion approaches based on time-of-flight-based algorithm. The extension to 
the Green's approach for this specific geometry has not been completed, but will be completed soon.

Please read the description of examples for better understanding the examples!
____________________________________________________________________________

SIMULATION OF ULTRASOUND DATA USING THE K-WAVE
____________________________________________________________________________
The pressure field used as the benchmark is simulated using the k-Wave toolbox.  www.k-Wave.org (v. 1.3.) [5].
The k-Wave toolbox and the functions for an off-grid interpolation is available in this project. 
..............................................................................
Notation: It was noticed that correction steps described in [8] must be applied on the k-Wave, either v. 1.3 or v. 1.4,
for getting accurate signals matching the analytic Greens' formula in homogenous media (only water).
Here, instead of applying the required corrections to the k-Wave, the pressure time series simulated by the k-Wave were kept
unchanged, and an inverse of the required correction steps was enforced on the emission pulse which is used as the input to
the Green's formula in homogeneous and heterogeneous media. Taking both approaches are the same, and it will affect only our 
assumption about the emission pulse. (please contact me if you need more information.)

____________________________________________________________________________

RUNNING EXAMPLES-EXPERIMENT
____________________________________________________________________________
Two real ultrasound experiment is also included in this version.
The time-of-flight-based approach is validated using experimental data.

1) .../experiment-real/examples-experiment/example_2d_tof_real_data_circle.m 
This example uses a time-of-flight-based image reconstruction for reconstructing a low resolution
sound speed image from an experimental data. The transducers are placed on the arc of a circular ring.
The ultrasound data (time series) can be downloaded via www.zenodo.org.

2) .../experiment-real/examples-experiment/example_2d_tof_real_data_plane.m 
This example uses a time-of-flight-based image reconstruction for reconstructing a low resolution
sound speed image from a experimental data, which is available online. The transducers are placed on 16 planar
arrays, which are placed on a ring. (So the detection ring is not actually a ring, but a number of line arrays 
which form an approximate ring.) Each array contains 16 transducers. For this specific geometry, the code is run a
bit slowly, but the time issue will be fixed soon. The reason is that for each emitteer, the ray linking 
mus be performed sepsrately for each linear array, which has its own equation.  
The ultrasound data (time series) are not available for me, but the computed time-of-flights can be downloaded via
the open-source link: https://github.com/ucl-bug/ust-sart.


____________________________________________________________________________
ABOUT IMPLEMENTATION OF THE GREENS APPROACH IN EXPERIMNTAL SETTING
____________________________________________________________________________
The Green's approaches have not been tested and have not been applied on any experimental data yet. One reason was my specific
circustance in the last two years of my position at UCL.

Currently, the Green's approach is applied in 2D simulation setting and using point transducers which
are positioned on a ring. But an extension to experimental setting is my next plan.

An extension to experimental setting (the next plan) can also be done by the users of this toolbox, but that needs reading and understanding the details of
the inversion approaches in the mentioned papers and going through the details of the codes. 

If the users are interested in extending the approach to experimental setting before me, the following points should be considered.
1) Because the quantitative recosntruction of the sound speed is also dependent on the harware, a full-wave approach must be used 
as the benchmark.
 
2) The emitters must be changed from points (assumed in the simulation studies) to real transducers with finite area (2D) or volume (3D).

In practice, the area (2D) or volume (3D) of the transducers can be approximated using the ultrasound data set measured in water (only-water data set).
3) The frequency response of the transducers may also be needed. 
4)...

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
