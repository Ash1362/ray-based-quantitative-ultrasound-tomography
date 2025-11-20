____________________________________________________________________________

       r-Wave:
An open-source MATLAB package for quantitative ultrasound tomography
     via ray-Born inversion with in vitro and in vivo validation

____________________________________________________________________________

VERSION INFORMATION
____________________________________________________________________________
Version 1.1, Released 3rd September 2023
Version 1.2, Released 20 November 2025
Written by Ashkan Javaherian

Please report bugs and suggestions to the email: ajavaherian62@gmail.com
The toolbox can be downloaded from https://github.com/Ash1362/ray-based-quantitative-ultrasound-tomography/


_______________________________________________________________________

PUBLICATIONS/PREPRINTS
_______________________________________________________________________

A) SIMULATON STUDIES:
_______________________
The simulation studies are based on the papers:
1 - A. Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution quantitative ultrasound tomography❞, 2022,
https://arxiv.org/abs/2211.00316/.
2 - A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering for biomedical ultrasound tomography❞, Inverse Problems vol. 37, no.11, 115003, 2021.
https://iopscience.iop.org/article/10.1088/1361-6420/ac28ed/.
3- A. Javaherian et al., ❝Refraction-corrected ray-based inversion for three-dimensional ultrasound tomography of the breast❞, Inverse Problems, 36 125010.
https://iopscience.iop.org/article/10.1088/1361-6420/abc0fc/.
4- A. Javaherian, ❝An open-source MATLAB package for quantitative ultrasound tomography via ray-Born inversion with in vitro and in vivo validaton❞, 2025,
https://arxiv.org/abs/


The references [1] and [2] correspond to ray-Born inverson approaches proposed for high-resolution and quantitative reconstruction of the sound speed
images from synthetic transmision ultrasound datasets simulated using the k-Wave toolbox (www.k-wave.org/). The transducers are assumed points 
which are placed on a 2D ring, and the image reconstruction is done in 2D, but the approaches can be extended to 3D.

The reference [3] corresponds to a full-3D and refraction-corrected ray-based algorithm for reconstructing a qunatitatively accurate and
low-resolution image of the sound speed from time-of-flight data. This study was sucessfully implemented on in-vivo ultrasound data
measured from the full-3D pammoth system. The application of the reference [3] on the Pammoth system shows significant contrast between the
bent-ray approach proposed in [3] and straight-ray approach for a fast and quantitative reconstruction of the sound speed from time-of-flight data.

The reference [4] corresponds to the numerical validation of the four ray tracing algorithms developed in this toolbox. 
Using the "Maxwell Fish-eye lens" phantom, the accuracy of ray tracing algorithms are quantified in terms of:
a) Mean Radiaus Deviation from the expected analytic circular or spherical paths for 2d and 3D cases. respectively.
b) Mean Deviation from the expected analytic acoustic lengths.



_______________________________________________________________________

B) EXPERIMENTAL STUDIES:
_______________________
The experimental studies studies are based on the paper:
1 - A. Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution quantitative ultrasound tomography❞, 2022,
https://arxiv.org/abs/2211.00316/.
5- A. Javaherian, ❝The first in vitro and in vivo validation of the hessian-free ray-Born inversion for quantitative ultarsound tomography❞, 2025.

The reference [5] corresponds to the implementation of the Hessian-free ray-Born image reconstruction approach, proposed in reference [1], on in vitro and in vivo
datsets released by the University of Rochester Medical Center [9,12].

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
Download the folder "data_ust_kWave_transmission" via the link:
https://zenodo.org/records/8330926 [7]. 

Add the folder "data_ust_kWave_transmission" to the path:
.../r-Wave/data/simulation/...

This folder includes the ultrasound data simulated by the k-Wave toolbox. Alternatively, the user can set "data_sim = true;" in the
simulation example scripts. By setting this parameter true, the ultrasound data are simulated and are stored in the associated path.

By setting "data_sim=true;", the following point must be considered. 
For the results in the references [1] and [2], the transducers are assumed off-grid points which are placed on a 2D ring.
In the reference [3], which fully corresponds to an image reconstruction using time-of-flight data in full-3D geometry,
the interpolation is done using a neighboring approach.

It must be reminded that that in [1] and [2], early iterations of a time-of-flight-based image reconstruction algorithm from the extracted 
first arrival of the simulated time series are used for providing an initial guess for the proposed ray-Born image reconstruction approaches.


2) A digital breast phantom developed by the team of Professor Mark Anastasio is used in this project. This phantom must be downloaded via the link:...
https://anastasio.bioengineering.illinois.edu/downloadable-content/oa-breast-database/ [11].
The folder 'Neg_47_Left' is used for simulation of UST data for all studies in this toolbox.
The folder 'Neg_47_Left' must be added to the directory:
".../data/simulation/phantom/OA-BREAST/Neg_47_Left/..."

3) Run "startup_simulation_ust.m"

4) Run examples in the path "...r-Wave/simulation/examples-simulation/...."

____________________________________________________________________________

RUNNING EXAMPLES-SIMULATION
____________________________________________________________________________
1) .../simulation/examples-simulation/example_2D_validate_greens_circle.m  (reference [1])
This example validates ray approximation to heterogeneous Green's function in approximating phases and amplitudes via
a comparison with thes phase and amplitudes simulated by the k-Wave. The transducers are placed on a circular ring.
For the k-wave simulations, the incusion of the source has been adapted to a volumtric source \( s \), according
to section 8.1 of the reference [6]. 

2) .../simulation/examples-simulation/example_2D_image_tof_greens_circle.m (reference [1])
The transducers are placed on a circular ring. This example reconstructs image of the sound speed maps using inversion approaches based on ray 
approaximation to heterogeneous Green's function [1,2].
Early iterations of an image reconstruction from time-of-flight data will be used for providing initial guess for the Green's inversion approaches.
Note that the Green's inversion approach used in [2] was used in [1] as the benchmark. In [1], uing both inverion approache, computing the geomterical
portion of the amplitudes of the Green's function was done using the paraxial ray tracing system, as described in the reference [1].
Therefore, the difference betweeen running the Hessian-based and the Hessian-free ray-Born image reconstruction approaches using this example script
arises only in the types of the ray-Born inversion approaches, as described in the references [2] and [1], respectively.

All the developed ray tracing algorithms in this tolbox can be used for time-of-flight-based reconstruction of the sound speed,
but only 'Runge-kutta-2nd' ray tracing algorithm has been adapted to the Green's image reconstruction approaches. 


3) .../simulation/examples-simulation/example_3D_image_tof_sphere.m (reference [3])
This example reconstructs a quantitatively accurate, low resolution and volumetric (full-3D) image of the sound speed
from time-of-flight of ultrasound data simulated on a hemispherical detection surface. The user can compare the images
reconstructed using the bent-ray inversion approach [3] with the same image reconstructed using straight rays.
In this example script, he straigt and bent ray can be chosen by setting:  
'matrix_construction_method= 'bent-ray', or ('matrix_construction_method= 'straight-ray'.


4) .../simulation/examples-simulation/example_validate_ray_fish_eye_phantom.m (reference [4])
This example quantifies accuracy of the developed ray tracing algorithms in terms of two criteria:
a) Mean Radiaus Deviation from the expected analytic circular (or spherical) paths (criterion '2'),
b) Mean Deviation from the expected acoustic lengths (criterion '1').

The ray tracing algorithms are performed on the "Maxwell's Fish-eye lens" phantom, for which the rays' trajectories and
and the accumulated acoustic lengths along the rays can be computed analytically.
All results for the four cases—two criteria each in 2D and 3D—are reported in reference [4].

5) .../simulation/examples-simulation/example_2D_image_tof_plane.m (reference [1])
This example reconstructs image of the sound speed maps using a time-of-flight-based inversion algorithm for transmission ultrasound data simulated on 
a number of rotating linear ultrasund arrays. The extension to the Green's inversion approaches for this specific geometry has not been completed yet.


____________________________________________________________________________

SIMULATION OF ULTRASOUND DATA USING THE K-WAVE
____________________________________________________________________________
The pressure field used as the benchmark is simulated using the k-Wave toolbox.  www.k-Wave.org (Version 1.3.) [10].
The k-Wave toolbox is available in this project [10]. 
..............................................................................
Notation: In the references [1] and [2], the sources are assumed as point sources in terms of \( s\), so 
inclusion of the source in the the k-wave toolbox should be modified to adapt this specific assumption [8]. 






____________________________________________________________________________

B) EXPERIMENTAL STUDIES
____________________________________________________________________________

1) Get access to the UST in-vitro and in-vivo datasets
Download via the link:
https://github.com/rehmanali1994/WaveformInversionUST [12]. 

2) Add the folder "example_rochester" to the path:
.../r-Wave/data/example_rochester/...within our toolbox.

3) Run "startup_simulation_ust.m"

4) Run examples in the path "......./experiment-real/example_urmc/...."   


____________________________________________________________________________

RUNNING EXAMPLES-EXPERIMENT
____________________________________________________________________________
The Green’s image reconstruction approach based on the Hessian-free ray_born inversion has now been successfully
applied to the open-source transmission-ultrasound datasets released by the University of Rochester Medical Center (URMC) [9,12].
This data package, publicly available at the link provided in reference [12], included the following four datasets:
1) VSX_Yezitronix_Phantom1.mat,
2) VSX_Yezitronix_Phantom2.mat,
3) BenignCyst.mat, and
4) Malignancy.mat.

The datasets should be downloaded from the link provided in the reference [12], and added to the path: "..../data/experment_rochester/...."
In the line plots in the reference [5], the images reconstructed using a full-wave apparoach based on the frequency-domain Helholtz equation
are used as the benchmark [9].
For reproducing the line plots in the reference [5], the user must run the script "display_images_urmc_data.m ". Running this script needs 
the images reconstructed using the full-wave approach to be stored in the "..../results/experment_rochester/full-wave/...." before running the
script. The images reconstructed using the full-wave approach can be reproduced through running the scripts in the Github link provided in the 
reference [12].


These examples are included in this toolbox.
1) ..../experiment-real/example_urmc/example_2d_real_data_vsx.m (reference [5])
This example employs the Hessian-free ray-Born inversion approach, proposed in the reference [1], for reconsructing
high-resolution images of the sound speed from transmission ultrasound datasets released by the University of 
Rochester Medical Center [9,12].
The transducers are placed on a ring with a radius of approximately 11cm. For further details, the users are referred to:
Section "IV. RESULTS, C. In-Vitro Experiments" in reference [9].


2) ..../experiment-real/example_urmc/example_2d_real_data_benigncyst.m (rerence [5]) and
3) ..../experiment-real/example_urmc/example_2d_real_data_malignancy.m (reference [5]).
The example scripts 2 and 3 employ the Hessian-free ray-Born inversion approach, proposed in the reference [1], for reconstructing
high-resolution images of the sound speed from transmission ultrasound datasets released by the University of 
Rochester Medical Center [9,12].
The transducers are placed on a ring with a radius of approximately 11cm. For further details, the users are referred to:
Section "IV. RESULTS, D. In-Vivo Demonstration" in reference [9].

In all these three example scripts, the image reconstruction follows the process in simulation experiments, i.e., the example script ''example_2D_image_tof_greens_circle.m''. 
Correspondingly, early iterations (linearizations) of a time-of-based algorithm is used to provide a low-resolution and low-contrast sound speed image
as the initial guess for the main Hessian-free ray-Born image reconstruction algorithm.	

____________________________________________________________________________

ACKNOWLEDGEMENT
__________________________________________________________________________

I would like to express my sincere gratitude to Professor Nebojsa Duric and his team at URMC, as well as to
Delphinus Medical Technologies, for publicly releasing these invaluable datasets. I also gratefully acknowledge
Professor Mohammad Mehrmohammadi, Dr. Rehman Ali, and Mr. Gaofei Jin for their helpful advice and for sharing their
experience with experimental ultrasound data. 

____________________________________________________________________________

REFERENCES
__________________________________________________________________________
Papers/Preprints:
1 - A. Javaherian, ❝Hessian-inversion-free ray-born inversion for high-resolution quantitative ultrasound tomography❞, 2022. https://arxiv.org/abs/2211.00316/ 
2 - A. Javaherian and B. Cox, ❝Ray-based inversion accounting for scattering for biomedical ultrasound tomography❞, Inverse Problems vol. 37, no.11, 115003, 2021. https://iopscience.iop.org/article/10.1088/1361-6420/ac28ed/ 
3 - A. Javaherian, F. Lucka and B. T. Cox, ❝Refraction-corrected ray-based inversion for three-dimensional ultrasound tomography of the breast❞, Inverse Problems, 36 125010. https://iopscience.iop.org/article/10.1088/1361-6420/abc0fc/  
4 - A. Javaherian, ❝An open-source MATLAB package for quantitative ultrasound tomography via ray-Born inversion with in vitro and in vivo validaton❞, 2025.
5 - A. Javaherian, ❝The first in vitro and in vivo validation of the hessian-free ray-Born inversion for quantitative ultarsound tomography❞, 2025.
6 - A. Javaherian and S.K. Setarehdan, ❝Full-waveform approximation of finite-Sized acoustic apertures: forward and adjoint wavefields❝, https://arxiv.org/abs/2212.04466/
7 - A. Javaherian, 2023, ❝Transmission ultrasound data simulated using the k-Wave toolbox as a benchmark for biomedical quantitative ultrasound tomography using a ray approximation to Green's function❞ (1.1) [Data set]. Zenodo. 
https://zenodo.org/records/8330926 
8 - E. S. Wise, B. T. Cox, J. Jaros, B. E. Treeby, ❝Representing arbitrary acoustic source and sensor distributions in Fourier collocation methods❞, J. Acoust. Soc. of Am., vol. 146, no. 1, pp. 278-288, 2019.
9 - R. Ali et al., ❝2-D Slicewise Waveform Inversion of Sound Speed and Acoustic Attenuation for Ring Array Ultrasound Tomography Based on a 
Block LU Solver,❞ in IEEE Transactions on Medical Imaging, vol. 43, no. 8, pp. 2988-3000, Aug. 2024.
Toolboxes/Datasets:
10 - B. E. Treeby and B. T. Cox, ❝k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave fields❞, J. Biomed. Opt. vol. 15, no. 2, 021314, 2010. http://www.k-wave.org/ 
11 - Y. Lou, W. Zhou, T. P. Matthews, C. M. Appleton and M. A. Anastasio, ❝Generation of anatomically realistic numerical phantoms for photoacoustic and ultrasonic breast imaging❞, J. Biomed. Opt., vol. 22, no. 4, pp. 041015, 2017. 
https://anastasio.bioengineering.illinois.edu/downloadable-content/oa-breast-database/
12 - R. Ali, https://github.com/rehmanali1994/WaveformInversionUST (DOI: https://zenodo.org/badge/latestdoi/684631232 )


____________________________________________________________________________
RELEASE NOTES
____________________________________________________________________________
Version 1.1, Released 3rd September 2023
Version 1.2, Released 20 November 2025
Written by Ashkan Javaherian

_______________________________________________________________

LICENSE
____________________________________________________________________________

r-Wave (c) 2022 Ashkan Javaherian
____________________________________________________________________________