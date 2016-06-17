# Sublabel-Accurate Convex Relaxations

## Installation
To run the code, first install the convex optimization framework [prost](https://github.com/tum-vision/prost) following the instructions presented there. 

It is required to add some additional proximal and linear operators to `prost`. To do so, create the file `CustomSources.cmake` in the directory `prost/cmake/` with the following contents:

```
set(PROST_CUSTOM_SOURCES
  "relative_path_to_sublabel_relax"/cvpr2016/prost/block_dataterm_sublabel.cu
  "relative_path_to_sublabel_relax"/cvpr2016/prost/prox_ind_epi_polyhedral_1d.cu
  "relative_path_to_sublabel_relax"/cvpr2016/prost/prox_ind_epi_conjquad_1d.cu
  
  "relative_path_to_sublabel_relax"/cvpr2016/prost/block_dataterm_sublabel.hpp
  "relative_path_to_sublabel_relax"/cvpr2016/prost/prox_ind_epi_polyhedral_1d.hpp
  "relative_path_to_sublabel_relax"/cvpr2016/prost/prox_ind_epi_conjquad_1d.hpp
  )
  
set(MATLAB_CUSTOM_SOURCES
  "relative_path_to_sublabel_relax"/cvpr2016/prost/custom.cpp
  )
```

Replace `"relative_path_to_sublabel_relax"` with the relative path to go from the directory `prost/cmake` to the directory where you cloned this repository into (e.g., `../../sublabel_relax`).

After adding this file, recompile `prost` again, e.g., run in the directory `prost/build`
```
cmake ..
make -j16
```

Finally, before running the MATLAB scripts the mex-File for computing convex envelopes has to be compiled. In the directory of this repository, run from within MATLAB the following command:
```
mex compute_convex_conjugate.cpp
```

Finally, add the folder `prost/matlab` to your MATLAB path, e.g. by writing
```
addpath('~/Documents/Projects/prost/matlab')
```
from within MATLAB.
## Usage

The code for reproducing individual numerical experiments from the [paper](https://vision.in.tum.de/_media/spezial/bib/moellenhoff_laude_cvpr_16.pdf) are organized in different files.

**ROF denoising, Figure 5:**
 * `rof_direct.m` Direct optimization of the ROF model without functional lifting
 * `rof_baseline.m` Optimization of the ROF model with the baseline approach
 * `rof_sublabel.m` Sublabel-accurate optimization of the ROF model

**Denoising with robust truncated quadratic dataterm, Figure 6:**
 * `truncrof_baseline.m` Optimzation using the baseline approach
 * `truncrof_sublabel.m` Sublabel-accurate version
 
**Stereo matching, Figure 9:**
  * `stereo_baseline.m` Baseline approach for stereo matching
  * `stereo_sublabel.m` Sublabel-accurate approach
  
Running the sublabel-stereo experiment (with `L=4`) should then produce the following expected output
``` 
>> stereo_sublabel
[compute_convex_conjugate] Total number of slopes: 5025269
[compute_convex_conjugate] Original points: 50017500 (Reduction factor 1.00e-01).
prost v0.2-build-2016-06-17
Running on device number 0: GeForce GTX TITAN X (12.0 GB, 3072 Cores).
# primal variables: 2223000
# dual variables: 4816500
Memory requirements: 255MB (11858/12287MB available).
It     1: Feas_p=6.09e+02, Eps_p=2.80e-02, Feas_d=0.00e+00, Eps_d=1.49e-02; 
It  1043: Feas_p=3.21e+00, Eps_p=4.09e-02, Feas_d=1.76e-02, Eps_d=1.49e-02; 
It  2085: Feas_p=7.98e-01, Eps_p=4.09e-02, Feas_d=1.68e-02, Eps_d=1.49e-02; 
It  3126: Feas_p=2.25e-01, Eps_p=4.09e-02, Feas_d=1.63e-02, Eps_d=1.49e-02; 
It  4168: Feas_p=7.54e-02, Eps_p=4.09e-02, Feas_d=1.51e-02, Eps_d=1.49e-02; 
It  4825: Feas_p=3.97e-02, Eps_p=4.09e-02, Feas_d=1.49e-02, Eps_d=1.49e-02; 
Reached convergence tolerance.
Elapsed time is 24.236633 seconds.
```

![alt text](https://github.com/tum-vision/sublabel_relax/raw/master/images/stereo_result.png "Expected output for 4 labels.")

## Publications

The following publications describe the approach:

 *   **Sublabel-Accurate Relaxation of Nonconvex Energies**
     (T. Möllenhoff, E. Laude, M. Moeller, J. Lellmann, D. Cremers),
     In IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2016.
