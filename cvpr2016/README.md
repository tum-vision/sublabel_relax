# Sublabel-Accurate Relaxation of Nonconvex Energies (CVPR16)

## Installation
Install prost as described [here](https://github.com/tum-vision/sublabel_relax/blob/master/README.md).

Before running the MATLAB scripts the mex-File for computing convex envelopes has to be compiled. In the directory of this repository, run from within MATLAB the following command:
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
  
Running the sublabel-stereo experiment with 4 labels should then produce the following expected output
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

![Sublabel stereo result](https://github.com/tum-vision/sublabel_relax/raw/master/cvpr2016/images/stereo_result.png)

## Publications
 *   **Sublabel-Accurate Relaxation of Nonconvex Energies**
     (T. Möllenhoff, E. Laude, M. Moeller, J. Lellmann, D. Cremers),
     In IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2016.

 *   **Sublabel-Accurate Convex Relaxation of Vectorial Multilabel Energies**
     (E. Laude, T. Möllenhoff, M. Moeller, J. Lellmann, D. Cremers),
     In European Conference on Computer Vision and Pattern Recognition (ECCV), 2016.



