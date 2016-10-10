# Sublabel-Accurate Convex Relaxation of Vectorial Multilabel Energies (ECCV16)

## Installation
Install prost as described [here](https://github.com/tum-vision/sublabel_relax/blob/master/README.md).

Add the folder `prost/matlab` to your MATLAB path, e.g. by writing
```
addpath('~/Documents/Projects/prost/matlab')
```
from within MATLAB.

For optical flow you need to compile the file `compute_cost_flow_piecw_conv.cpp` that first computes the optical flow matching cost given a pair of images and second, convexifies the cost on each triangle. To compile execute the command within MATLAB:
```
mex compute_cost_flow_piecw_conv.cpp quickhull/QuickHull.cpp
```
The code uses the Quickhull algorithm to compute convex hulls of three-dimensional point-clouds. The implementation is borrowed from [https://github.com/akuukka/quickhull](https://github.com/akuukka/quickhull).

For visualization of the optical flow results you need to download the [flow-code-matlab.zip](http://vision.middlebury.edu/flow/code/flow-code-matlab.zip). Extraxt the content and add the folder to your path via
```
addpath('~/flow-code-matlab')
```
Finally download the Middlebury optical flow benchmark data [other-color-allframes.zip](http://vision.middlebury.edu/flow/data/comp/zip/other-color-allframes.zip).
## Usage
The following files contain the code to reproduce the numerical experiments as described in the paper [Sublabel-Accurate Convex Relaxation of Vectorial Multilabel Energies](https://vision.in.tum.de/_media/spezial/bib/laude16eccv.pdf). Remark: Only the Lellmann'13 baselines are included.

* **ROF denoising of a vector-valued signal, Figure 4:** `main_spirale_rof_2d.m` (Direct + Sublabel)

* **ROF denoising of a color-image, Figure 5:** `main_rof_3d.m` (Direct + Sublabel)

* **Denoising of a color-image with a robust non-convex dataterm, Figure 6:** `main_sublabel_robust_rof_3d.m`, `main_baseline_robust_rof_3d.m`

* **Optical flow, Figure 7:** `main_sublabel_flow.m`, `main_baseline_flow.m`

## Publications
 *   **Sublabel-Accurate Relaxation of Nonconvex Energies**
     (T. Möllenhoff, E. Laude, M. Moeller, J. Lellmann, D. Cremers),
     In IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2016.

 *   **Sublabel-Accurate Convex Relaxation of Vectorial Multilabel Energies**
     (E. Laude, T. Möllenhoff, M. Moeller, J. Lellmann, D. Cremers),
     In European Conference on Computer Vision and Pattern Recognition (ECCV), 2016.