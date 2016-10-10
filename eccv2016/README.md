# Sublabel-Accurate Convex Relaxation of Vectorial Multilabel Energies (ECCV16)

## Installation
Add the folder `prost/matlab` to your MATLAB path, e.g. by writing
```
addpath('~/Documents/Projects/prost/matlab')
```
from within MATLAB.

For optical flow you need to compile the file `compute_cost_flow_piecw_conv.cpp` that first computes the optical flow matching cost given a pair of images and second, convexifies the cost on each triangle. To compile execute the command within MATLAB:
```
mex compute_cost_flow_piecw_conv.cpp quickhull/QuickHull.cpp
```
The code uses the Quickhull algorithm to compute convex Hulls of a three-dimensional point-cloud. The implementation is borrowed from [https://github.com/akuukka/quickhull](https://github.com/akuukka/quickhull)

## Usage

* **ROF denoising of a vector-valued signal, Figure 4:** `main_spirale_rof_2d.m` (Direct + Sublabel)

* **ROF denoising of a color-image, Figure 5:** `main_rof_3d.m` (Direct + Sublabel)

* **Denoising of a color-image with a robust non-convex dataterm, Figure 6:** `main_sublabel_robust_rof_3d.m`

* **Optical flow, Figure 7:** `main_sublabel_flow.m`

