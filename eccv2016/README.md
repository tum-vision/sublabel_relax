# Sublabel-Accurate Convex Relaxations

## Installation
Add the folder `prost/matlab` to your MATLAB path, e.g. by writing
```
addpath('~/Documents/Projects/prost/matlab')
```
from within MATLAB.
## Usage

The code for reproducing individual numerical experiments from the [paper](https://vision.in.tum.de/_media/spezial/bib/moellenhoff_laude_cvpr_16.pdf) are organized in different files.

**ROF denoising of a vector-valued signal, Figure 4:** `main_spirale_rof_2d.m` (Direct + Sublabel)

**ROF denoising of a color-image, Figure 5:** `main_rof_3d.m` (Direct + Sublabel)

**Denoising of a color-image with a robust non-convex dataterm, Figure 6:** `main_sublabel_robust_rof_3d.m`

**Optical flow, Figure 6:** `main_sublabel_flow.m`

