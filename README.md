# Sublabel-Accurate Convex Relaxations

## Installation
To run the code, first install the convex optimization framework [prost](https://github.com/tum-vision/prost) following the instructions presented there. 

It is required to add some additional proximal and linear operators to [prost](https://github.com/tum-vision/prost). To do so, create the file `CustomSources.cmake` in the directory `prost/cmake/` with the following contents:

```
set(PROST_CUSTOM_SOURCES
  "relative_path_to_sublabel_relax"/cvpr2016/prost/block_dataterm_sublabel.cu
  "relative_path_to_sublabel_relax"/cvpr2016/prost/prox_ind_epi_polyhedral_1d.cu
  "relative_path_to_sublabel_relax"/cvpr2016/prost/prox_ind_epi_conjquad_1d.cu
  "relative_path_to_sublabel_relax"/eccv2016/prost/prox_ind_epi_polyhedral.cu
  
  "relative_path_to_sublabel_relax"/cvpr2016/prost/block_dataterm_sublabel.hpp
  "relative_path_to_sublabel_relax"/cvpr2016/prost/prox_ind_epi_polyhedral_1d.hpp
  "relative_path_to_sublabel_relax"/cvpr2016/prost/prox_ind_epi_conjquad_1d.hpp
  "relative_path_to_sublabel_relax"/eccv2016/prost/prox_ind_epi_polyhedral.hpp
  )
  
set(MATLAB_CUSTOM_SOURCES
  "relative_path_to_sublabel_relax"/cvpr2016/prost/custom.cpp
  "relative_path_to_sublabel_relax"/eccv2016/prost/custom.cpp
  )
```

Replace `"relative_path_to_sublabel_relax"` with the relative path to go from the directory `prost/cmake` to the directory where you cloned this repository into (e.g., `../../sublabel_relax`).

After adding this file, recompile prost again, e.g., run in the directory `prost/build`
```
cmake ..
make -j16
```


## One-dimensional ranges (CVPR16)
Further instructions for reproducing the individual numerical experiments from the paper [Sublabel-Accurate Relaxation of Nonconvex Energies](https://vision.in.tum.de/_media/spezial/bib/moellenhoff_laude_cvpr_16.pdf) can be found [here](tree/master/cvpr2016/README.md).
## Multi-dimensional ranges (ECCV16)
Further instructions for reproducing the individual numerical experiments from the paper [Sublabel-Accurate Convex Relaxation of Vectorial Multilabel Energies](https://vision.in.tum.de/_media/spezial/bib/laude16eccv.pdf) can be found [here](tree/master/eccv2016/README.md).


## Publications
 *   **Sublabel-Accurate Relaxation of Nonconvex Energies**
     (T. Möllenhoff, E. Laude, M. Moeller, J. Lellmann, D. Cremers),
     In IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2016.

 *   **Sublabel-Accurate Convex Relaxation of Vectorial Multilabel Energies**
     (E. Laude, T. Möllenhoff, M. Moeller, J. Lellmann, D. Cremers),
     In European Conference on Computer Vision and Pattern Recognition (ECCV), 2016.
