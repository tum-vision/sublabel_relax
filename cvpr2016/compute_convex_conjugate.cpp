/**
* This file is part of sublabel_relax.
*
* Copyright 2016 Thomas Möllenhoff <thomas dot moellenhoff at in dot tum dot de> 
* and Emanuel Laude <emanuel dot laude at in dot tum dot de> (Technical University of Munich)
*
* sublabel_relax is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* prost is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with sublabel_relax. If not, see <http://www.gnu.org/licenses/>.
*/

#include "mex.h"
#include <memory.h>
#include <algorithm>
#include <iostream>
#include <list>
#include <cmath>

using std::cout;
using std::endl;

size_t num_wrong = 0;

template<typename T>
std::list<double>
linspace(T start_in, T end_in, int num_in)
{
  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);
  double delta = (end - start) / (num - 1);

  std::list<double> linspaced; 
  for(int i = 0; i < num - 1; ++i) {
    linspaced.push_back(start + delta * i);
  }
  linspaced.push_back(end);
  
  return linspaced;
}

template<typename T>
inline int
Turn(const T px,
     const T py,
     const T qx,
     const T qy,
     const T rx,
     const T ry)
{
  const T det = (qx - px) * (ry - py) - (rx - px) * (qy - py);

  if(det < 1e-6)
    return -1;
  else if(det > 1e-6)
    return 1;

  return 0;
}

template<typename T>
inline T
Dist(const T px, const T py, const T qx, const T qy)
{
  const T dx = qx - px;
  const T dy = qy - py;

  return dx * dx + dy * dy;
}

/// \brief Returns the convex envelope of the function described by the 
///        N points in_x, in_y. Returns the number of points in the envelope,
///        resulting points are writting in cx, cy.
///
///        Assumes in_x[0], in_y[0] is the left most, and in_x[N-1], in_y[N-1]
///        is the rightmost point on the function.
///
template<typename T>
size_t
ConvexEnvelope(T *cx,
               T *cy,
               const T *in_x,
               const T *in_y,
               size_t N)
{
  static const T eps = 1e-12;

  T px, py, epx, epy;
  size_t i = 0;  
  px = in_x[0];
  py = in_y[0];

  do {
    cx[i] = px;
    cy[i] = py;

    // check if we just added the right-most point.
    // note we only need to compare x coordinates, since we are
    // dealing with functions
    if(std::abs(cx[i] - in_x[N - 1]) < eps)
    {
      i++;
      break;
    }

    epx = in_x[0];
    epy = in_y[0];

    for(size_t j = 1; j < N; j++) {
      int t = Turn<double>(cx[i], cy[i], epx, epy, in_x[j], in_y[j]);
      
      if((std::abs(epx - px) < eps) ||
         (t == -1) ||
        ((t == 0) && (Dist<double>(cx[i], cy[i], in_x[j], in_y[j]) >
          Dist<double>(cx[i], cy[i], epx, epy))))
      {
        epx = in_x[j];
        epy = in_y[j];
      }
    }
    
    i++;
    px = epx;
    py = epy;
    
  } while(std::abs(epx - cx[0]) > eps);

  return i;
}

/// \brief Computes the Legendre-Fenchel conjugate of the piecewise linear
///        function described by the points pt_x, pt_y. 
///
template<typename T>
size_t
ConvexConjugate(
  T *out_x,
  T *out_y,
  const T *pt_x,
  const T *pt_y,
  size_t num_pts)
{
  static double eps = 1e-6;

  T new_x, new_y;

  for(size_t i = 0; i < num_pts - 1; i++) {
    T slope = (pt_y[i + 1] - pt_y[i]) / (pt_x[i + 1] - pt_x[i]);

    out_x[i] = slope;
    out_y[i] = slope * pt_x[i] - pt_y[i];
  }

  if(num_pts > 2) {

    // sanity check for conjugate
    for(int i = 0; i < num_pts - 2; i++) {
      T slope_conj = (out_y[i + 1] - out_y[i]) / (out_x[i + 1] - out_x[i]);

      if(std::abs(slope_conj - pt_x[i + 1]) > eps) {
        num_wrong ++;
      }
    }
  }

  return num_pts - 1;
}

/// \brief Given a cost volume of size nx * ny * nz this function first computes
///        the piecewise convex envelope at every pixel. The number of pieces
///        is defined as L - 1. Then it computes the convex conjugate of this
///        piecewise linear function (which is a piecewise linear function again).
///
///        The function then outputs an array containing the points (x, y) on the
///        convex conjugate and an array of size nx * ny * (L - 1) which contains
///        the start index and counts of the convex conjugate points for that label.
///
/// \param the cost volume
/// \param number of desired convex pieces + 1
/// \param sub-label space
/// \param label space
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  const mwSize *dims = mxGetDimensions(prhs[0]);
  double *vals = mxGetPr(prhs[0]);
  size_t L = (size_t) mxGetScalar(prhs[1]);

  size_t ny = dims[0];
  size_t nx = dims[1];
  size_t nz = dims[2];

  plhs[2] = mxCreateNumericMatrix(nx * ny * (L - 1), 1, mxINT32_CLASS, mxREAL);
  plhs[3] = mxCreateNumericMatrix(nx * ny * (L - 1), 1, mxINT32_CLASS, mxREAL);
  plhs[4] = mxCreateDoubleMatrix(L, 1, mxREAL); 

  int *indices = (int*)mxGetPr(plhs[2]);
  int *counts = (int*)mxGetPr(plhs[3]);

  double *sublabels = (double *) mxGetPr(prhs[2]);
  double *labels = (double *) mxGetPr(prhs[3]);

  double *points_x = new double[nz + L];
  double *points_y = new double[nz + L];
  double *out_x = new double[nx*ny*(nz+L)];
  double *out_y = new double[nx*ny*(nz+L)];
  double *conj_x = new double[nx*ny*(nz+L)];
  double *conj_y = new double[nx*ny*(nz+L)];

  size_t num_pts_total = 0;
  for(int x=0;x<nx;x++) {
    for(int y=0;y<ny;y++) {
      
      // loop over all convex envelopes at that pixel

      int j = 0;
      int k = 0;
      for(int i = 0; i < L - 1; i++) {

        while((sublabels[j] < labels[i + 1]) && (j < nz)) {
          points_x[k] = sublabels[j];
          points_y[k] = vals[y + ny * x + j * ny * nx]; // non-convex energy

          k++;
          j++;
        } 

        // add additional point for label 
        if(j < nz)
        {
          if(std::abs(points_x[k] - points_x[k-1]) < 1e-5) { // merge
            points_x[k - 1] = labels[i + 1]; 
          }
          else // add additional label
          {
            points_x[k] = labels[i + 1];

            double alpha = (labels[i + 1] - sublabels[j-1]) / (sublabels[j] - sublabels[j - 1]);
            points_y[k] = (1 - alpha) * vals[y + ny * x + (j - 1) * ny * nx] + 
              alpha * vals[y + ny * x + j * ny * nx];
            
            k++;
          }
        }

        // compute convex envelope
        size_t num_pts;
        num_pts = ConvexEnvelope<double>(&out_x[num_pts_total], &out_y[num_pts_total], points_x, points_y, k);

        // compute convex conjugate of piecewise linear function
        num_pts = ConvexConjugate<double>(&conj_x[num_pts_total],
          &conj_y[num_pts_total],
          &out_x[num_pts_total], 
          &out_y[num_pts_total], 
          num_pts);

        // ordered label first, write indices and counts
        indices[i + y * (L-1) + x * ny * (L-1)] = num_pts_total;
        counts[i + y * (L-1) + x * ny * (L-1)] = num_pts;
        
        num_pts_total += num_pts;

        points_x[0] = points_x[k - 1];
        points_y[0] = points_y[k - 1];
        k = 1;
      }
    }
  }

  cout << "[compute_convex_conjugate] Total number of slopes: " << num_pts_total << endl;
  cout << "[compute_convex_conjugate] Original points: " << nx * ny * nz << " (Reduction factor " << ((double)num_pts_total / (double)(nx * ny * nz)) << ")." << endl;
  if(num_wrong > 0) {
    cout << "[compute_convex_conjugate] Wrong slopes: " << num_wrong << " (" << num_wrong / (double)num_pts_total << ")" << endl;
  }

  plhs[0] = mxCreateDoubleMatrix(num_pts_total, 1, mxREAL); 
  plhs[1] = mxCreateDoubleMatrix(num_pts_total, 1, mxREAL);

  // copy result
  double *pts_x = (double*)mxGetPr(plhs[0]);
  double *pts_y = (double*)mxGetPr(plhs[1]);
  memcpy(pts_x, conj_x, sizeof(double) * num_pts_total);
  memcpy(pts_y, conj_y, sizeof(double) * num_pts_total);

  delete [] points_x;
  delete [] points_y;
  delete [] out_x;
  delete [] out_y;
  delete [] conj_x;
  delete [] conj_y;
}
