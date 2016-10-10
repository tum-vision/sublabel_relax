#include "mex.h"
#include "quickhull/QuickHull.hpp"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>


using namespace quickhull;

inline size_t compute_convex_envelope_2d(const std::vector<Vector3<double>>& points, std::vector<double>& verts, std::vector<double>& vals) {
    QuickHull<double> qh;

    ConvexHull<double> h = qh.getConvexHull(points, true, true);
    std::vector<size_t>& indexBuffer = h.getIndexBuffer();
    
    // remember vertices already visited when looping over triangle list of convex hull
    std::vector<bool> visited(points.size(), false);
    //VertexDataSource<double> vertexBuffer = h.getVertexBuffer();
    // mexPrintf("Convex hull contains %d points\n", vertexBuffer.size());
    
    // cut convex hull at the plane defined by
    // the 3 vertices of the current triangle.
    // By defintion those are the first 3 vertices in points.
    const Vector3<double>& t1 = points[0];
    const Vector3<double>& t2 = points[1];
    const Vector3<double>& t3 = points[2];
    
    // construct plane
    Vector3<double> w1 = t3 - t1;
    Vector3<double> w2 = t2 - t1;
    
    Vector3<double> normal = w1.crossProduct(w2).getNormalized();
    if(normal.z > 0)
        normal = -normal;
    Plane<double> plane(normal, t1);
    
    verts.push_back(t1.x);
    verts.push_back(t1.y);
    vals.push_back(t1.z);
    
    verts.push_back(t2.x);
    verts.push_back(t2.y);
    vals.push_back(t2.z);
    
    verts.push_back(t3.x);
    verts.push_back(t3.y);
    vals.push_back(t3.z);
    
    visited[0] = true;
    visited[1] = true;
    visited[2] = true;
    
    size_t count = 3;
    
    
    for(size_t v : indexBuffer) {
        if(visited.at(v)) {
            continue;
        }
        
        if(plane.isPointOnPositiveSide(points.at(v))) {
            verts.push_back(points.at(v).x);
            verts.push_back(points.at(v).y);
            vals.push_back(points.at(v).z);

            count++;
        }
        visited.at(v) = true;
    }
    return count;
    return count;
}

inline double interp_bilinear(double v1, double v2, double v3, double v4, double alpha, double beta) {
    double w1 = alpha * v1 + (1-alpha) * v2;
    double w2 = alpha * v3 + (1-alpha) * v4;
    return beta * w1 + (1-beta) * w2;
}

inline double compute_cost_flow(std::vector<double>& im0, std::vector<double>& im1, int nx, int ny, int j, int i, int width, double x, double y) {
  double cost = 0.;

  // for a window of size 31 we have
  // x=0..30
  // x-(w-1)/2.=-15..15
  double px = j + x - (width-1)/2.;
  double py = i + y - (width-1)/2.;

  int px_low = std::min(std::max((int)std::floor(px), 0), nx-1);
  int py_low = std::min(std::max((int)std::floor(py), 0), ny-1);
  int px_high = std::min(std::max((int)std::ceil(px), 0), nx-1);
  int py_high = std::min(std::max((int)std::ceil(py), 0), ny-1);


  for(int k = 0; k < 3; k++) 
  {
    double c0 = 0.;
    double c1 = 0.;
    try {
      c0 = im0.at(k*nx*ny + ny*j + i);
      c1 = interp_bilinear(
        im1.at(k*nx*ny + ny*px_low + py_low), 
        im1.at(k*nx*ny + ny*px_low + py_high), 
        im1.at(k*nx*ny + ny*px_high + py_low),
        im1.at(k*nx*ny + ny*px_high + py_high),
        std::ceil(y) - y,
        std::ceil(x) - x
        );
    } catch (const std::out_of_range& oor) {
      std::cerr << "Out of Range error: " << oor.what() << '\n';
    }
    //cost += std::abs(c0 - c1);
    cost += (c0 - c1) * (c0 - c1);
  }
 
  //return cost;
  return sqrt(cost);
}


inline void triangulate(double* vert, int* tri, size_t l, double delta) {
    for(size_t k = 0; k < l; k++) {
        for(size_t m = 0; m < l; m++) {
            vert[0*l*l + k*l + m] = k * delta;
            vert[1*l*l + k*l + m] = m * delta;
        }
    }

    for(size_t k = 0; k < l-1; k++) {
        for(size_t m = 0; m < l-1; m++) {
            size_t idx = k*l+m+1;

            // lower
            tri[k*2*(l-1) + 2*m]                   = idx;
            tri[(l-1)*(l-1)*2 + k*2*(l-1) + 2*m]   = idx+1;
            tri[2*(l-1)*(l-1)*2 + k*2*(l-1) + 2*m] = idx+l+1;
            
            // upper
            tri[k*2*(l-1) + 2*m + 1]                   = idx+l+1;
            tri[(l-1)*(l-1)*2 + k*2*(l-1) + 2*m + 1]   = idx+l;
            tri[2*(l-1)*(l-1)*2 + k*2*(l-1) + 2*m + 1] = idx;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mwSize *dims = mxGetDimensions(prhs[0]);
    
    size_t ny = dims[0];
    size_t nx = dims[1];
    size_t c = dims[2];
    
    std::vector<double> im0(mxGetPr(prhs[0]), mxGetPr(prhs[0])+nx*ny*c);
    std::vector<double> im1(mxGetPr(prhs[1]), mxGetPr(prhs[1])+nx*ny*c);

    // width of the flow vector search window
    size_t width = (size_t) mxGetScalar(prhs[2]);
    
    // number of labels
    size_t l = (size_t) mxGetScalar(prhs[3]);
    size_t L = l*l;
    
    // number of sublabels 
    size_t s  = (size_t) mxGetScalar(prhs[4]);
    
    // vert
    plhs[0] = mxCreateDoubleMatrix(L, 2, mxREAL);
    double *vert = (double*)mxGetPr(plhs[0]);
    
    double delta = (width-1.) / (l-1.);
    
    // tri
    plhs[1] = mxCreateNumericMatrix(2*(l-1)*(l-1), 3, mxINT32_CLASS, mxREAL);
    int *tri = (int*)mxGetPr(plhs[1]);
    triangulate(vert, tri, l, delta);

    bool output_costvolume = false;

    if(nrhs > 5)
      output_costvolume =  (bool) mxGetScalar(prhs[5]);

    // output a volume of size ny * nx * l * l 
    if(output_costvolume)
    {
      std::cout << "Computing cost volume... ";
      std::cout.flush();
      
      mwSize dims[4];
      dims[0] = ny;
      dims[1] = nx;
      dims[2] = l;
      dims[3] = l;

      plhs[2] = mxCreateNumericArray(4, dims, mxSINGLE_CLASS, mxREAL);

      float *cost_volume = (float *)mxGetPr(plhs[2]);
      for(int l2 = 0; l2 < l; l2++)
      {
        for(int l1 = 0; l1 < l; l1++)
        {
          double dx = vert[1*l*l + l2 * l + l1];
          double dy = vert[0*l*l + l2 * l + l1];
          
          for(int j = 0; j < nx; j++)
          {
            for(int i = 0; i < ny; i++)
            {
              cost_volume[i + j * ny + l1 * ny * nx + l2 * ny * nx * l] =
                compute_cost_flow(im0, im1, nx, ny, j, i, width, dx, dy);
            }
          }
        }
      }

      std::cout << " done!\n";
      
      return;
    }
                    
    size_t N = nx*ny;
    
    // index
    plhs[2] = mxCreateNumericMatrix((l-1)*(l-1)*2*N, 1, mxINT32_CLASS, mxREAL);
    int *index = (int*)mxGetPr(plhs[2]);
    
    // counts
    plhs[3] = mxCreateNumericMatrix((l-1)*(l-1)*2*N, 1, mxINT32_CLASS, mxREAL);
    int *counts = (int*)mxGetPr(plhs[3]);

    size_t total_sublabels = 0;
    size_t total_triangles = 0; 
    
    std::vector<double> sublabels;
    std::vector<double> cost;
    
    std::cout << "Convexifiying data... ";
    std::cout.flush();
    
    // loop over all squares
    for(size_t k = 0; k < l-1; k++) {

        for(size_t m = 0; m < l-1; m++) {

            // loop over all squares and collect sublabels belonging to
            // either one or both of the triangles lower and upper
            std::vector<Vector3<double>> curr_triangle[2];

            for(int a = 0; a < 2; a++) {
                for(int b = 0; b < 2; b++) {
                    double x = (k+a)*delta;
                    double y = (m+b)*delta;
                    Vector3<double> v(x, y, 0.);
                    if(b >= a) {
                        curr_triangle[0].push_back(v);
                    }
                    if(b <= a) {
                        curr_triangle[1].push_back(v);
                    }
                }
            }




            for(int a = 0; a < s; a++) {
                for(int b = 0; b < s; b++) {
                    if((a == 0 && b == 0) || (a == s-1 && b == s-1) || (a == 0 && b == s-1) || (a == s-1 && b == 0))
                        continue;
                    
                    double x = delta*(k + a / (s-1.));
                    double y = delta*(m + b / (s-1.));
                    Vector3<double> v(x, y, 0.);


                    if(b >= a) {    
                        curr_triangle[0].push_back(v);
                    }
                    if(b <= a) {
                        curr_triangle[1].push_back(v);
                    }


                }
            }


            // compute convex envelope on lower and upper triangle at each pixel
            for(int g = 0; g < 2; g++) {
                for(size_t j = 0; j < nx; j++) {
                    for(size_t i = 0; i < ny; i++) {
                        // compute costs
                        for(auto& vec : curr_triangle[g]) {
                            vec.z = compute_cost_flow(im0, im1, nx, ny, j, i, width, vec.x, vec.y);
                        }

                        size_t count = compute_convex_envelope_2d(curr_triangle[g], sublabels, cost);

                        index[total_triangles] = total_sublabels;
                        total_sublabels += count;
                        counts[total_triangles] = count;
                        total_triangles++;
                    }
                }
            }
        }
    }
    std::cout << " done!\n";
    
    // sublabels
    plhs[4] = mxCreateDoubleMatrix(sublabels.size(), 1, mxREAL);
    double *pts_x = (double*)mxGetPr(plhs[4]);
    
    std::copy(sublabels.begin(), sublabels.end(), pts_x);
    
    // cost
    plhs[5] = mxCreateDoubleMatrix(cost.size(), 1, mxREAL);
    double *pts_y = (double*)mxGetPr(plhs[5]);
    
    std::copy(cost.begin(), cost.end(), pts_y);
}
