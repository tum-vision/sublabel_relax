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

#include "block_dataterm_sublabel.hpp"

namespace prost {

template<typename T>
__global__
void BlockDatatermSublabelKernel(T *d_res,
				 const T *d_rhs,
				 size_t nx,
				 size_t ny,
				 size_t L,
				 T t_min,
				 T t_max)
{
  size_t x = threadIdx.x + blockDim.x * blockIdx.x;
  size_t y_tilde = threadIdx.y + blockDim.y * blockIdx.y;
  size_t l = y_tilde % (L-1);
  size_t y = y_tilde / (L-1);
  
  if(x >= nx || y >= ny || l >= L-1)
    return;
  
  size_t idx = (L-1)*ny*x + y*(L-1) + l;
  
  T delta_t = (t_max - t_min) / (L-1);
  T t = t_min + l * delta_t;
  
  T add_v = -t*d_rhs[idx];
  for(size_t i = 1; i < L-1-l; i++)
  {
    add_v += delta_t*d_rhs[idx+i];
  }

  d_res[idx] += add_v;
}

template<typename T>
__global__
void BlockDatatermSublabelAdjointKernel(T *d_res,
					const T *d_rhs,
					size_t nx,
					size_t ny,
					size_t L,
					T t_min,
					T t_max)
{
  size_t x = threadIdx.x + blockDim.x * blockIdx.x;
  size_t y_tilde = threadIdx.y + blockDim.y * blockIdx.y;
  size_t l = y_tilde % (L-1);
  size_t y = y_tilde / (L-1);

  if(x >= nx || y >= ny || l >= L-1)
    return;

  size_t idx = (L-1)*ny*x + y*(L-1) + l;
  
  T delta_t = (t_max - t_min) / (L-1);
  T t = t_min + l * delta_t;
  
  T add_s = -t*d_rhs[idx];
  for(size_t i = 1; i <= l; i++)
    add_s += delta_t*d_rhs[idx-i];

  d_res[idx] += add_s;
}
  
template<typename T>
BlockDatatermSublabel<T>::BlockDatatermSublabel(size_t row, 
						size_t col, 
						size_t nx, 
						size_t ny, 
						size_t L, 
						T left, 
						T right)
  : Block<T>(row,col,nx*ny*(L-1),nx*ny*(L-1)), nx_(nx), ny_(ny), L_(L), t_min_(left), t_max_(right)
{
}

template<typename T>
T BlockDatatermSublabel<T>::row_sum(size_t row, T alpha) const
{
  size_t l = row % (L_-1);
  T delta_t = (t_max_ - t_min_) / (L_-1);
  T t = t_min_ + l * delta_t;
  return t + (L_-2-l)*delta_t;  
}

template<typename T>
T BlockDatatermSublabel<T>::col_sum(size_t col, T alpha) const
{
  size_t l = col % (L_-1);
  T delta_t = (t_max_ - t_min_) / (L_-1);
  T t = t_min_ + l * delta_t;

  return t + l*delta_t;
}

template<typename T>
void BlockDatatermSublabel<T>::EvalLocalAdd(
    const typename device_vector<T>::iterator& res_begin,
    const typename device_vector<T>::iterator& res_end,
    const typename device_vector<T>::const_iterator& rhs_begin,
    const typename device_vector<T>::const_iterator& rhs_end)
{
  dim3 block(1, 128, 1);
  dim3 grid((nx_ + block.x - 1) / block.x,
            (ny_*(L_-1) + block.y - 1) / block.y,
            1);

  BlockDatatermSublabelKernel<T>
    <<<grid, block>>>(thrust::raw_pointer_cast(&(*res_begin)),
		      thrust::raw_pointer_cast(&(*rhs_begin)),
		      nx_,
		      ny_,
		      L_,
		      t_min_,
		      t_max_);
}

template<typename T>
void BlockDatatermSublabel<T>::EvalAdjointLocalAdd(
    const typename device_vector<T>::iterator& res_begin,
    const typename device_vector<T>::iterator& res_end,
    const typename device_vector<T>::const_iterator& rhs_begin,
    const typename device_vector<T>::const_iterator& rhs_end)
{
  dim3 block(1, 128, 1);
  dim3 grid((nx_ + block.x - 1) / block.x,
            (ny_*(L_-1) + block.y - 1) / block.y,
            1);

  BlockDatatermSublabelAdjointKernel<T>
    <<<grid, block>>>(thrust::raw_pointer_cast(&(*res_begin)),
		      thrust::raw_pointer_cast(&(*rhs_begin)),
		      nx_,
		      ny_,
		      L_,
		      t_min_,
		      t_max_);
}

// Explicit template instantiation
template class BlockDatatermSublabel<float>;
template class BlockDatatermSublabel<double>;

} // namespace prost