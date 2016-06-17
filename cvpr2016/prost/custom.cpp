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

#include <map>
#include <string>

#include "factory.hpp"

#include "prost/prox/prox_elem_operation.hpp"
#include "prox_ind_epi_polyhedral_1d.hpp"
#include "prox_ind_epi_conjquad_1d.hpp"
#include "block_dataterm_sublabel.hpp"

namespace matlab
{

using namespace prost;

ProxElemOperation<real, ElemOperationIndEpiConjQuad1D<real> >*
CreateProxIndEpiConjQuad1D(size_t idx, size_t size, bool diagsteps, const mxArray *data)
{
  size_t count = (size_t) mxGetScalar(mxGetCell(data, 0));
  size_t dim = (size_t) mxGetScalar(mxGetCell(data, 1));
  bool interleaved = (bool) mxGetScalar(mxGetCell(data, 2));
  const mxArray *mx_coeffs = mxGetCell(data, 3);

  std::array<std::vector<real>, 5> coeffs;

  for(int i = 0; i < 5; i++) {
    const mwSize *dims = mxGetDimensions(mxGetCell(mx_coeffs, i));
    double *val = mxGetPr(mxGetCell(mx_coeffs, i));
    coeffs[i] = std::vector<real>(val, val + dims[0]);
  }

  return new ProxElemOperation<real, ElemOperationIndEpiConjQuad1D<real> >(
    idx, count, dim, interleaved, diagsteps, coeffs);
}

ProxIndEpiPolyhedral1D<real>*
CreateProxIndEpiPolyhedral1D(size_t idx, size_t size, bool diagsteps, const mxArray *data)
{
  size_t count = (size_t) mxGetScalar(mxGetCell(data, 0));
  size_t dim = (size_t) mxGetScalar(mxGetCell(data, 1));
  bool interleaved = (bool) mxGetScalar(mxGetCell(data, 2));
  const mxArray *mx_coeffs = mxGetCell(data, 3);

  std::array<std::vector<real>, 4> coeffs_xyab;
  std::array<std::vector<size_t>, 2> coeffs_ci;

  for(int i = 0; i < 4; i++)
  {
    const mwSize *dims = mxGetDimensions(mxGetCell(mx_coeffs, i));
    double *val = mxGetPr(mxGetCell(mx_coeffs, i));
    coeffs_xyab[i] = std::vector<real>(val, val + dims[0]);
  }

  for(int i = 0; i < 2; i++)
  {
    const mwSize *dims = mxGetDimensions(mxGetCell(mx_coeffs, 4 + i));
    double *val = mxGetPr(mxGetCell(mx_coeffs, 4 + i));
    coeffs_ci[i] = std::vector<size_t>(val, val + dims[0]);
  }
  
  return new ProxIndEpiPolyhedral1D<real>(idx, count, interleaved,
					  coeffs_xyab[0],
					  coeffs_xyab[1],
					  coeffs_xyab[2],
					  coeffs_xyab[3],
					  coeffs_ci[0],
					  coeffs_ci[1]);
}

BlockDatatermSublabel<real>*
CreateBlockDatatermSublabel(size_t row, size_t col, const mxArray *data)
{
  size_t nx = (size_t) mxGetScalar(mxGetCell(data, 0));
  size_t ny = (size_t) mxGetScalar(mxGetCell(data, 1));
  size_t L = (size_t) mxGetScalar(mxGetCell(data, 2));
  real left = (real) mxGetScalar(mxGetCell(data, 3));
  real right = (real) mxGetScalar(mxGetCell(data, 4));

  return new BlockDatatermSublabel<real>(row, col, nx, ny, L, left, right);    
}

static map<string, function<Prox<real>*(size_t, size_t, bool, const mxArray*)>> custom_prox_reg =
{
  { "ind_epi_conjquad_1d",   CreateProxIndEpiConjQuad1D   },
  { "ind_epi_polyhedral_1d", CreateProxIndEpiPolyhedral1D },
};
 
static map<string, function<Block<real>*(size_t, size_t, const mxArray*)>> custom_block_reg =
{
  { "dataterm_sublabel", CreateBlockDatatermSublabel },
};

struct RegisterCustom {
  RegisterCustom() {
    get_prox_reg().insert(custom_prox_reg.begin(), custom_prox_reg.end());
    get_block_reg().insert(custom_block_reg.begin(), custom_block_reg.end());
  }
};

static RegisterCustom registerCustom;

} // namespace matlab
