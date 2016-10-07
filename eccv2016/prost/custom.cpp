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

#include "prox_ind_epi_polyhedral.hpp"

namespace matlab
{

using namespace prost;

template<typename T>
std::vector<T> GetVector(const mxArray *p)
{
  const mwSize *dims = mxGetDimensions(p);
  
  if(dims[1] != 1 && dims[0] != 1)
    throw Exception("Vector has to be Nx1 or 1xN.");

  if(dims[0] == 0 || dims[1] == 0)
    throw Exception("Empty vector passed.");

  if(mxIsDouble(p))
  {
    double *val = mxGetPr(p);
    
    if(dims[1] == 1)
      return std::vector<T>(val, val + dims[0]);
    else
      return std::vector<T>(val, val + dims[1]);
  }
  else if(mxIsSingle(p))
  {
    float *val = (float *)mxGetPr(p);
    
    if(dims[1] == 1)
      return std::vector<T>(val, val + dims[0]);
    else
      return std::vector<T>(val, val + dims[1]);
  }
  else
    throw Exception("Argument has to be passed as a vector of type single or double.");
}

ProxIndEpiPolyhedral<real>*
CreateProxIndEpiPolyhedral(size_t idx, size_t size, bool diagsteps, const mxArray *data)
{
  size_t count = (size_t) mxGetScalar(mxGetCell(data, 0));
  size_t dim = (size_t) mxGetScalar(mxGetCell(data, 1));
  bool interleaved = (bool) mxGetScalar(mxGetCell(data, 2));
  const mxArray *mx_coeffs = mxGetCell(data, 3);

  std::vector<real> coeffs_a = GetVector<real>(mxGetCell(mx_coeffs, 0));
  std::vector<real> coeffs_b = GetVector<real>(mxGetCell(mx_coeffs, 1));
  std::vector<uint32_t> count_vec = GetVector<uint32_t>(mxGetCell(mx_coeffs, 2));
  std::vector<uint32_t> index_vec = GetVector<uint32_t>(mxGetCell(mx_coeffs, 3));

  return new ProxIndEpiPolyhedral<real>(idx, count, dim, interleaved,
    coeffs_a, coeffs_b, count_vec, index_vec);
}

static map<string, function<Prox<real>*(size_t, size_t, bool, const mxArray*)>> custom_prox_reg =
{
  { "ind_epi_polyhedral", CreateProxIndEpiPolyhedral }
};
 
static map<string, function<Block<real>*(size_t, size_t, const mxArray*)>> custom_block_reg = {};

struct RegisterCustom {
  RegisterCustom() {
    get_prox_reg().insert(custom_prox_reg.begin(), custom_prox_reg.end());
    get_block_reg().insert(custom_block_reg.begin(), custom_block_reg.end());
  }
};

static RegisterCustom registerCustom;

} // namespace matlab
