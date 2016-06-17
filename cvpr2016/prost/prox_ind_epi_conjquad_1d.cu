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

#include "prox_ind_epi_conjquad_1d.hpp"

#include "prost/prox/prox_elem_operation.hpp"
#include "prost/prox/prox_elem_operation.inl"

namespace prost {

template class ProxElemOperation<float, ElemOperationIndEpiConjQuad1D<float>>;
template class ProxElemOperation<double, ElemOperationIndEpiConjQuad1D<double>>;
  
} // namespace prost
