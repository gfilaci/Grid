    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_trace_pt.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: Gianluca Filaci <g.filaci@ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_MATH_TRACE_PT_H
#define GRID_MATH_TRACE_PT_H
namespace Grid {

//////////////////////////////////////////////////////////////////
// Traces: both all indices and a specific index. Indices must be
// either scalar or matrix
/////////////////////////////////////////////////////////////////

template<class vtype,int N>
  inline auto trace(const iPert<vtype,N> &arg) -> iPert<decltype(trace(arg._internal[0])),N>
{
  iPert<decltype(trace(arg._internal[0])),N> ret;
  for(int i=0;i<N;i++){
    ret._internal[i]=trace(arg._internal[i]);
  }
  return ret;
}


}
#endif
