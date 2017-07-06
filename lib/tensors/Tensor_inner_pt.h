    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_inner_pt.h

    Copyright (C) 2015

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
#ifndef GRID_MATH_INNER_PT_H
#define GRID_MATH_INNER_PT_H
namespace Grid {
  ///////////////////////////////////////////////////////////////////////////////////////
  // innerProduct Perturbative series x Perturbative series -> Scalar
  ///////////////////////////////////////////////////////////////////////////////////////

  template<class l,class r,int N> inline
  auto innerProductD (const iPert<l,N>& lhs,const iPert<r,N>& rhs) -> iScalar<decltype(innerProductD(lhs._internal[0],rhs._internal[0]))>
  {
    typedef decltype(innerProductD(lhs._internal[0],rhs._internal[0])) ret_t;
    iScalar<ret_t> ret;
    ret=zero;
    for(int c1=0;c1<N;c1++){
      ret._internal += innerProductD(lhs._internal[c1],rhs._internal[c1]);
    }
    return ret;
  }
  
  //////////////////////
  // Keep same precison
  //////////////////////
  template<class l,class r,int N> inline
  auto innerProduct (const iPert<l,N>& lhs,const iPert<r,N>& rhs) -> iScalar<decltype(innerProduct(lhs._internal[0],rhs._internal[0]))>
  {
    typedef decltype(innerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
    iScalar<ret_t> ret;
    ret=zero;
    for(int c1=0;c1<N;c1++){
      ret._internal += innerProduct(lhs._internal[c1],rhs._internal[c1]);
    }
    return ret;
  }
  
}
#endif
