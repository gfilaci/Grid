    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_arith_scalar_pt.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef GRID_MATH_ARITH_SCALAR_PT_H
#define GRID_MATH_ARITH_SCALAR_PT_H

namespace Grid {


//////////////////////////////////////////////////////////////////////////////////////////
// Must support native C++ types Integer, Complex, Real
//////////////////////////////////////////////////////////////////////////////////////////

// multiplication by fundamental scalar type
template<class l,int N> strong_inline iPert<l,N> operator * (const iPert<l,N>& lhs,const typename iScalar<l>::scalar_type rhs) 
{
  typename iPert<l,N>::tensor_reduced srhs; srhs=rhs;
  return lhs*srhs;
}
template<class l,int N> strong_inline iPert<l,N> operator * (const typename iScalar<l>::scalar_type lhs,const iPert<l,N>& rhs) {  return rhs*lhs; }

////////////////////////////////////////////////////////////////////
// Double support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////

template<class l,int N> strong_inline iPert<l,N> operator * (const iPert<l,N>& lhs,double rhs) 
{
  typename iScalar<l>::scalar_type t;t=rhs;
  typename iScalar<l>::tensor_reduced srhs;srhs=t;
  return lhs*srhs;
}
template<class l,int N> strong_inline iPert<l,N> operator * (double lhs,const iPert<l,N>& rhs) {  return rhs*lhs; }

////////////////////////////////////////////////////////////////////
// Complex support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////

template<class l,int N> strong_inline iPert<l,N> operator * (const iPert<l,N>& lhs,ComplexD rhs) 
{
  typename iScalar<l>::scalar_type t;t=rhs;
  typename iScalar<l>::tensor_reduced srhs;srhs=t;
  return lhs*srhs;
}
template<class l,int N> strong_inline iPert<l,N> operator * (ComplexD lhs,const iPert<l,N>& rhs) {  return rhs*lhs; }

////////////////////////////////////////////////////////////////////
// Integer support; cast to "scalar_type" through constructor
////////////////////////////////////////////////////////////////////

template<class l,int N> strong_inline iPert<l,N> operator * (const iPert<l,N>& lhs,Integer rhs) 
{
  typename iScalar<l>::scalar_type t;t=rhs;
  typename iScalar<l>::tensor_reduced srhs;srhs=t;
  return lhs*srhs;
}
template<class l,int N> strong_inline iPert<l,N> operator * (Integer lhs,const iPert<l,N>& rhs) {  return rhs*lhs; }

}
#endif
