    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_arith_mul_pt.h

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
#ifndef GRID_MATH_ARITH_MUL_PT_H
#define GRID_MATH_ARITH_MUL_PT_H

namespace Grid {


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// MUL         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    
template<class rtype,class vtype,class mtype,int N>
strong_inline void mult(iPert<rtype,N> * __restrict__ ret,
                 const iScalar<mtype>   * __restrict__ lhs,
                 const iPert<vtype,N> * __restrict__ rhs){
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1],&lhs->_internal,&rhs->_internal[c1]);
    }
}
template<class rtype,class vtype,class mtype,int N>
strong_inline void mult(iPert<rtype,N> * __restrict__ ret,
                 const iPert<vtype,N> * __restrict__ rhs,
                 const iScalar<mtype> * __restrict__ lhs){
    mult(ret,lhs,rhs);
}
template<class rtype,class vtype,class mtype,int N>
strong_inline void mult(iPert<rtype,N> * __restrict__ ret,const iPert<mtype,N> * __restrict__ lhs,const iPert<vtype,N> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1],&lhs->_internal[0],&rhs->_internal[c1]);
        for(int c2=1;c2<=c1;c2++){
            mac(&ret->_internal[c1],&lhs->_internal[c2],&rhs->_internal[c1-c2]);
        }
    }
    return;
}

//////////////////////////////////////////////////////////////////
// Divide by scalar
//////////////////////////////////////////////////////////////////
template<class rtype,class vtype,int N> strong_inline
iPert<rtype,N> operator / (const iPert<rtype,N>& lhs,const iScalar<vtype>& rhs)
{
    iPert<rtype,N> ret;
    for(int i=0;i<N;i++){
      ret._internal[i] = lhs._internal[i]/rhs._internal;
    }
    return ret;
}
    
    //////////////////////////////////////////////////////////////////
    // Glue operators to mult routines. Must resolve return type cleverly from typeof(internal)
    // since nesting matrix<scalar> x matrix<matrix>-> matrix<matrix>
    // while         matrix<scalar> x matrix<scalar>-> matrix<scalar>
    // so return type depends on argument types in nasty way.
    //////////////////////////////////////////////////////////////////
    // scal x pert = pert
    // pert x scal = pert
    // pert x pert = pert
    //
template<class l,class r,int N> strong_inline
auto operator * (const iScalar<l>& lhs,const iPert<r,N>& rhs) -> iPert<decltype(lhs._internal*rhs._internal[0]),N>
{
    typedef decltype(lhs._internal*rhs._internal[0]) ret_t;
    iPert<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
        mult(&ret._internal[c1],&lhs._internal,&rhs._internal[c1]);
    }
    return ret;
}
template<class l,class r,int N> strong_inline
auto operator * (const iPert<l,N>& lhs,const iScalar<r>& rhs) -> iPert<decltype(lhs._internal[0]*rhs._internal),N>
{
    typedef decltype(lhs._internal[0]*rhs._internal) ret_t;
    iPert<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
        mult(&ret._internal[c1],&lhs._internal[c1],&rhs._internal);
    }
    return ret;
}
template<class l,class r,int N> strong_inline
auto operator * (const iPert<l,N>& lhs,const iPert<r,N>& rhs) -> iPert<decltype(lhs._internal[0]*rhs._internal[0]),N>
{
    typedef decltype(lhs._internal[0]*rhs._internal[0]) ret_t;
    iPert<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
        mult(&ret._internal[c1],&lhs._internal[0],&rhs._internal[c1]);
         for(int c2=1;c2<=c1;c2++){
            mac(&ret._internal[c1],&lhs._internal[c2],&rhs._internal[c1-c2]);
        }
    }
    return ret;
}

}
#endif
