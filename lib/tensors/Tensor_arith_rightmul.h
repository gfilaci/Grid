    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_arith_rightmul.h

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
#ifndef GRID_MATH_ARITH_RIGHTMUL_H
#define GRID_MATH_ARITH_RIGHTMUL_H

namespace Grid {
    
    //////////////////////////////////////////////////////////////////
    // Matrix right multiplies vector (missing in non _pt files)
    //////////////////////////////////////////////////////////////////
    // vec x mat = vec

template<class rrtype,class ltype,class rtype,int N>
strong_inline void mac(iVector<rrtype,N> * __restrict__ ret,const iVector<ltype,N> * __restrict__ lhs,const iMatrix<rtype,N> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        mac(&ret->_internal[c1],&lhs->_internal[c2],&rhs->_internal[c2][c1]);
    }}
    return;
}

template<class rtype,class vtype,class mtype,int N>
strong_inline void mult(iVector<rtype,N> * __restrict__ ret,const iVector<mtype,N> * __restrict__ lhs,const iMatrix<vtype,N> * __restrict__ rhs)
{
    for(int c1=0;c1<N;c1++){
        mult(&ret->_internal[c1],&lhs->_internal[0],&rhs->_internal[0][c1]);
        for(int c2=1;c2<N;c2++){
            mac(&ret->_internal[c1],&lhs->_internal[c2],&rhs->_internal[c2][c1]);
        }
    }
    return;
}

template<class l,class r,int N> strong_inline
auto operator * (const iVector<l,N>& lhs,const iMatrix<r,N>& rhs) -> iVector<decltype(lhs._internal[0]*rhs._internal[0][0]),N>
{
    typedef decltype(lhs._internal[0]*rhs._internal[0][0]) ret_t;
    iVector<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
        mult(&ret._internal[c1],&lhs._internal[0],&rhs._internal[0][c1]);
        for(int c2=1;c2<N;c2++){
            mac(&ret._internal[c1],&lhs._internal[c2],&rhs._internal[c2][c1]);
        }
    }
    return ret;
}

}
#endif
